/*
 * sim.c  --  the simulation engine (Srikripa: correctness / base kernel,
 *                                    Gautami : opt kernel / SIMD)
 *
 * Model (per UNBURNT cell c with burning neighbours n):
 *     p_n = p0 * (1 + fuel[c]) * exp(wind . d_n) * exp(cs * tan(theta_n->c))
 *     P   = 1 - prod_n (1 - p_n)              (combined probability)
 *     c ignites if hash_rand(seed,i,j,t) < P
 * BURNING cells count down a timer, then become BURNED.
 *
 * PULL-STYLE UPDATE: a cell only READS neighbours from R and only WRITES its
 * own entry in W (and its own timer).  Therefore there are no data races and
 * no atomics / critical sections are needed.
 */
#include "wildfire.h"
#include <math.h>
#include <omp.h>

/* ---------------------------------------------------------------------- */
/* Cell-level logic, shared by ALL kernels so results are bit-identical.   */
/* ---------------------------------------------------------------------- */

/* Does UNBURNT cell (i,j) ignite at step t ? */
static inline int ignites(const Sim *s, const uint8_t *R, int idx, int i, int j, int t)
{
    double p_none = 1.0;
    int any = 0;
    const double base = s->p0 * (1.0 + (double)s->fuel[idx]);
    const double e_c  = (double)s->elev[idx];

    for (int d = 0; d < 8; d++) {
        int n = idx + s->off[d];
        if (R[n] == BURNING) {
            /* tan(theta): positive when the fire climbs uphill (n -> c) */
            double tan_th = (e_c - (double)s->elev[n]) * s->inv_dist[d];
            double p = base * s->wind_fac[d] * exp(s->cs * tan_th);
            if (p > 0.999) p = 0.999;
            p_none *= (1.0 - p);
            any = 1;
        }
    }
    if (!any) return 0;                       /* no random number needed */
    return hash_rand(s->seed, i, j, t) < (1.0 - p_none);
}

static inline void update_unburnt(const Sim *s, const uint8_t *R, uint8_t *W,
                                  uint8_t *timer, int idx, int i, int j, int t)
{
    if (ignites(s, R, idx, i, j, t)) {
        W[idx] = BURNING;
        timer[idx] = (uint8_t)(1 + (int)(s->fuel[idx] * (float)s->burn_time));
    } else {
        W[idx] = UNBURNT;
    }
}

static inline void update_burning(uint8_t *W, uint8_t *timer, int idx)
{
    uint8_t tm = (uint8_t)(timer[idx] - 1);
    timer[idx] = tm;
    W[idx] = tm ? BURNING : BURNED;
}

static inline void update_cell(const Sim *s, const uint8_t *R, uint8_t *W,
                               uint8_t *timer, int i, int j, int t)
{
    int idx = i * s->M + j;
    uint8_t st = R[idx];
    if (st == BURNED)       W[idx] = BURNED;
    else if (st == BURNING) update_burning(W, timer, idx);
    else                    update_unburnt(s, R, W, timer, idx, i, j, t);
}

static inline void swap_buffers(Sim *s)
{
    uint8_t *tmp = s->R; s->R = s->W; s->W = tmp;      /* pointer swap, no copy */
}

static inline void record_busy(Sim *s, double t0)
{
    int tid = omp_get_thread_num();
    if (tid < MAX_THREADS) s->busy[tid * BUSY_STRIDE] += omp_get_wtime() - t0;
}

/* ---------------------------------------------------------------------- */
/* 1) SERIAL reference -- no OpenMP pragmas. Used for correctness + T_serial */
/* ---------------------------------------------------------------------- */
void step_serial(Sim *s, int t)
{
    const int N = s->N, M = s->M;
    const uint8_t *R = s->R;
    uint8_t *W = s->W, *tm = s->timer;

    for (int i = 1; i < N - 1; i++)
        for (int j = 1; j < M - 1; j++)
            update_cell(s, R, W, tm, i, j, t);

    swap_buffers(s);
}

/* ---------------------------------------------------------------------- */
/* 2) BASE kernel -- straightforward OpenMP (Srikripa)                      */
/*    collapse(2) merges both loops; schedule(runtime) lets us pick          */
/*    static / dynamic / guided from the command line or OMP_SCHEDULE.       */
/*    NOTE: with collapse(2) the chunk size counts CELLS, not rows.          */
/* ---------------------------------------------------------------------- */
void step_base(Sim *s, int t)
{
    const int N = s->N, M = s->M;
    const uint8_t *R = s->R;
    uint8_t *W = s->W, *tm = s->timer;

    #pragma omp parallel
    {
        double t0 = omp_get_wtime();

        #pragma omp for collapse(2) schedule(runtime) nowait
        for (int i = 1; i < N - 1; i++)
            for (int j = 1; j < M - 1; j++)
                update_cell(s, R, W, tm, i, j, t);

        record_busy(s, t0);          /* time until THIS thread ran out of work */
    }                                /* implicit barrier here                  */

    swap_buffers(s);
}

/* ---------------------------------------------------------------------- */
/* 3) OPT kernel -- row-wise + vectorised pre-pass (Gautami)                */
/*    Pass 1 (SIMD): flag every UNBURNT cell that has a burning neighbour.   */
/*    Pass 2       : only flagged cells run the expensive exp()/RNG code.    */
/*    restrict tells the compiler the arrays do not overlap.                 */
/*    Chunk size counts ROWS here.                                           */
/* ---------------------------------------------------------------------- */
void step_opt(Sim *s, int t)
{
    const int N = s->N, M = s->M;
    const uint8_t *restrict R = s->R;
    uint8_t *restrict W  = s->W;
    uint8_t *restrict tm = s->timer;
    uint8_t *restrict act = s->act;

    #pragma omp parallel
    {
        double t0 = omp_get_wtime();

        #pragma omp for schedule(runtime) nowait
        for (int i = 1; i < N - 1; i++) {
            const uint8_t *restrict r0 = R + (size_t)(i - 1) * M;
            const uint8_t *restrict r1 = R + (size_t)i * M;
            const uint8_t *restrict r2 = R + (size_t)(i + 1) * M;
            uint8_t *restrict a = act + (size_t)i * M;

            /* pass 1: branch-free neighbour test, auto/explicit SIMD */
            OMP_SIMD
            for (int j = 1; j < M - 1; j++) {
                unsigned nb = (r0[j-1] == BURNING) + (r0[j] == BURNING) + (r0[j+1] == BURNING)
                            + (r1[j-1] == BURNING)                       + (r1[j+1] == BURNING)
                            + (r2[j-1] == BURNING) + (r2[j] == BURNING) + (r2[j+1] == BURNING);
                a[j] = (uint8_t)((r1[j] == UNBURNT) & (nb > 0));
            }

            /* pass 2: scalar work only where something can happen */
            for (int j = 1; j < M - 1; j++) {
                int idx = i * M + j;
                uint8_t st = r1[j];
                if (a[j])                 update_unburnt(s, R, W, tm, idx, i, j, t);
                else if (st == BURNING)   update_burning(W, tm, idx);
                else                      W[idx] = st;     /* UNBURNT or BURNED */
            }
        }

        record_busy(s, t0);
    }

    swap_buffers(s);
}
