/*
 * wildfire.h  --  shared definitions for the Parallel Wildfire CA project
 * (write this together on Day 1: it is the contract between sim.c, init.c, main.c)
 *
 * Grid layout : 1D flat row-major array, idx = i*M + j   (i = row, j = column)
 * Cell states : UNBURNT -> BURNING -> BURNED
 * Borders     : the outermost ring of cells is never updated (stays UNBURNT).
 */
#ifndef WILDFIRE_H
#define WILDFIRE_H

#include <stdint.h>
#include <stddef.h>

/* ---- states ---------------------------------------------------------- */
enum { UNBURNT = 0, BURNING = 1, BURNED = 2 };

/* ---- SIMD pragma can be switched off at compile time (-DNO_SIMD) ----- */
#ifdef NO_SIMD
#define OMP_SIMD
#else
#define OMP_SIMD _Pragma("omp simd")
#endif

#define MAX_THREADS 512
#define BUSY_STRIDE 8            /* doubles per slot -> one cache line each */

/* ---- simulation state ------------------------------------------------- */
typedef struct {
    int N, M;                    /* rows, columns                            */
    uint8_t *R, *W;              /* state double buffer (read / write)       */
    uint8_t *timer;              /* remaining burn steps per cell            */
    uint8_t *act;                /* scratch: "may ignite" flags (opt kernel) */
    float   *fuel;               /* K_fuel in [0,1]                          */
    float   *elev;               /* elevation (DEM), 1 unit = 1 cell width   */

    /* model parameters */
    double   p0;                 /* base ignition probability                */
    double   cs;                 /* slope sensitivity                        */
    int      burn_time;          /* max burn duration (steps)                */
    uint64_t seed;               /* RNG seed                                 */

    /* tables precomputed from grid width and wind (see sim_set_wind) */
    int    off[8];               /* neighbour index offsets                  */
    double inv_dist[8];          /* 1 / distance to neighbour (1 or 1/sqrt2) */
    double wind_fac[8];          /* exp(V_wind . d_n) for each direction     */

    /* per-thread busy time, used for the load-imbalance metric */
    double *busy;
} Sim;

/* ---- RNG: counter-based, so results do NOT depend on thread count ------ */
static inline uint64_t splitmix64(uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

/* uniform double in [0,1) that depends only on (seed, i, j, t) */
static inline double hash_rand(uint64_t seed, int i, int j, int t)
{
    uint64_t k = splitmix64(seed ^ ((uint64_t)(uint32_t)i * 0x9E3779B97F4A7C15ULL));
    k = splitmix64(k ^ ((uint64_t)(uint32_t)j * 0xC2B2AE3D27D4EB4FULL));
    k = splitmix64(k ^ ((uint64_t)(uint32_t)t * 0x165667B19E3779F9ULL));
    return (double)(k >> 11) * (1.0 / 9007199254740992.0);   /* 2^-53 */
}

/* ======================= init.c  (Kumud) ================================ */
Sim *sim_create(int N, int M, uint64_t seed);          /* allocate + defaults */
void sim_free(Sim *s);
void sim_set_wind(Sim *s, double speed, double angle_deg);
void init_terrain(Sim *s, const char *kind);           /* "flat" | "hilly"    */
void init_ignition(Sim *s, int k);                     /* k ignition points   */
int  write_ppm(const Sim *s, const char *path);        /* colour image        */
void print_ascii(const Sim *s, int max_cols);          /* quick terminal view */

/* ======================= sim.c   (Srikripa / Gautami) =================== */
void step_serial(Sim *s, int t);   /* reference: no OpenMP at all             */
void step_base(Sim *s, int t);     /* OpenMP: collapse(2) + schedule(runtime) */
void step_opt(Sim *s, int t);      /* OpenMP: row-wise, SIMD pre-pass         */

#endif
