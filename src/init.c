/*
 * init.c  --  setup and visualisation (Kumud)
 *   - memory allocation
 *   - wind lookup tables
 *   - synthetic terrain (elevation) and fuel map
 *   - ignition points
 *   - PPM image output and ASCII preview
 *
 * No real DEM data is needed: terrain is generated from sine waves and a
 * Gaussian hill, fuel from smooth patches plus hash noise.  Everything is
 * deterministic given the seed.
 */
#include "wildfire.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979323846

/* aligned allocation, zero-filled; size rounded up to a multiple of 64 */
static void *xalloc(size_t bytes)
{
    size_t n = (bytes + 63) & ~(size_t)63;
    void *p = aligned_alloc(64, n);
    if (!p) { fprintf(stderr, "out of memory (%zu bytes)\n", n); exit(1); }
    memset(p, 0, n);
    return p;
}

Sim *sim_create(int N, int M, uint64_t seed)
{
    Sim *s = calloc(1, sizeof(Sim));
    if (!s) exit(1);
    size_t cells = (size_t)N * (size_t)M;

    s->N = N; s->M = M; s->seed = seed;
    s->R     = xalloc(cells);
    s->W     = xalloc(cells);
    s->timer = xalloc(cells);
    s->act   = xalloc(cells);
    s->fuel  = xalloc(cells * sizeof(float));
    s->elev  = xalloc(cells * sizeof(float));
    s->busy  = xalloc((size_t)MAX_THREADS * BUSY_STRIDE * sizeof(double));

    /* default model parameters (can be overridden from the command line) */
    s->p0 = 0.10;
    s->cs = 2.0;
    s->burn_time = 8;

    /* neighbour offsets: (di,dj) pairs, row-major index */
    static const int di[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
    static const int dj[8] = {-1,  0,  1, -1, 1, -1, 0, 1};
    for (int d = 0; d < 8; d++) {
        s->off[d] = di[d] * M + dj[d];
        s->inv_dist[d] = (di[d] != 0 && dj[d] != 0) ? 1.0 / sqrt(2.0) : 1.0;
    }
    sim_set_wind(s, 0.0, 0.0);
    return s;
}

void sim_free(Sim *s)
{
    if (!s) return;
    free(s->R); free(s->W); free(s->timer); free(s->act);
    free(s->fuel); free(s->elev); free(s->busy); free(s);
}

/*
 * Wind blows TOWARD angle_deg:  0 = east (+j), 90 = south (+i / down in the
 * image), 180 = west, 270 = north.  For a neighbour n at offset (di,dj) the
 * fire travels n -> cell, i.e. in direction (-dj, -di) (x = column, y = row).
 * wind_fac = exp( speed * cos(angle between wind and that direction) ).
 */
void sim_set_wind(Sim *s, double speed, double angle_deg)
{
    static const int di[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
    static const int dj[8] = {-1,  0,  1, -1, 1, -1, 0, 1};
    double a = angle_deg * PI / 180.0;
    double wx = speed * cos(a), wy = speed * sin(a);
    for (int d = 0; d < 8; d++) {
        double ux = -dj[d], uy = -di[d];
        double len = sqrt(ux * ux + uy * uy);
        s->wind_fac[d] = exp((wx * ux + wy * uy) / len);
    }
}

static inline double clamp01(double x) { return x < 0 ? 0 : (x > 1 ? 1 : x); }

void init_terrain(Sim *s, const char *kind)
{
    const int N = s->N, M = s->M;

    if (strcmp(kind, "flat") == 0) {            /* uniform: even workload */
        for (size_t k = 0; k < (size_t)N * M; k++) { s->fuel[k] = 0.5f; s->elev[k] = 0.0f; }
        return;
    }

    /* "hilly": rolling terrain, patchy fuel -> irregular, uneven workload */
    const double lam = (N > M ? N : M) / 3.0;   /* wavelength of the hills   */
    const double amp = 0.04 * lam;              /* => max slope about 0.25   */
    const double hi = 0.7 * N, hj = 0.3 * M, sig = 0.08 * N;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            size_t idx = (size_t)i * M + j;

            double e = amp * sin(2 * PI * i / lam + 0.7) * cos(2 * PI * j / lam + 0.3);
            double dx = i - hi, dy = j - hj;
            e += 0.5 * amp * exp(-(dx * dx + dy * dy) / (2 * sig * sig));
            s->elev[idx] = (float)e;

            double patch = 0.5 + 0.5 * sin(2 * PI * j / (0.7 * lam) + 1.3)
                                     * sin(2 * PI * i / (0.9 * lam));
            double noise = hash_rand(s->seed, i, j, 0x7fffffff);
            s->fuel[idx] = (float)clamp01(0.15 + 0.6 * patch + 0.25 * noise);
        }
    }
}

static void ignite_cell(Sim *s, int i, int j)
{
    size_t idx = (size_t)i * s->M + j;
    s->R[idx] = s->W[idx] = BURNING;
    s->timer[idx] = (uint8_t)(1 + (int)(s->fuel[idx] * (float)s->burn_time));
}

/* k ignition points: the first is a 3x3 block in the centre, the rest are
 * placed pseudo-randomly (deterministic from the seed). */
void init_ignition(Sim *s, int k)
{
    int ci = s->N / 2, cj = s->M / 2;
    for (int di = -1; di <= 1; di++)
        for (int dj = -1; dj <= 1; dj++)
            ignite_cell(s, ci + di, cj + dj);

    for (int p = 1; p < k; p++) {
        int i = 3 + (int)(hash_rand(s->seed, p, 1, -7) * (s->N - 6));
        int j = 3 + (int)(hash_rand(s->seed, p, 2, -7) * (s->M - 6));
        ignite_cell(s, i, j);
    }
}

/* PPM (P6) image: green = unburnt (brighter = more fuel), orange = burning,
 * dark grey = burned.  View with viz.py (Pillow) or any image viewer. */
int write_ppm(const Sim *s, const char *path)
{
    FILE *f = fopen(path, "wb");
    if (!f) { perror(path); return -1; }
    fprintf(f, "P6\n%d %d\n255\n", s->M, s->N);

    unsigned char *row = malloc((size_t)s->M * 3);
    if (!row) { fclose(f); return -1; }
    for (int i = 0; i < s->N; i++) {
        for (int j = 0; j < s->M; j++) {
            size_t idx = (size_t)i * s->M + j;
            unsigned char r, g, b;
            switch (s->R[idx]) {
                case BURNING: r = 255; g = 110; b = 0;  break;
                case BURNED:  r = 45;  g = 45;  b = 45; break;
                default:      r = 25;  g = (unsigned char)(70 + 120 * s->fuel[idx]); b = 25; break;
            }
            row[3 * j] = r; row[3 * j + 1] = g; row[3 * j + 2] = b;
        }
        fwrite(row, 1, (size_t)s->M * 3, f);
    }
    free(row);
    fclose(f);
    return 0;
}

/* '.' unburnt   '*' burning   '#' burned   (downsampled to <= max_cols) */
void print_ascii(const Sim *s, int max_cols)
{
    int f = (s->M + max_cols - 1) / max_cols;
    if (f < 1) f = 1;
    for (int i = 0; i < s->N; i += 2 * f) {          /* rows are ~2x taller */
        for (int j = 0; j < s->M; j += f) {
            int burning = 0, burned = 0;
            for (int a = i; a < i + 2 * f && a < s->N; a++)
                for (int b = j; b < j + f && b < s->M; b++) {
                    uint8_t st = s->R[(size_t)a * s->M + b];
                    burning |= (st == BURNING);
                    burned  |= (st == BURNED);
                }
            putchar(burning ? '*' : (burned ? '#' : '.'));
        }
        putchar('\n');
    }
}
