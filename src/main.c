/*
 * main.c  --  driver (Gautami): command-line, timing, checksum, verification,
 *             metrics, CSV output.
 *
 * stdout : ONE csv line per run (scripts append it to the results CSV files)
 * stderr : human readable summary (hidden with --quiet)
 */
#include "wildfire.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

typedef struct {
    int n, m, steps, threads, chunk, ignitions, burn_time;
    int verify, verify_every, snap_every, ascii, quiet;
    char terrain[16], kernel[16], sched[16], out[256], tag[32];
    double wind_speed, wind_angle, p0, cs;
    uint64_t seed;
} Args;

static void usage(const char *p)
{
    fprintf(stderr,
"Usage: %s [options]\n"
"  --n N              grid rows (default 1024)       --m M   columns (default = n)\n"
"  --steps S          time steps (default 200)       --seed X  RNG seed (default 42)\n"
"  --kernel K         serial | base | opt   (default opt)\n"
"  --threads T        number of OpenMP threads (default: all)\n"
"  --sched S          static | dynamic | guided | auto   (default: OMP_SCHEDULE / static)\n"
"  --chunk C          schedule chunk size (base: cells, opt: rows; 0 = default)\n"
"  --terrain T        flat | hilly (default hilly)   --ignitions K (default 1)\n"
"  --wind-speed W     (default 1.0)   --wind-angle A  degrees, 0=E 90=S (default 45)\n"
"  --p0 P  --cs C  --burn-time B      model parameters (0.10 / 2.0 / 8)\n"
"  --verify           also run the serial reference and compare every step\n"
"  --verify-every K   compare every K steps (default 1)\n"
"  --out FILE.ppm     write final image;   --snap-every K  also write every K steps\n"
"  --ascii            print a small ASCII picture at the end\n"
"  --tag LABEL        free text label stored in the CSV (e.g. 'novec')\n"
"  --quiet            no summary on stderr      --csv-header  print CSV header and exit\n",
        p);
}

#define CSV_HEADER "kernel,terrain,n,m,steps,threads,schedule,chunk,seed,ignitions," \
                   "time_s,mcups,imbalance,burned_frac,checksum,verify,tag"

/* ---- FNV-1a 64-bit checksum ------------------------------------------- */
static uint64_t fnv1a(const uint8_t *p, size_t n, uint64_t h)
{
    for (size_t i = 0; i < n; i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}
static uint64_t checksum_state(const Sim *s)
{
    return fnv1a(s->R, (size_t)s->N * s->M, 1469598103934665603ULL);
}
static uint64_t checksum_full(const Sim *s)
{
    size_t n = (size_t)s->N * s->M;
    return fnv1a(s->timer, n, fnv1a(s->R, n, 1469598103934665603ULL));
}

static Sim *make_sim(const Args *a)
{
    Sim *s = sim_create(a->n, a->m, a->seed);
    s->p0 = a->p0; s->cs = a->cs; s->burn_time = a->burn_time;
    sim_set_wind(s, a->wind_speed, a->wind_angle);
    init_terrain(s, a->terrain);
    init_ignition(s, a->ignitions);
    return s;
}

static void snapshot_name(char *dst, size_t len, const char *out, int t)
{
    char base[256];
    snprintf(base, sizeof base, "%s", out);
    size_t L = strlen(base);
    if (L > 4 && strcmp(base + L - 4, ".ppm") == 0) base[L - 4] = '\0';
    snprintf(dst, len, "%s_t%05d.ppm", base, t);
}

int main(int argc, char **argv)
{
    Args a = { .n = 1024, .m = 0, .steps = 200, .threads = 0, .chunk = 0,
               .ignitions = 1, .burn_time = 8, .verify = 0, .verify_every = 1,
               .snap_every = 0, .ascii = 0, .quiet = 0,
               .terrain = "hilly", .kernel = "opt", .sched = "", .out = "", .tag = "default",
               .wind_speed = 1.0, .wind_angle = 45.0, .p0 = 0.10, .cs = 2.0, .seed = 42 };

    for (int k = 1; k < argc; k++) {
        const char *o = argv[k];
        #define NEXT (k + 1 < argc ? argv[++k] : (usage(argv[0]), exit(1), (char *)0))
        if      (!strcmp(o, "--n"))           a.n = atoi(NEXT);
        else if (!strcmp(o, "--m"))           a.m = atoi(NEXT);
        else if (!strcmp(o, "--steps"))       a.steps = atoi(NEXT);
        else if (!strcmp(o, "--threads"))     a.threads = atoi(NEXT);
        else if (!strcmp(o, "--chunk"))       a.chunk = atoi(NEXT);
        else if (!strcmp(o, "--ignitions"))   a.ignitions = atoi(NEXT);
        else if (!strcmp(o, "--burn-time"))   a.burn_time = atoi(NEXT);
        else if (!strcmp(o, "--verify-every"))a.verify_every = atoi(NEXT);
        else if (!strcmp(o, "--snap-every"))  a.snap_every = atoi(NEXT);
        else if (!strcmp(o, "--seed"))        a.seed = strtoull(NEXT, NULL, 10);
        else if (!strcmp(o, "--wind-speed"))  a.wind_speed = atof(NEXT);
        else if (!strcmp(o, "--wind-angle"))  a.wind_angle = atof(NEXT);
        else if (!strcmp(o, "--p0"))          a.p0 = atof(NEXT);
        else if (!strcmp(o, "--cs"))          a.cs = atof(NEXT);
        else if (!strcmp(o, "--kernel"))      snprintf(a.kernel, sizeof a.kernel, "%s", NEXT);
        else if (!strcmp(o, "--terrain"))     snprintf(a.terrain, sizeof a.terrain, "%s", NEXT);
        else if (!strcmp(o, "--sched"))       snprintf(a.sched, sizeof a.sched, "%s", NEXT);
        else if (!strcmp(o, "--out"))         snprintf(a.out, sizeof a.out, "%s", NEXT);
        else if (!strcmp(o, "--tag"))         snprintf(a.tag, sizeof a.tag, "%s", NEXT);
        else if (!strcmp(o, "--verify"))      a.verify = 1;
        else if (!strcmp(o, "--ascii"))       a.ascii = 1;
        else if (!strcmp(o, "--quiet"))       a.quiet = 1;
        else if (!strcmp(o, "--csv-header"))  { puts(CSV_HEADER); return 0; }
        else if (!strcmp(o, "--help") || !strcmp(o, "-h")) { usage(argv[0]); return 0; }
        else { fprintf(stderr, "unknown option %s\n", o); usage(argv[0]); return 1; }
    }
    if (a.m <= 0) a.m = a.n;
    if (a.n < 5 || a.m < 5 || a.steps < 1) { fprintf(stderr, "grid too small\n"); return 1; }
    if (a.verify_every < 1) a.verify_every = 1;
    if (a.burn_time < 1 || a.burn_time > 200) { fprintf(stderr, "burn-time must be 1..200\n"); return 1; }
    if (strcmp(a.terrain, "flat") && strcmp(a.terrain, "hilly")) { fprintf(stderr, "terrain: flat|hilly\n"); return 1; }

    /* ---- pick kernel ---- */
    void (*step)(Sim *, int);
    if      (!strcmp(a.kernel, "serial")) { step = step_serial; a.threads = 1; }
    else if (!strcmp(a.kernel, "base"))     step = step_base;
    else if (!strcmp(a.kernel, "opt"))      step = step_opt;
    else { fprintf(stderr, "kernel: serial|base|opt\n"); return 1; }

    /* ---- OpenMP setup ---- */
    if (a.threads > 0) omp_set_num_threads(a.threads);
    int nth = omp_get_max_threads();
    if (nth > MAX_THREADS) { nth = MAX_THREADS; omp_set_num_threads(nth); }
    if (a.sched[0]) {
        omp_sched_t kind = omp_sched_static;
        if      (!strcmp(a.sched, "static"))  kind = omp_sched_static;
        else if (!strcmp(a.sched, "dynamic")) kind = omp_sched_dynamic;
        else if (!strcmp(a.sched, "guided"))  kind = omp_sched_guided;
        else if (!strcmp(a.sched, "auto"))    kind = omp_sched_auto;
        else { fprintf(stderr, "sched: static|dynamic|guided|auto\n"); return 1; }
        omp_set_schedule(kind, a.chunk);
    }
    else if (!getenv("OMP_SCHEDULE"))
        omp_set_schedule(omp_sched_static, 0);   /* libgomp's own default is dynamic,1 !  */
    omp_sched_t sk; int sc;
    omp_get_schedule(&sk, &sc);
    const char *sname = sk == omp_sched_static ? "static" : sk == omp_sched_dynamic ? "dynamic"
                      : sk == omp_sched_guided ? "guided" : "auto";

    /* ---- build the world ---- */
    Sim *s = make_sim(&a);
    Sim *ref = a.verify ? make_sim(&a) : NULL;

    char snap[300];
    int verify_ok = 1, fail_step = -1;
    double elapsed = 0.0;

    if (a.snap_every > 0 && a.out[0]) { snapshot_name(snap, sizeof snap, a.out, 0); write_ppm(s, snap); }

    for (int t = 0; t < a.steps; t++) {
        double t0 = omp_get_wtime();
        step(s, t);
        elapsed += omp_get_wtime() - t0;            /* only the kernel is timed */

        if (ref) {
            step_serial(ref, t);
            if (((t + 1) % a.verify_every == 0 || t == a.steps - 1) &&
                checksum_state(s) != checksum_state(ref)) { verify_ok = 0; fail_step = t; break; }
        }
        if (a.snap_every > 0 && a.out[0] && (t + 1) % a.snap_every == 0) {
            snapshot_name(snap, sizeof snap, a.out, t + 1);
            write_ppm(s, snap);
        }
    }
    if (ref && verify_ok && checksum_full(s) != checksum_full(ref)) { verify_ok = 0; fail_step = a.steps - 1; }

    /* ---- metrics ---- */
    size_t cells = (size_t)s->N * s->M, burned = 0, burning = 0;
    for (size_t k = 0; k < cells; k++) { burned += (s->R[k] == BURNED); burning += (s->R[k] == BURNING); }
    double burned_frac = (double)(burned + burning) / (double)cells;

    double updates = (double)(a.n - 2) * (double)(a.m - 2) * (double)a.steps;
    double mcups = updates / elapsed / 1e6;

    double bmax = 0, bsum = 0;
    for (int k = 0; k < nth; k++) {
        double b = s->busy[k * BUSY_STRIDE];
        if (b > bmax) bmax = b;
        bsum += b;
    }
    double imbalance = (bsum > 0) ? bmax / (bsum / nth) : 1.0;   /* 1.0 = perfectly balanced */

    uint64_t cs = checksum_full(s);
    const char *vtxt = !a.verify ? "NA" : (verify_ok ? "PASS" : "FAIL");

    printf("%s,%s,%d,%d,%d,%d,%s,%d,%llu,%d,%.6f,%.3f,%.3f,%.5f,%016llx,%s,%s\n",
           a.kernel, a.terrain, a.n, a.m, a.steps, nth, sname, sc,
           (unsigned long long)a.seed, a.ignitions, elapsed, mcups, imbalance, burned_frac,
           (unsigned long long)cs, vtxt, a.tag);

    if (!a.quiet) {
        fprintf(stderr, "\n=== Wildfire CA (%s kernel) ===\n", a.kernel);
        fprintf(stderr, "grid %d x %d, %d steps, terrain=%s, threads=%d, schedule=%s,%d\n",
                a.n, a.m, a.steps, a.terrain, nth, sname, sc);
        fprintf(stderr, "time        : %.4f s\n", elapsed);
        fprintf(stderr, "throughput  : %.1f MCUPS (million cell updates / s)\n", mcups);
        fprintf(stderr, "imbalance   : %.3f (max thread busy time / mean; 1.0 = balanced)\n", imbalance);
        fprintf(stderr, "burnt area  : %.2f %% of grid\n", 100.0 * burned_frac);
        fprintf(stderr, "checksum    : %016llx\n", (unsigned long long)cs);
        if (a.verify) {
            if (verify_ok) fprintf(stderr, "verification: PASS (identical to serial at every checked step)\n");
            else           fprintf(stderr, "verification: FAIL at step %d\n", fail_step);
        }
    }
    if (a.out[0]) write_ppm(s, a.out);
    if (a.ascii) print_ascii(s, 78);

    sim_free(s); sim_free(ref);
    return (a.verify && !verify_ok) ? 2 : 0;
}
