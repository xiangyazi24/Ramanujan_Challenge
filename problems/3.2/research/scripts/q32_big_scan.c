/*
 * Large-scale scan of the Apery zero sets and the conjecture quantity.
 *
 *   Z_p = { r < p : p | b_r },   b_r = Apery numbers for zeta(3).
 *
 * Division-free recurrence: with B_r = b_r (r!)^3,
 *     B_{r+1} = P(r) B_r - r^6 B_{r-1},   P(r) = 34r^3+51r^2+27r+5,
 * and for r < p, p | b_r <=> p | B_r (since (r!)^3 is a unit).  No modular inverses.
 *
 * Phase 1 (threaded): for every prime p <= PMAX collect the pairs (p, r), r in Z_p.
 * Phase 2: for every n <= NMAX accumulate
 *     R(n) = sum of log p over p <= n with (n mod p) in Z_p        [the crux sum]
 *     K(n) = the number of such p with p > n/2                     [top-window targets]
 *
 * Build: gcc -O3 -march=native -pthread -o big_scan big_scan.c -lm
 * Run:   ./big_scan PMAX NMAX [threads]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

typedef unsigned long long u64;
typedef unsigned int u32;

static char *sieve(u32 limit) {
    char *comp = calloc((size_t)limit + 1, 1);
    for (u64 i = 2; i * i <= limit; i++)
        if (!comp[i])
            for (u64 j = i * i; j <= limit; j += i) comp[j] = 1;
    return comp;
}

static u32 *primes;
static u32 nprimes;
static char *comp;

/* collected pairs */
typedef struct { u32 p, r; } pair_t;
static pair_t *pairs;
static size_t npairs, cap_pairs;
static pthread_mutex_t plock = PTHREAD_MUTEX_INITIALIZER;

static u32 PMAX, NMAX, NTHREADS;
static volatile u32 next_index;
static pthread_mutex_t ilock = PTHREAD_MUTEX_INITIALIZER;

static u64 zeros_total, maxZ, maxZ_p;

static void *worker(void *arg) {
    (void)arg;
    pair_t *local = malloc(sizeof(pair_t) * 4096);
    size_t nlocal = 0;
    u64 local_zeros = 0, local_maxZ = 0, local_maxZ_p = 0;
    for (;;) {
        u32 idx;
        pthread_mutex_lock(&ilock);
        idx = next_index++;
        pthread_mutex_unlock(&ilock);
        if (idx >= nprimes) break;
        u32 p = primes[idx];
        if (p < 5) continue;
        u64 P = p;
        /* Barrett: mu = floor(2^64 / P); mulmod via 128-bit high product */
        const u64 mu = (u64)(((((__uint128_t)1) << 64) - 1) / P);
        #define MULMOD(a, b) ({ u64 _x = (a) * (b); \
            u64 _q = (u64)(((__uint128_t)_x * mu) >> 64); \
            u64 _r = _x - _q * P; if (_r >= P) _r -= P; if (_r >= P) _r -= P; _r; })
        #define ADDM(a, b) ({ u64 _s = (a) + (b); if (_s >= P) _s -= P; _s; })
        u64 b0 = 1 % P, b1 = 5 % P;      /* B_0 = b_0, B_1 = b_1 (1!)^3 = 5 */
        u64 zc = 0;
        if (b0 == 0) { local[nlocal++] = (pair_t){p, 0}; zc++; }
        if (b1 == 0) { local[nlocal++] = (pair_t){p, 1}; zc++; }
        /* incremental cubics: track j, j^2, j^3 and P(j) by finite differences */
        u64 jm = 1 % P, j2 = 1 % P, j3 = 1 % P;
        for (u64 j = 1; j + 1 < P; j++) {
            u64 j6 = MULMOD(j3, j3);
            /* 34 j3 + 51 j2 + 27 j + 5 < 112 P < 2^39: reduce once */
            u64 praw = 34 * j3 + 51 * j2 + 27 * jm + 5;
            u64 pq = (u64)(((__uint128_t)praw * mu) >> 64);
            u64 Pv = praw - pq * P; if (Pv >= P) Pv -= P; if (Pv >= P) Pv -= P;
            u64 t1 = MULMOD(Pv, b1), t2 = MULMOD(j6, b0);
            u64 nxt = (t1 >= t2) ? (t1 - t2) : (t1 + P - t2);
            if (nxt == 0) {
                if (nlocal == 4096) {
                    pthread_mutex_lock(&plock);
                    if (npairs + nlocal > cap_pairs) {
                        cap_pairs = (npairs + nlocal) * 2;
                        pairs = realloc(pairs, cap_pairs * sizeof(pair_t));
                    }
                    memcpy(pairs + npairs, local, nlocal * sizeof(pair_t));
                    npairs += nlocal;
                    pthread_mutex_unlock(&plock);
                    nlocal = 0;
                }
                local[nlocal++] = (pair_t){p, (u32)(j + 1)};
                zc++;
            }
            b0 = b1; b1 = nxt;
            /* j -> j+1 */
            u64 jn = ADDM(jm, 1 % P);
            j3 = ADDM(ADDM(j3, MULMOD(3 % P, j2)), ADDM(MULMOD(3 % P, jm), 1 % P));
            j2 = ADDM(ADDM(j2, ADDM(jm, jm)), 1 % P);
            jm = jn;
        }
        #undef MULMOD
        #undef ADDM
        local_zeros += zc;
        if (zc > local_maxZ) { local_maxZ = zc; local_maxZ_p = p; }
    }
    pthread_mutex_lock(&plock);
    if (npairs + nlocal > cap_pairs) {
        cap_pairs = (npairs + nlocal) * 2 + 1024;
        pairs = realloc(pairs, cap_pairs * sizeof(pair_t));
    }
    memcpy(pairs + npairs, local, nlocal * sizeof(pair_t));
    npairs += nlocal;
    zeros_total += local_zeros;
    if (local_maxZ > maxZ) { maxZ = local_maxZ; maxZ_p = local_maxZ_p; }
    pthread_mutex_unlock(&plock);
    free(local);
    return NULL;
}

int main(int argc, char **argv) {
    PMAX = (argc > 1) ? (u32)atof(argv[1]) : 1000000;
    NMAX = (argc > 2) ? (u32)atof(argv[2]) : PMAX;
    NTHREADS = (argc > 3) ? atoi(argv[3]) : 8;
    clock_t t0 = clock();
    time_t w0 = time(NULL);

    comp = sieve(PMAX);
    nprimes = 0;
    for (u32 i = 2; i <= PMAX; i++) if (!comp[i]) nprimes++;
    primes = malloc(sizeof(u32) * nprimes);
    { u32 k = 0; for (u32 i = 2; i <= PMAX; i++) if (!comp[i]) primes[k++] = i; }
    fprintf(stderr, "sieved: %u primes <= %u  [%lds]\n", nprimes, PMAX, time(NULL) - w0);

    cap_pairs = 1 << 20;
    pairs = malloc(cap_pairs * sizeof(pair_t));
    pthread_t *th = malloc(sizeof(pthread_t) * NTHREADS);
    for (u32 i = 0; i < NTHREADS; i++) pthread_create(&th[i], NULL, worker, NULL);
    for (u32 i = 0; i < NTHREADS; i++) pthread_join(th[i], NULL);
    fprintf(stderr, "phase 1 done: %zu zero pairs, %lds wall\n", npairs, time(NULL) - w0);

    printf("PMAX=%u NMAX=%u\n", PMAX, NMAX);
    printf("sum |Z_p| = %llu over %u primes -> mean %.4f\n",
           zeros_total, nprimes, (double)zeros_total / nprimes);
    printf("max |Z_p| = %llu at p = %llu\n", maxZ, maxZ_p);

    /* optional dump of (p,r) pairs */
    { const char *dp = getenv("DUMP_PAIRS");
      if (dp) { FILE *f = fopen(dp, "w");
        for (size_t i = 0; i < npairs; i++) fprintf(f, "%u %u\n", pairs[i].p, pairs[i].r);
        fclose(f); fprintf(stderr, "dumped %zu pairs to %s\n", npairs, dp); } }

    /* phase 2: accumulate R(n), K(n) */
    float *R = calloc((size_t)NMAX + 1, sizeof(float));
    unsigned char *K = calloc((size_t)NMAX + 1, 1);
    for (size_t i = 0; i < npairs; i++) {
        u32 p = pairs[i].p, r = pairs[i].r;
        float lp = (float)log((double)p);
        for (u64 n = (u64)p + r; n <= NMAX; n += p) {
            R[n] += lp;
            if ((u64)2 * p > n && K[n] < 255) K[n]++;
        }
    }
    fprintf(stderr, "phase 2 done [%lds wall]\n", time(NULL) - w0);

    double best = 0; u64 bestn = 0; int maxK = 0; u64 maxKn = 0;
    for (u64 n = 2; n <= NMAX; n++) {
        double v = R[n] / (double)n;
        if (v > best) { best = v; bestn = n; }
        if (K[n] > maxK) { maxK = K[n]; maxKn = n; }
    }
    printf("max R(n)/n = %.7f at n = %llu   (R = %.2f)\n", best, bestn, (double)R[bestn]);
    printf("max K(n) [top-window targets] = %d at n = %llu\n", maxK, maxKn);
    for (u64 t = 1000; t <= NMAX; t *= 10) {
        double b2 = 0; u64 b2n = 0; int k2 = 0;
        for (u64 n = t / 10; n <= t && n <= NMAX; n++) {
            double v = R[n] / (double)n;
            if (v > b2) { b2 = v; b2n = n; }
            if (K[n] > k2) k2 = K[n];
        }
        printf("  n in [%llu,%llu]: max R/n = %.7f at n=%llu (R=%.1f), max K = %d\n",
               t / 10, t, b2, b2n, (double)R[b2n], k2);
    }
    printf("cpu %.0fs wall %lds\n", (double)(clock() - t0) / CLOCKS_PER_SEC,
           time(NULL) - w0);
    return 0;
}
