/*
 * Fast Z(p) computation using renormalized Apéry recurrence (NO division).
 * B_{n+1} = P(n) B_n - n^6 B_{n-1}  where P(n) = 34n³+51n²+27n+5
 * B_n = (n!)³ b_n, so b_n ≡ 0 mod p iff B_n ≡ 0 mod p (for n < p).
 *
 * Compile: gcc -O3 -march=native -o zp_fast2 zp_fast2.c -lm
 * Usage:   ./zp_fast2 [PMAX]    (default 50000000)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

typedef unsigned long long u64;
typedef long long i64;

static char *sieve;
static int sieve_limit;

static void init_sieve(int limit) {
    sieve_limit = limit;
    sieve = calloc(limit + 1, 1);
    if (!sieve) { fprintf(stderr, "OOM sieve\n"); exit(1); }
    for (i64 i = 2; i * i <= limit; i++)
        if (!sieve[i])
            for (i64 j = i * i; j <= limit; j += i)
                sieve[j] = 1;
}

static inline i64 mulmod(i64 a, i64 b, i64 p) {
    return ((__int128)a * b) % p;
}

static int compute_Zp(i64 p) {
    if (p < 5) return 0;
    int zcount = 0;

    // B_0 = 0!^3 * b_0 = 1
    // B_1 = 1!^3 * b_1 = 5
    i64 Bprev = 1 % p;
    i64 Bcur = 5 % p;
    if (Bprev == 0) zcount++;  // b_0 = 0?
    if (Bcur == 0) zcount++;   // b_1 = 0?

    for (i64 n = 1; n < p - 1; n++) {
        // B_{n+1} = P(n) B_n - n^6 B_{n-1}  mod p
        // P(n) = 34n³ + 51n² + 27n + 5
        i64 n2 = mulmod(n, n, p);
        i64 n3 = mulmod(n2, n, p);
        i64 Pn = (mulmod(34, n3, p) + mulmod(51, n2, p) + mulmod(27, n, p) + 5) % p;
        i64 n6 = mulmod(n3, n3, p);

        i64 term1 = mulmod(Pn, Bcur, p);
        i64 term2 = mulmod(n6, Bprev, p);
        i64 Bnext = (term1 - term2 + p) % p;

        Bprev = Bcur;
        Bcur = Bnext;
        if (Bnext == 0) zcount++;
    }
    return zcount;
}

int main(int argc, char **argv) {
    int pmax = 50000000;
    if (argc > 1) pmax = atoi(argv[1]);

    printf("Computing Z(p) for primes 5 <= p <= %d (no-division version)...\n", pmax);
    fflush(stdout);

    time_t t0 = time(NULL);
    init_sieve(pmax);
    printf("  Sieve done in %.0fs.\n", difftime(time(NULL), t0));
    fflush(stdout);

    int maxZ = 0;
    int maxZ_p = 0;
    double sumZ = 0;
    int count = 0;
    int hist[64] = {0};
    int record_progression_p[64];
    int record_progression_z[64];
    int record_count = 0;
    int non_ordinary_count = 0;

    for (int p = 5; p <= pmax; p++) {
        if (sieve[p]) continue;
        int z = compute_Zp(p);
        count++;
        sumZ += z;

        if (z < 64) hist[z]++;
        if (z > maxZ) {
            maxZ = z;
            maxZ_p = p;
            if (record_count < 64) {
                record_progression_p[record_count] = p;
                record_progression_z[record_count] = z;
                record_count++;
            }
            printf("  NEW MAX Z(p) = %d at p = %d\n", z, p);
            fflush(stdout);
        }
        if (z % 2 == 1) {
            non_ordinary_count++;
            if (non_ordinary_count <= 20)
                printf("  NON-ORDINARY: p=%d, Z(p)=%d\n", p, z);
        }

        if (count % 100000 == 0) {
            time_t now = time(NULL);
            double elapsed = difftime(now, t0);
            printf("  [%d primes] p=%d, avg Z=%.4f, max Z=%d (p=%d), "
                   "elapsed %.0fs\n",
                   count, p, sumZ / count, maxZ, maxZ_p, elapsed);
            fflush(stdout);
        }
    }

    time_t tf = time(NULL);
    double total = difftime(tf, t0);

    printf("\n=== FINAL RESULTS ===\n");
    printf("  Primes tested: %d (from 5 to %d)\n", count, pmax);
    printf("  Max Z(p) = %d at p = %d\n", maxZ, maxZ_p);
    printf("  Mean Z(p) = %.6f\n", sumZ / count);
    printf("  Non-ordinary primes (odd Z): %d\n", non_ordinary_count);

    printf("\n  Histogram:\n");
    for (int z = 0; z < 64; z++) {
        if (hist[z] > 0)
            printf("    Z(p) = %2d: %8d (%.4f%%)\n",
                   z, hist[z], 100.0 * hist[z] / count);
    }

    printf("\n  Poisson(1/2) comparison for pair count K = Z/2:\n");
    double lambda = 0.5;
    double poi[32];
    poi[0] = exp(-lambda);
    for (int k = 1; k < 32; k++) poi[k] = poi[k-1] * lambda / k;
    for (int k = 0; k < 8; k++) {
        int z = 2*k;
        double obs = (double)hist[z] / count;
        printf("    K=%d (Z=%2d): observed=%.5f, Poisson(1/2)=%.5f, ratio=%.3f\n",
               k, z, obs, poi[k], obs / (poi[k] > 0 ? poi[k] : 1));
    }

    printf("\n  Record progression:\n");
    for (int j = 0; j < record_count; j++)
        printf("    Z(p) = %d at p = %d\n", record_progression_z[j], record_progression_p[j]);

    printf("\n  Total time: %.0f seconds (%.1f minutes)\n", total, total/60);

    free(sieve);
    return 0;
}
