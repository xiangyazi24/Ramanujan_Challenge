/*
 * Fast computation of Z(p) = #{j < p : b_j ≡ 0 mod p} for all primes p up to PMAX.
 * Uses the Apéry recurrence b_{n+1} = [P(n) b_n - n³ b_{n-1}] / (n+1)³ mod p.
 * P(n) = 34n³ + 51n² + 27n + 5.
 *
 * Compile: gcc -O3 -march=native -o zp_fast zp_fast.c -lm
 * Usage:   ./zp_fast [PMAX]    (default 50000000)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

typedef unsigned long long u64;
typedef long long i64;

static int *sieve_primes(int limit, int *count) {
    char *is_composite = calloc(limit + 1, 1);
    if (!is_composite) { fprintf(stderr, "OOM sieve\n"); exit(1); }
    for (int i = 2; (i64)i * i <= limit; i++)
        if (!is_composite[i])
            for (int j = i * i; j <= limit; j += i)
                is_composite[j] = 1;
    int n = 0;
    for (int i = 2; i <= limit; i++)
        if (!is_composite[i]) n++;
    int *primes = malloc(n * sizeof(int));
    if (!primes) { fprintf(stderr, "OOM primes\n"); exit(1); }
    int k = 0;
    for (int i = 2; i <= limit; i++)
        if (!is_composite[i]) primes[k++] = i;
    free(is_composite);
    *count = n;
    return primes;
}

static inline i64 mod(i64 a, i64 p) {
    a %= p;
    return a < 0 ? a + p : a;
}

static inline i64 modinv(i64 a, i64 p) {
    i64 g = p, x = 0, y = 1;
    i64 aa = a;
    while (aa != 0) {
        i64 q = g / aa;
        i64 t = g - q * aa; g = aa; aa = t;
        t = x - q * y; x = y; y = t;
    }
    return mod(x, p);
}

static inline i64 mulmod(i64 a, i64 b, i64 p) {
    return ((__int128)a * b) % p;
}

static int compute_Zp(int p) {
    if (p < 5) return 0;
    i64 P = p;
    int zcount = 0;

    i64 bprev = 1; // b_0
    i64 bcur = 5;  // b_1
    if (bprev == 0) zcount++;
    if (bcur == 0) zcount++;

    for (int n = 1; n < p - 1; n++) {
        // b_{n+1} = [P(n) * b_n - n³ * b_{n-1}] / (n+1)³
        i64 nn = n;
        i64 Pn = mod(((34 * nn + 51) * nn + 27) * nn + 5, P);
        i64 n3 = mulmod(mulmod(nn, nn, P), nn, P);
        i64 np1 = nn + 1;
        i64 np1_3 = mulmod(mulmod(np1, np1, P), np1, P);
        i64 np1_3_inv = modinv(np1_3, P);

        i64 term1 = mulmod(Pn, bcur, P);
        i64 term2 = mulmod(n3, bprev, P);
        i64 bnext = mulmod(mod(term1 - term2, P), np1_3_inv, P);

        bprev = bcur;
        bcur = bnext;
        if (bnext == 0) zcount++;
    }
    return zcount;
}

int main(int argc, char **argv) {
    int pmax = 50000000;
    if (argc > 1) pmax = atoi(argv[1]);

    printf("Computing Z(p) for primes 5 <= p <= %d...\n", pmax);
    fflush(stdout);

    int nprimes;
    int *primes = sieve_primes(pmax, &nprimes);
    printf("  %d primes total, sieve done.\n", nprimes);
    fflush(stdout);

    // Skip p=2,3
    int start = 0;
    while (start < nprimes && primes[start] < 5) start++;

    int maxZ = 0;
    int maxZ_p = 0;
    double sumZ = 0;
    int hist[64] = {0};
    int count = 0;
    int record_progression[32];
    int record_count = 0;
    int non_ordinary_count = 0;

    time_t t0 = time(NULL);

    for (int i = start; i < nprimes; i++) {
        int p = primes[i];
        int z = compute_Zp(p);
        count++;
        sumZ += z;

        if (z < 64) hist[z]++;
        if (z > maxZ) {
            maxZ = z;
            maxZ_p = p;
            if (record_count < 32) {
                record_progression[record_count++] = p;
            }
            printf("  NEW MAX Z(p) = %d at p = %d\n", z, p);
            fflush(stdout);
        }
        if (z % 2 == 1) {
            non_ordinary_count++;
            if (non_ordinary_count <= 20)
                printf("  NON-ORDINARY: p=%d, Z(p)=%d\n", p, z);
        }

        if (count % 100000 == 0 || i == nprimes - 1) {
            time_t now = time(NULL);
            double elapsed = difftime(now, t0);
            double frac = (double)(i - start + 1) / (nprimes - start);
            double eta = frac > 0 ? elapsed / frac * (1 - frac) : 0;
            printf("  [%5.1f%%] p=%d, avg Z=%.4f, max Z=%d (p=%d), "
                   "elapsed %.0fs, ETA %.0fs\n",
                   frac * 100, p, sumZ / count, maxZ, maxZ_p,
                   elapsed, eta);
            fflush(stdout);
        }
    }

    printf("\n=== FINAL RESULTS ===\n");
    printf("  Primes tested: %d (from 5 to %d)\n", count, primes[nprimes-1]);
    printf("  Max Z(p) = %d at p = %d\n", maxZ, maxZ_p);
    printf("  Mean Z(p) = %.6f\n", sumZ / count);
    printf("  Non-ordinary primes (odd Z): %d\n", non_ordinary_count);

    printf("\n  Histogram:\n");
    for (int z = 0; z < 64; z++) {
        if (hist[z] > 0)
            printf("    Z(p) = %2d: %d (%.4f%%)\n",
                   z, hist[z], 100.0 * hist[z] / count);
    }

    printf("\n  Record progression (primes where max Z increased):\n");
    for (int j = 0; j < record_count; j++)
        printf("    p = %d\n", record_progression[j]);

    double elapsed = difftime(time(NULL), t0);
    printf("\n  Total time: %.0f seconds (%.1f minutes)\n", elapsed, elapsed/60);

    free(primes);
    return 0;
}
