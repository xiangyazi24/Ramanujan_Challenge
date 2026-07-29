/*
 * p32_zp_fast.c — Compute Z(p) for all primes up to PMAX.
 * Z(p) = #{j in [0,p) : b_j ≡ 0 (mod p)} where b_n are Apery numbers.
 *
 * Compile: gcc -O3 -o p32_zp_fast p32_zp_fast.c -lm
 * Usage:   ./p32_zp_fast [PMAX]    (default 1000000)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef unsigned long long u64;
typedef long long i64;

/* Modular exponentiation */
static inline u64 powmod(u64 base, u64 exp, u64 mod) {
    u64 result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = (__uint128_t)result * base % mod;
        exp >>= 1;
        base = (__uint128_t)base * base % mod;
    }
    return result;
}

/* Modular inverse via Fermat */
static inline u64 modinv(u64 a, u64 p) {
    return powmod(a, p - 2, p);
}

/* Sieve of Eratosthenes */
static int *sieve_primes(int N, int *count) {
    char *is_prime = calloc(N + 1, 1);
    for (int i = 2; i <= N; i++) is_prime[i] = 1;
    for (int i = 2; (i64)i * i <= N; i++)
        if (is_prime[i])
            for (int j = i * i; j <= N; j += i)
                is_prime[j] = 0;
    *count = 0;
    for (int i = 2; i <= N; i++) if (is_prime[i]) (*count)++;
    int *primes = malloc(*count * sizeof(int));
    int k = 0;
    for (int i = 2; i <= N; i++) if (is_prime[i]) primes[k++] = i;
    free(is_prime);
    return primes;
}

/* Count Z(p) using recurrence mod p */
static int count_zeros(int p) {
    if (p < 5) return 0;
    u64 P = (u64)p;
    int z = 0;

    u64 b_prev = 1;  /* b[0] */
    u64 b_curr = 5 % P;  /* b[1] */
    if (b_prev == 0) z++;
    if (b_curr == 0) z++;

    for (int n = 1; n < p - 1; n++) {
        u64 N = (u64)n;
        /* coeff = 34n^3 + 51n^2 + 27n + 5 */
        u64 n2 = (__uint128_t)N * N % P;
        u64 n3 = (__uint128_t)n2 * N % P;
        u64 coeff = ((__uint128_t)34 * n3 + (__uint128_t)51 * n2 + 27 * N + 5) % P;
        u64 prev_n3 = n3;
        u64 np1 = N + 1;
        u64 np1_2 = (__uint128_t)np1 * np1 % P;
        u64 np1_3 = (__uint128_t)np1_2 * np1 % P;
        if (np1_3 == 0) break;  /* shouldn't happen for n < p-1 */
        u64 inv = modinv(np1_3, P);

        u64 num = ((__uint128_t)coeff * b_curr + P - (__uint128_t)prev_n3 * b_prev % P) % P;
        u64 b_next = (__uint128_t)num * inv % P;
        b_prev = b_curr;
        b_curr = b_next;
        if (b_curr == 0) z++;
    }
    return z;
}

int main(int argc, char **argv) {
    int PMAX = 1000000;
    if (argc > 1) PMAX = atoi(argv[1]);

    printf("Computing Z(p) for primes 5 <= p <= %d...\n", PMAX);
    fflush(stdout);

    int nprimes;
    int *primes = sieve_primes(PMAX, &nprimes);
    printf("  %d primes total\n", nprimes);
    fflush(stdout);

    /* Filter to p >= 5 */
    int start = 0;
    while (start < nprimes && primes[start] < 5) start++;
    int count5 = nprimes - start;
    printf("  %d primes >= 5\n", count5);
    fflush(stdout);

    int *zp = malloc(count5 * sizeof(int));
    int hist[256] = {0};
    long long total_z = 0;
    int max_z = 0, max_z_p = 0;
    int report_interval = count5 / 100;
    if (report_interval < 1) report_interval = 1;

    time_t t0 = time(NULL);
    clock_t c0 = clock();

    for (int i = 0; i < count5; i++) {
        int p = primes[start + i];
        int z = count_zeros(p);
        zp[i] = z;
        if (z < 256) hist[z]++;
        total_z += z;
        if (z > max_z) {
            max_z = z;
            max_z_p = p;
            printf("  NEW MAX Z(p) = %d at p = %d\n", z, p);
            fflush(stdout);
        }
        if ((i + 1) % report_interval == 0) {
            double pct = 100.0 * (i + 1) / count5;
            double elapsed = (double)(clock() - c0) / CLOCKS_PER_SEC;
            double avg = (double)total_z / (i + 1);
            double rate = (i + 1) / elapsed;
            double eta = (count5 - i - 1) / rate;
            printf("  [%5.1f%%] p=%d, avg Z=%.3f, max Z=%d (p=%d), "
                   "elapsed %.0fs, ETA %.0fs\n",
                   pct, p, avg, max_z, max_z_p, elapsed, eta);
            fflush(stdout);
        }
    }

    double elapsed = (double)(clock() - c0) / CLOCKS_PER_SEC;
    printf("\nDone in %.1fs (%d primes)\n", elapsed, count5);

    /* Histogram */
    printf("\n=== Z(p) histogram (p <= %d) ===\n", PMAX);
    for (int z = 0; z < 256; z++) {
        if (hist[z] > 0) {
            double pct = 100.0 * hist[z] / count5;
            printf("  Z(p) = %3d: %6d primes (%6.2f%%)\n", z, hist[z], pct);
        }
    }

    /* Summary */
    printf("\n=== Summary ===\n");
    printf("  Total primes (>= 5): %d\n", count5);
    printf("  Mean Z(p): %.4f\n", (double)total_z / count5);
    printf("  Max Z(p): %d at p = %d\n", max_z, max_z_p);
    printf("  Primes with Z(p)=0: %d (%.1f%%)\n",
           hist[0], 100.0 * hist[0] / count5);

    /* Average by range */
    printf("\n=== Average Z(p) by range ===\n");
    int boundaries[] = {5, 100, 1000, 10000, 100000, 500000, 1000000, 5000000, 10000000};
    int nb = sizeof(boundaries) / sizeof(boundaries[0]);
    for (int b = 0; b < nb - 1; b++) {
        int lo = boundaries[b], hi = boundaries[b + 1];
        if (lo >= PMAX) break;
        if (hi > PMAX + 1) hi = PMAX + 1;
        long long s = 0;
        int cnt = 0, mx = 0;
        for (int i = 0; i < count5; i++) {
            int p = primes[start + i];
            if (p >= lo && p < hi) {
                s += zp[i];
                cnt++;
                if (zp[i] > mx) mx = zp[i];
            }
        }
        if (cnt > 0)
            printf("  p in [%7d, %7d): avg=%.3f, max=%d, n=%d\n",
                   lo, hi, (double)s / cnt, mx, cnt);
    }

    /* Primes with Z(p) >= 6 */
    printf("\n=== Primes with Z(p) >= 6 ===\n");
    int notable_count = 0;
    for (int i = 0; i < count5; i++) {
        if (zp[i] >= 6) {
            int p = primes[start + i];
            printf("  p=%d: Z(p)=%d\n", p, zp[i]);
            notable_count++;
            if (notable_count >= 100) {
                printf("  ... (showing first 100)\n");
                break;
            }
        }
    }

    /* Non-ordinary primes (odd Z(p)) */
    printf("\n=== Non-ordinary primes (odd Z(p)) ===\n");
    for (int i = 0; i < count5; i++) {
        if (zp[i] % 2 == 1) {
            int p = primes[start + i];
            printf("  p=%d: Z(p)=%d\n", p, zp[i]);
        }
    }

    free(primes);
    free(zp);
    return 0;
}
