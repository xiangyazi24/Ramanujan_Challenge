/*
 * The full conjecture quantity.  By the Apery-Lucas property b_n = prod_i b_{d_i} (mod p),
 * where d_i are the base-p digits of n, we have
 *     p | b_n   <=>   some base-p digit of n lies in Z_p = {r < p : p | b_r}.
 * This computes, for every n <= NMAX,
 *     R(n) = sum over primes p <= n with p | b_n of log p      (= log rad_{p<=n}(b_n)),
 *     W(n) = the number of such primes,
 * which is what the conjecture requires to be o(n) (up to the bounded multiplicities
 * e_p <= 3 of the denominator defect).
 *
 * gcc -O3 -o rad_scan rad_scan.c -lm ; ./rad_scan [NMAX]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef unsigned long long u64;

static char *sieve(int limit) {
    char *comp = calloc(limit + 1, 1);
    for (int i = 2; (long long)i * i <= limit; i++)
        if (!comp[i])
            for (int j = i * i; j <= limit; j += i) comp[j] = 1;
    return comp;
}

int main(int argc, char **argv) {
    int NMAX = (argc > 1) ? atoi(argv[1]) : 100000;
    char *comp = sieve(NMAX);
    double *R = calloc((size_t)NMAX + 1, sizeof(double));
    unsigned short *W = calloc((size_t)NMAX + 1, sizeof(unsigned short));
    u64 *inv = malloc(sizeof(u64) * ((size_t)NMAX + 2));
    char *inZ = malloc((size_t)NMAX + 2);
    clock_t t0 = clock();

    for (int p = 5; p <= NMAX; p++) {
        if (comp[p]) continue;
        inv[1] = 1;
        for (int k = 2; k < p; k++)
            inv[k] = (u64)(p - (p / k)) * inv[p % k] % (u64)p;
        memset(inZ, 0, p);
        u64 b0 = 1 % p, b1 = 5 % p;
        if (b0 == 0) inZ[0] = 1;
        if (b1 == 0) inZ[1] = 1;
        for (int j = 1; j <= p - 2; j++) {
            u64 jj = (u64)j % p;
            u64 P = (34ULL * jj % p * jj % p * jj + 51ULL * jj % p * jj
                     + 27ULL * jj + 5ULL) % p;
            u64 j3 = jj * jj % p * jj % p;
            u64 ik = inv[j + 1];
            u64 lead = ik * ik % p * ik % p;
            u64 nxt = (P * b1 % p + (u64)p * p - j3 * b0 % p) % p * lead % p;
            if (nxt == 0) inZ[j + 1] = 1;
            b0 = b1; b1 = nxt;
        }
        double lp = log((double)p);
        /* mark every n <= NMAX having a base-p digit in Z_p */
        for (int n = p; n <= NMAX; n++) {
            int m = n, hit = 0;
            while (m) { if (inZ[m % p]) { hit = 1; break; } m /= p; }
            if (hit) { R[n] += lp; W[n]++; }
        }
        if (p < 200 || p % 20011 == 0)
            fprintf(stderr, "  p=%d %.0fs\n", p, (double)(clock() - t0) / CLOCKS_PER_SEC);
    }

    double best = 0; int bestn = 0, maxW = 0, maxWn = 0;
    for (int n = 2; n <= NMAX; n++) {
        if (R[n] / n > best) { best = R[n] / n; bestn = n; }
        if (W[n] > maxW) { maxW = W[n]; maxWn = n; }
    }
    printf("NMAX=%d\n", NMAX);
    printf("max R(n)/n = %.6f at n=%d   (R = %.2f)\n", best, bestn, R[bestn]);
    printf("max #{p<=n : p|b_n} = %d at n=%d\n", maxW, maxWn);
    for (int t = 1000; t <= NMAX; t *= 10) {
        double b2 = 0; int b2n = 0;
        for (int n = t / 10; n <= t && n <= NMAX; n++)
            if (R[n] / n > b2) { b2 = R[n] / n; b2n = n; }
        printf("  n in [%d,%d]: max R(n)/n = %.6f at n=%d\n", t / 10, t, b2, b2n);
    }
    return 0;
}
