/*
 * K(n) = #{ p prime, n/2 < p <= n : p | b_n }, via Apery-Lucas  p | b_n <=> p | b_{n-p}.
 *
 * For each prime p <= NMAX we iterate the Apery recurrence mod p, record every r < p with
 * b_r == 0 (mod p), and bump K(p+r) (note r < p forces p > n/2 automatically).
 *
 * Compile: gcc -O3 -march=native -o k_scan k_scan.c
 * Usage:   ./k_scan [NMAX]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

typedef unsigned long long u64;

static char *sieve(int limit) {
    char *comp = calloc(limit + 1, 1);
    for (int i = 2; (long long)i * i <= limit; i++)
        if (!comp[i])
            for (int j = i * i; j <= limit; j += i) comp[j] = 1;
    return comp;
}

int main(int argc, char **argv) {
    int NMAX = (argc > 1) ? atoi(argv[1]) : 1000000;
    char *comp = sieve(NMAX);
    unsigned char *K = calloc((size_t)NMAX + 1, 1);
    u64 *inv = malloc(sizeof(u64) * ((size_t)NMAX + 2));
    long long totalzeros = 0;
    int maxZp = 0, maxZp_p = 0;
    clock_t t0 = clock();

    for (int p = 5; p <= NMAX; p++) {
        if (comp[p]) continue;
        /* inverses 1..p-1 mod p */
        inv[1] = 1;
        for (int k = 2; k < p; k++)
            inv[k] = (u64)(p - (p / k)) * inv[p % k] % (u64)p;

        u64 b0 = 1 % p, b1 = 5 % p;
        int zp = 0;
        if (b0 == 0) { if ((long long)p + 0 <= NMAX) K[p]++; zp++; totalzeros++; }
        if (b1 == 0) { if ((long long)p + 1 <= NMAX) K[p + 1]++; zp++; totalzeros++; }
        for (int j = 1; j <= p - 2; j++) {
            u64 jj = (u64)j % p;
            u64 P = (34ULL * jj % p * jj % p * jj
                     + 51ULL * jj % p * jj
                     + 27ULL * jj + 5ULL) % p;
            u64 j3 = jj * jj % p * jj % p;
            u64 ik = inv[j + 1];
            u64 lead = ik * ik % p * ik % p;
            u64 nxt = (P * b1 % p + (u64)p * p - j3 * b0 % p) % p * lead % p;
            if (nxt == 0) {
                long long n = (long long)p + (j + 1);
                if (n <= NMAX) K[n]++;
                zp++;
                totalzeros++;
            }
            b0 = b1; b1 = nxt;
        }
        if (zp > maxZp) { maxZp = zp; maxZp_p = p; }
        if (p % 100000 < 2)
            fprintf(stderr, "  p=%d  %.0fs\n", p, (double)(clock() - t0) / CLOCKS_PER_SEC);
    }

    long long cnt[8] = {0};
    int best = 0, bestn = 0;
    for (int n = 2; n <= NMAX; n++) {
        int k = K[n];
        if (k > best) { best = k; bestn = n; printf("new max K=%d at n=%d\n", k, n); }
        if (k < 8) cnt[k]++;
    }
    printf("NMAX=%d  max K(n)=%d at n=%d\n", NMAX, best, bestn);
    for (int k = 0; k < 6; k++) printf("  #{n : K(n)=%d} = %lld\n", k, cnt[k]);
    printf("total zeros (sum |Z_p|) = %lld over primes <= %d ; max |Z_p| = %d at p=%d\n",
           totalzeros, NMAX, maxZp, maxZp_p);
    printf("mean |Z_p| ~ %.4f\n", (double)totalzeros / (NMAX / (double)(log((double)NMAX) - 1)));
    return 0;
}
