/* Q6190 scratch-only first-Witt/companion census.
 *
 * For primes p <= LIMIT this computes C_n=(n!)^3 b_n modulo p^2 and the
 * normalized companion D_n=(n!)^3 c_n modulo p, both with
 *
 *   X_{n+1}=P(n)X_n-n^6 X_{n-1},  P(n)=34n^3+51n^2+27n+5.
 *
 * At an actual zero p|b_r, eta_r=(b_r/p) mod p is recovered from C_r/p.
 * The original companion c_r is recovered from D_r/(r!)^3.  For every
 * consecutive zero pair x<y we record the exact transfer unit
 *
 *   lambda = -c_y/c_x mod p
 *
 * and the divided-transfer/Bockstein residue
 *
 *   B = eta_y + lambda eta_x mod p.
 *
 * For every reflected pair y=p-1-x we also record
 *   g_x = eta_x-eta_y mod p,
 * the first-Witt reflection jet.  Legendre classes are evidence only.
 */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ZEROS 256

typedef unsigned __int128 u128;

typedef struct {
    uint32_t r;
    uint32_t eta;
    uint32_t comp;
} Zero;

static uint64_t modpow64(uint64_t a, uint64_t e, uint64_t m) {
    uint64_t out = 1 % m;
    while (e) {
        if (e & 1) out = (uint64_t)((u128)out * a % m);
        a = (uint64_t)((u128)a * a % m);
        e >>= 1;
    }
    return out;
}

static int color(uint32_t a, uint32_t p) {
    if (a == 0) return 0;
    uint32_t q = (uint32_t)modpow64(a, (p - 1) / 2, p);
    if (q == 1) return 1;
    if (q == p - 1) return 2;
    fprintf(stderr, "bad Legendre value p=%u a=%u q=%u\n", p, a, q);
    exit(2);
}

static uint32_t subp(uint32_t a, uint32_t b, uint32_t p) {
    return a >= b ? a - b : a + p - b;
}

static uint32_t addp(uint32_t a, uint32_t b, uint32_t p) {
    uint64_t s = (uint64_t)a + b;
    return (uint32_t)(s >= p ? s - p : s);
}

static uint64_t polyP(uint64_t n, uint64_t mod) {
    return (uint64_t)((((u128)34 * n + 51) * n + 27) * n + 5) % mod;
}

static uint64_t pow6(uint64_t n, uint64_t mod) {
    uint64_t n2 = (uint64_t)((u128)n * n % mod);
    uint64_t n3 = (uint64_t)((u128)n2 * n % mod);
    return (uint64_t)((u128)n3 * n3 % mod);
}

static uint64_t linrec_next(uint64_t Pn, uint64_t cur, uint64_t n6,
                            uint64_t prev, uint64_t mod) {
    uint64_t a = (uint64_t)((u128)Pn * cur % mod);
    uint64_t b = (uint64_t)((u128)n6 * prev % mod);
    return a >= b ? a - b : a + mod - b;
}

static uint8_t *sieve(uint32_t limit) {
    uint8_t *prime = calloc((size_t)limit + 1, 1);
    if (!prime) exit(2);
    memset(prime, 1, (size_t)limit + 1);
    prime[0] = prime[1] = 0;
    for (uint32_t d = 2; (uint64_t)d * d <= limit; ++d)
        if (prime[d])
            for (uint32_t m = d * d; m <= limit; m += d)
                prime[m] = 0;
    return prime;
}

int main(int argc, char **argv) {
    uint32_t limit = 30000;
    const char *csv_path = "q6190_witt_rows.csv";
    if (argc > 1) limit = (uint32_t)strtoul(argv[1], NULL, 10);
    if (argc > 2) csv_path = argv[2];
    if (limit < 5 || limit > 1000000) {
        fprintf(stderr, "LIMIT must be in [5,1000000]\n");
        return 2;
    }
    FILE *csv = fopen(csv_path, "w");
    if (!csv) { perror("fopen"); return 2; }
    fputs("p,r,eta,comp\n", csv);

    uint8_t *isprime = sieve(limit);
    uint64_t primes = 0, total_zeros = 0, reflection_pairs = 0;
    uint64_t consecutive_pairs = 0, midpoint_zeros = 0;
    uint64_t gcolor[3] = {0,0,0};
    uint64_t lcolor[3] = {0,0,0};
    uint64_t bcolor[3] = {0,0,0};
    uint32_t maxz = 0, maxz_p = 0;
    uint32_t midpoint_list[64]; size_t nmid = 0;

    for (uint32_t p = 5; p <= limit; ++p) {
        if (!isprime[p]) continue;
        ++primes;
        uint64_t p2 = (uint64_t)p * p;
        Zero zs[MAX_ZEROS];
        uint32_t nz = 0;

        uint64_t bprev = 1 % p2, bcur = 5 % p2;
        uint32_t cprev = 0, ccur = 1 % p;
        uint32_t fact = 1;

        if (bcur % p == 0) {
            if (nz >= MAX_ZEROS) return 3;
            uint32_t eta = (uint32_t)((bcur / p) % p);
            zs[nz++] = (Zero){1, eta, ccur};
        }

        for (uint32_t n = 1; n < p - 1; ++n) {
            uint32_t r = n + 1;
            uint64_t P2 = polyP(n, p2);
            uint64_t N62 = pow6(n, p2);
            uint64_t bnext = linrec_next(P2, bcur, N62, bprev, p2);

            uint32_t P1 = (uint32_t)(P2 % p);
            uint32_t N61 = (uint32_t)(N62 % p);
            uint32_t cnext = (uint32_t)linrec_next(P1, ccur, N61, cprev, p);

            fact = (uint32_t)((uint64_t)fact * r % p);
            if (bnext % p == 0) {
                if (nz >= MAX_ZEROS) {
                    fprintf(stderr, "MAX_ZEROS exceeded p=%u\n", p);
                    return 3;
                }
                uint32_t fact2 = (uint32_t)((uint64_t)fact * fact % p);
                uint32_t fact3 = (uint32_t)((uint64_t)fact2 * fact % p);
                uint32_t invfact3 = (uint32_t)modpow64(fact3, p - 2, p);
                uint32_t eta_norm = (uint32_t)((bnext / p) % p);
                uint32_t eta = (uint32_t)((uint64_t)eta_norm * invfact3 % p);
                uint32_t comp = (uint32_t)((uint64_t)cnext * invfact3 % p);
                zs[nz++] = (Zero){r, eta, comp};
            }
            bprev = bcur; bcur = bnext;
            cprev = ccur; ccur = cnext;
        }

        total_zeros += nz;
        if (nz > maxz) { maxz = nz; maxz_p = p; }
        if ((nz & 1U) != 0) {
            uint32_t mid = (p - 1) / 2;
            int found = 0;
            for (uint32_t i = 0; i < nz; ++i) if (zs[i].r == mid) found = 1;
            if (!found) { fprintf(stderr, "parity/midpoint failure p=%u\n", p); return 4; }
        }
        for (uint32_t i = 0; i < nz; ++i) {
            fprintf(csv, "%u,%u,%u,%u\n", p, zs[i].r, zs[i].eta, zs[i].comp);
            if (zs[i].r == (p - 1) / 2) {
                ++midpoint_zeros;
                if (nmid < 64) midpoint_list[nmid++] = p;
            }
            uint32_t j = nz - 1 - i;
            if (zs[i].r + zs[j].r != p - 1) {
                fprintf(stderr, "reflection failure p=%u i=%u r=%u mirror=%u\n",
                        p, i, zs[i].r, zs[j].r);
                return 5;
            }
        }
        for (uint32_t i = 0; i < nz / 2; ++i) {
            uint32_t j = nz - 1 - i;
            uint32_t g = subp(zs[i].eta, zs[j].eta, p);
            ++gcolor[color(g,p)];
            ++reflection_pairs;
        }
        for (uint32_t i = 0; i + 1 < nz; ++i) {
            Zero x = zs[i], y = zs[i+1];
            if (x.comp == 0 || y.comp == 0) {
                fprintf(stderr, "companion vanishes at zero p=%u x=%u y=%u\n", p,x.r,y.r);
                return 6;
            }
            uint32_t invcx = (uint32_t)modpow64(x.comp, p - 2, p);
            uint32_t lambda = (uint32_t)((uint64_t)y.comp * invcx % p);
            lambda = lambda ? p - lambda : 0; /* -c_y/c_x */
            uint32_t bock = addp(y.eta, (uint32_t)((uint64_t)lambda * x.eta % p), p);
            ++lcolor[color(lambda,p)];
            ++bcolor[color(bock,p)];
            ++consecutive_pairs;
        }
    }
    fclose(csv);
    free(isprime);

    printf("Q6190_WITT_CENSUS limit=%u primes=%" PRIu64 " total_zeros=%" PRIu64
           " maxZ=%u maxZ_p=%u\n", limit, primes, total_zeros, maxz, maxz_p);
    printf("midpoint_zeros=%" PRIu64 " midpoint_primes=", midpoint_zeros);
    for (size_t i=0;i<nmid;++i) printf("%s%u", i?",":"", midpoint_list[i]);
    printf("\n");
    printf("reflection_pairs=%" PRIu64 " g_color_zero=%" PRIu64
           " g_color_square=%" PRIu64 " g_color_nonsquare=%" PRIu64 "\n",
           reflection_pairs,gcolor[0],gcolor[1],gcolor[2]);
    printf("consecutive_zero_pairs=%" PRIu64 " lambda_zero=%" PRIu64
           " lambda_square=%" PRIu64 " lambda_nonsquare=%" PRIu64 "\n",
           consecutive_pairs,lcolor[0],lcolor[1],lcolor[2]);
    printf("bockstein_zero=%" PRIu64 " bockstein_square=%" PRIu64
           " bockstein_nonsquare=%" PRIu64 "\n",
           bcolor[0],bcolor[1],bcolor[2]);
    printf("rows_csv=%s\n", csv_path);
    return 0;
}
