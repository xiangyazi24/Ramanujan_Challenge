/*
 * bn_bigscan.c -- large-scale top-half Apéry bad-prime scan.
 *
 * For every prime p >= 7, this program finds
 *
 *     Z_p = { r : 0 <= r < p and b_r == 0 (mod p) }
 *
 * and scatters every (p,r) to n=p+r.  Thus the scatter count at n is
 * exactly #{p : n/2 < p <= n and p divides b_{n-p}}.
 *
 * The division-free state A_m=(m!)^3 b_m obeys
 *
 *     A_{m+1} = (34m^3+51m^2+27m+5) A_m - m^6 A_{m-1} (mod p).
 *
 * P(m) and m^6 are advanced by finite differences.  The two remaining
 * modular products per recurrence step use 32-bit Montgomery arithmetic.
 * Primes are dynamically distributed, largest first, among pthreads.
 *
 * Compile:
 *     cc -O2 -o bn_bigscan bn_bigscan.c -lpthread
 *
 * Usage:
 *     ./bn_bigscan [N_MAX [P_MAX [THREADS [PAIR_FILE]]]]
 *
 * Defaults: N_MAX=P_MAX=2000000, THREADS=online CPUs,
 *           PAIR_FILE=data_zp_pairs.bin.
 *
 * PAIR_FILE is a headerless stream of 8-byte little-endian records:
 * uint32 p followed by uint32 r.  Records are sorted by (p,r), and include
 * every zero for p <= P_MAX, even when p+r > N_MAX.
 */

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define DEFAULT_MAX 2000000U
#define FIRST_PRIME 7U
#define DYADIC_START 1024U
#define WINDOW_START (1U << 17)
#define LN2_CONST 0.693147180559945309417232121458176568

typedef struct {
    uint32_t p;
    uint32_t r;
} Pair;

typedef struct {
    Pair *data;
    size_t len;
    size_t cap;
} PairVec;

typedef struct {
    uint32_t *primes;
    size_t nprimes;
    uint32_t *zcounts;
    atomic_size_t next_job;
    atomic_size_t primes_done;
    atomic_uint_fast64_t steps_done;
    atomic_int workers_done;
    atomic_int failed;
} Shared;

typedef struct {
    Shared *shared;
    PairVec pairs;
} Worker;

static double monotonic_seconds(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
        perror("clock_gettime");
        exit(EXIT_FAILURE);
    }
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int pair_push(PairVec *v, uint32_t p, uint32_t r) {
    if (v->len == v->cap) {
        size_t new_cap = v->cap ? 2 * v->cap : 1024;
        if (new_cap < v->cap || new_cap > SIZE_MAX / sizeof(*v->data))
            return 0;
        Pair *new_data = (Pair *)realloc(v->data,
                                         new_cap * sizeof(*v->data));
        if (!new_data)
            return 0;
        v->data = new_data;
        v->cap = new_cap;
    }
    v->data[v->len++] = (Pair){p, r};
    return 1;
}

static uint32_t *sieve_primes(uint32_t limit, size_t *count_out,
                              double *seconds_out) {
    double t0 = monotonic_seconds();
    uint8_t *composite = (uint8_t *)calloc((size_t)limit + 1, 1);
    if (!composite) {
        fprintf(stderr, "out of memory in sieve\n");
        return NULL;
    }
    for (uint32_t i = 2; (uint64_t)i * i <= limit; ++i) {
        if (!composite[i]) {
            for (uint32_t j = i * i; j <= limit; j += i)
                composite[j] = 1;
        }
    }
    size_t count = 0;
    for (uint32_t p = FIRST_PRIME; p <= limit; ++p)
        if (!composite[p])
            ++count;
    uint32_t *primes = (uint32_t *)malloc(count * sizeof(*primes));
    if (!primes) {
        fprintf(stderr, "out of memory for prime list\n");
        free(composite);
        return NULL;
    }
    size_t k = 0;
    for (uint32_t p = FIRST_PRIME; p <= limit; ++p)
        if (!composite[p])
            primes[k++] = p;
    free(composite);
    *count_out = count;
    *seconds_out = monotonic_seconds() - t0;
    return primes;
}

/* -p^{-1} modulo 2^32, for odd p. */
static inline uint32_t montgomery_nprime(uint32_t p) {
    uint32_t x = 1;
    x *= 2U - p * x;
    x *= 2U - p * x;
    x *= 2U - p * x;
    x *= 2U - p * x;
    x *= 2U - p * x;
    return 0U - x;
}

/* Inputs and output are Montgomery residues with R=2^32. */
static inline uint32_t montgomery_mul(uint32_t a, uint32_t b,
                                      uint32_t p, uint32_t nprime) {
    uint64_t t = (uint64_t)a * b;
    uint32_t m = (uint32_t)t * nprime;
    uint64_t u = (t + (uint64_t)m * p) >> 32;
    uint32_t r = (uint32_t)u;
    return r >= p ? r - p : r;
}

static inline uint32_t to_montgomery(uint32_t x, uint32_t p) {
    return (uint32_t)(((uint64_t)(x % p) << 32) % p);
}

static inline uint32_t add_mod(uint32_t a, uint32_t b, uint32_t p) {
    uint32_t s = a + b;
    return s >= p ? s - p : s;
}

static inline uint32_t sub_mod(uint32_t a, uint32_t b, uint32_t p) {
    return a >= b ? a - b : a + p - b;
}

/* Advance a forward-difference row in place, low order first. */
static inline void advance_differences(uint32_t *d, int degree, uint32_t p) {
    for (int k = 0; k < degree; ++k)
        d[k] = add_mod(d[k], d[k + 1], p);
}

/*
 * Append all zero positions for one prime.  At m=1, the forward-difference
 * rows are
 *   P(m): [117, 418, 510, 204]
 *   m^6 : [1, 63, 602, 2100, 3360, 2520, 720].
 */
static int scan_prime(uint32_t p, PairVec *pairs, uint32_t *zcount_out) {
    static const uint32_t p_init[4] = {117, 418, 510, 204};
    static const uint32_t sixth_init[7] =
        {1, 63, 602, 2100, 3360, 2520, 720};
    uint32_t pd[4], sixth[7];
    for (int k = 0; k < 4; ++k)
        pd[k] = to_montgomery(p_init[k], p);
    for (int k = 0; k < 7; ++k)
        sixth[k] = to_montgomery(sixth_init[k], p);

    uint32_t nprime = montgomery_nprime(p);
    uint32_t prev = to_montgomery(1, p);
    uint32_t cur = to_montgomery(5, p);
    uint32_t zcount = 0;

    for (uint32_t m = 1; m < p - 1; ++m) {
        uint32_t term1 = montgomery_mul(pd[0], cur, p, nprime);
        uint32_t term2 = montgomery_mul(sixth[0], prev, p, nprime);
        uint32_t next = sub_mod(term1, term2, p);
        uint32_t r = m + 1;
        if (next == 0) {
            if (!pair_push(pairs, p, r))
                return 0;
            ++zcount;
        }
        prev = cur;
        cur = next;
        advance_differences(pd, 3, p);
        advance_differences(sixth, 6, p);
    }
    *zcount_out = zcount;
    return 1;
}

/* Independent, deliberately simple reference kernel for startup tests. */
static int scan_prime_reference(uint32_t p, PairVec *pairs) {
    uint64_t prev = 1;
    uint64_t cur = 5 % p;
    for (uint32_t m = 1; m < p - 1; ++m) {
        uint64_t n2 = (uint64_t)m * m % p;
        uint64_t n3 = n2 * m % p;
        uint64_t coeff = (34 * n3 + 51 * n2 + 27ULL * m + 5) % p;
        uint64_t n6 = n3 * n3 % p;
        uint64_t next = (coeff * cur + p - n6 * prev % p) % p;
        if (next == 0 && !pair_push(pairs, p, m + 1))
            return 0;
        prev = cur;
        cur = next;
    }
    return 1;
}

static int pair_cmp(const void *va, const void *vb) {
    const Pair *a = (const Pair *)va;
    const Pair *b = (const Pair *)vb;
    if (a->p != b->p)
        return a->p < b->p ? -1 : 1;
    if (a->r != b->r)
        return a->r < b->r ? -1 : 1;
    return 0;
}

static int startup_self_test(void) {
    static const uint32_t tests[] = {7, 11, 17, 31, 181, 379, 3137};
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i) {
        PairVec fast = {0}, reference = {0};
        uint32_t zcount = 0;
        int ok = scan_prime(tests[i], &fast, &zcount) &&
                 scan_prime_reference(tests[i], &reference) &&
                 fast.len == reference.len && zcount == fast.len &&
                 (fast.len == 0 ||
                  memcmp(fast.data, reference.data,
                         fast.len * sizeof(*fast.data)) == 0);
        free(fast.data);
        free(reference.data);
        if (!ok) {
            fprintf(stderr, "startup self-test failed at p=%" PRIu32 "\n",
                    tests[i]);
            return 0;
        }
    }
    printf("Startup kernel self-test: PASS (7 selected primes through 3137)\n");
    return 1;
}

static void *worker_main(void *arg) {
    Worker *worker = (Worker *)arg;
    Shared *s = worker->shared;
    for (;;) {
        if (atomic_load_explicit(&s->failed, memory_order_relaxed))
            break;
        size_t job = atomic_fetch_add_explicit(&s->next_job, 1,
                                               memory_order_relaxed);
        if (job >= s->nprimes)
            break;
        size_t index = s->nprimes - 1 - job; /* expensive primes first */
        uint32_t p = s->primes[index];
        uint32_t zcount = 0;
        if (!scan_prime(p, &worker->pairs, &zcount)) {
            atomic_store_explicit(&s->failed, 1, memory_order_relaxed);
            break;
        }
        s->zcounts[index] = zcount;
        atomic_fetch_add_explicit(&s->steps_done, p - 2,
                                  memory_order_relaxed);
        atomic_fetch_add_explicit(&s->primes_done, 1, memory_order_relaxed);
    }
    atomic_fetch_add_explicit(&s->workers_done, 1, memory_order_release);
    return NULL;
}

static int write_pairs_little_endian(const char *path,
                                     const Pair *pairs, size_t count) {
    FILE *fp = fopen(path, "wb");
    if (!fp) {
        fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
        return 0;
    }
    enum { RECORDS_PER_BLOCK = 4096 };
    uint8_t *buf = (uint8_t *)malloc(8 * RECORDS_PER_BLOCK);
    if (!buf) {
        fprintf(stderr, "out of memory for output buffer\n");
        fclose(fp);
        return 0;
    }
    size_t offset = 0;
    while (offset < count) {
        size_t nr = count - offset;
        if (nr > RECORDS_PER_BLOCK)
            nr = RECORDS_PER_BLOCK;
        for (size_t j = 0; j < nr; ++j) {
            uint32_t p = pairs[offset + j].p;
            uint32_t r = pairs[offset + j].r;
            uint8_t *dst = buf + 8 * j;
            dst[0] = (uint8_t)p;
            dst[1] = (uint8_t)(p >> 8);
            dst[2] = (uint8_t)(p >> 16);
            dst[3] = (uint8_t)(p >> 24);
            dst[4] = (uint8_t)r;
            dst[5] = (uint8_t)(r >> 8);
            dst[6] = (uint8_t)(r >> 16);
            dst[7] = (uint8_t)(r >> 24);
        }
        if (fwrite(buf, 8, nr, fp) != nr) {
            fprintf(stderr, "write failed for %s: %s\n", path,
                    strerror(errno));
            free(buf);
            fclose(fp);
            return 0;
        }
        offset += nr;
    }
    free(buf);
    if (fclose(fp) != 0) {
        fprintf(stderr, "close failed for %s: %s\n", path, strerror(errno));
        return 0;
    }
    return 1;
}

static int validate_prefix_200000(const uint32_t *b, uint32_t nmax,
                                  uint32_t pmax) {
    if (nmax < 200000 || pmax < 200000)
        return 1;
    static const uint64_t expected[4] = {187494, 12094, 400, 12};
    uint64_t hist[4] = {0, 0, 0, 0};
    uint32_t max_b = 0;
    for (uint32_t n = 1; n <= 200000; ++n) {
        if (b[n] > max_b)
            max_b = b[n];
        if (b[n] < 4)
            ++hist[b[n]];
    }
    int ok = max_b == 3;
    for (int k = 0; k < 4; ++k)
        ok = ok && hist[k] == expected[k];
    printf("\n=== Mandatory n <= 200000 gate ===\n");
    printf("B=0: %" PRIu64 " (%.3f%%)\n", hist[0], hist[0] / 2000.0);
    printf("B=1: %" PRIu64 " (%.3f%%)\n", hist[1], hist[1] / 2000.0);
    printf("B=2: %" PRIu64 " (%.3f%%)\n", hist[2], hist[2] / 2000.0);
    printf("B=3: %" PRIu64 " (%.3f%%)\n", hist[3], hist[3] / 2000.0);
    printf("max B: %" PRIu32 "\n", max_b);
    printf("Validation gate: %s\n", ok ? "PASS" : "FAIL");
    if (!ok)
        fprintf(stderr, "mandatory 200000-prefix validation failed\n");
    return ok;
}

static int validate_prime_prefix_1000000(const uint32_t *primes,
                                         const uint32_t *zcounts,
                                         size_t nprimes, uint32_t pmax) {
    if (pmax < 1000000)
        return 1;
    static const uint64_t expected[13] = {
        47632, 2, 23729, 0, 6045, 0, 951,
        0, 123, 0, 12, 0, 1
    };
    uint64_t hist[13] = {0};
    uint64_t pairs = 0;
    uint32_t max_z = 0, max_p = 0;
    size_t prime_count = 0;
    int ok = 1;
    for (size_t i = 0; i < nprimes && primes[i] <= 1000000; ++i) {
        uint32_t z = zcounts[i];
        ++prime_count;
        pairs += z;
        if (z > max_z) {
            max_z = z;
            max_p = primes[i];
        }
        if (z >= 13)
            ok = 0;
        else
            ++hist[z];
    }
    ok = ok && prime_count == 78495 && pairs == 78462 &&
         max_z == 12 && max_p == 159977;
    for (int z = 0; z <= 12; ++z)
        ok = ok && hist[z] == expected[z];
    printf("Known p <= 1000000 Z(p) cross-check: %s "
           "(primes=%zu, pairs=%" PRIu64 ", max=%" PRIu32
           " at p=%" PRIu32 ")\n",
           ok ? "PASS" : "FAIL", prime_count, pairs, max_z, max_p);
    if (!ok)
        fprintf(stderr, "known million-prime prefix validation failed\n");
    return ok;
}

static void print_integer_list_for_value(const uint32_t *b, uint32_t nmax,
                                         uint32_t value) {
    int first = 1;
    for (uint32_t n = 1; n <= nmax; ++n) {
        if (b[n] == value) {
            printf("%s%" PRIu32, first ? "" : ", ", n);
            first = 0;
        }
    }
    if (first)
        printf("none");
    putchar('\n');
}

static void print_ge_four_list(const uint32_t *b, uint32_t nmax) {
    int first = 1;
    for (uint32_t n = 1; n <= nmax; ++n) {
        if (b[n] >= 4) {
            printf("%s%" PRIu32 ":%" PRIu32,
                   first ? "" : ", ", n, b[n]);
            first = 0;
        }
    }
    if (first)
        printf("none");
    putchar('\n');
}

static void print_global_statistics(const uint32_t *b, uint32_t nmax) {
    uint32_t max_b = 0;
    for (uint32_t n = 1; n <= nmax; ++n)
        if (b[n] > max_b)
            max_b = b[n];
    uint64_t *hist = (uint64_t *)calloc((size_t)max_b + 1, sizeof(*hist));
    if (!hist) {
        fprintf(stderr, "out of memory for B histogram\n");
        return;
    }
    uint64_t argmax_count = 0, ge_four_count = 0;
    long double sum = 0, sumsq = 0;
    for (uint32_t n = 1; n <= nmax; ++n) {
        ++hist[b[n]];
        sum += b[n];
        sumsq += (long double)b[n] * b[n];
        if (b[n] == max_b)
            ++argmax_count;
        if (b[n] >= 4)
            ++ge_four_count;
    }
    long double mean = sum / nmax;
    long double variance = sumsq / nmax - mean * mean;
    printf("\n=== Global B(n), 1 <= n <= %" PRIu32 " ===\n", nmax);
    printf("HEADLINE max B(n) = %" PRIu32 " (argmax count=%" PRIu64 ")\n",
           max_b, argmax_count);
    printf("Argmax n: ");
    print_integer_list_for_value(b, nmax, max_b);
    printf("All n with B(n)>=4 (count=%" PRIu64 "): ", ge_four_count);
    print_ge_four_list(b, nmax);
    printf("Histogram:\n");
    for (uint32_t k = 0; k <= max_b; ++k)
        printf("  B=%" PRIu32 ": %" PRIu64 " (%.6f%%)\n",
               k, hist[k], 100.0 * (double)hist[k] / nmax);
    printf("mean=%.12Lf variance(population)=%.12Lf Var/Mean=%.12Lf\n",
           mean, variance, mean > 0 ? variance / mean : 0);
    free(hist);
}

static void print_z_statistics(const uint32_t *primes,
                               const uint32_t *zcounts, size_t nprimes,
                               size_t pair_count) {
    uint32_t max_z = 0;
    for (size_t i = 0; i < nprimes; ++i)
        if (zcounts[i] > max_z)
            max_z = zcounts[i];
    uint64_t *hist = (uint64_t *)calloc((size_t)max_z + 1, sizeof(*hist));
    if (!hist) {
        fprintf(stderr, "out of memory for Z histogram\n");
        return;
    }
    uint64_t sum_z = 0;
    uint64_t odd_count = 0;
    for (size_t i = 0; i < nprimes; ++i) {
        ++hist[zcounts[i]];
        sum_z += zcounts[i];
        odd_count += zcounts[i] & 1U;
    }
    printf("\n=== Per-prime Z(p) statistics ===\n");
    printf("primes=%zu pair records=%zu mean Z=%.12f max Z=%" PRIu32 "\n",
           nprimes, pair_count, nprimes ? (double)sum_z / nprimes : 0,
           max_z);
    printf("Histogram:\n");
    for (uint32_t z = 0; z <= max_z; ++z) {
        if (hist[z])
            printf("  Z=%" PRIu32 ": %" PRIu64 " (%.6f%%)\n",
                   z, hist[z], 100.0 * (double)hist[z] / nprimes);
    }
    printf("Prime(s) attaining max Z: ");
    int first = 1;
    for (size_t i = 0; i < nprimes; ++i) {
        if (zcounts[i] == max_z) {
            printf("%s%" PRIu32, first ? "" : ", ", primes[i]);
            first = 0;
        }
    }
    putchar('\n');
    uint64_t ordinary_count = nprimes - odd_count;
    printf("Odd-Z exceptional primes (%" PRIu64 "): ", odd_count);
    int odd_first = 1;
    for (size_t i = 0; i < nprimes; ++i) {
        if (zcounts[i] & 1U) {
            printf("%s%" PRIu32 " (Z=%" PRIu32 ")",
                   odd_first ? "" : ", ", primes[i], zcounts[i]);
            odd_first = 0;
        }
    }
    if (odd_first)
        printf("none");
    putchar('\n');
    printf("Poisson(1/2) pair-count comparison "
           "(denominator=%" PRIu64 " ordinary even-Z primes):\n",
           ordinary_count);
    double probability = exp(-0.5);
    for (uint32_t k = 0; 2 * k <= max_z; ++k) {
        uint32_t z = 2 * k;
        printf("  K=%" PRIu32 " (Z=%" PRIu32 "): obs=%.9f pred=%.9f\n",
               k, z, ordinary_count ? (double)hist[z] / ordinary_count : 0,
               probability);
        probability *= 0.5 / (k + 1);
    }
    if (sum_z != pair_count)
        fprintf(stderr, "internal error: sum Z != pair count\n");
    free(hist);
}

static void print_dyadic_statistics(const uint32_t *b, const double *w,
                                    uint32_t nmax) {
    printf("\n=== Dyadic B histogram vs Poisson(log(2)/log(n)) mixture ===\n");
    printf("N hi size partial lambda_bar "
           "obs0 pred0 obs1 pred1 obs2 pred2 obs3 pred3 obs4+ pred4+\n");
    for (uint64_t lower = DYADIC_START; lower < nmax; lower *= 2) {
        uint64_t upper = 2 * lower;
        int partial = upper > nmax;
        if (partial)
            upper = nmax;
        uint64_t size = upper - lower;
        uint64_t observed[5] = {0, 0, 0, 0, 0};
        long double predicted[5] = {0, 0, 0, 0, 0};
        long double lambda_sum = 0;
        for (uint64_t n = lower + 1; n <= upper; ++n) {
            uint32_t bin = b[n] < 4 ? b[n] : 4;
            ++observed[bin];
            double lambda = LN2_CONST / log((double)n);
            double pk = exp(-lambda);
            double used = 0;
            for (int k = 0; k < 4; ++k) {
                predicted[k] += pk;
                used += pk;
                pk *= lambda / (k + 1);
            }
            predicted[4] += 1.0 - used;
            lambda_sum += lambda;
        }
        printf("%" PRIu64 " %" PRIu64 " %" PRIu64 " %s %.9Lf",
               lower, upper, size, partial ? "yes" : "no",
               lambda_sum / size);
        for (int k = 0; k < 5; ++k)
            printf(" %.9f %.9Lf",
                   (double)observed[k] / size, predicted[k] / size);
        putchar('\n');
        if (partial)
            break;
    }

    printf("\n=== Dyadic moments and max W_top(n)/n ===\n");
    printf("N hi size partial mean_B variance_B Var/Mean "
           "max_W_over_n argmax_n W_top B_at_n\n");
    for (uint64_t lower = DYADIC_START; lower < nmax; lower *= 2) {
        uint64_t upper = 2 * lower;
        int partial = upper > nmax;
        if (partial)
            upper = nmax;
        uint64_t size = upper - lower;
        long double sum = 0, sumsq = 0;
        double max_ratio = -1;
        uint32_t argmax = (uint32_t)(lower + 1);
        for (uint64_t n = lower + 1; n <= upper; ++n) {
            sum += b[n];
            sumsq += (long double)b[n] * b[n];
            double ratio = w[n] / (double)n;
            if (ratio > max_ratio) {
                max_ratio = ratio;
                argmax = (uint32_t)n;
            }
        }
        long double mean = sum / size;
        long double variance = sumsq / size - mean * mean;
        printf("%" PRIu64 " %" PRIu64 " %" PRIu64 " %s "
               "%.12Lf %.12Lf %.12Lf %.12g %" PRIu32 " %.12g %" PRIu32 "\n",
               lower, upper, size, partial ? "yes" : "no",
               mean, variance, mean > 0 ? variance / mean : 0,
               max_ratio, argmax, w[argmax], b[argmax]);
        if (partial)
            break;
    }
}

static void print_window_statistics(const uint32_t *b, uint32_t nmax) {
    printf("\n=== 64-window localized dispersion (complete shells only) ===\n");
    printf("N window_length mean_window_sum model_mean_window_sum "
           "sample_variance variance/observed_mean variance/model_mean "
           "max_sum max_window\n");
    for (uint64_t lower = WINDOW_START; 2 * lower <= nmax; lower *= 2) {
        uint64_t length = lower / 64;
        uint64_t sums[64] = {0};
        long double model_total = 0;
        for (int j = 0; j < 64; ++j) {
            uint64_t lo = lower + (uint64_t)j * length;
            uint64_t hi = lo + length;
            for (uint64_t n = lo + 1; n <= hi; ++n) {
                sums[j] += b[n];
                model_total += LN2_CONST / log((double)n);
            }
        }
        uint64_t total = 0, max_sum = 0;
        int max_j = 0;
        for (int j = 0; j < 64; ++j) {
            total += sums[j];
            if (sums[j] > max_sum) {
                max_sum = sums[j];
                max_j = j;
            }
        }
        long double mean = (long double)total / 64;
        long double variance = 0;
        for (int j = 0; j < 64; ++j) {
            long double d = sums[j] - mean;
            variance += d * d;
        }
        variance /= 63; /* unbiased: iid Poisson prediction is the mean */
        uint64_t max_lo = lower + (uint64_t)max_j * length;
        uint64_t max_hi = max_lo + length;
        long double model_mean = model_total / 64;
        printf("%" PRIu64 " %" PRIu64 " %.9Lf %.9Lf %.9Lf %.9Lf "
               "%.9Lf %" PRIu64 " (%" PRIu64 ",%" PRIu64 "]\n",
               lower, length, mean, model_mean, variance,
               mean > 0 ? variance / mean : 0,
               model_mean > 0 ? variance / model_mean : 0,
               max_sum, max_lo, max_hi);
    }
}

static int parse_u32(const char *text, uint32_t *value) {
    char *end = NULL;
    errno = 0;
    unsigned long long x = strtoull(text, &end, 10);
    if (errno || !end || *end || x > UINT32_MAX)
        return 0;
    *value = (uint32_t)x;
    return 1;
}

int main(int argc, char **argv) {
    uint32_t nmax = DEFAULT_MAX;
    uint32_t pmax = DEFAULT_MAX;
    long detected = sysconf(_SC_NPROCESSORS_ONLN);
    uint32_t nthreads = detected > 0 ? (uint32_t)detected : 1;
    const char *pair_path = "data_zp_pairs.bin";

    if (argc > 5) {
        fprintf(stderr,
                "usage: %s [N_MAX [P_MAX [THREADS [PAIR_FILE]]]]\n",
                argv[0]);
        return EXIT_FAILURE;
    }
    if (argc >= 2 && !parse_u32(argv[1], &nmax)) {
        fprintf(stderr, "invalid N_MAX: %s\n", argv[1]);
        return EXIT_FAILURE;
    }
    pmax = nmax;
    if (argc >= 3 && !parse_u32(argv[2], &pmax)) {
        fprintf(stderr, "invalid P_MAX: %s\n", argv[2]);
        return EXIT_FAILURE;
    }
    if (argc >= 4 && !parse_u32(argv[3], &nthreads)) {
        fprintf(stderr, "invalid THREADS: %s\n", argv[3]);
        return EXIT_FAILURE;
    }
    if (argc >= 5)
        pair_path = argv[4];
    if (nmax < 1 || pmax < FIRST_PRIME || nthreads < 1) {
        fprintf(stderr, "require N_MAX>=1, P_MAX>=7, THREADS>=1\n");
        return EXIT_FAILURE;
    }
    if (pmax < nmax) {
        fprintf(stderr,
                "require P_MAX>=N_MAX for a complete top-half B(n) scan\n");
        return EXIT_FAILURE;
    }

    setvbuf(stdout, NULL, _IOLBF, 0);
    printf("bn_bigscan: N_MAX=%" PRIu32 " P_MAX=%" PRIu32
           " threads=%" PRIu32 " (online CPUs=%ld)\n",
           nmax, pmax, nthreads, detected);
    printf("pair output: %s\n", pair_path);
    if (!startup_self_test())
        return EXIT_FAILURE;

    double program_t0 = monotonic_seconds();
    size_t nprimes = 0;
    double sieve_seconds = 0;
    uint32_t *primes = sieve_primes(pmax, &nprimes, &sieve_seconds);
    if (!primes)
        return EXIT_FAILURE;
    uint64_t total_steps = 0;
    for (size_t i = 0; i < nprimes; ++i)
        total_steps += primes[i] - 2;
    printf("sieve: %zu primes in [7,%" PRIu32 "] in %.3f s\n",
           nprimes, pmax, sieve_seconds);
    printf("recurrence steps: %" PRIu64 "\n", total_steps);
    if (nthreads > nprimes)
        nthreads = (uint32_t)nprimes;

    uint32_t *zcounts = (uint32_t *)calloc(nprimes, sizeof(*zcounts));
    Worker *workers = (Worker *)calloc(nthreads, sizeof(*workers));
    pthread_t *threads = (pthread_t *)malloc(nthreads * sizeof(*threads));
    if (!zcounts || !workers || !threads) {
        fprintf(stderr, "out of memory for scan state\n");
        free(primes);
        free(zcounts);
        free(workers);
        free(threads);
        return EXIT_FAILURE;
    }
    Shared shared = {
        .primes = primes,
        .nprimes = nprimes,
        .zcounts = zcounts
    };
    atomic_init(&shared.next_job, 0);
    atomic_init(&shared.primes_done, 0);
    atomic_init(&shared.steps_done, 0);
    atomic_init(&shared.workers_done, 0);
    atomic_init(&shared.failed, 0);
    double scan_t0 = monotonic_seconds();
    uint32_t created = 0;
    for (; created < nthreads; ++created) {
        workers[created].shared = &shared;
        int err = pthread_create(&threads[created], NULL, worker_main,
                                 &workers[created]);
        if (err) {
            fprintf(stderr, "pthread_create: %s\n", strerror(err));
            atomic_store(&shared.failed, 1);
            break;
        }
    }
    if (created != nthreads) {
        for (uint32_t t = 0; t < created; ++t)
            pthread_join(threads[t], NULL);
        for (uint32_t t = 0; t < nthreads; ++t)
            free(workers[t].pairs.data);
        free(primes); free(zcounts); free(workers); free(threads);
        return EXIT_FAILURE;
    }

    int last_percent = -1;
    while (atomic_load_explicit(&shared.workers_done, memory_order_acquire) <
           (int)nthreads) {
        uint64_t done = atomic_load_explicit(&shared.steps_done,
                                             memory_order_relaxed);
        int percent = total_steps ? (int)(100 * done / total_steps) : 100;
        if (percent != last_percent && (percent % 2 == 0 || percent == 100)) {
            double elapsed = monotonic_seconds() - scan_t0;
            double eta = done ? elapsed * (double)(total_steps - done) / done : 0;
            size_t pdone = atomic_load_explicit(&shared.primes_done,
                                                memory_order_relaxed);
            printf("progress: %3d%%, primes=%zu/%zu, steps=%" PRIu64
                   ", elapsed=%.1fs, ETA=%.1fs\n",
                   percent, pdone, nprimes, done, elapsed, eta);
            last_percent = percent;
        }
        struct timespec pause = {.tv_sec = 1, .tv_nsec = 0};
        nanosleep(&pause, NULL);
    }
    for (uint32_t t = 0; t < nthreads; ++t)
        pthread_join(threads[t], NULL);
    double scan_seconds = monotonic_seconds() - scan_t0;
    if (atomic_load(&shared.failed)) {
        fprintf(stderr, "scan failed (probably out of memory)\n");
        for (uint32_t t = 0; t < nthreads; ++t)
            free(workers[t].pairs.data);
        free(primes); free(zcounts); free(workers); free(threads);
        return EXIT_FAILURE;
    }
    printf("kernel scan complete in %.3f s (%.3f billion steps/s)\n",
           scan_seconds, scan_seconds > 0 ? total_steps / scan_seconds / 1e9 : 0);

    size_t pair_count = 0;
    for (uint32_t t = 0; t < nthreads; ++t) {
        if (SIZE_MAX - pair_count < workers[t].pairs.len) {
            fprintf(stderr, "pair count overflow\n");
            return EXIT_FAILURE;
        }
        pair_count += workers[t].pairs.len;
    }
    Pair *pairs = (Pair *)malloc(pair_count * sizeof(*pairs));
    if (!pairs && pair_count) {
        fprintf(stderr, "out of memory merging pairs\n");
        return EXIT_FAILURE;
    }
    size_t offset = 0;
    for (uint32_t t = 0; t < nthreads; ++t) {
        memcpy(pairs + offset, workers[t].pairs.data,
               workers[t].pairs.len * sizeof(*pairs));
        offset += workers[t].pairs.len;
        free(workers[t].pairs.data);
    }
    qsort(pairs, pair_count, sizeof(*pairs), pair_cmp);
    for (size_t i = 0; i < pair_count; ++i) {
        if (pairs[i].r >= pairs[i].p ||
            (i && pair_cmp(&pairs[i - 1], &pairs[i]) >= 0)) {
            fprintf(stderr, "invalid or duplicate pair at record %zu\n", i);
            return EXIT_FAILURE;
        }
    }
    if (!write_pairs_little_endian(pair_path, pairs, pair_count))
        return EXIT_FAILURE;
    printf("wrote %zu sorted pair records (%zu bytes)\n",
           pair_count, 8 * pair_count);

    uint32_t *b = (uint32_t *)calloc((size_t)nmax + 1, sizeof(*b));
    double *w = (double *)calloc((size_t)nmax + 1, sizeof(*w));
    if (!b || !w) {
        fprintf(stderr, "out of memory for scatter arrays\n");
        return EXIT_FAILURE;
    }
    size_t scattered = 0;
    for (size_t i = 0; i < pair_count; ++i) {
        uint64_t n = (uint64_t)pairs[i].p + pairs[i].r;
        if (n <= nmax) {
            ++b[n];
            w[n] += log((double)pairs[i].p);
            ++scattered;
        }
    }
    printf("scattered %zu/%zu records to n <= %" PRIu32 "\n",
           scattered, pair_count, nmax);

    int validations_ok = validate_prefix_200000(b, nmax, pmax) &&
        validate_prime_prefix_1000000(primes, zcounts, nprimes, pmax);
    if (!validations_ok) {
        fprintf(stderr, "validation failed; suppressing statistical claims\n");
        return EXIT_FAILURE;
    }
    print_z_statistics(primes, zcounts, nprimes, pair_count);
    print_global_statistics(b, nmax);
    print_dyadic_statistics(b, w, nmax);
    print_window_statistics(b, nmax);
    printf("\nTotal program wall time: %.3f s\n",
           monotonic_seconds() - program_t0);

    free(b);
    free(w);
    free(pairs);
    free(primes);
    free(zcounts);
    free(workers);
    free(threads);
    return EXIT_SUCCESS;
}
