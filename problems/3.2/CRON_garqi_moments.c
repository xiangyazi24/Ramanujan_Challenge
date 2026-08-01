/*
 * CRON_garqi_moments.c -- dyadic factorial moments for Apéry zero rows.
 *
 * For each prime p <= 2^20, compute
 *
 *   Z_p = {0 <= r < p : b_r == 0 (mod p)},
 *
 * where b_0=1, b_1=5 and
 *
 *   (n+1)^3 b_{n+1}
 *     = (34n^3 + 51n^2 + 27n + 5)b_n - n^3 b_{n-1}.
 *
 * The hot loop uses c_n=(n!)^3 b_n.  Since n! is a unit modulo p for
 * n<p, c_n and b_n have exactly the same zero set, while c satisfies the
 * division-free recurrence
 *
 *   c_{n+1} = (34n^3 + 51n^2 + 27n + 5)c_n - n^6 c_{n-1}.
 *
 * Cubic and sixth-power coefficients are advanced by finite differences,
 * and products are reduced with 32-bit Montgomery arithmetic.  An
 * independent startup implementation of the original b recurrence, using
 * a table of modular inverses, checks the optimized kernel.
 *
 * For X in {4000, 8000, 2^13, ..., 2^19}, every zero pair is scattered to
 *
 *   H_X(n) = #{p in (X/2,2X] : n=p+r, r in Z_p},  X<n<=2X.
 *
 * The X=4000 and X=8000 results are mandatory exact sanity gates.  A second
 * mandatory gate compares 20 deterministic pseudorandom H_X(n) values at
 * X=8192 with direct tests of 5*b_{n-p} modulo p.
 *
 * Compile and run (the optional argument is the thread count):
 *
 *   cc -O3 -march=native -Wall -Wextra -Wpedantic \
 *      -o /tmp/CRON_garqi_moments CRON_garqi_moments.c -lpthread -lm
 *   /tmp/CRON_garqi_moments 8
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

#define DEFAULT_THREADS 8U
#define MAX_X 524288U
#define PRIME_LIMIT (2U * MAX_X + 1U)
#define MIDPOINT_LIMIT 1000000U
#define TARGET_COUNT 9U
#define RANDOM_CHECK_COUNT 20U

typedef struct {
    uint32_t p;
    uint32_t z;
} ZeroPair;

typedef struct {
    ZeroPair *data;
    size_t len;
    size_t cap;
} PairVec;

typedef struct {
    const uint32_t *primes;
    size_t nprimes;
    atomic_size_t next_job;
    atomic_size_t primes_done;
    atomic_uint_fast64_t steps_done;
    atomic_int failed;
} SharedScan;

typedef struct {
    SharedScan *shared;
    PairVec zeros;
} Worker;

typedef struct {
    uint32_t X;
    uint16_t *hits;
    uint64_t S;
    uint64_t moments[7];
    long double ratios[7];
    uint32_t max_h;
    uint32_t first_argmax;
    uint32_t argmax_count;
} Target;

static const uint32_t target_sizes[TARGET_COUNT] = {
    4000U, 8000U, 8192U, 16384U, 32768U,
    65536U, 131072U, 262144U, 524288U
};

static double monotonic_seconds(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
        perror("clock_gettime");
        exit(EXIT_FAILURE);
    }
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int parse_threads(const char *text, uint32_t *value) {
    char *end = NULL;
    errno = 0;
    unsigned long parsed = strtoul(text, &end, 10);
    if (errno != 0 || end == text || end == NULL || *end != '\0' ||
        parsed == 0 || parsed > UINT32_MAX)
        return 0;
    *value = (uint32_t)parsed;
    return 1;
}

static int pair_push(PairVec *v, uint32_t p, uint32_t z) {
    if (v->len == v->cap) {
        size_t next_cap = v->cap == 0 ? 1024U : 2U * v->cap;
        if (next_cap < v->cap || next_cap > SIZE_MAX / sizeof(*v->data))
            return 0;
        ZeroPair *next =
            (ZeroPair *)realloc(v->data, next_cap * sizeof(*v->data));
        if (next == NULL)
            return 0;
        v->data = next;
        v->cap = next_cap;
    }
    v->data[v->len++] = (ZeroPair){p, z};
    return 1;
}

static uint32_t *sieve_primes(size_t *count_out) {
    uint8_t *composite = (uint8_t *)calloc(PRIME_LIMIT, 1);
    if (composite == NULL) {
        fprintf(stderr, "fatal: cannot allocate sieve\n");
        return NULL;
    }
    for (uint32_t d = 2; (uint64_t)d * d < PRIME_LIMIT; ++d) {
        if (composite[d] != 0)
            continue;
        for (uint32_t m = d * d; m < PRIME_LIMIT; m += d)
            composite[m] = 1;
    }

    size_t count = 0;
    for (uint32_t p = 7; p < PRIME_LIMIT; ++p)
        if (composite[p] == 0)
            ++count;

    uint32_t *primes = (uint32_t *)malloc(count * sizeof(*primes));
    if (primes == NULL) {
        fprintf(stderr, "fatal: cannot allocate prime list\n");
        free(composite);
        return NULL;
    }
    size_t out = 0;
    for (uint32_t p = 7; p < PRIME_LIMIT; ++p)
        if (composite[p] == 0)
            primes[out++] = p;
    free(composite);
    *count_out = count;
    return primes;
}

/* Return -p^{-1} modulo 2^32 for odd p. */
static inline uint32_t montgomery_nprime(uint32_t p) {
    uint32_t inverse = 1;
    inverse *= 2U - p * inverse;
    inverse *= 2U - p * inverse;
    inverse *= 2U - p * inverse;
    inverse *= 2U - p * inverse;
    inverse *= 2U - p * inverse;
    return 0U - inverse;
}

/* Inputs and output are Montgomery residues, with R=2^32. */
static inline uint32_t montgomery_mul(uint32_t a, uint32_t b,
                                      uint32_t p, uint32_t nprime) {
    uint64_t product = (uint64_t)a * b;
    uint32_t correction = (uint32_t)product * nprime;
    uint64_t reduced =
        (product + (uint64_t)correction * p) >> 32;
    uint32_t result = (uint32_t)reduced;
    return result >= p ? result - p : result;
}

static inline uint32_t to_montgomery(uint32_t value, uint32_t p) {
    return (uint32_t)(((uint64_t)(value % p) << 32) % p);
}

static inline uint32_t add_mod(uint32_t a, uint32_t b, uint32_t p) {
    uint32_t sum = a + b;
    return sum >= p ? sum - p : sum;
}

static inline uint32_t sub_mod(uint32_t a, uint32_t b, uint32_t p) {
    return a >= b ? a - b : a + p - b;
}

static inline void advance_differences(uint32_t *row, int degree,
                                       uint32_t p) {
    for (int order = 0; order < degree; ++order)
        row[order] = add_mod(row[order], row[order + 1], p);
}

/*
 * Forward-difference rows at n=1, in increasing difference order:
 *
 *   34n^3+51n^2+27n+5 : 117, 418, 510, 204
 *   n^6                  :   1,  63, 602, 2100, 3360, 2520, 720.
 */
static int scan_prime_fast(uint32_t p, PairVec *zeros) {
    static const uint32_t cubic_initial[4] = {117U, 418U, 510U, 204U};
    static const uint32_t sixth_initial[7] = {
        1U, 63U, 602U, 2100U, 3360U, 2520U, 720U
    };
    uint32_t cubic[4];
    uint32_t sixth[7];
    for (int i = 0; i < 4; ++i)
        cubic[i] = to_montgomery(cubic_initial[i], p);
    for (int i = 0; i < 7; ++i)
        sixth[i] = to_montgomery(sixth_initial[i], p);

    uint32_t nprime = montgomery_nprime(p);
    uint32_t previous = to_montgomery(1U, p);
    uint32_t current = to_montgomery(5U, p);

    for (uint32_t n = 1; n < p - 1; ++n) {
        uint32_t first =
            montgomery_mul(cubic[0], current, p, nprime);
        uint32_t second =
            montgomery_mul(sixth[0], previous, p, nprime);
        uint32_t next = sub_mod(first, second, p);
        if (next == 0U && !pair_push(zeros, p, n + 1U))
            return 0;
        previous = current;
        current = next;
        advance_differences(cubic, 3, p);
        advance_differences(sixth, 6, p);
    }
    return 1;
}

/* Independent reference kernel: the original b recurrence with inv[1..p-1]. */
static int scan_prime_reference(uint32_t p, PairVec *zeros) {
    uint32_t *inverse = (uint32_t *)malloc((size_t)p * sizeof(*inverse));
    if (inverse == NULL)
        return 0;
    inverse[1] = 1U;
    for (uint32_t k = 2; k < p; ++k) {
        uint64_t correction =
            (uint64_t)(p / k) * inverse[p % k] % p;
        inverse[k] = p - (uint32_t)correction;
    }

    uint32_t previous = 1U % p;
    uint32_t current = 5U % p;
    for (uint32_t n = 1; n < p - 1; ++n) {
        uint64_t n2 = (uint64_t)n * n % p;
        uint64_t n3 = n2 * n % p;
        uint64_t coefficient =
            (34U * n3 + 51U * n2 + 27ULL * n + 5U) % p;
        uint32_t numerator = sub_mod(
            (uint32_t)(coefficient * current % p),
            (uint32_t)(n3 * previous % p), p);
        uint64_t inv = inverse[n + 1U];
        uint64_t inv_cube = inv * inv % p * inv % p;
        uint32_t next = (uint32_t)((uint64_t)numerator * inv_cube % p);
        if (next == 0U && !pair_push(zeros, p, n + 1U)) {
            free(inverse);
            return 0;
        }
        previous = current;
        current = next;
    }
    free(inverse);
    return 1;
}

static int startup_kernel_gate(void) {
    static const uint32_t test_primes[] = {
        7U, 11U, 17U, 31U, 181U, 379U, 3137U
    };
    for (size_t i = 0;
         i < sizeof(test_primes) / sizeof(test_primes[0]); ++i) {
        PairVec fast = {0};
        PairVec reference = {0};
        uint32_t p = test_primes[i];
        int ok = scan_prime_fast(p, &fast) &&
                 scan_prime_reference(p, &reference) &&
                 fast.len == reference.len &&
                 (fast.len == 0 ||
                  memcmp(fast.data, reference.data,
                         fast.len * sizeof(*fast.data)) == 0);
        free(fast.data);
        free(reference.data);
        if (!ok) {
            fprintf(stderr,
                    "SANITY ABORT: optimized/reference Z_p mismatch at p=%"
                    PRIu32 "\n", p);
            return 0;
        }
    }
    printf("SANITY kernel: PASS (original recurrence vs optimized kernel, "
           "7 primes through 3137)\n");
    return 1;
}

static void *worker_main(void *argument) {
    Worker *worker = (Worker *)argument;
    SharedScan *shared = worker->shared;
    for (;;) {
        if (atomic_load_explicit(&shared->failed, memory_order_relaxed))
            break;
        size_t job = atomic_fetch_add_explicit(
            &shared->next_job, 1U, memory_order_relaxed);
        if (job >= shared->nprimes)
            break;
        size_t prime_index = shared->nprimes - 1U - job;
        uint32_t p = shared->primes[prime_index];
        if (!scan_prime_fast(p, &worker->zeros)) {
            atomic_store_explicit(&shared->failed, 1, memory_order_relaxed);
            break;
        }
        atomic_fetch_add_explicit(
            &shared->steps_done, (uint_fast64_t)(p - 2U),
            memory_order_relaxed);
        atomic_fetch_add_explicit(
            &shared->primes_done, 1U, memory_order_relaxed);
    }
    return NULL;
}

static int pair_compare(const void *left, const void *right) {
    const ZeroPair *a = (const ZeroPair *)left;
    const ZeroPair *b = (const ZeroPair *)right;
    if (a->p != b->p)
        return a->p < b->p ? -1 : 1;
    if (a->z != b->z)
        return a->z < b->z ? -1 : 1;
    return 0;
}

static int checked_add_u64(uint64_t *accumulator, uint64_t value) {
    if (UINT64_MAX - *accumulator < value)
        return 0;
    *accumulator += value;
    return 1;
}

static int analyze_targets(const ZeroPair *pairs, size_t npairs,
                           Target targets[TARGET_COUNT]) {
    for (size_t t = 0; t < TARGET_COUNT; ++t) {
        targets[t].X = target_sizes[t];
        targets[t].hits =
            (uint16_t *)calloc(targets[t].X, sizeof(*targets[t].hits));
        if (targets[t].hits == NULL) {
            fprintf(stderr, "fatal: cannot allocate H array for X=%" PRIu32
                    "\n", targets[t].X);
            return 0;
        }
    }

    for (size_t i = 0; i < npairs; ++i) {
        uint32_t p = pairs[i].p;
        uint64_t n = (uint64_t)p + pairs[i].z;
        for (size_t t = 0; t < TARGET_COUNT; ++t) {
            uint32_t X = targets[t].X;
            if (2ULL * p <= X || p > 2U * X || n <= X || n > 2U * X)
                continue;
            size_t index = (size_t)(n - X - 1U);
            if (targets[t].hits[index] == UINT16_MAX) {
                fprintf(stderr, "fatal: H overflow at X=%" PRIu32
                        " n=%" PRIu64 "\n", X, n);
                return 0;
            }
            ++targets[t].hits[index];
        }
    }

    for (size_t t = 0; t < TARGET_COUNT; ++t) {
        Target *target = &targets[t];
        target->first_argmax = target->X + 1U;
        for (uint32_t offset = 0; offset < target->X; ++offset) {
            uint32_t h = target->hits[offset];
            if (!checked_add_u64(&target->S, h)) {
                fprintf(stderr, "fatal: S_X overflow at X=%" PRIu32 "\n",
                        target->X);
                return 0;
            }
            if (h > target->max_h) {
                target->max_h = h;
                target->first_argmax = target->X + 1U + offset;
                target->argmax_count = 1U;
            } else if (h == target->max_h) {
                ++target->argmax_count;
            }

            uint64_t falling = 1U;
            for (uint32_t k = 1; k <= 6; ++k) {
                if (h < k) {
                    falling = 0U;
                } else {
                    uint32_t factor = h - k + 1U;
                    if (factor != 0U && falling > UINT64_MAX / factor) {
                        fprintf(stderr,
                                "fatal: falling-factorial overflow at X=%"
                                PRIu32 "\n", target->X);
                        return 0;
                    }
                    falling *= factor;
                }
                if (!checked_add_u64(&target->moments[k], falling)) {
                    fprintf(stderr, "fatal: M_%" PRIu32
                            " overflow at X=%" PRIu32 "\n",
                            k, target->X);
                    return 0;
                }
            }
        }
        if (target->moments[1] != target->S) {
            fprintf(stderr, "SANITY ABORT: M_1 != S_X at X=%" PRIu32 "\n",
                    target->X);
            return 0;
        }
        long double lambda = (long double)target->S / target->X;
        long double lambda_power = 1.0L;
        for (uint32_t k = 1; k <= 6; ++k) {
            lambda_power *= lambda;
            long double denominator = target->X * lambda_power;
            target->ratios[k] = denominator > 0.0L
                ? target->moments[k] / denominator
                : NAN;
        }
    }
    return 1;
}

static Target *find_target(Target targets[TARGET_COUNT], uint32_t X) {
    for (size_t t = 0; t < TARGET_COUNT; ++t)
        if (targets[t].X == X)
            return &targets[t];
    return NULL;
}

static int exact_reference_gates(Target targets[TARGET_COUNT]) {
    Target *x4000 = find_target(targets, 4000U);
    Target *x8000 = find_target(targets, 8000U);
    if (x4000 == NULL || x8000 == NULL)
        return 0;

    if (x4000->S != 311U || x4000->max_h != 2U ||
        x4000->moments[2] != 18U) {
        fprintf(stderr,
                "SANITY ABORT X=4000: got S=%" PRIu64 " maxH=%" PRIu32
                " M2=%" PRIu64 "; expected 311,2,18\n",
                x4000->S, x4000->max_h, x4000->moments[2]);
        return 0;
    }
    if (x8000->S != 622U || x8000->max_h != 3U ||
        x8000->moments[2] != 54U || x8000->moments[3] != 6U) {
        fprintf(stderr,
                "SANITY ABORT X=8000: got S=%" PRIu64 " maxH=%" PRIu32
                " M2=%" PRIu64 " M3=%" PRIu64
                "; expected 622,3,54,6\n",
                x8000->S, x8000->max_h,
                x8000->moments[2], x8000->moments[3]);
        return 0;
    }
    printf("SANITY X=4000: PASS (S=311, maxH=2, M2=18)\n");
    printf("SANITY X=8000: PASS (S=622, maxH=3, M2=54, M3=6)\n");
    return 1;
}

static uint64_t splitmix64(uint64_t *state) {
    uint64_t z = (*state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

static int u32_compare(const void *left, const void *right) {
    uint32_t a = *(const uint32_t *)left;
    uint32_t b = *(const uint32_t *)right;
    return a < b ? -1 : a > b ? 1 : 0;
}

static int random_direct_gate(const uint32_t *primes, size_t nprimes,
                              Target *target) {
    uint32_t samples[RANDOM_CHECK_COUNT];
    uint32_t direct[RANDOM_CHECK_COUNT] = {0};
    uint64_t rng = UINT64_C(0x47524151492d3230);
    size_t sample_count = 0;
    while (sample_count < RANDOM_CHECK_COUNT) {
        uint32_t n = target->X + 1U +
                     (uint32_t)(splitmix64(&rng) % target->X);
        int duplicate = 0;
        for (size_t i = 0; i < sample_count; ++i)
            if (samples[i] == n)
                duplicate = 1;
        if (!duplicate)
            samples[sample_count++] = n;
    }
    qsort(samples, RANDOM_CHECK_COUNT, sizeof(*samples), u32_compare);

    uint32_t *inverse =
        (uint32_t *)malloc((target->X + 1U) * sizeof(*inverse));
    uint32_t *brow =
        (uint32_t *)malloc((target->X + 1U) * sizeof(*brow));
    if (inverse == NULL || brow == NULL) {
        fprintf(stderr, "fatal: direct-gate allocation failed\n");
        free(inverse);
        free(brow);
        return 0;
    }

    for (size_t prime_index = 0; prime_index < nprimes; ++prime_index) {
        uint32_t p = primes[prime_index];
        uint32_t max_r = 0;
        int used = 0;
        for (size_t i = 0; i < RANDOM_CHECK_COUNT; ++i) {
            uint32_t n = samples[i];
            if (2ULL * p > n && p <= n) {
                uint32_t r = n - p;
                if (!used || r > max_r)
                    max_r = r;
                used = 1;
            }
        }
        if (!used)
            continue;

        if (max_r >= p) {
            fprintf(stderr, "SANITY ABORT: direct gate reached r>=p\n");
            free(inverse);
            free(brow);
            return 0;
        }
        if (max_r >= 1U) {
            inverse[1] = 1U;
            for (uint32_t k = 2; k <= max_r; ++k) {
                uint64_t correction =
                    (uint64_t)(p / k) * inverse[p % k] % p;
                inverse[k] = p - (uint32_t)correction;
            }
        }
        brow[0] = 1U % p;
        if (max_r >= 1U)
            brow[1] = 5U % p;
        for (uint32_t n = 1; n < max_r; ++n) {
            uint64_t n2 = (uint64_t)n * n % p;
            uint64_t n3 = n2 * n % p;
            uint64_t coefficient =
                (34U * n3 + 51U * n2 + 27ULL * n + 5U) % p;
            uint32_t numerator = sub_mod(
                (uint32_t)(coefficient * brow[n] % p),
                (uint32_t)(n3 * brow[n - 1U] % p), p);
            uint64_t inv = inverse[n + 1U];
            uint64_t inv_cube = inv * inv % p * inv % p;
            brow[n + 1U] =
                (uint32_t)((uint64_t)numerator * inv_cube % p);
        }
        for (size_t i = 0; i < RANDOM_CHECK_COUNT; ++i) {
            uint32_t n = samples[i];
            if (2ULL * p > n && p <= n) {
                uint32_t r = n - p;
                if ((uint64_t)5U * brow[r] % p == 0U)
                    ++direct[i];
            }
        }
    }

    printf("SANITY X=8192 direct Gessel checks (seed=0x47524151492d3230):\n");
    for (size_t i = 0; i < RANDOM_CHECK_COUNT; ++i) {
        uint32_t n = samples[i];
        uint32_t row = target->hits[n - target->X - 1U];
        printf("  n=%" PRIu32 " row_H=%" PRIu32 " direct_H=%" PRIu32
               " %s\n", n, row, direct[i],
               row == direct[i] ? "PASS" : "MISMATCH");
        if (row != direct[i]) {
            fprintf(stderr,
                    "SANITY ABORT: random direct mismatch at n=%" PRIu32
                    "\n", n);
            free(inverse);
            free(brow);
            return 0;
        }
    }
    printf("SANITY X=8192 random direct checks: PASS (20/20)\n");
    free(inverse);
    free(brow);
    return 1;
}

static void print_midpoint_primes(const ZeroPair *pairs, size_t npairs) {
    size_t count = 0;
    printf("MIDPOINT_ZERO_PRIMES_BELOW_1000000:");
    for (size_t i = 0; i < npairs; ++i) {
        uint32_t p = pairs[i].p;
        if (p >= MIDPOINT_LIMIT)
            break;
        if (pairs[i].z == (p - 1U) / 2U) {
            printf(" %" PRIu32, p);
            ++count;
        }
    }
    printf("\nMIDPOINT_ZERO_COUNT=%zu\n", count);
}

static void print_results(const Target targets[TARGET_COUNT]) {
    printf("\n=== DYADIC GARQI FACTORIAL MOMENTS ===\n");
    printf("Columns: X S lambda lambda_logX maxH first_argmax argmax_count\n");
    for (size_t t = 2; t < TARGET_COUNT; ++t) {
        const Target *target = &targets[t];
        long double lambda = (long double)target->S / target->X;
        printf("SUMMARY X=%" PRIu32 " S=%" PRIu64
               " lambda=%.15Lf lambda_logX=%.15Lf maxH=%" PRIu32
               " first_argmax=%" PRIu32 " argmax_count=%" PRIu32 "\n",
               target->X, target->S, lambda,
               lambda * logl((long double)target->X),
               target->max_h, target->first_argmax,
               target->argmax_count);
        printf("  k M_k ratio_to_X_lambda^k\n");
        for (uint32_t k = 1; k <= 6; ++k)
            printf("  %" PRIu32 " %" PRIu64 " %.15Lf\n",
                   k, target->moments[k], target->ratios[k]);
    }
}

static void free_targets(Target targets[TARGET_COUNT]) {
    for (size_t t = 0; t < TARGET_COUNT; ++t) {
        free(targets[t].hits);
        targets[t].hits = NULL;
    }
}

int main(int argc, char **argv) {
    uint32_t nthreads = DEFAULT_THREADS;
    if (argc > 2) {
        fprintf(stderr, "usage: %s [THREADS]\n", argv[0]);
        return EXIT_FAILURE;
    }
    if (argc == 2 && !parse_threads(argv[1], &nthreads)) {
        fprintf(stderr, "invalid thread count: %s\n", argv[1]);
        return EXIT_FAILURE;
    }

    setvbuf(stdout, NULL, _IOLBF, 0);
    printf("CRON_garqi_moments prime_limit=%u threads=%" PRIu32 "\n",
           PRIME_LIMIT - 1U, nthreads);
    if (!startup_kernel_gate())
        return EXIT_FAILURE;

    double total_start = monotonic_seconds();
    size_t nprimes = 0;
    uint32_t *primes = sieve_primes(&nprimes);
    if (primes == NULL)
        return EXIT_FAILURE;
    if (nthreads > nprimes)
        nthreads = (uint32_t)nprimes;

    uint64_t expected_steps = 0;
    for (size_t i = 0; i < nprimes; ++i)
        expected_steps += primes[i] - 2U;
    printf("primes=%zu expected_recurrence_steps=%" PRIu64 "\n",
           nprimes, expected_steps);

    Worker *workers = (Worker *)calloc(nthreads, sizeof(*workers));
    pthread_t *threads =
        (pthread_t *)malloc((size_t)nthreads * sizeof(*threads));
    if (workers == NULL || threads == NULL) {
        fprintf(stderr, "fatal: cannot allocate worker state\n");
        free(workers);
        free(threads);
        free(primes);
        return EXIT_FAILURE;
    }

    SharedScan shared = {
        .primes = primes,
        .nprimes = nprimes
    };
    atomic_init(&shared.next_job, 0U);
    atomic_init(&shared.primes_done, 0U);
    atomic_init(&shared.steps_done, 0U);
    atomic_init(&shared.failed, 0);

    double scan_start = monotonic_seconds();
    uint32_t created = 0;
    for (; created < nthreads; ++created) {
        workers[created].shared = &shared;
        int error = pthread_create(
            &threads[created], NULL, worker_main, &workers[created]);
        if (error != 0) {
            fprintf(stderr, "pthread_create: %s\n", strerror(error));
            atomic_store_explicit(&shared.failed, 1, memory_order_relaxed);
            break;
        }
    }

    double last_progress = scan_start - 5.0;
    while (!atomic_load_explicit(&shared.failed, memory_order_relaxed) &&
           atomic_load_explicit(&shared.primes_done, memory_order_relaxed) <
               nprimes) {
        struct timespec pause = {.tv_sec = 0, .tv_nsec = 200000000L};
        (void)nanosleep(&pause, NULL);
        double now = monotonic_seconds();
        if (now - last_progress >= 5.0) {
            uint64_t done = (uint64_t)atomic_load_explicit(
                &shared.steps_done, memory_order_relaxed);
            size_t finished = atomic_load_explicit(
                &shared.primes_done, memory_order_relaxed);
            double elapsed = now - scan_start;
            double fraction = expected_steps == 0
                ? 1.0 : (double)done / expected_steps;
            double eta = done == 0
                ? INFINITY : elapsed * (expected_steps - done) / done;
            fprintf(stderr,
                    "[progress] primes=%zu/%zu steps=%" PRIu64 "/%" PRIu64
                    " (%.2f%%) elapsed=%.1fs ETA=%.1fs\n",
                    finished, nprimes, done, expected_steps,
                    100.0 * fraction, elapsed, eta);
            last_progress = now;
        }
    }

    for (uint32_t i = 0; i < created; ++i) {
        int error = pthread_join(threads[i], NULL);
        if (error != 0) {
            fprintf(stderr, "pthread_join: %s\n", strerror(error));
            atomic_store_explicit(&shared.failed, 1, memory_order_relaxed);
        }
    }
    double scan_seconds = monotonic_seconds() - scan_start;
    uint64_t steps_done = (uint64_t)atomic_load_explicit(
        &shared.steps_done, memory_order_relaxed);
    size_t primes_done = atomic_load_explicit(
        &shared.primes_done, memory_order_relaxed);
    if (created != nthreads ||
        atomic_load_explicit(&shared.failed, memory_order_relaxed) ||
        primes_done != nprimes || steps_done != expected_steps) {
        fprintf(stderr,
                "fatal: incomplete scan: threads=%" PRIu32 "/%" PRIu32
                " primes=%zu/%zu steps=%" PRIu64 "/%" PRIu64 "\n",
                created, nthreads, primes_done, nprimes,
                steps_done, expected_steps);
        for (uint32_t i = 0; i < nthreads; ++i)
            free(workers[i].zeros.data);
        free(workers);
        free(threads);
        free(primes);
        return EXIT_FAILURE;
    }

    size_t npairs = 0;
    for (uint32_t i = 0; i < nthreads; ++i) {
        if (SIZE_MAX - npairs < workers[i].zeros.len) {
            fprintf(stderr, "fatal: zero-pair count overflow\n");
            return EXIT_FAILURE;
        }
        npairs += workers[i].zeros.len;
    }
    ZeroPair *pairs =
        npairs == 0 ? NULL : (ZeroPair *)malloc(npairs * sizeof(*pairs));
    if (pairs == NULL && npairs != 0) {
        fprintf(stderr, "fatal: cannot merge zero pairs\n");
        return EXIT_FAILURE;
    }
    size_t offset = 0;
    for (uint32_t i = 0; i < nthreads; ++i) {
        memcpy(pairs + offset, workers[i].zeros.data,
               workers[i].zeros.len * sizeof(*pairs));
        offset += workers[i].zeros.len;
        free(workers[i].zeros.data);
    }
    qsort(pairs, npairs, sizeof(*pairs), pair_compare);
    for (size_t i = 0; i < npairs; ++i) {
        if (pairs[i].p >= PRIME_LIMIT || pairs[i].z >= pairs[i].p ||
            (i > 0 && pair_compare(&pairs[i - 1U], &pairs[i]) >= 0)) {
            fprintf(stderr,
                    "fatal: invalid or duplicate zero pair at index %zu\n", i);
            free(pairs);
            free(workers);
            free(threads);
            free(primes);
            return EXIT_FAILURE;
        }
    }
    printf("scan_seconds=%.6f throughput_steps_per_second=%.6e "
           "zero_pairs=%zu\n",
           scan_seconds, steps_done / scan_seconds, npairs);

    Target targets[TARGET_COUNT] = {0};
    int ok = analyze_targets(pairs, npairs, targets) &&
             exact_reference_gates(targets);
    Target *x8192 = find_target(targets, 8192U);
    if (ok && x8192 != NULL)
        ok = random_direct_gate(primes, nprimes, x8192);
    else if (x8192 == NULL)
        ok = 0;

    if (ok) {
        print_midpoint_primes(pairs, npairs);
        print_results(targets);
        printf("\nALL_SANITY_GATES=PASS\n");
        printf("total_seconds=%.6f\n", monotonic_seconds() - total_start);
    }

    free_targets(targets);
    free(pairs);
    free(workers);
    free(threads);
    free(primes);
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
