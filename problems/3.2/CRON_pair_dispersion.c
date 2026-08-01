/*
 * CRON_pair_dispersion.c -- shell decomposition of centered Apéry pair sums.
 *
 * For X in {4000, 8000, 2^13, ..., 2^19}, this program constructs
 *
 *   A_{p,X} = {n in (X,2X] : n-p in Z_p},
 *   Z_p     = {0 <= r < p : b_r == 0 (mod p)},
 *
 * for primes X/2 < p <= 2X.  It reports the ordered-pair quantities
 *
 *   C_j = sum_{p != q, 2^j <= |p-q| < 2^(j+1)} |A_p intersect A_q|,
 *   E_j = (1/X) sum_{same p,q} |A_p| |A_q|,
 *   D_j = C_j - E_j.
 *
 * The E_j numerator is accumulated from prefix sums over sorted primes, so
 * no all-prime-pairs loop is used.  The hot zero-set scan is the transformed
 * c_n=(n!)^3 b_n Montgomery kernel from CRON_garqi_moments.c.
 *
 * Compile and run from problems/3.2:
 *
 *   cc -O3 -march=native -Wall -Wextra -Wpedantic \
 *      -o /tmp/CRON_pair_dispersion CRON_pair_dispersion.c -lpthread -lm
 *   /tmp/CRON_pair_dispersion 8
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
#define TARGET_COUNT 9U
#define FIRST_REPORT_TARGET 2U
#define MAX_SHELLS 32U
#define SPOTCHECK_COUNT 5U

typedef struct {
    uint32_t p;
    uint32_t r;
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
    uint64_t actual;
    uint64_t expected_numerator;
    uint64_t reflection_actual;
    uint64_t reflection_expected_numerator;
    uint64_t central_actual;
    uint64_t folded_match_actual;
} Shell;

typedef struct {
    uint32_t p_index;
    uint32_t q_index;
    uint32_t witness_n;
    uint32_t hit_count;
    uint32_t direct_count;
} Spotcheck;

typedef struct {
    uint32_t X;
    uint16_t *a;
    uint16_t *a_reflection;
    uint16_t *a_central;
    uint16_t *hit_count;
    size_t *hit_offset;
    uint32_t *hit_primes;
    uint64_t S;
    uint64_t M2;
    uint64_t expected_numerator;
    uint64_t reflection_expected_numerator;
    uint32_t nshells;
    Shell shells[MAX_SHELLS];
    Spotcheck spots[SPOTCHECK_COUNT];
    uint32_t nspots;
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

static int pair_push(PairVec *v, uint32_t p, uint32_t r) {
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
    v->data[v->len++] = (ZeroPair){p, r};
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
        uint32_t first = montgomery_mul(cubic[0], current, p, nprime);
        uint32_t second = montgomery_mul(sixth[0], previous, p, nprime);
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

/* Independent startup path: the original b recurrence and inverse table. */
static int scan_prime_reference(uint32_t p, PairVec *zeros) {
    uint32_t *inverse = (uint32_t *)malloc((size_t)p * sizeof(*inverse));
    if (inverse == NULL)
        return 0;
    inverse[1] = 1U;
    for (uint32_t k = 2; k < p; ++k) {
        uint64_t correction = (uint64_t)(p / k) * inverse[p % k] % p;
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
                    "SANITY ABORT: kernel mismatch at p=%" PRIu32 "\n", p);
            return 0;
        }
    }
    printf("SANITY kernel: PASS (original vs transformed recurrence)\n");
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
    if (a->r != b->r)
        return a->r < b->r ? -1 : 1;
    return 0;
}

static uint32_t shell_index(uint32_t difference) {
    uint32_t shell = 0;
    while (difference >>= 1U)
        ++shell;
    return shell;
}

static size_t lower_bound_prime(const uint32_t *primes, size_t nprimes,
                                uint64_t value) {
    size_t lo = 0;
    size_t hi = nprimes;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2U;
        if ((uint64_t)primes[mid] < value)
            lo = mid + 1U;
        else
            hi = mid;
    }
    return lo;
}

static int zero_contains(const ZeroPair *pairs, const size_t *zero_offset,
                         size_t prime_index, uint32_t residue) {
    size_t lo = zero_offset[prime_index];
    size_t hi = zero_offset[prime_index + 1U];
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2U;
        if (pairs[mid].r < residue)
            lo = mid + 1U;
        else
            hi = mid;
    }
    return lo < zero_offset[prime_index + 1U] && pairs[lo].r == residue;
}

static int validate_zero_pairs(const ZeroPair *pairs, size_t npairs,
                               const uint32_t *primes, size_t nprimes,
                               size_t *zero_offset) {
    size_t cursor = 0;
    for (size_t i = 0; i < nprimes; ++i) {
        zero_offset[i] = cursor;
        while (cursor < npairs && pairs[cursor].p == primes[i]) {
            if (pairs[cursor].r >= primes[i] ||
                (cursor > zero_offset[i] &&
                 pairs[cursor - 1U].r >= pairs[cursor].r)) {
                fprintf(stderr,
                        "SANITY ABORT: invalid zero record p=%" PRIu32
                        " r=%" PRIu32 "\n",
                        pairs[cursor].p, pairs[cursor].r);
                return 0;
            }
            ++cursor;
        }
    }
    zero_offset[nprimes] = cursor;
    if (cursor != npairs) {
        fprintf(stderr, "SANITY ABORT: zero record has unknown prime\n");
        return 0;
    }
    for (size_t i = 0; i < nprimes; ++i) {
        uint32_t p = primes[i];
        for (size_t k = zero_offset[i]; k < zero_offset[i + 1U]; ++k) {
            uint32_t reflected = p - 1U - pairs[k].r;
            if (!zero_contains(pairs, zero_offset, i, reflected)) {
                fprintf(stderr,
                        "SANITY ABORT: reflection missing at p=%" PRIu32
                        " r=%" PRIu32 "\n", p, pairs[k].r);
                return 0;
            }
        }
    }
    printf("SANITY zero records: PASS (sorted, unique, in range, reflection)\n");
    return 1;
}

static int allocate_targets(Target targets[TARGET_COUNT], size_t nprimes) {
    for (size_t t = 0; t < TARGET_COUNT; ++t) {
        Target *target = &targets[t];
        target->X = target_sizes[t];
        target->a = (uint16_t *)calloc(nprimes, sizeof(*target->a));
        target->a_reflection =
            (uint16_t *)calloc(nprimes, sizeof(*target->a_reflection));
        target->a_central =
            (uint16_t *)calloc(nprimes, sizeof(*target->a_central));
        target->hit_count =
            (uint16_t *)calloc(target->X, sizeof(*target->hit_count));
        target->hit_offset =
            (size_t *)calloc((size_t)target->X + 1U,
                             sizeof(*target->hit_offset));
        if (target->a == NULL || target->a_reflection == NULL ||
            target->a_central == NULL || target->hit_count == NULL ||
            target->hit_offset == NULL) {
            fprintf(stderr, "fatal: target allocation failed at X=%" PRIu32
                    "\n", target->X);
            return 0;
        }
    }
    return 1;
}

static int build_hit_lists(const ZeroPair *pairs, size_t npairs,
                           const uint32_t *prime_index,
                           size_t nprimes, Target targets[TARGET_COUNT]) {
    (void)nprimes;
    for (size_t k = 0; k < npairs; ++k) {
        uint32_t p = pairs[k].p;
        uint32_t r = pairs[k].r;
        uint32_t pi = prime_index[p];
        uint64_t n = (uint64_t)p + r;
        for (size_t t = 0; t < TARGET_COUNT; ++t) {
            Target *target = &targets[t];
            uint32_t X = target->X;
            if (2ULL * p <= X || p > 2U * X || n <= X || n > 2U * X)
                continue;
            uint32_t offset = (uint32_t)(n - X - 1U);
            if (target->a[pi] == UINT16_MAX ||
                target->hit_count[offset] == UINT16_MAX) {
                fprintf(stderr, "fatal: hit counter overflow\n");
                return 0;
            }
            ++target->a[pi];
            ++target->hit_count[offset];
            ++target->S;
            if (2U * r == p - 1U)
                ++target->a_central[pi];
            else
                ++target->a_reflection[pi];
        }
    }

    for (size_t t = 0; t < TARGET_COUNT; ++t) {
        Target *target = &targets[t];
        for (uint32_t offset = 0; offset < target->X; ++offset) {
            target->hit_offset[offset + 1U] =
                target->hit_offset[offset] + target->hit_count[offset];
            uint64_t h = target->hit_count[offset];
            target->M2 += h * (h - 1U);
        }
        if (target->hit_offset[target->X] != target->S) {
            fprintf(stderr, "SANITY ABORT: hit prefix mismatch\n");
            return 0;
        }
        target->hit_primes =
            (uint32_t *)malloc((size_t)target->S * sizeof(*target->hit_primes));
        if (target->S != 0 && target->hit_primes == NULL) {
            fprintf(stderr, "fatal: hit-list allocation failed\n");
            return 0;
        }
    }

    size_t *cursor[TARGET_COUNT];
    for (size_t t = 0; t < TARGET_COUNT; ++t) {
        cursor[t] = (size_t *)malloc(
            (size_t)targets[t].X * sizeof(*cursor[t]));
        if (cursor[t] == NULL) {
            fprintf(stderr, "fatal: hit cursor allocation failed\n");
            for (size_t u = 0; u < t; ++u)
                free(cursor[u]);
            return 0;
        }
        memcpy(cursor[t], targets[t].hit_offset,
               (size_t)targets[t].X * sizeof(*cursor[t]));
    }

    for (size_t k = 0; k < npairs; ++k) {
        uint32_t p = pairs[k].p;
        uint32_t pi = prime_index[p];
        uint64_t n = (uint64_t)p + pairs[k].r;
        for (size_t t = 0; t < TARGET_COUNT; ++t) {
            Target *target = &targets[t];
            uint32_t X = target->X;
            if (2ULL * p <= X || p > 2U * X || n <= X || n > 2U * X)
                continue;
            uint32_t offset = (uint32_t)(n - X - 1U);
            target->hit_primes[cursor[t][offset]++] = pi;
        }
    }
    for (size_t t = 0; t < TARGET_COUNT; ++t)
        free(cursor[t]);
    return 1;
}

static uint32_t residue_for_hit(uint32_t n, uint32_t p) {
    return n - p;
}

static int same_spot_pair(const Spotcheck *spot, uint32_t a, uint32_t b) {
    return spot->p_index == a && spot->q_index == b;
}

static void maybe_add_spot(Target *target, uint32_t a, uint32_t b,
                           uint32_t n) {
    if (target->nspots >= SPOTCHECK_COUNT)
        return;
    for (uint32_t i = 0; i < target->nspots; ++i)
        if (same_spot_pair(&target->spots[i], a, b))
            return;
    target->spots[target->nspots++] = (Spotcheck){a, b, n, 0U, 0U};
}

static int enumerate_actual(Target *target, const uint32_t *primes) {
    for (uint32_t offset = 0; offset < target->X; ++offset) {
        size_t begin = target->hit_offset[offset];
        size_t end = target->hit_offset[offset + 1U];
        uint32_t n = target->X + 1U + offset;
        for (size_t i = begin; i < end; ++i) {
            for (size_t k = i + 1U; k < end; ++k) {
                uint32_t ai = target->hit_primes[i];
                uint32_t bi = target->hit_primes[k];
                if (ai == bi) {
                    fprintf(stderr, "SANITY ABORT: duplicate prime in H(n)\n");
                    return 0;
                }
                uint32_t a = ai < bi ? ai : bi;
                uint32_t b = ai < bi ? bi : ai;
                uint32_t p = primes[a];
                uint32_t q = primes[b];
                uint32_t r = residue_for_hit(n, p);
                uint32_t s = residue_for_hit(n, q);
                uint32_t shell = shell_index(q - p);
                if (shell >= MAX_SHELLS) {
                    fprintf(stderr, "fatal: shell index overflow\n");
                    return 0;
                }
                if (shell + 1U > target->nshells)
                    target->nshells = shell + 1U;
                target->shells[shell].actual += 2U;
                int r_noncentral = 2U * r != p - 1U;
                int s_noncentral = 2U * s != q - 1U;
                if (r_noncentral && s_noncentral)
                    target->shells[shell].reflection_actual += 2U;
                if (!r_noncentral && !s_noncentral)
                    target->shells[shell].central_actual += 2U;
                uint32_t r_fold = r < p - 1U - r ? r : p - 1U - r;
                uint32_t s_fold = s < q - 1U - s ? s : q - 1U - s;
                if (r_fold == s_fold)
                    target->shells[shell].folded_match_actual += 2U;
                if (target->X == 8192U)
                    maybe_add_spot(target, a, b, n);
            }
        }
    }
    uint64_t actual_total = 0;
    for (uint32_t j = 0; j < target->nshells; ++j)
        actual_total += target->shells[j].actual;
    if (actual_total != target->M2) {
        fprintf(stderr,
                "SANITY ABORT X=%" PRIu32 ": shell actual=%" PRIu64
                " != M2=%" PRIu64 "\n",
                target->X, actual_total, target->M2);
        return 0;
    }
    return 1;
}

static int add_expected_shells(Target *target, const uint32_t *primes,
                               size_t nprimes) {
    uint64_t *prefix =
        (uint64_t *)calloc(nprimes + 1U, sizeof(*prefix));
    uint64_t *reflection_prefix =
        (uint64_t *)calloc(nprimes + 1U, sizeof(*reflection_prefix));
    if (prefix == NULL || reflection_prefix == NULL) {
        fprintf(stderr, "fatal: expected-value prefix allocation failed\n");
        free(prefix);
        free(reflection_prefix);
        return 0;
    }
    for (size_t i = 0; i < nprimes; ++i) {
        prefix[i + 1U] = prefix[i] + target->a[i];
        reflection_prefix[i + 1U] =
            reflection_prefix[i] + target->a_reflection[i];
    }

    size_t first = lower_bound_prime(primes, nprimes,
                                     (uint64_t)target->X / 2U + 1U);
    size_t last = lower_bound_prime(primes, nprimes,
                                    2ULL * target->X + 1U);
    for (size_t i = first; i < last; ++i) {
        if (target->a[i] == 0U && target->a_reflection[i] == 0U)
            continue;
        uint64_t p = primes[i];
        for (uint32_t j = 0; j < MAX_SHELLS; ++j) {
            uint64_t low_distance = UINT64_C(1) << j;
            uint64_t high_distance = UINT64_C(1) << (j + 1U);
            size_t lo = lower_bound_prime(primes, nprimes,
                                          p + low_distance);
            size_t hi = lower_bound_prime(primes, nprimes,
                                          p + high_distance);
            if (lo < i + 1U)
                lo = i + 1U;
            if (hi > last)
                hi = last;
            if (lo >= hi)
                continue;
            uint64_t weight = prefix[hi] - prefix[lo];
            uint64_t reflection_weight =
                reflection_prefix[hi] - reflection_prefix[lo];
            target->shells[j].expected_numerator +=
                2U * (uint64_t)target->a[i] * weight;
            target->shells[j].reflection_expected_numerator +=
                2U * (uint64_t)target->a_reflection[i] * reflection_weight;
            if (j + 1U > target->nshells &&
                (weight != 0U || reflection_weight != 0U))
                target->nshells = j + 1U;
        }
    }
    free(prefix);
    free(reflection_prefix);

    uint64_t expected_total = 0;
    uint64_t reflection_expected_total = 0;
    for (uint32_t j = 0; j < target->nshells; ++j) {
        expected_total += target->shells[j].expected_numerator;
        reflection_expected_total +=
            target->shells[j].reflection_expected_numerator;
    }
    target->expected_numerator = expected_total;
    target->reflection_expected_numerator = reflection_expected_total;

    uint64_t sum_a = 0;
    uint64_t sum_a2 = 0;
    uint64_t sum_ref = 0;
    uint64_t sum_ref2 = 0;
    for (size_t i = first; i < last; ++i) {
        uint64_t a = target->a[i];
        uint64_t ar = target->a_reflection[i];
        sum_a += a;
        sum_a2 += a * a;
        sum_ref += ar;
        sum_ref2 += ar * ar;
    }
    if (sum_a != target->S ||
        expected_total != sum_a * sum_a - sum_a2 ||
        reflection_expected_total != sum_ref * sum_ref - sum_ref2) {
        fprintf(stderr,
                "SANITY ABORT X=%" PRIu32 ": benchmark prefix mismatch\n",
                target->X);
        return 0;
    }

    int64_t centered_numerator_sum = 0;
    for (uint32_t j = 0; j < target->nshells; ++j) {
        uint64_t actual_scaled = target->shells[j].actual * target->X;
        uint64_t expected = target->shells[j].expected_numerator;
        if (actual_scaled > (uint64_t)INT64_MAX ||
            expected > (uint64_t)INT64_MAX) {
            fprintf(stderr, "fatal: centered numerator overflow\n");
            return 0;
        }
        centered_numerator_sum +=
            (int64_t)actual_scaled - (int64_t)expected;
    }
    if (centered_numerator_sum + (int64_t)expected_total !=
        (int64_t)(target->M2 * target->X)) {
        fprintf(stderr,
                "SANITY ABORT X=%" PRIu32 ": centered identity failed\n",
                target->X);
        return 0;
    }
    printf("SANITY X=%" PRIu32
           " shell identity: PASS (sum D_j + E = M2)\n", target->X);
    return 1;
}

static uint32_t count_pair_from_hits(const Target *target,
                                     uint32_t p_index, uint32_t q_index) {
    uint32_t count = 0;
    for (uint32_t offset = 0; offset < target->X; ++offset) {
        int has_p = 0;
        int has_q = 0;
        for (size_t k = target->hit_offset[offset];
             k < target->hit_offset[offset + 1U]; ++k) {
            has_p |= target->hit_primes[k] == p_index;
            has_q |= target->hit_primes[k] == q_index;
        }
        if (has_p && has_q)
            ++count;
    }
    return count;
}

static uint32_t count_pair_direct(const Target *target,
                                  const ZeroPair *pairs,
                                  const size_t *zero_offset,
                                  const uint32_t *primes,
                                  uint32_t p_index, uint32_t q_index) {
    uint32_t p = primes[p_index];
    uint32_t q = primes[q_index];
    uint32_t count = 0;
    for (size_t k = zero_offset[p_index];
         k < zero_offset[p_index + 1U]; ++k) {
        uint64_t n = (uint64_t)p + pairs[k].r;
        if (n <= target->X || n > 2ULL * target->X || n < q)
            continue;
        uint32_t s = (uint32_t)n - q;
        if (s < q && zero_contains(pairs, zero_offset, q_index, s))
            ++count;
    }
    return count;
}

static int run_spotchecks(Target *target, const ZeroPair *pairs,
                          const size_t *zero_offset,
                          const uint32_t *primes) {
    if (target->nspots != SPOTCHECK_COUNT) {
        fprintf(stderr,
                "SANITY ABORT: only %" PRIu32 " spotcheck pairs found\n",
                target->nspots);
        return 0;
    }
    printf("SANITY X=%" PRIu32 " direct pair spotchecks:\n", target->X);
    for (uint32_t i = 0; i < target->nspots; ++i) {
        Spotcheck *spot = &target->spots[i];
        spot->hit_count = count_pair_from_hits(
            target, spot->p_index, spot->q_index);
        spot->direct_count = count_pair_direct(
            target, pairs, zero_offset, primes,
            spot->p_index, spot->q_index);
        printf("  p=%" PRIu32 " q=%" PRIu32 " witness_n=%" PRIu32
               " hit_count=%" PRIu32 " direct_count=%" PRIu32 " %s\n",
               primes[spot->p_index], primes[spot->q_index], spot->witness_n,
               spot->hit_count, spot->direct_count,
               spot->hit_count == spot->direct_count &&
               spot->direct_count > 0U ? "PASS" : "MISMATCH");
        if (spot->hit_count != spot->direct_count ||
            spot->direct_count == 0U)
            return 0;
    }
    return 1;
}

static int exact_m2_gates(const Target targets[TARGET_COUNT]) {
    static const uint64_t expected_m2[TARGET_COUNT] = {
        18U, 54U, 58U, 100U, 154U, 238U, 462U, 742U, 1458U
    };
    for (size_t t = 0; t < TARGET_COUNT; ++t) {
        if (targets[t].M2 != expected_m2[t]) {
            fprintf(stderr,
                    "SANITY ABORT X=%" PRIu32 ": M2=%" PRIu64
                    ", expected=%" PRIu64 "\n",
                    targets[t].X, targets[t].M2, expected_m2[t]);
            return 0;
        }
    }
    printf("SANITY M2: PASS (4000=18, 8000=54, 2^13=58, 2^19=1458; "
           "all dyadic values matched)\n");
    return 1;
}

static void print_results(const Target targets[TARGET_COUNT]) {
    printf("\n=== CENTERED PAIR DISPERSION ===\n");
    for (size_t t = FIRST_REPORT_TARGET; t < TARGET_COUNT; ++t) {
        const Target *target = &targets[t];
        long double E =
            (long double)target->expected_numerator / target->X;
        long double D = (long double)target->M2 - E;
        long double lambda = (long double)target->S / target->X;
        long double normalization = target->X * lambda * lambda;
        printf("SUMMARY X=%" PRIu32 " S=%" PRIu64 " lambda=%.12Lf"
               " M2=%" PRIu64 " E=%.12Lf D=%.12Lf Dnorm=%.12Lf\n",
               target->X, target->S, lambda, target->M2, E, D,
               normalization == 0.0L ? NAN : D / normalization);
        printf("  j low high C_j W_j E_j D_j absD_over_sqrtC "
               "R_actual R_E R_D central folded_match\n");
        for (uint32_t j = 0; j < target->nshells; ++j) {
            const Shell *shell = &target->shells[j];
            long double shell_E =
                (long double)shell->expected_numerator / target->X;
            long double shell_D = (long double)shell->actual - shell_E;
            long double z = shell->actual == 0U
                ? NAN : fabsl(shell_D) / sqrtl((long double)shell->actual);
            long double reflection_E =
                (long double)shell->reflection_expected_numerator /
                target->X;
            long double reflection_D =
                (long double)shell->reflection_actual - reflection_E;
            printf("  %" PRIu32 " %" PRIu64 " %" PRIu64
                   " %" PRIu64 " %" PRIu64 " %.12Lf %.12Lf %.12Lf"
                   " %" PRIu64 " %.12Lf %.12Lf %" PRIu64
                   " %" PRIu64 "\n",
                   j, UINT64_C(1) << j, UINT64_C(1) << (j + 1U),
                   shell->actual, shell->expected_numerator,
                   shell_E, shell_D, z,
                   shell->reflection_actual, reflection_E, reflection_D,
                   shell->central_actual, shell->folded_match_actual);
        }
    }
}

static void free_targets(Target targets[TARGET_COUNT]) {
    for (size_t t = 0; t < TARGET_COUNT; ++t) {
        free(targets[t].a);
        free(targets[t].a_reflection);
        free(targets[t].a_central);
        free(targets[t].hit_count);
        free(targets[t].hit_offset);
        free(targets[t].hit_primes);
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
    printf("CRON_pair_dispersion prime_limit=%u threads=%" PRIu32 "\n",
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

    uint32_t *prime_index =
        (uint32_t *)malloc(PRIME_LIMIT * sizeof(*prime_index));
    if (prime_index == NULL) {
        fprintf(stderr, "fatal: cannot allocate prime index\n");
        free(primes);
        return EXIT_FAILURE;
    }
    for (uint32_t p = 0; p < PRIME_LIMIT; ++p)
        prime_index[p] = UINT32_MAX;
    for (size_t i = 0; i < nprimes; ++i)
        prime_index[primes[i]] = (uint32_t)i;

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
        free(prime_index);
        free(primes);
        return EXIT_FAILURE;
    }
    SharedScan shared = {.primes = primes, .nprimes = nprimes};
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
            double eta = done == 0 ? INFINITY :
                elapsed * (expected_steps - done) / done;
            fprintf(stderr,
                    "[progress] primes=%zu/%zu steps=%" PRIu64 "/%" PRIu64
                    " (%.2f%%) elapsed=%.1fs ETA=%.1fs\n",
                    finished, nprimes, done, expected_steps,
                    expected_steps == 0 ? 100.0 :
                        100.0 * (double)done / expected_steps,
                    elapsed, eta);
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
    if (atomic_load_explicit(&shared.failed, memory_order_relaxed) ||
        created != nthreads || primes_done != nprimes ||
        steps_done != expected_steps) {
        fprintf(stderr,
                "SANITY ABORT: incomplete scan primes=%zu/%zu steps=%" PRIu64
                "/%" PRIu64 "\n",
                primes_done, nprimes, steps_done, expected_steps);
        for (uint32_t i = 0; i < nthreads; ++i)
            free(workers[i].zeros.data);
        free(workers);
        free(threads);
        free(prime_index);
        free(primes);
        return EXIT_FAILURE;
    }

    size_t npairs = 0;
    for (uint32_t i = 0; i < nthreads; ++i)
        npairs += workers[i].zeros.len;
    ZeroPair *pairs =
        npairs == 0 ? NULL : (ZeroPair *)malloc(npairs * sizeof(*pairs));
    if (npairs != 0 && pairs == NULL) {
        fprintf(stderr, "fatal: cannot merge zero pairs\n");
        return EXIT_FAILURE;
    }
    size_t merged = 0;
    for (uint32_t i = 0; i < nthreads; ++i) {
        memcpy(pairs + merged, workers[i].zeros.data,
               workers[i].zeros.len * sizeof(*pairs));
        merged += workers[i].zeros.len;
        free(workers[i].zeros.data);
    }
    qsort(pairs, npairs, sizeof(*pairs), pair_compare);
    printf("scan complete: zero_pairs=%zu scan_seconds=%.3f\n",
           npairs, scan_seconds);

    size_t *zero_offset =
        (size_t *)malloc((nprimes + 1U) * sizeof(*zero_offset));
    Target targets[TARGET_COUNT] = {0};
    int ok = zero_offset != NULL &&
             validate_zero_pairs(pairs, npairs, primes, nprimes, zero_offset) &&
             allocate_targets(targets, nprimes) &&
             build_hit_lists(pairs, npairs, prime_index, nprimes, targets);
    if (ok) {
        for (size_t t = 0; t < TARGET_COUNT && ok; ++t)
            ok = enumerate_actual(&targets[t], primes) &&
                 add_expected_shells(&targets[t], primes, nprimes);
    }
    if (ok)
        ok = exact_m2_gates(targets);
    if (ok)
        ok = run_spotchecks(&targets[2], pairs, zero_offset, primes);
    if (ok)
        print_results(targets);

    double total_seconds = monotonic_seconds() - total_start;
    if (ok)
        printf("TOTAL seconds=%.3f scan_seconds=%.3f throughput=%.6g steps/s\n",
               total_seconds, scan_seconds,
               scan_seconds == 0.0 ? 0.0 : expected_steps / scan_seconds);

    free_targets(targets);
    free(zero_offset);
    free(pairs);
    free(workers);
    free(threads);
    free(prime_index);
    free(primes);
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
