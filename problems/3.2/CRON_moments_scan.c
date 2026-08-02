/*
 * CRON_moments_scan.c -- factorial moments of the top-half Apéry bad-prime
 * count H(n).
 *
 * For every prime 7 <= p < N, compute
 *
 *   Z_p = {z : 0 <= z < p and b_z == 0 (mod p)}
 *
 * using c_z = (z!)^3 b_z and the division-free recurrence
 *
 *   c_{r+1} = (34r^3 + 51r^2 + 27r + 5)c_r - r^6 c_{r-1} (mod p).
 *
 * Each (p,z) is scattered to n=p+z, giving
 *
 *   H(n) = #{p in (n/2,n] : n-p in Z_p}.
 *
 * The hot loop advances the cubic and sixth-power coefficients by finite
 * differences and uses two 32-bit Montgomery products.  Primes are assigned
 * largest first to pthread workers.
 *
 * Compile:
 *   cc -O3 -march=native -Wall -Wextra -Wpedantic \
 *      -o /tmp/CRON_moments_scan CRON_moments_scan.c -lpthread -lm
 *
 * Usage:
 *   /tmp/CRON_moments_scan [N [START [THREADS]]]
 *
 * Defaults: N=2000000, START=10000, THREADS=online CPU count.  N and START
 * are exclusive and inclusive, respectively.  For the Python-reference gate
 * use N=30000 and START=1000.
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

#define DEFAULT_N 2000000U
#define DEFAULT_START 10000U
#define FIRST_PRIME 7U
#define TOP_COUNT 20U
#define LN2_CONST 0.693147180559945309417232121458176568L

typedef struct {
    uint32_t p;
    uint32_t z;
} Pair;

typedef struct {
    Pair *data;
    size_t len;
    size_t cap;
} PairVec;

typedef struct {
    uint32_t n;
    uint32_t h;
} TopEntry;

typedef struct {
    const uint32_t *primes;
    size_t nprimes;
    atomic_size_t next_job;
    atomic_uint_fast64_t steps_done;
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

static int parse_u32(const char *text, uint32_t *value) {
    char *end = NULL;
    errno = 0;
    unsigned long long x = strtoull(text, &end, 10);
    if (errno || end == text || !end || *end != '\0' || x > UINT32_MAX)
        return 0;
    *value = (uint32_t)x;
    return 1;
}

static int pair_push(PairVec *v, uint32_t p, uint32_t z) {
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
    v->data[v->len++] = (Pair){p, z};
    return 1;
}

static uint32_t *sieve_primes(uint32_t limit, size_t *count_out) {
    uint8_t *composite = (uint8_t *)calloc(limit, 1);
    if (!composite) {
        fprintf(stderr, "out of memory in sieve\n");
        return NULL;
    }
    for (uint32_t i = 2; (uint64_t)i * i < limit; ++i) {
        if (!composite[i]) {
            for (uint32_t j = i * i; j < limit; j += i)
                composite[j] = 1;
        }
    }
    size_t count = 0;
    for (uint32_t p = FIRST_PRIME; p < limit; ++p)
        if (!composite[p])
            ++count;
    uint32_t *primes = (uint32_t *)malloc(count * sizeof(*primes));
    if (!primes) {
        fprintf(stderr, "out of memory for prime list\n");
        free(composite);
        return NULL;
    }
    size_t k = 0;
    for (uint32_t p = FIRST_PRIME; p < limit; ++p)
        if (!composite[p])
            primes[k++] = p;
    free(composite);
    *count_out = count;
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

static inline void advance_differences(uint32_t *d, int degree,
                                       uint32_t p) {
    for (int k = 0; k < degree; ++k)
        d[k] = add_mod(d[k], d[k + 1], p);
}

/*
 * At r=1 the forward-difference rows (low order first) are
 *   34r^3+51r^2+27r+5: [117, 418, 510, 204]
 *   r^6:                 [1, 63, 602, 2100, 3360, 2520, 720].
 */
static int scan_prime(uint32_t p, PairVec *pairs) {
    static const uint32_t coeff_init[4] = {117, 418, 510, 204};
    static const uint32_t sixth_init[7] =
        {1, 63, 602, 2100, 3360, 2520, 720};
    uint32_t coeff[4], sixth[7];
    for (int k = 0; k < 4; ++k)
        coeff[k] = to_montgomery(coeff_init[k], p);
    for (int k = 0; k < 7; ++k)
        sixth[k] = to_montgomery(sixth_init[k], p);

    uint32_t nprime = montgomery_nprime(p);
    uint32_t prev = to_montgomery(1, p);
    uint32_t cur = to_montgomery(5, p);
    for (uint32_t r = 1; r < p - 1; ++r) {
        uint32_t term1 = montgomery_mul(coeff[0], cur, p, nprime);
        uint32_t term2 = montgomery_mul(sixth[0], prev, p, nprime);
        uint32_t next = sub_mod(term1, term2, p);
        if (next == 0 && !pair_push(pairs, p, r + 1))
            return 0;
        prev = cur;
        cur = next;
        advance_differences(coeff, 3, p);
        advance_differences(sixth, 6, p);
    }
    return 1;
}

/* Independent slow kernel, used only by startup tests. */
static int scan_prime_reference(uint32_t p, PairVec *pairs) {
    uint64_t prev = 1;
    uint64_t cur = 5 % p;
    for (uint32_t r = 1; r < p - 1; ++r) {
        uint64_t r2 = (uint64_t)r * r % p;
        uint64_t r3 = r2 * r % p;
        uint64_t coeff =
            (34 * r3 + 51 * r2 + 27ULL * r + 5) % p;
        uint64_t r6 = r3 * r3 % p;
        uint64_t next = (coeff * cur + p - r6 * prev % p) % p;
        if (next == 0 && !pair_push(pairs, p, r + 1))
            return 0;
        prev = cur;
        cur = next;
    }
    return 1;
}

static int startup_self_test(void) {
    static const uint32_t tests[] = {7, 11, 17, 31, 181, 379, 3137};
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i) {
        PairVec fast = {0}, reference = {0};
        int ok = scan_prime(tests[i], &fast) &&
                 scan_prime_reference(tests[i], &reference) &&
                 fast.len == reference.len &&
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
    printf("Startup kernel self-test: PASS (7 primes through 3137)\n");
    return 1;
}

static void *worker_main(void *arg) {
    Worker *worker = (Worker *)arg;
    Shared *shared = worker->shared;
    for (;;) {
        if (atomic_load_explicit(&shared->failed, memory_order_relaxed))
            break;
        size_t job = atomic_fetch_add_explicit(&shared->next_job, 1,
                                               memory_order_relaxed);
        if (job >= shared->nprimes)
            break;
        size_t index = shared->nprimes - 1 - job;
        uint32_t p = shared->primes[index];
        if (!scan_prime(p, &worker->pairs)) {
            atomic_store_explicit(&shared->failed, 1, memory_order_relaxed);
            break;
        }
        atomic_fetch_add_explicit(&shared->steps_done, p - 2,
                                  memory_order_relaxed);
    }
    return NULL;
}

static int pair_cmp(const void *va, const void *vb) {
    const Pair *a = (const Pair *)va;
    const Pair *b = (const Pair *)vb;
    if (a->p != b->p)
        return a->p < b->p ? -1 : 1;
    if (a->z != b->z)
        return a->z < b->z ? -1 : 1;
    return 0;
}

static uint32_t gcd_u32(uint32_t a, uint32_t b) {
    while (b != 0) {
        uint32_t r = a % b;
        a = b;
        b = r;
    }
    return a;
}

static int top_better(TopEntry a, TopEntry b) {
    return a.h > b.h || (a.h == b.h && a.n < b.n);
}

static void top_insert(TopEntry *top, size_t *len, TopEntry value) {
    size_t pos = 0;
    while (pos < *len && !top_better(value, top[pos]))
        ++pos;
    if (pos >= TOP_COUNT)
        return;
    size_t new_len = *len < TOP_COUNT ? *len + 1 : *len;
    for (size_t j = new_len - 1; j > pos; --j)
        top[j] = top[j - 1];
    top[pos] = value;
    *len = new_len;
}

static void print_structure_probes(const uint32_t *h, uint32_t start,
                                   uint32_t limit, const Pair *pairs,
                                   size_t npairs) {
    size_t nrecords = 0;
    printf("\n=== Structure probes for every n with H(n)>=4 ===\n");
    for (uint32_t n = start; n < limit; ++n) {
        if (h[n] < 4)
            continue;
        ++nrecords;
        printf("n=%" PRIu32 " H=%" PRIu32 "\n", n, h[n]);
        for (size_t i = 0; i < npairs; ++i) {
            uint64_t hit = (uint64_t)pairs[i].p + pairs[i].z;
            if (hit != n)
                continue;
            uint32_t g = gcd_u32(pairs[i].p - 1, n - 1);
            printf("  p=%" PRIu32 " z=%" PRIu32
                   " z/p=%.12f p_mod_24=%" PRIu32
                   " gcd(p-1,n-1)=%" PRIu32 " order=%" PRIu32 "\n",
                   pairs[i].p, pairs[i].z,
                   (double)pairs[i].z / pairs[i].p,
                   pairs[i].p % 24, g, (pairs[i].p - 1) / g);
        }
    }
    if (nrecords == 0)
        printf("none\n");
}

static int analyze(const Pair *pairs, size_t npairs, uint32_t limit,
                   uint32_t start) {
    uint32_t *h = (uint32_t *)calloc(limit, sizeof(*h));
    if (!h) {
        fprintf(stderr, "out of memory for H array\n");
        return 0;
    }
    size_t used_pairs = 0;
    for (size_t i = 0; i < npairs; ++i) {
        uint64_t n = (uint64_t)pairs[i].p + pairs[i].z;
        if (n < limit) {
            ++h[n];
            ++used_pairs;
        }
    }

    uint32_t max_h = 0;
    for (uint32_t n = start; n < limit; ++n)
        if (h[n] > max_h)
            max_h = h[n];
    uint64_t *hist =
        (uint64_t *)calloc((size_t)max_h + 1, sizeof(*hist));
    if (!hist) {
        fprintf(stderr, "out of memory for histogram\n");
        free(h);
        return 0;
    }

    long double moment_sum[7] = {0};
    long double predicted_sum[7] = {0};
    TopEntry top[TOP_COUNT];
    size_t top_len = 0;
    uint64_t count = (uint64_t)limit - start;
    for (uint32_t n = start; n < limit; ++n) {
        ++hist[h[n]];
        top_insert(top, &top_len, (TopEntry){n, h[n]});
        long double falling = 1;
        long double lambda = LN2_CONST / logl((long double)n);
        long double lambda_power = 1;
        for (int k = 1; k <= 6; ++k) {
            falling *= (long double)((int64_t)h[n] - (k - 1));
            if (h[n] < (uint32_t)k)
                falling = 0;
            lambda_power *= lambda;
            moment_sum[k] += falling;
            predicted_sum[k] += lambda_power;
        }
    }

    printf("\n=== H(n) statistics on [%" PRIu32 ",%" PRIu32 ") ===\n",
           start, limit);
    printf("sample_count=%" PRIu64 " used_zero_pairs=%zu max_H=%" PRIu32
           "\n", count, used_pairs, max_h);
    printf("Histogram:\n");
    for (uint32_t value = 0; value <= max_h; ++value)
        printf("  H=%" PRIu32 ": %" PRIu64 "\n", value, hist[value]);

    printf("Factorial moments and Poisson-mixture predictions:\n");
    printf("  k observed_Mk predicted_mean_lambda_pow_k ratio\n");
    for (int k = 1; k <= 6; ++k) {
        long double observed = moment_sum[k] / count;
        long double predicted = predicted_sum[k] / count;
        printf("  %d %.15Le %.15Le %.15Le\n", k, observed, predicted,
               predicted > 0 ? observed / predicted : 0);
    }

    printf("Top %u (ordered by H descending, then n ascending):\n",
           TOP_COUNT);
    for (size_t i = 0; i < top_len; ++i)
        printf("  rank=%zu n=%" PRIu32 " H=%" PRIu32 "\n",
               i + 1, top[i].n, top[i].h);

    print_structure_probes(h, start, limit, pairs, npairs);
    free(hist);
    free(h);
    return 1;
}

int main(int argc, char **argv) {
    uint32_t limit = DEFAULT_N;
    uint32_t start = DEFAULT_START;
    long online = sysconf(_SC_NPROCESSORS_ONLN);
    uint32_t nthreads = online > 0 ? (uint32_t)online : 1;
    if (argc > 4) {
        fprintf(stderr, "usage: %s [N [START [THREADS]]]\n", argv[0]);
        return EXIT_FAILURE;
    }
    if (argc >= 2 && !parse_u32(argv[1], &limit)) {
        fprintf(stderr, "invalid N: %s\n", argv[1]);
        return EXIT_FAILURE;
    }
    if (argc >= 3 && !parse_u32(argv[2], &start)) {
        fprintf(stderr, "invalid START: %s\n", argv[2]);
        return EXIT_FAILURE;
    }
    if (argc >= 4 && !parse_u32(argv[3], &nthreads)) {
        fprintf(stderr, "invalid THREADS: %s\n", argv[3]);
        return EXIT_FAILURE;
    }
    if (limit <= FIRST_PRIME || start < 2 || start >= limit ||
        nthreads == 0) {
        fprintf(stderr, "require N>7, 2<=START<N, and THREADS>=1\n");
        return EXIT_FAILURE;
    }

    setvbuf(stdout, NULL, _IOLBF, 0);
    printf("CRON_moments_scan N=%" PRIu32 " START=%" PRIu32
           " threads=%" PRIu32 " online_cpus=%ld\n",
           limit, start, nthreads, online);
    if (!startup_self_test())
        return EXIT_FAILURE;

    double total_t0 = monotonic_seconds();
    size_t nprimes = 0;
    uint32_t *primes = sieve_primes(limit, &nprimes);
    if (!primes)
        return EXIT_FAILURE;
    if (nthreads > nprimes)
        nthreads = (uint32_t)nprimes;
    uint64_t expected_steps = 0;
    for (size_t i = 0; i < nprimes; ++i)
        expected_steps += primes[i] - 2;
    printf("primes=%zu recurrence_steps=%" PRIu64 "\n",
           nprimes, expected_steps);

    Worker *workers = (Worker *)calloc(nthreads, sizeof(*workers));
    pthread_t *threads = (pthread_t *)malloc(nthreads * sizeof(*threads));
    if (!workers || !threads) {
        fprintf(stderr, "out of memory for worker state\n");
        free(primes);
        free(workers);
        free(threads);
        return EXIT_FAILURE;
    }
    Shared shared = {.primes = primes, .nprimes = nprimes};
    atomic_init(&shared.next_job, 0);
    atomic_init(&shared.steps_done, 0);
    atomic_init(&shared.failed, 0);

    double scan_t0 = monotonic_seconds();
    uint32_t created = 0;
    for (; created < nthreads; ++created) {
        workers[created].shared = &shared;
        int err = pthread_create(&threads[created], NULL, worker_main,
                                 &workers[created]);
        if (err != 0) {
            fprintf(stderr, "pthread_create: %s\n", strerror(err));
            atomic_store(&shared.failed, 1);
            break;
        }
    }
    for (uint32_t i = 0; i < created; ++i) {
        int err = pthread_join(threads[i], NULL);
        if (err != 0) {
            fprintf(stderr, "pthread_join: %s\n", strerror(err));
            atomic_store(&shared.failed, 1);
        }
    }
    double scan_seconds = monotonic_seconds() - scan_t0;
    uint64_t steps_done = atomic_load(&shared.steps_done);
    if (atomic_load(&shared.failed) || created != nthreads ||
        steps_done != expected_steps) {
        fprintf(stderr, "scan failed or incomplete: steps=%" PRIu64
                " expected=%" PRIu64 "\n", steps_done, expected_steps);
        for (uint32_t i = 0; i < nthreads; ++i)
            free(workers[i].pairs.data);
        free(threads);
        free(workers);
        free(primes);
        return EXIT_FAILURE;
    }

    size_t npairs = 0;
    for (uint32_t i = 0; i < nthreads; ++i) {
        if (SIZE_MAX - npairs < workers[i].pairs.len) {
            fprintf(stderr, "pair count overflow\n");
            return EXIT_FAILURE;
        }
        npairs += workers[i].pairs.len;
    }
    Pair *pairs = (Pair *)malloc(npairs * sizeof(*pairs));
    if (!pairs && npairs != 0) {
        fprintf(stderr, "out of memory merging pairs\n");
        return EXIT_FAILURE;
    }
    size_t offset = 0;
    for (uint32_t i = 0; i < nthreads; ++i) {
        memcpy(pairs + offset, workers[i].pairs.data,
               workers[i].pairs.len * sizeof(*pairs));
        offset += workers[i].pairs.len;
        free(workers[i].pairs.data);
    }
    qsort(pairs, npairs, sizeof(*pairs), pair_cmp);
    for (size_t i = 0; i < npairs; ++i) {
        if (pairs[i].p >= limit || pairs[i].z >= pairs[i].p ||
            (i > 0 && pair_cmp(&pairs[i - 1], &pairs[i]) >= 0)) {
            fprintf(stderr, "invalid or duplicate sorted pair at index %zu\n",
                    i);
            free(pairs);
            free(threads);
            free(workers);
            free(primes);
            return EXIT_FAILURE;
        }
    }
    printf("scan_seconds=%.6f throughput_steps_per_second=%.6e "
           "zero_pairs=%zu\n",
           scan_seconds, steps_done / scan_seconds, npairs);

    int ok = analyze(pairs, npairs, limit, start);
    printf("\ntotal_seconds=%.6f\n", monotonic_seconds() - total_t0);
    free(pairs);
    free(threads);
    free(workers);
    free(primes);
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
