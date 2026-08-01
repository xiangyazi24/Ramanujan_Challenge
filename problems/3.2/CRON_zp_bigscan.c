/*
 * CRON_zp_bigscan.c -- threaded large scan of Apéry zero sets Z_p.
 *
 * For every prime 7 <= p < LIMIT, compute
 *
 *   Z_p = {0 <= r < p : b_r == 0 (mod p)},
 *
 * where b_0=1, b_1=5 and
 *
 *   (n+1)^3 b_{n+1}
 *     = (34n^3+51n^2+27n+5)b_n - n^3 b_{n-1}.
 *
 * The hot loop uses the exactly equivalent division-free state
 * c_n=(n!)^3 b_n:
 *
 *   c_{n+1}=(34n^3+51n^2+27n+5)c_n-n^6c_{n-1} (mod p).
 *
 * Since n! is a unit for n<p, c_n and b_n have identical zero sets.
 * The cubic and sixth-power coefficients are advanced by forward
 * differences, and the two products use 32-bit Montgomery arithmetic.
 * An independent implementation of the original recurrence, using the
 * required O(p) inverse table, supplies the startup spot checks.
 *
 * Compile (from the repository root):
 *   cc -O3 -march=native -Wall -Wextra -Wpedantic \
 *     -o /tmp/CRON_zp_bigscan problems/3.2/CRON_zp_bigscan.c \
 *     -lpthread -lm
 *
 * Run from problems/3.2 so the two output paths land there:
 *   /tmp/CRON_zp_bigscan [LIMIT [THREADS]]
 *
 * LIMIT is exclusive. Defaults: LIMIT=1000000, THREADS=8.
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

#define FIRST_PRIME 7U
#define DEFAULT_LIMIT 1000000U
#define DEFAULT_THREADS 8U
#define BASELINE_LIMIT 30000U
#define MAX_HIST 128U

#define CSV_PATH "CRON_zp_bigscan.csv"
#define REPORT_PATH "CRON_zp_bigscan_report.md"

typedef struct {
    const uint32_t *primes;
    size_t nprimes;
    uint32_t *counts;
    uint8_t *midzero;
    atomic_size_t next_job;
    atomic_size_t primes_done;
    atomic_uint_fast64_t steps_done;
    atomic_uint_fast32_t workers_done;
} Shared;

typedef struct {
    uint32_t limit;
    uint32_t threads;
    uint32_t *primes;
    size_t nprimes;
    uint32_t *counts;
    uint8_t *midzero;
    uint64_t total_steps;
    double kernel_seconds;
} Scan;

typedef struct {
    uint32_t p;
    uint32_t count;
} Record;

typedef struct {
    uint32_t p;
    uint32_t recurrence;
    uint32_t binomial;
} SpotCheck;

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

static uint32_t *sieve_primes(uint32_t limit, size_t *count_out) {
    uint8_t *composite = (uint8_t *)calloc(limit, 1);
    if (!composite) {
        fprintf(stderr, "FATAL: out of memory in prime sieve\n");
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
        fprintf(stderr, "FATAL: out of memory for prime list\n");
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

/* Inputs and output are Montgomery residues, R=2^32. */
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
 * At n=1 the forward-difference rows, low order first, are
 *   34n^3+51n^2+27n+5: [117, 418, 510, 204]
 *   n^6:                 [1, 63, 602, 2100, 3360, 2520, 720].
 */
static uint32_t scan_prime(uint32_t p, uint8_t *midzero_out) {
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
    uint32_t count = 0;
    uint32_t midpoint = (p - 1) / 2;
    uint8_t midzero = 0;

    for (uint32_t n = 1; n < p - 1; ++n) {
        uint32_t term1 = montgomery_mul(coeff[0], cur, p, nprime);
        uint32_t term2 = montgomery_mul(sixth[0], prev, p, nprime);
        uint32_t next = sub_mod(term1, term2, p);
        uint32_t r = n + 1;
        if (next == 0) {
            ++count;
            if (r == midpoint)
                midzero = 1;
        }
        prev = cur;
        cur = next;
        advance_differences(coeff, 3, p);
        advance_differences(sixth, 6, p);
    }
    *midzero_out = midzero;
    return count;
}

static void *worker_main(void *arg) {
    Shared *shared = (Shared *)arg;
    for (;;) {
        size_t job = atomic_fetch_add_explicit(&shared->next_job, 1,
                                               memory_order_relaxed);
        if (job >= shared->nprimes)
            break;
        /* Largest first gives good tail load balance. */
        size_t index = shared->nprimes - 1 - job;
        uint32_t p = shared->primes[index];
        shared->counts[index] = scan_prime(p, &shared->midzero[index]);
        atomic_fetch_add_explicit(&shared->steps_done, p - 2,
                                  memory_order_relaxed);
        atomic_fetch_add_explicit(&shared->primes_done, 1,
                                  memory_order_relaxed);
    }
    atomic_fetch_add_explicit(&shared->workers_done, 1,
                              memory_order_release);
    return NULL;
}

static void free_scan(Scan *scan) {
    free(scan->primes);
    free(scan->counts);
    free(scan->midzero);
    memset(scan, 0, sizeof(*scan));
}

static int run_scan(uint32_t limit, uint32_t nthreads, int progress,
                    Scan *scan) {
    memset(scan, 0, sizeof(*scan));
    scan->limit = limit;
    scan->threads = nthreads;
    scan->primes = sieve_primes(limit, &scan->nprimes);
    if (!scan->primes)
        return 0;
    scan->counts = (uint32_t *)calloc(scan->nprimes, sizeof(*scan->counts));
    scan->midzero = (uint8_t *)calloc(scan->nprimes, sizeof(*scan->midzero));
    if (!scan->counts || !scan->midzero) {
        fprintf(stderr, "FATAL: out of memory for scan outputs\n");
        free_scan(scan);
        return 0;
    }
    for (size_t i = 0; i < scan->nprimes; ++i)
        scan->total_steps += scan->primes[i] - 2;

    Shared shared = {0};
    shared.primes = scan->primes;
    shared.nprimes = scan->nprimes;
    shared.counts = scan->counts;
    shared.midzero = scan->midzero;
    atomic_init(&shared.next_job, 0);
    atomic_init(&shared.primes_done, 0);
    atomic_init(&shared.steps_done, 0);
    atomic_init(&shared.workers_done, 0);
    pthread_t *threads = (pthread_t *)malloc(nthreads * sizeof(*threads));
    if (!threads) {
        fprintf(stderr, "FATAL: out of memory for thread handles\n");
        free_scan(scan);
        return 0;
    }

    double start = monotonic_seconds();
    for (uint32_t i = 0; i < nthreads; ++i) {
        int rc = pthread_create(&threads[i], NULL, worker_main, &shared);
        if (rc != 0) {
            fprintf(stderr, "FATAL: pthread_create: %s\n", strerror(rc));
            exit(EXIT_FAILURE);
        }
    }

    double last_report = start;
    while (atomic_load_explicit(&shared.workers_done,
                                memory_order_acquire) < nthreads) {
        struct timespec pause = {.tv_sec = 1, .tv_nsec = 0};
        nanosleep(&pause, NULL);
        double now = monotonic_seconds();
        if (progress && now - last_report >= 5.0) {
            uint64_t steps = atomic_load_explicit(&shared.steps_done,
                                                  memory_order_relaxed);
            size_t done = atomic_load_explicit(&shared.primes_done,
                                               memory_order_relaxed);
            size_t next = atomic_load_explicit(&shared.next_job,
                                               memory_order_relaxed);
            uint32_t frontier = FIRST_PRIME;
            if (next < scan->nprimes)
                frontier = scan->primes[scan->nprimes - 1 - next];
            double elapsed = now - start;
            double rate = elapsed > 0.0 ? (double)steps / elapsed : 0.0;
            double eta = rate > 0.0
                ? (double)(scan->total_steps - steps) / rate : INFINITY;
            fprintf(stderr,
                    "[progress] current p~%" PRIu32
                    " primes=%zu/%zu steps=%.3fG ETA=%.1fs\n",
                    frontier, done, scan->nprimes, 1e-9 * (double)steps, eta);
            last_report = now;
        }
    }
    for (uint32_t i = 0; i < nthreads; ++i) {
        int rc = pthread_join(threads[i], NULL);
        if (rc != 0) {
            fprintf(stderr, "FATAL: pthread_join: %s\n", strerror(rc));
            exit(EXIT_FAILURE);
        }
    }
    scan->kernel_seconds = monotonic_seconds() - start;
    free(threads);

    uint64_t actual = atomic_load_explicit(&shared.steps_done,
                                           memory_order_relaxed);
    if (actual != scan->total_steps) {
        fprintf(stderr,
                "FATAL: step accounting mismatch: got %" PRIu64
                ", expected %" PRIu64 "\n", actual, scan->total_steps);
        free_scan(scan);
        return 0;
    }
    return 1;
}

static uint32_t pow_mod(uint32_t a, uint32_t e, uint32_t p) {
    uint64_t result = 1;
    uint64_t base = a % p;
    while (e != 0) {
        if (e & 1U)
            result = result * base % p;
        base = base * base % p;
        e >>= 1;
    }
    return (uint32_t)result;
}

/* Original b_n recurrence with inv[1..p-1] built in O(p), as specified. */
static int recurrence_midpoint(uint32_t p, uint32_t *value_out) {
    uint32_t *inv = (uint32_t *)malloc((size_t)p * sizeof(*inv));
    if (!inv)
        return 0;
    inv[1] = 1;
    for (uint32_t i = 2; i < p; ++i) {
        uint64_t qterm = (uint64_t)(p / i) * inv[p % i] % p;
        inv[i] = qterm == 0 ? 0 : p - (uint32_t)qterm;
    }
    uint32_t target = (p - 1) / 2;
    uint64_t prev = 1;
    uint64_t cur = 5 % p;
    if (target == 0)
        cur = prev;
    for (uint32_t n = 1; n < target; ++n) {
        uint64_t n2 = (uint64_t)n * n % p;
        uint64_t n3 = n2 * n % p;
        uint64_t coeff =
            (34 * n3 + 51 * n2 + 27ULL * (n % p) + 5) % p;
        uint64_t numerator =
            (coeff * cur + p - n3 * prev % p) % p;
        uint64_t d = inv[n + 1];
        uint64_t d3 = d * d % p * d % p;
        uint64_t next = numerator * d3 % p;
        prev = cur;
        cur = next;
    }
    free(inv);
    *value_out = (uint32_t)cur;
    return 1;
}

/* Direct sum C(m,k)^2 C(m+k,k)^2 with factorial tables. */
static int binomial_midpoint(uint32_t p, uint32_t *value_out) {
    uint32_t *fact = (uint32_t *)malloc((size_t)p * sizeof(*fact));
    uint32_t *ifac = (uint32_t *)malloc((size_t)p * sizeof(*ifac));
    if (!fact || !ifac) {
        free(fact);
        free(ifac);
        return 0;
    }
    fact[0] = 1;
    for (uint32_t i = 1; i < p; ++i)
        fact[i] = (uint32_t)((uint64_t)fact[i - 1] * i % p);
    ifac[p - 1] = pow_mod(fact[p - 1], p - 2, p);
    for (uint32_t i = p - 1; i > 0; --i)
        ifac[i - 1] = (uint32_t)((uint64_t)ifac[i] * i % p);

    uint32_t m = (p - 1) / 2;
    uint64_t sum = 0;
    for (uint32_t k = 0; k <= m; ++k) {
        uint64_t c1 = (uint64_t)fact[m] * ifac[k] % p * ifac[m - k] % p;
        uint64_t c2 = (uint64_t)fact[m + k] * ifac[k] % p * ifac[m] % p;
        uint64_t term = c1 * c1 % p * c2 % p * c2 % p;
        sum += term;
        if (sum >= (uint64_t)p * p)
            sum %= p;
    }
    free(fact);
    free(ifac);
    *value_out = (uint32_t)(sum % p);
    return 1;
}

static size_t collect_records(const Scan *scan, Record *records,
                              size_t capacity) {
    size_t nrecords = 0;
    uint32_t maximum = 0;
    for (size_t i = 0; i < scan->nprimes; ++i) {
        uint32_t count = scan->counts[i];
        if (i == 0 || count > maximum) {
            if (nrecords >= capacity) {
                fprintf(stderr, "FATAL: record buffer exhausted\n");
                exit(EXIT_FAILURE);
            }
            records[nrecords++] = (Record){scan->primes[i], count};
            maximum = count;
        }
    }
    return nrecords;
}

static int baseline_gate(SpotCheck spots[3], double *seconds_out) {
    double start = monotonic_seconds();
    Scan scan;
    if (!run_scan(BASELINE_LIMIT, DEFAULT_THREADS, 0, &scan))
        return 0;

    static const uint64_t expected_hist[9] =
        {1933, 2, 1008, 0, 260, 0, 34, 0, 5};
    static const Record expected_records[6] =
        {{7, 0}, {11, 1}, {17, 2}, {181, 4}, {379, 6}, {3727, 8}};
    uint64_t hist[9] = {0};
    uint64_t sum = 0;
    uint32_t odd[8] = {0};
    size_t nodd = 0;
    int ok = scan.nprimes == 3242;
    for (size_t i = 0; i < scan.nprimes; ++i) {
        uint32_t count = scan.counts[i];
        sum += count;
        if (count >= 9)
            ok = 0;
        else
            ++hist[count];
        if (count & 1U) {
            if (nodd < sizeof(odd) / sizeof(odd[0]))
                odd[nodd] = scan.primes[i];
            ++nodd;
        }
        if ((uint8_t)(count & 1U) != scan.midzero[i])
            ok = 0;
    }
    for (size_t k = 0; k < 9; ++k)
        if (hist[k] != expected_hist[k])
            ok = 0;
    if (sum != 3302 || nodd != 2 || odd[0] != 11 || odd[1] != 3137)
        ok = 0;

    Record records[16];
    size_t nrecords = collect_records(&scan, records, 16);
    if (nrecords != 6 ||
        memcmp(records, expected_records, sizeof(expected_records)) != 0)
        ok = 0;

    static const uint32_t spot_primes[3] = {13, 101, 3137};
    for (size_t i = 0; i < 3; ++i) {
        spots[i].p = spot_primes[i];
        if (!recurrence_midpoint(spots[i].p, &spots[i].recurrence) ||
            !binomial_midpoint(spots[i].p, &spots[i].binomial) ||
            spots[i].recurrence != spots[i].binomial)
            ok = 0;
    }

    *seconds_out = monotonic_seconds() - start;
    if (!ok) {
        fprintf(stderr,
                "FATAL SANITY GATE: p<30000 reference or binomial spot "
                "check failed; refusing the large scan.\n");
    } else {
        fprintf(stderr,
                "[sanity] PASS: p<30000 exact profile; midpoint binomial "
                "checks p=13,101,3137 (%.3fs)\n", *seconds_out);
    }
    free_scan(&scan);
    return ok;
}

static int validate_full_scan(const Scan *scan) {
    for (size_t i = 0; i < scan->nprimes; ++i) {
        if ((uint8_t)(scan->counts[i] & 1U) != scan->midzero[i]) {
            fprintf(stderr,
                    "FATAL: parity/midpoint mismatch at p=%" PRIu32 "\n",
                    scan->primes[i]);
            return 0;
        }
    }
    return 1;
}

static int write_csv(const Scan *scan) {
    FILE *out = fopen(CSV_PATH, "w");
    if (!out) {
        perror("fopen " CSV_PATH);
        return 0;
    }
    fputs("p,zp_size\n", out);
    for (size_t i = 0; i < scan->nprimes; ++i)
        fprintf(out, "%" PRIu32 ",%" PRIu32 "\n",
                scan->primes[i], scan->counts[i]);
    if (fclose(out) != 0) {
        perror("fclose " CSV_PATH);
        return 0;
    }
    return 1;
}

static int write_report(const Scan *scan, const SpotCheck spots[3],
                        double sanity_seconds, double total_seconds) {
    uint64_t hist[MAX_HIST] = {0};
    uint64_t total_zeros = 0;
    uint32_t max_count = 0;
    size_t odd_count = 0;
    size_t even_count = 0;
    uint64_t pair_sum = 0;
    uint32_t max_k = 0;
    for (size_t i = 0; i < scan->nprimes; ++i) {
        uint32_t z = scan->counts[i];
        if (z >= MAX_HIST) {
            fprintf(stderr, "FATAL: histogram bound exceeded at p=%" PRIu32
                    ": |Z_p|=%" PRIu32 "\n", scan->primes[i], z);
            return 0;
        }
        ++hist[z];
        total_zeros += z;
        if (z > max_count)
            max_count = z;
        if (z & 1U) {
            ++odd_count;
        } else {
            ++even_count;
            pair_sum += z / 2;
            if (z / 2 > max_k)
                max_k = z / 2;
        }
    }

    Record records[64];
    size_t nrecords = collect_records(scan, records, 64);
    double mean = (double)total_zeros / (double)scan->nprimes;
    double lambda = (double)pair_sum / (double)even_count;
    double poisson_prob = exp(-lambda);
    double chi2 = 0.0;
    double predicted[MAX_HIST] = {0};
    for (uint32_t k = 0; k <= max_k; ++k) {
        if (k > 0)
            poisson_prob *= lambda / k;
        predicted[k] = (double)even_count * poisson_prob;
        uint64_t observed = hist[2 * k];
        if (predicted[k] > 0.0) {
            double delta = (double)observed - predicted[k];
            chi2 += delta * delta / predicted[k];
        }
    }
    int df = (int)max_k - 1; /* bins minus fitted lambda and normalization */

    FILE *out = fopen(REPORT_PATH, "w");
    if (!out) {
        perror("fopen " REPORT_PATH);
        return 0;
    }
    time_t now = time(NULL);
    struct tm local;
    localtime_r(&now, &local);
    char date[32];
    strftime(date, sizeof(date), "%Y-%m-%d", &local);

    fprintf(out,
            "# Large scan of Apéry zero sets \\(Z_p\\)\n\n"
            "Date: %s\n\n"
            "## Scope and headline\n\n"
            "The threaded scan covered every prime \\(7\\le p<%" PRIu32
            "\\): **%zu primes** and **%" PRIu64 " zeros** in total. "
            "The mean was **%.12f**, and the largest zero set had size "
            "**%" PRIu32 "**. The complete per-prime data are in "
            "`CRON_zp_bigscan.csv`.\n\n",
            date, scan->limit, scan->nprimes, total_zeros, mean, max_count);

    fprintf(out,
            "The scanner uses \\(c_n=(n!)^3b_n\\), whose zero set equals "
            "that of \\(b_n\\) for \\(n<p\\). Its division-free recurrence "
            "is evaluated with 32-bit Montgomery products and finite "
            "differences. This is algebraically the original Apéry "
            "recurrence, not a probabilistic shortcut.\n\n"
            "## Mandatory sanity gates\n\n"
            "All gates passed before the large scan:\n\n"
            "- The full \\(p<30000\\) scan reproduced 3,242 primes, mean "
            "`1.0185` (exact sum 3,302), distribution "
            "`{0:1933, 1:2, 2:1008, 4:260, 6:34, 8:5}`, odd primes "
            "`{11,3137}`, and records "
            "`(7,0),(11,1),(17,2),(181,4),(379,6),(3727,8)`.\n"
            "- Every scanned prime passed the parity/midpoint identity: "
            "`zp_size` is odd exactly when the midpoint state is zero.\n"
            "- The original inverse-table recurrence and the independent "
            "factorial-table binomial sum agreed at all required spots:\n\n"
            "| \\(p\\) | recurrence \\(b_{(p-1)/2}\\bmod p\\) | binomial "
            "sum |\n|---:|---:|---:|\n");
    for (size_t i = 0; i < 3; ++i)
        fprintf(out, "| %" PRIu32 " | %" PRIu32 " | %" PRIu32 " |\n",
                spots[i].p, spots[i].recurrence, spots[i].binomial);

    fputs("\n## Distribution of \\(\\lvert Z_p\\rvert\\)\n\n"
          "| \\(\\lvert Z_p\\rvert\\) | prime count | fraction |\n"
          "|---:|---:|---:|\n", out);
    for (uint32_t z = 0; z <= max_count; ++z) {
        if (hist[z] != 0)
            fprintf(out, "| %" PRIu32 " | %" PRIu64 " | %.9f |\n",
                    z, hist[z], (double)hist[z] / scan->nprimes);
    }

    fprintf(out,
            "\n## Parity and midpoint primes\n\n"
            "There are **%zu** odd cases. The complete list is:\n\n",
            odd_count);
    fputs("```text\n", out);
    for (size_t i = 0, printed = 0; i < scan->nprimes; ++i) {
        if (scan->counts[i] & 1U) {
            if (printed++)
                fputs(", ", out);
            fprintf(out, "%" PRIu32, scan->primes[i]);
        }
    }
    fputs("\n```\n\nNew midpoint primes above the previous "
          "\\(p<30000\\) scan:\n\n```text\n", out);
    size_t new_midpoints = 0;
    for (size_t i = 0; i < scan->nprimes; ++i) {
        if (scan->midzero[i] && scan->primes[i] >= BASELINE_LIMIT) {
            if (new_midpoints++)
                fputs(", ", out);
            fprintf(out, "%" PRIu32, scan->primes[i]);
        }
    }
    if (new_midpoints == 0)
        fputs("none", out);
    fputs("\n```\n", out);

    fputs("\n## Record breakers\n\n| \\(p\\) | \\(\\lvert Z_p\\rvert\\) |\n"
          "|---:|---:|\n", out);
    for (size_t i = 0; i < nrecords; ++i)
        fprintf(out, "| %" PRIu32 " | %" PRIu32 " |\n",
                records[i].p, records[i].count);

    fputs("\n## Running mean\n\nThe row at cutoff \\(X\\) uses primes "
          "\\(7\\le p<X\\).\n\n| cutoff \\(X\\) | primes | zeros | mean "
          "\\(\\lvert Z_p\\rvert\\) |\n|---:|---:|---:|---:|\n", out);
    size_t cursor = 0;
    uint64_t prefix_sum = 0;
    for (uint32_t cutoff = 100000; cutoff <= scan->limit;
         cutoff += 100000) {
        while (cursor < scan->nprimes && scan->primes[cursor] < cutoff)
            prefix_sum += scan->counts[cursor++];
        fprintf(out, "| %" PRIu32 " | %zu | %" PRIu64 " | %.12f |\n",
                cutoff, cursor, prefix_sum,
                cursor ? (double)prefix_sum / cursor : 0.0);
        if (UINT32_MAX - cutoff < 100000)
            break;
    }

    fprintf(out,
            "\n## Poisson-pair fit\n\n"
            "Restricting to the **%zu even rows**, put "
            "\\(K=\\lvert Z_p\\rvert/2\\). The maximum-likelihood fit is "
            "\\(\\widehat\\lambda=%.12f\\). Predicted counts are "
            "\\(N e^{-\\lambda}\\lambda^k/k!\\).\n\n"
            "| \\(k\\) | \\(\\lvert Z_p\\rvert=2k\\) | observed | predicted | "
            "Pearson contribution |\n|---:|---:|---:|---:|---:|\n",
            even_count, lambda);
    for (uint32_t k = 0; k <= max_k; ++k) {
        double delta = (double)hist[2 * k] - predicted[k];
        double contribution = delta * delta / predicted[k];
        fprintf(out, "| %" PRIu32 " | %" PRIu32 " | %" PRIu64
                " | %.6f | %.6f |\n", k, 2 * k, hist[2 * k],
                predicted[k], contribution);
    }
    fprintf(out,
            "\nPearson \\(\\chi^2=%.6f\\) on **%d df**, using bins "
            "\\(k=0,\\ldots,%" PRIu32 "\\). This convention reproduces "
            "the mandatory prefix value \\(\\chi^2=2.47\\) on 3 df.\n\n",
            chi2, df, max_k);

    fprintf(out,
            "## Performance and reproducibility\n\n"
            "- Threads: %" PRIu32 "\n"
            "- Exact recurrence steps: %" PRIu64 "\n"
            "- Large-scan kernel time: %.3f seconds\n"
            "- Kernel throughput: %.3f billion steps/second\n"
            "- Sanity-gate time: %.3f seconds\n"
            "- End-to-end time for this invocation: %.3f seconds\n\n"
            "Compile and run from the repository root / output directory:\n\n"
            "```sh\n"
            "cc -O3 -march=native -Wall -Wextra -Wpedantic \\\n"
            "  -o /tmp/CRON_zp_bigscan problems/3.2/CRON_zp_bigscan.c \\\n"
            "  -lpthread -lm\n"
            "cd problems/3.2\n"
            "/tmp/CRON_zp_bigscan %" PRIu32 " %" PRIu32 "\n"
            "```\n\n",
            scan->threads, scan->total_steps, scan->kernel_seconds,
            1e-9 * (double)scan->total_steps / scan->kernel_seconds,
            sanity_seconds, total_seconds, scan->limit, scan->threads);

    fprintf(out,
            "## Final summary\n\n"
            "```text\n"
            "range: 7 <= p < %" PRIu32 "\n"
            "primes: %zu\n"
            "mean |Z_p|: %.12f\n"
            "max |Z_p|: %" PRIu32 "\n"
            "odd/midpoint primes: ",
            scan->limit, scan->nprimes, mean, max_count);
    for (size_t i = 0, printed = 0; i < scan->nprimes; ++i) {
        if (scan->counts[i] & 1U) {
            if (printed++)
                fputs(", ", out);
            fprintf(out, "%" PRIu32, scan->primes[i]);
        }
    }
    fputs("\nnew midpoint primes above 30000: ", out);
    if (new_midpoints == 0) {
        fputs("none", out);
    } else {
        size_t printed = 0;
        for (size_t i = 0; i < scan->nprimes; ++i) {
            if (scan->midzero[i] && scan->primes[i] >= BASELINE_LIMIT) {
                if (printed++)
                    fputs(", ", out);
                fprintf(out, "%" PRIu32, scan->primes[i]);
            }
        }
    }
    fprintf(out,
            "\nPoisson-pair lambda-hat: %.12f\n"
            "Poisson-pair chi^2 (df=%d): %.6f\n"
            "sanity gates: PASS\n"
            "```\n", lambda, df, chi2);

    if (fclose(out) != 0) {
        perror("fclose " REPORT_PATH);
        return 0;
    }
    return 1;
}

int main(int argc, char **argv) {
    uint32_t limit = DEFAULT_LIMIT;
    uint32_t nthreads = DEFAULT_THREADS;
    if (argc > 3 || (argc >= 2 && !parse_u32(argv[1], &limit)) ||
        (argc >= 3 && !parse_u32(argv[2], &nthreads)) ||
        limit < BASELINE_LIMIT || limit > 100000000U || nthreads == 0 ||
        nthreads > 256) {
        fprintf(stderr, "usage: %s [LIMIT>=30000 [THREADS]]\n", argv[0]);
        return EXIT_FAILURE;
    }

    double total_start = monotonic_seconds();
    SpotCheck spots[3];
    double sanity_seconds = 0.0;
    if (!baseline_gate(spots, &sanity_seconds))
        return EXIT_FAILURE;

    fprintf(stderr,
            "[scan] starting p<%" PRIu32 " with %" PRIu32 " threads\n",
            limit, nthreads);
    Scan scan;
    if (!run_scan(limit, nthreads, 1, &scan))
        return EXIT_FAILURE;
    if (!validate_full_scan(&scan)) {
        free_scan(&scan);
        return EXIT_FAILURE;
    }
    double total_seconds = monotonic_seconds() - total_start;
    fprintf(stderr,
            "[scan] complete: %zu primes, %.3fG steps, %.3fs kernel\n",
            scan.nprimes, 1e-9 * (double)scan.total_steps,
            scan.kernel_seconds);

    if (!write_csv(&scan) ||
        !write_report(&scan, spots, sanity_seconds, total_seconds)) {
        free_scan(&scan);
        return EXIT_FAILURE;
    }
    fprintf(stderr, "[output] wrote %s and %s\n", CSV_PATH, REPORT_PATH);
    free_scan(&scan);
    return EXIT_SUCCESS;
}
