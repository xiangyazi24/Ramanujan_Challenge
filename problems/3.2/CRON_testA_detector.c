/*
 * Test A: exact zero-detector complexity for Apéry rows modulo primes.
 *
 * Build: cc -O2 CRON_testA_detector.c -o CRON_testA_detector
 * Run:   ./CRON_testA_detector > results.csv
 *
 * Exact checks use integer arithmetic only.  Floating point is used solely
 * for the two displayed ratios after every exact check has passed.
 */

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { PRIME_MIN = 500, PRIME_MAX = 4000, PROGRESS_STEP = 50 };

typedef struct {
    int p;
    int zero_count;
    int image_size;
    uint64_t energy;
    int detector_degree;
    int zero_function;
} PrimeStats;

static void fail(const char *message, int p)
{
    fprintf(stderr, "FAIL p=%d: %s\n", p, message);
    exit(EXIT_FAILURE);
}

static int is_prime(int n)
{
    int d;

    if (n < 2) return 0;
    if ((n & 1) == 0) return n == 2;
    for (d = 3; (int64_t)d * d <= n; d += 2) {
        if (n % d == 0) return 0;
    }
    return 1;
}

static uint32_t mul_mod(uint32_t a, uint32_t b, uint32_t p)
{
    return (uint32_t)(((uint64_t)a * b) % p);
}

static uint32_t pow_mod(uint32_t base, uint32_t exponent, uint32_t p)
{
    uint32_t result = 1;

    while (exponent != 0) {
        if (exponent & 1U) result = mul_mod(result, base, p);
        base = mul_mod(base, base, p);
        exponent >>= 1;
    }
    return result;
}

/* Fermat inverse; every caller has already established a != 0 mod p. */
static uint32_t inverse_mod(uint32_t a, uint32_t p)
{
    if (a == 0) fail("attempted to invert zero", (int)p);
    return pow_mod(a, p - 2, p);
}

static void compute_apery_row(int p, uint32_t *row)
{
    int n;

    row[0] = 1U;
    row[1] = 5U % (uint32_t)p;
    for (n = 1; n <= p - 3; ++n) {
        uint32_t nn = (uint32_t)n;
        uint32_t n2 = mul_mod(nn, nn, (uint32_t)p);
        uint32_t n3 = mul_mod(n2, nn, (uint32_t)p);
        uint32_t coefficient =
            (uint32_t)((34ULL * n3 + 51ULL * n2 + 27ULL * nn + 5ULL) %
                       (uint32_t)p);
        uint32_t next = (uint32_t)(n + 1);
        uint32_t next2 = mul_mod(next, next, (uint32_t)p);
        uint32_t denominator = mul_mod(next2, next, (uint32_t)p);
        uint32_t positive = mul_mod(coefficient, row[n], (uint32_t)p);
        uint32_t negative = mul_mod(n3, row[n - 1], (uint32_t)p);
        uint32_t numerator =
            positive >= negative ? positive - negative
                                 : positive + (uint32_t)p - negative;

        row[n + 1] = mul_mod(numerator,
                             inverse_mod(denominator, (uint32_t)p),
                             (uint32_t)p);
    }
}

/* Horner evaluation of coeff[0] + ... + coeff[degree] x^degree. */
static uint32_t evaluate_polynomial(const uint32_t *coeff, int degree,
                                    uint32_t x, uint32_t p)
{
    uint32_t value = coeff[degree];
    int j;

    for (j = degree; j-- > 0;) {
        value = (uint32_t)(((uint64_t)value * x + coeff[j]) % p);
    }
    return value;
}

/*
 * Build the Lagrange basis polynomial at the value-set node 0:
 * Q(x) = product_{v in image, v != 0} (1 - x/v).
 * The coefficient array is updated in place, from high degree to low.
 */
static int build_nonzero_root_detector(int p, const uint32_t *multiplicity,
                                       uint32_t *coeff)
{
    int degree = 0;
    int v;

    memset(coeff, 0, (size_t)(p + 1) * sizeof(*coeff));
    coeff[0] = 1U;
    for (v = 1; v < p; ++v) {
        int j;
        uint32_t inv;
        uint32_t alpha;

        if (multiplicity[v] == 0) continue;
        inv = inverse_mod((uint32_t)v, (uint32_t)p);
        alpha = (uint32_t)p - inv; /* -1/v mod p; nonzero. */
        coeff[degree + 1] = mul_mod(alpha, coeff[degree], (uint32_t)p);
        for (j = degree; j >= 1; --j) {
            coeff[j] = (uint32_t)
                (((uint64_t)coeff[j] + mul_mod(alpha, coeff[j - 1],
                                                (uint32_t)p)) %
                 (uint32_t)p);
        }
        ++degree;
    }
    return degree;
}

static PrimeStats analyze_prime(int p, int print_zero_set)
{
    uint32_t *row = calloc((size_t)(p - 1), sizeof(*row));
    uint32_t *multiplicity = calloc((size_t)p, sizeof(*multiplicity));
    uint32_t *coeff = calloc((size_t)(p + 1), sizeof(*coeff));
    uint32_t *detector_value = calloc((size_t)p, sizeof(*detector_value));
    PrimeStats stats = { p, 0, 0, 0, 0, 0 };
    int polynomial_degree;
    int r;
    int v;

    if (row == NULL || multiplicity == NULL || coeff == NULL ||
        detector_value == NULL) {
        fail("allocation failed", p);
    }

    compute_apery_row(p, row);
    if (p == 13 || p == 17 || p == 29) {
        for (r = 0; r <= p - 2; ++r) {
            int expected_zero = p == 17 && (r == 3 || r == 13);
            if ((row[r] == 0) != expected_zero)
                fail("spot-check zero set disagrees with ground truth", p);
        }
    }
    /* This exact pairing explains the effective half-row random model. */
    for (r = 1; r <= p - 2; ++r) {
        if (row[r] != row[p - 1 - r])
            fail("row pairing T(r)=T(p-1-r) failed", p);
    }
    for (r = 0; r <= p - 2; ++r) {
        ++multiplicity[row[r]];
        if (row[r] == 0) ++stats.zero_count;
    }
    for (v = 0; v < p; ++v) {
        if (multiplicity[v] != 0) ++stats.image_size;
        stats.energy += (uint64_t)multiplicity[v] * multiplicity[v];
    }

    /* Exact Cauchy--Schwarz: |S| E >= (sum_v N(v))^2 = (p-1)^2. */
    if ((uint64_t)stats.image_size * stats.energy <
        (uint64_t)(p - 1) * (uint64_t)(p - 1)) {
        fail("exact Cauchy--Schwarz check failed", p);
    }

    if (multiplicity[0] != 0) {
        polynomial_degree =
            build_nonzero_root_detector(p, multiplicity, coeff);
        stats.detector_degree = stats.image_size - 1;
        if (polynomial_degree != stats.detector_degree)
            fail("detector degree disagrees with |S|-1", p);
        if (coeff[0] != 1U || coeff[polynomial_degree] == 0U)
            fail("detector normalization or leading coefficient failed", p);
    } else {
        /*
         * The indicator is identically zero on the image.  Per the test
         * specification, record the zero function with degree convention
         * |S|, while its coefficient representation is the zero polynomial.
         */
        polynomial_degree = 0;
        stats.detector_degree = stats.image_size;
        stats.zero_function = 1;
        memset(coeff, 0, (size_t)(p + 1) * sizeof(*coeff));
    }

    /* Evaluate at every distinct image value, then check every row entry. */
    for (v = 0; v < p; ++v) {
        uint32_t expected;
        if (multiplicity[v] == 0) continue;
        detector_value[v] = evaluate_polynomial(coeff, polynomial_degree,
                                                 (uint32_t)v, (uint32_t)p);
        expected = (v == 0) ? 1U : 0U;
        if (detector_value[v] != expected)
            fail("detector failed on an image value", p);
    }
    for (r = 0; r <= p - 2; ++r) {
        uint32_t expected = row[r] == 0 ? 1U : 0U;
        if (detector_value[row[r]] != expected)
            fail("detector identity failed on a row entry", p);
    }

    if (print_zero_set) {
        int first = 1;
        printf("# spot p=%d Z={", p);
        for (r = 0; r <= p - 2; ++r) {
            if (row[r] != 0) continue;
            if (!first) putchar(',');
            printf("%d", r);
            first = 0;
        }
        printf("}\n");
    }

    free(detector_value);
    free(coeff);
    free(multiplicity);
    free(row);
    return stats;
}

static void run_spot_checks(void)
{
    PrimeStats s13 = analyze_prime(13, 1);
    PrimeStats s17 = analyze_prime(17, 1);
    PrimeStats s29 = analyze_prime(29, 1);

    if (s13.zero_count != 0) fail("expected empty Z_13", 13);
    if (s17.zero_count != 2) fail("expected |Z_17|=2", 17);
    if (s29.zero_count != 0) fail("expected empty Z_29", 29);
}

/* Newton iteration, used only to display E/p^(5/3), avoiding libm. */
static double cube_root_positive(double a)
{
    double x = a;
    int i;
    for (i = 0; i < 32; ++i) x = (2.0 * x + a / (x * x)) / 3.0;
    return x;
}

int main(void)
{
    int p;
    int prime_count = 0;

    run_spot_checks();
    printf("p,Z_size,S_size,E,detector_degree,zero_function,S_over_p,"
           "E_over_p_5_3,degree_over_p\n");

    for (p = PRIME_MIN; p <= PRIME_MAX; ++p) {
        PrimeStats stats;
        double root;
        double energy_scale;

        if (!is_prime(p)) continue;
        stats = analyze_prime(p, 0);
        root = cube_root_positive((double)p);
        energy_scale = (double)p * root * root;
        printf("%d,%d,%d,%" PRIu64 ",%d,%d,%.12f,%.12f,%.12f\n",
               stats.p, stats.zero_count, stats.image_size, stats.energy,
               stats.detector_degree, stats.zero_function,
               (double)stats.image_size / (double)p,
               (double)stats.energy / energy_scale,
               (double)stats.detector_degree / (double)p);
        ++prime_count;
        if (prime_count % PROGRESS_STEP == 0) {
            fprintf(stderr, "progress: %d primes complete (last p=%d)\n",
                    prime_count, p);
        }
    }

    fprintf(stderr, "complete: %d primes in [%d,%d]\n", prime_count,
            PRIME_MIN, PRIME_MAX);
    return EXIT_SUCCESS;
}
