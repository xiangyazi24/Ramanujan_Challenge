/*
 * CRON_testF_dispersion.c
 *
 * Test F for Q6420: dyadic cross-prime dispersion of centered Apéry
 * observables.  The primary statistic is the explicit fallback requested in
 * CODEX_SPEC_CRON_testF_dispersion.md.  The shifted statistic reconstructs
 * the linkage in Q6420 Sections 2.3--2.5 with d = p-q and the hard admissible
 * overlap window.
 *
 * Compile:
 *   cc -O2 -std=c11 -Wall -Wextra -pedantic \
 *      CRON_testF_dispersion.c -lm -o CRON_testF_dispersion
 *
 * The program writes machine-readable CSV records to stdout and progress
 * messages to stderr.  It has no non-system dependencies.
 */

#define _POSIX_C_SOURCE 200809L

#include <complex.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define PRIME_MIN 500
#define PRIME_MAX 2000
#define MAX_BINS 16
#define FCOUNT 4

static const long double PI_L =
    3.141592653589793238462643383279502884L;
static const char *const FNAMES[FCOUNT] = {
    "indicator", "T", "T2", "T3"
};

typedef struct {
    int p;
    int n;                  /* p-1 */
    uint32_t *t;
    unsigned char *zero;
    int *zero_pos;
    int zcount;
    int midpoint_hit;
    long double *f[FCOUNT];
    long double sumsq[FCOUNT];
    long double *no_mid;
    long double no_mid_sumsq;
    int red_n;
    long double *reduced;
    long double reduced_sumsq;
} Row;

typedef struct {
    uint64_t count;
    uint64_t support;
    long double sum;
    long double sumsq;
    long double varsum;
} Agg;

typedef struct {
    int present;
    int ip;
    int iq;
    long long cost;
} SamplePair;

static double monotonic_seconds(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
        perror("clock_gettime");
        exit(EXIT_FAILURE);
    }
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static void die(const char *message) {
    fprintf(stderr, "fatal: %s\n", message);
    exit(EXIT_FAILURE);
}

static void *xcalloc(size_t count, size_t size) {
    void *ptr = calloc(count, size);
    if (!ptr)
        die("out of memory");
    return ptr;
}

static int *prime_list(int low, int high, int *count_out) {
    unsigned char *composite = xcalloc((size_t)high + 1, 1);
    for (int i = 2; (long long)i * i <= high; ++i) {
        if (!composite[i]) {
            for (int j = i * i; j <= high; j += i)
                composite[j] = 1;
        }
    }
    int count = 0;
    for (int p = low; p <= high; ++p)
        if (p >= 2 && !composite[p])
            ++count;
    int *primes = xcalloc((size_t)count, sizeof(*primes));
    int k = 0;
    for (int p = low; p <= high; ++p)
        if (p >= 2 && !composite[p])
            primes[k++] = p;
    free(composite);
    *count_out = count;
    return primes;
}

static uint32_t mul_mod(uint32_t a, uint32_t b, uint32_t p) {
    return (uint32_t)((uint64_t)a * b % p);
}

static void build_apery_values(int p, uint32_t *t) {
    int nrow = p - 1;
    uint32_t *inv = xcalloc((size_t)p, sizeof(*inv));
    inv[1] = 1;
    for (int i = 2; i < p; ++i)
        inv[i] = (uint32_t)(p - (uint64_t)(p / i) * inv[p % i] % p);

    t[0] = 1U;
    if (nrow > 1)
        t[1] = 5U % (uint32_t)p;
    for (int n = 1; n <= p - 3; ++n) {
        uint32_t nn = (uint32_t)n;
        uint32_t n2 = mul_mod(nn, nn, (uint32_t)p);
        uint32_t n3 = mul_mod(n2, nn, (uint32_t)p);
        uint32_t coeff = (uint32_t)(
            (34ULL * n3 + 51ULL * n2 + 27ULL * nn + 5ULL) % (uint32_t)p
        );
        uint32_t left = mul_mod(coeff, t[n], (uint32_t)p);
        uint32_t right = mul_mod(n3, t[n - 1], (uint32_t)p);
        uint32_t numerator = left >= right ? left - right : left + p - right;
        uint32_t den_inv = inv[n + 1];
        den_inv = mul_mod(mul_mod(den_inv, den_inv, (uint32_t)p),
                          den_inv, (uint32_t)p);
        t[n + 1] = mul_mod(numerator, den_inv, (uint32_t)p);
    }
    free(inv);
}

static long double observable_value(uint32_t value, int which) {
    long double x = (long double)value;
    if (which == 0)
        return value == 0 ? 1.0L : 0.0L;
    if (which == 1)
        return x;
    if (which == 2)
        return x * x;
    return x * x * x;
}

static void build_row(Row *row, int p) {
    memset(row, 0, sizeof(*row));
    row->p = p;
    row->n = p - 1;
    row->red_n = row->n / 2;
    row->t = xcalloc((size_t)row->n, sizeof(*row->t));
    row->zero = xcalloc((size_t)row->n, sizeof(*row->zero));
    build_apery_values(p, row->t);

    for (int r = 0; r < row->n; ++r) {
        if (row->t[r] == 0) {
            row->zero[r] = 1;
            ++row->zcount;
        }
    }
    row->zero_pos = xcalloc((size_t)(row->zcount ? row->zcount : 1),
                            sizeof(*row->zero_pos));
    int zi = 0;
    for (int r = 0; r < row->n; ++r)
        if (row->zero[r])
            row->zero_pos[zi++] = r;
    row->midpoint_hit = row->zero[row->n / 2] != 0;

    for (int which = 0; which < FCOUNT; ++which) {
        row->f[which] = xcalloc((size_t)row->n, sizeof(*row->f[which]));
        long double mean = 0.0L;
        for (int r = 0; r < row->n; ++r)
            mean += observable_value(row->t[r], which);
        mean /= (long double)row->n;
        for (int r = 0; r < row->n; ++r) {
            long double x = observable_value(row->t[r], which) - mean;
            row->f[which][r] = x;
            row->sumsq[which] += x * x;
        }
    }

    /* Remove the exceptional midpoint zero, then center on the full cycle. */
    row->no_mid = xcalloc((size_t)row->n, sizeof(*row->no_mid));
    int generic_zeros = row->zcount - row->midpoint_hit;
    long double generic_mean = (long double)generic_zeros / row->n;
    for (int r = 0; r < row->n; ++r) {
        long double indicator =
            row->zero[r] && r != row->n / 2 ? 1.0L : 0.0L;
        row->no_mid[r] = indicator - generic_mean;
        row->no_mid_sumsq += row->no_mid[r] * row->no_mid[r];
    }

    /*
     * Project the midpoint-deleted row to reflection-even functions and use
     * the representatives 0 <= r < (p-1)/2.  The projection is explicit
     * even though the exact Apéry row is checked to be already even.
     */
    row->reduced = xcalloc((size_t)row->red_n, sizeof(*row->reduced));
    long double projected_mean = 0.0L;
    for (int r = 0; r < row->red_n; ++r) {
        int reflected = (row->n - r) % row->n;
        long double projected =
            0.5L * (row->no_mid[r] + row->no_mid[reflected]);
        row->reduced[r] = projected;
        projected_mean += projected;
    }
    projected_mean /= row->red_n;
    for (int r = 0; r < row->red_n; ++r) {
        row->reduced[r] -= projected_mean;
        row->reduced_sumsq += row->reduced[r] * row->reduced[r];
    }
}

static void free_row(Row *row) {
    free(row->t);
    free(row->zero);
    free(row->zero_pos);
    for (int which = 0; which < FCOUNT; ++which)
        free(row->f[which]);
    free(row->no_mid);
    free(row->reduced);
}

static int dyadic_bin(int d) {
    int k = 0;
    while ((1U << (k + 1)) <= (unsigned)d)
        ++k;
    return k;
}

static void agg_add(Agg *agg, long double value, long double variance,
                    int supported) {
    ++agg->count;
    if (supported)
        ++agg->support;
    agg->sum += value;
    agg->sumsq += value * value;
    agg->varsum += variance;
}

/* Exact random-injection variance for a centered full first row. */
static long double full_variance(long double qx, long double qy, int m) {
    if (m <= 1)
        return 0.0L;
    return qx * qy / (long double)(m - 1);
}

/*
 * Exact random-injection variance for a fixed subset of the first row:
 * Q_y/(m-1) * (sum x_i^2 - (sum x_i)^2/m).
 */
static long double subset_variance(const long double *x, int begin, int end,
                                   long double qy, int m) {
    if (m <= 1 || begin >= end)
        return 0.0L;
    long double sum = 0.0L;
    long double sumsq = 0.0L;
    for (int r = begin; r < end; ++r) {
        sum += x[r];
        sumsq += x[r] * x[r];
    }
    long double adjusted = sumsq - sum * sum / (long double)m;
    if (adjusted < 0.0L && fabsl(adjusted) < 1e-18L * (1.0L + sumsq))
        adjusted = 0.0L;
    return qy * adjusted / (long double)(m - 1);
}

static long double dot_prefix(const long double *x, const long double *y,
                              int n) {
    long double sum = 0.0L;
    for (int r = 0; r < n; ++r)
        sum += x[r] * y[r];
    return sum;
}

static long double dot_shifted_overlap(const long double *x,
                                       const long double *y,
                                       int n, int d) {
    if (d >= n)
        return 0.0L;
    long double sum = 0.0L;
    for (int r = d; r < n; ++r)
        sum += x[r] * y[r - d];
    return sum;
}

static long double complex dft_indicator(const Row *row, int signed_u) {
    if (signed_u == 0)
        return 0.0L + 0.0L * I;
    long double complex sum = 0.0L + 0.0L * I;
    for (int i = 0; i < row->zcount; ++i) {
        long double angle = -2.0L * PI_L * signed_u * row->zero_pos[i]
                            / (long double)row->n;
        sum += cosl(angle) + I * sinl(angle);
    }
    return sum;
}

static long double complex interval_kernel(int begin, int end,
                                           long double theta) {
    int length = end - begin;
    if (length <= 0)
        return 0.0L + 0.0L * I;
    theta -= nearbyintl(theta);
    long double denom = sinl(PI_L * theta);
    if (fabsl(denom) < 1e-15L) {
        long double angle = 2.0L * PI_L * theta * begin;
        return length * (cosl(angle) + I * sinl(angle));
    }
    long double scale = sinl(PI_L * theta * length) / denom;
    long double angle = PI_L * theta * (begin + end - 1);
    return scale * (cosl(angle) + I * sinl(angle));
}

static long double near_resonant_shifted(const Row *p, const Row *q, int d) {
    if (d >= p->n)
        return 0.0L;
    long double xscale = 0.5L * (p->p + q->p);
    int cutoff = (int)floorl(xscale / d);
    int signed_cap_p = (p->n - 1) / 2;
    int signed_cap_q = (q->n - 1) / 2;
    if (cutoff > signed_cap_p)
        cutoff = signed_cap_p;
    if (cutoff > signed_cap_q)
        cutoff = signed_cap_q;
    long double complex sum = 0.0L + 0.0L * I;
    for (int u = 1; u <= cutoff; ++u) {
        long double complex hp = dft_indicator(p, u);
        long double complex hq_minus = conjl(dft_indicator(q, u));
        long double phase_angle = 2.0L * PI_L * u * d / q->n;
        long double complex phase =
            cosl(phase_angle) + I * sinl(phase_angle);
        long double theta =
            (long double)u / p->n - (long double)u / q->n;
        long double complex kernel = interval_kernel(d, p->n, theta);
        sum += 2.0L * creall(hp * hq_minus * phase * kernel);
    }
    return creall(sum) / ((long double)p->n * q->n);
}

static long double cyclic_indicator_at_shift(const Row *p, const Row *q,
                                             int delta) {
    int matches = 0;
    for (int i = 0; i < p->zcount; ++i) {
        int s = (p->zero_pos[i] + delta) % q->n;
        if (s < 0)
            s += q->n;
        matches += q->zero[s] != 0;
    }
    int qzeros_in_window = 0;
    for (int i = 0; i < q->zcount; ++i) {
        int r = (q->zero_pos[i] - delta) % q->n;
        if (r < 0)
            r += q->n;
        qzeros_in_window += r < p->n;
    }
    long double ap = (long double)p->zcount / p->n;
    long double aq = (long double)q->zcount / q->n;
    return matches - ap * qzeros_in_window - aq * p->zcount
           + (long double)p->n * ap * aq;
}

static long double fourier_reconstruct(const Row *p, const Row *q,
                                       int delta, int begin, int end) {
    long double complex *hp = xcalloc((size_t)p->n, sizeof(*hp));
    long double complex *hq = xcalloc((size_t)q->n, sizeof(*hq));
    for (int u = 1; u < p->n; ++u)
        hp[u] = dft_indicator(p, u);
    for (int v = 1; v < q->n; ++v)
        hq[v] = dft_indicator(q, v);

    long double complex total = 0.0L + 0.0L * I;
    for (int u = 1; u < p->n; ++u) {
        for (int v = 1; v < q->n; ++v) {
            long double phase_angle =
                2.0L * PI_L * v * delta / q->n;
            long double complex phase =
                cosl(phase_angle) + I * sinl(phase_angle);
            long double theta =
                (long double)u / p->n + (long double)v / q->n;
            long double complex kernel = interval_kernel(begin, end, theta);
            total += hp[u] * hq[v] * phase * kernel;
        }
    }
    free(hp);
    free(hq);
    return creall(total) / ((long double)p->n * q->n);
}

static void print_agg(const char *family, int bin, const char *observable,
                      const Agg *agg) {
    long double rms = agg->count
        ? sqrtl(agg->sumsq / (long double)agg->count) : 0.0L;
    long double bench = agg->count
        ? sqrtl(agg->varsum / (long double)agg->count) : 0.0L;
    long double ratio = agg->varsum > 0.0L
        ? sqrtl(agg->sumsq / agg->varsum) : 0.0L;
    long double zscore = agg->varsum > 0.0L
        ? agg->sum / sqrtl(agg->varsum) : 0.0L;
    printf("STAT,%s,%d,%" PRIu64 ",%" PRIu64 ",%s,"
           "%.18Le,%.18Le,%.18Le,%.18Le,%.18Le,%.18Le,%.18Le\n",
           family, bin, agg->count, agg->support, observable,
           agg->sum, rms, bench, ratio, zscore, agg->sumsq, agg->varsum);
}

static void select_sample(SamplePair *sample, int ip, int iq,
                          const Row *rows) {
    long long cost = (long long)rows[ip].n * rows[iq].n;
    if (!sample->present || cost < sample->cost) {
        sample->present = 1;
        sample->ip = ip;
        sample->iq = iq;
        sample->cost = cost;
    }
}

static void validate_p17(void) {
    Row row;
    build_row(&row, 17);
    if (row.zcount != 2 || row.zero_pos[0] != 3 || row.zero_pos[1] != 13)
        die("p=17 zero set is not {3,13}");
    free_row(&row);
}

int main(void) {
    double started = monotonic_seconds();
    double last_progress = started;
    fprintf(stderr, "[testF] start: primes in [%d,%d]\n",
            PRIME_MIN, PRIME_MAX);

    validate_p17();
    fprintf(stderr, "[testF] ground truth p=17: Z_17={3,13}\n");

    int nprimes = 0;
    int *primes = prime_list(PRIME_MIN, PRIME_MAX, &nprimes);
    Row *rows = xcalloc((size_t)nprimes, sizeof(*rows));
    int symmetry_failures = 0;
    int parity_failures_all = 0;
    int parity_failures_sample20 = 0;
    int midpoint_hits = 0;
    for (int i = 0; i < nprimes; ++i) {
        build_row(&rows[i], primes[i]);
        for (int r = 0; r < rows[i].n; ++r) {
            int reflected = (rows[i].n - r) % rows[i].n;
            if (rows[i].t[r] != rows[i].t[reflected]) {
                ++symmetry_failures;
                break;
            }
        }
        int parity_ok = ((rows[i].zcount & 1) != 0) == rows[i].midpoint_hit;
        if (!parity_ok)
            ++parity_failures_all;
        if (i < 20 && !parity_ok)
            ++parity_failures_sample20;
        midpoint_hits += rows[i].midpoint_hit;
        double now = monotonic_seconds();
        if (now - last_progress >= 30.0) {
            fprintf(stderr, "[testF] rows: %d/%d (%.1fs)\n",
                    i + 1, nprimes, now - started);
            last_progress = now;
        }
    }
    if (symmetry_failures || parity_failures_all || parity_failures_sample20)
        die("reflection/parity validation failed");
    fprintf(stderr,
            "[testF] built %d rows; reflection exact; parity law passed "
            "sample 20 and all %d rows\n", nprimes, nprimes);

    Agg aligned[MAX_BINS][FCOUNT] = {{{0}}};
    Agg shifted[MAX_BINS][FCOUNT] = {{{0}}};
    Agg no_mid[MAX_BINS] = {{0}};
    Agg reduced[MAX_BINS] = {{0}};
    Agg near_shifted[MAX_BINS] = {{0}};
    Agg cyclic_exact[MAX_BINS] = {{0}};
    Agg cyclic_profile[MAX_BINS] = {{0}};
    uint64_t profile_arg_exact[MAX_BINS] = {0};
    SamplePair aligned_samples[MAX_BINS] = {{0}};
    SamplePair shifted_samples[MAX_BINS] = {{0}};

    uint64_t total_pairs = (uint64_t)nprimes * (nprimes - 1) / 2;
    uint64_t pairs_done = 0;
    for (int ip = 0; ip < nprimes; ++ip) {
        const Row *p = &rows[ip];
        for (int iq = ip + 1; iq < nprimes; ++iq) {
            const Row *q = &rows[iq];
            int d = q->p - p->p;
            int bin = dyadic_bin(d);
            if (bin >= MAX_BINS)
                die("dyadic bin overflow");
            if (p->zcount > 0 && q->zcount > 0)
                select_sample(&aligned_samples[bin], ip, iq, rows);
            if (d < p->n && p->zcount > 0 && q->zcount > 0)
                select_sample(&shifted_samples[bin], ip, iq, rows);

            for (int which = 0; which < FCOUNT; ++which) {
                long double c = dot_prefix(p->f[which], q->f[which], p->n);
                long double var = full_variance(
                    p->sumsq[which], q->sumsq[which], q->n
                );
                agg_add(&aligned[bin][which], c, var, var > 0.0L);

                long double cs = dot_shifted_overlap(
                    p->f[which], q->f[which], p->n, d
                );
                long double vars = subset_variance(
                    p->f[which], d, p->n, q->sumsq[which], q->n
                );
                agg_add(&shifted[bin][which], cs, vars, d < p->n);
            }

            long double cnm = dot_prefix(p->no_mid, q->no_mid, p->n);
            long double vnm = full_variance(
                p->no_mid_sumsq, q->no_mid_sumsq, q->n
            );
            agg_add(&no_mid[bin], cnm, vnm, vnm > 0.0L);

            long double cr = dot_prefix(p->reduced, q->reduced, p->red_n);
            long double vr = full_variance(
                p->reduced_sumsq, q->reduced_sumsq, q->red_n
            );
            agg_add(&reduced[bin], cr, vr, vr > 0.0L);

            long double near = near_resonant_shifted(p, q, d);
            agg_add(&near_shifted[bin], near, 0.0L, d < p->n);

            long double c_exact = cyclic_indicator_at_shift(p, q, -d);
            long double vcyc = full_variance(
                p->sumsq[0], q->sumsq[0], q->n
            );
            long double maxabs = -1.0L;
            for (int delta = -d; delta <= d; ++delta) {
                long double value = cyclic_indicator_at_shift(p, q, delta);
                long double av = fabsl(value);
                if (av > maxabs)
                    maxabs = av;
            }
            agg_add(&cyclic_exact[bin], c_exact, vcyc, vcyc > 0.0L);
            agg_add(&cyclic_profile[bin], maxabs, vcyc, vcyc > 0.0L);
            long double plus_endpoint = cyclic_indicator_at_shift(p, q, d);
            long double endpoint_max =
                fmaxl(fabsl(c_exact), fabsl(plus_endpoint));
            if (vcyc > 0.0L
                && maxabs - endpoint_max
                   <= 1e-18L * (1.0L + maxabs))
                ++profile_arg_exact[bin];

            ++pairs_done;
            double now = monotonic_seconds();
            if (now - last_progress >= 30.0) {
                fprintf(stderr,
                        "[testF] pairs: %" PRIu64 "/%" PRIu64
                        " (%.1fs)\n",
                        pairs_done, total_pairs, now - started);
                last_progress = now;
            }
        }
    }
    fprintf(stderr, "[testF] pair loop complete; Fourier checks next\n");

    printf("META,prime_min,%d\n", PRIME_MIN);
    printf("META,prime_max,%d\n", PRIME_MAX);
    printf("META,prime_count,%d\n", nprimes);
    printf("META,pair_count,%" PRIu64 "\n", total_pairs);
    printf("META,midpoint_hit_count,%d\n", midpoint_hits);
    printf("VALID,p17,1,3;13\n");
    printf("VALID,parity_sample20,%d,0\n", 20);
    printf("VALID,parity_all,%d,0\n", nprimes);
    printf("VALID,reflection_all,%d,0\n", nprimes);

    for (int bin = 0; bin < MAX_BINS; ++bin) {
        if (!aligned[bin][0].count)
            continue;
        for (int which = 0; which < FCOUNT; ++which) {
            print_agg("aligned", bin, FNAMES[which], &aligned[bin][which]);
            print_agg("shifted_overlap", bin, FNAMES[which],
                      &shifted[bin][which]);
        }
        print_agg("indicator_no_mid", bin, "indicator", &no_mid[bin]);
        print_agg("indicator_reflection_reduced", bin, "indicator",
                  &reduced[bin]);
        print_agg("near_resonant_shifted", bin, "indicator",
                  &near_shifted[bin]);
        print_agg("indicator_cyclic_exact_shift", bin, "indicator",
                  &cyclic_exact[bin]);
        print_agg("indicator_cyclic_profile_max", bin, "indicator",
                  &cyclic_profile[bin]);
        printf("PROFILE,%d,%" PRIu64 ",%" PRIu64 "\n", bin,
               profile_arg_exact[bin], cyclic_profile[bin].support);
    }

    for (int bin = 0; bin < MAX_BINS; ++bin) {
        if (aligned_samples[bin].present) {
            const Row *p = &rows[aligned_samples[bin].ip];
            const Row *q = &rows[aligned_samples[bin].iq];
            long double direct = dot_prefix(p->f[0], q->f[0], p->n);
            long double recon = fourier_reconstruct(p, q, 0, 0, p->n);
            printf("FOURIER,aligned,%d,%d,%d,%.18Le,%.18Le,%.18Le\n",
                   bin, p->p, q->p, direct, recon, fabsl(direct - recon));
        }
        if (shifted_samples[bin].present) {
            const Row *p = &rows[shifted_samples[bin].ip];
            const Row *q = &rows[shifted_samples[bin].iq];
            int d = q->p - p->p;
            long double direct = dot_shifted_overlap(
                p->f[0], q->f[0], p->n, d
            );
            long double recon = fourier_reconstruct(
                p, q, -d, d, p->n
            );
            printf("FOURIER,shifted_overlap,%d,%d,%d,"
                   "%.18Le,%.18Le,%.18Le\n",
                   bin, p->p, q->p, direct, recon,
                   fabsl(direct - recon));
        }
        double now = monotonic_seconds();
        if (now - last_progress >= 30.0) {
            fprintf(stderr, "[testF] Fourier bins through k=%d (%.1fs)\n",
                    bin, now - started);
            last_progress = now;
        }
    }

    printf("META,elapsed_seconds,%.6f\n", monotonic_seconds() - started);
    fprintf(stderr, "[testF] complete in %.2fs\n",
            monotonic_seconds() - started);

    for (int i = 0; i < nprimes; ++i)
        free_row(&rows[i]);
    free(rows);
    free(primes);
    return EXIT_SUCCESS;
}
