// Q7320: exact anti-diagonal Apéry-continuant audit.
//
// This program checks, independently of the proof in the accompanying answer:
//   * the polynomial reciprocity N_h(-h-1-X)=(-1)^(h-1)N_h(X);
//   * the exact factor (2X+h+1) for every even h in the polynomial range;
//   * divided/cleared Apéry recurrences, reflection for both basis solutions,
//     the Casoratian, and direct binomial values;
//   * both orientations of every top-shell zero for primes up to --prime-max;
//   * the quotient identity modulo p linking N_h/p to the first reflection
//     defect modulo p^2.
//
// Build:
//   g++ -std=c++17 -O3 -Wall -Wextra -Wpedantic \
//     q7320_antidiagonal_check.cpp -o q7320_check
// Run:
//   ./q7320_check --poly-max 20 --prime-max 10000

#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;
using Poly = std::vector<cpp_int>; // increasing coefficient order

static void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

static cpp_int abs_int(cpp_int x) { return x < 0 ? -x : x; }

static cpp_int gcd_int(cpp_int a, cpp_int b) {
    a = abs_int(a);
    b = abs_int(b);
    while (b != 0) {
        cpp_int r = a % b;
        a = b;
        b = r;
    }
    return a;
}

static void trim(Poly &f) {
    while (f.size() > 1 && f.back() == 0) f.pop_back();
}

static Poly poly_add(const Poly &a, const Poly &b) {
    Poly out(std::max(a.size(), b.size()), 0);
    for (std::size_t i = 0; i < a.size(); ++i) out[i] += a[i];
    for (std::size_t i = 0; i < b.size(); ++i) out[i] += b[i];
    trim(out);
    return out;
}

static Poly poly_sub(const Poly &a, const Poly &b) {
    Poly out(std::max(a.size(), b.size()), 0);
    for (std::size_t i = 0; i < a.size(); ++i) out[i] += a[i];
    for (std::size_t i = 0; i < b.size(); ++i) out[i] -= b[i];
    trim(out);
    return out;
}

static Poly poly_scale(Poly a, const cpp_int &c) {
    for (cpp_int &x : a) x *= c;
    trim(a);
    return a;
}

static Poly poly_mul(const Poly &a, const Poly &b) {
    Poly out(a.size() + b.size() - 1, 0);
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j)
            out[i + j] += a[i] * b[j];
    trim(out);
    return out;
}

static Poly poly_pow(Poly base, unsigned exponent) {
    Poly result{1};
    while (exponent != 0) {
        if (exponent & 1U) result = poly_mul(result, base);
        exponent >>= 1U;
        if (exponent != 0) base = poly_mul(base, base);
    }
    return result;
}

static Poly affine_compose(const Poly &f, const cpp_int &a,
                           const cpp_int &b) {
    const Poly linear{a, b};
    Poly power{1};
    Poly out{0};
    for (const cpp_int &coefficient : f) {
        out = poly_add(out, poly_scale(power, coefficient));
        power = poly_mul(power, linear);
    }
    trim(out);
    return out;
}

static Poly P_shift(u32 shift) {
    const Poly t{cpp_int(shift), cpp_int(1)};
    const Poly first = poly_add(poly_scale(t, 2), Poly{1});
    Poly second = poly_scale(poly_mul(t, t), 17);
    second = poly_add(second, poly_scale(t, 17));
    second = poly_add(second, Poly{5});
    return poly_mul(first, second);
}

static Poly linear_power(u32 shift, unsigned exponent) {
    return poly_pow(Poly{cpp_int(shift), cpp_int(1)}, exponent);
}

static Poly divide_by_antidiagonal(const Poly &f, u32 h) {
    require(h % 2 == 0, "anti-diagonal division requested for odd h");
    require(f.size() >= 2, "constant polynomial cannot have linear factor");
    const cpp_int b = h + 1;
    const std::size_t degree = f.size() - 1;
    Poly q(degree, 0);
    require(f[degree] % 2 == 0, "nonintegral leading quotient coefficient");
    q[degree - 1] = f[degree] / 2;
    for (std::size_t k = degree - 1; k >= 1; --k) {
        const cpp_int numerator = f[k] - b * q[k];
        require(numerator % 2 == 0, "nonintegral quotient coefficient");
        q[k - 1] = numerator / 2;
        if (k == 1) break;
    }
    require(f[0] == b * q[0], "nonzero remainder after anti-diagonal division");
    trim(q);
    return q;
}

static cpp_int content(const Poly &f) {
    cpp_int out = 0;
    for (const cpp_int &coefficient : f) out = gcd_int(out, coefficient);
    return out;
}

static unsigned bit_length(const cpp_int &x0) {
    const cpp_int x = abs_int(x0);
    if (x == 0) return 0;
    return boost::multiprecision::msb(x) + 1;
}

static std::vector<Poly> build_N_polynomials(u32 poly_max) {
    std::vector<Poly> N(poly_max + 1);
    N[0] = Poly{0};
    N[1] = Poly{1};
    cpp_int lead_prev = 0;
    cpp_int lead_cur = 1;
    for (u32 h = 1; h < poly_max; ++h) {
        N[h + 1] = poly_sub(poly_mul(P_shift(h), N[h]),
                            poly_mul(linear_power(h, 6), N[h - 1]));
        const std::size_t expected_degree = 3ULL * h;
        require(N[h + 1].size() - 1 == expected_degree,
                "unexpected continuant degree");
        const cpp_int lead_next = 34 * lead_cur - lead_prev;
        require(N[h + 1].back() == lead_next,
                "leading coefficient recurrence failed");
        lead_prev = lead_cur;
        lead_cur = lead_next;
    }
    return N;
}

static void polynomial_audit(u32 poly_max) {
    const std::vector<Poly> N = build_N_polynomials(poly_max);
    std::cout << "POLYNOMIAL_AUDIT poly_max=" << poly_max << "\n";
    std::cout << "h degree leading content quotient_degree quotient_content max_coeff_bits\n";
    for (u32 h = 1; h <= poly_max; ++h) {
        const Poly reflected = affine_compose(N[h], -cpp_int(h) - 1, -1);
        const Poly expected = poly_scale(N[h], (h % 2 == 1) ? 1 : -1);
        require(reflected == expected, "polynomial reflection reciprocity failed");

        unsigned max_bits = 0;
        for (const cpp_int &coefficient : N[h])
            max_bits = std::max(max_bits, bit_length(coefficient));

        std::cout << h << ' ' << (N[h].size() - 1) << ' '
                  << N[h].back() << ' ' << content(N[h]) << ' ';
        if (h % 2 == 0) {
            const Poly M = divide_by_antidiagonal(N[h], h);
            const Poly reflected_M = affine_compose(M, -cpp_int(h) - 1, -1);
            require(reflected_M == M, "quotient reflection symmetry failed");
            std::cout << (M.size() - 1) << ' ' << content(M) << ' ';
        } else {
            std::cout << "- - ";
        }
        std::cout << max_bits << "\n";
    }
    std::cout << "POLYNOMIAL_AUDIT PASS\n";
}

static inline u64 mul_mod(u64 a, u64 b, u64 modulus) {
    return static_cast<u64>((static_cast<u128>(a) * b) % modulus);
}

static u64 pow_mod(u64 a, u64 exponent, u64 modulus) {
    u64 out = 1 % modulus;
    a %= modulus;
    while (exponent != 0) {
        if (exponent & 1U) out = mul_mod(out, a, modulus);
        a = mul_mod(a, a, modulus);
        exponent >>= 1U;
    }
    return out;
}

static u64 sub_mod(u64 a, u64 b, u64 modulus) {
    return a >= b ? a - b : a + modulus - b;
}

static u64 P_mod(u64 m, u64 modulus) {
    const u64 x = m % modulus;
    const u64 x2 = mul_mod(x, x, modulus);
    const u64 x3 = mul_mod(x2, x, modulus);
    u64 out = mul_mod(34 % modulus, x3, modulus);
    out = (out + mul_mod(51 % modulus, x2, modulus)) % modulus;
    out = (out + mul_mod(27 % modulus, x, modulus)) % modulus;
    return (out + 5) % modulus;
}

static u64 cube_mod(u64 x, u64 modulus) {
    return mul_mod(mul_mod(x % modulus, x % modulus, modulus),
                   x % modulus, modulus);
}

static u64 sixth_mod(u64 x, u64 modulus) {
    const u64 y = cube_mod(x, modulus);
    return mul_mod(y, y, modulus);
}

static std::int64_t extended_gcd(std::int64_t a, std::int64_t b,
                                 std::int64_t &x, std::int64_t &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    std::int64_t x1 = 0, y1 = 0;
    const std::int64_t g = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

static u64 inverse_mod(u64 a, u64 modulus) {
    require(modulus <= static_cast<u64>(std::numeric_limits<std::int64_t>::max()),
            "modulus too large for reference inverse");
    std::int64_t x = 0, y = 0;
    const std::int64_t g = extended_gcd(static_cast<std::int64_t>(a),
                                        static_cast<std::int64_t>(modulus),
                                        x, y);
    require(g == 1, "noninvertible recurrence denominator");
    x %= static_cast<std::int64_t>(modulus);
    if (x < 0) x += static_cast<std::int64_t>(modulus);
    return static_cast<u64>(x);
}

static std::vector<u32> primes_upto(u32 limit) {
    std::vector<bool> is_prime(limit + 1, true);
    if (limit >= 0) is_prime[0] = false;
    if (limit >= 1) is_prime[1] = false;
    for (u32 p = 2; static_cast<u64>(p) * p <= limit; ++p) {
        if (!is_prime[p]) continue;
        for (u32 q = p * p; q <= limit; q += p) is_prime[q] = false;
    }
    std::vector<u32> out;
    for (u32 p = 2; p <= limit; ++p)
        if (is_prime[p]) out.push_back(p);
    return out;
}

struct BasisData {
    std::vector<u64> b2;
    std::vector<u64> c2;
};

static BasisData build_basis_mod_p2(u32 p) {
    const u64 modulus = static_cast<u64>(p) * p;
    BasisData out;
    out.b2.assign(p, 0);
    out.c2.assign(p, 0);
    out.b2[0] = 1;
    out.b2[1] = 5;
    out.c2[0] = 0;
    out.c2[1] = 1;
    for (u32 m = 1; m + 1 < p; ++m) {
        const u64 inv = inverse_mod(m + 1, modulus);
        const u64 inv_cube = cube_mod(inv, modulus);
        const auto next = [&](const std::vector<u64> &v) {
            const u64 numerator = sub_mod(
                mul_mod(P_mod(m, modulus), v[m], modulus),
                mul_mod(cube_mod(m, modulus), v[m - 1], modulus),
                modulus);
            return mul_mod(numerator, inv_cube, modulus);
        };
        out.b2[m + 1] = next(out.b2);
        out.c2[m + 1] = next(out.c2);
    }
    return out;
}

static u64 direct_apery_mod(u32 n, u32 p,
                            const std::vector<u64> &inverse) {
    u64 choose_n_j = 1;
    u64 choose_n_plus_j_j = 1;
    u64 sum = 0;
    for (u32 j = 0; j <= n; ++j) {
        const u64 a2 = mul_mod(choose_n_j, choose_n_j, p);
        const u64 c2 = mul_mod(choose_n_plus_j_j,
                               choose_n_plus_j_j, p);
        sum = (sum + mul_mod(a2, c2, p)) % p;
        if (j == n) break;
        choose_n_j = mul_mod(choose_n_j, n - j, p);
        choose_n_j = mul_mod(choose_n_j, inverse[j + 1], p);
        choose_n_plus_j_j = mul_mod(choose_n_plus_j_j,
                                    n + j + 1, p);
        choose_n_plus_j_j = mul_mod(choose_n_plus_j_j,
                                    inverse[j + 1], p);
    }
    return sum;
}

static u64 N_mod(u32 u, u32 h, u64 modulus) {
    if (h == 0) return 0;
    u64 previous = 0;
    u64 current = 1 % modulus;
    for (u32 step = 1; step < h; ++step) {
        const u64 x = static_cast<u64>(u) + step;
        const u64 next = sub_mod(
            mul_mod(P_mod(x, modulus), current, modulus),
            mul_mod(sixth_mod(x, modulus), previous, modulus), modulus);
        previous = current;
        current = next;
    }
    return current;
}

static u64 interval_factor_mod(u32 u, u32 v, u32 p) {
    u64 out = 1;
    for (u32 x = u + 1; x <= v; ++x)
        out = mul_mod(out, cube_mod(x, p), p);
    return out;
}

struct ScanCounts {
    u64 primes = 0;
    u64 recurrence_positions = 0;
    u64 reflection_positions = 0;
    u64 direct_positions = 0;
    u64 anti_diagonal_all_u = 0;
    u64 target_orientations = 0;
    u64 target_pairs = 0;
    u64 lower_orientations = 0;
    u64 upper_orientations = 0;
    u64 central_orientations = 0;
    u64 p2_extra_orientations = 0;
    u64 p2_extra_pairs = 0;
};

struct Example {
    u32 p = 0;
    u32 n = 0;
    u32 j = 0;
    u32 k = 0;
    u32 u = 0;
    u32 h = 0;
    u32 Mmod = 0;
    u32 defect = 0;
};

static void prime_audit(u32 prime_max) {
    ScanCounts counts;
    std::vector<Example> target_examples;
    std::vector<Example> extra_examples;

    for (u32 p : primes_upto(prime_max)) {
        if (p < 5) continue;
        ++counts.primes;
        const u64 modulus = static_cast<u64>(p) * p;
        const BasisData basis = build_basis_mod_p2(p);

        std::vector<u64> inverse(p, 0);
        inverse[1] = 1;
        for (u32 m = 2; m < p; ++m)
            inverse[m] = p - mul_mod(p / m, inverse[p % m], p);

        u64 factorial = 1;
        u64 Bprev = 1 % p, Bcur = 5 % p;
        u64 Dprev = 0, Dcur = 1 % p;
        for (u32 m = 0; m < p; ++m) {
            const u64 b = basis.b2[m] % p;
            const u64 c = basis.c2[m] % p;
            const u64 factorial_cube = cube_mod(factorial, p);
            const u64 B = mul_mod(factorial_cube, b, p);
            const u64 D = mul_mod(factorial_cube, c, p);
            if (m == 0) {
                require(B == Bprev && D == Dprev,
                        "cleared initial state mismatch");
            } else if (m == 1) {
                require(B == Bcur && D == Dcur,
                        "cleared first state mismatch");
            } else {
                require(B == Bcur && D == Dcur,
                        "divided/cleared recurrence mismatch");
            }
            if (m + 1 < p) {
                const u64 factorial_sixth = sixth_mod(factorial, p);
                const u64 next_factorial = mul_mod(factorial, m + 1, p);
                const u64 next_cube = cube_mod(next_factorial, p);
                const u64 Bnext_direct = mul_mod(next_cube,
                                                 basis.b2[m + 1] % p, p);
                const u64 Dnext_direct = mul_mod(next_cube,
                                                 basis.c2[m + 1] % p, p);
                const u64 W = sub_mod(mul_mod(B, Dnext_direct, p),
                                      mul_mod(Bnext_direct, D, p), p);
                require(W == factorial_sixth, "cleared Wronskian failed");
            }
            ++counts.recurrence_positions;

            const u32 reflected = p - 1 - m;
            require(b == basis.b2[reflected] % p,
                    "Apéry reflection failed");
            require(c == basis.c2[reflected] % p,
                    "companion reflection failed");
            ++counts.reflection_positions;

            if (m <= std::min<u32>(12, p - 1)) {
                require(direct_apery_mod(m, p, inverse) == b,
                        "direct binomial check failed");
                ++counts.direct_positions;
            }

            if (m + 1 < p) {
                const u64 polynomial = P_mod(m + 1, p);
                const u64 m1six = sixth_mod(m + 1, p);
                const u64 Bnext = sub_mod(mul_mod(polynomial, Bcur, p),
                                          mul_mod(m1six, Bprev, p), p);
                const u64 Dnext = sub_mod(mul_mod(polynomial, Dcur, p),
                                          mul_mod(m1six, Dprev, p), p);
                Bprev = Bcur;
                Bcur = Bnext;
                Dprev = Dcur;
                Dcur = Dnext;
                factorial = next_factorial;
            }
        }

        // Exhaustive anti-diagonal quotient/defect identity for small primes.
        if (p <= 251) {
            for (u32 u = 0; 2 * u + 1 < p; ++u) {
                const u32 v = p - 1 - u;
                const u32 h = v - u;
                require(h % 2 == 0, "odd anti-diagonal gap");
                const u64 N2 = N_mod(u, h, modulus);
                require(N2 % p == 0,
                        "universal anti-diagonal factor failed");
                const u64 Mmod = (N2 / p) % p;
                const u64 F = interval_factor_mod(u, v, p);
                const u64 db_num = (basis.b2[v] + modulus - basis.b2[u]) % modulus;
                const u64 dc_num = (basis.c2[v] + modulus - basis.c2[u]) % modulus;
                require(db_num % p == 0 && dc_num % p == 0,
                        "reflection defect is not integral modulo p");
                const u64 db = db_num / p;
                const u64 dc = dc_num / p;
                const u64 determinant_jet = sub_mod(
                    mul_mod(basis.b2[u] % p, dc, p),
                    mul_mod(basis.c2[u] % p, db, p), p);
                require(Mmod == mul_mod(F, determinant_jet, p),
                        "general quotient/reflection-jet identity failed");
                ++counts.anti_diagonal_all_u;
            }
        }

        std::vector<bool> pair_seen((p + 1) / 2, false);
        std::vector<bool> extra_pair_seen((p + 1) / 2, false);
        for (u32 j = 0; j < p; ++j) {
            if (basis.b2[j] % p != 0) continue;
            ++counts.target_orientations;
            const u32 n = p + j;
            const u32 k = p - 1 - j;
            require(p > n / 2, "top-shell inequality failed");
            require(basis.b2[k] % p == 0,
                    "reflected target endpoint is not a zero");
            if (j < k) ++counts.lower_orientations;
            else if (j > k) ++counts.upper_orientations;
            else ++counts.central_orientations;

            const u32 u = std::min(j, k);
            if (!pair_seen[u]) {
                pair_seen[u] = true;
                ++counts.target_pairs;
            }
            if (j == k) continue;
            const u32 v = std::max(j, k);
            const u32 h = v - u;
            require(h % 2 == 0, "target anti-diagonal gap is odd");
            const u64 N2 = N_mod(u, h, modulus);
            require(N2 % p == 0, "target continuant is not divisible by p");
            const u64 Mmod = (N2 / p) % p;
            const u64 F = interval_factor_mod(u, v, p);
            const u64 db_num = (basis.b2[v] + modulus - basis.b2[u]) % modulus;
            require(db_num % p == 0,
                    "target reflection defect is not integral");
            const u64 db = db_num / p;
            const u64 predicted = sub_mod(
                0, mul_mod(F, mul_mod(basis.c2[u] % p, db, p), p), p);
            require(Mmod == predicted,
                    "zero-fiber quotient/reflection-defect identity failed");
            require((Mmod == 0) == (db == 0),
                    "extra p factor is not equivalent to p^2 reflection");

            const Example example{p, n, j, k, u, h,
                                  static_cast<u32>(Mmod),
                                  static_cast<u32>(db)};
            if (target_examples.size() < 16) target_examples.push_back(example);
            if (Mmod == 0) {
                ++counts.p2_extra_orientations;
                if (!extra_pair_seen[u]) {
                    extra_pair_seen[u] = true;
                    ++counts.p2_extra_pairs;
                    if (extra_examples.size() < 32) extra_examples.push_back(example);
                }
            }
        }
    }

    std::cout << "PRIME_AUDIT prime_max=" << prime_max << "\n";
    std::cout << "primes=" << counts.primes
              << " recurrence_positions=" << counts.recurrence_positions
              << " reflection_positions=" << counts.reflection_positions
              << " direct_binomial_positions=" << counts.direct_positions
              << " all_u_jet_checks=" << counts.anti_diagonal_all_u << "\n";
    std::cout << "target_orientations=" << counts.target_orientations
              << " target_pairs=" << counts.target_pairs
              << " j_lt_k=" << counts.lower_orientations
              << " j_gt_k=" << counts.upper_orientations
              << " j_eq_k=" << counts.central_orientations << "\n";
    std::cout << "p2_extra_orientations=" << counts.p2_extra_orientations
              << " p2_extra_pairs=" << counts.p2_extra_pairs << "\n";

    std::cout << "TARGET_EXAMPLES p n j k u h Mmod defect\n";
    for (const Example &e : target_examples)
        std::cout << e.p << ' ' << e.n << ' ' << e.j << ' ' << e.k << ' '
                  << e.u << ' ' << e.h << ' ' << e.Mmod << ' '
                  << e.defect << "\n";
    std::cout << "P2_EXTRA_EXAMPLES p n j k u h Mmod defect\n";
    for (const Example &e : extra_examples)
        std::cout << e.p << ' ' << e.n << ' ' << e.j << ' ' << e.k << ' '
                  << e.u << ' ' << e.h << ' ' << e.Mmod << ' '
                  << e.defect << "\n";
    std::cout << "PRIME_AUDIT PASS\n";
}

static u32 parse_u32(const char *text) {
    const unsigned long value = std::stoul(text);
    require(value <= std::numeric_limits<u32>::max(), "argument too large");
    return static_cast<u32>(value);
}

int main(int argc, char **argv) {
    try {
        u32 poly_max = 20;
        u32 prime_max = 10000;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--poly-max") {
                require(i + 1 < argc, "missing --poly-max value");
                poly_max = parse_u32(argv[++i]);
            } else if (arg == "--prime-max") {
                require(i + 1 < argc, "missing --prime-max value");
                prime_max = parse_u32(argv[++i]);
            } else {
                throw std::runtime_error("unknown argument: " + arg);
            }
        }
        require(poly_max >= 2, "poly_max must be at least 2");
        require(prime_max >= 5, "prime_max must be at least 5");
        polynomial_audit(poly_max);
        prime_audit(prime_max);
        std::cout << "Q7320_AUDIT PASS\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Q7320_AUDIT FAIL: " << error.what() << "\n";
        return 1;
    }
}
