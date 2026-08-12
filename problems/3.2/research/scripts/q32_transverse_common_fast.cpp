// Exact high-prime scanner for Apéry transverse common pairs.
//
// Build (inside a Sage environment, or with system FLINT):
//   g++ -O3 -DNDEBUG -std=c++17 q32_transverse_common_fast.cpp \
//       -lflint -lgmp -lmpfr -o q32_transverse_common_fast
//
// The scanner is deliberately prime-by-prime.  For one prime p it first
// computes b_0,...,b_{p-1} in O(p) word operations.  If no b_r vanishes,
// there cannot be a common pair and no power-series arithmetic is done.
// Otherwise it truncates at one past the last zero of b and lets FLINT compute
//
//   g = F^{-2} (1 - 34 t + t^2)^{-1/2}
//
// modulo p using one inverse_series and two mullow operations.  Decisions are
// exact modulo p; floating point is used only for timing output.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if __has_include(<flint/nmod_poly.h>)
#include <flint/nmod_poly.h>
#elif __has_include(<nmod_poly.h>)
#include <nmod_poly.h>
#else
#error "FLINT nmod_poly.h not found"
#endif

namespace fs = std::filesystem;
using u64 = std::uint64_t;
using u128 = __uint128_t;

static constexpr const char *ALGO_VERSION = "q32-transverse-common-fast-v1";

static inline u64 add_mod(u64 a, u64 b, u64 p) {
    u64 c = a + b;
    return c >= p ? c - p : c;
}

static inline u64 sub_mod(u64 a, u64 b, u64 p) {
    return a >= b ? a - b : a + p - b;
}

static inline u64 mul_mod(u64 a, u64 b, u64 p) {
    return static_cast<u64>((static_cast<u128>(a) * b) % p);
}

static inline u64 cube_mod(u64 a, u64 p) {
    return mul_mod(mul_mod(a, a, p), a, p);
}

static bool is_prime(u64 n) {
    if (n < 2) return false;
    if ((n & 1U) == 0) return n == 2;
    for (u64 d = 3; d * d <= n; d += 2)
        if (n % d == 0) return false;
    return true;
}

static std::vector<u64> primes_up_to(u64 lo, u64 hi) {
    if (hi < 2 || hi < lo) return {};
    std::vector<bool> sieve(static_cast<std::size_t>(hi + 1), true);
    sieve[0] = false;
    if (hi >= 1) sieve[1] = false;
    for (u64 q = 2; q * q <= hi; ++q) {
        if (!sieve[static_cast<std::size_t>(q)]) continue;
        for (u64 m = q * q; m <= hi; m += q)
            sieve[static_cast<std::size_t>(m)] = false;
    }
    std::vector<u64> out;
    for (u64 p = std::max<u64>(2, lo); p <= hi; ++p)
        if (sieve[static_cast<std::size_t>(p)]) out.push_back(p);
    return out;
}

static std::vector<u64> primes_from_file(const fs::path &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open prime file: " + path.string());
    std::vector<u64> out;
    u64 p;
    while (in >> p) {
        if (!is_prime(p))
            throw std::runtime_error("non-prime in prime file: " + std::to_string(p));
        out.push_back(p);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

struct PrimeRecord {
    u64 p = 0;
    u64 r_limit = 0;
    u64 b_zero_count = 0;
    u64 truncation = 0;
    std::vector<u64> pairs;
    double seconds = 0.0;
};

static fs::path checkpoint_path(const fs::path &dir, u64 p) {
    std::ostringstream name;
    name << "p_";
    name.width(6);
    name.fill('0');
    name << p << ".json";
    return dir / name.str();
}

static void write_checkpoint_atomic(const fs::path &dir, const PrimeRecord &r) {
    fs::create_directories(dir);
    const fs::path final = checkpoint_path(dir, r.p);
    const fs::path temp = final.string() + ".tmp";
    {
        std::ofstream out(temp, std::ios::trunc);
        if (!out) throw std::runtime_error("cannot write checkpoint: " + temp.string());
        out << "{\"version\":\"" << ALGO_VERSION << "\","
            << "\"p\":" << r.p << ','
            << "\"r_limit\":" << r.r_limit << ','
            << "\"b_zero_count\":" << r.b_zero_count << ','
            << "\"truncation\":" << r.truncation << ','
            << "\"pairs\":[";
        for (std::size_t i = 0; i < r.pairs.size(); ++i) {
            if (i) out << ',';
            out << r.pairs[i];
        }
        out << "],\"seconds\":" << r.seconds << "}\n";
        out.flush();
        if (!out) throw std::runtime_error("failed while writing checkpoint");
    }
    if (fs::exists(final)) fs::remove(temp);
    else fs::rename(temp, final);
}

// b_0,...,b_{p-1}.  All denominators (n+1)^3 are units because n+1 < p.
static std::vector<u64> apery_mod_prime(u64 p, std::vector<u64> &inv) {
    const std::size_t N = static_cast<std::size_t>(p);
    inv.assign(N, 0);
    if (p > 1) inv[1] = 1;
    for (u64 n = 2; n < p; ++n) {
        // p is prime and n < p, so p % n != 0 unless n=1.
        const u64 term = mul_mod(p / n, inv[static_cast<std::size_t>(p % n)], p);
        inv[static_cast<std::size_t>(n)] = term == 0 ? 0 : p - term;
    }

    std::vector<u64> b(N, 0);
    b[0] = 1 % p;
    if (p == 1) return b;
    b[1] = 5 % p;
    for (u64 n = 1; n + 1 < p; ++n) {
        const u64 nn = n % p;
        const u64 n2 = mul_mod(nn, nn, p);
        const u64 n3 = mul_mod(n2, nn, p);
        u64 coeff = mul_mod(34 % p, n3, p);
        coeff = add_mod(coeff, mul_mod(51 % p, n2, p), p);
        coeff = add_mod(coeff, mul_mod(27 % p, nn, p), p);
        coeff = add_mod(coeff, 5 % p, p);
        const u64 num = sub_mod(
            mul_mod(coeff, b[static_cast<std::size_t>(n)], p),
            mul_mod(n3, b[static_cast<std::size_t>(n - 1)], p), p);
        const u64 iv = inv[static_cast<std::size_t>(n + 1)];
        const u64 iv3 = cube_mod(iv, p);
        b[static_cast<std::size_t>(n + 1)] = mul_mod(num, iv3, p);
    }
    return b;
}

// q_n for Q(t)=(1-34t+t^2)^(-1/2):
// (n+1) q_{n+1} = (34n+17) q_n - n q_{n-1}.
static std::vector<u64> inv_sqrt_D_coeffs(
    u64 p, std::size_t L, const std::vector<u64> &inv) {
    std::vector<u64> q(L, 0);
    if (L == 0) return q;
    q[0] = 1 % p;
    if (L == 1) return q;
    q[1] = 17 % p;
    for (u64 n = 1; n + 1 < L; ++n) {
        const u64 a = (34 * (n % p) + 17) % p;
        const u64 num = sub_mod(
            mul_mod(a, q[static_cast<std::size_t>(n)], p),
            mul_mod(n % p, q[static_cast<std::size_t>(n - 1)], p), p);
        q[static_cast<std::size_t>(n + 1)] =
            mul_mod(num, inv[static_cast<std::size_t>(n + 1)], p);
    }
    return q;
}

static PrimeRecord scan_prime(u64 p) {
    const auto start = std::chrono::steady_clock::now();
    PrimeRecord rec;
    rec.p = p;
    rec.r_limit = p - 1;

    if (p == 2) {
        // b_1 = 1 mod 2, so there is no r<p with b_r=0.
        rec.seconds = 0.0;
        return rec;
    }

    std::vector<u64> inv;
    std::vector<u64> b = apery_mod_prime(p, inv);
    std::vector<u64> zero_positions;
    for (u64 r = 1; r < p; ++r)
        if (b[static_cast<std::size_t>(r)] == 0) zero_positions.push_back(r);
    rec.b_zero_count = static_cast<u64>(zero_positions.size());

    if (zero_positions.empty()) {
        const auto stop = std::chrono::steady_clock::now();
        rec.seconds = std::chrono::duration<double>(stop - start).count();
        return rec;
    }

    const u64 max_zero = zero_positions.back();
    const std::size_t L = static_cast<std::size_t>(max_zero + 1);
    rec.truncation = static_cast<u64>(L);
    std::vector<u64> q = inv_sqrt_D_coeffs(p, L, inv);

    nmod_poly_t F, Q, Finv, Finv2, G;
    nmod_poly_init2(F, static_cast<ulong>(p), static_cast<slong>(L));
    nmod_poly_init2(Q, static_cast<ulong>(p), static_cast<slong>(L));
    nmod_poly_init2(Finv, static_cast<ulong>(p), static_cast<slong>(L));
    nmod_poly_init2(Finv2, static_cast<ulong>(p), static_cast<slong>(L));
    nmod_poly_init2(G, static_cast<ulong>(p), static_cast<slong>(L));

    for (std::size_t i = 0; i < L; ++i) {
        nmod_poly_set_coeff_ui(F, static_cast<slong>(i), static_cast<ulong>(b[i]));
        nmod_poly_set_coeff_ui(Q, static_cast<slong>(i), static_cast<ulong>(q[i]));
    }

    // F(0)=Q(0)=1, so all series inversions below are unconditionally valid.
    nmod_poly_inv_series(Finv, F, static_cast<slong>(L));
    nmod_poly_mullow(Finv2, Finv, Finv, static_cast<slong>(L));
    nmod_poly_mullow(G, Finv2, Q, static_cast<slong>(L));

    // Xi_0=-1 and Xi_r-Xi_{r-1}=-5 g_r b_{r-1}.
    u64 xi = p - 1;
    std::size_t zindex = 0;
    for (u64 m = 1; m <= max_zero; ++m) {
        const u64 gm = static_cast<u64>(
            nmod_poly_get_coeff_ui(G, static_cast<slong>(m)));
        const u64 delta = mul_mod(5 % p,
            mul_mod(gm, b[static_cast<std::size_t>(m - 1)], p), p);
        xi = sub_mod(xi, delta, p);
        if (zindex < zero_positions.size() && zero_positions[zindex] == m) {
            if (xi == 0) rec.pairs.push_back(m);
            ++zindex;
        }
    }

    nmod_poly_clear(F);
    nmod_poly_clear(Q);
    nmod_poly_clear(Finv);
    nmod_poly_clear(Finv2);
    nmod_poly_clear(G);

    const auto stop = std::chrono::steady_clock::now();
    rec.seconds = std::chrono::duration<double>(stop - start).count();
    return rec;
}

struct Options {
    u64 pmin = 2;
    u64 pmax = 100000;
    fs::path prime_file;
    fs::path checkpoint_dir = "q32_transverse_fast_checkpoints";
    bool force = false;
};

static Options parse_options(int argc, char **argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const char *flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(std::string("missing value after ") + flag);
            return argv[++i];
        };
        if (a == "--pmin") o.pmin = std::stoull(need("--pmin"));
        else if (a == "--pmax") o.pmax = std::stoull(need("--pmax"));
        else if (a == "--prime-file") o.prime_file = need("--prime-file");
        else if (a == "--checkpoint-dir") o.checkpoint_dir = need("--checkpoint-dir");
        else if (a == "--force") o.force = true;
        else if (a == "--help" || a == "-h") {
            std::cout
                << "usage: q32_transverse_common_fast [--pmin N] [--pmax N]\n"
                << "       [--prime-file FILE] [--checkpoint-dir DIR] [--force]\n";
            std::exit(0);
        } else throw std::runtime_error("unknown argument: " + a);
    }
    if (o.pmax > 100000000ULL)
        throw std::runtime_error("refusing pmax>1e8: this scanner is optimized for word primes");
    return o;
}

int main(int argc, char **argv) {
    try {
        const Options opt = parse_options(argc, argv);
        fs::create_directories(opt.checkpoint_dir);
        std::vector<u64> primes = opt.prime_file.empty()
            ? primes_up_to(opt.pmin, opt.pmax)
            : primes_from_file(opt.prime_file);

        std::size_t done = 0;
        for (u64 p : primes) {
            if (p < opt.pmin || p > opt.pmax) continue;
            const fs::path cp = checkpoint_path(opt.checkpoint_dir, p);
            if (!opt.force && fs::exists(cp)) {
                ++done;
                continue;
            }
            PrimeRecord rec = scan_prime(p);
            write_checkpoint_atomic(opt.checkpoint_dir, rec);
            ++done;
            if (!rec.pairs.empty()) {
                for (u64 r : rec.pairs)
                    std::cout << "PAIR p=" << p << " r=" << r
                              << " ratio=" << static_cast<double>(p) / r << "\n";
                std::cout.flush();
            }
            if ((done % 100) == 0) {
                std::cerr << "processed " << done << "/" << primes.size()
                          << " primes; last p=" << p << "\n";
            }
        }
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "fatal: " << e.what() << "\n";
        return 2;
    }
}
