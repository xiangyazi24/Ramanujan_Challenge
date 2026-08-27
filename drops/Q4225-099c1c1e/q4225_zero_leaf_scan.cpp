// Optimized exact Apéry zero and two-row selected-leaf scanner for Q4225.
//
// All prime, zero, raw-window, reflection, and minus-first decisions are exact.
// The scan computes only 0 <= r <= (p-1)/2 and uses the proved Apéry
// reflection b_r == b_(p-1-r) (mod p) to materialize the other half.  An
// independent full divided-recurrence computation checks every p <= 500.

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

namespace {

struct Options {
    int pmax = 200000;
    int threads = 0;
    std::string out_dir = "q4225-scan";
};

Options parse_options(int argc, char **argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto value = [&](const char *name) {
            if (++i >= argc) throw std::runtime_error(std::string("missing value after ") + name);
            return std::string(argv[i]);
        };
        if (arg == "--pmax") opt.pmax = std::stoi(value("--pmax"));
        else if (arg == "--threads") opt.threads = std::stoi(value("--threads"));
        else if (arg == "--out") opt.out_dir = value("--out");
        else if (arg == "--help") {
            std::cout << "usage: q4225_zero_leaf_scan [--pmax N] [--threads N] [--out DIR]\n";
            std::exit(0);
        } else throw std::runtime_error("unknown option: " + arg);
    }
    if (opt.pmax < 5000 || opt.pmax > 500000) {
        throw std::runtime_error("this implementation is certified for 5000 <= pmax <= 500000");
    }
    return opt;
}

std::vector<int> primes_upto(int limit) {
    std::vector<bool> a(limit + 1, true);
    a[0] = a[1] = false;
    for (int p = 2; 1LL * p * p <= limit; ++p) if (a[p]) {
        for (long long n = 1LL * p * p; n <= limit; n += p) a[(int)n] = false;
    }
    std::vector<int> out;
    for (int p = 2; p <= limit; ++p) if (a[p]) out.push_back(p);
    return out;
}

struct FastMod {
    uint32_t p;
    uint64_t reciprocal;
    uint64_t p2;
    explicit FastMod(uint32_t modulus)
        : p(modulus), reciprocal(std::numeric_limits<uint64_t>::max() / modulus),
          p2(uint64_t(modulus) * modulus) {}
    uint32_t reduce(uint64_t x) const {
        uint64_t q = uint64_t((__uint128_t(x) * reciprocal) >> 64);
        uint64_t r = x - q * p;
        while (r >= p) r -= p;
        return uint32_t(r);
    }
};

uint64_t P_exact(uint64_t n) {
    return (2 * n + 1) * (17 * n * n + 17 * n + 5);
}

uint64_t mod_pow(uint64_t a, uint64_t e, uint64_t p) {
    uint64_t r = 1;
    while (e) {
        if (e & 1) r = uint64_t((__uint128_t(r) * a) % p);
        a = uint64_t((__uint128_t(a) * a) % p);
        e >>= 1;
    }
    return r;
}

struct PrimeZeros {
    int p = 0;
    std::vector<int> zeros;
};

PrimeZeros zero_set_reflected_cleared(int p) {
    PrimeZeros out;
    out.p = p;
    if (p <= 5) return out;
    FastMod mod(uint32_t(p));
    uint32_t y0 = 1 % p;
    uint32_t y1 = 5 % p;
    std::vector<int> low;
    if (y0 == 0) low.push_back(0);
    if (y1 == 0) low.push_back(1);
    const int center = (p - 1) / 2;
    for (int n = 1; n < center; ++n) {
        // y_n=(n!)^3 b_n and y_(n+1)=P(n)y_n-n^6 y_(n-1).
        // n <= 250000 here, so n^3 and P(n) fit in uint64_t.
        uint64_t n3 = uint64_t(n) * uint64_t(n) * uint64_t(n);
        uint32_t n3m = mod.reduce(n3);
        uint32_t n6m = mod.reduce(uint64_t(n3m) * n3m);
        uint32_t pm = mod.reduce(P_exact(uint64_t(n)));
        uint64_t expression = uint64_t(pm) * y1 + mod.p2 - uint64_t(n6m) * y0;
        uint32_t y2 = mod.reduce(expression);
        if (y2 == 0) low.push_back(n + 1);
        y0 = y1;
        y1 = y2;
    }
    out.zeros.reserve(low.size() * 2);
    for (int r : low) {
        out.zeros.push_back(r);
        int reflected = p - 1 - r;
        if (reflected != r) out.zeros.push_back(reflected);
    }
    std::sort(out.zeros.begin(), out.zeros.end());
    out.zeros.erase(std::unique(out.zeros.begin(), out.zeros.end()), out.zeros.end());
    return out;
}

std::vector<int> zero_set_full_divided(int p) {
    std::vector<int> out;
    uint64_t b0 = 1 % p, b1 = 5 % p;
    if (b0 == 0) out.push_back(0);
    if (b1 == 0) out.push_back(1);
    for (int n = 1; n <= p - 2; ++n) {
        uint64_t nn = uint64_t(n);
        uint64_t P = P_exact(nn) % p;
        uint64_t n3 = nn * nn % p * nn % p;
        uint64_t rhs = (P * b1 + p - n3 * b0 % p) % p;
        uint64_t den = uint64_t(n + 1) % p;
        den = den * den % p * uint64_t(n + 1) % p;
        uint64_t b2 = uint64_t((__uint128_t(rhs) * mod_pow(den, p - 2, p)) % p);
        if (b2 == 0) out.push_back(n + 1);
        b0 = b1;
        b1 = b2;
    }
    return out;
}

bool contains(const PrimeZeros &d, int r) {
    return std::binary_search(d.zeros.begin(), d.zeros.end(), r);
}

std::string join(const std::vector<int> &v) {
    std::ostringstream os;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) os << ';';
        os << v[i];
    }
    return os.str();
}

struct Leaf {
    char sign;
    int q, t, p, rho, alpha, j6, js;
};

} // namespace

int main(int argc, char **argv) {
    try {
        Options opt = parse_options(argc, argv);
#ifdef _OPENMP
        if (opt.threads > 0) omp_set_num_threads(opt.threads);
        int threads = omp_get_max_threads();
#else
        int threads = 1;
#endif
        fs::create_directories(opt.out_dir);
        auto start = std::chrono::steady_clock::now();
        std::vector<int> primes = primes_upto(opt.pmax);
        std::vector<int> index(opt.pmax + 1, -1);
        for (size_t i = 0; i < primes.size(); ++i) index[primes[i]] = int(i);

        // Independent exact reduction check.
        uint64_t reduction_cases = 0, reduction_failures = 0;
        for (int p : primes) {
            if (p > 10000) break;
            FastMod mod(uint32_t(p));
            std::array<uint64_t, 11> samples = {
                0ULL, 1ULL, uint64_t(p - 1), uint64_t(p), uint64_t(p + 1),
                uint64_t(p) * p - 1, uint64_t(p) * p,
                2ULL * uint64_t(p) * p - 1,
                P_exact(uint64_t((p - 1) / 2)),
                uint64_t((p - 1) / 2) * ((p - 1) / 2) * ((p - 1) / 2),
                std::numeric_limits<uint64_t>::max() / 4
            };
            for (uint64_t x : samples) {
                ++reduction_cases;
                if (mod.reduce(x) != x % p) ++reduction_failures;
            }
        }
        if (reduction_failures) throw std::runtime_error("FastMod correctness gate failed");

        std::vector<PrimeZeros> data(primes.size());
        std::atomic<size_t> completed{0};
#pragma omp parallel for schedule(dynamic, 1) if(threads > 1)
        for (long long i = 0; i < (long long)primes.size(); ++i) {
            data[(size_t)i] = zero_set_reflected_cleared(primes[(size_t)i]);
            size_t done = ++completed;
            if (done % 2000 == 0) {
#pragma omp critical
                std::cerr << "ZERO_PROGRESS " << done << '/' << primes.size() << '\n';
            }
        }

        uint64_t recurrence_checks = 0, recurrence_failures = 0;
        uint64_t reflection_checks = 0, reflection_failures = 0;
        uint64_t consecutive_checks = 0, consecutive_failures = 0;
        uint64_t total_zero_records = 0;
        for (size_t i = 0; i < data.size(); ++i) {
            int p = data[i].p;
            const auto &z = data[i].zeros;
            total_zero_records += z.size();
            if (p >= 7 && p <= 500) {
                ++recurrence_checks;
                std::vector<int> full = zero_set_full_divided(p);
                if (z != full) ++recurrence_failures;
            }
            for (int r : z) {
                ++reflection_checks;
                if (!std::binary_search(z.begin(), z.end(), p - 1 - r)) ++reflection_failures;
            }
            for (size_t j = 1; j < z.size(); ++j) {
                ++consecutive_checks;
                if (z[j] == z[j - 1] + 1) ++consecutive_failures;
            }
        }
        if (recurrence_failures || reflection_failures || consecutive_failures) {
            throw std::runtime_error("zero-set validation gate failed");
        }

        std::unordered_map<long long, std::vector<int>> row6;
        row6.reserve(total_zero_records * 2 + 1);
        for (const auto &d : data) if (d.p >= 7) {
            for (int r : d.zeros) row6[6LL * d.p + r].push_back(d.p);
        }
        for (auto &[row, ps] : row6) {
            std::sort(ps.begin(), ps.end());
            ps.erase(std::unique(ps.begin(), ps.end()), ps.end());
        }

        std::vector<std::pair<int,int>> states;
        std::vector<Leaf> leaves;
        uint64_t selected_le_5000 = 0;
        uint64_t plus_raw_count = 0, minus_raw_count = 0, overlap_count = 0;
        for (const auto &qd : data) {
            int q = qd.p;
            if (q < 17) continue;
            for (int t : qd.zeros) {
                states.emplace_back(q,t);
                if (q <= 5000) ++selected_le_5000;
                int n = 6 * q + t, N = 12 * q + t, s = q - 1 - t;
                auto it = row6.find(n);
                if (it == row6.end()) continue;
                std::set<int> minus, plus_raw;
                for (int p : it->second) {
                    if (p < 17 || p >= q || n >= 7 * p) continue;
                    int pi = index[p];
                    if (pi < 0) continue;
                    int rho = n - 6 * p;
                    int am = s - p;
                    if (0 <= am && am < p && contains(data[(size_t)pi], am)) minus.insert(p);
                    int ap = N - 13 * p;
                    if (0 <= ap && ap < p && contains(data[(size_t)pi], ap)) plus_raw.insert(p);
                }
                minus_raw_count += minus.size();
                plus_raw_count += plus_raw.size();
                for (int p : minus) if (plus_raw.count(p)) ++overlap_count;
                for (int p : minus) {
                    int rho = n - 6 * p, alpha = s - p;
                    leaves.push_back({'-',q,t,p,rho,alpha,
                        std::min(rho,p-1-rho),std::min(alpha,p-1-alpha)});
                }
                for (int p : plus_raw) if (!minus.count(p)) {
                    int rho = n - 6 * p, alpha = N - 13 * p;
                    leaves.push_back({'+',q,t,p,rho,alpha,
                        std::min(rho,p-1-rho),std::min(alpha,p-1-alpha)});
                }
            }
        }
        if (selected_le_5000 != 605) {
            throw std::runtime_error("selected-state checkpoint failed: expected 605 through q=5000");
        }

        {
            std::ofstream out(fs::path(opt.out_dir) / "zero_sets.csv");
            out << "p,zero_count,zeros\n";
            for (const auto &d : data) out << d.p << ',' << d.zeros.size() << ',' << join(d.zeros) << '\n';
        }
        {
            std::ofstream out(fs::path(opt.out_dir) / "selected_states.csv");
            out << "q,t\n";
            for (auto [q,t] : states) out << q << ',' << t << '\n';
        }
        {
            std::ofstream out(fs::path(opt.out_dir) / "raw_leaves.csv");
            out << "sign,q,t,p,rho,alpha,j6,jsigma\n";
            for (const auto &x : leaves) out << x.sign << ',' << x.q << ',' << x.t << ',' << x.p << ','
                << x.rho << ',' << x.alpha << ',' << x.j6 << ',' << x.js << '\n';
        }
        {
            std::ofstream out(fs::path(opt.out_dir) / "witnesses.txt");
            if (leaves.empty()) {
                out << "NO RAW TWO-ROW SELECTED LEAF FOUND THROUGH p,q <= " << opt.pmax << "\n";
                out << "The inherited residual/depth/primitive filters are therefore vacuous on this range.\n";
            } else {
                out << "RAW TWO-ROW LEAF WITNESSES; apply any extra inherited residual mask before DPRM:\n";
                for (const auto &x : leaves) out << x.sign << " q=" << x.q << " t=" << x.t << " p=" << x.p
                    << " rho=" << x.rho << " alpha=" << x.alpha << "\n";
            }
        }

        double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        {
            std::ofstream out(fs::path(opt.out_dir) / "scan_manifest.json");
            size_t plus = std::count_if(leaves.begin(),leaves.end(),[](const Leaf&x){return x.sign=='+';});
            size_t minus = leaves.size() - plus;
            out << "{\n"
                << "  \"pmax\": " << opt.pmax << ",\n"
                << "  \"threads\": " << threads << ",\n"
                << "  \"elapsed_seconds\": " << std::setprecision(12) << elapsed << ",\n"
                << "  \"prime_count\": " << primes.size() << ",\n"
                << "  \"total_zero_records\": " << total_zero_records << ",\n"
                << "  \"selected_states\": " << states.size() << ",\n"
                << "  \"selected_states_le_5000\": " << selected_le_5000 << ",\n"
                << "  \"raw_plus_before_minus_first\": " << plus_raw_count << ",\n"
                << "  \"raw_minus\": " << minus_raw_count << ",\n"
                << "  \"raw_sign_overlap\": " << overlap_count << ",\n"
                << "  \"actual_plus_after_minus_first\": " << plus << ",\n"
                << "  \"actual_minus\": " << minus << ",\n"
                << "  \"checks\": {\n"
                << "    \"reduction_cases\": " << reduction_cases << ", \"reduction_failures\": " << reduction_failures << ",\n"
                << "    \"full_recurrence_crosschecks\": " << recurrence_checks << ", \"recurrence_failures\": " << recurrence_failures << ",\n"
                << "    \"reflection_checks\": " << reflection_checks << ", \"reflection_failures\": " << reflection_failures << ",\n"
                << "    \"consecutive_checks\": " << consecutive_checks << ", \"consecutive_failures\": " << consecutive_failures << ",\n"
                << "    \"q5000_checkpoint_expected\": 605, \"q5000_checkpoint_ok\": true\n"
                << "  }\n"
                << "}\n";
        }
        std::cout << "Q4225_ZERO_LEAF_SCAN_COMPLETE pmax=" << opt.pmax
                  << " primes=" << primes.size() << " zeros=" << total_zero_records
                  << " states=" << states.size() << " leaves=" << leaves.size()
                  << " elapsed=" << elapsed << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Q4225_FATAL " << e.what() << '\n';
        return 2;
    }
}
