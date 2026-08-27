// Q4225 exact canonical popular-edge ledger runner.
//
// This program uses only integer arithmetic for every prime, zero, leaf,
// cell, matching, fold, and divisibility decision.  Floating point is used
// only to turn the user-specified beta into the integer scale parameters and
// to report logarithms.

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>

using boost::multiprecision::cpp_int;
namespace fs = std::filesystem;

namespace {

constexpr int TYPE_DENOMINATOR = 512; // explicit finite-type qualifying constant
constexpr std::array<int, 2> FIXED_SMALL_EXCEPTIONS = {5, 17};

struct Options {
    int pmax = 200000;
    int threads = 0;
    int exact_b_limit = 1000;
    std::string out_dir = "q4225-output";
    std::vector<int> Ts = {1000, 2000, 5000, 10000, 20000, 40000, 80000, 100000};
    std::vector<double> betas = {0.02, 0.04, 0.06, 0.08, 0.09, 0.095, 0.099};
};

std::vector<std::string> split(const std::string &text, char sep) {
    std::vector<std::string> out;
    std::string part;
    std::stringstream ss(text);
    while (std::getline(ss, part, sep)) {
        if (!part.empty()) out.push_back(part);
    }
    return out;
}

Options parse_options(int argc, char **argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto require_value = [&](const char *name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error(std::string("missing value after ") + name);
            }
            return argv[++i];
        };
        if (arg == "--pmax") opt.pmax = std::stoi(require_value("--pmax"));
        else if (arg == "--threads") opt.threads = std::stoi(require_value("--threads"));
        else if (arg == "--exact-b-limit") opt.exact_b_limit = std::stoi(require_value("--exact-b-limit"));
        else if (arg == "--out") opt.out_dir = require_value("--out");
        else if (arg == "--Ts") {
            opt.Ts.clear();
            for (const auto &x : split(require_value("--Ts"), ',')) opt.Ts.push_back(std::stoi(x));
        } else if (arg == "--betas") {
            opt.betas.clear();
            for (const auto &x : split(require_value("--betas"), ',')) opt.betas.push_back(std::stod(x));
        } else if (arg == "--help") {
            std::cout << "usage: q4225_popular_ledger [--pmax N] [--threads N] "
                         "[--Ts a,b,c] [--betas a,b,c] [--out DIR]\n";
            std::exit(0);
        } else {
            throw std::runtime_error("unknown option: " + arg);
        }
    }
    if (opt.pmax < 100) throw std::runtime_error("pmax must be at least 100");
    std::sort(opt.Ts.begin(), opt.Ts.end());
    opt.Ts.erase(std::unique(opt.Ts.begin(), opt.Ts.end()), opt.Ts.end());
    std::sort(opt.betas.begin(), opt.betas.end());
    return opt;
}

std::vector<int> primes_upto(int limit) {
    std::vector<bool> is_prime(limit + 1, true);
    if (limit >= 0) is_prime[0] = false;
    if (limit >= 1) is_prime[1] = false;
    for (int p = 2; 1LL * p * p <= limit; ++p) {
        if (!is_prime[p]) continue;
        for (long long n = 1LL * p * p; n <= limit; n += p) is_prime[(int)n] = false;
    }
    std::vector<int> primes;
    for (int p = 2; p <= limit; ++p) if (is_prime[p]) primes.push_back(p);
    return primes;
}

uint64_t poly_P_exact(uint64_t n) {
    return (2 * n + 1) * (17 * n * n + 17 * n + 5);
}

uint64_t pow6_exact(uint64_t n) {
    uint64_t n2 = n * n;
    uint64_t n3 = n2 * n;
    return n3 * n3;
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

template <size_t N, class F>
std::array<uint32_t, N> initial_forward_differences(uint32_t p, F value) {
    std::array<uint64_t, N> row{};
    std::array<uint32_t, N> out{};
    for (size_t i = 0; i < N; ++i) row[i] = value(uint64_t(i + 1)) % p;
    for (size_t order = 0; order < N; ++order) {
        out[order] = uint32_t(row[0] % p);
        for (size_t i = 0; i + 1 < N - order; ++i) {
            row[i] = (row[i + 1] + p - row[i]) % p;
        }
    }
    return out;
}

template <size_t N>
inline void advance_differences(std::array<uint32_t, N> &d, uint32_t p) {
    for (size_t i = 0; i + 1 < N; ++i) {
        uint32_t x = d[i] + d[i + 1];
        if (x >= p) x -= p;
        d[i] = x;
    }
}

struct PrimeZeros {
    int p = 0;
    std::vector<int> zeros;
};

PrimeZeros zero_set_cleared(int p) {
    PrimeZeros result;
    result.p = p;
    if (p <= 5) return result;

    FastMod mod(uint32_t(p));
    auto pd = initial_forward_differences<4>(uint32_t(p), poly_P_exact);
    auto n6d = initial_forward_differences<7>(uint32_t(p), pow6_exact);

    uint32_t y_prev = 1 % p;
    uint32_t y = 5 % p;
    if (y_prev == 0) result.zeros.push_back(0);
    if (y == 0) result.zeros.push_back(1);

    // y_n=(n!)^3 b_n and y_{n+1}=P(n)y_n-n^6 y_{n-1}.
    // For n<p, n! is a unit, so y_n and b_n have exactly the same zeros.
    for (int n = 1; n <= p - 2; ++n) {
        uint64_t x = uint64_t(pd[0]) * y + mod.p2 - uint64_t(n6d[0]) * y_prev;
        uint32_t y_next = mod.reduce(x);
        if (y_next == 0) result.zeros.push_back(n + 1);
        y_prev = y;
        y = y_next;
        advance_differences(pd, uint32_t(p));
        advance_differences(n6d, uint32_t(p));
    }
    return result;
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

std::vector<int> zero_set_divided(int p) {
    std::vector<int> zeros;
    uint64_t b0 = 1, b1 = 5 % p;
    if (b0 == 0) zeros.push_back(0);
    if (b1 == 0) zeros.push_back(1);
    for (int n = 1; n <= p - 2; ++n) {
        uint64_t nn = uint64_t(n) % p;
        uint64_t P = poly_P_exact(nn) % p;
        uint64_t rhs = (P * b1 + p - uint64_t((__uint128_t(nn * nn % p) * nn) % p) * b0 % p) % p;
        uint64_t den = uint64_t(n + 1) % p;
        den = uint64_t((__uint128_t(den) * den) % p);
        den = uint64_t((__uint128_t(den) * (n + 1)) % p);
        uint64_t b2 = uint64_t((__uint128_t(rhs) * mod_pow(den, p - 2, p)) % p);
        if (b2 == 0) zeros.push_back(n + 1);
        b0 = b1;
        b1 = b2;
    }
    return zeros;
}

bool contains_zero(const PrimeZeros &data, int r) {
    return std::binary_search(data.zeros.begin(), data.zeros.end(), r);
}

struct Checks {
    uint64_t recurrence_crosschecks = 0;
    uint64_t recurrence_failures = 0;
    uint64_t reflection_checks = 0;
    uint64_t reflection_failures = 0;
    uint64_t consecutive_checks = 0;
    uint64_t consecutive_failures = 0;
    uint64_t barrett_checks = 0;
    uint64_t barrett_failures = 0;
    uint64_t selected_states_le_5000 = 0;
    uint64_t expected_selected_states_le_5000 = 605;
    bool selected_state_checkpoint_ok = false;
};

void check_barrett(const std::vector<int> &primes, Checks &checks) {
    for (int p : primes) {
        if (p < 7 || p > 10000) continue;
        FastMod fm(uint32_t(p));
        const std::array<uint64_t, 9> seeds = {
            0ULL, 1ULL, uint64_t(p - 1), uint64_t(p), uint64_t(p) + 1,
            uint64_t(p) * (p - 1), uint64_t(p) * p - 1,
            uint64_t(p) * p, 2ULL * uint64_t(p) * p - 1
        };
        for (uint64_t x : seeds) {
            ++checks.barrett_checks;
            if (fm.reduce(x) != x % p) ++checks.barrett_failures;
        }
    }
}

struct Leaf {
    char sign = '?'; // '+' or '-'
    int q = 0;
    int t = 0;
    int p = 0;
    int rho = 0;
    int alpha = 0;
    int j6 = 0;
    int js = 0;
};

int folded_depth(int M, int c, int p) {
    int r = M - c * p;
    return std::min(r, p - 1 - r);
}

struct VertexType {
    int s6 = 0;
    int ss = 0;
    int order = 0; // -1 means j6<js, +1 means j6>js
    bool valid = false;

    auto key() const { return std::tuple<int,int,int>(s6, ss, order); }
};

VertexType vertex_type(const Leaf &leaf) {
    const int c = leaf.sign == '+' ? 13 : 1;
    auto slope = [](int r, int p, int c) -> std::optional<int> {
        long long cmp = 2LL * r - (p - 1);
        if (cmp < 0) return -c;
        if (cmp > 0) return c + 1;
        return std::nullopt; // midpoint transition
    };
    auto a = slope(leaf.rho, leaf.p, 6);
    auto b = slope(leaf.alpha, leaf.p, c);
    if (!a || !b || leaf.j6 == leaf.js) return {};
    VertexType result;
    result.s6 = *a;
    result.ss = *b;
    result.order = leaf.j6 < leaf.js ? -1 : 1;
    result.valid = true;
    return result;
}

struct EdgeType {
    int s6 = 0;
    int ss = 0;
    int order = 0;
    int lambda = 0;
    bool stable = false;

    auto key() const { return std::tuple<int,int,int>(s6, ss, order); }
    std::string text() const {
        if (!stable) return "transition";
        std::ostringstream os;
        os << "s6=" << s6 << ";ss=" << ss << ";order=" << (order < 0 ? "6<sigma" : "sigma<6")
           << ";lambda=" << lambda;
        return os.str();
    }
};

EdgeType edge_type(const Leaf &left, const Leaf &right) {
    VertexType a = vertex_type(left), b = vertex_type(right);
    if (!a.valid || !b.valid || a.key() != b.key()) return {};
    EdgeType result;
    result.s6 = a.s6;
    result.ss = a.ss;
    result.order = a.order;
    result.lambda = std::abs(a.ss - a.s6);
    result.stable = true;
    return result;
}

struct MatchEdge {
    const Leaf *left = nullptr;
    const Leaf *right = nullptr;
    EdgeType type;
};

struct CanonicalResult {
    bool qualifies = false;
    int delta = 0;
    int R = 0;
    int qualifying_min = 0;
    int delta_max = 0;
    EdgeType type;
    std::vector<MatchEdge> edges;
};

CanonicalResult canonical_matching(const std::vector<const Leaf*> &deep, int Lambda) {
    CanonicalResult best;
    const int K = int(deep.size());
    if (K < 2 || Lambda <= 0) return best;
    best.qualifying_min = std::max(1, int((1LL * K * K + 1LL * TYPE_DENOMINATOR * Lambda - 1) /
                                         (1LL * TYPE_DENOMINATOR * Lambda)));
    best.delta_max = std::max(1, int((4LL * Lambda + K - 1) / K));

    std::set<int> gaps;
    for (int i = 0; i + 1 < K; ++i) {
        int d = deep[i + 1]->p - deep[i]->p;
        if (d > 0 && d <= best.delta_max) gaps.insert(d);
    }

    for (int d : gaps) {
        std::vector<MatchEdge> greedy;
        int i = 0;
        while (i + 1 < K) {
            if (deep[i + 1]->p - deep[i]->p == d) {
                MatchEdge edge;
                edge.left = deep[i];
                edge.right = deep[i + 1];
                edge.type = edge_type(*deep[i], *deep[i + 1]);
                greedy.push_back(edge);
                i += 2;
            } else {
                ++i;
            }
        }

        std::map<std::tuple<int,int,int>, std::vector<MatchEdge>> classes;
        for (const auto &edge : greedy) if (edge.type.stable) classes[edge.type.key()].push_back(edge);
        if (classes.empty()) continue;

        auto chosen = classes.begin();
        for (auto it = classes.begin(); it != classes.end(); ++it) {
            if (it->second.size() > chosen->second.size()) chosen = it;
        }
        int R = int(chosen->second.size());
        if (R < best.qualifying_min) continue;

        best.qualifies = true;
        best.delta = d;
        best.R = R;
        best.edges = chosen->second;
        best.type = chosen->second.front().type;
        return best; // lexicographically least qualifying gap
    }
    return best;
}

std::vector<cpp_int> exact_apery_row(int limit) {
    std::vector<cpp_int> b;
    b.reserve(limit + 1);
    b.push_back(1);
    if (limit == 0) return b;
    b.push_back(5);
    for (int n = 1; n < limit; ++n) {
        cpp_int P = cpp_int(34) * n * n * n + cpp_int(51) * n * n + cpp_int(27) * n + 5;
        cpp_int num = P * b[n] - cpp_int(n) * n * n * b[n - 1];
        cpp_int den = cpp_int(n + 1) * (n + 1) * (n + 1);
        if (num % den != 0) throw std::runtime_error("nonintegral exact Apéry recurrence");
        b.push_back(num / den);
    }
    return b;
}

cpp_int gcd_cpp(cpp_int a, cpp_int b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        cpp_int r = a % b;
        a = b;
        b = r;
    }
    return a;
}

void strip_cpp(cpp_int &x, int p) {
    while (x != 0 && x % p == 0) x /= p;
}

std::string join_ints(const std::vector<int> &values, char sep=';') {
    std::ostringstream os;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) os << sep;
        os << values[i];
    }
    return os.str();
}

std::string beta_text(double beta) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(3) << beta;
    return os.str();
}

long long floor_power(long long base, double exponent) {
    long double x = std::pow((long double)base, (long double)exponent);
    long long candidate = std::max<long long>(1, (long long)std::floor(x + 1e-12L));
    // Correct rare floating boundary slips by comparison in logarithms.
    while (std::log((long double)(candidate + 1)) <= exponent * std::log((long double)base) + 1e-15L) ++candidate;
    while (candidate > 1 && std::log((long double)candidate) > exponent * std::log((long double)base) + 1e-15L) --candidate;
    return candidate;
}

struct DprmRecord {
    int T = 0;
    double beta = 0;
    int Y = 0, Lambda = 0, H = 0, Jstar = 0, tau = 0;
    char sign = '?';
    int q = 0, t = 0, cell = 0;
    int Kraw = 0, Kdeep = 0;
    int delta = 0, R = 0, qualifying_min = 0, delta_max = 0;
    std::string type;
    double R_normalized = 0;
    int content_bits = 0;
    std::string content_value;
};

struct LedgerRow {
    int T = 0;
    double beta = 0;
    int Y = 0, Lambda = 0, H = 0, Jstar = 0;
    char sign = '?';
    long long M = 0;
    double rich_mass = 0;
    double deep_rich_mass = 0;
    double dprm_mass = 0;
    double Epop = 0;
    long long rich_cells = 0;
    long long dprm_cells = 0;
};

uint64_t state_key(int q, int t) {
    return (uint64_t(uint32_t(q)) << 32) | uint32_t(t);
}

struct StateLeaves {
    int q = 0, t = 0;
    std::vector<const Leaf*> plus;
    std::vector<const Leaf*> minus;
};

struct HistogramKey {
    int T;
    std::string beta;
    char sign;
    std::string key;
    bool operator<(const HistogramKey &o) const {
        return std::tie(T,beta,sign,key) < std::tie(o.T,o.beta,o.sign,o.key);
    }
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
        auto started = std::chrono::steady_clock::now();

        std::vector<int> primes = primes_upto(opt.pmax);
        std::vector<int> prime_index(opt.pmax + 1, -1);
        for (size_t i = 0; i < primes.size(); ++i) prime_index[primes[i]] = int(i);

        Checks checks;
        check_barrett(primes, checks);
        if (checks.barrett_failures) throw std::runtime_error("Barrett reduction self-check failed");

        std::vector<PrimeZeros> zero_data(primes.size());
        std::atomic<size_t> done{0};
#pragma omp parallel for schedule(dynamic, 1) if(threads > 1)
        for (long long i = 0; i < (long long)primes.size(); ++i) {
            zero_data[(size_t)i] = zero_set_cleared(primes[(size_t)i]);
            size_t d = ++done;
            if (d % 2000 == 0) {
#pragma omp critical
                std::cerr << "ZERO_PROGRESS " << d << "/" << primes.size() << "\n";
            }
        }

        for (size_t i = 0; i < primes.size(); ++i) {
            int p = primes[i];
            if (p < 7) continue;
            const auto &z = zero_data[i].zeros;
            if (p <= 500) {
                ++checks.recurrence_crosschecks;
                if (z != zero_set_divided(p)) ++checks.recurrence_failures;
            }
            for (int r : z) {
                ++checks.reflection_checks;
                if (!std::binary_search(z.begin(), z.end(), p - 1 - r)) ++checks.reflection_failures;
            }
            for (size_t j = 1; j < z.size(); ++j) {
                ++checks.consecutive_checks;
                if (z[j] == z[j - 1] + 1) ++checks.consecutive_failures;
            }
        }
        if (checks.recurrence_failures || checks.reflection_failures || checks.consecutive_failures) {
            throw std::runtime_error("zero-set correctness gate failed");
        }

        std::unordered_map<long long, std::vector<int>> row6;
        size_t total_zero_records = 0;
        for (const auto &d : zero_data) total_zero_records += d.zeros.size();
        row6.reserve(total_zero_records * 2 + 1);
        for (const auto &d : zero_data) {
            if (d.p < 7) continue;
            for (int r : d.zeros) row6[6LL * d.p + r].push_back(d.p);
        }
        for (auto &[row, ps] : row6) {
            std::sort(ps.begin(), ps.end());
            ps.erase(std::unique(ps.begin(), ps.end()), ps.end());
        }

        std::vector<Leaf> leaves;
        std::vector<std::pair<int,int>> selected_states;
        uint64_t plus_raw_count = 0, minus_count = 0, overlap_count = 0;
        for (const auto &qd : zero_data) {
            int q = qd.p;
            if (q < 17) continue;
            for (int t : qd.zeros) {
                selected_states.emplace_back(q,t);
                if (q <= 5000) ++checks.selected_states_le_5000;
                int n = 6 * q + t;
                int N = 12 * q + t;
                int s = q - 1 - t;
                auto hit = row6.find(n);
                if (hit == row6.end()) continue;
                std::set<int> minus, plus_raw;
                for (int p : hit->second) {
                    if (p < 17 || p >= q || n >= 7 * p) continue;
                    int idx = p <= opt.pmax ? prime_index[p] : -1;
                    if (idx < 0) continue;
                    int rho = n - 6 * p;
                    int am = s - p;
                    if (0 <= am && am < p && contains_zero(zero_data[(size_t)idx], am)) minus.insert(p);
                    int ap = N - 13 * p;
                    if (0 <= ap && ap < p && contains_zero(zero_data[(size_t)idx], ap)) plus_raw.insert(p);
                }
                minus_count += minus.size();
                plus_raw_count += plus_raw.size();
                for (int p : minus) if (plus_raw.count(p)) ++overlap_count;
                for (int p : minus) {
                    int rho = n - 6 * p, alpha = s - p;
                    leaves.push_back(Leaf{'-',q,t,p,rho,alpha,
                        std::min(rho,p-1-rho), std::min(alpha,p-1-alpha)});
                }
                for (int p : plus_raw) if (!minus.count(p)) {
                    int rho = n - 6 * p, alpha = N - 13 * p;
                    leaves.push_back(Leaf{'+',q,t,p,rho,alpha,
                        std::min(rho,p-1-rho), std::min(alpha,p-1-alpha)});
                }
            }
        }
        checks.selected_state_checkpoint_ok =
            opt.pmax < 5000 || checks.selected_states_le_5000 == checks.expected_selected_states_le_5000;
        if (!checks.selected_state_checkpoint_ok) throw std::runtime_error("q<=5000 selected-state checkpoint failed");

        // Build exact state index.
        std::map<std::pair<int,int>, StateLeaves> states;
        for (const Leaf &leaf : leaves) {
            auto &state = states[{leaf.q, leaf.t}];
            state.q = leaf.q; state.t = leaf.t;
            (leaf.sign == '+' ? state.plus : state.minus).push_back(&leaf);
        }
        for (auto &[key, state] : states) {
            auto cmp = [](const Leaf *a, const Leaf *b){ return a->p < b->p; };
            std::sort(state.plus.begin(), state.plus.end(), cmp);
            std::sort(state.minus.begin(), state.minus.end(), cmp);
        }

        // Exact Apéry row used only for small-index primitive-content diagnostics.
        std::vector<cpp_int> exact_b = exact_apery_row(opt.exact_b_limit);

        std::vector<LedgerRow> ledger_rows;
        std::vector<DprmRecord> dprm_records;
        std::map<HistogramKey,long long> gap_hist, R_hist, type_hist, residue_hist;

        for (int T : opt.Ts) {
            if (2LL * T > opt.pmax + 1LL) continue;
            for (double beta : opt.betas) {
                int Y = int(floor_power(T, beta));
                int Lambda = int(floor_power(Y, 5.0/6.0));
                int H = int(floor_power(Y, 1.0/2.0));
                int Jstar = int(floor_power(Y, 1.0/6.0));
                std::string btxt = beta_text(beta);

                for (char sign : std::array<char,2>{'+','-'}) {
                    long long M = 0;
                    for (const auto &[key, state] : states) {
                        if (!(T <= state.q && state.q < 2*T)) continue;
                        const auto &v = sign == '+' ? state.plus : state.minus;
                        for (const Leaf *leaf : v) if (T <= leaf->p && leaf->p < 2*T) ++M;
                    }

                    double rich_mass_sum = 0, deep_mass_sum = 0, dprm_mass_sum = 0, epop_sum = 0;
                    long long rich_cells_sum = 0, dprm_cells_sum = 0;

                    for (int tau = 0; tau < Lambda; ++tau) {
                        for (const auto &[key, state] : states) {
                            if (!(T <= state.q && state.q < 2*T)) continue;
                            const auto &source = sign == '+' ? state.plus : state.minus;
                            std::map<int,std::vector<const Leaf*>> cells;
                            for (const Leaf *leaf : source) {
                                if (!(T <= leaf->p && leaf->p < 2*T)) continue;
                                int cell = (leaf->p - tau) / Lambda;
                                cells[cell].push_back(leaf);
                            }
                            for (auto &[cell_index, cell] : cells) {
                                std::sort(cell.begin(), cell.end(), [](auto *a, auto *b){return a->p<b->p;});
                                int Kraw = int(cell.size());
                                if (Kraw <= H) continue;
                                rich_mass_sum += Kraw;
                                ++rich_cells_sum;
                                std::vector<const Leaf*> deep;
                                for (const Leaf *leaf : cell) if (leaf->j6 > Jstar && leaf->js > Jstar) deep.push_back(leaf);
                                if ((int)deep.size() <= H) continue;
                                deep_mass_sum += deep.size();
                                CanonicalResult canon = canonical_matching(deep, Lambda);
                                if (!canon.qualifies) continue;
                                dprm_mass_sum += deep.size();
                                epop_sum += canon.R;
                                ++dprm_cells_sum;

                                DprmRecord rec;
                                rec.T=T; rec.beta=beta; rec.Y=Y; rec.Lambda=Lambda; rec.H=H; rec.Jstar=Jstar;
                                rec.tau=tau; rec.sign=sign; rec.q=state.q; rec.t=state.t; rec.cell=cell_index;
                                rec.Kraw=Kraw; rec.Kdeep=int(deep.size()); rec.delta=canon.delta; rec.R=canon.R;
                                rec.qualifying_min=canon.qualifying_min; rec.delta_max=canon.delta_max;
                                rec.type=canon.type.text(); rec.R_normalized=double(canon.R)*Lambda/(double(deep.size())*deep.size());

                                cpp_int common = 0;
                                bool content_available = !canon.edges.empty();
                                for (const auto &edge : canon.edges) {
                                    const Leaf &a=*edge.left, &b=*edge.right;
                                    int maxidx=std::max({a.rho,a.alpha,b.rho,b.alpha});
                                    if (maxidx>opt.exact_b_limit) {content_available=false; break;}
                                    cpp_int carrier=exact_b[a.rho]*exact_b[b.alpha]+exact_b[a.alpha]*exact_b[b.rho];
                                    common = common==0 ? carrier : gcd_cpp(common,carrier);
                                }
                                if (content_available && common!=0) {
                                    strip_cpp(common,state.q); strip_cpp(common,5); strip_cpp(common,17);
                                    rec.content_bits = common==0 ? 0 : int(boost::multiprecision::msb(common)+1);
                                    if (rec.content_bits <= 1024) rec.content_value=common.convert_to<std::string>();
                                    else rec.content_value="too-large";
                                } else rec.content_value="not-evaluated";
                                dprm_records.push_back(rec);

                                gap_hist[{T,btxt,sign,std::to_string(canon.delta)}]++;
                                R_hist[{T,btxt,sign,std::to_string(canon.R)}]++;
                                type_hist[{T,btxt,sign,canon.type.text()}]++;
                                for (const auto &edge : canon.edges) {
                                    for (int mod : {3,4,5,7,8,9,11,13}) {
                                        std::ostringstream keytxt;
                                        keytxt << "mod" << mod << "=" << edge.left->p%mod;
                                        residue_hist[{T,btxt,sign,keytxt.str()}]++;
                                    }
                                }
                            }
                        }
                    }

                    LedgerRow row;
                    row.T=T; row.beta=beta; row.Y=Y; row.Lambda=Lambda; row.H=H; row.Jstar=Jstar; row.sign=sign; row.M=M;
                    row.rich_mass=rich_mass_sum/Lambda; row.deep_rich_mass=deep_mass_sum/Lambda;
                    row.dprm_mass=dprm_mass_sum/Lambda; row.Epop=epop_sum/Lambda;
                    row.rich_cells=rich_cells_sum; row.dprm_cells=dprm_cells_sum;
                    ledger_rows.push_back(row);
                }
            }
        }

        // Raw zero sets.
        {
            std::ofstream out(fs::path(opt.out_dir)/"zero_sets.csv");
            out << "p,zero_count,zeros\n";
            for (const auto &d : zero_data) out << d.p << ',' << d.zeros.size() << ',' << join_ints(d.zeros) << '\n';
        }
        {
            std::ofstream out(fs::path(opt.out_dir)/"selected_states.csv");
            out << "q,t\n";
            for (auto [q,t] : selected_states) out << q << ',' << t << '\n';
        }
        {
            std::ofstream out(fs::path(opt.out_dir)/"raw_leaves.csv");
            out << "sign,q,t,p,rho,alpha,j6,jsigma\n";
            for (const auto &x : leaves) out << x.sign << ',' << x.q << ',' << x.t << ',' << x.p << ','
                << x.rho << ',' << x.alpha << ',' << x.j6 << ',' << x.js << '\n';
        }
        {
            std::ofstream out(fs::path(opt.out_dir)/"ledger_summary.csv");
            out << "T,beta,Y,Lambda,H,Jstar,sign,M,rich_mass,deep_rich_mass,dprm_mass,Epop,Epop_over_M,rich_cells_raw_sum,dprm_cells_raw_sum,dprm_fraction_of_rich\n";
            out << std::setprecision(17);
            for (const auto &r : ledger_rows) {
                out << r.T << ',' << r.beta << ',' << r.Y << ',' << r.Lambda << ',' << r.H << ',' << r.Jstar << ','
                    << r.sign << ',' << r.M << ',' << r.rich_mass << ',' << r.deep_rich_mass << ',' << r.dprm_mass << ',' << r.Epop << ',';
                if (r.M) out << r.Epop/r.M;
                else out << "NA";
                out << ',' << r.rich_cells << ',' << r.dprm_cells << ',';
                if (r.rich_mass) out << r.dprm_mass/r.rich_mass;
                else out << "NA";
                out << '\n';
            }
        }
        {
            std::ofstream out(fs::path(opt.out_dir)/"dprm_cells.csv");
            out << "T,beta,Y,Lambda,H,Jstar,tau,sign,q,t,cell,Kraw,Kdeep,delta,R,qualifying_min,delta_max,type,R_Lambda_over_K2,primitive_content_bits,primitive_content\n";
            out << std::setprecision(17);
            for (const auto &r : dprm_records) {
                out << r.T << ',' << r.beta << ',' << r.Y << ',' << r.Lambda << ',' << r.H << ',' << r.Jstar << ','
                    << r.tau << ',' << r.sign << ',' << r.q << ',' << r.t << ',' << r.cell << ',' << r.Kraw << ',' << r.Kdeep << ','
                    << r.delta << ',' << r.R << ',' << r.qualifying_min << ',' << r.delta_max << ',"' << r.type << '",' << r.R_normalized
                    << ',' << r.content_bits << ',"' << r.content_value << '"' << '\n';
            }
        }
        auto write_hist = [&](const char *name, const std::map<HistogramKey,long long> &hist) {
            std::ofstream out(fs::path(opt.out_dir)/name);
            out << "T,beta,sign,key,count\n";
            for (const auto &[k,v] : hist) out << k.T << ',' << k.beta << ',' << k.sign << ',"' << k.key << '",' << v << '\n';
        };
        write_hist("gap_hist.csv",gap_hist);
        write_hist("R_hist.csv",R_hist);
        write_hist("type_hist.csv",type_hist);
        write_hist("residue_hist.csv",residue_hist);

        auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        {
            std::ofstream out(fs::path(opt.out_dir)/"manifest.json");
            out << "{\n";
            out << "  \"pmax\": " << opt.pmax << ",\n";
            out << "  \"threads\": " << threads << ",\n";
            out << "  \"elapsed_seconds\": " << std::setprecision(12) << elapsed << ",\n";
            out << "  \"prime_count\": " << primes.size() << ",\n";
            out << "  \"total_zero_records\": " << total_zero_records << ",\n";
            out << "  \"selected_states\": " << selected_states.size() << ",\n";
            out << "  \"selected_states_le_5000\": " << checks.selected_states_le_5000 << ",\n";
            out << "  \"raw_plus_before_minus_first\": " << plus_raw_count << ",\n";
            out << "  \"raw_minus\": " << minus_count << ",\n";
            out << "  \"raw_sign_overlap\": " << overlap_count << ",\n";
            out << "  \"actual_plus_after_minus_first\": " << std::count_if(leaves.begin(),leaves.end(),[](const Leaf&x){return x.sign=='+';} ) << ",\n";
            out << "  \"actual_minus\": " << std::count_if(leaves.begin(),leaves.end(),[](const Leaf&x){return x.sign=='-';} ) << ",\n";
            out << "  \"dprm_cell_records\": " << dprm_records.size() << ",\n";
            out << "  \"scale_rounding\": \"Y=floor(T^beta), Lambda=floor(Y^(5/6)), H=floor(Y^(1/2)), Jstar=floor(Y^(1/6)), each clamped to >=1\",\n";
            out << "  \"shell_convention\": \"T <= p,q < 2T\",\n";
            out << "  \"canonical_type_denominator\": " << TYPE_DENOMINATOR << ",\n";
            out << "  \"checks\": {\n";
            out << "    \"barrett_cases\": " << checks.barrett_checks << ", \"barrett_failures\": " << checks.barrett_failures << ",\n";
            out << "    \"recurrence_crosschecks\": " << checks.recurrence_crosschecks << ", \"recurrence_failures\": " << checks.recurrence_failures << ",\n";
            out << "    \"reflection_checks\": " << checks.reflection_checks << ", \"reflection_failures\": " << checks.reflection_failures << ",\n";
            out << "    \"consecutive_checks\": " << checks.consecutive_checks << ", \"consecutive_failures\": " << checks.consecutive_failures << ",\n";
            out << "    \"q5000_checkpoint_expected\": 605, \"q5000_checkpoint_ok\": " << (checks.selected_state_checkpoint_ok?"true":"false") << "\n";
            out << "  }\n";
            out << "}\n";
        }
        {
            std::ofstream out(fs::path(opt.out_dir)/"witnesses.txt");
            if (leaves.empty()) {
                out << "NO RAW TWO-ROW SELECTED LEAF FOUND THROUGH p,q <= " << opt.pmax << "\n";
                out << "Consequently every inherited residual/depth/product filter is vacuous on this range.\n";
            } else {
                out << "RAW LEAF WITNESSES (must be checked against any additional inherited residual mask):\n";
                for (const auto &x : leaves) out << x.sign << " q=" << x.q << " t=" << x.t << " p=" << x.p
                    << " rho=" << x.rho << " alpha=" << x.alpha << "\n";
            }
        }

        std::cout << "Q4225_SCAN_COMPLETE pmax=" << opt.pmax << " primes=" << primes.size()
                  << " zeros=" << total_zero_records << " states=" << selected_states.size()
                  << " plus=" << std::count_if(leaves.begin(),leaves.end(),[](const Leaf&x){return x.sign=='+';})
                  << " minus=" << std::count_if(leaves.begin(),leaves.end(),[](const Leaf&x){return x.sign=='-';})
                  << " dprm_cells=" << dprm_records.size() << " elapsed=" << elapsed << "\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Q4225_FATAL " << e.what() << "\n";
        return 2;
    }
}
