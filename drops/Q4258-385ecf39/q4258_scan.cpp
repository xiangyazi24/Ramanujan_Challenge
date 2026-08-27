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
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_set>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;
namespace fs = std::filesystem;

struct FastMod {
    u64 d;
    u64 reciprocal;
    explicit FastMod(u32 modulus)
        : d(modulus), reciprocal(std::numeric_limits<u64>::max() / modulus) {}
    u32 reduce(u64 x) const {
        const u64 q = static_cast<u64>((static_cast<u128>(x) * reciprocal) >> 64);
        u64 r = x - q * d;
        if (r >= d) r -= d;
        if (r >= d) r %= d;
        return static_cast<u32>(r);
    }
    u32 mul(u32 a, u32 b) const { return reduce(static_cast<u64>(a) * b); }
};

static u32 add_mod(u32 a, u32 b, u32 p) {
    const u64 c = static_cast<u64>(a) + b;
    return static_cast<u32>(c >= p ? c - p : c);
}
static u32 sub_mod(u32 a, u32 b, u32 p) {
    return a >= b ? a - b : static_cast<u32>(static_cast<u64>(a) + p - b);
}
template <std::size_t N>
static void advance_differences(std::array<u32, N>& d, u32 p) {
    for (std::size_t k = 0; k + 1 < N; ++k) d[k] = add_mod(d[k], d[k + 1], p);
}

static std::vector<u32> prime_sieve(u32 n, std::vector<std::uint8_t>& is_prime) {
    is_prime.assign(static_cast<std::size_t>(n) + 1, 1);
    is_prime[0] = 0;
    if (n >= 1) is_prime[1] = 0;
    for (u64 i = 2; i * i <= n; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) continue;
        for (u64 j = i * i; j <= n; j += i) is_prime[static_cast<std::size_t>(j)] = 0;
    }
    std::vector<u32> primes;
    for (u32 i = 2; i <= n; ++i) if (is_prime[i]) primes.push_back(i);
    return primes;
}

static std::vector<u32> apery_zeros(u32 p, u64& recurrence_steps) {
    if (p < 7) return {};
    FastMod mod{p};
    std::array<u32, 4> dp{
        static_cast<u32>(117 % p), static_cast<u32>(418 % p),
        static_cast<u32>(510 % p), static_cast<u32>(204 % p)
    };
    std::array<u32, 7> d6{
        static_cast<u32>(1 % p), static_cast<u32>(63 % p),
        static_cast<u32>(602 % p), static_cast<u32>(2100 % p),
        static_cast<u32>(3360 % p), static_cast<u32>(2520 % p),
        static_cast<u32>(720 % p)
    };
    const u32 midpoint = (p - 1) / 2;
    u32 previous = 1 % p;
    u32 current = 5 % p;
    std::vector<u32> half;
    if (previous == 0) half.push_back(0);
    if (midpoint >= 1 && current == 0) half.push_back(1);
    for (u32 n = 1; n < midpoint; ++n) {
        const u32 next = sub_mod(mod.mul(dp[0], current), mod.mul(d6[0], previous), p);
        ++recurrence_steps;
        if (next == 0) half.push_back(n + 1);
        previous = current;
        current = next;
        advance_differences(dp, p);
        advance_differences(d6, p);
    }
    std::vector<u32> zeros;
    zeros.reserve(half.size() * 2);
    for (u32 z : half) {
        zeros.push_back(z);
        const u32 reflected = p - 1 - z;
        if (reflected != z) zeros.push_back(reflected);
    }
    std::sort(zeros.begin(), zeros.end());
    return zeros;
}

static u64 apery_small(unsigned target) {
    std::vector<u128> b(target + 1);
    b[0] = 1;
    if (target >= 1) b[1] = 5;
    for (unsigned n = 1; n < target; ++n) {
        const u128 P = 34u * static_cast<u128>(n) * n * n
                     + 51u * static_cast<u128>(n) * n
                     + 27u * static_cast<u128>(n) + 5u;
        const u128 numerator = P * b[n] - static_cast<u128>(n) * n * n * b[n - 1];
        const u128 denominator = static_cast<u128>(n + 1) * (n + 1) * (n + 1);
        if (numerator % denominator != 0) throw std::runtime_error("nonintegral small Apéry recurrence");
        b[n + 1] = numerator / denominator;
    }
    if (b[target] > std::numeric_limits<u64>::max()) throw std::runtime_error("small Apéry value overflow");
    return static_cast<u64>(b[target]);
}

enum class Type : std::uint8_t { Plus = 0, MPlus = 1, MMinus = 2 };
static const char* type_name(Type t) {
    switch (t) {
        case Type::Plus: return "P";
        case Type::MPlus: return "M+";
        case Type::MMinus: return "M-";
    }
    return "?";
}

struct Candidate {
    u32 q{};
    u32 t{};
    u32 p{};
    u32 x{};
    u32 y{};
    u32 h{};
    u32 zp_size{};
    Type type{};
    char sign() const { return type == Type::Plus ? '+' : '-'; }
};

static bool candidate_less(const Candidate& a, const Candidate& b) {
    return std::tie(a.q, a.t, a.p, a.type, a.x, a.y)
         < std::tie(b.q, b.t, b.p, b.type, b.x, b.y);
}
static bool logical_same(const Candidate& a, const Candidate& b) {
    return a.q == b.q && a.t == b.t && a.p == b.p && a.sign() == b.sign();
}
struct LeafKey {
    u32 q, t, p;
    bool operator==(const LeafKey& o) const { return q == o.q && t == o.t && p == o.p; }
};
struct LeafHash {
    std::size_t operator()(const LeafKey& k) const {
        u64 x = static_cast<u64>(k.q) * 0x9E3779B185EBCA87ULL;
        x ^= static_cast<u64>(k.t) + 0x9E3779B97F4A7C15ULL + (x << 6) + (x >> 2);
        x ^= static_cast<u64>(k.p) + 0xC2B2AE3D27D4EB4FULL + (x << 6) + (x >> 2);
        return static_cast<std::size_t>(x);
    }
};

static void verify_candidate(const Candidate& c, const std::vector<u32>& zeros) {
    const std::int64_t p = c.p, q = c.q, t = c.t;
    const std::int64_t z = 7 * p - 6 * q - t - 1;
    if (z < 0 || z >= p) throw std::runtime_error("candidate primary digit out of range");
    auto has = [&](u32 a) { return std::binary_search(zeros.begin(), zeros.end(), a); };
    if (!has(static_cast<u32>(z))) throw std::runtime_error("candidate primary is not a zero");
    if (c.type == Type::Plus) {
        const std::int64_t a = 12 * q + t - 13 * p;
        if (a < 0 || a >= p || !has(static_cast<u32>(a))) throw std::runtime_error("bad plus secondary");
        if (static_cast<u32>(z) != c.x || static_cast<u32>(a) != c.p - 1 - c.y)
            throw std::runtime_error("plus normalization mismatch");
    } else {
        const std::int64_t a = q - p - t - 1;
        if (a < 0 || a >= p || !has(static_cast<u32>(a))) throw std::runtime_error("bad minus secondary");
        if (c.type == Type::MPlus) {
            if (static_cast<u32>(a) != c.x || static_cast<u32>(z) != c.y)
                throw std::runtime_error("M+ normalization mismatch");
        } else if (static_cast<u32>(z) != c.x || static_cast<u32>(a) != c.y) {
            throw std::runtime_error("M- normalization mismatch");
        }
    }
}

static void generate_candidates(
    u32 p, const std::vector<u32>& zeros,
    const std::vector<std::uint8_t>& is_prime,
    u32 q_min_exclusive, u32 q_max, u64 b6, u64 b13,
    std::vector<Candidate>& out, u64& admissible_zero_pairs
) {
    if (zeros.size() < 3) return;
    const bool plus_coeff_safe = p > 13 && b6 % p != 0 && b13 % p != 0;
    const bool minus_coeff_safe = p != 7 && b6 % p != 0 && 5 % p != 0;
    for (std::size_t i = 0; i < zeros.size(); ++i) {
        for (std::size_t j = i + 1; j < zeros.size(); ++j) {
            const u32 x = zeros[i], y = zeros[j], h = y - x;
            if ((h & 1U) == 0 || h < 3) continue;
            ++admissible_zero_pairs;
            const auto emit = [&](u32 q, u32 t, Type type) {
                if (q <= q_min_exclusive || q > q_max || t >= q || !is_prime[q]) return;
                Candidate c{q, t, p, x, y, h, static_cast<u32>(zeros.size()), type};
                verify_candidate(c, zeros);
                out.push_back(c);
            };

            if (plus_coeff_safe && p > h && (p - h) % 12 == 0) {
                const u32 d = (p - h) / 6;
                if ((d & 1U) == 0 && d > 0 && x < std::min(h, p - h)) emit(p + d, h - 1 - x, Type::Plus);
            }
            if (minus_coeff_safe && p > h && (p - h) % 14 == 0) {
                const u32 d = (p - h) / 7;
                if ((d & 1U) == 0 && d > 0 && x < d) emit(p + d, d - 1 - x, Type::MPlus);
            }
            if (minus_coeff_safe && (static_cast<u64>(p) + h) % 14 == 0) {
                const u32 d = (p + h) / 7;
                if ((d & 1U) == 0 && d > h && x < d - h) emit(p + d, d - h - 1 - x, Type::MMinus);
            }
        }
    }
}

static void fastmod_self_test() {
    for (u32 p : {7U, 11U, 101U, 1009U, 999983U, 9999991U}) {
        FastMod mod{p};
        u64 state = 0x123456789ABCDEF0ULL ^ p;
        for (int i = 0; i < 100000; ++i) {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            const u32 a = static_cast<u32>(state % p);
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            const u32 b = static_cast<u32>(state % p);
            if (mod.mul(a, b) != static_cast<u32>((static_cast<u64>(a) * b) % p))
                throw std::runtime_error("FastMod self-test failed");
        }
    }
}

struct Evaluated {
    Candidate c;
    u32 zq_size{};
    std::int64_t nearest_zero{-1};
    std::int64_t linear_distance{-1};
    std::int64_t cyclic_distance{-1};
    bool selected{};
};

struct Options {
    u32 q_min_exclusive = 1000000;
    u32 q_max = 2000000;
    unsigned shard_index = 0;
    unsigned shard_count = 1;
    unsigned threads = std::max(1U, std::thread::hardware_concurrency());
    fs::path out_dir = ".";
};

static Options parse_options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&]() -> std::string {
            if (++i >= argc) throw std::runtime_error("missing value after " + a);
            return argv[i];
        };
        if (a == "--q-min-exclusive") o.q_min_exclusive = static_cast<u32>(std::stoul(need()));
        else if (a == "--q-max") o.q_max = static_cast<u32>(std::stoul(need()));
        else if (a == "--shard-index") o.shard_index = static_cast<unsigned>(std::stoul(need()));
        else if (a == "--shard-count") o.shard_count = static_cast<unsigned>(std::stoul(need()));
        else if (a == "--threads") o.threads = static_cast<unsigned>(std::stoul(need()));
        else if (a == "--out-dir") o.out_dir = need();
        else throw std::runtime_error("unknown argument: " + a);
    }
    if (o.q_max <= o.q_min_exclusive || o.shard_count == 0 || o.shard_index >= o.shard_count || o.threads == 0)
        throw std::runtime_error("invalid range/shard/thread options");
    return o;
}

static unsigned shard_for_prime(u32 p, u32 p_min, u32 p_max_exclusive, unsigned count) {
    const u128 p2 = static_cast<u128>(p) * p;
    const u128 lo2 = static_cast<u128>(p_min) * p_min;
    const u128 hi2 = static_cast<u128>(p_max_exclusive) * p_max_exclusive;
    if (hi2 <= lo2) return 0;
    u128 bucket = (p2 - lo2) * count / (hi2 - lo2);
    if (bucket >= count) bucket = count - 1;
    return static_cast<unsigned>(bucket);
}

int main(int argc, char** argv) {
    try {
        const Options opt = parse_options(argc, argv);
        fastmod_self_test();
        const u64 b6 = apery_small(6), b13 = apery_small(13);
        fs::create_directories(opt.out_dir);
        const auto start = std::chrono::steady_clock::now();

        std::vector<std::uint8_t> is_prime;
        const auto primes = prime_sieve(opt.q_max, is_prime);
        const u32 p_min = std::max<u32>(7, static_cast<u32>((static_cast<u64>(6) * opt.q_min_exclusive) / 7 > 2
            ? (static_cast<u64>(6) * opt.q_min_exclusive) / 7 - 2 : 7));

        std::vector<u32> lower;
        for (u32 p : primes) {
            if (p < p_min || p >= opt.q_max) continue;
            if (shard_for_prime(p, p_min, opt.q_max, opt.shard_count) == opt.shard_index) lower.push_back(p);
        }

        std::atomic<std::size_t> next{0};
        std::atomic<u64> done{0};
        std::mutex progress_mu;
        std::vector<std::vector<Candidate>> local_candidates(opt.threads);
        std::vector<u64> local_steps(opt.threads, 0), local_pairs(opt.threads, 0), local_primes(opt.threads, 0);

        auto lower_worker = [&](unsigned tid) {
            while (true) {
                const std::size_t idx = next.fetch_add(1, std::memory_order_relaxed);
                if (idx >= lower.size()) break;
                const u32 p = lower[idx];
                auto zeros = apery_zeros(p, local_steps[tid]);
                ++local_primes[tid];
                generate_candidates(p, zeros, is_prime, opt.q_min_exclusive, opt.q_max,
                                    b6, b13, local_candidates[tid], local_pairs[tid]);
                const u64 n = done.fetch_add(1, std::memory_order_relaxed) + 1;
                if (n % 2000 == 0) {
                    std::lock_guard<std::mutex> lock(progress_mu);
                    std::cerr << "lower-progress shard=" << opt.shard_index << " " << n << "/" << lower.size() << " p=" << p << "\n";
                }
            }
        };
        std::vector<std::thread> pool;
        for (unsigned t = 0; t < opt.threads; ++t) pool.emplace_back(lower_worker, t);
        for (auto& th : pool) th.join();

        std::vector<Candidate> generated;
        u64 lower_steps = 0, admissible_pairs = 0, lower_primes_scanned = 0;
        for (unsigned t = 0; t < opt.threads; ++t) {
            lower_steps += local_steps[t];
            admissible_pairs += local_pairs[t];
            lower_primes_scanned += local_primes[t];
            generated.insert(generated.end(), local_candidates[t].begin(), local_candidates[t].end());
        }
        std::sort(generated.begin(), generated.end(), candidate_less);
        std::vector<Candidate> logical;
        for (const Candidate& c : generated) if (logical.empty() || !logical_same(logical.back(), c)) logical.push_back(c);

        std::unordered_set<LeafKey, LeafHash> minus_keys;
        for (const Candidate& c : logical) if (c.sign() == '-') minus_keys.insert({c.q, c.t, c.p});
        std::vector<Candidate> assigned;
        for (const Candidate& c : logical) {
            if (c.sign() == '+' && minus_keys.count({c.q, c.t, c.p})) continue;
            assigned.push_back(c);
        }
        std::sort(assigned.begin(), assigned.end(), candidate_less);

        struct Group { std::size_t begin, end; u32 q; };
        std::vector<Group> groups;
        for (std::size_t i = 0; i < assigned.size();) {
            std::size_t j = i + 1;
            while (j < assigned.size() && assigned[j].q == assigned[i].q) ++j;
            groups.push_back({i, j, assigned[i].q});
            i = j;
        }

        std::vector<Evaluated> evaluated(assigned.size());
        std::atomic<std::size_t> qnext{0};
        std::vector<u64> upper_steps_by_thread(opt.threads, 0);
        auto upper_worker = [&](unsigned tid) {
            while (true) {
                const std::size_t gi = qnext.fetch_add(1, std::memory_order_relaxed);
                if (gi >= groups.size()) break;
                const Group g = groups[gi];
                const auto zq = apery_zeros(g.q, upper_steps_by_thread[tid]);
                for (std::size_t i = g.begin; i < g.end; ++i) {
                    const Candidate& c = assigned[i];
                    Evaluated e;
                    e.c = c;
                    e.zq_size = static_cast<u32>(zq.size());
                    e.selected = std::binary_search(zq.begin(), zq.end(), c.t);
                    if (!zq.empty()) {
                        auto it = std::lower_bound(zq.begin(), zq.end(), c.t);
                        u32 best = zq.front();
                        u64 bestd = std::numeric_limits<u64>::max();
                        auto consider = [&](u32 z) {
                            const u64 d = z > c.t ? static_cast<u64>(z - c.t) : static_cast<u64>(c.t - z);
                            if (d < bestd || (d == bestd && z < best)) { bestd = d; best = z; }
                        };
                        if (it != zq.end()) consider(*it);
                        if (it != zq.begin()) consider(*std::prev(it));
                        consider(zq.front());
                        consider(zq.back());
                        e.nearest_zero = best;
                        e.linear_distance = static_cast<std::int64_t>(bestd);
                        e.cyclic_distance = static_cast<std::int64_t>(std::min<u64>(bestd, static_cast<u64>(g.q) - bestd));
                    }
                    if (e.zq_size == 1 && c.sign() == '-' && e.selected)
                        throw std::runtime_error("upper singleton minus restriction violated");
                    evaluated[i] = e;
                }
            }
        };
        pool.clear();
        for (unsigned t = 0; t < opt.threads; ++t) pool.emplace_back(upper_worker, t);
        for (auto& th : pool) th.join();
        u64 upper_steps = 0;
        for (u64 x : upper_steps_by_thread) upper_steps += x;

        const fs::path csv_path = opt.out_dir / ("candidates-" + std::to_string(opt.shard_index) + ".csv");
        std::ofstream csv(csv_path);
        if (!csv) throw std::runtime_error("cannot create candidate CSV");
        csv << "q,t,p,sign,type,x,y,h,zp_size,zq_size,selected,nearest_zero,linear_distance,cyclic_distance\n";
        u64 selected_plus = 0, selected_minus = 0, assigned_plus = 0, assigned_minus = 0;
        for (const Evaluated& e : evaluated) {
            if (e.c.sign() == '+') ++assigned_plus; else ++assigned_minus;
            if (e.selected) { if (e.c.sign() == '+') ++selected_plus; else ++selected_minus; }
            csv << e.c.q << ',' << e.c.t << ',' << e.c.p << ',' << e.c.sign() << ',' << type_name(e.c.type)
                << ',' << e.c.x << ',' << e.c.y << ',' << e.c.h << ',' << e.c.zp_size << ',' << e.zq_size
                << ',' << (e.selected ? 1 : 0) << ',' << e.nearest_zero << ',' << e.linear_distance << ',' << e.cyclic_distance << '\n';
        }
        csv.close();

        const auto finish = std::chrono::steady_clock::now();
        const double wall = std::chrono::duration<double>(finish - start).count();
        const fs::path summary_path = opt.out_dir / ("summary-" + std::to_string(opt.shard_index) + ".json");
        std::ofstream js(summary_path);
        js << "{\n"
           << "  \"q_min_exclusive\": " << opt.q_min_exclusive << ",\n"
           << "  \"q_max\": " << opt.q_max << ",\n"
           << "  \"shard_index\": " << opt.shard_index << ",\n"
           << "  \"shard_count\": " << opt.shard_count << ",\n"
           << "  \"threads\": " << opt.threads << ",\n"
           << "  \"p_min_global\": " << p_min << ",\n"
           << "  \"lower_primes_scanned\": " << lower_primes_scanned << ",\n"
           << "  \"lower_recurrence_steps\": " << lower_steps << ",\n"
           << "  \"admissible_opposite_parity_zero_pairs\": " << admissible_pairs << ",\n"
           << "  \"generated_records\": " << generated.size() << ",\n"
           << "  \"logical_raw_records\": " << logical.size() << ",\n"
           << "  \"assigned_plus\": " << assigned_plus << ",\n"
           << "  \"assigned_minus\": " << assigned_minus << ",\n"
           << "  \"candidate_upper_primes\": " << groups.size() << ",\n"
           << "  \"upper_recurrence_steps\": " << upper_steps << ",\n"
           << "  \"selected_plus\": " << selected_plus << ",\n"
           << "  \"selected_minus\": " << selected_minus << ",\n"
           << "  \"wall_seconds\": " << std::fixed << std::setprecision(6) << wall << "\n"
           << "}\n";
        js.close();

        std::cout << summary_path << '\n' << csv_path << '\n';
        return (selected_plus || selected_minus) ? 10 : 0;
    } catch (const std::exception& e) {
        std::cerr << "fatal: " << e.what() << '\n';
        return 2;
    }
}
