#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
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

    u32 mul(u32 a, u32 b) const {
        return reduce(static_cast<u64>(a) * b);
    }
};

static u32 add_mod(u32 a, u32 b, u32 p) {
    u32 c = a + b;
    if (c >= p) c -= p;
    return c;
}

static u32 sub_mod(u32 a, u32 b, u32 p) {
    return a >= b ? a - b : static_cast<u32>(static_cast<u64>(a) + p - b);
}

template <std::size_t N>
static void advance_differences(std::array<u32, N>& d, u32 p) {
    for (std::size_t k = 0; k + 1 < N; ++k) {
        d[k] = add_mod(d[k], d[k + 1], p);
    }
}

static std::vector<u32> prime_sieve(u32 n, std::vector<std::uint8_t>& is_prime) {
    is_prime.assign(static_cast<std::size_t>(n) + 1, 1);
    if (n >= 0) is_prime[0] = 0;
    if (n >= 1) is_prime[1] = 0;
    for (u64 i = 2; i * i <= n; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) continue;
        for (u64 j = i * i; j <= n; j += i) {
            is_prime[static_cast<std::size_t>(j)] = 0;
        }
    }
    std::vector<u32> primes;
    for (u32 i = 2; i <= n; ++i) {
        if (is_prime[i]) primes.push_back(i);
    }
    return primes;
}

static std::vector<u32> apery_zeros(u32 p, u64& recurrence_steps) {
    if (p < 7) return {};
    FastMod mod{p};

    std::array<u32, 4> dp{
        static_cast<u32>(117 % p),
        static_cast<u32>(418 % p),
        static_cast<u32>(510 % p),
        static_cast<u32>(204 % p)
    };
    std::array<u32, 7> d6{
        static_cast<u32>(1 % p),
        static_cast<u32>(63 % p),
        static_cast<u32>(602 % p),
        static_cast<u32>(2100 % p),
        static_cast<u32>(3360 % p),
        static_cast<u32>(2520 % p),
        static_cast<u32>(720 % p)
    };

    const u32 midpoint = (p - 1) / 2;
    u32 previous = 1 % p;
    u32 current = 5 % p;
    std::vector<u32> half;
    if (previous == 0) half.push_back(0);
    if (midpoint >= 1 && current == 0) half.push_back(1);

    for (u32 n = 1; n < midpoint; ++n) {
        const u32 term1 = mod.mul(dp[0], current);
        const u32 term2 = mod.mul(d6[0], previous);
        const u32 next = sub_mod(term1, term2, p);
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

struct Candidate {
    u32 q{};
    u32 t{};
    u32 p{};
    u32 x{};
    u32 y{};
    u32 h{};
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
        if (a < 0 || a >= p) throw std::runtime_error("plus secondary out of range");
        if (!has(static_cast<u32>(a))) throw std::runtime_error("plus secondary is not a zero");
        if (static_cast<u32>(z) != c.x || static_cast<u32>(a) != c.p - 1 - c.y)
            throw std::runtime_error("plus inverse normalization mismatch");
    } else {
        const std::int64_t a = q - p - t - 1;
        if (a < 0 || a >= p) throw std::runtime_error("minus secondary out of range");
        if (!has(static_cast<u32>(a))) throw std::runtime_error("minus secondary is not a zero");
        if (c.type == Type::MPlus) {
            if (static_cast<u32>(a) != c.x || static_cast<u32>(z) != c.y)
                throw std::runtime_error("M+ inverse normalization mismatch");
        } else {
            if (static_cast<u32>(z) != c.x || static_cast<u32>(a) != c.y)
                throw std::runtime_error("M- inverse normalization mismatch");
        }
    }
}

static void generate_candidates(
    u32 p,
    const std::vector<u32>& zeros,
    const std::vector<std::uint8_t>& is_prime,
    u32 q_max,
    u64 b6,
    u64 b13,
    std::vector<Candidate>& out,
    u64& zero_pairs
) {
    if (zeros.size() < 3) return;
    const bool plus_coeff_safe = p > 13 && b6 % p != 0 && b13 % p != 0;
    const bool minus_coeff_safe = p != 7 && b6 % p != 0 && 5 % p != 0;

    for (std::size_t i = 0; i < zeros.size(); ++i) {
        for (std::size_t j = i + 1; j < zeros.size(); ++j) {
            const u32 x = zeros[i];
            const u32 y = zeros[j];
            const u32 h = y - x;
            if ((h & 1U) == 0 || h < 3) continue;
            ++zero_pairs;

            if (plus_coeff_safe && p > h && (p - h) % 12 == 0) {
                const u32 d = (p - h) / 6;
                if ((d & 1U) == 0 && d > 0 && x < std::min(h, p - h)) {
                    const u64 q64 = static_cast<u64>(p) + d;
                    const u32 t = h - 1 - x;
                    if (q64 <= q_max && t < q64 && is_prime[static_cast<std::size_t>(q64)]) {
                        Candidate c{static_cast<u32>(q64), t, p, x, y, h, Type::Plus};
                        verify_candidate(c, zeros);
                        out.push_back(c);
                    }
                }
            }

            if (minus_coeff_safe && p > h && (p - h) % 14 == 0) {
                const u32 d = (p - h) / 7;
                if ((d & 1U) == 0 && d > 0 && x < d) {
                    const u64 q64 = static_cast<u64>(p) + d;
                    const u32 t = d - 1 - x;
                    if (q64 <= q_max && t < q64 && is_prime[static_cast<std::size_t>(q64)]) {
                        Candidate c{static_cast<u32>(q64), t, p, x, y, h, Type::MPlus};
                        verify_candidate(c, zeros);
                        out.push_back(c);
                    }
                }
            }

            if (minus_coeff_safe && (static_cast<u64>(p) + h) % 14 == 0) {
                const u32 d = (p + h) / 7;
                if ((d & 1U) == 0 && d > h && x < d - h) {
                    const u64 q64 = static_cast<u64>(p) + d;
                    const u32 t = d - h - 1 - x;
                    if (q64 <= q_max && t < q64 && is_prime[static_cast<std::size_t>(q64)]) {
                        Candidate c{static_cast<u32>(q64), t, p, x, y, h, Type::MMinus};
                        verify_candidate(c, zeros);
                        out.push_back(c);
                    }
                }
            }
        }
    }
}

static void fastmod_self_test() {
    for (u32 p : {7U, 11U, 101U, 1009U, 999983U}) {
        FastMod mod{p};
        u64 state = 0x123456789ABCDEF0ULL ^ p;
        for (int i = 0; i < 100000; ++i) {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            const u32 a = static_cast<u32>(state % p);
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            const u32 b = static_cast<u32>(state % p);
            const u32 got = mod.mul(a, b);
            const u32 want = static_cast<u32>((static_cast<u64>(a) * b) % p);
            if (got != want) throw std::runtime_error("FastMod self-test failed");
        }
    }
}

int main(int argc, char** argv) {
    try {
        u32 q_max = 1000000;
        unsigned threads = std::max(1U, std::thread::hardware_concurrency());
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--q-max" && i + 1 < argc) q_max = static_cast<u32>(std::stoul(argv[++i]));
            else if (a == "--threads" && i + 1 < argc) threads = static_cast<unsigned>(std::stoul(argv[++i]));
            else throw std::runtime_error("unknown or incomplete argument: " + a);
        }

        fastmod_self_test();
        const u64 b6 = apery_small(6);
        const u64 b13 = apery_small(13);
        std::cerr << "b6=" << b6 << " b13=" << b13 << "\n";

        const auto wall0 = std::chrono::steady_clock::now();
        std::vector<std::uint8_t> is_prime;
        const auto primes = prime_sieve(q_max, is_prime);
        std::vector<u32> lower;
        for (u32 p : primes) if (p >= 7 && p < q_max) lower.push_back(p);

        std::atomic<std::size_t> next{0};
        std::vector<std::vector<Candidate>> local_candidates(threads);
        std::vector<u64> local_steps(threads, 0), local_pairs(threads, 0), local_primes(threads, 0);
        std::atomic<u64> done{0};
        std::mutex progress_mu;

        auto worker = [&](unsigned tid) {
            auto& out = local_candidates[tid];
            while (true) {
                const std::size_t idx = next.fetch_add(1, std::memory_order_relaxed);
                if (idx >= lower.size()) break;
                const u32 p = lower[idx];
                auto zeros = apery_zeros(p, local_steps[tid]);
                ++local_primes[tid];
                generate_candidates(p, zeros, is_prime, q_max, b6, b13, out, local_pairs[tid]);
                const u64 n = done.fetch_add(1, std::memory_order_relaxed) + 1;
                if (n % 5000 == 0) {
                    std::lock_guard<std::mutex> lock(progress_mu);
                    std::cerr << "progress=" << n << "/" << lower.size() << " p=" << p << "\n";
                }
            }
        };

        std::vector<std::thread> pool;
        for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker, t);
        for (auto& th : pool) th.join();

        std::vector<Candidate> candidates;
        u64 recurrence_steps = 0, zero_pairs = 0, scanned_primes = 0;
        for (unsigned t = 0; t < threads; ++t) {
            recurrence_steps += local_steps[t];
            zero_pairs += local_pairs[t];
            scanned_primes += local_primes[t];
            candidates.insert(candidates.end(), local_candidates[t].begin(), local_candidates[t].end());
        }
        std::sort(candidates.begin(), candidates.end(), candidate_less);

        std::vector<Candidate> logical;
        for (const Candidate& c : candidates) {
            if (logical.empty() || !logical_same(logical.back(), c)) logical.push_back(c);
        }

        std::unordered_set<LeafKey, LeafHash> minus_keys;
        for (const Candidate& c : logical) if (c.sign() == '-') minus_keys.insert({c.q, c.t, c.p});

        std::vector<Candidate> assigned;
        for (const Candidate& c : logical) {
            if (c.sign() == '+' && minus_keys.count({c.q, c.t, c.p})) continue;
            assigned.push_back(c);
        }

        u64 raw_plus = 0, raw_minus = 0, assigned_plus = 0, assigned_minus = 0;
        for (const Candidate& c : logical) (c.sign() == '+' ? raw_plus : raw_minus)++;
        for (const Candidate& c : assigned) (c.sign() == '+' ? assigned_plus : assigned_minus)++;

        u64 selected_plus = 0, selected_minus = 0;
        for (const Candidate& c : assigned) {
            u64 dummy = 0;
            const auto zq = apery_zeros(c.q, dummy);
            const bool selected = std::binary_search(zq.begin(), zq.end(), c.t);
            if (selected) (c.sign() == '+' ? selected_plus : selected_minus)++;
        }

        const auto wall1 = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(wall1 - wall0).count();

        std::cout << "{\n";
        std::cout << "  \"q_max\": " << q_max << ",\n";
        std::cout << "  \"threads\": " << threads << ",\n";
        std::cout << "  \"prime_count\": " << primes.size() << ",\n";
        std::cout << "  \"lower_primes_scanned\": " << scanned_primes << ",\n";
        std::cout << "  \"lower_recurrence_steps\": " << recurrence_steps << ",\n";
        std::cout << "  \"opposite_parity_zero_pairs\": " << zero_pairs << ",\n";
        std::cout << "  \"generated_records\": " << candidates.size() << ",\n";
        std::cout << "  \"raw_plus\": " << raw_plus << ",\n";
        std::cout << "  \"raw_minus\": " << raw_minus << ",\n";
        std::cout << "  \"assigned_plus\": " << assigned_plus << ",\n";
        std::cout << "  \"assigned_minus\": " << assigned_minus << ",\n";
        std::cout << "  \"selected_plus\": " << selected_plus << ",\n";
        std::cout << "  \"selected_minus\": " << selected_minus << ",\n";
        std::cout << "  \"wall_seconds\": " << std::fixed << std::setprecision(6) << seconds << "\n";
        std::cout << "}\n";

        if (q_max == 1000000 && (assigned_plus != 255 || assigned_minus != 86)) {
            std::cerr << "KNOWN-LEDGER MISMATCH: expected assigned_plus=255 assigned_minus=86\n";
            return 3;
        }
        if (selected_plus != 0 || selected_minus != 0) {
            std::cerr << "ACTUAL LEAF FOUND IN VALIDATION RANGE\n";
            return 4;
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "fatal: " << e.what() << "\n";
        return 2;
    }
}
