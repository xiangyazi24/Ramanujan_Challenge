// Q7174 exact projective-state collision scan.
//
// For each prime p, compute the two recurrence solutions
//   B_0=1, B_1=5 and C_0=0, C_1=1
// modulo p, form q_n=[B_n:C_n] represented by B_n/C_n or p for infinity,
// pack each projective fiber into an ordered occurrence list, and inspect only
// sliding windows of four consecutive occurrences.  Thus the scan is O(p)
// per prime, not O(p^2).
//
// A window z0<z1<z2<z3 is counted exactly when:
//   * all three gaps are at least 2;
//   * it is primitive (automatic from consecutive fiber occurrences);
//   * the gap triple is not all equal;
//   * no adjacent pair sums to p-1 (off-center);
//   * (z3-z0)^2 <= p;
//   * p does not divide U_s = product_{j<=s} j! B_j V_{j+1},
//     V_0=0,V_1=1,V_{j+1}=34V_j-V_{j-1}.
//
// Every recurrence solution is palindromic modulo p, so q_{p-1-n}=q_n.
// The reflection of (z0,z1,z2,z3) is
//   (p-1-z3,p-1-z2,p-1-z1,p-1-z0).
// Raw counts and reflection-orbit counts are both reported.  Since centered
// adjacent pairs are excluded, no valid chain is fixed by reflection, and the
// raw count is exactly twice the orbit count.

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

int add_mod(i64 x, i64 y, int p) {
    i64 z = (x + y) % p;
    if (z < 0) z += p;
    return static_cast<int>(z);
}

int mul_mod(i64 x, i64 y, int p) {
    return static_cast<int>((x * y) % p);
}

std::vector<int> primes_up_to(int limit) {
    std::vector<unsigned char> sieve(static_cast<std::size_t>(limit) + 1, 1);
    if (limit >= 0) sieve[0] = 0;
    if (limit >= 1) sieve[1] = 0;
    for (int p = 2; static_cast<i64>(p) * p <= limit; ++p) {
        if (!sieve[p]) continue;
        for (int m = p * p; m <= limit; m += p) sieve[m] = 0;
    }
    std::vector<int> primes;
    for (int p = 7; p <= limit; ++p) {
        if (sieve[p]) primes.push_back(p);
    }
    return primes;
}

struct Chain {
    std::array<int, 4> z{};
};

bool operator==(const Chain& a, const Chain& b) { return a.z == b.z; }

Chain reflected(const Chain& chain, int p) {
    return Chain{{p - 1 - chain.z[3], p - 1 - chain.z[2],
                  p - 1 - chain.z[1], p - 1 - chain.z[0]}};
}

std::array<int, 3> gaps(const Chain& chain) {
    return {chain.z[1] - chain.z[0], chain.z[2] - chain.z[1],
            chain.z[3] - chain.z[2]};
}

std::string vector_text(const std::vector<int>& values) {
    std::ostringstream out;
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << ',';
        out << values[i];
    }
    out << ']';
    return out.str();
}

std::string chain_text(const Chain& chain) {
    const auto g = gaps(chain);
    std::ostringstream out;
    out << "start=" << chain.z[0] << " gaps=" << g[0] << ',' << g[1]
        << ',' << g[2] << " zeros=" << chain.z[0] << ',' << chain.z[1]
        << ',' << chain.z[2] << ',' << chain.z[3];
    return out.str();
}

struct Witness {
    int prime = 0;
    int state = 0;
    std::vector<int> occurrences;
    std::vector<int> actual_zeros;
    std::vector<Chain> canonical_chains;
};

bool earlier(const Witness& a, const Witness& b) {
    if (b.prime == 0) return a.prime != 0;
    if (a.prime == 0) return false;
    if (a.prime != b.prime) return a.prime < b.prime;
    return a.state < b.state;
}

u64 fnv_mix(u64 hash, u64 value) {
    constexpr u64 prime = 1099511628211ULL;
    for (int byte = 0; byte < 8; ++byte) {
        hash ^= (value >> (8 * byte)) & 0xffULL;
        hash *= prime;
    }
    return hash;
}

struct PrimeResult {
    int prime = 0;
    std::uint64_t raw_total = 0;
    std::uint64_t orbit_total = 0;
    int max_raw = 0;
    int max_raw_state = -1;
    int max_orbits = 0;
    int max_orbit_state = -1;
    int actual_raw = 0;
    int actual_orbits = 0;
    u64 checksum = 1469598103934665603ULL;
    std::optional<Witness> first_raw_duplicate;
    std::optional<Witness> first_nonreflection_duplicate;
};

PrimeResult scan_prime(int p) {
    PrimeResult result;
    result.prime = p;

    std::vector<int> inverse(p, 0);
    inverse[1] = 1;
    for (int v = 2; v < p; ++v) {
        inverse[v] = static_cast<int>(
            p - (static_cast<i64>(p / v) * inverse[p % v]) % p);
        if (inverse[v] == p) inverse[v] = 0;
    }

    std::vector<int> apery(p, 0), companion(p, 0);
    apery[0] = 1;
    apery[1] = 5 % p;
    companion[0] = 0;
    companion[1] = 1;

    for (int n = 1; n < p - 1; ++n) {
        const int n2 = mul_mod(n, n, p);
        const int n3 = mul_mod(n2, n, p);
        int coefficient = 5 % p;
        coefficient = add_mod(coefficient, static_cast<i64>(27) * n, p);
        coefficient = add_mod(coefficient, static_cast<i64>(51) * n2, p);
        coefficient = add_mod(coefficient, static_cast<i64>(34) * n3, p);
        const int inv_den = mul_mod(mul_mod(inverse[n + 1], inverse[n + 1], p),
                                    inverse[n + 1], p);

        auto advance = [&](const std::vector<int>& u) {
            int numerator = mul_mod(coefficient, u[n], p);
            numerator = add_mod(numerator, -static_cast<i64>(n3) * u[n - 1], p);
            return mul_mod(numerator, inv_den, p);
        };
        apery[n + 1] = advance(apery);
        companion[n + 1] = advance(companion);
    }

    // This is both an implementation invariant and the key reason raw chains
    // occur in reflection pairs inside the same projective state.
    for (int n = 0; n < p; ++n) {
        if (apery[n] != apery[p - 1 - n] ||
            companion[n] != companion[p - 1 - n]) {
            throw std::runtime_error("palindrome failure at p=" +
                                     std::to_string(p));
        }
        if (apery[n] == 0 && companion[n] == 0) {
            throw std::runtime_error("projective basis vanished at p=" +
                                     std::to_string(p));
        }
    }

    std::vector<int> state_of(p, 0);
    std::vector<int> counts(p + 1, 0);
    std::vector<int> actual_zeros;
    for (int n = 0; n < p; ++n) {
        const int q = companion[n] == 0
                          ? p
                          : mul_mod(apery[n], inverse[companion[n]], p);
        state_of[n] = q;
        ++counts[q];
        if (apery[n] == 0) actual_zeros.push_back(n);
        if ((q == 0) != (apery[n] == 0)) {
            throw std::runtime_error("q=0 mismatch at p=" + std::to_string(p));
        }
        if (q != state_of[p - 1 - n] && p - 1 - n < n) {
            throw std::runtime_error("state reflection failure at p=" +
                                     std::to_string(p));
        }
    }

    std::vector<unsigned char> carrier_bad(p, 0);
    int v_previous = 0;
    int v_current = 1;
    for (int s = 1; s < p; ++s) {
        const int v_next = add_mod(static_cast<i64>(34) * v_current,
                                   -v_previous, p);
        v_previous = v_current;
        v_current = v_next;
        carrier_bad[s] = static_cast<unsigned char>(
            carrier_bad[s - 1] || apery[s] == 0 || v_current == 0);
    }

    std::vector<int> offset(p + 2, 0);
    for (int q = 0; q <= p; ++q) offset[q + 1] = offset[q] + counts[q];
    if (offset[p + 1] != p) throw std::runtime_error("packed size mismatch");
    std::vector<int> cursor = offset;
    std::vector<int> packed(p, 0);
    for (int n = 0; n < p; ++n) packed[cursor[state_of[n]]++] = n;

    for (int q = 0; q <= p; ++q) {
        const int begin = offset[q];
        const int end = offset[q + 1];
        const int count = end - begin;
        if (count < 4) continue;

        for (int j = 0; j < count; ++j) {
            if (packed[begin + j] + packed[end - 1 - j] != p - 1) {
                throw std::runtime_error("fiber palindrome failure at p=" +
                                         std::to_string(p));
            }
        }

        std::vector<Chain> raw_chains;
        std::vector<Chain> canonical_chains;
        raw_chains.reserve(static_cast<std::size_t>(count - 3));
        canonical_chains.reserve(static_cast<std::size_t>((count - 3 + 1) / 2));

        for (int j = 0; j + 3 < count; ++j) {
            Chain chain{{packed[begin + j], packed[begin + j + 1],
                         packed[begin + j + 2], packed[begin + j + 3]}};
            const auto g = gaps(chain);
            const int span = chain.z[3] - chain.z[0];
            const bool valid_gaps = g[0] >= 2 && g[1] >= 2 && g[2] >= 2;
            const bool non_all_equal = !(g[0] == g[1] && g[1] == g[2]);
            const bool off_center =
                chain.z[0] + chain.z[1] != p - 1 &&
                chain.z[1] + chain.z[2] != p - 1 &&
                chain.z[2] + chain.z[3] != p - 1;
            const bool short_span = static_cast<i64>(span) * span <= p;
            const bool carrier_clean = !carrier_bad[span];
            if (!(valid_gaps && non_all_equal && off_center && short_span &&
                  carrier_clean)) {
                continue;
            }

            raw_chains.push_back(chain);
            const Chain mirror = reflected(chain, p);
            if (chain.z < mirror.z) canonical_chains.push_back(chain);
        }

        if (raw_chains.empty()) continue;
        if (raw_chains.size() != 2 * canonical_chains.size()) {
            throw std::runtime_error("raw/orbit parity failure at p=" +
                                     std::to_string(p) + " q=" +
                                     std::to_string(q));
        }
        for (std::size_t j = 0; j < raw_chains.size(); ++j) {
            if (!(raw_chains[j] ==
                  reflected(raw_chains[raw_chains.size() - 1 - j], p))) {
                throw std::runtime_error("reflection pairing failure at p=" +
                                         std::to_string(p) + " q=" +
                                         std::to_string(q));
            }
        }

        const int raw_count = static_cast<int>(raw_chains.size());
        const int orbit_count = static_cast<int>(canonical_chains.size());
        result.raw_total += raw_count;
        result.orbit_total += orbit_count;
        if (raw_count > result.max_raw ||
            (raw_count == result.max_raw && q < result.max_raw_state)) {
            result.max_raw = raw_count;
            result.max_raw_state = q;
        }
        if (orbit_count > result.max_orbits ||
            (orbit_count == result.max_orbits && q < result.max_orbit_state)) {
            result.max_orbits = orbit_count;
            result.max_orbit_state = q;
        }
        if (q == 0) {
            result.actual_raw = raw_count;
            result.actual_orbits = orbit_count;
        }

        for (const Chain& chain : canonical_chains) {
            result.checksum = fnv_mix(result.checksum, static_cast<u64>(p));
            result.checksum = fnv_mix(result.checksum, static_cast<u64>(q));
            for (int z : chain.z) {
                result.checksum = fnv_mix(result.checksum, static_cast<u64>(z));
            }
        }

        std::vector<int> occurrences(packed.begin() + begin, packed.begin() + end);
        if (raw_count >= 2 && !result.first_raw_duplicate) {
            result.first_raw_duplicate =
                Witness{p, q, occurrences, actual_zeros, canonical_chains};
        }
        if (orbit_count >= 2 && !result.first_nonreflection_duplicate) {
            result.first_nonreflection_duplicate =
                Witness{p, q, occurrences, actual_zeros, canonical_chains};
        }
    }

    return result;
}

struct Options {
    int limit = 100000;
    int shard = 0;
    int shards = 1;
    int threads = 1;
};

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&]() -> int {
            if (++i >= argc) throw std::runtime_error("missing option value");
            return std::stoi(argv[i]);
        };
        if (arg == "--limit") {
            options.limit = require_value();
        } else if (arg == "--shard") {
            options.shard = require_value();
        } else if (arg == "--shards") {
            options.shards = require_value();
        } else if (arg == "--threads") {
            options.threads = require_value();
        } else {
            throw std::runtime_error("unknown option: " + arg);
        }
    }
    if (options.limit < 7 || options.shards < 1 || options.shard < 0 ||
        options.shard >= options.shards || options.threads < 1) {
        throw std::runtime_error("invalid options");
    }
    return options;
}

struct Aggregate {
    std::uint64_t primes = 0;
    std::uint64_t raw_total = 0;
    std::uint64_t orbit_total = 0;
    int max_raw = 0;
    int max_raw_prime = 0;
    int max_raw_state = -1;
    int max_orbits = 0;
    int max_orbit_prime = 0;
    int max_orbit_state = -1;
    u64 checksum = 1469598103934665603ULL;
    std::optional<Witness> first_raw;
    std::optional<Witness> first_nonreflection;
};

void merge_prime(const PrimeResult& prime_result, Aggregate& aggregate) {
    ++aggregate.primes;
    aggregate.raw_total += prime_result.raw_total;
    aggregate.orbit_total += prime_result.orbit_total;
    aggregate.checksum = fnv_mix(aggregate.checksum,
                                 static_cast<u64>(prime_result.prime));
    aggregate.checksum = fnv_mix(aggregate.checksum, prime_result.checksum);

    if (prime_result.max_raw > aggregate.max_raw ||
        (prime_result.max_raw == aggregate.max_raw &&
         std::pair{prime_result.prime, prime_result.max_raw_state} <
             std::pair{aggregate.max_raw_prime, aggregate.max_raw_state})) {
        aggregate.max_raw = prime_result.max_raw;
        aggregate.max_raw_prime = prime_result.prime;
        aggregate.max_raw_state = prime_result.max_raw_state;
    }
    if (prime_result.max_orbits > aggregate.max_orbits ||
        (prime_result.max_orbits == aggregate.max_orbits &&
         std::pair{prime_result.prime, prime_result.max_orbit_state} <
             std::pair{aggregate.max_orbit_prime, aggregate.max_orbit_state})) {
        aggregate.max_orbits = prime_result.max_orbits;
        aggregate.max_orbit_prime = prime_result.prime;
        aggregate.max_orbit_state = prime_result.max_orbit_state;
    }

    if (prime_result.first_raw_duplicate &&
        (!aggregate.first_raw ||
         earlier(*prime_result.first_raw_duplicate, *aggregate.first_raw))) {
        aggregate.first_raw = prime_result.first_raw_duplicate;
    }
    if (prime_result.first_nonreflection_duplicate &&
        (!aggregate.first_nonreflection ||
         earlier(*prime_result.first_nonreflection_duplicate,
                 *aggregate.first_nonreflection))) {
        aggregate.first_nonreflection =
            prime_result.first_nonreflection_duplicate;
    }
}

void print_witness(const char* label, const std::optional<Witness>& witness) {
    if (!witness) {
        std::cout << label << " none\n";
        return;
    }
    const Witness& w = *witness;
    std::cout << label << " p=" << w.prime << " q=";
    if (w.state == w.prime) {
        std::cout << "inf";
    } else {
        std::cout << w.state;
    }
    std::cout << " actual=" << (w.state == 0 ? 1 : 0)
              << " occurrences=" << vector_text(w.occurrences)
              << " actual_zero_set=" << vector_text(w.actual_zeros)
              << " reflection_orbits=" << w.canonical_chains.size() << '\n';
    for (std::size_t i = 0; i < w.canonical_chains.size(); ++i) {
        const Chain& chain = w.canonical_chains[i];
        std::cout << label << " orbit=" << i + 1 << " canonical "
                  << chain_text(chain) << '\n';
        std::cout << label << " orbit=" << i + 1 << " reflected "
                  << chain_text(reflected(chain, w.prime)) << '\n';
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        const std::vector<int> all_primes = primes_up_to(options.limit);
        std::vector<int> assigned;
        for (std::size_t i = 0; i < all_primes.size(); ++i) {
            if (static_cast<int>(i % options.shards) == options.shard) {
                assigned.push_back(all_primes[i]);
            }
        }
        // Descending order balances dynamic worker allocation and makes large
        // cases run while all threads are active.  Earliest witnesses are
        // selected by exact (p,q) comparison, so order cannot affect output.
        std::reverse(assigned.begin(), assigned.end());

        std::atomic<std::size_t> next{0};
        std::mutex mutex;
        Aggregate aggregate;
        std::exception_ptr failure;

        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t i = next.fetch_add(1);
                    if (i >= assigned.size()) break;
                    const PrimeResult result = scan_prime(assigned[i]);
                    std::lock_guard<std::mutex> lock(mutex);
                    merge_prime(result, aggregate);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(mutex);
                if (!failure) failure = std::current_exception();
            }
        };

        std::vector<std::thread> workers;
        for (int i = 0; i < options.threads; ++i) workers.emplace_back(worker);
        for (auto& thread : workers) thread.join();
        if (failure) std::rethrow_exception(failure);

        std::cout << "Q7174_SUMMARY limit=" << options.limit
                  << " shard=" << options.shard << '/' << options.shards
                  << " threads=" << options.threads
                  << " primes=" << aggregate.primes
                  << " raw_total=" << aggregate.raw_total
                  << " orbit_total=" << aggregate.orbit_total
                  << " max_raw=" << aggregate.max_raw
                  << " max_raw_at=" << aggregate.max_raw_prime << ':'
                  << aggregate.max_raw_state
                  << " max_reflection_quotient=" << aggregate.max_orbits
                  << " max_quotient_at=" << aggregate.max_orbit_prime << ':'
                  << aggregate.max_orbit_state
                  << " checksum=0x" << std::hex << aggregate.checksum << std::dec
                  << '\n';
        print_witness("Q7174_FIRST_RAW", aggregate.first_raw);
        print_witness("Q7174_FIRST_NONREFLECTION", aggregate.first_nonreflection);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Q7174_FAIL " << error.what() << '\n';
        return 2;
    }
}
