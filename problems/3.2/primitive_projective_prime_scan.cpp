// Exact prime-first census of primitive projective Apéry return chains.
//
// For every prime p <= limit, compute two independent solutions of the Apéry
// recurrence modulo p and group the indices 0,...,p-1 by their common point
// [b_n:c_n] in P^1(F_p).  Four consecutive occurrences of one point give an
// exact primitive return chain.  The scan removes centered adjacent pairs,
// the all-equal gap slice, and primes supported on U_s.  It therefore audits
// QPRS directly, without a resultant superset or a factorization cutoff.
// It also restricts to span <= floor(sqrt(p)), computes the raw projective
// chain variance and its fixed-point-free reflection quotient, and reports
// every state supporting at least two quotient chains.
//
// Example:
//   clang++ -std=c++20 -O3 -Wall -Wextra -Wpedantic -Werror \
//     primitive_projective_prime_scan.cpp -o /tmp/primitive_projective_scan
//   /tmp/primitive_projective_scan --limit 200000

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using i128 = __int128_t;

i64 mul_mod(i64 left, i64 right, i64 prime) {
    return static_cast<i64>((static_cast<i128>(left) * right) % prime);
}

i64 add_mod(i64 value, i64 term, i64 prime) {
    value += term;
    value %= prime;
    if (value < 0) value += prime;
    return value;
}

std::vector<int> primes_up_to(int limit) {
    std::vector<bool> sieve(limit + 1, true);
    if (limit >= 0) sieve[0] = false;
    if (limit >= 1) sieve[1] = false;
    for (int prime = 2; static_cast<i64>(prime) * prime <= limit; ++prime) {
        if (!sieve[prime]) continue;
        for (int multiple = prime * prime; multiple <= limit; multiple += prime) {
            sieve[multiple] = false;
        }
    }
    std::vector<int> primes;
    for (int value = 7; value <= limit; ++value) {
        if (sieve[value]) primes.push_back(value);
    }
    return primes;
}

struct Record {
    int prime;
    int start;
    int a;
    int b;
    int c;
    int state;
    bool actual;
};

struct DplsBin {
    std::uint64_t prime_count = 0;
    std::uint64_t raw_mass = 0;
    std::uint64_t raw_energy = 0;
    std::uint64_t quotient_mass = 0;
    std::uint64_t quotient_energy = 0;
    std::uint64_t quotient_overlap_energy = 0;
    std::uint64_t quotient_separated_energy = 0;
    std::uint64_t raw_collision_states = 0;
    std::uint64_t quotient_collision_states = 0;
    std::uint64_t actual_raw = 0;
    long double raw_variance = 0.0L;
    long double quotient_variance = 0.0L;
    int max_raw_count = 0;
    int max_quotient_count = 0;
};

struct Summary {
    std::uint64_t primitive_off_center = 0;
    std::uint64_t ratio_above_one = 0;
    std::uint64_t ratio_above_two = 0;
    int best_prime = 0;
    int best_span = 1;
    Record best_record{};
    std::vector<Record> records_above_one;
    bool saw_1297_forward = false;
    bool saw_1297_reverse = false;
    bool saw_7411_forward = false;
    bool saw_7411_reverse = false;
    bool saw_128047_forward = false;
    bool saw_128047_reverse = false;
    bool saw_dpls_1297_pair = false;
    std::vector<DplsBin> dpls_bins = std::vector<DplsBin>(32);
    std::vector<Record> quotient_collision_records;
};

void consider_record(const Record& record, Summary& summary) {
    const int span = record.a + record.b + record.c;
    ++summary.primitive_off_center;
    const i64 span_squared = static_cast<i64>(span) * span;
    if (record.prime > span_squared) {
        ++summary.ratio_above_one;
        summary.records_above_one.push_back(record);
    }
    if (static_cast<i128>(record.prime) > 2 * static_cast<i128>(span_squared)) {
        ++summary.ratio_above_two;
    }

    if (summary.best_prime == 0 ||
        static_cast<i128>(record.prime) * summary.best_span * summary.best_span >
            static_cast<i128>(summary.best_prime) * span_squared) {
        summary.best_prime = record.prime;
        summary.best_span = span;
        summary.best_record = record;
    }

    summary.saw_1297_forward |=
        record.prime == 1297 && record.start == 360 && record.a == 5 &&
        record.b == 20 && record.c == 10;
    summary.saw_1297_reverse |=
        record.prime == 1297 && record.start == 901 && record.a == 10 &&
        record.b == 20 && record.c == 5;
    summary.saw_7411_forward |=
        record.prime == 7411 && record.start == 4681 && record.a == 3 &&
        record.b == 40 && record.c == 31;
    summary.saw_7411_reverse |=
        record.prime == 7411 && record.start == 2655 && record.a == 31 &&
        record.b == 40 && record.c == 3;
    summary.saw_128047_forward |=
        record.prime == 128047 && record.start == 42375 && record.a == 41 &&
        record.b == 86 && record.c == 37;
    summary.saw_128047_reverse |=
        record.prime == 128047 && record.start == 85507 && record.a == 37 &&
        record.b == 86 && record.c == 41;
}

int integer_square_root(int value) {
    int root = static_cast<int>(std::sqrt(static_cast<double>(value)));
    while (static_cast<i64>(root + 1) * (root + 1) <= value) ++root;
    while (static_cast<i64>(root) * root > value) --root;
    return root;
}

int dyadic_index(int value) {
    int index = 0;
    while (value > 1) {
        value >>= 1;
        ++index;
    }
    return index;
}

bool closed_intervals_meet(int left_a, int right_a,
                           int left_b, int right_b) {
    return std::max(left_a, left_b) <= std::min(right_a, right_b);
}

bool reflection_orbits_overlap(const Record& left, const Record& right) {
    const int prime = left.prime;
    if (right.prime != prime) {
        throw std::runtime_error("compared records from different primes");
    }
    const int left_end = left.start + left.a + left.b + left.c;
    const int right_end = right.start + right.a + right.b + right.c;
    const int reflected_left_start = prime - 1 - left_end;
    const int reflected_left_end = prime - 1 - left.start;
    const int reflected_right_start = prime - 1 - right_end;
    const int reflected_right_end = prime - 1 - right.start;
    return
        closed_intervals_meet(
            left.start, left_end, right.start, right_end) ||
        closed_intervals_meet(
            left.start, left_end,
            reflected_right_start, reflected_right_end) ||
        closed_intervals_meet(
            reflected_left_start, reflected_left_end,
            right.start, right_end) ||
        closed_intervals_meet(
            reflected_left_start, reflected_left_end,
            reflected_right_start, reflected_right_end);
}

void audit_reflection_orbit_overlap() {
    const Record left{101, 10, 2, 2, 2, 5, false};
    const Record overlapping{101, 15, 2, 2, 2, 5, false};
    const Record separated{101, 30, 2, 2, 2, 5, false};
    if (!reflection_orbits_overlap(left, overlapping)) {
        throw std::runtime_error("synthetic quotient-overlap regression failed");
    }
    if (reflection_orbits_overlap(left, separated)) {
        throw std::runtime_error("synthetic quotient-separation regression failed");
    }
}

void scan_prime(int prime, Summary& summary) {
    std::vector<int> inverse(prime, 0);
    inverse[1] = 1;
    for (int value = 2; value < prime; ++value) {
        inverse[value] = static_cast<int>(
            prime - mul_mod(prime / value, inverse[prime % value], prime));
    }

    std::vector<int> apery(prime, 0);
    std::vector<int> companion(prime, 0);
    apery[0] = 1;
    apery[1] = 5 % prime;
    companion[0] = 0;
    companion[1] = 1;
    for (int index = 1; index < prime - 1; ++index) {
        const i64 n = index;
        const i64 n2 = mul_mod(n, n, prime);
        const i64 n3 = mul_mod(n2, n, prime);
        i64 coefficient = 5;
        coefficient = add_mod(coefficient, 27 * n, prime);
        coefficient = add_mod(coefficient, 51 * n2, prime);
        coefficient = add_mod(coefficient, 34 * n3, prime);
        const i64 denominator_inverse = mul_mod(
            mul_mod(inverse[index + 1], inverse[index + 1], prime),
            inverse[index + 1], prime);

        auto advance = [&](const std::vector<int>& values) {
            i64 numerator = mul_mod(coefficient, values[index], prime);
            numerator = add_mod(
                numerator, -mul_mod(n3, values[index - 1], prime), prime);
            return static_cast<int>(mul_mod(numerator, denominator_inverse, prime));
        };
        apery[index + 1] = advance(apery);
        companion[index + 1] = advance(companion);
    }

    // p does not divide any factorial j! with j < p.  Thus the prefix below
    // is exactly the p-support of U_s = product_{j<=s} j! b_j V_{j+1}.
    std::vector<unsigned char> carrier_bad(prime, 0);
    i64 v_previous = 0;  // V_0
    i64 v_current = 1;   // V_1
    for (int span = 1; span < prime; ++span) {
        const i64 v_next = add_mod(34 * v_current, -v_previous, prime);
        v_previous = v_current;
        v_current = v_next;  // V_{span+1}
        carrier_bad[span] = static_cast<unsigned char>(
            carrier_bad[span - 1] || apery[span] == 0 || v_current == 0);
    }

    // The value prime represents the point at infinity.  For each projective
    // state retain its three most recent occurrences; a new occurrence then
    // exposes exactly one sliding primitive four-return chain.
    std::vector<unsigned char> count(prime + 1, 0);
    std::vector<int> last0(prime + 1, 0);
    std::vector<int> last1(prime + 1, 0);
    std::vector<int> last2(prime + 1, 0);
    std::vector<int> short_count(prime + 1, 0);
    std::vector<Record> short_records;
    const int short_cutoff = integer_square_root(prime);

    for (int index = 0; index < prime; ++index) {
        if (apery[index] == 0 && companion[index] == 0) {
            throw std::runtime_error("projective solution pair vanished");
        }
        const int state = companion[index] == 0
            ? prime
            : static_cast<int>(mul_mod(
                  apery[index], inverse[companion[index]], prime));
        const unsigned char seen = count[state];
        if (seen >= 3) {
            const int x0 = last0[state];
            const int x1 = last1[state];
            const int x2 = last2[state];
            const int a = x1 - x0;
            const int b = x2 - x1;
            const int c = index - x2;
            const int span = index - x0;
            const bool valid_gaps = a >= 2 && b >= 2 && c >= 2;
            const bool non_progression = !(a == b && b == c);
            const bool off_center =
                x0 + x1 != prime - 1 &&
                x1 + x2 != prime - 1 &&
                x2 + index != prime - 1;
            if (valid_gaps && non_progression && off_center &&
                !carrier_bad[span]) {
                const Record record{
                    prime, x0, a, b, c, state, state == 0
                };
                consider_record(record, summary);
                if (span <= short_cutoff) {
                    ++short_count[state];
                    short_records.push_back(record);
                }
            }
        }

        if (seen == 0) {
            last0[state] = index;
            count[state] = 1;
        } else if (seen == 1) {
            last1[state] = index;
            count[state] = 2;
        } else if (seen == 2) {
            last2[state] = index;
            count[state] = 3;
        } else {
            last0[state] = last1[state];
            last1[state] = last2[state];
            last2[state] = index;
        }
    }

    std::uint64_t raw_mass = 0;
    std::uint64_t raw_energy = 0;
    std::uint64_t quotient_mass = 0;
    std::uint64_t quotient_energy = 0;
    std::uint64_t raw_collision_states = 0;
    std::uint64_t quotient_collision_states = 0;
    std::uint64_t quotient_overlap_energy = 0;
    std::uint64_t quotient_separated_energy = 0;
    int max_raw_count = 0;
    int max_quotient_count = 0;

    // Store exactly one representative of each reflection orbit.  Freeness
    // implies that the two reflected starts are distinct, so the smaller
    // start is a canonical choice.  Linked lists group representatives by
    // state without allocating one vector object per projective point.
    std::vector<Record> canonical_records;
    canonical_records.reserve(short_records.size() / 2);
    std::vector<int> canonical_head(prime + 1, -1);
    std::vector<int> canonical_next;
    canonical_next.reserve(short_records.size() / 2);
    for (const Record& record : short_records) {
        const int span = record.a + record.b + record.c;
        const int reflected_start = prime - 1 - record.start - span;
        if (record.start >= reflected_start) continue;
        const int record_index = static_cast<int>(canonical_records.size());
        canonical_records.push_back(record);
        canonical_next.push_back(canonical_head[record.state]);
        canonical_head[record.state] = record_index;
    }

    for (int state = 0; state <= prime; ++state) {
        const int raw_count = short_count[state];
        if (raw_count % 2 != 0) {
            throw std::runtime_error(
                "off-center short-chain count is not reflection-paired");
        }
        const int quotient_count = raw_count / 2;
        raw_mass += static_cast<std::uint64_t>(raw_count);
        raw_energy += static_cast<std::uint64_t>(raw_count) * raw_count;
        quotient_mass += static_cast<std::uint64_t>(quotient_count);
        quotient_energy +=
            static_cast<std::uint64_t>(quotient_count) * quotient_count;
        if (raw_count >= 2) ++raw_collision_states;
        if (quotient_count >= 2) ++quotient_collision_states;
        if (raw_count > max_raw_count) max_raw_count = raw_count;
        if (quotient_count > max_quotient_count) {
            max_quotient_count = quotient_count;
        }
        if (quotient_count >= 2 &&
            summary.quotient_collision_records.size() < 100) {
            for (const Record& record : short_records) {
                if (record.state == state &&
                    summary.quotient_collision_records.size() < 100) {
                    summary.quotient_collision_records.push_back(record);
                }
            }
        }

        for (int left = canonical_head[state]; left != -1;
             left = canonical_next[left]) {
            for (int right = canonical_next[left]; right != -1;
                 right = canonical_next[right]) {
                // Energies count ordered distinct pairs.
                if (reflection_orbits_overlap(
                        canonical_records[left], canonical_records[right])) {
                    quotient_overlap_energy += 2;
                } else {
                    quotient_separated_energy += 2;
                }
            }
        }
    }
    if (canonical_records.size() != quotient_mass) {
        throw std::runtime_error("reflection-orbit canonicalization failed");
    }
    if (quotient_energy != quotient_mass + quotient_overlap_energy +
            quotient_separated_energy) {
        throw std::runtime_error("quotient energy decomposition failed");
    }
    if (raw_mass > static_cast<std::uint64_t>(prime)) {
        throw std::runtime_error("short-chain mass exceeds number of starts");
    }
    const int bin_index = dyadic_index(prime);
    DplsBin& bin = summary.dpls_bins[bin_index];
    ++bin.prime_count;
    bin.raw_mass += raw_mass;
    bin.raw_energy += raw_energy;
    bin.quotient_mass += quotient_mass;
    bin.quotient_energy += quotient_energy;
    bin.quotient_overlap_energy += quotient_overlap_energy;
    bin.quotient_separated_energy += quotient_separated_energy;
    bin.raw_collision_states += raw_collision_states;
    bin.quotient_collision_states += quotient_collision_states;
    bin.actual_raw += static_cast<std::uint64_t>(short_count[0]);
    bin.raw_variance += static_cast<long double>(raw_energy) -
        static_cast<long double>(raw_mass) * raw_mass / (prime + 1);
    bin.quotient_variance += static_cast<long double>(quotient_energy) -
        static_cast<long double>(quotient_mass) * quotient_mass / (prime + 1);
    if (max_raw_count > bin.max_raw_count) {
        bin.max_raw_count = max_raw_count;
    }
    if (max_quotient_count > bin.max_quotient_count) {
        bin.max_quotient_count = max_quotient_count;
    }
    if (prime == 1297) {
        summary.saw_dpls_1297_pair =
            short_count[454] == 2 && max_raw_count == 2 &&
            max_quotient_count == 1;
    }
}

int parse_limit(int argc, char** argv) {
    int limit = 10000;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--limit" && index + 1 < argc) {
            limit = std::stoi(argv[++index]);
        } else {
            throw std::runtime_error("usage: primitive_projective_prime_scan --limit N");
        }
    }
    if (limit < 7) throw std::runtime_error("limit must be at least 7");
    return limit;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        audit_reflection_orbit_overlap();
        const int limit = parse_limit(argc, argv);
        const std::vector<int> primes = primes_up_to(limit);
        Summary summary;
        for (std::size_t index = 0; index < primes.size(); ++index) {
            scan_prime(primes[index], summary);
            if ((index + 1) % 1000 == 0) {
                std::cerr << "scanned " << index + 1 << "/" << primes.size()
                          << " primes through " << primes[index] << '\n';
            }
        }

        if (limit >= 1297 &&
            !(summary.saw_1297_forward && summary.saw_1297_reverse)) {
            throw std::runtime_error("missing the canonical p=1297 regression pair");
        }
        if (limit >= 1297 && !summary.saw_dpls_1297_pair) {
            throw std::runtime_error(
                "p=1297 raw/quotient DPLS regression failed");
        }
        if (limit >= 7411 &&
            !(summary.saw_7411_forward && summary.saw_7411_reverse)) {
            throw std::runtime_error("missing the canonical p=7411 regression pair");
        }
        if (limit >= 128047 &&
            !(summary.saw_128047_forward && summary.saw_128047_reverse)) {
            throw std::runtime_error("missing the canonical p=128047 regression pair");
        }

        std::cout << "limit=" << limit << " primes=" << primes.size()
                  << " primitive_off_center=" << summary.primitive_off_center
                  << " ratio_above_one=" << summary.ratio_above_one
                  << " ratio_above_two=" << summary.ratio_above_two << '\n';
        if (summary.best_prime != 0) {
            const Record& record = summary.best_record;
            const int span = record.a + record.b + record.c;
            const i64 span_squared = static_cast<i64>(span) * span;
            std::cout << "max_p_over_span_squared=" << summary.best_prime << "/"
                      << span_squared << " approx="
                      << static_cast<double>(summary.best_prime) / span_squared
                      << " record=(" << record.prime << ',' << record.start << ','
                      << record.a << ',' << record.b << ',' << record.c << ','
                      << record.state << ',' << (record.actual ? "actual" : "phantom")
                      << ")\n";
        }
        for (const Record& record : summary.records_above_one) {
            const int span = record.a + record.b + record.c;
            std::cout << "record p=" << record.prime << " x=" << record.start
                      << " gaps=" << record.a << ',' << record.b << ',' << record.c
                      << " span=" << span << " state=" << record.state
                      << " orbit=" << (record.actual ? "actual" : "phantom") << '\n';
        }
        std::cout << std::setprecision(12);
        DplsBin aggregate;
        for (const DplsBin& bin : summary.dpls_bins) {
            aggregate.prime_count += bin.prime_count;
            aggregate.raw_mass += bin.raw_mass;
            aggregate.raw_energy += bin.raw_energy;
            aggregate.quotient_mass += bin.quotient_mass;
            aggregate.quotient_energy += bin.quotient_energy;
            aggregate.quotient_overlap_energy +=
                bin.quotient_overlap_energy;
            aggregate.quotient_separated_energy +=
                bin.quotient_separated_energy;
            aggregate.raw_collision_states += bin.raw_collision_states;
            aggregate.quotient_collision_states +=
                bin.quotient_collision_states;
            aggregate.actual_raw += bin.actual_raw;
            aggregate.raw_variance += bin.raw_variance;
            aggregate.quotient_variance += bin.quotient_variance;
            if (bin.max_raw_count > aggregate.max_raw_count) {
                aggregate.max_raw_count = bin.max_raw_count;
            }
            if (bin.max_quotient_count > aggregate.max_quotient_count) {
                aggregate.max_quotient_count = bin.max_quotient_count;
            }
        }
        std::cout << "dpls_aggregate"
                  << " primes=" << aggregate.prime_count
                  << " raw_M=" << aggregate.raw_mass
                  << " raw_E=" << aggregate.raw_energy
                  << " raw_V="
                  << static_cast<double>(aggregate.raw_variance)
                  << " raw_maxC=" << aggregate.max_raw_count
                  << " raw_collision_states="
                  << aggregate.raw_collision_states
                  << " quotient_M=" << aggregate.quotient_mass
                  << " quotient_E=" << aggregate.quotient_energy
                  << " quotient_E_overlap="
                  << aggregate.quotient_overlap_energy
                  << " quotient_E_separated="
                  << aggregate.quotient_separated_energy
                  << " quotient_V="
                  << static_cast<double>(aggregate.quotient_variance)
                  << " quotient_maxC="
                  << aggregate.max_quotient_count
                  << " quotient_collision_states="
                  << aggregate.quotient_collision_states
                  << " actual_raw=" << aggregate.actual_raw << '\n';
        for (std::size_t index = 0; index < summary.dpls_bins.size(); ++index) {
            const DplsBin& bin = summary.dpls_bins[index];
            if (bin.prime_count == 0) continue;
            const std::uint64_t lower = std::uint64_t{1} << index;
            const std::uint64_t upper =
                (std::uint64_t{1} << (index + 1)) - 1;
            const long double raw_ratio = bin.raw_mass == 0
                ? 0.0L : bin.raw_variance / bin.raw_mass;
            const long double quotient_ratio = bin.quotient_mass == 0
                ? 0.0L : bin.quotient_variance / bin.quotient_mass;
            std::cout
                << "dpls range=[" << lower << ',' << upper << "]"
                << " primes=" << bin.prime_count
                << " raw_M=" << bin.raw_mass
                << " raw_E=" << bin.raw_energy
                << " raw_V=" << static_cast<double>(bin.raw_variance)
                << " raw_V_over_M=" << static_cast<double>(raw_ratio)
                << " raw_maxC=" << bin.max_raw_count
                << " raw_collision_states=" << bin.raw_collision_states
                << " quotient_M=" << bin.quotient_mass
                << " quotient_E=" << bin.quotient_energy
                << " quotient_E_overlap=" << bin.quotient_overlap_energy
                << " quotient_E_separated="
                << bin.quotient_separated_energy
                << " quotient_V="
                << static_cast<double>(bin.quotient_variance)
                << " quotient_V_over_M="
                << static_cast<double>(quotient_ratio)
                << " quotient_maxC=" << bin.max_quotient_count
                << " quotient_collision_states="
                << bin.quotient_collision_states
                << " actual_raw=" << bin.actual_raw << '\n';
        }
        for (const Record& record : summary.quotient_collision_records) {
            std::cout << "dpls_quotient_collision p=" << record.prime
                      << " state=" << record.state
                      << " x=" << record.start
                      << " gaps=" << record.a << ',' << record.b << ','
                      << record.c
                      << " span=" << record.a + record.b + record.c
                      << " orbit="
                      << (record.actual ? "actual" : "phantom") << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
