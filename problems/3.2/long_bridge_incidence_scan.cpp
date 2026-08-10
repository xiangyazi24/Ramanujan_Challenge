// Exact dyadic scan of two-equation Apéry long-bridge incidences.
//
// For every prime X < p <= 2X, group the projective orbit positions by
// state.  In one occurrence list r_0 < ... < r_m, a tuple
//
//     x = r_j,  y = r_i,  z = r_k
//
// with y-x=s and z-y=G is exactly an F_p-root witness for
//
//     N_s(x) = N_G(x+s) = 0.
//
// The scan counts all such tuples with 2 <= s <= floor(sqrt(p)), then the
// subset whose first four occurrences are consecutive with gaps >= 2, and
// finally the subset whose first four-occurrence chain passes the actual
// off-center/non-AP/carrier filters.  The latter two still do not require a
// second selected chain, so they are explicit upper-bound incidences rather
// than the canonical separated energy itself.  For selected first chains,
// the scan also removes every later point in that chain's own reflection
// orbit; a distinct separated quotient chain cannot start at such a point.
// A final column retains only external points that themselves begin a
// four-consecutive-occurrence window of span at most floor(sqrt(p)), a
// necessary condition for the second selected chain.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using i128 = __int128_t;
using u64 = std::uint64_t;

i64 mul_mod(i64 left, i64 right, i64 prime) {
    return static_cast<i64>((static_cast<i128>(left) * right) % prime);
}

i64 add_mod(i64 value, i64 term, i64 prime) {
    value += term;
    value %= prime;
    if (value < 0) value += prime;
    return value;
}

int integer_square_root(int value) {
    int root = static_cast<int>(std::sqrt(static_cast<double>(value)));
    while (static_cast<i64>(root + 1) * (root + 1) <= value) ++root;
    while (static_cast<i64>(root) * root > value) --root;
    return root;
}

std::vector<int> primes_in_dyadic_interval(int scale) {
    const int limit = 2 * scale;
    std::vector<bool> sieve(limit + 1, true);
    sieve[0] = false;
    sieve[1] = false;
    for (int prime = 2; static_cast<i64>(prime) * prime <= limit; ++prime) {
        if (!sieve[prime]) continue;
        for (int multiple = prime * prime; multiple <= limit; multiple += prime) {
            sieve[multiple] = false;
        }
    }
    std::vector<int> primes;
    for (int value = std::max(7, scale + 1); value <= limit; ++value) {
        if (sieve[value]) primes.push_back(value);
    }
    return primes;
}

struct BridgeBin {
    int lower;
    int upper;
    u64 all_short = 0;
    u64 four_gap_valid = 0;
    u64 selected_first = 0;
    u64 selected_external = 0;
    u64 selected_to_external_short_four = 0;
};

struct PrimeStats {
    u64 short_return_pairs = 0;
    u64 four_gap_valid_chains = 0;
    u64 selected_first_chains = 0;
    std::vector<BridgeBin> bins;
};

struct Options {
    int scale = 5000;
    int max_primes = 0;
    bool verbose = false;
};

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--x" && index + 1 < argc) {
            options.scale = std::stoi(argv[++index]);
        } else if (argument == "--max-primes" && index + 1 < argc) {
            options.max_primes = std::stoi(argv[++index]);
        } else if (argument == "--verbose") {
            options.verbose = true;
        } else {
            throw std::runtime_error(
                "usage: long_bridge_incidence_scan --x X "
                "[--max-primes N] [--verbose]");
        }
    }
    if (options.scale < 7) throw std::runtime_error("X must be at least 7");
    if (options.scale > 1000000000 / 2) {
        throw std::runtime_error("X is too large for 32-bit indexing");
    }
    if (options.max_primes < 0) {
        throw std::runtime_error("max-primes must be nonnegative");
    }
    return options;
}

std::vector<BridgeBin> make_bins(int scale) {
    const int maximum = 2 * scale;
    const int height = integer_square_root(maximum);
    const int near_cutoff = std::max(
        2,
        static_cast<int>(height / std::sqrt(std::log(static_cast<double>(scale)))));
    std::vector<BridgeBin> bins;
    bins.push_back(BridgeBin{2, near_cutoff});
    if (near_cutoff < height) {
        bins.push_back(BridgeBin{near_cutoff + 1, height});
    }
    int lower = height + 1;
    while (lower <= maximum) {
        const i64 doubled = 2 * static_cast<i64>(lower) - 1;
        const int upper = static_cast<int>(
            std::min<i64>(maximum, doubled));
        bins.push_back(BridgeBin{lower, upper});
        lower = upper + 1;
    }
    return bins;
}

PrimeStats scan_prime(int prime, const std::vector<BridgeBin>& template_bins) {
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

    const int short_cutoff = integer_square_root(prime);
    std::vector<unsigned char> carrier_bad(short_cutoff + 1, 0);
    i64 v_previous = 0;
    i64 v_current = 1;
    for (int span = 1; span <= short_cutoff; ++span) {
        const i64 v_next = add_mod(34 * v_current, -v_previous, prime);
        v_previous = v_current;
        v_current = v_next;
        carrier_bad[span] = static_cast<unsigned char>(
            carrier_bad[span - 1] || apery[span] == 0 || v_current == 0);
    }

    std::vector<int> states(prime, 0);
    std::vector<std::vector<int>> occurrences(prime + 1);
    for (int index = 0; index < prime; ++index) {
        if (apery[index] == 0 && companion[index] == 0) {
            throw std::runtime_error("projective solution pair vanished");
        }
        const int state = companion[index] == 0
            ? prime
            : static_cast<int>(mul_mod(
                  apery[index], inverse[companion[index]], prime));
        states[index] = state;
        occurrences[state].push_back(index);
    }
    for (int index = 0; index < prime; ++index) {
        if (states[index] != states[prime - 1 - index]) {
            throw std::runtime_error("projective reflection regression failed");
        }
        if (index + 1 < prime && states[index] == states[index + 1]) {
            throw std::runtime_error("adjacent projective return regression failed");
        }
    }

    PrimeStats stats;
    stats.bins = template_bins;
    for (const std::vector<int>& list : occurrences) {
        std::vector<unsigned char> short_four_start(list.size(), 0);
        std::vector<u64> short_four_prefix(list.size() + 1, 0);
        for (std::size_t index = 0; index < list.size(); ++index) {
            short_four_start[index] = static_cast<unsigned char>(
                index + 3 < list.size() &&
                list[index + 3] - list[index] <= short_cutoff);
            short_four_prefix[index + 1] =
                short_four_prefix[index] + short_four_start[index];
        }
        for (std::size_t middle = 0; middle < list.size(); ++middle) {
            const int y = list[middle];
            const auto prior_begin = std::lower_bound(
                list.begin(), list.begin() + static_cast<std::ptrdiff_t>(middle),
                y - short_cutoff);
            const auto prior_end = std::upper_bound(
                prior_begin, list.begin() + static_cast<std::ptrdiff_t>(middle),
                y - 2);
            const u64 prior_count = static_cast<u64>(prior_end - prior_begin);
            stats.short_return_pairs += prior_count;

            bool four_gap_valid = false;
            bool selected_first = false;
            std::array<int, 4> reflected_support{};
            if (middle >= 3) {
                const int x0 = list[middle - 3];
                const int x1 = list[middle - 2];
                const int x2 = list[middle - 1];
                const int a = x1 - x0;
                const int b = x2 - x1;
                const int c = y - x2;
                const int span = y - x0;
                four_gap_valid =
                    span <= short_cutoff && a >= 2 && b >= 2 && c >= 2;
                selected_first = four_gap_valid && !(a == b && b == c) &&
                    x0 + x1 != prime - 1 && x1 + x2 != prime - 1 &&
                    x2 + y != prime - 1 && !carrier_bad[span];
                reflected_support = {
                    prime - 1 - y,
                    prime - 1 - x2,
                    prime - 1 - x1,
                    prime - 1 - x0,
                };
            }
            stats.four_gap_valid_chains += static_cast<u64>(four_gap_valid);
            stats.selected_first_chains += static_cast<u64>(selected_first);

            for (BridgeBin& bin : stats.bins) {
                const auto later_begin = std::lower_bound(
                    list.begin() + static_cast<std::ptrdiff_t>(middle + 1),
                    list.end(), y + bin.lower);
                const auto later_end = std::upper_bound(
                    later_begin, list.end(), y + bin.upper);
                const u64 later_count = static_cast<u64>(later_end - later_begin);
                bin.all_short += prior_count * later_count;
                if (four_gap_valid) bin.four_gap_valid += later_count;
                if (selected_first) {
                    bin.selected_first += later_count;
                    u64 external_count = later_count;
                    const std::size_t later_begin_index = static_cast<std::size_t>(
                        later_begin - list.begin());
                    const std::size_t later_end_index = static_cast<std::size_t>(
                        later_end - list.begin());
                    u64 external_short_four =
                        short_four_prefix[later_end_index] -
                        short_four_prefix[later_begin_index];
                    for (const int reflected : reflected_support) {
                        const int gap = reflected - y;
                        if (bin.lower <= gap && gap <= bin.upper) {
                            if (external_count == 0) {
                                throw std::runtime_error(
                                    "reflected-support subtraction underflow");
                            }
                            --external_count;
                            const auto reflected_position = std::lower_bound(
                                list.begin(), list.end(), reflected);
                            if (reflected_position == list.end() ||
                                *reflected_position != reflected) {
                                throw std::runtime_error(
                                    "reflected support left its projective fiber");
                            }
                            const std::size_t reflected_index =
                                static_cast<std::size_t>(
                                    reflected_position - list.begin());
                            if (short_four_start[reflected_index]) {
                                if (external_short_four == 0) {
                                    throw std::runtime_error(
                                        "short-four subtraction underflow");
                                }
                                --external_short_four;
                            }
                        }
                    }
                    bin.selected_external += external_count;
                    bin.selected_to_external_short_four += external_short_four;
                }
            }
        }
    }
    return stats;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        std::vector<int> primes = primes_in_dyadic_interval(options.scale);
        if (options.max_primes > 0 &&
            static_cast<int>(primes.size()) > options.max_primes) {
            primes.resize(static_cast<std::size_t>(options.max_primes));
        }
        const std::vector<BridgeBin> template_bins = make_bins(options.scale);
        PrimeStats aggregate;
        aggregate.bins = template_bins;

        for (std::size_t index = 0; index < primes.size(); ++index) {
            const PrimeStats stats = scan_prime(primes[index], template_bins);
            aggregate.short_return_pairs += stats.short_return_pairs;
            aggregate.four_gap_valid_chains += stats.four_gap_valid_chains;
            aggregate.selected_first_chains += stats.selected_first_chains;
            for (std::size_t bin = 0; bin < aggregate.bins.size(); ++bin) {
                aggregate.bins[bin].all_short += stats.bins[bin].all_short;
                aggregate.bins[bin].four_gap_valid +=
                    stats.bins[bin].four_gap_valid;
                aggregate.bins[bin].selected_first +=
                    stats.bins[bin].selected_first;
                aggregate.bins[bin].selected_external +=
                    stats.bins[bin].selected_external;
                aggregate.bins[bin].selected_to_external_short_four +=
                    stats.bins[bin].selected_to_external_short_four;
            }
            if (options.verbose &&
                (stats.selected_first_chains > 0 || index < 10)) {
                std::cout << "prime=" << primes[index]
                          << " short_pairs=" << stats.short_return_pairs
                          << " four_chains=" << stats.four_gap_valid_chains
                          << " selected_chains=" << stats.selected_first_chains
                          << '\n';
            }
            if ((index + 1) % 500 == 0) {
                std::cerr << "scanned " << index + 1 << "/" << primes.size()
                          << " primes through " << primes[index] << '\n';
            }
        }

        const int height = integer_square_root(2 * options.scale);
        const int near_cutoff = template_bins.front().upper;
        std::cout << "LONG_BRIDGE_INCIDENCE"
                  << " X=" << options.scale
                  << " primes=" << primes.size()
                  << " H_global=" << height
                  << " K=" << near_cutoff
                  << " short_pairs=" << aggregate.short_return_pairs
                  << " four_chains=" << aggregate.four_gap_valid_chains
                  << " selected_chains=" << aggregate.selected_first_chains
                  << '\n';
        for (const BridgeBin& bin : aggregate.bins) {
            std::cout << "bridge_bin=[" << bin.lower << ',' << bin.upper << ']'
                      << " all_short=" << bin.all_short
                      << " four_gap_valid=" << bin.four_gap_valid
                      << " selected_first=" << bin.selected_first
                      << " selected_external=" << bin.selected_external
                      << " selected_to_external_short_four="
                      << bin.selected_to_external_short_four << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
