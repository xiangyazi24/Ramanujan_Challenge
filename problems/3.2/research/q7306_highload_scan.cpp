// Q7306: exact high-load invariant scanner for Apéry zero sets.
//
// Primary enumeration: pair-CRT records followed by fail-closed clique
// decoding.  The X=256 run is independently checked by direct hit scatter.
// Output: deterministic JSON on stdout and a compact exact summary on stderr.
//
// Build:
//   g++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic \
//       problems/3.2/research/q7306_highload_scan.cpp -o q7306_scan
// Run:
//   ./q7306_scan 256 512 1024 2048 > q7306-output.json \
//       2> q7306-summary.txt

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using u8 = std::uint8_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

static void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

static std::string u128_decimal(u128 value) {
    if (value == 0) return "0";
    std::string out;
    while (value != 0) {
        out.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

static inline u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

static u64 pow_mod(u64 a, u64 exponent, u64 mod) {
    u64 result = 1 % mod;
    a %= mod;
    while (exponent != 0) {
        if (exponent & 1U) result = mul_mod(result, a, mod);
        a = mul_mod(a, a, mod);
        exponent >>= 1U;
    }
    return result;
}

static inline u64 sub_mod(u64 a, u64 b, u64 mod) {
    return a >= b ? a - b : a + mod - b;
}

static u64 P_mod(u64 m, u64 mod) {
    const u64 x = m % mod;
    const u64 x2 = mul_mod(x, x, mod);
    const u64 x3 = mul_mod(x2, x, mod);
    u64 value = mul_mod(34 % mod, x3, mod);
    value = (value + mul_mod(51 % mod, x2, mod)) % mod;
    value = (value + mul_mod(27 % mod, x, mod)) % mod;
    return (value + 5) % mod;
}

static inline u64 cube_mod(u64 value, u64 mod) {
    return mul_mod(mul_mod(value % mod, value % mod, mod), value % mod, mod);
}

static inline u64 sixth_mod(u64 value, u64 mod) {
    const u64 cube = cube_mod(value, mod);
    return mul_mod(cube, cube, mod);
}

static std::vector<u32> primes_between(u32 X) {
    const u32 limit = 2 * X;
    std::vector<bool> is_prime(limit + 1, true);
    if (limit >= 0) is_prime[0] = false;
    if (limit >= 1) is_prime[1] = false;
    for (u32 p = 2; static_cast<u64>(p) * p <= limit; ++p) {
        if (!is_prime[p]) continue;
        for (u32 multiple = p * p; multiple <= limit; multiple += p)
            is_prime[multiple] = false;
    }
    std::vector<u32> primes;
    for (u32 p = X + 1; p <= limit; ++p)
        if (is_prime[p]) primes.push_back(p);
    return primes;
}

struct BuildChecks {
    u64 primes_checked = 0;
    u64 cleared_divided_positions = 0;
    u64 direct_binomial_positions = 0;
    u64 reflection_positions = 0;
    u64 reflection_zero_positions = 0;
    u64 wronskian_positions = 0;
    u64 nonadjacent_zero_positions = 0;
};

struct PrimeData {
    u32 p = 0;
    std::vector<u32> B;
    std::vector<u32> D;
    std::vector<u32> b;
    std::vector<u32> factorial;
    std::vector<u32> zeros;
};

static u32 apery_direct_mod(u32 n, u32 p,
                            const std::vector<u32> &inverse) {
    u64 choose_n_j = 1;
    u64 choose_n_plus_j_j = 1;
    u64 sum = 0;
    for (u32 j = 0; j <= n; ++j) {
        const u64 a2 = mul_mod(choose_n_j, choose_n_j, p);
        const u64 c2 = mul_mod(choose_n_plus_j_j,
                               choose_n_plus_j_j, p);
        sum = (sum + mul_mod(a2, c2, p)) % p;
        if (j == n) break;
        const u64 inv = inverse[j + 1];
        choose_n_j = mul_mod(choose_n_j, n - j, p);
        choose_n_j = mul_mod(choose_n_j, inv, p);
        choose_n_plus_j_j = mul_mod(choose_n_plus_j_j,
                                    n + j + 1, p);
        choose_n_plus_j_j = mul_mod(choose_n_plus_j_j, inv, p);
    }
    return static_cast<u32>(sum);
}

static PrimeData build_prime_data(u32 p, BuildChecks &checks) {
    require(p >= 7, "prime kernel expects p>=7");
    PrimeData out;
    out.p = p;
    out.B.assign(p, 0);
    out.D.assign(p, 0);
    out.b.assign(p, 0);
    out.factorial.assign(p, 0);

    out.B[0] = 1 % p;
    out.B[1] = 5 % p;
    out.D[0] = 0;
    out.D[1] = 1;
    for (u32 m = 1; m + 1 < p; ++m) {
        out.B[m + 1] = static_cast<u32>(sub_mod(
            mul_mod(P_mod(m, p), out.B[m], p),
            mul_mod(sixth_mod(m, p), out.B[m - 1], p), p));
        out.D[m + 1] = static_cast<u32>(sub_mod(
            mul_mod(P_mod(m, p), out.D[m], p),
            mul_mod(sixth_mod(m, p), out.D[m - 1], p), p));
    }

    std::vector<u32> inverse(p, 0);
    inverse[1] = 1;
    for (u32 m = 2; m < p; ++m) {
        inverse[m] = static_cast<u32>(
            p - mul_mod(p / m, inverse[p % m], p));
        if (inverse[m] == p) inverse[m] = 0;
        require(mul_mod(m, inverse[m], p) == 1,
                "modular inverse precomputation failed");
    }

    out.b[0] = 1 % p;
    out.b[1] = 5 % p;
    for (u32 m = 1; m + 1 < p; ++m) {
        const u64 numerator = sub_mod(
            mul_mod(P_mod(m, p), out.b[m], p),
            mul_mod(cube_mod(m, p), out.b[m - 1], p), p);
        const u64 inv_cube = cube_mod(inverse[m + 1], p);
        out.b[m + 1] = static_cast<u32>(mul_mod(numerator,
                                                        inv_cube, p));
    }

    out.factorial[0] = 1;
    for (u32 m = 1; m < p; ++m)
        out.factorial[m] = static_cast<u32>(
            mul_mod(out.factorial[m - 1], m, p));

    for (u32 m = 0; m < p; ++m) {
        const u64 factor_cube = cube_mod(out.factorial[m], p);
        const u64 expected_B = mul_mod(factor_cube, out.b[m], p);
        if (out.B[m] != expected_B) {
            std::ostringstream message;
            message << "cleared/divided recurrence mismatch p=" << p
                    << " m=" << m << " B=" << out.B[m]
                    << " expected=" << expected_B;
            throw std::runtime_error(message.str());
        }
        ++checks.cleared_divided_positions;
        require((out.B[m] == 0) == (out.b[m] == 0),
                "zero-set scaling mismatch");
        if (out.b[m] == 0) out.zeros.push_back(m);
    }

    const u32 direct_stop = std::min<u32>(12, p - 1);
    for (u32 n = 0; n <= direct_stop; ++n) {
        const u32 direct = apery_direct_mod(n, p, inverse);
        if (direct != out.b[n]) {
            std::ostringstream message;
            message << "direct binomial mismatch p=" << p
                    << " n=" << n << " direct=" << direct
                    << " recurrence=" << out.b[n];
            throw std::runtime_error(message.str());
        }
        ++checks.direct_binomial_positions;
    }

    for (u32 m = 0; m < p; ++m) {
        if (out.b[m] != out.b[p - 1 - m]) {
            std::ostringstream message;
            message << "reflection value mismatch p=" << p
                    << " m=" << m;
            throw std::runtime_error(message.str());
        }
        ++checks.reflection_positions;
    }

    for (u32 r : out.zeros) {
        require(std::binary_search(out.zeros.begin(), out.zeros.end(),
                                   p - 1 - r),
                "reflection zero-set closure failed");
        ++checks.reflection_zero_positions;
        if (r > 0) {
            require(out.b[r - 1] != 0, "adjacent Apéry zeros detected");
            ++checks.nonadjacent_zero_positions;
        }
        if (r + 1 < p) {
            require(out.b[r + 1] != 0, "adjacent Apéry zeros detected");
            ++checks.nonadjacent_zero_positions;
        }
    }

    for (u32 m = 0; m + 1 < p; ++m) {
        const u64 lhs = sub_mod(
            mul_mod(out.B[m], out.D[m + 1], p),
            mul_mod(out.B[m + 1], out.D[m], p), p);
        const u64 rhs = sixth_mod(out.factorial[m], p);
        if (lhs != rhs) {
            std::ostringstream message;
            message << "Wronskian mismatch p=" << p << " m=" << m
                    << " lhs=" << lhs << " rhs=" << rhs;
            throw std::runtime_error(message.str());
        }
        ++checks.wronskian_positions;
    }

    ++checks.primes_checked;
    return out;
}

struct PairRec {
    u64 n = 0;
    u32 a = 0;
    u32 b = 0;
};

struct RowBase {
    u64 n = 0;
    std::vector<u32> ids;
};

static std::vector<PairRec>
enumerate_pair_records(const std::vector<PrimeData> &prime_data,
                       u64 limit) {
    std::vector<PairRec> records;
    for (u32 i = 0; i < prime_data.size(); ++i) {
        if (prime_data[i].zeros.empty()) continue;
        for (u32 j = i + 1; j < prime_data.size(); ++j) {
            if (prime_data[j].zeros.empty()) continue;
            const u64 p = prime_data[i].p;
            const u64 q = prime_data[j].p;
            const u64 p_inverse = pow_mod(p % q, q - 2, q);
            require(mul_mod(p % q, p_inverse, q) == 1,
                    "CRT inverse check failed");
            for (u32 rp : prime_data[i].zeros) {
                for (u32 rq : prime_data[j].zeros) {
                    const u64 difference = (rq + q - (rp % q)) % q;
                    const u64 multiplier = mul_mod(difference,
                                                   p_inverse, q);
                    const u64 n = static_cast<u64>(rp) + p * multiplier;
                    require(n < p * q, "CRT representative out of range");
                    if (n < limit) records.push_back({n, i, j});
                }
            }
        }
    }
    std::sort(records.begin(), records.end(),
              [](const PairRec &left, const PairRec &right) {
                  return std::tie(left.n, left.a, left.b) <
                         std::tie(right.n, right.a, right.b);
              });
    return records;
}

static std::vector<u32>
validate_clique_edges(const std::vector<std::pair<u32, u32>> &raw_edges,
                      u32 vertex_bound) {
    std::set<std::pair<u32, u32>> edges;
    std::vector<u32> ids;
    for (const auto &edge : raw_edges) {
        require(edge.first < edge.second, "noncanonical clique edge");
        require(edge.second < vertex_bound, "clique vertex out of range");
        require(edges.insert(edge).second, "duplicate clique edge");
        ids.push_back(edge.first);
        ids.push_back(edge.second);
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    const u64 K = ids.size();
    require(edges.size() == K * (K - 1) / 2,
            "equal-n pair group is not a complete clique");
    for (u32 i = 0; i < ids.size(); ++i)
        for (u32 j = i + 1; j < ids.size(); ++j)
            require(edges.count({ids[i], ids[j]}) == 1,
                    "clique edge missing despite count check");
    return ids;
}

static bool synthetic_clique_self_check() {
    const std::vector<std::pair<u32, u32>> triangle = {
        {0, 1}, {0, 2}, {1, 2}};
    require(validate_clique_edges(triangle, 10) ==
                std::vector<u32>({0, 1, 2}),
            "synthetic triangle decoding failed");

    std::vector<std::pair<u32, u32>> K4;
    for (u32 i = 0; i < 4; ++i)
        for (u32 j = i + 1; j < 4; ++j) K4.push_back({i, j});
    require(validate_clique_edges(K4, 10) ==
                std::vector<u32>({0, 1, 2, 3}),
            "synthetic K4 decoding failed");

    bool missing_rejected = false;
    std::vector<std::pair<u32, u32>> missing = K4;
    missing.pop_back();
    try {
        (void)validate_clique_edges(missing, 10);
    } catch (const std::runtime_error &) {
        missing_rejected = true;
    }
    require(missing_rejected, "incomplete synthetic clique was accepted");

    bool duplicate_rejected = false;
    std::vector<std::pair<u32, u32>> duplicate = triangle;
    duplicate.push_back({0, 1});
    try {
        (void)validate_clique_edges(duplicate, 10);
    } catch (const std::runtime_error &) {
        duplicate_rejected = true;
    }
    require(duplicate_rejected, "duplicate synthetic edge was accepted");
    return true;
}

struct DecodeResult {
    std::vector<RowBase> rows;
    u64 groups_checked = 0;
    u64 row_memberships_checked = 0;
};

static DecodeResult decode_rows(const std::vector<PairRec> &records,
                                const std::vector<PrimeData> &prime_data) {
    DecodeResult result;
    for (std::size_t lo = 0; lo < records.size();) {
        std::size_t hi = lo + 1;
        while (hi < records.size() && records[hi].n == records[lo].n) ++hi;
        std::vector<std::pair<u32, u32>> edges;
        edges.reserve(hi - lo);
        for (std::size_t cursor = lo; cursor < hi; ++cursor)
            edges.push_back({records[cursor].a, records[cursor].b});
        std::vector<u32> ids = validate_clique_edges(
            edges, static_cast<u32>(prime_data.size()));
        ++result.groups_checked;
        for (u32 id : ids) {
            const u32 residue = static_cast<u32>(records[lo].n %
                                                  prime_data[id].p);
            require(std::binary_search(prime_data[id].zeros.begin(),
                                       prime_data[id].zeros.end(), residue),
                    "decoded row fails zero-set membership");
            ++result.row_memberships_checked;
        }
        if (ids.size() >= 3)
            result.rows.push_back({records[lo].n, std::move(ids)});
        lo = hi;
    }
    return result;
}

static inline u64 choose2(u64 K) {
    return K < 2 ? 0 : K * (K - 1) / 2;
}

static inline u64 choose3(u64 K) {
    return K < 3 ? 0 : K * (K - 1) * (K - 2) / 6;
}

struct ScatterOutcome {
    bool performed = false;
    bool passed = false;
    u64 hit_events = 0;
    u64 rows = 0;
    u64 pair_records = 0;
    u64 canonical_triples = 0;
};

static ScatterOutcome direct_scatter_check(
    u32 X, const std::vector<PrimeData> &prime_data,
    const std::vector<RowBase> &decoded_rows, u64 expected_pair_records) {
    ScatterOutcome out;
    out.performed = true;
    const u64 limit = static_cast<u64>(X) * X;
    std::vector<u8> load(limit, 0);
    for (const PrimeData &item : prime_data) {
        for (u32 residue : item.zeros) {
            for (u64 n = residue; n < limit; n += item.p) {
                require(load[n] != std::numeric_limits<u8>::max(),
                        "scatter load overflow");
                ++load[n];
                ++out.hit_events;
            }
        }
    }

    std::size_t row_cursor = 0;
    for (u64 n = 0; n < limit; ++n) {
        out.pair_records += choose2(load[n]);
        out.canonical_triples += choose3(load[n]);
        if (load[n] < 3) continue;
        require(row_cursor < decoded_rows.size(),
                "scatter found row absent from clique decoder");
        require(decoded_rows[row_cursor].n == n,
                "scatter/clique row index mismatch");
        require(decoded_rows[row_cursor].ids.size() == load[n],
                "scatter/clique row load mismatch");
        ++row_cursor;
        ++out.rows;
    }
    require(row_cursor == decoded_rows.size(),
            "clique decoder found row absent from scatter");
    require(out.pair_records == expected_pair_records,
            "scatter pair count mismatch");
    u64 decoded_triples = 0;
    for (const RowBase &row : decoded_rows)
        decoded_triples += choose3(row.ids.size());
    require(out.canonical_triples == decoded_triples,
            "scatter canonical triple count mismatch");
    out.passed = true;
    return out;
}

struct Baseline {
    u32 prime_count = 0;
    u32 active_count = 0;
    u32 zero_count = 0;
    u64 pair_records = 0;
    std::map<u32, u64> highload_histogram;
};

static bool baseline_for(u32 X, Baseline &out) {
    if (X == 256) {
        out = {43, 17, 46, 467, {{3, 10}}};
        return true;
    }
    if (X == 512) {
        out = {75, 30, 70, 1095, {{3, 25}}};
        return true;
    }
    if (X == 1024) {
        out = {137, 47, 120, 3447, {{3, 77}, {4, 1}}};
        return true;
    }
    if (X == 2048) {
        out = {255, 94, 231, 12405, {{3, 272}, {4, 8}}};
        return true;
    }
    return false;
}

static u32 continuant_mod(u32 x, u32 h, u32 p) {
    if (h == 0) return 0;
    if (h == 1) return 1 % p;
    u64 previous = 0;
    u64 current = 1 % p;
    for (u32 k = 1; k < h; ++k) {
        const u64 m = static_cast<u64>(x) + k;
        const u64 next = sub_mod(
            mul_mod(P_mod(m, p), current, p),
            mul_mod(sixth_mod(m, p), previous, p), p);
        previous = current;
        current = next;
    }
    return static_cast<u32>(current);
}

static u32 state_minor(const PrimeData &item, u32 s, u32 t) {
    return static_cast<u32>(sub_mod(
        mul_mod(item.B[s], item.D[t], item.p),
        mul_mod(item.D[s], item.B[t], item.p), item.p));
}

struct Hit {
    u32 id = 0;
    u32 p = 0;
    u32 q = 0;
    u32 r = 0;
    u32 rho = 0;
    u32 zero_count = 0;
};

struct Digit {
    u32 value = 0;
    std::vector<u32> q_origins;
    std::vector<u32> rho_origins;
};

static std::vector<std::string> edge_types(const Digit &left,
                                           const Digit &right) {
    std::vector<std::string> types;
    if (!left.q_origins.empty() && !right.q_origins.empty())
        types.push_back("QQ");
    if ((!left.q_origins.empty() && !right.rho_origins.empty()) ||
        (!left.rho_origins.empty() && !right.q_origins.empty()))
        types.push_back("QR");
    if (!left.rho_origins.empty() && !right.rho_origins.empty())
        types.push_back("RR");
    require(!types.empty(), "untyped digit edge");
    return types;
}

struct Collision {
    u32 hit_index = 0;
    u32 p = 0;
    u32 s = 0;
    u32 t = 0;
    u32 h = 0;
    std::vector<std::string> types;
    bool automatic_reflection = false;
};

struct ContinuantSupport {
    u32 s = 0;
    u32 t = 0;
    u32 h = 0;
    std::vector<std::string> types;
    std::vector<u32> primes;
    std::vector<u32> reflection_primes;
    std::vector<u32> nonreflection_primes;
    std::string prime_product;
    std::string nonreflection_prime_product;
};

struct GapSource {
    std::string kind;
    std::vector<u32> hit_indices;
};

struct GapInfo {
    u32 h = 0;
    std::vector<GapSource> sources;
};

struct ExtraReturn {
    u32 hit_index = 0;
    u32 p = 0;
    u32 h = 0;
    u32 x = 0;
    u32 y = 0;
    u32 other_zero = 0;
    std::string direction;
    std::vector<GapSource> sources;
};

struct RowResult {
    u64 n = 0;
    std::vector<Hit> hits;
    std::vector<Digit> digits;
    std::vector<std::vector<int>> quotient_zero_matrix;
    std::vector<std::vector<int>> residue_zero_matrix;
    std::vector<Collision> collisions;
    std::vector<ContinuantSupport> continuant_supports;
    std::vector<GapInfo> quotient_gaps;
    std::vector<ExtraReturn> extra_returns;
};

struct ScaleChecks {
    BuildChecks recurrence;
    bool baseline_available = false;
    bool baseline_passed = false;
    u64 clique_groups_checked = 0;
    u64 clique_memberships_checked = 0;
    u64 canonical_hit_checks = 0;
    u64 continuant_minor_checks = 0;
    u64 extra_return_zero_equivalence_checks = 0;
    u64 reflection_gap_tests_excluded = 0;
    ScatterOutcome scatter;
};

struct ScaleResult {
    u32 X = 0;
    u64 limit = 0;
    u32 prime_count = 0;
    u32 active_prime_count = 0;
    u32 total_zeros = 0;
    u64 pair_records = 0;
    u64 pair_groups = 0;
    std::map<u32, u64> highload_histogram;
    u32 max_K = 0;
    u64 canonical_triples = 0;
    u64 ordered_S3 = 0;

    u64 quotient_zero_ones = 0;
    u64 residue_zero_ones = 0;
    u64 residue_offdiagonal_raw = 0;
    u64 residue_offdiagonal_distinct_digit = 0;
    u64 collision_events = 0;
    u64 collision_type_labels = 0;
    u64 automatic_reflection_collisions = 0;
    u64 nonreflection_collisions = 0;
    u64 nonempty_continuant_supports = 0;
    u64 shared_continuant_supports = 0;
    u64 shared_nonreflection_continuant_supports = 0;
    u64 quotient_gap_count = 0;
    u64 quotient_generated_extra_returns = 0;

    ScaleChecks checks;
    std::vector<RowResult> rows;
};

struct SupportBuilder {
    std::vector<u32> primes;
    std::vector<u32> reflection_primes;
    std::vector<u32> nonreflection_primes;
};

static const Digit &find_digit(const std::vector<Digit> &digits, u32 value) {
    auto it = std::lower_bound(
        digits.begin(), digits.end(), value,
        [](const Digit &digit, u32 target) { return digit.value < target; });
    require(it != digits.end() && it->value == value,
            "digit lookup failed");
    return *it;
}

static RowResult enrich_row(u32 X, const RowBase &base,
                            const std::vector<PrimeData> &prime_data,
                            ScaleResult &scale) {
    RowResult row;
    row.n = base.n;
    const u32 K = static_cast<u32>(base.ids.size());
    row.hits.reserve(K);
    std::map<u32, Digit> digit_map;

    for (u32 hit_index = 0; hit_index < K; ++hit_index) {
        const u32 id = base.ids[hit_index];
        const PrimeData &item = prime_data[id];
        const u32 q = static_cast<u32>(base.n / item.p);
        const u32 r = static_cast<u32>(base.n % item.p);
        const u32 rho = std::min(r, item.p - 1 - r);
        require(q < X, "quotient digit is not below X");
        require(rho < X, "canonical residue digit is not below X");
        require(static_cast<u64>(q) * item.p + r == base.n,
                "quotient/residue reconstruction failed");
        require(item.B[r] == 0 && item.B[rho] == 0,
                "hit or canonical reflection is not an Apéry zero");
        ++scale.checks.canonical_hit_checks;
        row.hits.push_back({id, item.p, q, r, rho,
                            static_cast<u32>(item.zeros.size())});

        Digit &qdigit = digit_map[q];
        qdigit.value = q;
        qdigit.q_origins.push_back(hit_index);
        Digit &rdigit = digit_map[rho];
        rdigit.value = rho;
        rdigit.rho_origins.push_back(hit_index);
    }
    for (auto &entry : digit_map) row.digits.push_back(std::move(entry.second));

    row.quotient_zero_matrix.assign(K, std::vector<int>(K, 0));
    row.residue_zero_matrix.assign(K, std::vector<int>(K, 0));
    for (u32 i = 0; i < K; ++i) {
        const PrimeData &item = prime_data[row.hits[i].id];
        for (u32 j = 0; j < K; ++j) {
            const bool qzero = item.B[row.hits[j].q] == 0;
            const bool rzero = item.B[row.hits[j].rho] == 0;
            row.quotient_zero_matrix[i][j] = qzero ? 1 : 0;
            row.residue_zero_matrix[i][j] = rzero ? 1 : 0;
            if (qzero) ++scale.quotient_zero_ones;
            if (rzero) {
                ++scale.residue_zero_ones;
                if (i != j) {
                    ++scale.residue_offdiagonal_raw;
                    if (row.hits[i].rho != row.hits[j].rho)
                        ++scale.residue_offdiagonal_distinct_digit;
                }
            }
        }
        require(row.residue_zero_matrix[i][i] == 1,
                "residue matrix diagonal is not one");
    }

    std::map<std::pair<u32, u32>, SupportBuilder> support_map;
    for (u32 i = 0; i < K; ++i) {
        const PrimeData &item = prime_data[row.hits[i].id];
        for (std::size_t a = 0; a < row.digits.size(); ++a) {
            for (std::size_t b = a + 1; b < row.digits.size(); ++b) {
                const u32 s = row.digits[a].value;
                const u32 t = row.digits[b].value;
                const u32 h = t - s;
                require(t < item.p, "digit state wraps characteristic");
                require(!(item.B[s] == 0 && item.D[s] == 0),
                        "zero projective state at left endpoint");
                require(!(item.B[t] == 0 && item.D[t] == 0),
                        "zero projective state at right endpoint");
                const u32 minor = state_minor(item, s, t);
                const u32 continuant = continuant_mod(s, h, item.p);
                const u64 expected = mul_mod(
                    sixth_mod(item.factorial[s], item.p),
                    continuant, item.p);
                require(minor == expected,
                        "state-minor/continuant identity mismatch");
                ++scale.checks.continuant_minor_checks;
                if (minor != 0) continue;

                const bool reflection =
                    s == row.hits[i].rho &&
                    t == item.p - 1 - row.hits[i].rho;
                const std::vector<std::string> types =
                    edge_types(row.digits[a], row.digits[b]);
                row.collisions.push_back(
                    {i, item.p, s, t, h, types, reflection});
                ++scale.collision_events;
                scale.collision_type_labels += types.size();
                if (reflection)
                    ++scale.automatic_reflection_collisions;
                else
                    ++scale.nonreflection_collisions;

                SupportBuilder &builder = support_map[{s, t}];
                builder.primes.push_back(item.p);
                if (reflection)
                    builder.reflection_primes.push_back(item.p);
                else
                    builder.nonreflection_primes.push_back(item.p);
            }
        }
    }

    for (const auto &entry : support_map) {
        const u32 s = entry.first.first;
        const u32 t = entry.first.second;
        const Digit &left = find_digit(row.digits, s);
        const Digit &right = find_digit(row.digits, t);
        u128 product = 1;
        for (u32 p : entry.second.primes) product *= p;
        u128 nonreflection_product = 1;
        for (u32 p : entry.second.nonreflection_primes)
            nonreflection_product *= p;
        row.continuant_supports.push_back(
            {s, t, t - s, edge_types(left, right),
             entry.second.primes, entry.second.reflection_primes,
             entry.second.nonreflection_primes, u128_decimal(product),
             u128_decimal(nonreflection_product)});
        ++scale.nonempty_continuant_supports;
        if (entry.second.primes.size() >= 2)
            ++scale.shared_continuant_supports;
        if (entry.second.nonreflection_primes.size() >= 2)
            ++scale.shared_nonreflection_continuant_supports;
    }

    std::map<u32, std::vector<GapSource>> gap_map;
    for (u32 j = 0; j < K; ++j) {
        if (row.hits[j].q >= 2)
            gap_map[row.hits[j].q].push_back({"q", {j}});
    }
    for (u32 j = 0; j < K; ++j) {
        for (u32 k = j + 1; k < K; ++k) {
            const u32 left = row.hits[j].q;
            const u32 right = row.hits[k].q;
            const u32 difference = left > right ? left - right : right - left;
            if (difference >= 2)
                gap_map[difference].push_back(
                    {"abs_quotient_difference", {j, k}});
        }
    }
    for (auto &entry : gap_map) {
        auto &sources = entry.second;
        std::sort(sources.begin(), sources.end(),
                  [](const GapSource &left, const GapSource &right) {
                      const int left_rank = left.kind == "q" ? 0 : 1;
                      const int right_rank = right.kind == "q" ? 0 : 1;
                      if (left_rank != right_rank)
                          return left_rank < right_rank;
                      return left.hit_indices < right.hit_indices;
                  });
        row.quotient_gaps.push_back({entry.first, sources});
    }
    scale.quotient_gap_count += row.quotient_gaps.size();

    for (u32 i = 0; i < K; ++i) {
        const PrimeData &item = prime_data[row.hits[i].id];
        const u32 rho = row.hits[i].rho;
        const u32 reflection_gap = item.p - 1 - 2 * rho;
        for (const GapInfo &gap : row.quotient_gaps) {
            if (gap.h == reflection_gap) {
                ++scale.checks.reflection_gap_tests_excluded;
                continue;
            }
            std::vector<std::pair<u32, std::string>> starts;
            if (static_cast<u64>(rho) + gap.h < item.p)
                starts.push_back({rho, "up"});
            if (rho >= gap.h)
                starts.push_back({rho - gap.h, "down"});
            for (const auto &start : starts) {
                const u32 x = start.first;
                const u32 y = x + gap.h;
                require(x < y && y < item.p,
                        "invalid quotient-generated return interval");
                const u32 minor = state_minor(item, x, y);
                const u32 continuant = continuant_mod(x, gap.h, item.p);
                const u64 expected = mul_mod(
                    sixth_mod(item.factorial[x], item.p),
                    continuant, item.p);
                require(minor == expected,
                        "extra-return minor/continuant mismatch");
                ++scale.checks.continuant_minor_checks;
                const bool direct_zero_pair =
                    item.B[x] == 0 && item.B[y] == 0;
                require((minor == 0) == direct_zero_pair,
                        "extra-return state/zero-set equivalence mismatch");
                ++scale.checks.extra_return_zero_equivalence_checks;
                if (minor != 0) continue;
                const u32 other = x == rho ? y : x;
                require(other != item.p - 1 - rho,
                        "reflection return survived exclusion");
                row.extra_returns.push_back(
                    {i, item.p, gap.h, x, y, other,
                     start.second, gap.sources});
                ++scale.quotient_generated_extra_returns;
            }
        }
    }

    std::sort(row.collisions.begin(), row.collisions.end(),
              [](const Collision &left, const Collision &right) {
                  return std::tie(left.p, left.s, left.t) <
                         std::tie(right.p, right.s, right.t);
              });
    std::sort(row.extra_returns.begin(), row.extra_returns.end(),
              [](const ExtraReturn &left, const ExtraReturn &right) {
                  return std::tie(left.p, left.h, left.x) <
                         std::tie(right.p, right.h, right.x);
              });
    return row;
}

static ScaleResult analyze_scale(u32 X) {
    require(X >= 8, "X must be at least 8");
    require(static_cast<u64>(X) * X <=
                std::numeric_limits<u64>::max(),
            "X^2 overflow");
    ScaleResult result;
    result.X = X;
    result.limit = static_cast<u64>(X) * X;

    const std::vector<u32> primes = primes_between(X);
    result.prime_count = static_cast<u32>(primes.size());
    std::vector<PrimeData> prime_data;
    prime_data.reserve(primes.size());
    for (u32 p : primes)
        prime_data.push_back(build_prime_data(p,
                                              result.checks.recurrence));
    for (const PrimeData &item : prime_data) {
        if (!item.zeros.empty()) ++result.active_prime_count;
        result.total_zeros += static_cast<u32>(item.zeros.size());
    }

    const std::vector<PairRec> pair_records =
        enumerate_pair_records(prime_data, result.limit);
    result.pair_records = pair_records.size();
    const DecodeResult decoded = decode_rows(pair_records, prime_data);
    result.pair_groups = decoded.groups_checked;
    result.checks.clique_groups_checked = decoded.groups_checked;
    result.checks.clique_memberships_checked =
        decoded.row_memberships_checked;

    for (const RowBase &row : decoded.rows) {
        const u32 K = static_cast<u32>(row.ids.size());
        ++result.highload_histogram[K];
        result.max_K = std::max(result.max_K, K);
        result.canonical_triples += choose3(K);
    }
    result.ordered_S3 = 6 * result.canonical_triples;

    Baseline baseline;
    result.checks.baseline_available = baseline_for(X, baseline);
    if (result.checks.baseline_available) {
        require(result.prime_count == baseline.prime_count,
                "prime-count baseline mismatch");
        require(result.active_prime_count == baseline.active_count,
                "active-prime baseline mismatch");
        require(result.total_zeros == baseline.zero_count,
                "zero-count baseline mismatch");
        require(result.pair_records == baseline.pair_records,
                "pair-record baseline mismatch");
        require(result.highload_histogram == baseline.highload_histogram,
                "high-load histogram baseline mismatch");
        result.checks.baseline_passed = true;
    }

    if (X == 256) {
        result.checks.scatter = direct_scatter_check(
            X, prime_data, decoded.rows, result.pair_records);
    }

    result.rows.reserve(decoded.rows.size());
    for (const RowBase &row : decoded.rows)
        result.rows.push_back(enrich_row(X, row, prime_data, result));
    require(result.rows.size() == decoded.rows.size(),
            "enrichment lost rows");
    return result;
}

static std::string json_escape(const std::string &text) {
    std::ostringstream out;
    for (unsigned char character : text) {
        switch (character) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\b': out << "\\b"; break;
        case '\f': out << "\\f"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (character < 0x20) {
                out << "\\u" << std::hex << std::setw(4)
                    << std::setfill('0') << static_cast<int>(character)
                    << std::dec << std::setfill(' ');
            } else {
                out << character;
            }
        }
    }
    return out.str();
}

static void indent(std::ostream &out, int spaces) {
    for (int i = 0; i < spaces; ++i) out.put(' ');
}

template <class T>
static void print_number_array(std::ostream &out,
                               const std::vector<T> &values) {
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << ',';
        out << values[i];
    }
    out << ']';
}

static void print_string_array(std::ostream &out,
                               const std::vector<std::string> &values) {
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << ',';
        out << '"' << json_escape(values[i]) << '"';
    }
    out << ']';
}

static void print_matrix(std::ostream &out,
                         const std::vector<std::vector<int>> &matrix) {
    out << '[';
    for (std::size_t i = 0; i < matrix.size(); ++i) {
        if (i) out << ',';
        print_number_array(out, matrix[i]);
    }
    out << ']';
}

static void print_gap_sources(std::ostream &out,
                              const std::vector<GapSource> &sources) {
    out << '[';
    for (std::size_t i = 0; i < sources.size(); ++i) {
        if (i) out << ',';
        out << "{\"kind\":\"" << json_escape(sources[i].kind)
            << "\",\"hit_indices\":";
        print_number_array(out, sources[i].hit_indices);
        out << '}';
    }
    out << ']';
}

static void print_row_json(std::ostream &out, const RowResult &row,
                           int base_indent) {
    indent(out, base_indent); out << "{\n";
    indent(out, base_indent + 2); out << "\"n\":" << row.n << ",\n";
    indent(out, base_indent + 2); out << "\"K\":"
                                     << row.hits.size() << ",\n";
    indent(out, base_indent + 2); out << "\"hits\":[\n";
    for (std::size_t i = 0; i < row.hits.size(); ++i) {
        const Hit &hit = row.hits[i];
        indent(out, base_indent + 4);
        out << "{\"hit_index\":" << i
            << ",\"p\":" << hit.p
            << ",\"q\":" << hit.q
            << ",\"r\":" << hit.r
            << ",\"rho\":" << hit.rho
            << ",\"Zp\":" << hit.zero_count << '}';
        if (i + 1 != row.hits.size()) out << ',';
        out << '\n';
    }
    indent(out, base_indent + 2); out << "],\n";

    indent(out, base_indent + 2); out << "\"digit_pool\":[\n";
    for (std::size_t i = 0; i < row.digits.size(); ++i) {
        const Digit &digit = row.digits[i];
        indent(out, base_indent + 4);
        out << "{\"value\":" << digit.value << ",\"roles\":[";
        bool comma = false;
        if (!digit.q_origins.empty()) {
            out << "\"Q\"";
            comma = true;
        }
        if (!digit.rho_origins.empty()) {
            if (comma) out << ',';
            out << "\"R\"";
        }
        out << "],\"q_origins\":";
        print_number_array(out, digit.q_origins);
        out << ",\"rho_origins\":";
        print_number_array(out, digit.rho_origins);
        out << '}';
        if (i + 1 != row.digits.size()) out << ',';
        out << '\n';
    }
    indent(out, base_indent + 2); out << "],\n";

    indent(out, base_indent + 2);
    out << "\"quotient_zero_matrix\":{";
    out << "\"row_primes\":[";
    for (std::size_t i = 0; i < row.hits.size(); ++i) {
        if (i) out << ',';
        out << row.hits[i].p;
    }
    out << "],\"column_hit_indices\":[";
    for (std::size_t i = 0; i < row.hits.size(); ++i) {
        if (i) out << ',';
        out << i;
    }
    out << "],\"values\":";
    print_matrix(out, row.quotient_zero_matrix);
    out << "},\n";

    indent(out, base_indent + 2);
    out << "\"residue_cross_zero_matrix\":{";
    out << "\"row_primes\":[";
    for (std::size_t i = 0; i < row.hits.size(); ++i) {
        if (i) out << ',';
        out << row.hits[i].p;
    }
    out << "],\"column_hit_indices\":[";
    for (std::size_t i = 0; i < row.hits.size(); ++i) {
        if (i) out << ',';
        out << i;
    }
    out << "],\"values\":";
    print_matrix(out, row.residue_zero_matrix);
    out << "},\n";

    indent(out, base_indent + 2);
    out << "\"typed_projective_state_collisions\":[\n";
    for (std::size_t i = 0; i < row.collisions.size(); ++i) {
        const Collision &collision = row.collisions[i];
        indent(out, base_indent + 4);
        out << "{\"hit_index\":" << collision.hit_index
            << ",\"p\":" << collision.p
            << ",\"s\":" << collision.s
            << ",\"t\":" << collision.t
            << ",\"h\":" << collision.h
            << ",\"types\":";
        print_string_array(out, collision.types);
        out << ",\"minor\":0,\"automatic_reflection\":"
            << (collision.automatic_reflection ? "true" : "false")
            << '}';
        if (i + 1 != row.collisions.size()) out << ',';
        out << '\n';
    }
    indent(out, base_indent + 2); out << "],\n";

    indent(out, base_indent + 2); out << "\"continuant_supports\":[\n";
    for (std::size_t i = 0; i < row.continuant_supports.size(); ++i) {
        const ContinuantSupport &support = row.continuant_supports[i];
        indent(out, base_indent + 4);
        out << "{\"s\":" << support.s
            << ",\"t\":" << support.t
            << ",\"h\":" << support.h
            << ",\"types\":";
        print_string_array(out, support.types);
        out << ",\"primes\":";
        print_number_array(out, support.primes);
        out << ",\"reflection_primes\":";
        print_number_array(out, support.reflection_primes);
        out << ",\"nonreflection_primes\":";
        print_number_array(out, support.nonreflection_primes);
        out << ",\"prime_product\":\"" << support.prime_product
            << "\",\"nonreflection_prime_product\":\""
            << support.nonreflection_prime_product << "\"}";
        if (i + 1 != row.continuant_supports.size()) out << ',';
        out << '\n';
    }
    indent(out, base_indent + 2); out << "],\n";

    indent(out, base_indent + 2); out << "\"quotient_gap_catalog\":[\n";
    for (std::size_t i = 0; i < row.quotient_gaps.size(); ++i) {
        const GapInfo &gap = row.quotient_gaps[i];
        indent(out, base_indent + 4);
        out << "{\"h\":" << gap.h << ",\"sources\":";
        print_gap_sources(out, gap.sources);
        out << '}';
        if (i + 1 != row.quotient_gaps.size()) out << ',';
        out << '\n';
    }
    indent(out, base_indent + 2); out << "],\n";

    indent(out, base_indent + 2);
    out << "\"nonreflection_quotient_generated_extra_returns\":[\n";
    for (std::size_t i = 0; i < row.extra_returns.size(); ++i) {
        const ExtraReturn &event = row.extra_returns[i];
        indent(out, base_indent + 4);
        out << "{\"hit_index\":" << event.hit_index
            << ",\"p\":" << event.p
            << ",\"h\":" << event.h
            << ",\"x\":" << event.x
            << ",\"y\":" << event.y
            << ",\"other_zero\":" << event.other_zero
            << ",\"direction\":\"" << event.direction
            << "\",\"sources\":";
        print_gap_sources(out, event.sources);
        out << '}';
        if (i + 1 != row.extra_returns.size()) out << ',';
        out << '\n';
    }
    indent(out, base_indent + 2); out << "],\n";

    indent(out, base_indent + 2); out << "\"row_counts\":{";
    out << "\"digit_count\":" << row.digits.size()
        << ",\"collision_events\":" << row.collisions.size()
        << ",\"continuant_supports\":"
        << row.continuant_supports.size()
        << ",\"quotient_gaps\":" << row.quotient_gaps.size()
        << ",\"extra_returns\":" << row.extra_returns.size()
        << "}\n";
    indent(out, base_indent); out << '}';
}

static void print_scale_json(std::ostream &out, const ScaleResult &scale,
                             int base_indent) {
    indent(out, base_indent); out << "{\n";
    indent(out, base_indent + 2); out << "\"X\":" << scale.X << ",\n";
    indent(out, base_indent + 2); out << "\"n_limit\":"
                                     << scale.limit << ",\n";
    indent(out, base_indent + 2); out << "\"prime_count\":"
                                     << scale.prime_count << ",\n";
    indent(out, base_indent + 2); out << "\"active_prime_count\":"
                                     << scale.active_prime_count << ",\n";
    indent(out, base_indent + 2); out << "\"sum_Zp\":"
                                     << scale.total_zeros << ",\n";
    indent(out, base_indent + 2); out << "\"pair_records\":"
                                     << scale.pair_records << ",\n";
    indent(out, base_indent + 2); out << "\"pair_groups\":"
                                     << scale.pair_groups << ",\n";
    indent(out, base_indent + 2); out << "\"rows_K_ge_3\":"
                                     << scale.rows.size() << ",\n";
    indent(out, base_indent + 2); out << "\"K_histogram\":{";
    std::size_t hist_index = 0;
    for (const auto &entry : scale.highload_histogram) {
        if (hist_index++) out << ',';
        out << '"' << entry.first << "\":" << entry.second;
    }
    out << "},\n";
    indent(out, base_indent + 2); out << "\"max_K\":"
                                     << scale.max_K << ",\n";
    indent(out, base_indent + 2); out << "\"canonical_triples\":"
                                     << scale.canonical_triples << ",\n";
    indent(out, base_indent + 2); out << "\"ordered_S3\":"
                                     << scale.ordered_S3 << ",\n";

    indent(out, base_indent + 2); out << "\"aggregate_invariants\":{";
    out << "\"quotient_zero_ones\":" << scale.quotient_zero_ones
        << ",\"residue_zero_ones\":" << scale.residue_zero_ones
        << ",\"residue_offdiagonal_raw\":"
        << scale.residue_offdiagonal_raw
        << ",\"residue_offdiagonal_distinct_digit\":"
        << scale.residue_offdiagonal_distinct_digit
        << ",\"collision_events\":" << scale.collision_events
        << ",\"collision_type_labels\":"
        << scale.collision_type_labels
        << ",\"automatic_reflection_collisions\":"
        << scale.automatic_reflection_collisions
        << ",\"nonreflection_collisions\":"
        << scale.nonreflection_collisions
        << ",\"nonempty_continuant_supports\":"
        << scale.nonempty_continuant_supports
        << ",\"shared_continuant_supports\":"
        << scale.shared_continuant_supports
        << ",\"shared_nonreflection_continuant_supports\":"
        << scale.shared_nonreflection_continuant_supports
        << ",\"quotient_gap_count\":" << scale.quotient_gap_count
        << ",\"quotient_generated_extra_returns\":"
        << scale.quotient_generated_extra_returns << "},\n";

    const ScaleChecks &checks = scale.checks;
    indent(out, base_indent + 2); out << "\"self_checks\":{";
    out << "\"status\":\"PASS\""
        << ",\"baseline_available\":"
        << (checks.baseline_available ? "true" : "false")
        << ",\"baseline_passed\":"
        << (checks.baseline_passed ? "true" : "false")
        << ",\"recurrence_primes_checked\":"
        << checks.recurrence.primes_checked
        << ",\"cleared_divided_positions\":"
        << checks.recurrence.cleared_divided_positions
        << ",\"direct_binomial_positions\":"
        << checks.recurrence.direct_binomial_positions
        << ",\"reflection_positions\":"
        << checks.recurrence.reflection_positions
        << ",\"reflection_zero_positions\":"
        << checks.recurrence.reflection_zero_positions
        << ",\"wronskian_positions\":"
        << checks.recurrence.wronskian_positions
        << ",\"nonadjacent_zero_positions\":"
        << checks.recurrence.nonadjacent_zero_positions
        << ",\"clique_groups_checked\":"
        << checks.clique_groups_checked
        << ",\"clique_memberships_checked\":"
        << checks.clique_memberships_checked
        << ",\"canonical_hit_checks\":"
        << checks.canonical_hit_checks
        << ",\"continuant_minor_checks\":"
        << checks.continuant_minor_checks
        << ",\"extra_return_zero_equivalence_checks\":"
        << checks.extra_return_zero_equivalence_checks
        << ",\"reflection_gap_tests_excluded\":"
        << checks.reflection_gap_tests_excluded
        << ",\"direct_scatter_performed\":"
        << (checks.scatter.performed ? "true" : "false")
        << ",\"direct_scatter_passed\":"
        << (checks.scatter.passed ? "true" : "false")
        << ",\"direct_scatter_hit_events\":"
        << checks.scatter.hit_events
        << ",\"direct_scatter_rows\":" << checks.scatter.rows
        << ",\"direct_scatter_pair_records\":"
        << checks.scatter.pair_records
        << ",\"direct_scatter_canonical_triples\":"
        << checks.scatter.canonical_triples << "},\n";

    indent(out, base_indent + 2); out << "\"rows\":[\n";
    for (std::size_t i = 0; i < scale.rows.size(); ++i) {
        print_row_json(out, scale.rows[i], base_indent + 4);
        if (i + 1 != scale.rows.size()) out << ',';
        out << '\n';
    }
    indent(out, base_indent + 2); out << "]\n";
    indent(out, base_indent); out << '}';
}

static void print_json(std::ostream &out,
                       const std::vector<ScaleResult> &scales,
                       const std::string &source_sha256,
                       bool synthetic_clique_passed) {
    out << "{\n";
    out << "  \"schema\":\"q7306-highload-invariant-scan-v1\",\n";
    out << "  \"project\":\"Ramanujan_Challenge\",\n";
    out << "  \"sequence\":\"b_n=sum_j binom(n,j)^2 binom(n+j,j)^2\",\n";
    out << "  \"prime_interval\":\"X < p <= 2X\",\n";
    out << "  \"row_interval\":\"0 <= n < X^2\",\n";
    out << "  \"canonical_residue\":\"rho=min(r,p-1-r)\",\n";
    out << "  \"source_path\":\"problems/3.2/research/q7306_highload_scan.cpp\",\n";
    out << "  \"source_sha256\":\"" << json_escape(source_sha256)
        << "\",\n";
    out << "  \"compile_flags\":\"-std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic\",\n";
    out << "  \"global_self_checks\":{";
    out << "\"status\":\"PASS\",\"synthetic_clique_decoder\":"
        << (synthetic_clique_passed ? "true" : "false") << "},\n";
    out << "  \"scales\":[\n";
    for (std::size_t i = 0; i < scales.size(); ++i) {
        print_scale_json(out, scales[i], 4);
        if (i + 1 != scales.size()) out << ',';
        out << '\n';
    }
    out << "  ]\n";
    out << "}\n";
}

static u64 histogram_value(const ScaleResult &scale, u32 K) {
    auto it = scale.highload_histogram.find(K);
    return it == scale.highload_histogram.end() ? 0 : it->second;
}

static void print_summary(std::ostream &out,
                          const std::vector<ScaleResult> &scales,
                          const std::string &source_sha256,
                          bool synthetic_clique_passed) {
    out << "Q7306 exact high-load invariant scanner\n";
    out << "source_sha256 " << source_sha256 << '\n';
    out << "synthetic_clique_decoder "
        << (synthetic_clique_passed ? "PASS" : "FAIL") << '\n';
    out << "X primes active sumZ pairrec rows K3 K4 maxK Q1 Roff* "
           "coll(nonrefl/refl) sharedNR extra checks\n";
    for (const ScaleResult &scale : scales) {
        out << scale.X << ' '
            << scale.prime_count << ' '
            << scale.active_prime_count << ' '
            << scale.total_zeros << ' '
            << scale.pair_records << ' '
            << scale.rows.size() << ' '
            << histogram_value(scale, 3) << ' '
            << histogram_value(scale, 4) << ' '
            << scale.max_K << ' '
            << scale.quotient_zero_ones << ' '
            << scale.residue_offdiagonal_distinct_digit << ' '
            << scale.nonreflection_collisions << '/'
            << scale.automatic_reflection_collisions << ' '
            << scale.shared_nonreflection_continuant_supports << ' '
            << scale.quotient_generated_extra_returns << ' '
            << "PASS" << '\n';
    }
    out << "Roff* excludes off-diagonal columns that repeat the row prime's "
           "own canonical digit.\n";
    out << "All reported extra returns exclude the automatic reflection gap.\n";
}

static u32 parse_X(const std::string &text) {
    std::size_t used = 0;
    const unsigned long value = std::stoul(text, &used, 10);
    require(used == text.size(), "invalid X argument: " + text);
    require(value <= std::numeric_limits<u32>::max(),
            "X exceeds uint32 range");
    return static_cast<u32>(value);
}

int main(int argc, char **argv) {
    try {
        std::string source_sha256 = "unspecified";
        std::vector<u32> scales_requested;
        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            if (argument == "--source-sha256") {
                require(i + 1 < argc, "--source-sha256 requires a value");
                source_sha256 = argv[++i];
            } else if (argument == "--help") {
                std::cout << "usage: q7306_scan [--source-sha256 HEX] "
                             "[X ...]\n";
                return 0;
            } else {
                scales_requested.push_back(parse_X(argument));
            }
        }
        if (scales_requested.empty())
            scales_requested = {256, 512, 1024, 2048};
        std::sort(scales_requested.begin(), scales_requested.end());
        require(std::adjacent_find(scales_requested.begin(),
                                   scales_requested.end()) ==
                    scales_requested.end(),
                "duplicate X argument");

        const bool synthetic_clique_passed = synthetic_clique_self_check();
        std::vector<ScaleResult> results;
        results.reserve(scales_requested.size());
        for (u32 X : scales_requested) results.push_back(analyze_scale(X));

        print_summary(std::cerr, results, source_sha256,
                      synthetic_clique_passed);
        print_json(std::cout, results, source_sha256,
                   synthetic_clique_passed);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Q7306 scanner failure: " << error.what() << '\n';
        return 1;
    }
}
