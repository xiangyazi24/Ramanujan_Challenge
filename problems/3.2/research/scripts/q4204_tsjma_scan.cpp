// Q4204 rigorous computational TSJMA/SRT campaign.
// Exact arithmetic backend; C++20, no third-party libraries.
//
// The program exhaustively generates Apéry zero sets for every prime up to
// qmax, enumerates the exact q6 raw fibres by an inverse lower-row map,
// applies the raw minus-first assignment, and computes translated Q4190
// cell statistics at L=Y^(5/6), H=Y^(1/2).  It also audits folded-prefix
// depths, local stripped radical support, p-adic prefix valuations, Smith
// factors, and auxiliary-modulus fingerprints.

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using std::int64_t;
using std::uint64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

namespace {

struct Options {
    int qmax = 300000;
    int min_lower = 7;
    int clean_lower = 17;
    int threads = 0;
    int top_cells = 24;
    std::string betas = "0.10,0.15,0.30,0.50";
    std::string out_dir = "q4204_results";
};

struct Timer {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double seconds() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }
};

uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t m) {
    return static_cast<uint64_t>((static_cast<u128>(a) * b) % m);
}

uint64_t add_mod(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t r = a + b;
    if (r >= m || r < a) r %= m;
    return r;
}

uint64_t sub_mod(uint64_t a, uint64_t b, uint64_t m) {
    return a >= b ? a - b : m - (b - a) % m;
}

uint64_t pow_mod(uint64_t a, uint64_t e, uint64_t m) {
    uint64_t r = 1 % m;
    while (e) {
        if (e & 1) r = mul_mod(r, a, m);
        a = mul_mod(a, a, m);
        e >>= 1;
    }
    return r;
}

uint64_t inv_mod(uint64_t a, uint64_t m) {
    // Extended Euclid with 128-bit intermediates.  All callers ensure gcd=1.
    i128 t = 0, new_t = 1;
    i128 r = static_cast<i128>(m), new_r = static_cast<i128>(a % m);
    while (new_r != 0) {
        i128 q = r / new_r;
        i128 tmp_t = t - q * new_t;
        t = new_t;
        new_t = tmp_t;
        i128 tmp_r = r - q * new_r;
        r = new_r;
        new_r = tmp_r;
    }
    if (r != 1) throw std::runtime_error("nonunit modular inverse");
    t %= static_cast<i128>(m);
    if (t < 0) t += m;
    return static_cast<uint64_t>(t);
}

std::vector<int> sieve_primes(int n, std::vector<uint8_t>& is_prime) {
    is_prime.assign(n + 1, 1);
    if (n >= 0) is_prime[0] = 0;
    if (n >= 1) is_prime[1] = 0;
    for (int p = 2; static_cast<int64_t>(p) * p <= n; ++p) {
        if (!is_prime[p]) continue;
        for (int64_t k = static_cast<int64_t>(p) * p; k <= n; k += p) {
            is_prime[static_cast<size_t>(k)] = 0;
        }
    }
    std::vector<int> primes;
    primes.reserve(static_cast<size_t>(n / std::max(1.0, std::log(std::max(3, n)))) + 100);
    for (int p = 2; p <= n; ++p) if (is_prime[p]) primes.push_back(p);
    return primes;
}

bool has_zero(const std::vector<std::vector<int>>& zeros, int p, int r) {
    if (p < 0 || p >= static_cast<int>(zeros.size()) || r < 0 || r >= p) return false;
    const auto& z = zeros[p];
    return std::binary_search(z.begin(), z.end(), r);
}

std::vector<int> apery_zero_set_half(
    int p,
    std::vector<uint32_t>& inv,
    uint64_t& recurrence_steps
) {
    if (p < 3) return {};
    const int half = (p - 1) / 2;
    inv.assign(static_cast<size_t>(half + 1), 0);
    if (half >= 1) inv[1] = 1;
    for (int a = 2; a <= half; ++a) {
        inv[a] = static_cast<uint32_t>(
            (p - (static_cast<uint64_t>(p / a) * inv[p % a]) % p) % p
        );
    }

    std::vector<int> left;
    uint64_t previous = 1 % p;
    uint64_t current = 5 % p;
    if (half >= 1 && current == 0) left.push_back(1);

    for (int n = 1; n < half; ++n) {
        uint64_t nn = static_cast<uint64_t>(n);
        uint64_t n2 = (nn * nn) % p;
        uint64_t n3 = (n2 * nn) % p;
        uint64_t middle = ((2 * nn + 1) % p) *
            ((17 * n2 + 17 * nn + 5) % p) % p;
        uint64_t rhs = (middle * current) % p;
        uint64_t sub = (n3 * previous) % p;
        rhs = rhs >= sub ? rhs - sub : rhs + p - sub;
        uint64_t u = inv[n + 1];
        uint64_t next = rhs * u % p * u % p * u % p;
        previous = current;
        current = next;
        if (current == 0) left.push_back(n + 1);
        ++recurrence_steps;
    }

    std::vector<int> result = left;
    for (int r : left) {
        int reflected = p - 1 - r;
        if (reflected != r) result.push_back(reflected);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    if (std::binary_search(result.begin(), result.end(), 0)) {
        throw std::runtime_error("b_0 unexpectedly zero");
    }
    for (size_t i = 1; i < result.size(); ++i) {
        if (result[i] == result[i - 1] + 1) {
            std::ostringstream os;
            os << "consecutive Apéry zeros at p=" << p;
            throw std::runtime_error(os.str());
        }
    }
    return result;
}

struct RowKey {
    int q = 0;
    int t = 0;
    int p = 0;
    bool operator<(const RowKey& other) const {
        return std::tie(q, t, p) < std::tie(other.q, other.t, other.p);
    }
    bool operator==(const RowKey& other) const {
        return q == other.q && t == other.t && p == other.p;
    }
};

struct Row {
    int q = 0;
    int t = 0;
    int p = 0;
    int r = 0;
    int a = 0;
    char sign = '+';
    bool selected = false;
    bool clean = false;
    int zq = 0;

    RowKey key() const { return {q, t, p}; }
    int n() const { return 6 * q + t; }
    int m() const { return sign == '+' ? 12 * q + t : q - 1 - t; }
    int slope2() const { return sign == '+' ? 13 : 1; }
    int delta_state() const { return sign == '+' ? 6 * q + 7 * t : 7 * t + 6; }
    int j6() const { return std::min(r, p - 1 - r); }
    int js() const { return std::min(a, p - 1 - a); }
};

struct Enumeration {
    std::vector<Row> plus_candidate_raw;
    std::vector<Row> minus_candidate_raw;
    std::vector<Row> plus_actual_raw;
    std::vector<Row> minus_actual_raw;
    std::vector<Row> plus_candidate_marked;
    std::vector<Row> minus_candidate_marked;
    std::vector<Row> plus_actual_marked;
    std::vector<Row> minus_actual_marked;
    uint64_t plus_integral = 0;
    uint64_t minus_integral = 0;
    uint64_t plus_qprime_legal = 0;
    uint64_t minus_qprime_legal = 0;
};

bool coefficient_clean(const std::vector<std::vector<int>>& zeros, int p, char sign) {
    if (p < 17) return false;
    if (has_zero(zeros, p, 6)) return false;
    if (sign == '+' && has_zero(zeros, p, 13)) return false;
    if (sign == '-' && has_zero(zeros, p, 1)) return false;
    return true;
}

Enumeration enumerate_rows(
    int qmax,
    int min_lower,
    const std::vector<int>& primes,
    const std::vector<uint8_t>& is_prime,
    const std::vector<std::vector<int>>& zeros
) {
    Enumeration out;
    std::set<RowKey> plus_candidate_keys, minus_candidate_keys;
    std::set<RowKey> plus_actual_keys, minus_actual_keys;

    for (int p : primes) {
        if (p < min_lower || p > qmax) continue;
        const auto& zp = zeros[p];
        if (zp.empty()) continue;
        for (int r : zp) {
            for (int a : zp) {
                // Plus inverse: q=(7p+a-r)/6, t=2r-a-p.
                int64_t numerator_plus = static_cast<int64_t>(7) * p + a - r;
                if (numerator_plus % 6 == 0) {
                    ++out.plus_integral;
                    int q = static_cast<int>(numerator_plus / 6);
                    int t = 2 * r - a - p;
                    if (q > p && q <= qmax && q >= 0 && is_prime[q] && 0 <= t && t < q) {
                        ++out.plus_qprime_legal;
                        Row row{q, t, p, r, a, '+', false,
                                coefficient_clean(zeros, p, '+'),
                                static_cast<int>(zeros[q].size())};
                        if (row.n() - 6 * p != r || row.m() - 13 * p != a) {
                            throw std::runtime_error("plus inverse mismatch");
                        }
                        if (plus_candidate_keys.insert(row.key()).second) {
                            out.plus_candidate_raw.push_back(row);
                        }
                        if (has_zero(zeros, q, t)) {
                            row.selected = true;
                            if (plus_actual_keys.insert(row.key()).second) {
                                out.plus_actual_raw.push_back(row);
                            }
                        }
                    }
                }

                // Minus inverse: q=p+(r+a+1)/7, t=(r-6a-6)/7.
                int64_t numerator_minus = static_cast<int64_t>(r) + a + 1;
                if (numerator_minus % 7 == 0) {
                    ++out.minus_integral;
                    int q = p + static_cast<int>(numerator_minus / 7);
                    int64_t tnum = static_cast<int64_t>(r) - 6LL * a - 6;
                    if (tnum % 7 != 0) throw std::runtime_error("minus integrality mismatch");
                    int t = static_cast<int>(tnum / 7);
                    if (q > p && q <= qmax && is_prime[q] && 0 <= t && t < q) {
                        ++out.minus_qprime_legal;
                        Row row{q, t, p, r, a, '-', false,
                                coefficient_clean(zeros, p, '-'),
                                static_cast<int>(zeros[q].size())};
                        if (row.n() - 6 * p != r || row.m() - p != a) {
                            throw std::runtime_error("minus inverse mismatch");
                        }
                        if (minus_candidate_keys.insert(row.key()).second) {
                            out.minus_candidate_raw.push_back(row);
                        }
                        if (has_zero(zeros, q, t)) {
                            row.selected = true;
                            if (minus_actual_keys.insert(row.key()).second) {
                                out.minus_actual_raw.push_back(row);
                            }
                        }
                    }
                }
            }
        }
    }

    // Authoritative raw minus-first assignment, separately for candidate and actual ledgers.
    for (const Row& row : out.minus_candidate_raw) out.minus_candidate_marked.push_back(row);
    for (const Row& row : out.plus_candidate_raw) {
        if (!minus_candidate_keys.count(row.key())) out.plus_candidate_marked.push_back(row);
    }
    for (const Row& row : out.minus_actual_raw) out.minus_actual_marked.push_back(row);
    for (const Row& row : out.plus_actual_raw) {
        if (!minus_actual_keys.count(row.key())) out.plus_actual_marked.push_back(row);
    }

    auto order_rows = [](std::vector<Row>& rows) {
        std::sort(rows.begin(), rows.end(), [](const Row& x, const Row& y) {
            return std::tie(x.q, x.t, x.p, x.sign) < std::tie(y.q, y.t, y.p, y.sign);
        });
    };
    order_rows(out.plus_candidate_raw);
    order_rows(out.minus_candidate_raw);
    order_rows(out.plus_actual_raw);
    order_rows(out.minus_actual_raw);
    order_rows(out.plus_candidate_marked);
    order_rows(out.minus_candidate_marked);
    order_rows(out.plus_actual_marked);
    order_rows(out.minus_actual_marked);
    return out;
}

std::vector<double> parse_betas(const std::string& text) {
    std::vector<double> result;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) continue;
        double b = std::stod(item);
        if (!(b > 0.0 && b < 1.0)) throw std::runtime_error("beta must lie in (0,1)");
        result.push_back(b);
    }
    if (result.empty()) throw std::runtime_error("empty beta list");
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

int64_t floor_div(int64_t x, int64_t d) {
    if (d <= 0) throw std::runtime_error("nonpositive divisor");
    if (x >= 0) return x / d;
    return - ((-x + d - 1) / d);
}

struct StateKey {
    int q = 0;
    int t = 0;
    char sign = '+';
    bool operator<(const StateKey& other) const {
        return std::tie(q, t, sign) < std::tie(other.q, other.t, other.sign);
    }
};

using Groups = std::map<StateKey, std::vector<Row>>;

Groups group_rows(const std::vector<Row>& plus, const std::vector<Row>& minus, bool clean_only) {
    Groups groups;
    auto add = [&](const std::vector<Row>& rows) {
        for (const Row& row : rows) {
            if (clean_only && !row.clean) continue;
            groups[{row.q, row.t, row.sign}].push_back(row);
        }
    };
    add(plus);
    add(minus);
    for (auto& [key, rows] : groups) {
        std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.p < b.p; });
        rows.erase(std::unique(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
            return a.p == b.p;
        }), rows.end());
    }
    return groups;
}

const std::vector<int> AUX_PRIMES = {5,7,11,13,17,19,23,29,31,37,41,43,47};

std::vector<int> apery_values_mod_small(int ell) {
    std::vector<int> b(ell, 0), inv(ell, 0);
    b[0] = 1 % ell;
    if (ell > 1) b[1] = 5 % ell;
    inv[1] = 1;
    for (int a = 2; a < ell; ++a) {
        inv[a] = static_cast<int>((ell - static_cast<int64_t>(ell / a) * inv[ell % a] % ell) % ell);
    }
    for (int n = 1; n < ell - 1; ++n) {
        int64_t n2 = static_cast<int64_t>(n) * n % ell;
        int64_t n3 = n2 * n % ell;
        int64_t middle = (2LL * n + 1) % ell * ((17 * n2 + 17LL * n + 5) % ell) % ell;
        int64_t rhs = (middle * b[n] - n3 * b[n - 1]) % ell;
        if (rhs < 0) rhs += ell;
        int64_t u = inv[n + 1];
        b[n + 1] = static_cast<int>(rhs * u % ell * u % ell * u % ell);
    }
    return b;
}

int apery_lucas_mod(int64_t n, int ell, const std::vector<int>& digits) {
    if (n < 0) return 0;
    int64_t x = n;
    int result = 1 % ell;
    while (x) {
        int d = static_cast<int>(x % ell);
        result = static_cast<int>(static_cast<int64_t>(result) * digits[d] % ell);
        x /= ell;
    }
    return result;
}

int aux_divisor_count(int delta) {
    int count = 0;
    for (int ell : AUX_PRIMES) if (delta % ell == 0) ++count;
    return count;
}

int pair_aux_zero_count(
    const Row& base,
    const Row& second,
    const std::map<int, std::vector<int>>& aux_values
) {
    int h = second.p - base.p;
    int c = base.slope2();
    int count = 0;
    for (int ell : AUX_PRIMES) {
        const auto& vals = aux_values.at(ell);
        int br = apery_lucas_mod(base.r, ell, vals);
        int ba = apery_lucas_mod(base.a, ell, vals);
        int bshift_a = apery_lucas_mod(static_cast<int64_t>(base.a) - static_cast<int64_t>(c) * h, ell, vals);
        int bshift_r = apery_lucas_mod(static_cast<int64_t>(base.r) - 6LL * h, ell, vals);
        int value = static_cast<int>((static_cast<int64_t>(br) * bshift_a +
                                      static_cast<int64_t>(ba) * bshift_r) % ell);
        if (value == 0) ++count;
    }
    return count;
}

struct CellFeature {
    bool actual = false;
    char sign = '+';
    int sector = 0; // 1 singleton, 2 non-singleton
    int T = 0;
    double beta = 0;
    int Y = 0, L = 0, H = 0;
    int q = 0, t = 0, tau = 0;
    int64_t cell = 0;
    int K = 0, J6 = 0, Js = 0;
    int delta = 0;
    int delta_aux = 0;
    int max_gap_delta_gcd = 1;
    double mean_pair_aux_zeros = 0;
};

struct WorstCell {
    bool actual = false;
    char sign = '+';
    int sector = 0;
    int T = 0;
    double beta = 0;
    int Y = 0, L = 0, H = 0;
    int q = 0, t = 0, tau = 0;
    int64_t cell = 0;
    int K = 0, J6 = 0, Js = 0;
    std::vector<Row> rows;

    int64_t lo() const { return static_cast<int64_t>(tau) + cell * L; }
    int64_t hi() const { return lo() + L; }
};

struct Agg {
    uint64_t states = 0;
    uint64_t M = 0;
    long double triangular_num = 0; // numerator with denominator Y
    uint64_t energy_num = 0;       // denominator L
    uint64_t rich_num = 0;         // K>H, denominator L
    uint64_t sat_num = 0;          // K>H and J6,Js>H
    uint64_t rich2_num = 0;        // K>2H
    uint64_t sat2_num = 0;
    uint64_t rich4_num = 0;
    uint64_t sat4_num = 0;
    uint64_t prefix_thin_rich_num = 0;
    uint64_t multi_cells = 0;
    int maxK = 0;
    int max_prefix_ratio_num = 0;
    int max_prefix_ratio_den = 1;
};

struct MetricRow {
    std::string dataset;
    char sign = '+';
    int sector = 0;
    int T = 0, Y = 0, L = 0, H = 0;
    double beta = 0;
    Agg agg;
};

void maybe_store_worst(std::vector<WorstCell>& worst, WorstCell cell, int cap) {
    if (cell.K < 2) return;
    worst.push_back(std::move(cell));
    std::sort(worst.begin(), worst.end(), [](const WorstCell& a, const WorstCell& b) {
        if (a.K != b.K) return a.K > b.K;
        int amin = std::min(a.J6, a.Js), bmin = std::min(b.J6, b.Js);
        if (amin != bmin) return amin > bmin;
        return std::tie(a.actual, a.T, a.q, a.t, a.tau, a.cell) <
               std::tie(b.actual, b.T, b.q, b.t, b.tau, b.cell);
    });
    if (static_cast<int>(worst.size()) > cap) worst.resize(cap);
}

void add_group_metrics(
    const std::vector<Row>& rows,
    int T,
    double beta,
    int Y,
    int L,
    int H,
    bool actual,
    int sector,
    std::array<Agg, 3>& aggs,
    std::vector<WorstCell>& worst,
    std::vector<CellFeature>& features,
    int top_cap,
    const std::map<int, std::vector<int>>& aux_values
) {
    if (rows.empty()) return;
    auto apply_all = [&](auto fn) {
        fn(aggs[sector - 1]);
        fn(aggs[2]);
    };
    apply_all([&](Agg& a) {
        ++a.states;
        a.M += rows.size();
    });

    for (size_t i = 0; i < rows.size(); ++i) {
        for (size_t j = i + 1; j < rows.size(); ++j) {
            int d = rows[j].p - rows[i].p;
            if (d < Y) {
                apply_all([&](Agg& a) { a.triangular_num += (Y - d); });
            }
        }
    }

    if (rows.size() == 1) {
        apply_all([&](Agg& a) {
            a.energy_num += L;
            a.maxK = std::max(a.maxK, 1);
            int den = 1 + std::min(rows[0].j6(), rows[0].js());
            if (a.max_prefix_ratio_num * den < a.max_prefix_ratio_den) {
                a.max_prefix_ratio_num = 1;
                a.max_prefix_ratio_den = den;
            }
        });
        return;
    }

    for (int tau = 0; tau < L; ++tau) {
        std::map<int64_t, std::vector<Row>> cells;
        for (const Row& row : rows) {
            int64_t j = floor_div(static_cast<int64_t>(row.p) - tau, L);
            cells[j].push_back(row);
        }
        for (auto& [cell_index, cr] : cells) {
            int K = static_cast<int>(cr.size());
            int J6 = 0, Js = 0;
            for (const Row& row : cr) {
                J6 = std::max(J6, row.j6());
                Js = std::max(Js, row.js());
            }
            bool sat = J6 > H && Js > H;
            apply_all([&](Agg& a) {
                a.energy_num += static_cast<uint64_t>(K) * K;
                a.maxK = std::max(a.maxK, K);
                int den = 1 + std::min(J6, Js);
                if (static_cast<int64_t>(a.max_prefix_ratio_num) * den <
                    static_cast<int64_t>(K) * a.max_prefix_ratio_den) {
                    a.max_prefix_ratio_num = K;
                    a.max_prefix_ratio_den = den;
                }
                if (K >= 2) ++a.multi_cells;
                if (K > H) {
                    a.rich_num += K;
                    if (sat) a.sat_num += K;
                    if (!sat) a.prefix_thin_rich_num += K;
                }
                if (K > 2 * H) {
                    a.rich2_num += K;
                    if (sat) a.sat2_num += K;
                }
                if (K > 4 * H) {
                    a.rich4_num += K;
                    if (sat) a.sat4_num += K;
                }
            });

            if (K >= 2) {
                int delta = cr.front().delta_state();
                int maxg = 1;
                int auxsum = 0;
                int pairs = 0;
                for (size_t i = 0; i < cr.size(); ++i) {
                    for (size_t j = i + 1; j < cr.size(); ++j) {
                        int gap = cr[j].p - cr[i].p;
                        maxg = std::max(maxg, std::gcd(gap, std::abs(delta)));
                        auxsum += pair_aux_zero_count(cr[i], cr[j], aux_values);
                        ++pairs;
                    }
                }
                CellFeature f;
                f.actual = actual; f.sign = cr.front().sign; f.sector = sector;
                f.T = T; f.beta = beta; f.Y = Y; f.L = L; f.H = H;
                f.q = cr.front().q; f.t = cr.front().t; f.tau = tau;
                f.cell = cell_index; f.K = K; f.J6 = J6; f.Js = Js;
                f.delta = delta; f.delta_aux = aux_divisor_count(delta);
                f.max_gap_delta_gcd = maxg;
                f.mean_pair_aux_zeros = pairs ? static_cast<double>(auxsum) / pairs : 0.0;
                features.push_back(f);

                WorstCell w;
                w.actual = actual; w.sign = cr.front().sign; w.sector = sector;
                w.T = T; w.beta = beta; w.Y = Y; w.L = L; w.H = H;
                w.q = cr.front().q; w.t = cr.front().t; w.tau = tau;
                w.cell = cell_index; w.K = K; w.J6 = J6; w.Js = Js;
                w.rows = cr;
                maybe_store_worst(worst, std::move(w), top_cap);
            }
        }
    }
}

std::vector<MetricRow> compute_metrics(
    const std::string& dataset,
    const Groups& groups,
    int qmax,
    const std::vector<double>& betas,
    bool actual,
    std::vector<WorstCell>& worst,
    std::vector<CellFeature>& features,
    int top_cap,
    const std::map<int, std::vector<int>>& aux_values
) {
    std::vector<MetricRow> result;
    std::vector<int> scales;
    for (int T = 16; 2LL * T <= qmax; T *= 2) scales.push_back(T);

    for (double beta : betas) {
        for (int T : scales) {
            int Y = std::max(1, static_cast<int>(std::floor(std::pow(static_cast<long double>(T), beta))));
            int L = std::max(1, static_cast<int>(std::floor(std::pow(static_cast<long double>(Y), 5.0L / 6.0L))));
            int H = std::max(1, static_cast<int>(std::floor(std::sqrt(static_cast<long double>(Y)))));
            for (char sign : {'+', '-'}) {
                std::array<Agg, 3> aggs{};
                for (const auto& [key, rows] : groups) {
                    if (key.sign != sign || key.q < T || key.q >= 2 * T) continue;
                    int sector = rows.front().zq == 1 ? 1 : 2;
                    add_group_metrics(rows, T, beta, Y, L, H, actual, sector,
                                      aggs, worst, features, top_cap, aux_values);
                }
                for (int sector = 1; sector <= 3; ++sector) {
                    result.push_back({dataset, sign, sector, T, Y, L, H, beta, aggs[sector - 1]});
                }
            }
        }
    }
    return result;
}

uint64_t prefix_mod(uint64_t M, int R, uint64_t mod) {
    if (R < 0) return 0;
    uint64_t A = 1 % mod; // binom(M,k)
    uint64_t B = 1 % mod; // binom(M+k,k)
    uint64_t sum = 0;
    for (int k = 0; k <= R; ++k) {
        uint64_t term = mul_mod(mul_mod(A, A, mod), mul_mod(B, B, mod), mod);
        sum = add_mod(sum, term, mod);
        if (k == R) break;
        uint64_t den = static_cast<uint64_t>(k + 1);
        uint64_t inv = inv_mod(den, mod);
        uint64_t x1 = (M - static_cast<uint64_t>(k)) % mod;
        uint64_t x2 = (M + static_cast<uint64_t>(k) + 1) % mod;
        A = mul_mod(mul_mod(A, x1, mod), inv, mod);
        B = mul_mod(mul_mod(B, x2, mod), inv, mod);
    }
    return sum;
}

bool prefix_zero_mod_prime(
    int M,
    int J,
    int ell,
    const std::vector<std::vector<int>>& zeros
) {
    int r = M % ell;
    int folded = std::min(r, ell - 1 - r);
    if (J >= folded) return has_zero(zeros, ell, folded);
    if (J >= ell) throw std::runtime_error("prefix cutoff reached modulus");
    return prefix_mod(static_cast<uint64_t>(M), J, ell) == 0;
}

bool row_zero_clean(
    int M,
    int ell,
    const std::vector<std::vector<int>>& zeros
) {
    int c = M / ell;
    int r = M % ell;
    if (c >= ell) throw std::runtime_error("Gessel-Lucas quotient too large");
    if (has_zero(zeros, ell, c)) return false; // fixed coefficient exception stripped
    return has_zero(zeros, ell, r);
}

int prefix_valuation_cap3(int M, int J, int p) {
    u128 p3wide = static_cast<u128>(p) * p * p;
    if (p3wide > std::numeric_limits<uint64_t>::max()) return -1;
    uint64_t p2 = static_cast<uint64_t>(p) * p;
    uint64_t p3 = static_cast<uint64_t>(p3wide);
    uint64_t value = prefix_mod(static_cast<uint64_t>(M), J, p3);
    if (value % p != 0) return 0;
    if (value % p2 != 0) return 1;
    if (value % p3 != 0) return 2;
    return 3;
}

std::vector<std::pair<int,int>> factor_small(int n, const std::vector<int>& primes) {
    int x = std::abs(n);
    std::vector<std::pair<int,int>> f;
    for (int p : primes) {
        if (static_cast<int64_t>(p) * p > x) break;
        if (x % p) continue;
        int e = 0;
        while (x % p == 0) { x /= p; ++e; }
        f.push_back({p,e});
    }
    if (x > 1) f.push_back({x,1});
    return f;
}

std::string factor_string(int n, const std::vector<int>& primes) {
    auto f = factor_small(n, primes);
    if (f.empty()) return "1";
    std::ostringstream os;
    for (size_t i = 0; i < f.size(); ++i) {
        if (i) os << "*";
        os << f[i].first;
        if (f[i].second > 1) os << "^" << f[i].second;
    }
    return os.str();
}

struct RadicalDiag {
    WorstCell cell;
    std::vector<int> support;
    std::vector<int> phantoms;
    double log_radical = 0;
    double log_actual = 0;
    bool all_actual_supported = true;
    std::vector<std::tuple<int,int,int,int>> valuations; // p,v6,vs,min
};

RadicalDiag diagnose_radical(
    const WorstCell& cell,
    const std::vector<int>& primes,
    const std::vector<std::vector<int>>& zeros
) {
    RadicalDiag d;
    d.cell = cell;
    int n = 6 * cell.q + cell.t;
    int m = cell.sign == '+' ? 12 * cell.q + cell.t : cell.q - 1 - cell.t;
    int64_t lo = std::max<int64_t>(17, cell.lo());
    int64_t hi = cell.hi();
    auto it = std::lower_bound(primes.begin(), primes.end(), static_cast<int>(lo));
    for (; it != primes.end() && *it < hi; ++it) {
        int ell = *it;
        if (ell == cell.q) continue;
        if (!row_zero_clean(n, ell, zeros) || !row_zero_clean(m, ell, zeros)) continue;
        if (!prefix_zero_mod_prime(n, cell.J6, ell, zeros)) continue;
        if (!prefix_zero_mod_prime(m, cell.Js, ell, zeros)) continue;
        d.support.push_back(ell);
        d.log_radical += std::log(static_cast<double>(ell));
    }
    std::set<int> actual;
    for (const Row& row : cell.rows) {
        actual.insert(row.p);
        d.log_actual += std::log(static_cast<double>(row.p));
        int v6 = prefix_valuation_cap3(n, cell.J6, row.p);
        int vs = prefix_valuation_cap3(m, cell.Js, row.p);
        d.valuations.push_back({row.p, v6, vs, std::min(v6,vs)});
    }
    std::set<int> support(d.support.begin(), d.support.end());
    for (int p : actual) if (!support.count(p)) d.all_actual_supported = false;
    for (int ell : d.support) if (!actual.count(ell)) d.phantoms.push_back(ell);
    return d;
}

double pearson(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 3) return std::numeric_limits<double>::quiet_NaN();
    double mx = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    double my = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
    long double num = 0, dx = 0, dy = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        long double a = x[i] - mx, b = y[i] - my;
        num += a * b; dx += a * a; dy += b * b;
    }
    if (dx == 0 || dy == 0) return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(num / std::sqrt(dx * dy));
}

std::string sector_name(int sector) {
    if (sector == 1) return "singleton";
    if (sector == 2) return "non-singleton";
    return "all";
}

std::string format_ratio(uint64_t num, uint64_t den) {
    if (den == 0) return "NA";
    std::ostringstream os;
    os << std::setprecision(8) << static_cast<long double>(num) / den;
    return os.str();
}

std::string format_ld(long double x) {
    std::ostringstream os;
    os << std::setprecision(10) << x;
    return os.str();
}

std::optional<std::pair<double,double>> fit_power(
    const std::vector<MetricRow>& metrics,
    const std::string& dataset,
    char sign,
    int sector,
    double beta,
    bool saturated
) {
    std::vector<double> x,y;
    for (const auto& row : metrics) {
        if (row.dataset != dataset || row.sign != sign || row.sector != sector ||
            std::abs(row.beta - beta) > 1e-12 || row.agg.M == 0) continue;
        uint64_t num = saturated ? row.agg.sat_num : row.agg.rich_num;
        if (num == 0) continue;
        double ratio = static_cast<double>(num) / (row.L * static_cast<double>(row.agg.M));
        if (ratio <= 0) continue;
        x.push_back(std::log(static_cast<double>(row.T)));
        y.push_back(std::log(ratio));
    }
    if (x.size() < 2) return std::nullopt;
    double mx = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    double my = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
    double den = 0, num = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        den += (x[i]-mx)*(x[i]-mx);
        num += (x[i]-mx)*(y[i]-my);
    }
    if (den == 0) return std::nullopt;
    double slope = num / den;
    return std::make_pair(-slope, static_cast<double>(x.size()));
}

void write_csv_rows(const std::string& path, const std::vector<Row>& plus, const std::vector<Row>& minus) {
    std::ofstream out(path);
    out << "sign,q,t,zq,p,r,a,selected,clean,delta_state,j6,js\n";
    auto emit = [&](const std::vector<Row>& rows) {
        for (const Row& r : rows) {
            out << r.sign << ',' << r.q << ',' << r.t << ',' << r.zq << ',' << r.p << ','
                << r.r << ',' << r.a << ',' << r.selected << ',' << r.clean << ','
                << r.delta_state() << ',' << r.j6() << ',' << r.js() << '\n';
        }
    };
    emit(plus); emit(minus);
}

void write_metric_csv(const std::string& path, const std::vector<MetricRow>& metrics) {
    std::ofstream out(path);
    out << "dataset,sign,sector,T,beta,Y,L,H,states,M,triangular,energy_over_M,rich_mass,rich_over_M,sat_mass,sat_over_M,rich2,sat2,rich4,sat4,prefix_thin_rich,maxK,max_prefix_ratio,multi_cells\n";
    for (const auto& r : metrics) {
        long double tri = r.Y ? r.agg.triangular_num / r.Y : 0;
        long double energy = r.L ? static_cast<long double>(r.agg.energy_num) / r.L : 0;
        long double rich = r.L ? static_cast<long double>(r.agg.rich_num) / r.L : 0;
        long double sat = r.L ? static_cast<long double>(r.agg.sat_num) / r.L : 0;
        out << r.dataset << ',' << r.sign << ',' << sector_name(r.sector) << ','
            << r.T << ',' << r.beta << ',' << r.Y << ',' << r.L << ',' << r.H << ','
            << r.agg.states << ',' << r.agg.M << ',' << static_cast<double>(tri) << ','
            << (r.agg.M ? static_cast<double>(energy / r.agg.M) : 0) << ','
            << static_cast<double>(rich) << ',' << (r.agg.M ? static_cast<double>(rich/r.agg.M) : 0) << ','
            << static_cast<double>(sat) << ',' << (r.agg.M ? static_cast<double>(sat/r.agg.M) : 0) << ','
            << static_cast<double>(r.agg.rich2_num) / r.L << ','
            << static_cast<double>(r.agg.sat2_num) / r.L << ','
            << static_cast<double>(r.agg.rich4_num) / r.L << ','
            << static_cast<double>(r.agg.sat4_num) / r.L << ','
            << static_cast<double>(r.agg.prefix_thin_rich_num) / r.L << ','
            << r.agg.maxK << ','
            << static_cast<double>(r.agg.max_prefix_ratio_num) / r.agg.max_prefix_ratio_den << ','
            << r.agg.multi_cells << '\n';
    }
}

Options parse_options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + name);
            return argv[++i];
        };
        if (a == "--qmax") o.qmax = std::stoi(need("--qmax"));
        else if (a == "--min-lower") o.min_lower = std::stoi(need("--min-lower"));
        else if (a == "--clean-lower") o.clean_lower = std::stoi(need("--clean-lower"));
        else if (a == "--threads") o.threads = std::stoi(need("--threads"));
        else if (a == "--top-cells") o.top_cells = std::stoi(need("--top-cells"));
        else if (a == "--betas") o.betas = need("--betas");
        else if (a == "--out-dir") o.out_dir = need("--out-dir");
        else throw std::runtime_error("unknown option: " + a);
    }
    if (o.qmax < 1000) throw std::runtime_error("qmax must be at least 1000");
    return o;
}

} // namespace

int main(int argc, char** argv) {
    try {
        Options opt = parse_options(argc, argv);
        std::vector<double> betas = parse_betas(opt.betas);
        int threads = opt.threads > 0 ? opt.threads : static_cast<int>(std::thread::hardware_concurrency());
        threads = std::max(1, threads);
        Timer total_timer;

        std::cerr << "Q4204_START qmax=" << opt.qmax << " threads=" << threads
                  << " betas=" << opt.betas << "\n";

        std::vector<uint8_t> is_prime;
        std::vector<int> primes = sieve_primes(opt.qmax, is_prime);
        std::vector<std::vector<int>> zeros(static_cast<size_t>(opt.qmax + 1));

        std::vector<int> work;
        for (int p : primes) if (p >= 5) work.push_back(p);
        std::sort(work.rbegin(), work.rend());
        std::atomic<size_t> next{0};
        std::atomic<uint64_t> completed{0};
        std::atomic<uint64_t> global_steps{0};
        std::mutex print_mutex;
        Timer zero_timer;
        std::vector<std::thread> pool;
        for (int tid = 0; tid < threads; ++tid) {
            pool.emplace_back([&, tid]() {
                std::vector<uint32_t> inv;
                uint64_t local_steps = 0;
                while (true) {
                    size_t index = next.fetch_add(1);
                    if (index >= work.size()) break;
                    int p = work[index];
                    zeros[p] = apery_zero_set_half(p, inv, local_steps);
                    uint64_t done = completed.fetch_add(1) + 1;
                    if (done % 2500 == 0 || done == work.size()) {
                        std::lock_guard<std::mutex> lock(print_mutex);
                        std::cerr << "ZERO_PROGRESS " << done << '/' << work.size()
                                  << " elapsed=" << std::fixed << std::setprecision(2)
                                  << zero_timer.seconds() << "s\n";
                    }
                }
                global_steps.fetch_add(local_steps);
            });
        }
        for (auto& th : pool) th.join();
        double zero_seconds = zero_timer.seconds();
        std::cerr << "ZERO_DONE seconds=" << zero_seconds
                  << " recurrence_steps=" << global_steps.load() << "\n";

        // Exact zero-set census and banked checkpoints.
        uint64_t selected_states = 0, singleton_primes = 0, nonempty_primes = 0;
        int max_z = 0;
        std::vector<int> max_z_primes, singleton_list;
        for (int q : primes) {
            if (q < 7) continue;
            int z = static_cast<int>(zeros[q].size());
            selected_states += z;
            if (z) ++nonempty_primes;
            if (z == 1) { ++singleton_primes; singleton_list.push_back(q); }
            if (z > max_z) { max_z = z; max_z_primes = {q}; }
            else if (z == max_z && z) max_z_primes.push_back(q);
        }
        uint64_t states_1000 = 0, lower_columns_1000 = 0, states_17_5000 = 0;
        for (int q : primes) {
            if (q >= 7 && q <= 1000) states_1000 += zeros[q].size();
            if (q >= 7 && q <= 1000 && zeros[q].size() >= 2) ++lower_columns_1000;
            if (q >= 17 && q <= 5000) states_17_5000 += zeros[q].size();
        }

        Timer enum_timer;
        Enumeration enumeration = enumerate_rows(opt.qmax, opt.min_lower, primes, is_prime, zeros);
        double enum_seconds = enum_timer.seconds();
        std::cerr << "ENUM_DONE seconds=" << enum_seconds
                  << " actual_plus_raw=" << enumeration.plus_actual_raw.size()
                  << " actual_minus_raw=" << enumeration.minus_actual_raw.size() << "\n";

        auto clean_count = [](const std::vector<Row>& rows) {
            return std::count_if(rows.begin(), rows.end(), [](const Row& r){ return r.clean; });
        };
        auto count_q_le = [](const std::vector<Row>& rows, int bound, bool clean) {
            return std::count_if(rows.begin(), rows.end(), [&](const Row& r){
                return r.q <= bound && (!clean || r.clean);
            });
        };
        uint64_t prior_actual_5000 = count_q_le(enumeration.plus_actual_marked, 5000, true) +
                                     count_q_le(enumeration.minus_actual_marked, 5000, true);
        uint64_t prior_actual_100000 = count_q_le(enumeration.plus_actual_marked, 100000, true) +
                                       count_q_le(enumeration.minus_actual_marked, 100000, true);

        Groups actual_groups = group_rows(enumeration.plus_actual_marked,
                                          enumeration.minus_actual_marked, true);
        Groups ambient_groups = group_rows(enumeration.plus_candidate_marked,
                                           enumeration.minus_candidate_marked, true);

        std::map<int, std::vector<int>> aux_values;
        for (int ell : AUX_PRIMES) aux_values[ell] = apery_values_mod_small(ell);

        std::vector<WorstCell> worst_actual, worst_ambient;
        std::vector<CellFeature> features_actual, features_ambient;
        Timer metric_timer;
        std::vector<MetricRow> metrics_actual = compute_metrics(
            "actual", actual_groups, opt.qmax, betas, true,
            worst_actual, features_actual, opt.top_cells, aux_values);
        std::vector<MetricRow> metrics_ambient = compute_metrics(
            "unselected-overcarrier", ambient_groups, opt.qmax, betas, false,
            worst_ambient, features_ambient, opt.top_cells, aux_values);
        std::vector<MetricRow> metrics = metrics_actual;
        metrics.insert(metrics.end(), metrics_ambient.begin(), metrics_ambient.end());
        double metric_seconds = metric_timer.seconds();
        std::cerr << "METRIC_DONE seconds=" << metric_seconds
                  << " actual_multi_cells=" << features_actual.size()
                  << " ambient_multi_cells=" << features_ambient.size() << "\n";

        // Diagnose exact stripped local radicals for the worst cells.  Keep the
        // workload bounded but exact on every reported cell.
        std::vector<RadicalDiag> radical_actual, radical_ambient;
        Timer radical_timer;
        for (size_t i = 0; i < worst_actual.size() && i < static_cast<size_t>(opt.top_cells); ++i) {
            radical_actual.push_back(diagnose_radical(worst_actual[i], primes, zeros));
        }
        for (size_t i = 0; i < worst_ambient.size() && i < static_cast<size_t>(opt.top_cells); ++i) {
            radical_ambient.push_back(diagnose_radical(worst_ambient[i], primes, zeros));
        }
        double radical_seconds = radical_timer.seconds();
        std::cerr << "RADICAL_DONE seconds=" << radical_seconds << "\n";

        // Correlation diagnostics on all multi-leaf translated cells.
        auto correlations = [](const std::vector<CellFeature>& f) {
            std::map<std::string,double> result;
            std::vector<double> y, depth, delta_aux, gap_gcd, pair_aux;
            for (const auto& x : f) {
                y.push_back(std::log(static_cast<double>(x.K)));
                depth.push_back(std::log1p(static_cast<double>(std::min(x.J6,x.Js)) / std::max(1,x.H)));
                delta_aux.push_back(x.delta_aux);
                gap_gcd.push_back(std::log(static_cast<double>(x.max_gap_delta_gcd)));
                pair_aux.push_back(x.mean_pair_aux_zeros);
            }
            result["min_folded_depth_over_H"] = pearson(y, depth);
            result["small_aux_divisors_of_Delta"] = pearson(y, delta_aux);
            result["max_gcd(gap,Delta)"] = pearson(y, gap_gcd);
            result["primitive_bilinear_aux_zero_count"] = pearson(y, pair_aux);
            return result;
        };
        auto corr_actual = correlations(features_actual);
        auto corr_ambient = correlations(features_ambient);

        std::string mkdir_cmd = "mkdir -p '" + opt.out_dir + "'";
        if (std::system(mkdir_cmd.c_str()) != 0) throw std::runtime_error("mkdir failed");
        write_csv_rows(opt.out_dir + "/actual_rows.csv",
                       enumeration.plus_actual_marked, enumeration.minus_actual_marked);
        write_csv_rows(opt.out_dir + "/unselected_overcarrier_rows.csv",
                       enumeration.plus_candidate_marked, enumeration.minus_candidate_marked);
        write_metric_csv(opt.out_dir + "/translated_metrics.csv", metrics);

        std::ofstream report(opt.out_dir + "/report.md");
        report << "# Q4204 exact TSJMA/SRT computation\n\n";
        report << "## Execution ledger\n\n```text\n";
        report << "qmax                       = " << opt.qmax << "\n";
        report << "threads                    = " << threads << "\n";
        report << "primes <= qmax             = " << primes.size() << "\n";
        report << "half-orbit recurrence steps= " << global_steps.load() << "\n";
        report << "zero-set seconds           = " << zero_seconds << "\n";
        report << "enumeration seconds        = " << enum_seconds << "\n";
        report << "translated-metric seconds  = " << metric_seconds << "\n";
        report << "radical-diagnostic seconds = " << radical_seconds << "\n";
        report << "total seconds              = " << total_timer.seconds() << "\n";
        report << "betas                      = " << opt.betas << "\n```\n\n";

        report << "## Exact Apéry zero-set census\n\n```text\n";
        report << "selected states (q,t)      = " << selected_states << "\n";
        report << "nonempty zero-set primes   = " << nonempty_primes << "\n";
        report << "singleton zero-set primes  = " << singleton_primes << "\n";
        report << "maximum z_q                = " << max_z << "\n";
        report << "states q<=1000             = " << states_1000 << " (Q4186 checkpoint 163)\n";
        report << "lower z_p>=2, p<=1000      = " << lower_columns_1000 << " (Q4186 checkpoint 68)\n";
        if (opt.qmax >= 5000)
            report << "states 17<=q<=5000         = " << states_17_5000 << " (Q4181 checkpoint 605)\n";
        report << "singleton primes           =";
        for (int q : singleton_list) report << ' ' << q;
        report << "\nmax-z primes               =";
        for (size_t i = 0; i < max_z_primes.size() && i < 20; ++i) report << ' ' << max_z_primes[i];
        if (max_z_primes.size() > 20) report << " ...";
        report << "\n```\n\n";

        report << "## Inverse q6 row enumeration\n\n";
        report << "Every lower zero pair `(p,r,a)` is mapped exactly by\n\n";
        report << "```text\nplus:  q=(7p+a-r)/6,       t=2r-a-p\n";
        report << "minus: q=p+(r+a+1)/7,      t=(r-6a-6)/7\n```\n\n";
        report << "| ledger | plus raw | minus raw | plus minus-first | minus minus-first | clean plus | clean minus |\n";
        report << "|:---|---:|---:|---:|---:|---:|---:|\n";
        report << "| q-prime lower-complete overcarrier | "
               << enumeration.plus_candidate_raw.size() << " | " << enumeration.minus_candidate_raw.size()
               << " | " << enumeration.plus_candidate_marked.size() << " | "
               << enumeration.minus_candidate_marked.size() << " | "
               << clean_count(enumeration.plus_candidate_marked) << " | "
               << clean_count(enumeration.minus_candidate_marked) << " |\n";
        report << "| actual selected q|b_t | "
               << enumeration.plus_actual_raw.size() << " | " << enumeration.minus_actual_raw.size()
               << " | " << enumeration.plus_actual_marked.size() << " | "
               << enumeration.minus_actual_marked.size() << " | "
               << clean_count(enumeration.plus_actual_marked) << " | "
               << clean_count(enumeration.minus_actual_marked) << " |\n\n";
        report << "```text\n";
        report << "plus integral lower pairs      = " << enumeration.plus_integral << "\n";
        report << "minus integral lower pairs     = " << enumeration.minus_integral << "\n";
        report << "plus q-prime/legal candidates  = " << enumeration.plus_qprime_legal << "\n";
        report << "minus q-prime/legal candidates = " << enumeration.minus_qprime_legal << "\n";
        report << "clean actual rows q<=5000      = " << prior_actual_5000 << " (Q4181 checkpoint 0)\n";
        if (opt.qmax >= 100000)
            report << "clean actual rows q<=100000    = " << prior_actual_100000 << " (prior certified range 0)\n";
        report << "```\n\n";

        auto emit_rows = [&](const std::vector<Row>& rows, const char* title) {
            report << "### " << title << "\n\n";
            if (rows.empty()) { report << "None.\n\n"; return; }
            report << "| sign | q | t | z_q | p | r | a | clean | Delta_sigma | j6 | js |\n";
            report << "|:---:|---:|---:|---:|---:|---:|---:|:---:|---:|---:|---:|\n";
            size_t cap = std::min<size_t>(rows.size(), 100);
            for (size_t i = 0; i < cap; ++i) {
                const Row& r = rows[i];
                report << '|' << r.sign << '|' << r.q << '|' << r.t << '|' << r.zq << '|'
                       << r.p << '|' << r.r << '|' << r.a << '|' << (r.clean?"yes":"no") << '|'
                       << r.delta_state() << '|' << r.j6() << '|' << r.js() << "|\n";
            }
            if (rows.size() > cap) report << "\nFirst " << cap << " of " << rows.size() << " shown.\n";
            report << '\n';
        };
        std::vector<Row> all_actual = enumeration.plus_actual_marked;
        all_actual.insert(all_actual.end(), enumeration.minus_actual_marked.begin(), enumeration.minus_actual_marked.end());
        std::sort(all_actual.begin(), all_actual.end(), [](const Row& a,const Row& b){
            return std::tie(a.q,a.t,a.p,a.sign)<std::tie(b.q,b.t,b.p,b.sign);
        });
        emit_rows(all_actual, "Actual selected rows");

        report << "## Translated Q4190 metrics\n\n";
        report << "For each complete dyadic shell `T<=q<2T`, `Y=floor(T^beta)`, "
               << "`L=floor(Y^(5/6))`, and `H=floor(Y^(1/2))`. "
               << "Rich and saturated-rich masses are exact averages over all `tau mod L`.\n\n";
        report << "| dataset | sign | sector | T | beta | Y | L | H | states | M | A_L/M | rich/M | sat/M | max K | max K/(1+min J) |\n";
        report << "|:---|:---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n";
        for (const auto& r : metrics) {
            if (r.agg.M == 0 && r.dataset == "unselected-overcarrier") continue;
            long double energy = r.L ? static_cast<long double>(r.agg.energy_num) / r.L : 0;
            long double rich = r.L ? static_cast<long double>(r.agg.rich_num) / r.L : 0;
            long double sat = r.L ? static_cast<long double>(r.agg.sat_num) / r.L : 0;
            report << '|' << r.dataset << '|' << r.sign << '|' << sector_name(r.sector) << '|'
                   << r.T << '|' << r.beta << '|' << r.Y << '|' << r.L << '|' << r.H << '|'
                   << r.agg.states << '|' << r.agg.M << '|';
            if (r.agg.M) report << static_cast<double>(energy/r.agg.M) << '|'
                                << static_cast<double>(rich/r.agg.M) << '|'
                                << static_cast<double>(sat/r.agg.M) << '|';
            else report << "NA|NA|NA|";
            report << r.agg.maxK << '|'
                   << static_cast<double>(r.agg.max_prefix_ratio_num)/r.agg.max_prefix_ratio_den << "|\n";
        }
        report << '\n';

        report << "## Power-law fits (positive observations only)\n\n";
        report << "The fitted quantity is `(translated rich first mass)/M = C*T^{-kappa}`. "
               << "Zero tails are censored and never replaced by pseudocounts.\n\n";
        report << "| dataset | sign | sector | beta | ordinary rich kappa_hat | saturated rich kappa_hat |\n";
        report << "|:---|:---:|:---|---:|:---|:---|\n";
        for (const std::string& ds : {std::string("actual"), std::string("unselected-overcarrier")}) {
            for (char sign : {'+','-'}) for (int sector : {1,2,3}) for (double beta : betas) {
                auto a = fit_power(metrics, ds, sign, sector, beta, false);
                auto b = fit_power(metrics, ds, sign, sector, beta, true);
                report << '|' << ds << '|' << sign << '|' << sector_name(sector) << '|' << beta << '|';
                if (a) report << a->first << " (n=" << static_cast<int>(a->second) << ")|";
                else report << "NA|";
                if (b) report << b->first << " (n=" << static_cast<int>(b->second) << ")|\n";
                else report << "NA|\n";
            }
        }
        report << '\n';

        auto emit_worst = [&](const std::vector<WorstCell>& cells, const char* title) {
            report << "## " << title << "\n\n";
            if (cells.empty()) { report << "No translated cell has occupancy at least two.\n\n"; return; }
            report << "| sign | sector | T | beta | Y | L | H | q | t | tau | interval | K | J6 | Js | Delta | max gcd(gap,Delta) |\n";
            report << "|:---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|---:|---:|---:|---:|\n";
            for (const auto& c : cells) {
                int maxg = 1;
                int delta = c.rows.front().delta_state();
                for (size_t i=0;i<c.rows.size();++i) for(size_t j=i+1;j<c.rows.size();++j)
                    maxg=std::max(maxg,std::gcd(c.rows[j].p-c.rows[i].p,std::abs(delta)));
                report << '|' << c.sign << '|' << sector_name(c.sector) << '|' << c.T << '|'
                       << c.beta << '|' << c.Y << '|' << c.L << '|' << c.H << '|'
                       << c.q << '|' << c.t << '|' << c.tau << "|[" << c.lo() << ',' << c.hi() << ")|"
                       << c.K << '|' << c.J6 << '|' << c.Js << '|' << delta << '|' << maxg << "|\n";
            }
            report << '\n';
        };
        emit_worst(worst_actual, "Worst actual selected cells");
        emit_worst(worst_ambient, "Worst unselected lower-complete overcarrier cells");

        auto emit_radicals = [&](const std::vector<RadicalDiag>& diags, const char* title) {
            report << "## " << title << "\n\n";
            if (diags.empty()) { report << "No multi-leaf cell was available for a local radical audit.\n\n"; return; }
            report << "Each support is the exact set of primes `ell` in the displayed cell for which both "
                   << "folded prefixes and both full quotient rows vanish, after stripping `q` and coefficient exceptions.\n\n";
            for (size_t i = 0; i < diags.size(); ++i) {
                const auto& d = diags[i];
                const auto& c = d.cell;
                int delta = c.rows.front().delta_state();
                report << "### Cell " << (i+1) << ": sign " << c.sign << ", state (" << c.q << ',' << c.t
                       << "), interval [" << c.lo() << ',' << c.hi() << ")\n\n```text\n";
                report << "K/J6/Js/H                 = " << c.K << '/' << c.J6 << '/' << c.Js << '/' << c.H << "\n";
                report << "Delta_sigma               = " << delta << " = " << factor_string(delta, primes) << "\n";
                report << "local radical support size= " << d.support.size() << "\n";
                report << "phantom support size      = " << d.phantoms.size() << "\n";
                report << "log radical               = " << d.log_radical << "\n";
                report << "log actual cell product   = " << d.log_actual << "\n";
                report << "logR/(K log T)            = " << (c.K ? d.log_radical/(c.K*std::log(static_cast<double>(c.T))) : 0) << "\n";
                report << "logR/(H log T)            = " << d.log_radical/(std::max(1,c.H)*std::log(static_cast<double>(c.T))) << "\n";
                report << "all actual primes supported= " << (d.all_actual_supported?"yes":"NO") << "\n";
                report << "actual primes             =";
                for (const Row& r : c.rows) report << ' ' << r.p;
                report << "\nsupport primes            =";
                for (int p : d.support) report << ' ' << p;
                report << "\nphantoms                  =";
                for (int p : d.phantoms) report << ' ' << p;
                report << "\nprefix valuations (cap3) =";
                for (auto [p,v6,vs,vg] : d.valuations)
                    report << " p" << p << ":(" << v6 << ',' << vs << ";g=" << vg << ')';
                report << "\n```\n\n";
            }
        };
        emit_radicals(radical_actual, "Exact stripped local radicals: actual cells");
        emit_radicals(radical_ambient, "Exact stripped local radicals: overcarrier cells");

        auto emit_corr = [&](const std::map<std::string,double>& corr, size_t n, const char* title) {
            report << "## " << title << "\n\n";
            report << "Pearson correlations use `log K` over every translated multi-leaf cell (`n=" << n << "`).\n\n";
            report << "| fingerprint | correlation with log K |\n|:---|---:|\n";
            for (const auto& [name,value] : corr) {
                report << '|' << name << '|';
                if (std::isnan(value)) report << "NA|\n";
                else report << value << "|\n";
            }
            report << '\n';
        };
        emit_corr(corr_actual, features_actual.size(), "Auxiliary/invariant fingerprints: actual cells");
        emit_corr(corr_ambient, features_ambient.size(), "Auxiliary/invariant fingerprints: overcarrier cells");

        uint64_t clean_actual = clean_count(enumeration.plus_actual_marked) +
                                clean_count(enumeration.minus_actual_marked);
        bool actual_rich = std::any_of(metrics_actual.begin(), metrics_actual.end(), [](const MetricRow& r){
            return r.agg.rich_num > 0;
        });
        bool actual_sat = std::any_of(metrics_actual.begin(), metrics_actual.end(), [](const MetricRow& r){
            return r.agg.sat_num > 0;
        });
        report << "## Finite counterexample protocol\n\n```text\n";
        report << "clean actual marked leaves       = " << clean_actual << "\n";
        report << "actual Q4190 rich cell observed   = " << (actual_rich?"YES":"no") << "\n";
        report << "actual saturated-rich cell        = " << (actual_sat?"YES":"no") << "\n";
        report << "singleton saturated-rich witness  = "
               << (std::any_of(features_actual.begin(), features_actual.end(), [](const CellFeature& f){
                    return f.sector==1 && f.K>f.H && f.J6>f.H && f.Js>f.H;
                  }) ? "YES" : "no") << "\n";
        report << "prefix carrier support failures   = "
               << std::count_if(radical_actual.begin(), radical_actual.end(), [](const RadicalDiag& d){return !d.all_actual_supported;})
               << "\n```\n\n";
        report << "A finite saturated-rich witness is not by itself a counterexample to TSJMA/SRT, whose constants and exponents are asymptotic. "
               << "Conversely, an event-free range cannot identify an exponent. The report therefore distinguishes literal witnesses, censored zeros, and theorem violations.\n\n";

        report << "## Strongest data statement and proof boundary\n\n";
        if (clean_actual == 0) {
            report << "The exhaustive clean scan proves the finite statement\n\n";
            report << "```text\nFor every prime q <= " << opt.qmax
                   << ", every actual selected t in Z_q, every prime p>=17,\n"
                   << "and both exact raw q6 signs with minus-first assignment, the marked fibre is empty.\n```\n\n";
            report << "Thus the literal TSJMA and SRT numerators and denominators are both zero on every complete dyadic shell in the scan; no scaling exponent is identifiable. "
                   << "The unselected lower-complete overcarrier remains non-vacuous and diagnoses where the upper selected mark deletes candidate rows.\n\n";
        } else {
            report << "Nonempty actual fibres occur. The exact translated tables above are the finite law; exponent fits use only positive observations and remain evidence rather than proof.\n\n";
        }
        report << "The symbolic invariants tested are the common folded depth `min(J6,Js)`, the state determinant factor\n\n";
        report << "```text\nDelta_+ = 6q+7t,        Delta_- = 7t+6,\n"
               << "det(row_p,row_p')=(p-p')*Delta_sigma,\n"
               << "SNF pair=(d1, |p-p'|*|Delta_sigma|/d1), d1=gcd(r,a,p'-p),\n"
               << "and the auxiliary vanishing mask of the primitive bilinear carrier.\n```\n\n";
        report << "Among nonempty multi-leaf data, the correlation table identifies the most predictive finite fingerprint. "
               << "No such correlation is promoted to a theorem. The strongest precise asymptotic target remains the occurrence-preserving combined selected saturated rich tail `CSRT`, separately normalized by singleton and non-singleton actual masses.\n";
        report.close();

        // Compact machine-readable summary without external JSON dependency.
        std::ofstream json(opt.out_dir + "/summary.json");
        json << "{\n";
        json << "  \"qmax\": " << opt.qmax << ",\n";
        json << "  \"prime_count\": " << primes.size() << ",\n";
        json << "  \"recurrence_steps\": " << global_steps.load() << ",\n";
        json << "  \"selected_states\": " << selected_states << ",\n";
        json << "  \"singleton_primes\": " << singleton_primes << ",\n";
        json << "  \"max_z\": " << max_z << ",\n";
        json << "  \"candidate_plus_marked\": " << enumeration.plus_candidate_marked.size() << ",\n";
        json << "  \"candidate_minus_marked\": " << enumeration.minus_candidate_marked.size() << ",\n";
        json << "  \"actual_plus_marked\": " << enumeration.plus_actual_marked.size() << ",\n";
        json << "  \"actual_minus_marked\": " << enumeration.minus_actual_marked.size() << ",\n";
        json << "  \"clean_actual\": " << clean_actual << ",\n";
        json << "  \"actual_rich_observed\": " << (actual_rich?"true":"false") << ",\n";
        json << "  \"actual_saturated_rich_observed\": " << (actual_sat?"true":"false") << ",\n";
        json << "  \"zero_seconds\": " << zero_seconds << ",\n";
        json << "  \"total_seconds\": " << total_timer.seconds() << "\n";
        json << "}\n";
        json.close();

        std::cout << "Q4204_RESULT_BEGIN\n";
        std::ifstream in(opt.out_dir + "/report.md");
        std::cout << in.rdbuf();
        std::cout << "\nQ4204_RESULT_END\n";
        std::cerr << "Q4204_DONE total_seconds=" << total_timer.seconds() << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Q4204_FATAL " << e.what() << "\n";
        return 2;
    }
}
