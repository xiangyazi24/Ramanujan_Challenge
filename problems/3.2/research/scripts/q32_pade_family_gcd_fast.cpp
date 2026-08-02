// Exact growing-family Padé gcd scan using one fraction-free elimination.
//
// For fixed H, reverse the bottom Newton multiplication rows and form
//
//   K^+_{i,l} = M_{H-i,l},  0 <= i < H, 0 <= l <= H.
//
// If b=H-a, the denominator q_{H,a} is the primitive integer relation
// among the first b+1 columns of the first b rows.  A single Bareiss
// elimination of K^+ gives a unit upper matrix U.  The monic rational
// relation for every b is obtained by back substitution in
//
//   U[0:b,0:b+1] q = 0, q_b=1.
//
// Clearing denominators and coefficient gcd gives the exact primitive
// denominator.  This replaces a separate expansion of b+1 large cofactors
// for every numerator degree a by one O(H^3) fraction-free elimination plus
// O(A H^2) back substitutions.

#include <gmpxx.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Big = mpz_class;
using Rat = mpq_class;

static Big abs_big(const Big& x) { return x >= 0 ? x : -x; }

static Big gcd_big(const Big& x, const Big& y) {
  Big a = abs_big(x), b = abs_big(y), g;
  mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return g;
}

static Big lcm_big(const Big& x, const Big& y) {
  if (x == 0 || y == 0) return 0;
  Big a = abs_big(x), b = abs_big(y), l;
  mpz_lcm(l.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return l;
}

static std::size_t bit_length(const Big& x) {
  Big a = abs_big(x);
  if (a == 0) return 0;
  return mpz_sizeinbase(a.get_mpz_t(), 2);
}

static unsigned long mod_ui(const Big& x, unsigned long p) {
  return mpz_fdiv_ui(x.get_mpz_t(), p);
}

static std::vector<Big> apery_values(int maximum_index) {
  if (maximum_index < 1) return {Big(1)};
  std::vector<Big> values(maximum_index + 1);
  values[0] = 1;
  values[1] = 5;
  for (int n = 1; n < maximum_index; ++n) {
    const long long nn = n;
    const Big polynomial =
        Big(34) * nn * nn * nn + Big(51) * nn * nn + Big(27) * nn + 5;
    const Big numerator = polynomial * values[n] - Big(nn) * nn * nn * values[n - 1];
    const Big denominator = Big(n + 1) * (n + 1) * (n + 1);
    assert(numerator % denominator == 0);
    values[n + 1] = numerator / denominator;
  }
  return values;
}

static std::vector<std::vector<Big>> pascal_table(int height) {
  std::vector<std::vector<Big>> choose(
      height + 1, std::vector<Big>(height + 1, 0));
  choose[0][0] = 1;
  for (int n = 1; n <= height; ++n) {
    choose[n][0] = 1;
    choose[n][n] = 1;
    for (int k = 1; k < n; ++k) {
      choose[n][k] = choose[n - 1][k - 1] + choose[n - 1][k];
    }
  }
  return choose;
}

static std::vector<std::vector<Big>> difference_table(
    const std::vector<Big>& values, int height) {
  std::vector<std::vector<Big>> diff(height + 1);
  diff[0].assign(values.begin(), values.begin() + height + 1);
  for (int order = 1; order <= height; ++order) {
    diff[order].resize(height + 1 - order);
    for (int start = 0; start + order <= height; ++start) {
      diff[order][start] =
          diff[order - 1][start + 1] - diff[order - 1][start];
    }
  }
  return diff;
}

static std::vector<std::vector<Big>> multiplication_matrix(
    int height, const std::vector<std::vector<Big>>& choose,
    const std::vector<std::vector<Big>>& diff) {
  std::vector<std::vector<Big>> matrix(
      height + 1, std::vector<Big>(height + 1, 0));
  for (int k = 0; k <= height; ++k) {
    for (int ell = 0; ell <= k; ++ell) {
      matrix[k][ell] = choose[k][ell] * diff[k - ell][ell];
    }
  }
  return matrix;
}

static std::vector<Big> binomial_row(int n, int maximum_k) {
  std::vector<Big> result(maximum_k + 1);
  result[0] = 1;
  for (int k = 0; k < maximum_k; ++k) {
    const Big numerator = result[k] * (n - k);
    assert(numerator % (k + 1) == 0);
    result[k + 1] = numerator / (k + 1);
  }
  return result;
}

static void bareiss_eliminate(std::vector<std::vector<Big>>& matrix) {
  const int rows = static_cast<int>(matrix.size());
  if (rows == 0) return;
  const int columns = static_cast<int>(matrix[0].size());
  assert(columns == rows + 1);

  Big previous = 1;
  for (int pivot_index = 0; pivot_index < rows; ++pivot_index) {
    const Big pivot = matrix[pivot_index][pivot_index];
    if (pivot == 0) {
      std::ostringstream message;
      message << "zero nested pivot at index " << pivot_index;
      throw std::runtime_error(message.str());
    }
    if (pivot_index + 1 == rows) break;

    for (int row = pivot_index + 1; row < rows; ++row) {
      const Big below = matrix[row][pivot_index];
      for (int column = pivot_index + 1; column < columns; ++column) {
        const Big numerator =
            matrix[row][column] * pivot - below * matrix[pivot_index][column];
        assert(numerator % previous == 0);
        matrix[row][column] = numerator / previous;
      }
      matrix[row][pivot_index] = 0;
    }
    previous = pivot;
  }
}

static std::vector<Big> primitive_relation(
    const std::vector<std::vector<Big>>& eliminated, int b) {
  std::vector<Rat> rational(b + 1);
  rational[b] = 1;
  for (int row = b - 1; row >= 0; --row) {
    Rat sum = 0;
    for (int column = row + 1; column <= b; ++column) {
      Rat upper(eliminated[row][column], eliminated[row][row]);
      upper.canonicalize();
      sum += upper * rational[column];
    }
    rational[row] = -sum;
    rational[row].canonicalize();
  }

  Big denominator_lcm = 1;
  for (const Rat& value : rational) {
    denominator_lcm = lcm_big(denominator_lcm, value.get_den());
  }

  std::vector<Big> relation(b + 1);
  Big common = 0;
  for (int index = 0; index <= b; ++index) {
    relation[index] =
        rational[index].get_num() * (denominator_lcm / rational[index].get_den());
    common = gcd_big(common, relation[index]);
  }
  assert(common != 0);
  for (Big& value : relation) value /= common;
  if (relation.back() < 0) {
    for (Big& value : relation) value = -value;
  }

  Big check = 0;
  for (const Big& value : relation) check = gcd_big(check, value);
  assert(check == 1);
  return relation;
}

static Big evaluate_newton(const std::vector<Big>& coefficients, int x) {
  const int maximum = std::min<int>(x, coefficients.size() - 1);
  Big choose = 1;
  Big result = 0;
  for (int k = 0; k <= maximum; ++k) {
    result += coefficients[k] * choose;
    if (k < maximum) {
      const Big numerator = choose * (x - k);
      assert(numerator % (k + 1) == 0);
      choose = numerator / (k + 1);
    }
  }
  return result;
}

static int ceil_two_thirds(int height) {
  int cutoff = 0;
  const long long target = 1LL * height * height;
  while (1LL * cutoff * cutoff * cutoff < target) ++cutoff;
  return cutoff;
}

static std::vector<bool> prime_sieve(int limit) {
  std::vector<bool> prime(limit + 1, true);
  if (limit >= 0) prime[0] = false;
  if (limit >= 1) prime[1] = false;
  for (int divisor = 2; 1LL * divisor * divisor <= limit; ++divisor) {
    if (!prime[divisor]) continue;
    for (int multiple = divisor * divisor; multiple <= limit;
         multiple += divisor) {
      prime[multiple] = false;
    }
  }
  return prime;
}

static std::string join_ints(const std::vector<int>& values) {
  if (values.empty()) return "-";
  std::ostringstream out;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) out << ',';
    out << values[i];
  }
  return out.str();
}

struct HeightResult {
  int height;
  int cutoff;
  Big gcd;
  Big p01_gcd;
  Big apery_n;
  Big binomial_carrier;
  Big direct_core;
  Big gcd_with_apery_n;
  Big gcd_with_carrier;
  int maximum_zero_count;
  std::vector<int> target_primes;
  std::vector<int> candidate_gcd_primes;
  std::vector<int> expected_primes;
  std::vector<std::size_t> value_bits;
};

static HeightResult audit_height(int height) {
  const int n = 3 * height + 1;
  const int cutoff = std::min(height, ceil_two_thirds(height));
  const std::vector<Big> apery = apery_values(n);
  const auto choose = pascal_table(height);
  const auto diff = difference_table(apery, height);
  const auto multiplication = multiplication_matrix(height, choose, diff);
  const std::vector<Big> choose_n = binomial_row(n, height + 1);

  Big apery_lcm = 1;
  for (int s = 0; s <= height; ++s) {
    apery_lcm = lcm_big(apery_lcm, apery[s]);
  }

  std::vector<Big> evaluation_row(height + 1, 0);
  for (int ell = 0; ell <= height; ++ell) {
    for (int k = ell; k <= height; ++k) {
      evaluation_row[ell] += choose_n[k] * multiplication[k][ell];
    }
  }

  std::vector<std::vector<Big>> reversed(
      height, std::vector<Big>(height + 1, 0));
  for (int row = 0; row < height; ++row) {
    const int original_row = height - row;
    for (int column = 0; column <= height; ++column) {
      reversed[row][column] = multiplication[original_row][column];
    }
  }
  bareiss_eliminate(reversed);

  Big common = 0;
  Big value_a0 = 0;
  Big value_a1 = 0;
  std::vector<std::size_t> value_bits;
  value_bits.reserve(cutoff + 1);

  for (int a = 0; a <= cutoff; ++a) {
    const int b = height - a;
    const std::vector<Big> denominator = primitive_relation(reversed, b);

    std::vector<Big> numerator(a + 1, 0);
    for (int k = 0; k <= height; ++k) {
      Big coefficient = 0;
      for (int ell = 0; ell <= std::min(k, b); ++ell) {
        coefficient += multiplication[k][ell] * denominator[ell];
      }
      if (k <= a) {
        numerator[k] = coefficient;
      } else {
        assert(coefficient == 0);
      }
    }

    Big value_from_row = 0;
    for (int ell = 0; ell <= b; ++ell) {
      value_from_row += evaluation_row[ell] * denominator[ell];
    }
    Big value_from_numerator = 0;
    for (int k = 0; k <= a; ++k) {
      value_from_numerator += choose_n[k] * numerator[k];
    }
    assert(value_from_row == value_from_numerator);
    const Big value = value_from_row;
    assert(value != 0);

    if (a == 0 || a == 1 || a == cutoff) {
      for (int s = 0; s <= height; ++s) {
        assert(evaluate_newton(numerator, s) ==
               apery[s] * evaluate_newton(denominator, s));
      }
    }

    if (a == 0) {
      value_a0 = value;
      assert(abs_big(value_a0) == apery_lcm);
    }
    if (a == 1) value_a1 = value;

    common = gcd_big(common, value);
    value_bits.push_back(bit_length(value));
    std::cerr << "H=" << height << " a=" << a << "/" << cutoff
              << " value_bits=" << bit_length(value)
              << " gcd_bits=" << bit_length(common) << "\n";
  }

  // Independent Q850 check for the a=1 primitive value.
  if (cutoff >= 1) {
    Big X = 0, Y = 0;
    for (int s = 0; s <= height; ++s) {
      Big term = choose[height][s] * (apery_lcm / apery[s]);
      if (s & 1) term = -term;
      X += term;
      Y += Big(s) * term;
    }
    const Big h = gcd_big(X, Y);
    assert(h != 0);
    Big multiplier = 1;
    for (int s = 0; s <= height; ++s) {
      const Big numerator = Y - Big(s) * X;
      assert(numerator % h == 0);
      const Big z = numerator / h;
      multiplier = lcm_big(multiplier, apery[s] / gcd_big(apery[s], z));
    }
    const Big root_value = (Y - Big(n) * X) / h;
    const Big expected_a1 = multiplier * root_value;
    assert(abs_big(expected_a1) == abs_big(value_a1));
  }

  const Big p01_gcd = gcd_big(value_a0, value_a1);
  assert(common != 0);
  assert(apery_lcm % common == 0);
  assert(p01_gcd % common == 0);

  const std::vector<bool> prime = prime_sieve(n);
  int maximum_zero_count = 0;
  std::vector<int> targets;
  std::vector<int> candidate_support;
  std::vector<int> expected;
  for (int node = 0; node <= height; ++node) {
    const int p = n - node;
    if (!prime[p]) continue;
    int zero_count = 0;
    for (int s = 0; s <= height; ++s) {
      if (mod_ui(apery[s], p) == 0) ++zero_count;
    }
    maximum_zero_count = std::max(maximum_zero_count, zero_count);
    const bool target = mod_ui(apery[node], p) == 0;
    const bool predicted = target || zero_count > cutoff;
    const bool observed = mod_ui(common, p) == 0;
    if (target) targets.push_back(p);
    if (predicted) expected.push_back(p);
    if (observed) candidate_support.push_back(p);
    assert(observed == predicted);
  }
  assert(candidate_support == expected);
  if (cutoff >= maximum_zero_count) assert(candidate_support == targets);

  Big carrier;
  mpz_bin_uiui(carrier.get_mpz_t(), static_cast<unsigned long>(n),
               static_cast<unsigned long>(height + 1));
  const Big direct_core = gcd_big(apery[n], carrier);

  HeightResult result;
  result.height = height;
  result.cutoff = cutoff;
  result.gcd = abs_big(common);
  result.p01_gcd = p01_gcd;
  result.apery_n = apery[n];
  result.binomial_carrier = carrier;
  result.direct_core = direct_core;
  result.gcd_with_apery_n = gcd_big(common, apery[n]);
  result.gcd_with_carrier = gcd_big(common, carrier);
  result.maximum_zero_count = maximum_zero_count;
  result.target_primes = targets;
  result.candidate_gcd_primes = candidate_support;
  result.expected_primes = expected;
  result.value_bits = value_bits;
  return result;
}

static std::string join_sizes(const std::vector<std::size_t>& values) {
  std::ostringstream out;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) out << ',';
    out << values[i];
  }
  return out.str();
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " H [H ...]\n";
    return 2;
  }
  try {
    for (int argument = 1; argument < argc; ++argument) {
      const int height = std::stoi(argv[argument]);
      if (height < 2) throw std::runtime_error("height must be at least 2");
      const HeightResult result = audit_height(height);
      std::cout << "RESULT"
                << " H=" << result.height
                << " A=" << result.cutoff
                << " G=" << result.gcd
                << " G_bits=" << bit_length(result.gcd)
                << " P01=" << result.p01_gcd
                << " P01_bits=" << bit_length(result.p01_gcd)
                << " max_z=" << result.maximum_zero_count
                << " targets=" << join_ints(result.target_primes)
                << " candidate_support=" << join_ints(result.candidate_gcd_primes)
                << " expected_support=" << join_ints(result.expected_primes)
                << " direct_core=" << result.direct_core
                << " gcd_G_An=" << result.gcd_with_apery_n
                << " gcd_G_B=" << result.gcd_with_carrier
                << " value_bits=" << join_sizes(result.value_bits)
                << "\n";
    }
  } catch (const std::exception& error) {
    std::cerr << "fatal: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
