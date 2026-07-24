// Exact modular scan of the unavoidable a=1 multiplier pollution.
//
// For interpolation height H and a candidate prime 2H < p <= 3H+1,
// Q850 gives
//
//   m_H = lcm_{0<=s<=H} A_s / gcd(A_s,z_H(s)),
//
// with z_H primitive linear.  If A_s vanishes modulo p at two distinct
// nodes s<=H, then p divides m_H: a primitive linear polynomial modulo
// p cannot vanish at both nodes.  Conversely, p can divide m_H only if at
// least one A_s vanishes modulo p in the prefix.
//
// The scan records, at every H:
//   possible    = primes having at least one prefix zero (an upper bound);
//   unavoidable = primes having at least two prefix zeros (a lower bound);
//   target      = primes vanishing at the moving node 3H+1-p.
//
// It does not construct the enormous primitive Padé coefficients.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using i64 = std::int64_t;

static i64 mul(i64 left, i64 right, i64 prime) {
  return static_cast<i64>((__int128)left * right % prime);
}

static void report_dyadic(const std::string& name,
                          const std::vector<double>& mass,
                          const std::vector<int>& count, int height_max) {
  std::cout << name << "\n";
  for (int cutoff : {100, 1000, 10000, 50000, height_max}) {
    if (cutoff > height_max) continue;
    double best = 0.0;
    int at = 0;
    for (int height = std::max(2, cutoff / 2); height <= cutoff; ++height) {
      const double ratio = mass[height] / height;
      if (ratio > best) {
        best = ratio;
        at = height;
      }
    }
    std::cout << "dyadic<= " << cutoff << " best=" << std::setprecision(10)
              << best << " at=" << at
              << " count=" << (at ? count[at] : 0) << "\n";
  }
  for (int height : {100, 300, 1000, 3000, 10000, height_max}) {
    if (height > height_max) continue;
    std::cout << "H=" << height << " mass/H=" << std::setprecision(10)
              << mass[height] / height << " count=" << count[height] << "\n";
  }
}

int main(int argc, char** argv) {
  const int height_max = argc > 1 ? std::stoi(argv[1]) : 10000;
  if (height_max < 2) {
    std::cerr << "height_max must be at least 2\n";
    return 2;
  }
  const int prime_limit = 3 * height_max + 1;
  std::vector<bool> is_prime(prime_limit + 1, true);
  is_prime[0] = is_prime[1] = false;
  for (int divisor = 2; divisor * divisor <= prime_limit; ++divisor) {
    if (!is_prime[divisor]) continue;
    for (int multiple = divisor * divisor; multiple <= prime_limit;
         multiple += divisor) {
      is_prime[multiple] = false;
    }
  }

  std::vector<double> possible_mass(height_max + 1, 0.0);
  std::vector<double> unavoidable_mass(height_max + 1, 0.0);
  std::vector<double> target_mass(height_max + 1, 0.0);
  std::vector<int> possible_count(height_max + 1, 0);
  std::vector<int> unavoidable_count(height_max + 1, 0);
  std::vector<int> target_count(height_max + 1, 0);

  for (int prime = 3; prime <= prime_limit; ++prime) {
    if (!is_prime[prime]) continue;
    const int first_height = std::max(2, (prime + 1) / 3);
    const int last_height = std::min(height_max, (prime - 1) / 2);
    if (first_height > last_height) continue;

    std::vector<i64> inverse(last_height + 1);
    std::vector<i64> apery(last_height + 1);
    inverse[1] = 1;
    for (int index = 2; index <= last_height; ++index) {
      inverse[index] =
          prime - mul(prime / index, inverse[prime % index], prime);
    }
    apery[0] = 1;
    if (last_height >= 1) apery[1] = 5 % prime;
    for (i64 index = 1; index < last_height; ++index) {
      const i64 square = mul(index, index, prime);
      const i64 cube = mul(square, index, prime);
      const i64 polynomial =
          (34 * cube + 51 * square + 27 * index + 5) % prime;
      i64 right =
          (mul(polynomial, apery[index], prime) -
           mul(cube, apery[index - 1], prime)) %
          prime;
      if (right < 0) right += prime;
      const i64 reciprocal = inverse[index + 1];
      apery[index + 1] =
          mul(right, mul(reciprocal, mul(reciprocal, reciprocal, prime),
                         prime),
              prime);
    }

    int zero_count = 0;
    const double logarithm = std::log(static_cast<double>(prime));
    for (int height = 0; height <= last_height; ++height) {
      if (apery[height] == 0) ++zero_count;
      if (height < first_height) continue;
      if (zero_count >= 1) {
        possible_mass[height] += logarithm;
        ++possible_count[height];
      }
      if (zero_count >= 2) {
        unavoidable_mass[height] += logarithm;
        ++unavoidable_count[height];
      }
      const int moving_node = 3 * height + 1 - prime;
      if (apery[moving_node] == 0) {
        target_mass[height] += logarithm;
        ++target_count[height];
      }
    }
  }

  std::cout << "height_max=" << height_max
            << " prime_limit=" << prime_limit << "\n";
  report_dyadic("possible_multiplier_upper", possible_mass, possible_count,
                height_max);
  report_dyadic("unavoidable_multiplier_lower", unavoidable_mass,
                unavoidable_count, height_max);
  report_dyadic("direct_target", target_mass, target_count, height_max);
}
