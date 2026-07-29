#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

using i64 = std::int64_t;

static i64 mul(i64 a, i64 b, i64 p) {
  return static_cast<i64>((__int128)a * b % p);
}

int main(int argc, char** argv) {
  const int limit = argc > 1 ? std::stoi(argv[1]) : 100000;
  std::vector<bool> is_prime(limit + 1, true);
  is_prime[0] = is_prime[1] = false;
  for (int d = 2; d * d <= limit; ++d) {
    if (is_prime[d]) {
      for (int m = d * d; m <= limit; m += d) is_prime[m] = false;
    }
  }

  std::vector<double> lower(limit + 1, 0.0);
  std::vector<double> q1(limit + 1, 0.0);
  std::vector<double> q1_direct(limit + 1, 0.0);
  std::vector<double> q1_reflected(limit + 1, 0.0);
  std::vector<int> lower_count(limit + 1, 0);
  std::vector<int> q1_count(limit + 1, 0);
  std::vector<int> q1_direct_count(limit + 1, 0);
  std::vector<int> q1_reflected_count(limit + 1, 0);
  std::vector<std::vector<int>> q1_hits(limit + 1);
  std::size_t zero_pairs = 0;

  for (int p = 3; p <= limit; ++p) {
    if (!is_prime[p]) continue;
    std::vector<i64> inv(p), b(p);
    inv[1] = 1;
    for (int k = 2; k < p; ++k) {
      inv[k] = p - mul(p / k, inv[p % k], p);
    }
    b[0] = 1;
    if (p > 1) b[1] = 5 % p;
    for (i64 n = 1; n + 1 < p; ++n) {
      const i64 n2 = mul(n, n, p);
      const i64 n3 = mul(n2, n, p);
      i64 poly = (34 * n3 + 51 * n2 + 27 * n + 5) % p;
      i64 rhs = (mul(poly, b[n], p) - mul(n3, b[n - 1], p)) % p;
      if (rhs < 0) rhs += p;
      const i64 z = inv[n + 1];
      b[n + 1] = mul(rhs, mul(z, mul(z, z, p), p), p);
    }

    const double lp = std::log(static_cast<double>(p));
    for (int j = 0; j < p; ++j) {
      if (b[j] != 0) continue;
      ++zero_pairs;
      const int qmax = std::min(p - 1, (limit - j) / p);
      for (int q = 1; q <= qmax; ++q) {
        const int n = q * p + j;
        lower[n] += lp;
        ++lower_count[n];
        if (q == 1) {
          q1[n] += lp;
          ++q1_count[n];
          q1_hits[n].push_back(p);
          if (2 * j <= p - 1) {
            q1_direct[n] += lp;
            ++q1_direct_count[n];
          } else {
            q1_reflected[n] += lp;
            ++q1_reflected_count[n];
          }
        }
      }
    }
  }

  auto report = [&](const char* name, const std::vector<double>& mass,
                    const std::vector<int>& count) {
    std::vector<int> order(limit - 9);
    std::iota(order.begin(), order.end(), 10);
    std::partial_sort(
        order.begin(), order.begin() + std::min<int>(20, order.size()), order.end(),
        [&](int a, int b) { return mass[a] / a > mass[b] / b; });
    std::cout << name << " top mass/n\n";
    for (int k = 0; k < std::min<int>(20, order.size()); ++k) {
      const int n = order[k];
      std::cout << n << " count=" << count[n] << " mass=" << std::setprecision(8)
                << mass[n] << " ratio=" << mass[n] / n << "\n";
    }
    for (int cutoff : {100, 1000, 10000, 100000, limit}) {
      if (cutoff > limit) continue;
      double best = 0.0;
      int at = 0;
      for (int n = std::max(10, cutoff / 2); n <= cutoff; ++n) {
        if (mass[n] / n > best) {
          best = mass[n] / n;
          at = n;
        }
      }
      std::cout << "dyadic<= " << cutoff << " best=" << best << " at=" << at
                << " count=" << (at ? count[at] : 0) << "\n";
    }
    const int max_count =
        *std::max_element(count.begin() + std::min(10, limit), count.end());
    std::cout << "max_count=" << max_count << " at";
    for (int n = 10; n <= limit; ++n) {
      if (count[n] == max_count) std::cout << " " << n;
    }
    std::cout << "\n";
  };

  std::cout << "limit=" << limit << " zero_pairs=" << zero_pairs << "\n";
  report("lower", lower, lower_count);
  report("q1", q1, q1_count);
  report("q1_direct", q1_direct, q1_direct_count);
  report("q1_reflected", q1_reflected, q1_reflected_count);
  const int q1_max =
      *std::max_element(q1_count.begin() + std::min(10, limit), q1_count.end());
  std::cout << "q1 max-hit tuples\n";
  for (int n = 10; n <= limit; ++n) {
    if (q1_count[n] != q1_max) continue;
    std::cout << "n=" << n;
    for (int p : q1_hits[n]) {
      const int r = n - p;
      const int j = std::min(r, p - 1 - r);
      const int kernel = std::gcd(p - 1, n - 1);
      const int order = (p - 1) / kernel;
      std::cout << " (p=" << p << ",r=" << r << ",j=" << j
                << ",branch=" << (j == r ? "D" : "R")
                << ",g=" << kernel << ",d=" << order << ")";
    }
    std::cout << "\n";
  }
}
