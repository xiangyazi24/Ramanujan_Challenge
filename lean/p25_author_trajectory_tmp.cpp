#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <tuple>
#include <vector>

using C = std::complex<double>;

static C ipow(C z, int exponent) {
  if (exponent < 0) return C(1) / ipow(z, -exponent);
  C answer(1);
  while (exponent) {
    if (exponent & 1) answer *= z;
    z *= z;
    exponent >>= 1;
  }
  return answer;
}

// Return the three roots of a monic cubic.
static std::array<C, 3> cubic_roots(double a, double b, double c) {
  const C delta0 = a * a - 3 * b;
  const C delta1 = 2 * a * a * a - 9 * a * b + 27 * c;
  const C discr = delta1 * delta1 - 4.0 * delta0 * delta0 * delta0;
  C plus = (delta1 + std::sqrt(discr)) / 2.0;
  C minus = (delta1 - std::sqrt(discr)) / 2.0;
  C q = std::pow(std::abs(plus) >= std::abs(minus) ? plus : minus, 1.0 / 3.0);
  if (std::abs(q) < 1e-14) {
    // The only case relevant here is the triple root.
    return {C(-a / 3), C(-a / 3), C(-a / 3)};
  }
  const C omega(-0.5, std::sqrt(3.0) / 2.0);
  std::array<C, 3> roots;
  C u(1);
  for (int index = 0; index < 3; ++index) {
    roots[index] = -(a + u * q + delta0 / (u * q)) / 3.0;
    u *= omega;
  }
  return roots;
}

static std::array<double, 5> elementary(const std::vector<int>& values) {
  // Coefficients of prod (theta + value), in descending order.
  std::array<double, 5> result{1, 0, 0, 0, 0};
  int degree = 0;
  for (int value : values) {
    ++degree;
    for (int index = degree; index >= 1; --index)
      result[index] = result[index] + value * result[index - 1];
  }
  return result;
}

struct Hit {
  double score;
  std::array<int, 7> trajectory;
  double invariant1;
  double invariant2;
  std::array<C, 3> eigenvalues;
  bool operator<(const Hit& other) const { return score < other.score; }
};

int main(int argc, char** argv) {
  const int bound = argc > 1 ? std::stoi(argv[1]) : 8;
  const double target1 = 1225.0;
  const double target2 = 42875.0;
  std::priority_queue<Hit> best;
  unsigned long long tested = 0;
  for (int a0 = -bound; a0 <= bound; ++a0)
  for (int a1 = a0; a1 <= bound; ++a1)
  for (int a2 = a1; a2 <= bound; ++a2)
  for (int a3 = a2; a3 <= bound; ++a3) {
    const std::vector<int> av{a0, a1, a2, a3};
    const auto pa = elementary(av);
    for (int b0 = -bound; b0 <= bound; ++b0)
    for (int b1 = b0; b1 <= bound; ++b1)
    for (int b2 = b1; b2 <= bound; ++b2) {
      // This is the non-cancellation condition used by the author's search.
      if (std::find(av.begin(), av.end(), b0) != av.end() ||
          std::find(av.begin(), av.end(), b1) != av.end() ||
          std::find(av.begin(), av.end(), b2) != av.end()) continue;
      ++tested;
      const std::vector<int> bv{0, b0, b1, b2};
      const auto pb = elementary(bv);
      const double lead = pb[1] - pa[1];
      if (std::abs(lead) < 0.5) continue;
      const double qa = (pb[2] - pa[2]) / lead;
      const double qb = (pb[3] - pa[3]) / lead;
      const double qc = (pb[4] - pa[4]) / lead;
      const auto roots = cubic_roots(qa, qb, qc);
      std::array<C, 3> eigenvalues;
      bool valid = true;
      for (int root_index = 0; root_index < 3; ++root_index) {
        C value(1);
        for (int step : av)
          if (step) value *= ipow(C(1) + roots[root_index] / double(step), step);
        for (int step : std::array<int, 3>{b0, b1, b2})
          if (step) value *= ipow(C(1) + roots[root_index] / double(step), -step);
        eigenvalues[root_index] = value;
        valid &= std::isfinite(value.real()) && std::isfinite(value.imag());
      }
      if (!valid) continue;
      const C e1 = eigenvalues[0] + eigenvalues[1] + eigenvalues[2];
      const C e2 = eigenvalues[0] * eigenvalues[1] +
                   eigenvalues[0] * eigenvalues[2] +
                   eigenvalues[1] * eigenvalues[2];
      const C e3 = eigenvalues[0] * eigenvalues[1] * eigenvalues[2];
      if (std::abs(e3) < 1e-100) continue;
      const C inv1 = e1 * e2 / e3;
      const C inv2 = e1 * e1 * e1 / e3;
      if (std::abs(inv1.imag()) > 1e-5 * std::max(1.0, std::abs(inv1.real())) ||
          std::abs(inv2.imag()) > 1e-5 * std::max(1.0, std::abs(inv2.real())) ||
          inv1.real() <= 0 || inv2.real() <= 0) continue;
      const double score = std::abs(std::log(inv1.real() / target1)) +
                           std::abs(std::log(inv2.real() / target2));
      Hit hit{score, {a0, a1, a2, a3, b0, b1, b2}, inv1.real(),
              inv2.real(), eigenvalues};
      if (best.size() < 100) best.push(hit);
      else if (score < best.top().score) { best.pop(); best.push(hit); }
    }
  }
  std::vector<Hit> output;
  while (!best.empty()) { output.push_back(best.top()); best.pop(); }
  std::sort(output.begin(), output.end());
  std::cout << "tested " << tested << "\n" << std::setprecision(16);
  for (const auto& hit : output) {
    std::cout << hit.score << " (";
    for (int index = 0; index < 7; ++index)
      std::cout << (index ? "," : "") << hit.trajectory[index];
    std::cout << ") " << hit.invariant1 << " " << hit.invariant2 << " eig";
    for (const C value : hit.eigenvalues) std::cout << " " << value;
    std::cout << "\n";
  }
}
