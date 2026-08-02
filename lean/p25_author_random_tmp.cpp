#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
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

static std::array<C, 3> cubic_roots(double a, double b, double c) {
  const C delta0 = a * a - 3 * b;
  const C delta1 = 2 * a * a * a - 9 * a * b + 27 * c;
  const C discr = delta1 * delta1 - 4.0 * delta0 * delta0 * delta0;
  C plus = (delta1 + std::sqrt(discr)) / 2.0;
  C minus = (delta1 - std::sqrt(discr)) / 2.0;
  C q = std::pow(std::abs(plus) >= std::abs(minus) ? plus : minus,
                 1.0 / 3.0);
  if (std::abs(q) < 1e-14)
    return {C(-a / 3), C(-a / 3), C(-a / 3)};
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
  std::array<double, 5> result{1, 0, 0, 0, 0};
  int degree = 0;
  for (int value : values) {
    ++degree;
    for (int index = degree; index >= 1; --index)
      result[index] += value * result[index - 1];
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

static bool evaluate(std::array<int, 7> t, Hit& hit) {
  std::sort(t.begin(), t.begin() + 4);
  std::sort(t.begin() + 4, t.end());
  int divisor = 0;
  for (int value : t) divisor = std::gcd(divisor, std::abs(value));
  if (divisor > 1) for (int& value : t) value /= divisor;
  for (int i = 0; i < 4; ++i)
    for (int j = 4; j < 7; ++j)
      if (t[i] == t[j]) return false;
  const std::vector<int> av(t.begin(), t.begin() + 4);
  const std::vector<int> bv{0, t[4], t[5], t[6]};
  const auto pa = elementary(av);
  const auto pb = elementary(bv);
  const double lead = pb[1] - pa[1];
  if (std::abs(lead) < 0.5) return false;
  const auto roots = cubic_roots((pb[2] - pa[2]) / lead,
                                 (pb[3] - pa[3]) / lead,
                                 (pb[4] - pa[4]) / lead);
  std::array<C, 3> eigenvalues;
  for (int root_index = 0; root_index < 3; ++root_index) {
    C value(1);
    for (int step : av)
      if (step) value *= ipow(C(1) + roots[root_index] / double(step), step);
    for (int step : std::array<int, 3>{t[4], t[5], t[6]})
      if (step) value *= ipow(C(1) + roots[root_index] / double(step), -step);
    if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) return false;
    eigenvalues[root_index] = value;
  }
  const C e1 = eigenvalues[0] + eigenvalues[1] + eigenvalues[2];
  const C e2 = eigenvalues[0] * eigenvalues[1] +
               eigenvalues[0] * eigenvalues[2] +
               eigenvalues[1] * eigenvalues[2];
  const C e3 = eigenvalues[0] * eigenvalues[1] * eigenvalues[2];
  if (std::abs(e3) < 1e-100) return false;
  const C invariant1 = e1 * e2 / e3;
  const C invariant2 = e1 * e1 * e1 / e3;
  if (std::abs(invariant1.imag()) > 1e-5 * std::max(1.0, std::abs(invariant1.real())) ||
      std::abs(invariant2.imag()) > 1e-5 * std::max(1.0, std::abs(invariant2.real())) ||
      invariant1.real() <= 0 || invariant2.real() <= 0) return false;
  const double score = std::abs(std::log(invariant1.real() / 1225.0)) +
                       std::abs(std::log(invariant2.real() / 42875.0));
  hit = Hit{score, t, invariant1.real(), invariant2.real(), eigenvalues};
  return true;
}

int main(int argc, char** argv) {
  const std::uint64_t count = argc > 1 ? std::stoull(argv[1]) : 10000000;
  const int bound = argc > 2 ? std::stoi(argv[2]) : 100;
  std::mt19937_64 generator(25042026);
  std::uniform_int_distribution<int> coordinate(-bound, bound);
  std::priority_queue<Hit> best;
  for (std::uint64_t index = 0; index < count; ++index) {
    std::array<int, 7> trajectory;
    for (int& value : trajectory) value = coordinate(generator);
    Hit hit;
    if (!evaluate(trajectory, hit)) continue;
    if (best.size() < 100) best.push(hit);
    else if (hit.score < best.top().score) { best.pop(); best.push(hit); }
  }
  std::vector<Hit> output;
  while (!best.empty()) { output.push_back(best.top()); best.pop(); }
  std::sort(output.begin(), output.end());
  std::cout << std::setprecision(16);
  for (const Hit& hit : output) {
    std::cout << hit.score << " (";
    for (int index = 0; index < 7; ++index)
      std::cout << (index ? "," : "") << hit.trajectory[index];
    std::cout << ") " << hit.invariant1 << " " << hit.invariant2 << " eig";
    for (C value : hit.eigenvalues) std::cout << " " << value;
    std::cout << "\n";
  }
}
