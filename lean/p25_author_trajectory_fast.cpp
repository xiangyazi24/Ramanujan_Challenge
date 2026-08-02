#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

static double ipow(double value, int exponent) {
  if (exponent < 0) return 1.0 / ipow(value, -exponent);
  double answer = 1.0;
  while (exponent) {
    if (exponent & 1) answer *= value;
    value *= value;
    exponent >>= 1;
  }
  return answer;
}

static bool three_real_roots(double a, double b, double c,
                             std::array<double, 3>& roots) {
  const double p = b - a * a / 3.0;
  const double q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
  const double discriminant = q * q / 4.0 + p * p * p / 27.0;
  if (discriminant > 1e-11 * (1.0 + std::abs(q * q))) return false;
  if (std::abs(p) < 1e-14) {
    roots = {-a / 3.0, -a / 3.0, -a / 3.0};
    return true;
  }
  const double radius = 2.0 * std::sqrt(std::max(0.0, -p / 3.0));
  double cosine = -q / (2.0 * std::sqrt(std::max(0.0, -p * p * p / 27.0)));
  cosine = std::max(-1.0, std::min(1.0, cosine));
  const double angle = std::acos(cosine) / 3.0;
  constexpr double tau_over_three = 2.0943951023931954923;
  for (int i = 0; i < 3; ++i)
    roots[i] = radius * std::cos(angle - tau_over_three * i) - a / 3.0;
  return true;
}

static std::array<double, 5> elementary(const std::array<int, 4>& values) {
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
  std::array<double, 3> eigenvalues;
  bool operator<(const Hit& other) const { return score < other.score; }
};

int main(int argc, char** argv) {
  const int bound = argc > 1 ? std::stoi(argv[1]) : 12;
  const int workers = argc > 2 ? std::stoi(argv[2]) :
      std::max(1u, std::thread::hardware_concurrency());
  const double target1 = 1225.0, target2 = 42875.0;
  std::atomic<int> next_a0(-bound);
  std::atomic<unsigned long long> tested(0);
  std::mutex mutex;
  std::priority_queue<Hit> global_best;

  auto work = [&]() {
    std::priority_queue<Hit> best;
    while (true) {
      const int a0 = next_a0.fetch_add(1);
      if (a0 > bound) break;
      for (int a1 = a0; a1 <= bound; ++a1)
      for (int a2 = a1; a2 <= bound; ++a2)
      for (int a3 = a2; a3 <= bound; ++a3) {
        const std::array<int, 4> av{a0, a1, a2, a3};
        const auto pa = elementary(av);
        for (int b0 = -bound; b0 <= bound; ++b0)
        for (int b1 = b0; b1 <= bound; ++b1)
        for (int b2 = b1; b2 <= bound; ++b2) {
          if (std::find(av.begin(), av.end(), b0) != av.end() ||
              std::find(av.begin(), av.end(), b1) != av.end() ||
              std::find(av.begin(), av.end(), b2) != av.end()) continue;
          ++tested;
          const std::array<int, 4> bv{0, b0, b1, b2};
          const auto pb = elementary(bv);
          const double lead = pb[1] - pa[1];
          if (std::abs(lead) < 0.5) continue;
          std::array<double, 3> roots;
          if (!three_real_roots((pb[2] - pa[2]) / lead,
                                (pb[3] - pa[3]) / lead,
                                (pb[4] - pa[4]) / lead, roots)) continue;
          std::array<double, 3> eigenvalues;
          bool valid = true;
          for (int r = 0; r < 3; ++r) {
            double value = 1.0;
            for (int step : av)
              if (step) value *= ipow(1.0 + roots[r] / step, step);
            for (int step : std::array<int, 3>{b0, b1, b2})
              if (step) value *= ipow(1.0 + roots[r] / step, -step);
            eigenvalues[r] = value;
            valid &= std::isfinite(value);
          }
          if (!valid) continue;
          const double e1 = eigenvalues[0] + eigenvalues[1] + eigenvalues[2];
          const double e2 = eigenvalues[0] * eigenvalues[1] +
                            eigenvalues[0] * eigenvalues[2] +
                            eigenvalues[1] * eigenvalues[2];
          const double e3 = eigenvalues[0] * eigenvalues[1] * eigenvalues[2];
          if (!std::isfinite(e3) || std::abs(e3) < 1e-200) continue;
          const double invariant1 = e1 * e2 / e3;
          const double invariant2 = e1 * e1 * e1 / e3;
          if (!(invariant1 > 0 && invariant2 > 0) ||
              !std::isfinite(invariant1) || !std::isfinite(invariant2)) continue;
          const double score = std::abs(std::log(invariant1 / target1)) +
                               std::abs(std::log(invariant2 / target2));
          Hit hit{score, {a0, a1, a2, a3, b0, b1, b2},
                  invariant1, invariant2, eigenvalues};
          if (best.size() < 100) best.push(hit);
          else if (score < best.top().score) { best.pop(); best.push(hit); }
        }
      }
    }
    std::lock_guard<std::mutex> lock(mutex);
    while (!best.empty()) {
      Hit hit = best.top(); best.pop();
      if (global_best.size() < 100) global_best.push(hit);
      else if (hit.score < global_best.top().score) {
        global_best.pop(); global_best.push(hit);
      }
    }
  };

  std::vector<std::thread> threads;
  for (int i = 0; i < workers; ++i) threads.emplace_back(work);
  for (auto& thread : threads) thread.join();
  std::vector<Hit> output;
  while (!global_best.empty()) { output.push_back(global_best.top()); global_best.pop(); }
  std::sort(output.begin(), output.end());
  std::cout << "tested " << tested << "\n" << std::setprecision(16);
  for (const Hit& hit : output) {
    std::cout << hit.score << " (";
    for (int i = 0; i < 7; ++i) std::cout << (i ? "," : "") << hit.trajectory[i];
    std::cout << ") " << hit.invariant1 << " " << hit.invariant2 << " eig";
    for (double value : hit.eigenvalues) std::cout << " " << value;
    std::cout << "\n";
  }
}
