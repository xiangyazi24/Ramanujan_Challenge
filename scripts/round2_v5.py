#!/usr/bin/env python3
"""Round 2, V5: Monte Carlo zero counts for coefficients of random squares.

For each requested odd prime p, let d=(p-1)/2 and sample

    A(T) = a_0 + a_1 T + ... + a_{d-1} T^(d-1) + T^d

with the a_i independent and uniform in F_p.  We count zero coefficients of
A(T)^2 in degrees 1,...,2d-1, excluding both endpoint coefficients.

The convolution is evaluated by batched double-precision FFT, then every
integer coefficient of every sample is independently certified by two exact
number-theoretic transforms and CRT.  The CRT modulus is above 1e18, whereas
the largest possible integer coefficient in the requested experiment is below
6.3e10.  The first sample for each p also receives an exact int64 schoolbook
convolution check.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np


DEFAULT_PRIMES = (101, 503, 1009, 5003)
DEFAULT_SAMPLES = 1_000
DEFAULT_SEED = 20_260_715
DEFAULT_OUTPUT = Path("/tmp/round2_v5.txt")
NTT_PRIMES = ((998_244_353, 3), (1_004_535_809, 3))
CRT_MODULUS = NTT_PRIMES[0][0] * NTT_PRIMES[1][0]

_BIT_REVERSE_CACHE: dict[int, np.ndarray] = {}
_TWIDDLE_CACHE: dict[tuple[int, int, bool], np.ndarray] = {}


@dataclass(frozen=True)
class SimulationResult:
    p: int
    d: int
    samples: int
    mean: float
    population_variance: float
    sample_variance: float
    minimum: int
    maximum: int
    histogram: dict[int, int]
    max_fft_rounding_error: float
    schoolbook_check: bool
    all_sample_ntt_crt_check: bool
    counts_sha256: str
    elapsed: float
    counts: np.ndarray


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % q for q in range(3, math.isqrt(n) + 1, 2))


def next_power_of_two(n: int) -> int:
    return 1 << (n - 1).bit_length()


def bit_reverse_indices(n: int) -> np.ndarray:
    cached = _BIT_REVERSE_CACHE.get(n)
    if cached is not None:
        return cached
    source = np.arange(n, dtype=np.uint64)
    reverse = np.zeros(n, dtype=np.uint64)
    working = source.copy()
    for _ in range(n.bit_length() - 1):
        reverse = (reverse << 1) | (working & 1)
        working >>= 1
    result = reverse.astype(np.int64)
    _BIT_REVERSE_CACHE[n] = result
    return result


def twiddle_powers(modulus: int, length: int, inverse: bool) -> np.ndarray:
    key = (modulus, length, inverse)
    cached = _TWIDDLE_CACHE.get(key)
    if cached is not None:
        return cached
    primitive_root = next(root for prime, root in NTT_PRIMES if prime == modulus)
    root = pow(primitive_root, (modulus - 1) // length, modulus)
    if inverse:
        root = pow(root, modulus - 2, modulus)
    half = length // 2
    powers = np.empty(half, dtype=np.int64)
    powers[0] = 1
    for j in range(1, half):
        powers[j] = powers[j - 1] * root % modulus
    _TWIDDLE_CACHE[key] = powers
    return powers


def batched_ntt(values: np.ndarray, modulus: int, inverse: bool) -> np.ndarray:
    """Exact radix-2 NTT on the final axis of a two-dimensional int64 array."""
    batch, n = values.shape
    if n & (n - 1):
        raise ValueError("NTT length must be a power of two")
    if (modulus - 1) % n:
        raise ValueError(f"NTT length {n} is unsupported by modulus {modulus}")

    transformed = (values[:, bit_reverse_indices(n)] % modulus).copy()
    length = 2
    while length <= n:
        half = length // 2
        blocks = transformed.reshape(batch, n // length, length)
        left = blocks[:, :, :half].copy()
        right = (
            blocks[:, :, half:]
            * twiddle_powers(modulus, length, inverse)[None, None, :]
        ) % modulus
        blocks[:, :, :half] = (left + right) % modulus
        blocks[:, :, half:] = (left - right) % modulus
        length *= 2

    if inverse:
        transformed = transformed * pow(n, modulus - 2, modulus) % modulus
    return transformed


def exact_square_convolutions(polynomials: np.ndarray, nfft: int) -> np.ndarray:
    """Recover all integer square-convolution coefficients exactly by NTT/CRT."""
    batch, polynomial_length = polynomials.shape
    padded = np.zeros((batch, nfft), dtype=np.int64)
    padded[:, :polynomial_length] = polynomials
    residues = []
    for modulus, _ in NTT_PRIMES:
        transform = batched_ntt(padded, modulus, inverse=False)
        squared_transform = transform * transform % modulus
        residues.append(batched_ntt(squared_transform, modulus, inverse=True))

    q1 = NTT_PRIMES[0][0]
    q2 = NTT_PRIMES[1][0]
    inverse_q1_mod_q2 = pow(q1, q2 - 2, q2)
    correction = (
        ((residues[1] - residues[0]) % q2) * inverse_q1_mod_q2 % q2
    )
    # CRT_MODULUS < 2^63, so every operation here is exact in signed int64.
    return residues[0] + q1 * correction


def exact_expected_zero_count(p: int) -> float:
    d = (p - 1) // 2
    r = d // 2
    return (p - 1) / p - p ** (-(r + 1))


def simulate_prime(
    p: int, samples: int, base_seed: int, batch_size: int
) -> SimulationResult:
    if not is_prime(p) or p == 2:
        raise ValueError(f"p must be an odd prime, got {p}")
    if samples < 2:
        raise ValueError("at least two samples are needed for sample variance")

    start_time = time.perf_counter()
    d = (p - 1) // 2
    convolution_length = 2 * d + 1  # This equals p for the requested model.
    nfft = next_power_of_two(convolution_length)
    coefficient_bound = (d + 1) * (p - 1) ** 2
    if coefficient_bound >= CRT_MODULUS:
        raise ValueError(
            f"two-prime CRT is too small at p={p}: bound={coefficient_bound}"
        )
    rng = np.random.Generator(
        np.random.PCG64(np.random.SeedSequence([base_seed, p]))
    )

    counts = np.empty(samples, dtype=np.int64)
    max_rounding_error = 0.0
    schoolbook_check = False
    all_sample_ntt_crt_check = True

    for first in range(0, samples, batch_size):
        size = min(batch_size, samples - first)
        coefficients = np.zeros((size, nfft), dtype=np.float64)
        coefficients[:, :d] = rng.integers(
            0, p, size=(size, d), dtype=np.int64
        )
        coefficients[:, d] = 1.0

        transform = np.fft.rfft(coefficients, axis=1)
        raw_convolution = np.fft.irfft(
            transform * transform, n=nfft, axis=1
        )[:, :convolution_length]
        rounded_convolution = np.rint(raw_convolution).astype(np.int64)
        max_rounding_error = max(
            max_rounding_error,
            float(np.max(np.abs(raw_convolution - rounded_convolution))),
        )
        if max_rounding_error >= 0.5:
            raise ArithmeticError(
                f"unsafe FFT rounding at p={p}: error={max_rounding_error}"
            )

        integer_polynomials = coefficients[:, : d + 1].astype(np.int64)
        exact_convolution = exact_square_convolutions(integer_polynomials, nfft)[
            :, :convolution_length
        ]
        batch_exact = bool(np.array_equal(exact_convolution, rounded_convolution))
        all_sample_ntt_crt_check = all_sample_ntt_crt_check and batch_exact
        if not batch_exact:
            sample, degree = np.argwhere(
                exact_convolution != rounded_convolution
            )[0]
            raise ArithmeticError(
                f"FFT/NTT-CRT mismatch at p={p}, sample={first + sample}, "
                f"degree={degree}"
            )

        # Degrees 0 and 2d are deliberately excluded.  Counts are made from
        # the exact NTT/CRT result, not from the floating-point result.
        interior_mod_p = exact_convolution[:, 1:-1] % p
        counts[first : first + size] = np.count_nonzero(
            interior_mod_p == 0, axis=1
        )

        if first == 0:
            first_polynomial = coefficients[0, : d + 1].astype(np.int64)
            exact = np.convolve(first_polynomial, first_polynomial)
            schoolbook_check = bool(np.array_equal(exact, exact_convolution[0]))
            if not schoolbook_check:
                mismatch = np.flatnonzero(exact != exact_convolution[0])[0]
                raise ArithmeticError(
                    f"exact convolution check failed at p={p}, degree={mismatch}"
                )

    elapsed = time.perf_counter() - start_time
    histogram = dict(sorted(Counter(map(int, counts)).items()))
    return SimulationResult(
        p=p,
        d=d,
        samples=samples,
        mean=float(np.mean(counts)),
        population_variance=float(np.var(counts, ddof=0)),
        sample_variance=float(np.var(counts, ddof=1)),
        minimum=int(np.min(counts)),
        maximum=int(np.max(counts)),
        histogram=histogram,
        max_fft_rounding_error=max_rounding_error,
        schoolbook_check=schoolbook_check,
        all_sample_ntt_crt_check=all_sample_ntt_crt_check,
        counts_sha256=hashlib.sha256(counts.tobytes()).hexdigest(),
        elapsed=elapsed,
        counts=counts,
    )


def sieve_primes(limit: int) -> list[int]:
    flags = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        flags[0] = 0
    if limit >= 1:
        flags[1] = 0
    for q in range(2, math.isqrt(limit) + 1):
        if flags[q]:
            flags[q * q : limit + 1 : q] = b"\x00" * (
                (limit - q * q) // q + 1
            )
    return [q for q in range(2, limit + 1) if flags[q]]


def apery_zero_count(p: int) -> int:
    """Count zero Apéry coefficients b_j mod p for 0 <= j < p."""
    inverses = [0] * p
    inverses[1] = 1
    for n in range(2, p):
        inverses[n] = (p - (p // n) * inverses[p % n] % p) % p

    b_previous = 1
    b_current = 5 % p
    zero_count = int(b_current == 0)
    for n in range(1, p - 1):
        n3 = n * n % p * n % p
        recurrence_coefficient = (((34 * n + 51) * n + 27) * n + 5) % p
        inverse = inverses[n + 1]
        inverse_cube = inverse * inverse % p * inverse % p
        b_next = (
            (recurrence_coefficient * b_current - n3 * b_previous)
            * inverse_cube
            % p
        )
        b_previous, b_current = b_current, b_next
        zero_count += int(b_current == 0)
    return zero_count


def empirical_stats(values: np.ndarray) -> dict[str, object]:
    return {
        "n": len(values),
        "mean": float(np.mean(values)),
        "population_variance": float(np.var(values, ddof=0)),
        "sample_variance": float(np.var(values, ddof=1)),
        "minimum": int(np.min(values)),
        "maximum": int(np.max(values)),
        "histogram": dict(sorted(Counter(map(int, values)).items())),
    }


def poisson_probability(k: int, lam: float) -> float:
    return math.exp(-lam) * lam**k / math.factorial(k)


def twice_poisson_half_probability(k: int) -> float:
    if k % 2:
        return 0.0
    return poisson_probability(k // 2, 0.5)


def total_variation_distance(
    histogram: dict[int, int], n: int, probability: Callable[[int], float]
) -> float:
    # k=40 leaves a negligible tail for all distributions used here.
    return 0.5 * sum(
        abs(histogram.get(k, 0) / n - probability(k)) for k in range(41)
    )


def format_histogram(histogram: dict[int, int]) -> str:
    return "{" + ", ".join(f"{k}: {v}" for k, v in histogram.items()) + "}"


def comparison_rows(
    random_histogram: dict[int, int],
    random_n: int,
    apery_histogram: dict[int, int],
    apery_n: int,
    maximum: int,
) -> Iterable[str]:
    yield "  z   random-square   Apery observed   Poisson(1)   2*Poisson(1/2)"
    for z in range(maximum + 1):
        yield (
            f"  {z:>1d}   {random_histogram.get(z, 0) / random_n:>13.6f}"
            f"   {apery_histogram.get(z, 0) / apery_n:>14.6f}"
            f"   {poisson_probability(z, 1.0):>10.6f}"
            f"   {twice_poisson_half_probability(z):>16.6f}"
        )


def build_report(
    results: list[SimulationResult],
    base_seed: int,
    batch_size: int,
    observed_p_max: int,
    overall_start_time: float,
) -> str:
    pooled_counts = np.concatenate([result.counts for result in results])
    pooled = empirical_stats(pooled_counts)

    observed_primes = [p for p in sieve_primes(observed_p_max) if p >= 5]
    observed_counts = np.array(
        [apery_zero_count(p) for p in observed_primes], dtype=np.int64
    )
    observed = empirical_stats(observed_counts)
    total_elapsed = time.perf_counter() - overall_start_time

    lines = [
        "ROUND 2 V5: RANDOM-SQUARE CONVOLUTION ZERO COUNT",
        "=" * 78,
        "",
        "MODEL AND COUNTING CONVENTION",
        "-" * 78,
        "For each odd prime p, d=(p-1)/2 and",
        "  A(T)=a_0+a_1*T+...+a_(d-1)*T^(d-1)+T^d,",
        "where a_0,...,a_(d-1) are iid uniform in F_p.  Thus A is monic",
        "of degree exactly d; its constant coefficient is not conditioned to be nonzero.",
        "For A(T)^2, zero coefficients are counted only in degrees 1,...,2d-1.",
        "The endpoint degrees 0 and 2d are excluded.  There are 2d-1=p-2",
        "tested coefficients, exactly the same number of non-endpoint positions as in H_p.",
        f"Samples per prime: {results[0].samples}; base seed: {base_seed};",
        "independent PCG64 stream for p is seeded by SeedSequence([base_seed,p]).",
        f"FFT batch size: {batch_size}; NumPy version: {np.__version__}.",
        "",
        "EXACT EXPECTATION",
        "-" * 78,
        "Write A(T)^2=sum c_k*T^k and r=floor(d/2).  For d <= k <= 2d-1,",
        "c_k is affine with nonzero coefficient 2 in a_(k-d), conditional on the",
        "higher a_i, so Pr(c_k=0)=1/p.  For low even k=2m<d, c_(2m) is a",
        "nondegenerate odd-dimensional quadratic form, hence also Pr(c_(2m)=0)=1/p.",
        "For low odd k=2m-1<d,",
        "  c_(2m-1)=2*(a_0*a_(2m-1)+...+a_(m-1)*a_m),",
        "a split form 2*x dot y on F_p^m x F_p^m.  Counting x=0 separately gives",
        "  Pr(c_(2m-1)=0)=1/p+(p-1)/p^(m+1).",
        "Therefore, summing over all p-2 interior coefficients,",
        "  E[Z_square]=(p-2)/p + sum_(m=1)^r (p-1)/p^(m+1)",
        "             =(p-1)/p - p^(-(r+1)).",
        "The naive independent-uniform value (p-2)/p misses the split-form correction.",
        "",
        "MONTE CARLO RESULTS",
        "-" * 78,
        "Population variance divides by N=1000; sample variance is the unbiased N-1 value.",
        "Mean SE=sqrt(sample variance/N).",
        "  p      d    exact E     MC mean     diff       mean SE   pop var   sample var  min max  sec",
    ]
    for result in results:
        exact_mean = exact_expected_zero_count(result.p)
        mean_se = math.sqrt(result.sample_variance / result.samples)
        lines.append(
            f"  {result.p:<5d}  {result.d:<4d}  {exact_mean:>10.6f}"
            f"  {result.mean:>10.6f}  {result.mean - exact_mean:>9.6f}"
            f"  {mean_se:>9.6f}  {result.population_variance:>8.5f}"
            f"  {result.sample_variance:>10.5f}"
            f"  {result.minimum:>3d} {result.maximum:>3d}  {result.elapsed:>5.3f}"
        )
        lines.append(f"    histogram: {format_histogram(result.histogram)}")
        lines.append(
            f"    FFT max |x-round(x)|={result.max_fft_rounding_error:.3e}; "
            f"all-sample NTT/CRT check={'PASS' if result.all_sample_ntt_crt_check else 'FAIL'}; "
            f"first-sample schoolbook check={'PASS' if result.schoolbook_check else 'FAIL'}; "
            f"counts sha256={result.counts_sha256}"
        )

    random_histogram = pooled["histogram"]
    observed_histogram = observed["histogram"]
    assert isinstance(random_histogram, dict)
    assert isinstance(observed_histogram, dict)
    lines.extend(
        [
            "",
            "POOLED RANDOM-SQUARE RESULTS",
            "-" * 78,
            "Equal-weight average of the four exact expectations: "
            f"{sum(exact_expected_zero_count(result.p) for result in results) / len(results):.6f}.",
            f"N={pooled['n']}; mean={pooled['mean']:.6f}; "
            f"population variance={pooled['population_variance']:.6f}; "
            f"sample variance={pooled['sample_variance']:.6f}; "
            f"range=[{pooled['minimum']},{pooled['maximum']}].",
            f"Histogram: {format_histogram(random_histogram)}",
            "TV distance to Poisson(1): "
            f"{total_variation_distance(random_histogram, int(pooled['n']), lambda k: poisson_probability(k, 1.0)):.6f}.",
            "",
            "SAME-CURRENCY APERY COMPARISON",
            "-" * 78,
            f"Recomputed from the Apéry recurrence for every prime 5 <= p <= {observed_p_max} "
            f"({observed['n']} primes), counting 0 <= j < p.  The two endpoints are",
            "nonzero, so this is again a count over exactly p-2 eligible positions.",
            f"Mean Z(p)={observed['mean']:.6f}; "
            f"population variance={observed['population_variance']:.6f}; "
            f"sample variance={observed['sample_variance']:.6f}; "
            f"range=[{observed['minimum']},{observed['maximum']}].",
            f"Histogram: {format_histogram(observed_histogram)}",
            "",
        ]
    )
    maximum = max(
        int(pooled["maximum"]), int(observed["maximum"]), 8
    )
    lines.extend(
        comparison_rows(
            random_histogram,
            int(pooled["n"]),
            observed_histogram,
            int(observed["n"]),
            maximum,
        )
    )
    lines.extend(
        [
            "",
            "Distributional distances:",
            "  random square vs Poisson(1): "
            f"TV={total_variation_distance(random_histogram, int(pooled['n']), lambda k: poisson_probability(k, 1.0)):.6f}",
            "  Apery vs Poisson(1): "
            f"TV={total_variation_distance(observed_histogram, int(observed['n']), lambda k: poisson_probability(k, 1.0)):.6f}",
            "  Apery vs 2*Poisson(1/2): "
            f"TV={total_variation_distance(observed_histogram, int(observed['n']), twice_poisson_half_probability):.6f}",
            "",
            "INTERPRETATION",
            "-" * 78,
            "The exact random-square expectation tends to 1, and the simulation agrees within",
            "its Monte Carlo standard errors.  Thus random squaring matches the first-order",
            "scale Mean Z(p) about 1.  But its zero count is close to Poisson(1), including",
            "many odd counts and variance about 1.  The Apéry counts are almost all even",
            "because b_j=b_(p-1-j), and instead track 2*Poisson(1/2), whose variance is 2.",
            "Thus squareness explains the O(1) mean heuristic but, by itself, does not",
            "explain the observed paired distribution; reversal symmetry supplies that feature.",
            "The experiment follows the requested pure-square degree d=(p-1)/2 model and",
            "does not insert the extra quadratic S_p occurring for inert residue classes.",
            "It is also deliberately a generic-monic model, not the full law of the actual A_p.",
            "Since H_p(0)=S_p(0)=1, an actual factor has A_p(0)^2=1.  Since H_p and",
            "S_p are reciprocal, A_p^2 is reciprocal and hence A_p^*=+/-A_p in odd",
            "characteristic.  Thus its coefficients have a reversal constraint absent here.",
            "If one only fixes a_0 to a nonzero value (without imposing reciprocity), every",
            "coefficient below the middle degree becomes conditionally affine-uniform; the",
            "middle coefficient remains a quadratic level set.  This changes the small finite-p",
            "correction but not the limiting order-one expectation.",
            "",
            "VALIDATION AND RUNTIME",
            "-" * 78,
            "All coefficients of all 4,000 samples were recomputed by exact NTT modulo",
            f"{NTT_PRIMES[0][0]} and {NTT_PRIMES[1][0]}, then recovered by CRT modulo",
            f"{CRT_MODULUS}.  Maximum possible raw coefficient at p=5003 is",
            f"{(2502 * 5002**2):,}, so CRT recovery is unique.  The NTT/CRT integers matched",
            "the rounded floating-point FFT at every coefficient.  The reported rounding",
            "distance is only a diagnostic; the NTT/CRT comparison is the exact certificate.",
            "In addition, every p passed an int64 schoolbook check of all coefficients in",
            "sample 1.  Zero counts in the report are computed from the NTT/CRT integers.",
            f"Total script computation time (simulation + comparison): {total_elapsed:.3f} seconds.",
        ]
    )
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primes", type=int, nargs="+", default=DEFAULT_PRIMES)
    parser.add_argument("--samples", type=int, default=DEFAULT_SAMPLES)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--observed-p-max", type=int, default=5_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    start = time.perf_counter()
    results = []
    for p in args.primes:
        print(f"Simulating p={p}...", flush=True)
        result = simulate_prime(p, args.samples, args.seed, args.batch_size)
        results.append(result)
        print(
            f"  mean={result.mean:.5f}, popvar={result.population_variance:.5f}, "
            f"range=[{result.minimum},{result.maximum}], {result.elapsed:.3f}s",
            flush=True,
        )
    # Include the fast exact Apéry comparison in the reported total runtime.
    report = build_report(
        results=results,
        base_seed=args.seed,
        batch_size=args.batch_size,
        observed_p_max=args.observed_p_max,
        overall_start_time=start,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(report, encoding="utf-8")
    print(f"Wrote {args.output} ({len(report):,} bytes)")


if __name__ == "__main__":
    main()
