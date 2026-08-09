#!/usr/bin/env python3
"""Exact/rational and Fourier checks for Codex Task 002(A).

The exact count and error are computed with integer CRT and fractions.  The
Fourier reconstruction is independent numerical DFT evaluation of

    1/(p q ell) sum_{a,b,c} F_p(a) F_q(b) F_ell(c)
        sum_{0 <= m < L} exp(-2 pi i m(a/p+b/q+c/ell)).

For the orthogonality check, ``v^{-1}`` necessarily means that v ranges over
the units modulo p p'.  The exact check expands the Fourier transforms and
uses the prime Ramanujan sum; floating point is used only for an additional
direct spot check.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent


def primes_upto(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for prime in range(2, math.isqrt(limit) + 1):
        if sieve[prime]:
            start = prime * prime
            sieve[start : limit + 1 : prime] = b"\x00" * (
                (limit - start) // prime + 1
            )
    return [n for n in range(2, limit + 1) if sieve[n]]


def apery_zeros(prime: int) -> tuple[int, ...]:
    """Return {0 <= r < p : b_r = 0 mod p} without modular divisions."""

    previous, current = 1 % prime, 5 % prime
    zeros: list[int] = []
    if previous == 0:
        zeros.append(0)
    if current == 0:
        zeros.append(1)
    for index in range(1, prime - 1):
        polynomial = (
            34 * index**3 + 51 * index**2 + 27 * index + 5
        ) % prime
        following = (
            polynomial * current - pow(index, 6, prime) * previous
        ) % prime
        previous, current = current, following
        if current == 0:
            zeros.append(index + 1)
    return tuple(zeros)


def crt3(
    rp: int, p: int, rq: int, q: int, rs: int, ell: int
) -> int:
    """Least nonnegative CRT solution for three pairwise-coprime moduli."""

    modulus = p * q * ell
    result = 0
    for residue, prime in ((rp, p), (rq, q), (rs, ell)):
        cofactor = modulus // prime
        result += residue * cofactor * pow(cofactor, -1, prime)
    return result % modulus


def exact_triple_count(
    p: int,
    q: int,
    ell: int,
    zero_sets: dict[int, tuple[int, ...]],
    length: int,
) -> int:
    count = 0
    for rp, rq, rs in itertools.product(
        zero_sets[p], zero_sets[q], zero_sets[ell]
    ):
        count += crt3(rp, p, rq, q, rs, ell) < length
    return count


def positive_fourier(prime: int, zeros: Sequence[int]) -> np.ndarray:
    """F[a] = sum_{r in Z} exp(+2 pi i a r/p)."""

    indicator = np.zeros(prime, dtype=np.complex128)
    if zeros:
        indicator[np.asarray(zeros, dtype=np.int64)] = 1.0
    return np.conjugate(np.fft.fft(indicator))


def fourier_reconstruction(
    p: int,
    q: int,
    ell: int,
    zero_sets: dict[int, tuple[int, ...]],
    length: int,
) -> tuple[complex, float, float]:
    """Evaluate the full trilinear Fourier formula by one length-pqell FFT."""

    modulus = p * q * ell
    fp = positive_fourier(p, zero_sets[p])
    fq = positive_fourier(q, zero_sets[q])
    fell = positive_fourier(ell, zero_sets[ell])

    frequencies = np.arange(modulus, dtype=np.int64)
    a = (frequencies * pow(q * ell, -1, p)) % p
    b = (frequencies * pow(p * ell, -1, q)) % q
    c = (frequencies * pow(p * q, -1, ell)) % ell
    coefficient = fp[a] * fq[b] * fell[c]

    # Fourier inversion has the negative exponential, hence numpy.fft.fft.
    reconstructed = np.fft.fft(coefficient) / modulus
    shell_count = np.sum(reconstructed[:length])

    exact_indicator = np.fromiter(
        (
            1.0
            if (
                (m % p in zero_sets[p])
                and (m % q in zero_sets[q])
                and (m % ell in zero_sets[ell])
            )
            else 0.0
            for m in range(length)
        ),
        dtype=np.float64,
        count=length,
    )
    prefix_error = float(
        np.max(np.abs(reconstructed[:length].real - exact_indicator))
    )
    imaginary_error = float(np.max(np.abs(reconstructed[:length].imag)))
    return shell_count, prefix_error, imaginary_error


def exact_unit_sum(prime: int, zeros: Sequence[int], k: int) -> int:
    """Exact sum_{u != 0} F_p(k u), expanded via Ramanujan sums."""

    if not 0 < k < prime:
        raise ValueError("the claimed identity requires k nonzero modulo p")
    # sum_{u != 0} e_p(u k r) is p-1 if kr=0 and -1 otherwise.
    return sum(prime - 1 if (k * r) % prime == 0 else -1 for r in zeros)


def direct_unit_orthogonality_spot_check(
    p: int,
    q: int,
    zp: Sequence[int],
    zq: Sequence[int],
    k: int,
    kp: int,
) -> complex:
    fp = positive_fourier(p, zp)
    fq = positive_fourier(q, zq)
    total = 0.0j
    for v in range(p * q):
        if v % p == 0 or v % q == 0:
            continue
        ap = (k * pow(v, -1, p)) % p
        aq = (kp * pow(v, -1, q)) % q
        total += fp[ap] * np.conjugate(fq[aq])
    return total


def verify_orthogonality(
    primes: Sequence[int], zero_sets: dict[int, tuple[int, ...]]
) -> tuple[int, float]:
    exact_checks = 0
    max_numeric_error = 0.0
    for p, q in itertools.permutations(primes, 2):
        zp, zq = zero_sets[p], zero_sets[q]
        if 0 in zp or 0 in zq:
            raise AssertionError("b_0=1, so zero must never belong to Z_p")

        # This checks every nonzero k and k' exactly.  CRT makes the double
        # unit sum the product of the two displayed one-prime sums.
        left_factors = [exact_unit_sum(p, zp, k) for k in range(1, p)]
        right_factors = [exact_unit_sum(q, zq, k) for k in range(1, q)]
        if any(value != -len(zp) for value in left_factors):
            raise AssertionError((p, left_factors))
        if any(value != -len(zq) for value in right_factors):
            raise AssertionError((q, right_factors))
        target = len(zp) * len(zq)
        for left in left_factors:
            for right in right_factors:
                if left * right != target:
                    raise AssertionError((p, q, left, right, target))
                exact_checks += 1

        # Directly enumerate v for three representative frequency pairs.
        samples = ((1, 1), (min(2, p - 1), min(3, q - 1)), (p - 1, q - 1))
        for k, kp in samples:
            direct = direct_unit_orthogonality_spot_check(
                p, q, zp, zq, k, kp
            )
            max_numeric_error = max(max_numeric_error, abs(direct - target))
    return exact_checks, max_numeric_error


def write_csv(path: Path, rows: Iterable[dict[str, object]]) -> None:
    rows = list(rows)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--x", type=int, default=50)
    parser.add_argument(
        "--csv", type=Path, default=HERE / "task002_trilinear_X50.csv"
    )
    args = parser.parse_args()

    x = args.x
    length = x * x
    primes = [p for p in primes_upto(2 * x) if p > x]
    zero_sets = {p: apery_zeros(p) for p in primes}
    rows: list[dict[str, object]] = []
    max_count_error = 0.0
    max_prefix_error = 0.0
    max_imaginary_error = 0.0

    for p, q, ell in itertools.combinations(primes, 3):
        modulus = p * q * ell
        zproduct = len(zero_sets[p]) * len(zero_sets[q]) * len(zero_sets[ell])
        count = exact_triple_count(p, q, ell, zero_sets, length)
        main_term = Fraction(length * zproduct, modulus)
        exact_error = Fraction(count, 1) - main_term

        if zproduct:
            fourier_count, prefix_error, imaginary_error = fourier_reconstruction(
                p, q, ell, zero_sets, length
            )
        else:
            fourier_count, prefix_error, imaginary_error = 0.0j, 0.0, 0.0
        fourier_error = fourier_count - float(main_term)
        count_error = abs(fourier_count - count)
        max_count_error = max(max_count_error, count_error)
        max_prefix_error = max(max_prefix_error, prefix_error)
        max_imaginary_error = max(max_imaginary_error, imaginary_error)
        if count_error > 2e-8 or prefix_error > 2e-8 or imaginary_error > 2e-8:
            raise AssertionError(
                (p, q, ell, count_error, prefix_error, imaginary_error)
            )

        rows.append(
            {
                "p": p,
                "q": q,
                "ell": ell,
                "Zp": len(zero_sets[p]),
                "Zq": len(zero_sets[q]),
                "Zell": len(zero_sets[ell]),
                "count": count,
                "main_numerator": main_term.numerator,
                "main_denominator": main_term.denominator,
                "error_numerator": exact_error.numerator,
                "error_denominator": exact_error.denominator,
                "fourier_error_real": f"{fourier_error.real:.17g}",
                "fourier_error_imag": f"{fourier_error.imag:.3e}",
                "count_abs_error": f"{count_error:.3e}",
            }
        )

    exact_checks, orthogonality_error = verify_orthogonality(primes, zero_sets)
    write_csv(args.csv, rows)
    print(f"X={x}; primes={primes}")
    print("zero sets:")
    for p in primes:
        print(f"  p={p}: {zero_sets[p]}")
    print(f"triples checked: {len(rows)}")
    print(f"nonzero triples: {sum(int(row['Zp'])*int(row['Zq'])*int(row['Zell']) > 0 for row in rows)}")
    print(f"max |Fourier count - exact count|: {max_count_error:.3e}")
    print(f"max reconstructed-indicator error: {max_prefix_error:.3e}")
    print(f"max imaginary reconstruction error: {max_imaginary_error:.3e}")
    print(f"exact orthogonality (p,q,k,k') checks: {exact_checks}")
    print(f"max direct-v orthogonality error: {orthogonality_error:.3e}")
    print(f"CSV: {args.csv}")


if __name__ == "__main__":
    main()
