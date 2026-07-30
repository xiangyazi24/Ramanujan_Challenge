#!/usr/bin/env python3
"""Q5707: dependency-free exact shell/carrier audit for two Apéry models.

Models
------
Lambda = (1+x)(1+y)(1+z)((1+y)(1+z)+xyz)/(xyz)
Phi    = (x+y)(1+z)(x+y+z)(1+y+z)/(xyz)

For F^M = sum_v c_M(v) X^v, evaluate

    C_M^F(d) = sum_{v in d Z^3} c_M(v)

by residue-binomial sums, then form the primitive Newton carrier

    G_F(d0,L) = sum_{i=0}^L omega_i C_M^F(d0+i),
    omega_i = (-1)^i binom(d0+i,i) binom(d0+L+1,L-i).

Only the Python standard library is used.  All integer arithmetic is exact.
"""
from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from functools import lru_cache, reduce
from math import comb, gcd, isqrt, log, log10, prod
from typing import Dict, Iterable, List, Tuple

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

Monomial = Tuple[int, int, int]
Polynomial = Dict[Monomial, int]


@lru_cache(maxsize=None)
def B(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


def prime_sieve(n: int) -> List[int]:
    mark = bytearray(b"\x01") * (n + 1)
    if n >= 0:
        mark[0] = 0
    if n >= 1:
        mark[1] = 0
    for p in range(2, isqrt(n) + 1):
        if mark[p]:
            mark[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [p for p in range(2, n + 1) if mark[p]]


def apery_numbers(n: int) -> List[int]:
    if n == 0:
        return [1]
    out = [1, 5]
    for m in range(1, n):
        numerator = (34 * m**3 + 51 * m**2 + 27 * m + 5) * out[m] - m**3 * out[m - 1]
        denominator = (m + 1) ** 3
        q, r = divmod(numerator, denominator)
        assert r == 0
        out.append(q)
    return out


def valuation(n: int, p: int) -> int:
    if n == 0:
        raise ValueError("valuation of zero is infinite")
    n = abs(n)
    value = 0
    while n % p == 0:
        n //= p
        value += 1
    return value


def height(n: int) -> Dict[str, float | int]:
    n = abs(n)
    if n == 0:
        return {"bits": 0, "digits": 1, "log10": float("-inf"), "ln": float("-inf")}
    text = str(n)
    take = min(17, len(text))
    mantissa = int(text[:take]) / 10 ** (take - 1)
    lg10 = len(text) - 1 + log10(mantissa)
    return {"bits": n.bit_length(), "digits": len(text), "log10": lg10, "ln": lg10 * log(10.0)}


def residue_binomial_sum(n: int, residue: int, modulus: int) -> int:
    """sum_{0<=k<=n, k=residue (mod modulus)} binom(n,k)."""
    residue %= modulus
    return sum(B(n, k) for k in range(residue, n + 1, modulus))


def shell_lambda(M: int, d: int) -> int:
    """Exact C_M^Lambda(d).

    Symmetry k -> (2M-t)-k changes the natural residue M-t to M, giving

      sum_t binom(M,t) R(M;t,d) R(2M-t;M,d)^2.
    """
    return sum(
        B(M, t)
        * residue_binomial_sum(M, t, d)
        * residue_binomial_sum(2 * M - t, M, d) ** 2
        for t in range(M + 1)
    )


def shell_phi(M: int, d: int) -> int:
    """Exact C_M^Phi(d).

    Let R be the non-x degree in (x+y+z)^M, i the x-degree in (x+y)^M,
    beta the y-degree in the R-block, delta the y-degree in (1+y+z)^M,
    and k the combined z-degree from (1+z)^M and the remaining fourth
    factor.  The shell congruences are

      i = R                 (mod d),
      beta + delta = i      (mod d),
      k = M - R + beta      (mod d).
    """
    binM = [B(M, k) for k in range(M + 1)]
    z_cache: Dict[Tuple[int, int], int] = {}

    def z_sum(delta: int, residue: int) -> int:
        key = (delta, residue % d)
        if key not in z_cache:
            z_cache[key] = residue_binomial_sum(2 * M - delta, key[1], d)
        return z_cache[key]

    total = 0
    for R in range(M + 1):
        for i in range(R % d, M + 1, d):
            prefix = binM[R] * binM[i]
            for beta in range(R + 1):
                inner = 0
                delta0 = (i - beta) % d
                z_residue = (M - R + beta) % d
                for delta in range(delta0, M + 1, d):
                    inner += binM[delta] * z_sum(delta, z_residue)
                total += prefix * B(R, beta) * inner
    return total


def weights(d0: int, L: int) -> List[int]:
    result = [(-1) ** i * B(d0 + i, i) * B(d0 + L + 1, L - i) for i in range(L + 1)]
    assert sum(result) == 1
    assert reduce(gcd, (abs(w) for w in result), 0) == 1
    return result


def carrier(shells: List[int], shell_start: int, d0: int, L: int) -> int:
    ws = weights(d0, L)
    return sum(ws[i] * shells[d0 + i - shell_start] for i in range(L + 1))


def trial_factor(n: int, limit: int = 100_000) -> Tuple[Dict[int, int], int]:
    """Remove every prime factor <= limit; return exact factors and cofactor."""
    n = abs(n)
    factors: Dict[int, int] = {}
    for p in prime_sieve(limit):
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            factors[p] = exponent
    return factors, n


# ---------------------------------------------------------------------------
# Independent polynomial and mutation self-tests
# ---------------------------------------------------------------------------

def poly_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    out: Polynomial = {}
    for u, a in left.items():
        for v, b in right.items():
            exponent = (u[0] + v[0], u[1] + v[1], u[2] + v[2])
            out[exponent] = out.get(exponent, 0) + a * b
    return out


def poly_pow(base: Polynomial, exponent: int) -> Polynomial:
    out: Polynomial = {(0, 0, 0): 1}
    for _ in range(exponent):
        out = poly_mul(out, base)
    return out


def shifted_product(factors: Iterable[Polynomial], shift: Monomial) -> Polynomial:
    out: Polynomial = {(0, 0, 0): 1}
    for factor in factors:
        out = poly_mul(out, factor)
    return {(u[0] + shift[0], u[1] + shift[1], u[2] + shift[2]): c for u, c in out.items()}


def base_models() -> Tuple[Polynomial, Polynomial]:
    one_x = {(0, 0, 0): 1, (1, 0, 0): 1}
    one_y = {(0, 0, 0): 1, (0, 1, 0): 1}
    one_z = {(0, 0, 0): 1, (0, 0, 1): 1}
    inner = {(0, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1, (0, 1, 1): 1, (1, 1, 1): 1}
    lam = shifted_product((one_x, one_y, one_z, inner), (-1, -1, -1))

    x_plus_y = {(1, 0, 0): 1, (0, 1, 0): 1}
    x_plus_y_plus_z = {(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1}
    one_plus_y_plus_z = {(0, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1}
    phi = shifted_product((x_plus_y, one_z, x_plus_y_plus_z, one_plus_y_plus_z), (-1, -1, -1))
    return lam, phi


def shell_from_expansion(base: Polynomial, M: int, d: int) -> int:
    return sum(c for exponent, c in poly_pow(base, M).items() if all(e % d == 0 for e in exponent))


def eval_lambda(X: Fraction, Y: Fraction, Z: Fraction) -> Fraction:
    return (1 + X) * (1 + Y) * (1 + Z) * ((1 + Y) * (1 + Z) + X * Y * Z) / (X * Y * Z)


def eval_phi(x: Fraction, y: Fraction, z: Fraction) -> Fraction:
    return (x + y) * (1 + z) * (x + y + z) * (1 + y + z) / (x * y * z)


def self_test() -> None:
    lam, phi = base_models()
    for M in range(5):
        for d in range(1, 6):
            assert shell_lambda(M, d) == shell_from_expansion(lam, M, d)
            assert shell_phi(M, d) == shell_from_expansion(phi, M, d)

    for x, y, z in (
        (Fraction(2), Fraction(3), Fraction(5)),
        (Fraction(3, 2), Fraction(5, 3), Fraction(7, 2)),
        (Fraction(-2), Fraction(1, 3), Fraction(4)),
    ):
        X, Y, Z = (x + y) / z, 1 / (y + z), y / x
        assert eval_lambda(X, Y, Z) == eval_phi(x, y, z)


# ---------------------------------------------------------------------------
# Row and adjacent-family audits
# ---------------------------------------------------------------------------

def row_audit(n: int, L: int = 50) -> Dict[str, object]:
    M = n - 1
    d0 = M // 2 + 1
    ds = list(range(d0, d0 + L + 1))
    CL = [shell_lambda(M, d) for d in ds]
    CP = [shell_phi(M, d) for d in ds]
    GL = carrier(CL, d0, d0, L)
    GP = carrier(CP, d0, d0, L)
    D = GL - GP

    b = apery_numbers(M)
    prime_set = set(prime_sieve(n))
    candidates = [d + 1 for d in ds if d + 1 in prime_set]
    targets = [p for p in candidates if b[n - p] % p == 0]
    P = prod(candidates)
    R = prod(targets)
    support_GL = [p for p in candidates if GL % p == 0]
    support_GP = [p for p in candidates if GP % p == 0]
    support_D = [p for p in candidates if D % p == 0]

    # Both models have the same Cartier residue b_{n-p}; hence D contains
    # every candidate prime, not only the target primes.
    assert set(targets) <= set(support_GL)
    assert set(targets) <= set(support_GP)
    assert support_D == candidates

    return {
        "n": n,
        "M": M,
        "d0": d0,
        "L": L,
        "candidates": candidates,
        "targets": targets,
        "P": P,
        "R": R,
        "support_GL": support_GL,
        "support_GP": support_GP,
        "support_D": support_D,
        "G_lambda": GL,
        "G_phi": GP,
        "D": D,
        "gcd_models": gcd(abs(GL), abs(GP)),
        "gcd_GL_P": gcd(abs(GL), P),
        "gcd_GP_P": gcd(abs(GP), P),
        "gcd_D_P": gcd(abs(D), P),
    }


def emit_row(result: Dict[str, object]) -> None:
    print("=" * 100)
    print(f"ROW n={result['n']} M={result['M']} d0={result['d0']} L={result['L']}")
    for key in ("candidates", "targets", "P", "R", "support_GL", "support_GP", "support_D"):
        print(f"{key}={result[key]}")
    for key in ("G_lambda", "G_phi", "D", "gcd_models"):
        value = int(result[key])
        print(f"{key}={value}")
        print(f"{key}_height={height(value)}")
    factors, cofactor = trial_factor(int(result["D"]))
    print(f"D_factors_le_100000={factors}")
    print(f"D_cofactor={cofactor}")
    print(f"gcd_GL_P={result['gcd_GL_P']}")
    print(f"gcd_GP_P={result['gcd_GP_P']}")
    print(f"gcd_D_P={result['gcd_D_P']}")
    for p in (179, 193, 211):
        print(
            f"valuation p={p}: "
            f"GL={valuation(int(result['G_lambda']), p)} "
            f"GP={valuation(int(result['G_phi']), p)} "
            f"D={valuation(int(result['D']), p)} "
            f"gcd={valuation(int(result['gcd_models']), p)}"
        )


def family_audit(n: int) -> Dict[str, object]:
    """Three cover-preserving adjacent stencils.

    (d0,50), (d0,51), and (d0-1,51) all contain the original interval
    d0,...,d0+50, so every original target divides all six model carriers.
    """
    M = n - 1
    d0 = M // 2 + 1
    dmin, dmax = d0 - 1, d0 + 51
    CL = [shell_lambda(M, d) for d in range(dmin, dmax + 1)]
    CP = [shell_phi(M, d) for d in range(dmin, dmax + 1)]
    specs = ((d0, 50), (d0, 51), (d0 - 1, 51))
    GL = [carrier(CL, dmin, start, L) for start, L in specs]
    GP = [carrier(CP, dmin, start, L) for start, L in specs]
    DD = [a - b for a, b in zip(GL, GP)]
    return {
        "n": n,
        "specs": specs,
        "G_lambda": GL,
        "G_phi": GP,
        "D": DD,
        "gcd_lambda": reduce(gcd, (abs(x) for x in GL)),
        "gcd_phi": reduce(gcd, (abs(x) for x in GP)),
        "gcd_all_six": reduce(gcd, (abs(x) for x in GL + GP)),
        "gcd_D": reduce(gcd, (abs(x) for x in DD)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=int, nargs="*", default=[321])
    parser.add_argument("--family", action="store_true")
    parser.add_argument("--L", type=int, default=50)
    args = parser.parse_args()

    self_test()
    print("SELF_TEST=PASS")
    for n in args.rows:
        emit_row(row_audit(n, args.L))
        if args.family:
            print(f"family={family_audit(n)}")


if __name__ == "__main__":
    main()
