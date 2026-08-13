#!/usr/bin/env python3
"""Verifier for the reflected-polynomial no-go in Q8042.

This script uses only the Python standard library.  It does *not* model the
actual Apery shell coefficients.  It verifies two adversarial polynomial
families showing that degree/content/height/reflection/exact-window-support
alone cannot imply seven-prime nonalignment.

Default demonstration:
    python3 drops/Q8042-b50f802a-verify.py --X 40 --variant both

X must be even and at least 8.  The `monic` construction performs coefficient
CRT and is intended for moderate demonstration sizes (roughly X <= 80).
"""

from __future__ import annotations

import argparse
import math
from functools import reduce
from typing import Dict, Iterable, List, Sequence, Tuple


Poly = List[int]  # ascending coefficients


def sieve_primes(n: int) -> List[int]:
    if n < 2:
        return []
    mark = bytearray(b"\x01") * (n + 1)
    mark[0:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if mark[p]:
            mark[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [i for i in range(2, n + 1) if mark[i]]


def trim(a: Poly) -> Poly:
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def poly_add_int(a: Sequence[int], b: Sequence[int]) -> Poly:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
    return trim(out)


def poly_mul_int(a: Sequence[int], b: Sequence[int]) -> Poly:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def poly_mul_mod(a: Sequence[int], b: Sequence[int], p: int) -> Poly:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % p
    return trim(out)


def poly_eval(a: Sequence[int], x: int) -> int:
    y = 0
    for c in reversed(a):
        y = y * x + c
    return y


def poly_eval_mod(a: Sequence[int], x: int, p: int) -> int:
    y = 0
    for c in reversed(a):
        y = (y * x + c) % p
    return y


def poly_divmod_mod(a: Sequence[int], b: Sequence[int], p: int) -> Tuple[Poly, Poly]:
    aa = trim([x % p for x in a])
    bb = trim([x % p for x in b])
    if bb == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    if len(aa) < len(bb):
        return [0], aa
    q = [0] * (len(aa) - len(bb) + 1)
    inv = pow(bb[-1], -1, p)
    while len(aa) >= len(bb) and aa != [0]:
        k = len(aa) - len(bb)
        c = aa[-1] * inv % p
        q[k] = c
        for j in range(len(bb)):
            aa[k + j] = (aa[k + j] - c * bb[j]) % p
        trim(aa)
    return trim(q), trim(aa)


def poly_gcd_mod(a: Sequence[int], b: Sequence[int], p: int) -> Poly:
    aa = trim([x % p for x in a])
    bb = trim([x % p for x in b])
    while bb != [0]:
        _, rr = poly_divmod_mod(aa, bb, p)
        aa, bb = bb, rr
    if aa == [0]:
        return [0]
    inv = pow(aa[-1], -1, p)
    return [(x * inv) % p for x in aa]


def derivative(a: Sequence[int]) -> Poly:
    if len(a) <= 1:
        return [0]
    return trim([i * a[i] for i in range(1, len(a))])


def u_to_t(c_u: Sequence[int]) -> Poly:
    """Expand sum_j c_u[j] * (T(T+1))^j in the T basis."""
    d = len(c_u) - 1
    out = [0] * (2 * d + 1)
    for j, c in enumerate(c_u):
        if c == 0:
            continue
        # U^j = T^j (T+1)^j = sum_{k=0}^j binom(j,k) T^(j+k).
        for k in range(j + 1):
            out[j + k] += c * math.comb(j, k)
    return trim(out)


def balanced_crt(residues: Sequence[int], primes: Sequence[int], modulus: int) -> int:
    x = 0
    for a, p in zip(residues, primes):
        m = modulus // p
        x = (x + (a % p) * m * pow(m, -1, p)) % modulus
    return x if x <= modulus // 2 else x - modulus


def product(values: Iterable[int]) -> int:
    return reduce(lambda x, y: x * y, values, 1)


def falling(n: int, k: int) -> int:
    out = 1
    for j in range(k):
        out *= max(0, n - j)
    return out


def first_root_free_cubic(p: int) -> Poly:
    # x^3 - x is not injective (0 and 1 have the same image), hence it is
    # not surjective on F_p.  Choose a missing value c.
    image = {(x * x * x - x) % p for x in range(p)}
    for c in range(p):
        if c not in image:
            return [(-c) % p, (-1) % p, 0, 1]
    raise AssertionError("no missing value for x^3-x")


def root_free_squarefree_monic(p: int, degree: int) -> Poly:
    """A squarefree monic polynomial of given degree >= 2 with no F_p root."""
    if degree < 2:
        raise ValueError("degree must be at least 2")
    out: Poly = [1]
    remaining = degree
    if remaining % 2 == 1:
        out = poly_mul_mod(out, first_root_free_cubic(p), p)
        remaining -= 3
    need = remaining // 2
    nonresidues = [a for a in range(1, p) if pow(a, (p - 1) // 2, p) == p - 1]
    if len(nonresidues) < need:
        raise AssertionError("not enough distinct quadratic nonresidues")
    for a in nonresidues[:need]:
        out = poly_mul_mod(out, [(-a) % p, 0, 1], p)
    if len(out) - 1 != degree:
        raise AssertionError("cofactor degree mismatch")
    if any(poly_eval_mod(out, x, p) == 0 for x in range(p)):
        raise AssertionError("cofactor unexpectedly has an F_p root")
    if len(poly_gcd_mod(out, derivative(out), p)) != 1:
        raise AssertionError("cofactor unexpectedly is not squarefree")
    return out


def common_setup(X: int) -> Tuple[int, int, int, List[int], int]:
    if X < 8 or X % 2:
        raise ValueError("X must be even and at least 8")
    M = X * X
    m0 = 3 * M // 4
    u0 = m0 * (m0 + 1)
    primes = [p for p in sieve_primes(2 * X) if X < p <= 2 * X and (2 * m0 + 1) % p]
    if not primes:
        raise ValueError("empty good-prime window")
    P = product(primes)
    return M, m0, u0, primes, P


def rank3_polynomial_u(X: int, P: int, u0: int) -> Poly:
    D = X // 2
    q = [0] * (D + 1)
    q[0] = P - u0
    q[1] = 1
    q[D] += P
    return trim(q)


def monic_polynomial_u(X: int, primes: Sequence[int], P: int, u0: int) -> Tuple[Poly, Poly]:
    D = X // 2
    n = D - 1
    locals_by_p: Dict[int, Poly] = {p: root_free_squarefree_monic(p, n) for p in primes}
    H = []
    for j in range(n):
        H.append(balanced_crt([locals_by_p[p][j] for p in primes], primes, P))
    H.append(1)
    for p in primes:
        if [x % p for x in H] != locals_by_p[p]:
            raise AssertionError("coefficient CRT reconstruction failed")
    anchor = [P - u0, 1]  # U-u0+P
    Q = poly_mul_int(anchor, H)
    return Q, H


def verify_exact_support(
    name: str,
    X: int,
    Q_u: Sequence[int],
    primes: Sequence[int],
    P: int,
    M: int,
    m0: int,
    full_window: bool,
) -> None:
    Q_t = u_to_t(Q_u)
    if len(Q_t) - 1 != X:
        raise AssertionError(f"{name}: degree is {len(Q_t)-1}, expected {X}")
    content = reduce(math.gcd, (abs(c) for c in Q_t), 0)
    if content != 1:
        raise AssertionError(f"{name}: coefficient content is {content}, expected 1")

    # Reflection is exact because U(T)=T(T+1)=U(-T-1).  Check many integers too.
    for t in range(-2 * X, 2 * X + 1):
        u = t * (t + 1)
        ur = (-t - 1) * (-t)
        if poly_eval(Q_u, u) != poly_eval(Q_u, ur):
            raise AssertionError(f"{name}: reflection failed at T={t}")

    for p in primes:
        roots = {r for r in range(p) if poly_eval_mod(Q_t, r, p) == 0}
        expected = {m0 % p, (-m0 - 1) % p}
        if roots != expected:
            raise AssertionError(f"{name}: roots mod {p}: {sorted(roots)} != {sorted(expected)}")

    values = range(M) if full_window else list(range(min(M, 4 * X))) + [m0]
    for m in values:
        u = m * (m + 1)
        value = poly_eval(Q_u, u)
        expected_support = [p for p in primes if m % p in {m0 % p, (-m0 - 1) % p}]
        actual_support = [p for p in primes if value % p == 0]
        if actual_support != expected_support:
            raise AssertionError(f"{name}: support mismatch at m={m}")
        expected_gcd = product(expected_support)
        if math.gcd(abs(value), P) != expected_gcd:
            raise AssertionError(f"{name}: gcd mismatch at m={m}")

    anchor_value = poly_eval(Q_u, m0 * (m0 + 1))
    if anchor_value == 0 or math.gcd(abs(anchor_value), P) != P:
        raise AssertionError(f"{name}: common anchor failed")

    height = max(abs(c) for c in Q_t)
    lam = 2.0 * sum(1.0 / p for p in primes)
    star7 = falling(len(primes), 7)
    target_proxy = M * (lam**7)
    ratio = float("inf") if target_proxy == 0 else star7 / target_proxy
    print(f"[{name}] degree={len(Q_t)-1}, content={content}, leading={Q_t[-1]}")
    print(f"[{name}] good primes={len(primes)}, P bits={P.bit_length()}, log(height)={math.log(height):.6f}")
    print(f"[{name}] log(height)/(X log X)={math.log(height)/(X*math.log(X)):.6f}")
    print(f"[{name}] anchor K={len(primes)}, lambda={lam:.8f}, (K)_7/(M lambda^7)={ratio:.8e}")


def verify_monic_squarefree(X: int, Q_u: Sequence[int], primes: Sequence[int]) -> None:
    Q_t = u_to_t(Q_u)
    if Q_t[-1] != 1:
        raise AssertionError("monic construction is not monic")
    for p in primes:
        qmod = [c % p for c in Q_t]
        if len(poly_gcd_mod(qmod, derivative(qmod), p)) != 1:
            raise AssertionError(f"monic construction has bad squarefree reduction mod {p}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--X", type=int, default=40, help="even shell scale, default 40")
    parser.add_argument(
        "--variant",
        choices=("rank3", "monic", "both"),
        default="both",
        help="counterfamily to verify",
    )
    parser.add_argument(
        "--sample-only",
        action="store_true",
        help="do not scan every m in [0,X^2)",
    )
    args = parser.parse_args()

    M, m0, u0, primes, P = common_setup(args.X)
    excluded = [p for p in sieve_primes(2 * args.X) if args.X < p <= 2 * args.X and (2 * m0 + 1) % p == 0]
    print(f"X={args.X}, M={M}, m0={m0}, excluded discriminant primes={excluded}")

    if args.variant in ("rank3", "both"):
        q3 = rank3_polynomial_u(args.X, P, u0)
        verify_exact_support("rank3", args.X, q3, primes, P, M, m0, not args.sample_only)
        # Modulo every window prime this variant drops to U-u0; that feature is intentional.

    if args.variant in ("monic", "both"):
        qm, H = monic_polynomial_u(args.X, primes, P, u0)
        for p in primes:
            if poly_eval_mod(H, u0 % p, p) == 0:
                raise AssertionError(f"cofactor vanishes at anchor mod {p}")
        verify_monic_squarefree(args.X, qm, primes)
        verify_exact_support("monic-good-reduction", args.X, qm, primes, P, M, m0, not args.sample_only)

    print("ALL Q8042 COUNTERFAMILY CHECKS PASSED")


if __name__ == "__main__":
    main()
