#!/usr/bin/env python3
"""Independent exact verifier for the Q7694 transverse cross-row audit.

The script uses only Python's standard library.  It reconstructs the Apéry
sequence b, the series g = 1/(F^2 sqrt(1-34t+t^2)), Xi, the homogeneous
unit-Casoratian companion u, Phi, and the inhomogeneous coordinate kappa.

It then checks, in exact arithmetic,

    r^3 (b_{r-1} u_r - b_r u_{r-1}) = 1,
    kappa_r = Xi_r u_r + Phi_r b_r,
    kappa_{r-1} = Xi_r u_{r-1} + Phi_r b_{r-1},
    Xi_r = r^3 (b_{r-1} kappa_r - b_r kappa_{r-1}),
    Phi_r = r^3 (u_r kappa_{r-1} - u_{r-1} kappa_r),
    A_r kappa = -5 g_r  (r >= 2),

with the genuine boundary defect A_1 kappa = -36 = -5 g_1 - 1.

For block heights it writes kappa_r in lowest terms K_r/d_r and reports
raw K-height, maximal small-prime-saturated K-height (all q <= 2R removed),
and the actual high-prime transverse target.  No finite height table is promoted
to an asymptotic theorem.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd, isqrt, log


def P(n: int) -> int:
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def primes_up_to(limit: int) -> list[int]:
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (
                (limit - start) // p + 1
            )
    return [n for n, flag in enumerate(sieve) if flag]


def strip_primes_leq(value: int, bound: int, primes: list[int]) -> int:
    """Remove the full q-adic content for every prime q <= bound."""
    value = abs(value)
    if value == 0:
        return 0
    for q in primes:
        if q > bound:
            break
        while value % q == 0:
            value //= q
    return value


def reconstruct(N: int):
    # Apéry b_n.
    b = [0] * (N + 1)
    b[0] = 1
    if N >= 1:
        b[1] = 5
    for n in range(1, N):
        num = P(n) * b[n] - n**3 * b[n - 1]
        den = (n + 1) ** 3
        assert num % den == 0, ("b recurrence integrality", n)
        b[n + 1] = num // den

    # s = sqrt(D), D = 1 - 34t + t^2.  The coefficient equation is
    # 2 s_n + sum_{i=1}^{n-1} s_i s_{n-i} = [t^n]D.
    s = [0] * (N + 1)
    s[0] = 1
    for n in range(1, N + 1):
        dn = -34 if n == 1 else (1 if n == 2 else 0)
        cross = sum(s[i] * s[n - i] for i in range(1, n))
        num = dn - cross
        assert num % 2 == 0, ("sqrt(D) integrality", n, num)
        s[n] = num // 2

    # h = F^2 sqrt(D).
    f2 = [0] * (N + 1)
    for n in range(N + 1):
        f2[n] = sum(b[i] * b[n - i] for i in range(n + 1))
    h = [0] * (N + 1)
    for n in range(N + 1):
        h[n] = sum(f2[i] * s[n - i] for i in range(n + 1))
    assert h[0] == 1

    # g = 1/h.
    g = [0] * (N + 1)
    g[0] = 1
    for n in range(1, N + 1):
        g[n] = -sum(h[i] * g[n - i] for i in range(1, n + 1))

    Xi = [0] * (N + 1)
    Xi[0] = -1
    for r in range(1, N + 1):
        Xi[r] = Xi[r - 1] - 5 * g[r] * b[r - 1]

    # kappa is reconstructed independently from its inhomogeneous source law.
    # The r=1 initial value is fixed by the frame, and is the source boundary
    # defect that must not be silently replaced by -5*g_1.
    kappa = [Fraction(0) for _ in range(N + 1)]
    if N >= 1:
        kappa[1] = Fraction(-36)
    for r in range(2, N + 1):
        kappa[r] = Fraction(
            P(r - 1) * kappa[r - 1]
            - (r - 1) ** 3 * kappa[r - 2]
            - 5 * g[r],
            r**3,
        )

    return b, s, h, g, Xi, kappa


def verify_frame(b, g, Xi, kappa, frame_N: int) -> None:
    # Homogeneous companion u_0=0,u_1=1 and variation coordinate Phi.
    u = [Fraction(0) for _ in range(frame_N + 1)]
    if frame_N >= 1:
        u[1] = Fraction(1)
    for r in range(1, frame_N):
        u[r + 1] = Fraction(
            P(r) * u[r] - r**3 * u[r - 1], (r + 1) ** 3
        )

    Phi = [Fraction(0) for _ in range(frame_N + 1)]
    for r in range(1, frame_N + 1):
        Phi[r] = Phi[r - 1] + 5 * g[r] * u[r - 1]

    for r in range(1, frame_N + 1):
        cas = r**3 * (b[r - 1] * u[r] - b[r] * u[r - 1])
        assert cas == 1, ("unit Casoratian", r, cas)

        direct = Xi[r] * u[r] + Phi[r] * b[r]
        assert direct == kappa[r], ("frame", r, direct, kappa[r])

        sync = Xi[r] * u[r - 1] + Phi[r] * b[r - 1]
        assert sync == kappa[r - 1], ("synchronized previous row", r)

        xi_cross = r**3 * (
            b[r - 1] * kappa[r] - b[r] * kappa[r - 1]
        )
        assert xi_cross == Xi[r], ("cross-row Xi", r, xi_cross, Xi[r])

        phi_cross = r**3 * (
            u[r] * kappa[r - 1] - u[r - 1] * kappa[r]
        )
        assert phi_cross == Phi[r], ("cross-row Phi", r)

    # Summation-by-parts / Green identity.  It telescopes to the same adjacent
    # determinant, so the research note will not count it as a second height
    # improvement.
    for a in range(2, min(frame_N, 20) + 1):
        acc = Fraction(0)
        for m in range(a, frame_N + 1):
            acc += b[m - 1] * (-5 * g[m])
            assert acc == Xi[m] - Xi[a - 1], (
                "Green/telescoping identity",
                a,
                m,
            )


def verify_source(g, kappa, N: int) -> None:
    if N >= 1:
        A1 = kappa[1] - P(0) * kappa[0]
        assert A1 == -36
        assert A1 == -5 * g[1] - 1
        assert A1 != -5 * g[1]
    for r in range(2, N + 1):
        Ar = (
            r**3 * kappa[r]
            - P(r - 1) * kappa[r - 1]
            + (r - 1) ** 3 * kappa[r - 2]
        )
        assert Ar == -5 * g[r], ("source transport", r, Ar, -5 * g[r])


def denominator_supported_at_most_r(kappa, N: int) -> None:
    primes = primes_up_to(N)
    for r in range(1, N + 1):
        rem = kappa[r].denominator
        for q in primes:
            if q > r:
                break
            while rem % q == 0:
                rem //= q
        assert rem == 1, ("large prime in kappa denominator", r, rem)


def verify_locked_regressions(b, Xi, kappa) -> None:
    for p, r in ((17, 13), (2237, 492)):
        assert b[r] % p == 0, ("locked common b", p, r)
        assert Xi[r] % p == 0, ("locked common Xi", p, r)
        assert kappa[r].denominator % p != 0
        assert kappa[r].numerator % p == 0, ("locked common kappa", p, r)

    p, r = 11, 5
    assert b[r] % p == 0, ("locked nontransverse b", p, r)
    assert Xi[r] % p != 0, ("locked nontransverse Xi", p, r)
    assert kappa[r].denominator % p != 0
    assert kappa[r].numerator % p != 0, ("locked nontransverse kappa", p, r)


def verify_small_values(b, g, Xi, kappa) -> None:
    assert b[:6] == [1, 5, 73, 1445, 33001, 819005]
    assert g[:5] == [1, 7, 192, 5520, 165168]
    assert Xi[:4] == [-1, -36, -4836, -2019636]
    assert kappa[0] == 0
    assert kappa[1] == -36
    assert kappa[2] == Fraction(-1293, 2)
    assert kappa[3] == Fraction(-82931, 6)
    assert 8 * (5 * kappa[2] - 73 * kappa[1]) == Xi[2]
    assert 27 * (73 * kappa[3] - 1445 * kappa[2]) == Xi[3]


def log_height(values) -> float:
    return sum(log(max(1, abs(int(v)))) for v in values)


def block_height_table(b, Xi, kappa, Rs: list[int]) -> None:
    maxR = max(Rs)
    primes = primes_up_to(2 * maxR)
    print("\nBLOCK HEIGHT AUDIT")
    print(
        "R rawK/R^2 satK/R^2 target/R^2 rawXi/R^2 "
        "target_rows target_eq"
    )
    for R in Rs:
        raw_K = []
        sat_K = []
        raw_Xi = []
        target_Xi = []
        target_K = []
        target_rows = 0
        for r in range(R + 1, 2 * R + 1):
            K = kappa[r].numerator  # primitive numerator; Fraction is reduced
            raw_K.append(K)
            sat_K.append(strip_primes_leq(K, 2 * R, primes))
            raw_Xi.append(Xi[r])

            gx = gcd(abs(b[r]), abs(Xi[r]))
            gk = gcd(abs(b[r]), abs(K))
            hx = strip_primes_leq(gx, 2 * R, primes)
            hk = strip_primes_leq(gk, 2 * R, primes)
            assert hx == hk, ("high target mismatch", R, r, hx, hk)
            target_Xi.append(hx)
            target_K.append(hk)
            if hx > 1:
                target_rows += 1

        rawK_h = log_height(raw_K)
        satK_h = log_height(sat_K)
        target_h = log_height(target_K)
        rawXi_h = log_height(raw_Xi)
        assert target_K == target_Xi
        print(
            f"{R:4d} {rawK_h/R**2: .9f} {satK_h/R**2: .9f} "
            f"{target_h/R**2: .9f} {rawXi_h/R**2: .9f} "
            f"{target_rows:4d} yes"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=520)
    args = parser.parse_args()
    N = args.N
    if N < 492:
        raise SystemExit("--N must be at least 492 for the locked (2237,492) check")

    b, s, h, g, Xi, kappa = reconstruct(N)
    verify_small_values(b, g, Xi, kappa)
    verify_source(g, kappa, N)
    verify_frame(b, g, Xi, kappa, min(N, 80))
    denominator_supported_at_most_r(kappa, N)
    verify_locked_regressions(b, Xi, kappa)

    candidate_Rs = [8, 16, 32, 64, 128, 256]
    Rs = [R for R in candidate_Rs if 2 * R <= N]
    block_height_table(b, Xi, kappa, Rs)

    print("\nLOCKED REGRESSIONS")
    print("common: (17,13) PASS; (2237,492) PASS")
    print("nontransverse target: (11,5) PASS")
    print("fixed six-slope passport R=1024,p=4013: NOT USED (locked failure preserved)")
    print("reflected-depth laws: NOT USED")
    print("\nQ7694 VERIFIER: PASS")


if __name__ == "__main__":
    main()
