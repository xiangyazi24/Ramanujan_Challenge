#!/usr/bin/env python3
"""Self-contained exact verifier for Q8239.

This file uses only Python's standard library.  It reconstructs the scalar
Cartier shell, terminal Newton family, Pascal-normalized Turan carriers, the
first-jet alias, two natural 2x2 determinantal families, and the hostile
n=200 certificate discussed in

    ORACLE_COMM/chatgpt_turan_nonzero_certificate_attack.md

No floating-point arithmetic is used for any asserted identity.
"""

from __future__ import annotations

import argparse
from math import comb, gcd, isqrt, prod


def C(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


def primes_up_to(limit: int) -> list[int]:
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (
                (limit - p * p) // p + 1
            )
    return [p for p in range(2, limit + 1) if sieve[p]]


def valuation(value: int, prime: int) -> int:
    value = abs(value)
    exponent = 0
    while value and value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def apery_mod(index: int, prime: int) -> int:
    return sum(
        (C(index, k) * C(index + k, k)) ** 2
        for k in range(index + 1)
    ) % prime


def binomial_row(n: int) -> list[int]:
    row = [1]
    for k in range(n):
        row.append(row[-1] * (n - k) // (k + 1))
    return row


def shell_batch(moment: int, nodes) -> dict[int, int]:
    """Exact C_M(d) for the Section-48 Laurent shell.

    This is the dependency-free batched form of the one-fold coefficient
    formula used by q32_cartier_packet_audit.py.
    """

    nodes = tuple(dict.fromkeys(nodes))
    assert all(1 <= d <= moment for d in nodes)
    moment_row = binomial_row(moment)
    quotients = {d: moment // d for d in nodes}
    answer = {d: 0 for d in nodes}

    for t in range(moment + 1):
        upper = binomial_row(2 * moment - t)
        base = moment - t
        outer = moment_row[t]
        for d in nodes:
            q = quotients[d]
            x_packet = sum(
                moment_row[index]
                for u in range(-q, q + 1)
                if 0 <= (index := base + d * u) <= moment
            )
            yz_packet = sum(
                upper[index]
                for v in range(-q, q + 1)
                if 0 <= (index := base + d * v) < len(upper)
            )
            answer[d] += outer * x_packet * yz_packet * yz_packet
    return answer


def newton(values: dict[int, int], start: int, order: int) -> int:
    return sum(
        (-1) ** i
        * C(start + i, i)
        * C(start + order + 1, order - i)
        * values[start + i]
        for i in range(order + 1)
    )


def terminal_data(n: int):
    """Return (start,L0,F,P) with F_j=f_{L0-j}."""

    moment = n - 1
    start = moment // 2 + 1
    L0 = moment - start
    values = shell_batch(moment, range(start, moment + 1))
    f = [
        newton(values, moment - order, order)
        for order in range(L0 + 1)
    ]
    F = tuple(reversed(f))
    P = tuple(C(n, L0 - j) for j in range(L0 + 1))
    return start, L0, F, P


def egcd(a: int, b: int):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def centered_pair_bezout(a: int, b: int):
    """Return a centered Bezout solution x*a+y*b=gcd(a,b)."""

    g, x, _ = egcd(a, b)
    step = abs(b // g)
    if step:
        x %= step
        if 2 * x > step:
            x -= step
    y = (g - x * a) // b
    assert x * a + y * b == g
    return g, x, y


def exact_family(n: int, r: int):
    """Build F_0..F_{r+1}, E_1..E_r, D_0..D_r, R_1..R_r."""

    start, L0, F_all, P_all = terminal_data(n)
    assert 2 <= r <= L0 - 1
    F = F_all[: r + 2]
    P = P_all[: r + 1]

    D = tuple(F[j] - F[j + 1] for j in range(r + 1))
    q = tuple(gcd(P[j - 1], P[j]) for j in range(1, r + 1))
    E = []
    R = []

    for j in range(1, r + 1):
        numerator = F[j - 1] * F[j + 1] - F[j] * F[j]
        assert numerator % q[j - 1] == 0
        e = numerator // q[j - 1]
        rr_num = D[j - 1] - D[j]
        assert rr_num % q[j - 1] == 0
        rr = rr_num // q[j - 1]
        assert D[j - 1] * D[j] % q[j - 1] == 0
        assert e == F[j] * rr - D[j - 1] * D[j] // q[j - 1]
        E.append(e)
        R.append(rr)

        k = L0 - j
        gsmall = gcd(n - k, k + 1)
        A = (n - k) // gsmall
        B = (k + 1) // gsmall
        assert gcd(A, B) == 1
        assert P[j - 1] // q[j - 1] == A
        assert P[j] // q[j - 1] == B
        assert D[j - 1] % P[j - 1] == 0
        assert D[j] % P[j] == 0
        X_left = D[j - 1] // P[j - 1]
        X_right = D[j] // P[j]
        assert rr == A * X_left - B * X_right

    E = tuple(E)
    R = tuple(R)

    W = []
    for j in range(r - 1):
        w = R[j + 1] * E[j] - R[j] * E[j + 1]
        rhs_num = D[j + 1] * (
            D[j] * D[j + 2] - D[j + 1] * D[j + 1]
        )
        rhs_den = q[j] * q[j + 1]
        assert rhs_num % rhs_den == 0
        assert w == rhs_num // rhs_den
        W.append(w)

    # Columns are (F_j,E_j), j=1,...,r.
    FE_minors = []
    for i in range(1, r + 1):
        for j in range(i + 1, r + 1):
            FE_minors.append(E[i - 1] * F[j] - E[j - 1] * F[i])

    return {
        "n": n,
        "start": start,
        "L0": L0,
        "F": F,
        "P": P,
        "D": D,
        "q": q,
        "E": E,
        "R": R,
        "W": tuple(W),
        "FE_minors": tuple(FE_minors),
    }


def common_candidates(n: int, start: int, r: int) -> tuple[int, ...]:
    # F_0,...,F_{r+1} all contain p-1 exactly when p > start+r+1.
    # p=n is deliberately excluded: h=n-p=0 and b_0=1.
    return tuple(
        p for p in primes_up_to(n - 1) if start + r + 1 < p < n
    )


def audit_candidate_congruences(record):
    n = record["n"]
    start = record["start"]
    r = len(record["E"])
    F = record["F"]
    P = record["P"]
    q = record["q"]
    E = record["E"]
    R = record["R"]
    W = record["W"]

    candidates = common_candidates(n, start, r)
    targets = []
    for p in candidates:
        h = n - p
        assert 1 <= h <= p - 2
        b = apery_mod(h, p)
        for j in range(r + 2):
            assert F[j] % p == b
        for j in range(1, r + 1):
            assert valuation(P[j - 1], p) == 1
            assert valuation(P[j], p) == 1
            assert valuation(q[j - 1], p) == 1
            assert (E[j - 1] - b * R[j - 1]) % p == 0
        for w in W:
            assert w % p == 0
        if b == 0:
            targets.append(p)
            for e in E:
                assert e % p == 0

    # FE minors have the target-square law, but no universal candidate law.
    minors = record["FE_minors"]
    for p in candidates:
        b = apery_mod(n - p, p)
        if b == 0:
            assert all(value % (p * p) == 0 for value in minors)
        # Check the exact non-target first-jet reduction on every 2x2 minor.
        index = 0
        for i in range(1, r + 1):
            for j in range(i + 1, r + 1):
                expected = b * b * (R[i - 1] - R[j - 1])
                assert (minors[index] - expected) % p == 0
                index += 1

    if n in primes_up_to(n):
        assert apery_mod(0, n) == 1
    return candidates, tuple(targets)


def gcd_many(values) -> int:
    answer = 0
    for value in values:
        answer = gcd(answer, abs(value))
    return answer


def hostile_n200():
    rec = exact_family(200, 6)
    candidates, targets = audit_candidate_congruences(rec)
    expected_candidates = (
        109, 113, 127, 131, 137, 139, 149, 151, 157,
        163, 167, 173, 179, 181, 191, 193, 197, 199,
    )
    assert candidates == expected_candidates
    assert targets == (139, 181)

    # The first-jet Wronskian is universally candidate-divisible.  Its
    # primitive quotient loses the actual targets: there is no second target
    # digit hiding in this determinant.
    candidate_product = prod(candidates)
    wgcd = gcd_many(rec["W"])
    assert wgcd % candidate_product == 0
    wquot = wgcd // candidate_product
    assert wquot % 139 != 0
    assert wquot % 181 != 0
    for p in candidates:
        assert valuation(wgcd, p) == 1

    # The value--Turan second determinantal divisor does not acquire the
    # whole candidate primorial.  In this row it is exactly target^2*50.
    fe_delta2 = gcd_many(rec["FE_minors"])
    assert fe_delta2 == 2 * 5**2 * 139**2 * 181**2
    for p in targets:
        assert valuation(fe_delta2, p) == 2
    for p in candidates:
        if p not in targets:
            assert fe_delta2 % p != 0

    # A small gcd is not a bounded-coefficient Bezout certificate.
    e1, e2 = rec["E"][:2]
    g, x, y = centered_pair_bezout(e1, e2)
    assert g == 2**2 * 5 * 139 * 181
    stepx = abs(e2 // g)
    stepy = abs(e1 // g)
    assert 2 * abs(x) < stepx
    assert 2 * abs(y) < stepy
    # Hence t=0 is coordinatewise minimal in the full family
    # (x+t*e2/g, y-t*e1/g) of Bezout solutions.
    assert abs(x).bit_length() == 2092
    assert abs(y).bit_length() == 2091
    assert abs(e1).bit_length() == 2112
    assert abs(e2).bit_length() == 2112

    print("N200_CANDIDATES", candidates)
    print("N200_TARGETS", targets)
    print("N200_W_GCD_BITS", wgcd.bit_length())
    print("N200_W_PRIMITIVE", wquot)
    print("N200_FE_DELTA2", fe_delta2)
    print("N200_FE_DELTA2_FACTORS", "2 * 5^2 * 139^2 * 181^2")
    print("N200_E12_GCD", g)
    print("N200_MIN_PAIR_BEZOUT_BITS", (abs(x).bit_length(), abs(y).bit_length()))
    return rec


def extended_records():
    expected = {
        200: ((139, 181), 2 * 5**2 * 139**2 * 181**2),
        272: ((191, 233), 2 * 191**2 * 233**2),
        300: ((191, 227), 37 * 191**2 * 227**2),
        321: ((179, 193, 211), 179**2 * 193**2 * 211**2),
    }
    for n, (wanted_targets, wanted_delta) in expected.items():
        rec = exact_family(n, 6)
        _, targets = audit_candidate_congruences(rec)
        delta = gcd_many(rec["FE_minors"])
        assert targets == wanted_targets
        assert delta == wanted_delta
        print("RECORD", n, "TARGETS", targets, "FE_DELTA2", delta, "BITS", delta.bit_length())


def small_identity_sweep():
    total = 0
    for n in range(12, 61):
        moment = n - 1
        start = moment // 2 + 1
        L0 = moment - start
        if L0 < 3:
            continue
        r = min(5, L0 - 1)
        rec = exact_family(n, r)
        audit_candidate_congruences(rec)
        total += 1
    print("SMALL_IDENTITY_ROWS", total)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--extended", action="store_true", help="also run n=272,300,321 records")
    parser.add_argument("--small-sweep", action="store_true", help="audit all 12<=n<=60")
    args = parser.parse_args()

    hostile_n200()
    if args.extended:
        extended_records()
    if args.small_sweep:
        small_identity_sweep()
    print("Q8239_TURAN_NONZERO_CERTIFICATE_VERIFY=PASS")


if __name__ == "__main__":
    main()
