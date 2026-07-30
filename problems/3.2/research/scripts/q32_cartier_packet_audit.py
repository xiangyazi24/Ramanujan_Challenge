#!/usr/bin/env python3
"""Exact audits for the fixed-moment Cartier-packet reformulation.

The script checks the claims used in Sections 48--50 of
Q32_SEPARATION_ANALYSIS.md:

* the Newton polytope and its Ehrhart polynomial;
* normality/IDP of the support in small dilates;
* the closed coefficient formula for Lambda^n;
* the Cartier packet congruence
      c_{ap+r}(p mu) = b_r c_a(mu) (mod p);
* the prime-square affine shell law in the quotient parameter;
* the fact that the first p-adic shell jet is not target-selective.

Only Python's standard library is used.
"""

from collections import defaultdict
from math import comb, gcd, log


def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def multiply(left, right, modulus=None):
    out = defaultdict(int)
    for u, a in left.items():
        for v, b in right.items():
            w = tuple(u[i] + v[i] for i in range(3))
            out[w] += a * b
            if modulus:
                out[w] %= modulus
    return dict(out)


def lambda_polynomial():
    one = {(0, 0, 0): 1}
    x = {(0, 0, 0): 1, (1, 0, 0): 1}
    y = {(0, 0, 0): 1, (0, 1, 0): 1}
    z = {(0, 0, 0): 1, (0, 0, 1): 1}
    bracket = multiply(y, z)
    bracket[(1, 1, 1)] = bracket.get((1, 1, 1), 0) + 1
    numerator = one
    for factor in (x, y, z, bracket):
        numerator = multiply(numerator, factor)
    return {
        (u[0] - 1, u[1] - 1, u[2] - 1): coefficient
        for u, coefficient in numerator.items()
    }


LAMBDA = lambda_polynomial()


def poly_power(base, exponent, modulus=None):
    out = {(0, 0, 0): 1}
    for _ in range(exponent):
        out = multiply(out, base, modulus)
    return out


def polytope_points(a):
    """The lattice points of aP.

    P is cut out by
        -1 <= x,y,z <= 1, x-y <= 1, x-z <= 1.
    """

    return [
        (x, y, z)
        for x in range(-a, a + 1)
        for y in range(-a, a + 1)
        for z in range(-a, a + 1)
        if x - y <= a and x - z <= a
    ]


def ehrhart(a):
    return (38 * a**3 + 57 * a**2 + 31 * a + 6) // 6


def coefficient(n, u, v, w, modulus=None):
    """Return [x^u y^v z^w] Lambda^n by a one-fold sum."""

    out = 0
    for k in range(n + 1):
        out += (
            C(n, k)
            * C(n, k - u)
            * C(2 * n - k, n - v)
            * C(2 * n - k, n - w)
        )
        if modulus:
            out %= modulus
    return out if modulus is None else out % modulus


def apery(n):
    return sum(C(n, k) ** 2 * C(n + k, k) ** 2 for k in range(n + 1))


def shell(M, d, a, modulus=None):
    """Return C_M(d) when M=a*d+r with 0 <= r < d."""

    out = 0
    for u, v, w in polytope_points(a):
        out += coefficient(M, d * u, d * v, d * w, modulus)
        if modulus:
            out %= modulus
    return out if modulus is None else out % modulus


def shell_fast(M, d, modulus=None):
    """The factorized one-sum formula for C_M(d)."""

    a = M // d
    out = 0
    for t in range(M + 1):
        x_packet = sum(C(M, M - t + d * u) for u in range(-a, a + 1))
        yz_packet = sum(
            C(2 * M - t, M - t + d * v) for v in range(-a, a + 1)
        )
        out += C(M, t) * x_packet * yz_packet**2
        if modulus:
            out %= modulus
    return out if modulus is None else out % modulus


def primes_up_to(limit):
    return [
        p
        for p in range(2, limit + 1)
        if all(p % d for d in range(2, int(p**0.5) + 1))
    ]


def valuation(n, p):
    out = 0
    while n and n % p == 0:
        n //= p
        out += 1
    return out


def audit_polytope_and_coefficients():
    assert len(LAMBDA) == 22
    assert sum(LAMBDA.values()) == 40
    assert set(LAMBDA) == set(polytope_points(1))

    for a in range(7):
        assert len(polytope_points(a)) == ehrhart(a)

    # Positivity makes Supp(Lambda^a) the a-fold Minkowski sum.  Equality
    # with aP in these cases audits the constructive IDP proof in the note.
    for a in range(7):
        power = poly_power(LAMBDA, a)
        assert set(power) == set(polytope_points(a))

    for n in range(7):
        power = poly_power(LAMBDA, n)
        for u in polytope_points(n):
            assert coefficient(n, *u) == power.get(u, 0)
        assert power.get((0, 0, 0), 0) == apery(n)


def audit_cartier_packets():
    checks = 0
    for p in primes_up_to(19):
        if p < 5:
            continue
        for a in range(1, 3):
            n0 = a * p
            packet = polytope_points(a)
            for r in range(1, p):
                n = n0 + r
                br = apery(r) % p
                power_a = poly_power(LAMBDA, a, p)
                for mu in packet:
                    left = coefficient(
                        n, *(p * coordinate for coordinate in mu), modulus=p
                    )
                    right = br * power_a[mu] % p
                    assert left == right, (p, a, r, mu, left, right)
                    checks += 1
    return checks


def audit_shell_lifts():
    checks = 0
    for p in primes_up_to(17):
        if p < 5:
            continue
        d = p - 1
        modulus = p**3
        for r in range(1, d):
            values = [
                shell(a * d + r, d, a, modulus=modulus) for a in range(6)
            ]

            # E_a is constant modulo p.
            for a in range(5):
                assert (values[a + 1] - values[a]) % p == 0

            # At a=0 the second difference has p^2 only for r>=2.
            if r >= 2:
                assert (values[2] - 2 * values[1] + values[0]) % p**2 == 0

            # Starting at a=1, the third difference has p^3.
            for a in range(3):
                third = (
                    values[a + 3]
                    - 3 * values[a + 2]
                    + 3 * values[a + 1]
                    - values[a]
                )
                required = min(3, a * d + r)
                assert third % p**required == 0
            checks += 1

    # A concrete warning against reading the packet congruence p-adically.
    # b_3=5*17^2, but the (-1,-1,-1) packet coordinate for n=20 has
    # exactly one factor 17.  The shell itself also has only one factor.
    p, r, a = 17, 3, 1
    assert valuation(apery(r), p) == 2
    vertex = coefficient(a * p + r, -a * p, -a * p, -a * p)
    assert valuation(vertex, p) == 1
    first_shell = shell(a * (p - 1) + r, p - 1, a)
    assert valuation(first_shell, p) == 1

    # The shell-to-Apery congruence is genuinely prime-specific.
    assert shell(6, 5, 1) % 6 == 1
    assert apery(1) % 6 == 5
    return checks


def audit_grid_disjointness():
    checks = 0
    for a in range(1, 8):
        points = set(polytope_points(a))
        for d in range(2, 40):
            for h in range(1, 8):
                if d <= a * h:
                    continue
                e = d + h
                left = {
                    (d * x, d * y, d * z)
                    for x, y, z in points
                    if (x, y, z) != (0, 0, 0)
                }
                right = {
                    (e * x, e * y, e * z)
                    for x, y, z in points
                    if (x, y, z) != (0, 0, 0)
                }
                assert left.isdisjoint(right)
                checks += 1
    return checks


def newton_weight(d0, L, i):
    """Primitive Lagrange weight for evaluation at d=-1."""

    return (-1) ** i * C(d0 + i, i) * C(d0 + L + 1, L - i)


def audit_newton_carrier():
    # The hostile n=321 row has targets 179, 193, 211 in its a=1
    # quotient cell.  Interpolate the fixed M=320 shell on d=161,...,211.
    M, d0, L = 320, 161, 50
    values = [shell_fast(M, d0 + i) for i in range(L + 1)]
    weights = [newton_weight(d0, L, i) for i in range(L + 1)]
    carrier = sum(weight * value for weight, value in zip(weights, values))

    assert sum(weights) == 1
    content = 0
    for weight in weights:
        content = gcd(content, abs(weight))
    assert content == 1
    assert carrier
    for target in (179, 193, 211):
        assert carrier % target == 0

    # Exact height data: primitive normalization removes all universal
    # content, but this carrier still has linear exponential height.
    return len(str(abs(carrier))), log(abs(carrier)) / M


def main():
    audit_polytope_and_coefficients()
    packet_checks = audit_cartier_packets()
    lift_checks = audit_shell_lifts()
    grid_checks = audit_grid_disjointness()
    carrier_digits, carrier_rate = audit_newton_carrier()
    print(
        "PASS:",
        "Newton/IDP/coefficient formula;",
        f"{packet_checks} Cartier coordinates;",
        f"{lift_checks} shell-lift rows;",
        f"{grid_checks} disjoint-grid rows;",
        f"n=321 Newton carrier has {carrier_digits} digits",
        f"(log-height/M={carrier_rate:.6f})",
    )


if __name__ == "__main__":
    main()
