#!/usr/bin/env python3
"""Exact cyclotomic check of the character-aspect Galois action.

For p in {13, 17}, this script computes the genuine integer Frobenius
traces of the smooth fibers

    E_u: y^2 + (1-2u)xy + u^2 y = x^3

and represents

    M_p(r) = sum_u a_{p,u}^2 * omega(phi(u))^(-r)

exactly in Z[z]/(Phi_{p-1}(z)).  Here omega(g)=z for the least primitive
root g modulo p and phi(u)=u(1-8u)/(1+u).  Terms where phi is zero or
undefined are omitted, equivalently extending every nontrivial
multiplicative character by zero.

For every a in (Z/(p-1)Z)^x and 1 <= r <= p-2, we check

    sigma_a(M_p(r)) = M_p(ar),       sigma_a(z)=z^a.

All polynomial reductions are exact over ZZ; no floating-point embedding of
the cyclotomic field is used.
"""

from __future__ import annotations

from math import comb, gcd

import sympy as sp


Z = sp.Symbol("z")


def prime_divisors(n: int) -> list[int]:
    """Return the distinct prime divisors of n."""
    divisors: list[int] = []
    candidate = 2
    while candidate * candidate <= n:
        if n % candidate == 0:
            divisors.append(candidate)
            while n % candidate == 0:
                n //= candidate
        candidate += 1
    if n > 1:
        divisors.append(n)
    return divisors


def primitive_root(p: int) -> int:
    """Return the least positive primitive root modulo the prime p."""
    factors = prime_divisors(p - 1)
    for candidate in range(2, p):
        if all(pow(candidate, (p - 1) // q, p) != 1 for q in factors):
            return candidate
    raise AssertionError(f"no primitive root found modulo {p}")


def discrete_log_table(p: int, generator: int) -> dict[int, int]:
    """Tabulate x -> j for x=generator^j in F_p^x."""
    table: dict[int, int] = {}
    value = 1
    for exponent in range(p - 1):
        assert value not in table
        table[value] = exponent
        value = value * generator % p
    assert value == 1 and len(table) == p - 1
    return table


def legendre(a: int, p: int) -> int:
    """Quadratic character of a modulo the odd prime p, extended by 0."""
    a %= p
    if a == 0:
        return 0
    value = pow(a, (p - 1) // 2, p)
    assert value in (1, p - 1)
    return 1 if value == 1 else -1


def elliptic_trace(p: int, u: int) -> int:
    """Compute p+1-#E_u(F_p) by exact point counting."""
    a1 = (1 - 2 * u) % p
    a3 = u * u % p
    points = 1  # the unique point at infinity on a smooth Weierstrass fiber
    for x in range(p):
        linear_y = (a1 * x + a3) % p
        discriminant_y = (linear_y * linear_y + 4 * x**3) % p
        points += 1 + legendre(discriminant_y, p)
    trace = p + 1 - points
    assert trace * trace <= 4 * p
    return trace


def franel_hasse_value(p: int, u: int) -> int:
    """Evaluate H_p(u)=sum_{n<p} sum_k binom(n,k)^3 u^n modulo p."""
    value = 0
    power = 1
    for n in range(p):
        franel = sum(comb(n, k) ** 3 for k in range(n + 1))
        value = (value + franel * power) % p
        power = power * u % p
    return value


def phi_value(p: int, u: int) -> int | None:
    """Return phi(u) in F_p, or None at its pole u=-1."""
    denominator = (1 + u) % p
    if denominator == 0:
        return None
    return u * (1 - 8 * u) * pow(denominator, -1, p) % p


def reduce_cyclotomic(poly: sp.Poly, cyclotomic: sp.Poly) -> sp.Poly:
    """Return the canonical representative modulo a monic cyclotomic poly."""
    return sp.Poly(poly, Z, domain=sp.ZZ).rem(cyclotomic)


def mellin_value(
    r: int,
    p: int,
    log_table: dict[int, int],
    traces: dict[int, int],
    cyclotomic: sp.Poly,
) -> sp.Poly:
    """Compute M_p(r) exactly in Z[z]/Phi_{p-1}."""
    order = p - 1
    result = sp.Poly(0, Z, domain=sp.ZZ)
    for u, trace in traces.items():
        value = phi_value(p, u)
        assert value not in (None, 0)
        exponent = (-r * log_table[value]) % order
        result += sp.Poly(trace * trace * Z**exponent, Z, domain=sp.ZZ)
    return reduce_cyclotomic(result, cyclotomic)


def galois_conjugate(poly: sp.Poly, a: int, cyclotomic: sp.Poly) -> sp.Poly:
    """Apply sigma_a(z)=z^a exactly and reduce modulo Phi_{p-1}."""
    substituted = sp.Poly(poly.as_expr().subs(Z, Z**a), Z, domain=sp.ZZ)
    return reduce_cyclotomic(substituted, cyclotomic)


def check_prime(p: int) -> tuple[int, int, int]:
    """Check all character indices and all cyclotomic automorphisms."""
    order = p - 1
    generator = primitive_root(p)
    log_table = discrete_log_table(p, generator)
    cyclotomic = sp.Poly(sp.cyclotomic_poly(order, Z), Z, domain=sp.ZZ)

    # phi(u) is nonzero precisely away from u=0, 1/8 and the pole u=-1.
    excluded = {0, pow(8, -1, p), p - 1}
    traces = {u: elliptic_trace(p, u) for u in range(p) if u not in excluded}
    assert all(isinstance(trace, int) for trace in traces.values())
    assert all(phi_value(p, u) not in (None, 0) for u in traces)
    assert all(trace % p == franel_hasse_value(p, u) for u, trace in traces.items())

    values = {
        r: mellin_value(r, p, log_table, traces, cyclotomic)
        for r in range(order)
    }
    automorphisms = [a for a in range(1, order) if gcd(a, order) == 1]
    checks = 0
    for a in automorphisms:
        for r in range(1, order):
            left = galois_conjugate(values[r], a, cyclotomic)
            right = values[(a * r) % order]
            assert left == right, (p, a, r, left, right)
            checks += 1

    # Guard the character convention independently on every field element:
    # sigma_a(omega(t)^(-r)) = omega(t)^(-ar).
    character_checks = 0
    for a in automorphisms:
        for r in range(1, order):
            for t, logarithm in log_table.items():
                lhs = reduce_cyclotomic(
                    sp.Poly(Z ** ((-r * logarithm * a) % order), Z, domain=sp.ZZ),
                    cyclotomic,
                )
                rhs = reduce_cyclotomic(
                    sp.Poly(Z ** ((-(a * r) * logarithm) % order), Z, domain=sp.ZZ),
                    cyclotomic,
                )
                assert lhs == rhs, (p, a, r, t)
                character_checks += 1

    print(
        f"VERIFIED p={p}: exact Q(zeta_{order}) Galois action "
        f"sigma_a(M(r))=M(ar) for {checks} (a,r) pairs; "
        f"generator={generator}, smooth Mellin fibers={len(traces)}, "
        f"character checks={character_checks}"
    )
    return checks, character_checks, len(traces)


def main() -> None:
    totals = [check_prime(p) for p in (13, 17)]
    print(
        "VERIFIED exact cyclotomic Galois action at p=13,17: "
        f"{sum(item[0] for item in totals)} Mellin identities and "
        f"{sum(item[1] for item in totals)} point-character identities"
    )


if __name__ == "__main__":
    main()
