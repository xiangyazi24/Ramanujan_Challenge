#!/usr/bin/env python3
"""Audit the star-sum proof of the weight-five double congruence.

For p >= 11 the proof uses

  A4 = sum (-1)^(k-1) binom(p-1,k)/k^4
     = H_star(1,1,1,1;p-1),

the mod-p evaluation H(1,1,4)=-B_(p-3)^2/6, and Sun's
second-order Kummer interpolation on the p-5 Bernoulli branch.
Together they prove the H(1,4) lemma needed for H6.
"""

from __future__ import annotations

from functools import lru_cache
from math import comb, gcd

import sympy as sp


def primes_below(limit: int) -> list[int]:
    return list(sp.primerange(11, limit))


@lru_cache(maxsize=None)
def bernoulli(index: int) -> sp.Rational:
    return sp.bernoulli(index)


def rational_mod(value: sp.Rational, modulus: int) -> int:
    numerator, denominator = value.as_numer_denom()
    denominator_int = int(denominator)
    assert gcd(denominator_int, modulus) == 1
    return (
        int(numerator) % modulus
        * pow(denominator_int % modulus, -1, modulus)
        % modulus
    )


def strict_sums(p: int, modulus: int) -> dict[str, int]:
    """Return the strict sums used in the proof."""
    h1 = 0
    h2 = 0
    h3 = 0
    h11 = 0
    h14 = 0
    h23 = 0
    h32 = 0
    h114 = 0
    h141 = 0
    h15 = 0
    h24 = 0

    for k in range(1, p):
        inverse = pow(k, -1, modulus)
        inverse2 = inverse * inverse % modulus
        inverse3 = inverse2 * inverse % modulus
        inverse4 = inverse2 * inverse2 % modulus
        inverse5 = inverse4 * inverse % modulus

        # Every prefix on the right has upper index < k here.
        h23 = (h23 + h2 * inverse3) % modulus
        h32 = (h32 + h3 * inverse2) % modulus
        h114 = (h114 + h11 * inverse4) % modulus
        h141 = (h141 + h14 * inverse) % modulus
        h15 = (h15 + h1 * inverse5) % modulus
        h24 = (
            h24
            + sum(pow(j, -2, modulus) for j in range(1, k))
            * inverse4
        ) % modulus

        h11 = (h11 + h1 * inverse) % modulus
        h14 = (h14 + h1 * inverse4) % modulus
        h1 = (h1 + inverse) % modulus
        h2 = (h2 + inverse2) % modulus
        h3 = (h3 + inverse3) % modulus

    return {
        "H14": h14,
        "H23": h23,
        "H32": h32,
        "H114": h114,
        "H141": h141,
        "H15": h15,
        "H24": h24,
    }


def audit_prime(p: int) -> dict[str, bool]:
    p2 = p * p
    p3 = p2 * p
    inverses = [0] + [pow(k, -1, p3) for k in range(1, p)]
    single = {
        exponent: sum(
            pow(inverses[k], exponent, p3)
            for k in range(1, p)
        )
        % p3
        for exponent in range(1, 7)
    }

    strict_p = strict_sums(p, p)
    strict_p2 = strict_sums(p, p2)
    h14 = strict_p2["H14"]
    h23 = strict_p2["H23"]
    h32 = strict_p2["H32"]
    h114 = strict_p["H114"]
    h141 = strict_p["H141"]
    h15 = strict_p["H15"]
    h24 = strict_p["H24"]

    additive_g23 = sum(
        pow(u, -2, p2) * pow(v, -3, p2)
        for u in range(1, p)
        for v in range(1, p - u)
    ) % p2

    # V=sum_k H_k(1,1)/k^4 = H(1,1,4)+H(1,5).
    prefix1 = 0
    prefix11 = 0
    inclusive_v = 0
    for k in range(1, p):
        inverse = pow(k, -1, p)
        prefix11 = (prefix11 + prefix1 * inverse) % p
        prefix1 = (prefix1 + inverse) % p
        inclusive_v = (
            inclusive_v + prefix11 * pow(inverse, 4, p)
        ) % p

    a4 = sum(
        (-1 if k % 2 == 0 else 1)
        * comb(p - 1, k)
        * pow(inverses[k], 4, p3)
        for k in range(1, p)
    ) % p3
    star_newton = (
        pow(single[1], 4, p3)
        + 6 * single[1] * single[1] * single[2]
        + 3 * single[2] * single[2]
        + 8 * single[1] * single[3]
        + 6 * single[4]
    ) * pow(24, -1, p3) % p3
    star_reduced = (
        single[2] * single[2] * pow(8, -1, p3)
        + single[4] * pow(4, -1, p3)
    ) % p3
    binomial_expansion = (
        -single[4]
        + p * (h14 + single[5])
        - p2 * inclusive_v
    ) % p3

    b3_p = rational_mod(bernoulli(p - 3), p)
    b3_square_p = b3_p * b3_p % p
    beta0 = (
        rational_mod(bernoulli(p - 5), p2)
        * pow(p - 5, -1, p2)
        % p2
    )
    beta1 = (
        rational_mod(bernoulli(2 * p - 6), p2)
        * pow(2 * p - 6, -1, p2)
        % p2
    )

    s4_bridge = p * (-8 * beta0 + 4 * beta1) % p3
    h14_target = (
        single[5] * pow(2, -1, p2)
        - 10 * beta0
        + 5 * beta1
        - p * b3_square_p * pow(9, -1, p2)
    ) % p2

    # Inclusive D=sum_k H_k^(2)/k^3 for H6.
    prefix2 = 0
    inclusive_d = 0
    for k in range(1, p):
        inverse = inverses[k]
        prefix2 = (
            prefix2 + inverse * inverse
        ) % p2
        inclusive_d = (
            inclusive_d + prefix2 * pow(inverse, 3, p2)
        ) % p2
    h6 = (
        single[2] * single[2]
        - 5 * single[4]
        - 2 * p * inclusive_d
    ) % p3

    return {
        "star_exact": a4 == star_newton,
        "star_reduction": a4 == star_reduced,
        "binomial_expansion": a4 == binomial_expansion,
        "inclusive_split": inclusive_v == (h114 + h15) % p,
        "zhao_middle": h141
        == b3_square_p * pow(3, -1, p) % p,
        "stuffle_left": (
            2 * h114 + h141 + h24 + h15
        )
        % p
        == 0,
        "e1_partial_fraction": additive_g23
        == (
            6 * h14 + 3 * h23 + h32
        )
        % p2,
        "e1_additive_reflection": additive_g23
        == (
            -h23 - 3 * p * strict_p2["H24"]
        )
        % p2,
        "e1_stuffle": (
            h23 + h32 + single[5] - single[2] * single[3]
        )
        % p2
        == 0,
        "e1": (
            h23 + 2 * h14
        )
        % p2
        == 0,
        "triple_value": h114
        == -b3_square_p * pow(6, -1, p) % p,
        "s2_classical": single[2] % p2
        == 2 * p * b3_p * pow(3, -1, p2) % p2,
        "s5_vanishing": single[5] % p2 == 0,
        "s4_bridge": single[4] == s4_bridge,
        "h14_target": h14 == h14_target,
        "h6": h6 == 0,
    }


def symbolic_kummer_recombination() -> bool:
    p = sp.symbols("p")
    beta0, beta1 = sp.symbols("beta0 beta1")
    expression = 0
    for index, coefficient in ((3, 15), (4, -24), (5, 10)):
        exponent = p - 5 + index * (p - 1)
        beta_index = index * beta1 - (index - 1) * beta0
        expression += coefficient * exponent * beta_index
    return sp.expand(expression) == -8 * beta0 + 4 * beta1


def symbolic_partial_fraction() -> bool:
    u, v = sp.symbols("u v")
    total = u + v
    right = (
        3 / (u * total**4)
        + 1 / (u**2 * total**3)
        + 3 / (v * total**4)
        + 2 / (v**2 * total**3)
        + 1 / (v**3 * total**2)
    )
    return sp.cancel(right - 1 / (u**2 * v**3)) == 0


def main() -> None:
    primes = primes_below(300)
    counters: dict[str, int] = {}
    first_failure: dict[str, int] = {}
    for prime in primes:
        result = audit_prime(prime)
        for name, passed in result.items():
            counters[name] = counters.get(name, 0) + int(passed)
            if not passed and name not in first_failure:
                first_failure[name] = prime

    print(
        "symbolic_kummer_recombination="
        f"{int(symbolic_kummer_recombination())}/1"
    )
    print(
        "symbolic_partial_fraction="
        f"{int(symbolic_partial_fraction())}/1"
    )
    for name in sorted(counters):
        print(f"{name}={counters[name]}/{len(primes)}")
    print(f"first_failure={first_failure}")
    assert symbolic_kummer_recombination()
    assert symbolic_partial_fraction()
    assert not first_failure
    print("failures=0")


if __name__ == "__main__":
    main()
