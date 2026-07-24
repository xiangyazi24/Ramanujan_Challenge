#!/usr/bin/env python3
"""Audit the candidate-window support of the primitive a=1 Padé numerator.

Q850 gives, for interpolation height H,

    P_H^prim(x) = m_H z_H(x),

where z_H is primitive linear and

    m_H = lcm_{0 <= s <= H} A_s / gcd(A_s, z_H(s)).

This script computes the required local valuations without constructing
the degree H-1 denominator, using Fraction arithmetic, or forming the
global multiplier.  At n=3H+1 it classifies every prime 2H < p <= n
dividing P_H^prim(n) as a multiplier prime, a root prime, a genuine q=1
target, or good-prime pollution.

The computation is exact but is only an experiment.  Sparse support in a
finite range is not an asymptotic theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import gcd, log

from q32_pade_total_positivity import apery_values


def integer_lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def prime_sieve(limit: int) -> list[bool]:
    prime = [True] * (limit + 1)
    prime[:2] = [False, False]
    for divisor in range(2, int(limit**0.5) + 1):
        if not prime[divisor]:
            continue
        prime[divisor * divisor : limit + 1 : divisor] = [False] * (
            (limit - divisor * divisor) // divisor + 1
        )
    return prime


def prefix_lcms(values: list[int]) -> list[int]:
    result: list[int] = []
    current = 1
    for value in values:
        current = integer_lcm(current, value)
        result.append(current)
    return result


def cleared_moments(
    height: int,
    apery: list[int],
    apery_lcm: int,
) -> tuple[int, int]:
    """Return X_H,Y_H with z_H(x)=(Y_H-X_H*x)/gcd(X_H,Y_H)."""

    x_moment = 0
    y_moment = 0
    binomial = 1
    for node in range(height + 1):
        signed_binomial = -binomial if node % 2 else binomial
        cleared_reciprocal = apery_lcm // apery[node]
        term = signed_binomial * cleared_reciprocal
        x_moment += term
        y_moment += node * term
        if node < height:
            binomial = binomial * (height - node) // (node + 1)

    assert x_moment or y_moment
    return x_moment, y_moment


def valuation(value: int, prime: int) -> int:
    """Return v_prime(value), with a large sentinel for value zero."""

    if not value:
        return 10**18
    value = abs(value)
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


@dataclass(frozen=True)
class SupportPrime:
    prime: int
    node: int
    target: bool
    multiplier: bool
    root: bool


def support_at_height(
    height: int,
    apery: list[int],
    apery_lcm: int,
    prime_table: list[bool],
) -> tuple[list[SupportPrime], tuple[int, int]]:
    n = 3 * height + 1
    x_moment, y_moment = cleared_moments(height, apery, apery_lcm)
    support: list[SupportPrime] = []
    for node in range(height + 1):
        prime = n - node
        if not prime_table[prime]:
            continue
        common_valuation = min(
            valuation(x_moment, prime),
            valuation(y_moment, prime),
        )
        in_root = (
            valuation(y_moment - n * x_moment, prime)
            > common_valuation
        )
        in_multiplier = False
        for interpolation_node in range(height + 1):
            apery_valuation = valuation(apery[interpolation_node], prime)
            if not apery_valuation:
                continue
            z_valuation = (
                valuation(
                    y_moment - interpolation_node * x_moment,
                    prime,
                )
                - common_valuation
            )
            if apery_valuation > z_valuation:
                in_multiplier = True
                break
        is_target = apery[node] % prime == 0
        assert not is_target or in_multiplier or in_root
        if in_multiplier or in_root:
            support.append(
                SupportPrime(
                    prime=prime,
                    node=node,
                    target=is_target,
                    multiplier=in_multiplier,
                    root=in_root,
                )
            )
    return support, (x_moment, y_moment)


def log_weight(support: list[SupportPrime], predicate) -> float:
    return sum(log(item.prime) for item in support if predicate(item))


def format_support(support: list[SupportPrime]) -> str:
    fields = []
    for item in support:
        tags = (
            ("T" if item.target else "G")
            + ("M" if item.multiplier else "")
            + ("R" if item.root else "")
        )
        fields.append(f"{item.prime}@{item.node}:{tags}")
    return ",".join(fields) if fields else "-"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("maximum_height", nargs="?", type=int, default=160)
    parser.add_argument(
        "--all",
        action="store_true",
        help="print also heights with empty candidate-window support",
    )
    args = parser.parse_args()
    if args.maximum_height < 2:
        raise SystemExit("maximum_height must be at least 2")

    apery = apery_values(args.maximum_height + 1)
    lcms = prefix_lcms(apery)
    prime_table = prime_sieve(3 * args.maximum_height + 1)

    maximum_rate = (-1.0, -1)
    maximum_pollution_rate = (-1.0, -1)
    nonempty = 0
    for height in range(2, args.maximum_height + 1):
        support, (x_moment, y_moment) = support_at_height(
            height,
            apery,
            lcms[height],
            prime_table,
        )
        if support:
            nonempty += 1
        total_weight = log_weight(support, lambda _: True)
        pollution_weight = log_weight(support, lambda item: not item.target)
        if total_weight / height > maximum_rate[0]:
            maximum_rate = (total_weight / height, height)
        if pollution_weight / height > maximum_pollution_rate[0]:
            maximum_pollution_rate = (pollution_weight / height, height)
        if args.all or support:
            target_weight = log_weight(support, lambda item: item.target)
            multiplier_weight = log_weight(support, lambda item: item.multiplier)
            root_weight = log_weight(support, lambda item: item.root)
            print(
                f"H={height:4d} "
                f"logW/H={total_weight / height:.8f} "
                f"target/H={target_weight / height:.8f} "
                f"pollution/H={pollution_weight / height:.8f} "
                f"multiplier/H={multiplier_weight / height:.8f} "
                f"root/H={root_weight / height:.8f} "
                f"bits(X,Y)=({abs(x_moment).bit_length()},"
                f"{abs(y_moment).bit_length()}) "
                f"support={format_support(support)}"
            )

    print(
        f"scanned=2..{args.maximum_height} nonempty={nonempty} "
        f"max_logW/H={maximum_rate[0]:.8f}@H={maximum_rate[1]} "
        f"max_pollution/H={maximum_pollution_rate[0]:.8f}"
        f"@H={maximum_pollution_rate[1]}"
    )


if __name__ == "__main__":
    main()
