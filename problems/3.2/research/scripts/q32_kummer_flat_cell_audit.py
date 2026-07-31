#!/usr/bin/env python3
"""Audit the exact Kummer flat-cell mechanism for Apéry shell carriers.

For

    C_M(e) = sum_t binom(M,t)
             sum_u binom(M,M-t+e*u)
             (sum_v binom(2*M-t,M-t+e*v))^2,

the first cell is ``M/2 < e <= M``, so only ``u,v in {-1,0,1}``
can occur.  Let ``q`` be prime and write

    M = alpha*q + r,       e = (alpha-1)*q + sigma.

If

    r < sigma < q-r,

Lucas' theorem kills every noncentral packet term individually.
Consequently

    C_M(e) = b_alpha*b_r (mod q).

The script verifies the termwise statement, its two sharp witnesses,
the moving endpoint in the fixed 3-by-5 rectangle, and the reason why
the natural full-margin family exits the flat first cell.

This is an audit of coefficient-forced nuisance primes.  Failure of the
flat inequalities does *not* imply nondivisibility of an evaluated
cross-minor: targetness or an arithmetic boundary cancellation can
still make the minor vanish.
"""

from __future__ import annotations

from argparse import ArgumentParser
from dataclasses import dataclass
from math import comb, isqrt


def primes_upto(limit: int) -> list[int]:
    mark = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        mark[0] = 0
    if limit >= 1:
        mark[1] = 0
    for prime in range(2, isqrt(limit) + 1):
        if mark[prime]:
            start = prime * prime
            mark[start : limit + 1 : prime] = b"\x00" * (
                (limit - start) // prime + 1
            )
    return [prime for prime in range(2, limit + 1) if mark[prime]]


def lucas_binomial(upper: int, lower: int, prime: int) -> int:
    """Return ``binom(upper,lower) modulo prime`` by Lucas digits."""

    if lower < 0 or lower > upper:
        return 0
    out = 1
    while upper or lower:
        upper_digit = upper % prime
        lower_digit = lower % prime
        if lower_digit > upper_digit:
            return 0
        out = out * comb(upper_digit, lower_digit) % prime
        upper //= prime
        lower //= prime
    return out


def apery_mod(index: int, prime: int) -> int:
    """Return the Apéry number ``b_index`` modulo ``prime``."""

    return sum(
        lucas_binomial(index, k, prime) ** 2
        * lucas_binomial(index + k, k, prime) ** 2
        for k in range(index + 1)
    ) % prime


def shell_mod(moment: int, node: int, prime: int) -> int:
    """Return the exact shell ``C_moment(node)`` modulo ``prime``."""

    quotient = moment // node
    out = 0
    for t in range(moment + 1):
        outer = lucas_binomial(moment, t, prime)
        if not outer:
            continue
        base = moment - t
        x_packet = sum(
            lucas_binomial(moment, base + node * u, prime)
            for u in range(-quotient, quotient + 1)
        ) % prime
        yz_packet = sum(
            lucas_binomial(2 * moment - t, base + node * v, prime)
            for v in range(-quotient, quotient + 1)
        ) % prime
        out = (out + outer * x_packet * yz_packet**2) % prime
    return out


def is_first_cell(moment: int, node: int) -> bool:
    return moment < 2 * node <= 2 * moment


def flat_parameters(
    moment: int, node: int, prime: int, alpha: int
) -> tuple[int, int] | None:
    """Return ``(r,sigma)`` in the relevant cell, or ``None``."""

    r = moment - alpha * prime
    cell = alpha - 1
    sigma = node - cell * prime
    if not (0 <= r < prime):
        return None
    if not (0 < sigma < prime):
        return None
    if node // prime != cell:
        return None
    if not is_first_cell(moment, node):
        return None
    return r, sigma


def is_termwise_flat(
    moment: int, node: int, prime: int, alpha: int
) -> bool:
    parameters = flat_parameters(moment, node, prime, alpha)
    if parameters is None:
        return False
    r, sigma = parameters
    return r < sigma < prime - r


def audit_one_flat_shell(
    moment: int, node: int, prime: int, alpha: int
) -> int:
    """Check every Kummer carry and the resulting shell residue."""

    parameters = flat_parameters(moment, node, prime, alpha)
    assert parameters is not None
    r, sigma = parameters
    assert r < sigma < prime - r
    assert 2 * r < prime

    for t in range(moment + 1):
        outer = lucas_binomial(moment, t, prime)
        if not outer:
            continue

        # By symmetry, the two noncentral terms in the first packet are
        # binom(M,t+e) and binom(M,t-e).
        assert lucas_binomial(moment, t + node, prime) == 0
        assert lucas_binomial(moment, t - node, prime) == 0

        # The two noncentral terms in the second packet have lower
        # arguments M-t+e and M-t-e.
        upper = 2 * moment - t
        center = moment - t
        assert lucas_binomial(upper, center + node, prime) == 0
        assert lucas_binomial(upper, center - node, prime) == 0

    expected = apery_mod(alpha, prime) * apery_mod(r, prime) % prime
    assert shell_mod(moment, node, prime) == expected
    return expected


def audit_sharp_witness(
    moment: int, node: int, prime: int, alpha: int
) -> str:
    """Exhibit a noncentral Lucas unit when a flat inequality fails."""

    parameters = flat_parameters(moment, node, prime, alpha)
    assert parameters is not None
    r, sigma = parameters
    if sigma <= r:
        assert lucas_binomial(moment, 0, prime) != 0
        assert lucas_binomial(moment, node, prime) != 0
        return "left"
    if prime - sigma <= r:
        t = alpha * prime
        assert t <= moment
        assert lucas_binomial(moment, t, prime) != 0
        assert lucas_binomial(moment, t - node, prime) != 0
        return "right"
    raise AssertionError("the shell is flat, so no sharp witness exists")


@dataclass(frozen=True)
class Rectangle:
    n: int
    d: int
    length: int
    s_minus: int
    s_plus: int
    t_minus: int
    t_plus: int

    @property
    def node_minimum(self) -> int:
        # Z_x uses shells at x-1 and x+1.
        return self.d + self.s_minus - 1

    @property
    def node_maximum(self) -> int:
        # The right carrier is G_{D+1,K}; its last sequence index is
        # D+K+1, and Z at that index reaches the shell node D+K+2.
        return (
            self.d
            + self.s_plus
            + self.length
            + self.t_plus
            + 2
        )

    @property
    def moments(self) -> range:
        return range(self.n - 1, self.n + 5)

    def lies_strictly_in_first_cell(self) -> bool:
        return (
            2 * self.node_minimum > self.n + 4
            and self.node_maximum <= self.n - 1
        )


def rectangle_is_flat(
    rectangle: Rectangle, prime: int, alpha: int
) -> bool:
    """Endpoint test for flatness of every shell in the rectangle."""

    if not rectangle.lies_strictly_in_first_cell():
        return False
    if not (
        alpha * prime <= rectangle.n - 1
        and rectangle.n + 4 < (alpha + 1) * prime
    ):
        return False

    cell = alpha - 1
    if not (
        cell * prime < rectangle.node_minimum
        <= rectangle.node_maximum
        < (cell + 1) * prime
    ):
        return False

    largest_remainder = rectangle.n + 4 - alpha * prime
    sigma_minimum = rectangle.node_minimum - cell * prime
    sigma_maximum = rectangle.node_maximum - cell * prime
    return (
        largest_remainder < sigma_minimum
        and sigma_maximum < prime - largest_remainder
    )


def offset_bounds(rectangle: Rectangle) -> tuple[int, int]:
    """Return the exact alpha=1 and alpha=2 flat offset bounds.

    In the first branch ``q=n-k``.  In the second branch
    ``2*q=n-k``.  Stable-quotient and primality conditions are checked
    separately; these are the two sharp flat inequalities in integer
    form.
    """

    upper_bound = (
        rectangle.n - rectangle.node_maximum - 5
    ) // 2
    alpha_one = min(rectangle.node_minimum - 5, upper_bound)
    # q >= n+5-e_min and 2*q=n-k give
    # k <= 2*e_min-n-10.  The parity condition on k prevents replacing
    # this by the weaker off-by-one bound with -9.
    alpha_two = min(
        2 * rectangle.node_minimum - rectangle.n - 10,
        upper_bound,
    )
    return alpha_one, alpha_two


def forced_primes(
    rectangle: Rectangle, primes: set[int]
) -> tuple[list[int], list[int]]:
    alpha_one_bound, alpha_two_bound = offset_bounds(rectangle)
    near_n = [
        rectangle.n - k
        for k in range(1, max(0, alpha_one_bound) + 1)
        if rectangle.n - k in primes
        and rectangle_is_flat(rectangle, rectangle.n - k, 1)
    ]
    half_scale = [
        (rectangle.n - k) // 2
        for k in range(1, max(0, alpha_two_bound) + 1)
        if (rectangle.n - k) % 2 == 0
        and (rectangle.n - k) // 2 in primes
        and rectangle_is_flat(
            rectangle, (rectangle.n - k) // 2, 2
        )
    ]
    return near_n, half_scale


def audit_fixed_rectangle(exhaustive: bool) -> list[tuple[object, ...]]:
    """Reproduce the six exact moving-tail endpoints."""

    sample_indices = (236, 272, 300, 311, 350, 400)
    expected = {
        236: (13, 223),
        272: (16, 256),
        300: (17, 283),
        311: (18, 293),
        350: (21, 329),
        400: (25, 375),
    }
    prime_set = set(primes_upto(max(sample_indices)))
    rows = []

    for n in sample_indices:
        rectangle = Rectangle(
            n=n,
            d=(13 * n) // 20 - 2,
            length=n // 5 + 1,
            # Three translations and five neighboring lengths.
            s_minus=0,
            s_plus=2,
            t_minus=-2,
            t_plus=2,
        )
        assert rectangle.lies_strictly_in_first_cell()
        alpha_one_bound, alpha_two_bound = offset_bounds(rectangle)
        expected_bound, expected_endpoint = expected[n]
        assert alpha_one_bound == alpha_two_bound == expected_bound
        assert n - alpha_one_bound == expected_endpoint

        near_n, half_scale = forced_primes(rectangle, prime_set)
        for alpha, forced in ((1, near_n), (2, half_scale)):
            for prime in forced:
                nodes = (
                    range(
                        rectangle.node_minimum,
                        rectangle.node_maximum + 1,
                    )
                    if exhaustive
                    else (
                        rectangle.node_minimum,
                        (
                            rectangle.node_minimum
                            + rectangle.node_maximum
                        )
                        // 2,
                        rectangle.node_maximum,
                    )
                )
                for moment in rectangle.moments:
                    for node in nodes:
                        audit_one_flat_shell(
                            moment, node, prime, alpha
                        )

        rows.append(
            (
                n,
                rectangle.d,
                rectangle.length,
                rectangle.node_minimum,
                rectangle.node_maximum,
                alpha_one_bound,
                n - alpha_one_bound,
                tuple(near_n),
                tuple(half_scale),
            )
        )
    return rows


def full_margin_rectangle(n: int, D: int, N: int) -> tuple[int, Rectangle]:
    m = min(
        D - (n - 1) // 2,
        (n - 1) - D - N + 2,
    )
    assert m >= 2
    d0 = D - m + 1
    length = N + m - 2
    # This fixed-length family has translations 0,...,m-2 and no
    # independent length shifts.
    rectangle = Rectangle(
        n=n,
        d=d0,
        length=length,
        s_minus=0,
        s_plus=m - 2,
        t_minus=0,
        t_plus=0,
    )
    assert rectangle.node_minimum == D - m
    assert rectangle.node_maximum == D + N + m - 1
    return m, rectangle


def audit_full_margin_boundaries() -> list[tuple[object, ...]]:
    cases = (
        (200, 128, 63),
        (272, 180, 63),
        (300, 180, 57),
        (321, 168, 53),
    )
    rows = []
    for n, D, N in cases:
        m, rectangle = full_margin_rectangle(n, D, N)
        assert not rectangle.lies_strictly_in_first_cell()
        if 2 * rectangle.node_minimum <= n + 4:
            boundary = "half-cell"
        else:
            assert rectangle.node_maximum > n - 1
            boundary = "top-cell"
        rows.append(
            (
                n,
                D,
                N,
                m,
                rectangle.node_minimum,
                rectangle.node_maximum,
                boundary,
            )
        )
    return rows


def audit_small_flat_lemma(limit: int) -> tuple[int, int]:
    """Exhaust the flat lemma and both sharp failure witnesses."""

    flat_checks = 0
    witness_checks = 0
    for prime in primes_upto(limit):
        if prime < 5:
            continue
        for alpha in (1, 2):
            cell = alpha - 1
            for r in range((prime - 1) // 2 + 1):
                moment = alpha * prime + r
                for sigma in range(1, prime):
                    node = cell * prime + sigma
                    if not is_first_cell(moment, node):
                        continue
                    if r < sigma < prime - r:
                        audit_one_flat_shell(
                            moment, node, prime, alpha
                        )
                        flat_checks += 1
                    else:
                        audit_sharp_witness(
                            moment, node, prime, alpha
                        )
                        witness_checks += 1
    return flat_checks, witness_checks


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument(
        "--small-prime-limit",
        type=int,
        default=23,
    )
    parser.add_argument(
        "--exhaustive-rectangle",
        action="store_true",
        help="check every shell node in the six moving-tail rows",
    )
    args = parser.parse_args()

    flat_counts = audit_small_flat_lemma(args.small_prime_limit)
    fixed_rows = audit_fixed_rectangle(args.exhaustive_rectangle)
    margin_rows = audit_full_margin_boundaries()

    print("FLAT_AND_SHARP_CHECKS", flat_counts)
    print("FIXED_RECTANGLE_ROWS")
    for row in fixed_rows:
        print(row)
    print("FULL_MARGIN_BOUNDARIES")
    for row in margin_rows:
        print(row)
    print("PASS: exact Kummer flat-cell audit")


if __name__ == "__main__":
    main()
