#!/usr/bin/env python3
"""Exact audit of the shortest Newton carrier at the second p-adic layer.

For a top-half target ``n=p+r`` put ``M=n-1=p+s`` with ``s=r-1`` and
``Y_d=C_M(d)``.  The two length-one Newton carriers whose common node is
``p-1`` satisfy

    G_{p-2,1} - G_{p-1,1} = p * Delta^2 Y_{p-2}.

Since targetness is ``p | Y_{p-1}``, division by the single universal
node factor gives

    (G_{p-2,1} - G_{p-1,1}) / p
        = Y_{p-2} + Y_p                 (mod p).

Cartier reduction gives ``Y_p = 40*b_s (mod p)``.  This script audits
the resulting distinguished second-layer scalar for every top-half
target through a configurable bound.  Shells are evaluated modulo p
from the exact one-fold coefficient formula, with binomial coefficients
computed by Lucas' theorem; no coefficient symmetry is assumed.
"""

from __future__ import annotations

from argparse import ArgumentParser
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


def apery_values(limit: int) -> list[int]:
    if limit == 0:
        return [1]
    values = [1, 5]
    for m in range(1, limit):
        numerator = (
            (34 * m**3 + 51 * m**2 + 27 * m + 5) * values[m]
            - m**3 * values[m - 1]
        )
        denominator = (m + 1) ** 3
        value, remainder = divmod(numerator, denominator)
        assert remainder == 0
        values.append(value)
    return values[: limit + 1]


class LucasBinomial:
    """Binomial coefficients modulo one prime."""

    def __init__(self, prime: int) -> None:
        self.prime = prime
        factorial = [1] * prime
        for index in range(1, prime):
            factorial[index] = factorial[index - 1] * index % prime
        inverse_factorial = [1] * prime
        inverse_factorial[-1] = pow(factorial[-1], prime - 2, prime)
        for index in range(prime - 1, 0, -1):
            inverse_factorial[index - 1] = (
                inverse_factorial[index] * index % prime
            )
        self.factorial = factorial
        self.inverse_factorial = inverse_factorial

    def small(self, n: int, k: int) -> int:
        if not 0 <= k <= n < self.prime:
            return 0
        return (
            self.factorial[n]
            * self.inverse_factorial[k]
            * self.inverse_factorial[n - k]
            % self.prime
        )

    def __call__(self, n: int, k: int) -> int:
        if k < 0 or k > n:
            return 0
        out = 1
        while n or k:
            n_digit = n % self.prime
            k_digit = k % self.prime
            if k_digit > n_digit:
                return 0
            out = out * self.small(n_digit, k_digit) % self.prime
            n //= self.prime
            k //= self.prime
        return out


def shell_mod(moment: int, node: int, prime: int) -> int:
    """Return the exact shell ``C_moment(node)`` modulo ``prime``."""

    choose = LucasBinomial(prime)
    quotient = moment // node
    out = 0
    for t in range(moment + 1):
        outer = choose(moment, t)
        if not outer:
            continue
        base = moment - t
        x_packet = sum(
            choose(moment, base + node * u)
            for u in range(-quotient, quotient + 1)
        ) % prime
        yz_packet = sum(
            choose(2 * moment - t, base + node * v)
            for v in range(-quotient, quotient + 1)
        ) % prime
        out = (out + outer * x_packet * yz_packet**2) % prime
    return out


def shell_mod_direct(moment: int, node: int, prime: int) -> int:
    """Slow independent evaluation used only by the self-test."""

    quotient = moment // node
    out = 0
    for t in range(moment + 1):
        base = moment - t
        x_packet = sum(
            comb(moment, base + node * u)
            for u in range(-quotient, quotient + 1)
            if 0 <= base + node * u <= moment
        )
        yz_packet = sum(
            comb(2 * moment - t, base + node * v)
            for v in range(-quotient, quotient + 1)
            if 0 <= base + node * v <= 2 * moment - t
        )
        out += comb(moment, t) * x_packet * yz_packet**2
    return out % prime


def self_test() -> int:
    checks = 0
    for prime in primes_upto(31):
        if prime < 5:
            continue
        for s in range(prime - 2):
            moment = prime + s
            for node in (prime - 2, prime - 1, prime):
                assert shell_mod(moment, node, prime) == shell_mod_direct(
                    moment, node, prime
                )
                checks += 1
    return checks


def target_records(limit: int) -> list[tuple[int, int, int]]:
    values = apery_values(limit)
    records = []
    for prime in primes_upto(limit):
        if prime <= 5:
            continue
        for r in range(1, min(prime - 1, limit - prime) + 1):
            if values[r] % prime == 0:
                records.append((prime + r, prime, r))
    return records


def audit(limit: int) -> tuple[int, list[tuple[int, ...]]]:
    values = apery_values(limit)
    records = target_records(limit)
    exceptional = []
    for n, prime, r in records:
        s = r - 1
        moment = n - 1
        assert moment == prime + s

        left_shell = shell_mod(moment, prime - 2, prime)
        target_shell = shell_mod(moment, prime - 1, prime)
        right_shell = shell_mod(moment, prime, prime)

        assert target_shell == values[r] % prime == 0
        assert right_shell == 40 * values[s] % prime

        divided_difference = (left_shell + right_shell) % prime
        if divided_difference == 0:
            exceptional.append(
                (
                    n,
                    prime,
                    r,
                    left_shell,
                    right_shell,
                    divided_difference,
                )
            )
    return len(records), exceptional


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("--limit", type=int, default=5000)
    args = parser.parse_args()

    checks = self_test()
    target_count, exceptional = audit(args.limit)
    print("Q32_NEWTON_SECOND_LAYER_AUDIT=PASS")
    print("SHELL_SELF_TESTS", checks)
    print("LIMIT", args.limit)
    print("TOP_HALF_TARGETS_P_GT_5", target_count)
    print("VANISHING_DIVIDED_DIFFERENCES", len(exceptional))
    print("EXCEPTIONS", exceptional)


if __name__ == "__main__":
    main()
