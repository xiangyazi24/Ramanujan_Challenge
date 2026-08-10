#!/usr/bin/env python3
"""Check every shift and sign in bridge reflection reciprocity."""

from __future__ import annotations

import hashlib
import json
import math


EXPECTED_SHA256 = "4faf19b270b89c14f9e6e584fd2631bd79a242fc65704f902e9dfbcd78f98b07"


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, math.isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime :: prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [prime for prime in range(5, limit + 1) if sieve[prime]]


def middle(index: int, prime: int) -> int:
    return (
        34 * index**3 + 51 * index**2 + 27 * index + 5
    ) % prime


def continuants(prime: int) -> list[list[int]]:
    rows: list[list[int]] = []
    for start in range(prime):
        length = prime - start
        values = [0] * (length + 1)
        values[0] = 0
        if length >= 1:
            values[1] = 1
        if length >= 2:
            values[2] = middle(start + 1, prime)
        for gap in range(2, length):
            point = start + gap
            values[gap + 1] = (
                middle(point, prime) * values[gap]
                - pow(point, 6, prime) * values[gap - 1]
            ) % prime
        rows.append(values)
    return rows


def solution_rows(prime: int) -> tuple[list[int], list[int]]:
    apery = [0] * prime
    companion = [0] * prime
    apery[0], apery[1] = 1, 5 % prime
    companion[0], companion[1] = 0, 1
    for index in range(1, prime - 1):
        apery[index + 1] = (
            middle(index, prime) * apery[index]
            - pow(index, 6, prime) * apery[index - 1]
        ) % prime
        companion[index + 1] = (
            middle(index, prime) * companion[index]
            - pow(index, 6, prime) * companion[index - 1]
        ) % prime
    return apery, companion


def audit_prime(prime: int) -> dict[str, int]:
    gaps = continuants(prime)
    apery, companion = solution_rows(prime)
    factorial = 1
    rho = []
    for index in range(prime):
        if index:
            factorial = factorial * index % prime
        sign = -1 if (index + 1) % 2 else 1
        rho.append(sign * pow(factorial, -6, prime) % prime)

    for index in range(prime):
        reflected = prime - 1 - index
        assert apery[reflected] == rho[index] * apery[index] % prime
        assert companion[reflected] == rho[index] * companion[index] % prime

    return_count = 0
    for start in range(prime - 1):
        for gap in range(2, prime - start):
            if gaps[start][gap] != 0:
                continue
            return_count += 1
            end = start + gap
            multiplier = (
                -pow(start + 1, 6, prime) * gaps[start + 1][gap - 1]
            ) % prime
            assert multiplier
            assert apery[end] == multiplier * apery[start] % prime
            assert companion[end] == multiplier * companion[start] % prime

            mirror_start = prime - 1 - end
            assert gaps[mirror_start][gap] == 0
            mirror_multiplier = (
                -pow(mirror_start + 1, 6, prime)
                * gaps[mirror_start + 1][gap - 1]
            ) % prime
            expected = rho[start] * pow(rho[end], -1, prime) % prime
            assert multiplier * mirror_multiplier % prime == expected

    return {"prime": prime, "returns": return_count}


def audit_known_record() -> dict[str, int]:
    prime = 1297
    gaps = continuants(prime)
    apery, companion = solution_rows(prime)
    start = 360
    offsets = (5, 25, 35, 541)
    assert all(gaps[start][gap] == 0 for gap in offsets)
    state = apery[start] * pow(companion[start], -1, prime) % prime
    assert state == 454

    factorial = 1
    rho = []
    for index in range(prime):
        if index:
            factorial = factorial * index % prime
        sign = -1 if (index + 1) % 2 else 1
        rho.append(sign * pow(factorial, -6, prime) % prime)

    for gap in offsets:
        end = start + gap
        multiplier = (
            -pow(start + 1, 6, prime) * gaps[start + 1][gap - 1]
        ) % prime
        assert multiplier
        assert apery[end] == multiplier * apery[start] % prime
        assert companion[end] == multiplier * companion[start] % prime

        mirror_start = prime - 1 - end
        assert gaps[mirror_start][gap] == 0
        mirror_multiplier = (
            -pow(mirror_start + 1, 6, prime)
            * gaps[mirror_start + 1][gap - 1]
        ) % prime
        expected = rho[start] * pow(rho[end], -1, prime) % prime
        assert multiplier * mirror_multiplier % prime == expected

    return {
        "prime": prime,
        "start": start,
        "state": state,
        "bridges": len(offsets),
    }


def main() -> None:
    rows = [audit_prime(prime) for prime in primes_up_to(101)]
    known = audit_known_record()
    payload = {"rows": rows, "known": known}
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    assert digest == EXPECTED_SHA256
    print(
        "BRIDGE_REFLECTION_MULTIPLIER_VERIFY"
        f" primes={len(rows)} returns={sum(row['returns'] for row in rows)}"
        " known_prime=1297 known_start=360"
    )
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
