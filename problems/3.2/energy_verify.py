#!/usr/bin/env python3
"""Verify the finite checks in ``energy_result.tex``.

For every prime ``p <= 5000`` this script uses the continuant recurrence
pointwise over F_p.  It computes the roots of N_h for
``2 <= h <= floor(sqrt(p))`` and then checks the exact off-boundary
column/energy identity, the boundary inequalities, and the amplification
inequality.  No numerical result is used as a proof of a uniform estimate.

The implementation is pure Python and does not construct resultants.
"""

from __future__ import annotations

from dataclasses import astuple, dataclass
from hashlib import sha256
from math import isqrt
import sys


PRIME_LIMIT = 5000
EXPECTED_DIGEST = "0e134a06d7d253396cbbbbb46e32d991a37fb30bf1fc8e7f11ac5e39de6eaeca"


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for divisor in range(2, isqrt(limit) + 1):
        if not sieve[divisor]:
            continue
        start = divisor * divisor
        sieve[start : limit + 1 : divisor] = b"\x00" * (
            (limit - start) // divisor + 1
        )
    return [number for number, flag in enumerate(sieve) if flag]


def apery_coefficient(argument: int) -> int:
    return (2 * argument + 1) * (
        17 * argument * argument + 17 * argument + 5
    )


def apery_values(max_index: int) -> list[int]:
    values = [1, 5]
    for index in range(1, max_index):
        numerator = (
            apery_coefficient(index) * values[index]
            - index**3 * values[index - 1]
        )
        denominator = (index + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values[: max_index + 1]


def global_polluted_cuts(prime: int) -> list[int]:
    """Return m with 2 <= m < p and b_{m-1}=0 modulo p."""

    if prime <= 2:
        return []
    previous = 1 % prime
    current = 5 % prime
    cuts = [2] if current == 0 and 2 < prime else []
    for index in range(1, prime - 2):
        numerator = (
            apery_coefficient(index) * current - index**3 * previous
        ) % prime
        denominator = (index + 1) ** 3 % prime
        following = numerator * pow(denominator, -1, prime) % prime
        cut = index + 2
        if following == 0:
            cuts.append(cut)
        previous, current = current, following
    return cuts


def root_data(prime: int, height: int) -> tuple[list[list[int]], list[bytearray], list[int]]:
    """Return root lists, root masks, and column counts for N_2,...,N_H."""

    roots = [[] for _ in range(height + 1)]
    masks = [bytearray(prime) for _ in range(height + 1)]
    column_counts = [0] * prime
    if height < 2:
        return roots, masks, column_counts

    p_values: list[int] = []
    sixth_powers: list[int] = []
    for residue in range(prime):
        square = residue * residue % prime
        p_values.append(
            ((2 * residue + 1) * (17 * square + 17 * residue + 5)) % prime
        )
        sixth_powers.append(pow(residue, 6, prime))

    previous = [0] * prime  # N_0
    current = [1] * prime   # N_1
    for index in range(1, height):
        # At this step q=x+index.  Rotation avoids a modulus operation in
        # every entry of the list comprehension.
        p_shift = p_values[index:] + p_values[:index]
        sixth_shift = sixth_powers[index:] + sixth_powers[:index]
        following = [
            (coefficient * value - sixth * old_value) % prime
            for coefficient, value, sixth, old_value in zip(
                p_shift, current, sixth_shift, previous
            )
        ]
        level = index + 1
        level_roots = [x for x, value in enumerate(following) if value == 0]
        if prime >= 7:
            # Lemma nonvanish and the nominal degree give this independent
            # row check throughout the present range level < p.
            assert len(level_roots) <= 3 * (level - 1)
        mask = masks[level]
        for x in level_roots:
            mask[x] = 1
            column_counts[x] += 1
        roots[level] = level_roots
        previous, current = current, following

    return roots, masks, column_counts


@dataclass(frozen=True)
class PrimeRecord:
    prime: int
    height: int
    root_mass: int
    energy: int
    energy_nonboundary: int
    energy_boundary: int
    column_pairs_nonboundary: int
    column_pairs_all: int
    boundary_mass: int
    low_fiber_mass: int
    window_polluted_count: int
    polluted_mass: int
    polluted_pairs: int
    max_unpolluted_column: int


def prime_record(prime: int, apery: list[int]) -> PrimeRecord:
    height = isqrt(prime)
    roots, masks, column_counts = root_data(prime, height)

    boundary = bytearray(prime)
    for cut in range(2, height + 1):
        boundary[(-cut) % prime] = 1

    all_polluted_indices = global_polluted_cuts(prime)
    polluted_indices = [cut for cut in all_polluted_indices if cut <= height]
    assert polluted_indices == [
        cut
        for cut in range(2, height + 1)
        if apery[cut - 1] % prime == 0
    ]
    polluted = bytearray(prime)
    for cut in polluted_indices:
        polluted[(-cut) % prime] = 1
    globally_polluted = bytearray(prime)
    for cut in all_polluted_indices:
        globally_polluted[(-cut) % prime] = 1

    energy_nonboundary = 0
    energy_boundary = 0
    energy_by_witness = [0] * prime
    for first_gap in range(2, max(2, height - 1)):
        for second_gap in range(2, height - first_gap + 1):
            second_mask = masks[second_gap]
            for witness in roots[first_gap]:
                if not second_mask[(witness + first_gap) % prime]:
                    continue
                energy_by_witness[witness] += 1
                if boundary[witness]:
                    energy_boundary += 1
                else:
                    energy_nonboundary += 1

    column_pairs_nonboundary = sum(
        count * (count - 1) // 2
        for x, count in enumerate(column_counts)
        if not boundary[x]
    )
    column_pairs_all = sum(
        count * (count - 1) // 2 for count in column_counts
    )
    root_mass = sum(column_counts)
    boundary_mass = sum(
        count for x, count in enumerate(column_counts) if boundary[x]
    )
    polluted_mass = sum(
        count for x, count in enumerate(column_counts) if polluted[x]
    )
    polluted_pairs = sum(
        count * (count - 1) // 2
        for x, count in enumerate(column_counts)
        if polluted[x]
    )
    max_unpolluted_column = max(
        (
            count
            for x, count in enumerate(column_counts)
            if not globally_polluted[x]
        ),
        default=0,
    )

    threshold = isqrt(height)
    if threshold * threshold < height:
        threshold += 1
    low_fiber_mass = sum(
        count
        for x, count in enumerate(column_counts)
        if not boundary[x] and 0 < count < threshold
    )

    energy = energy_nonboundary + energy_boundary

    # Exact/sandwich relations proved in energy_result.tex.
    for x, count in enumerate(column_counts):
        column_pairs = count * (count - 1) // 2
        assert energy_by_witness[x] <= column_pairs
        if not boundary[x]:
            assert energy_by_witness[x] == column_pairs
    assert sum(energy_by_witness) == energy
    assert energy_nonboundary == column_pairs_nonboundary
    assert energy <= column_pairs_all
    assert root_mass <= 3 * height * (height - 1) // 2

    nonboundary_mass = root_mass - boundary_mass
    nonboundary_max = max(
        (count for x, count in enumerate(column_counts) if not boundary[x]),
        default=0,
    )
    assert 2 * column_pairs_nonboundary <= max(0, nonboundary_max - 1) * nonboundary_mass

    if threshold >= 2:
        high_fiber_mass = root_mass - boundary_mass - low_fiber_mass
        assert high_fiber_mass >= 0
        # This is the integer form of both the exact off-boundary estimate
        # and its weakened version with total energy.
        assert high_fiber_mass * (threshold - 1) <= 2 * column_pairs_nonboundary
        assert high_fiber_mass * (threshold - 1) <= 2 * energy

    return PrimeRecord(
        prime,
        height,
        root_mass,
        energy,
        energy_nonboundary,
        energy_boundary,
        column_pairs_nonboundary,
        column_pairs_all,
        boundary_mass,
        low_fiber_mass,
        len(polluted_indices),
        polluted_mass,
        polluted_pairs,
        max_unpolluted_column,
    )


def record_line(record: PrimeRecord) -> str:
    return ",".join(str(value) for value in astuple(record))


def check_combinatorial_formulas() -> None:
    # Shallow-strip weighted mass.
    for height in range(4, 80):
        for cutoff in range(2, height // 2 + 1):
            direct = sum(
                min(first, second) - 1
                for first in range(2, height - 1)
                for second in range(2, height - first + 1)
                if min(first, second) <= cutoff
            )
            formula = sum(
                (minimum - 1) * (2 * height - 4 * minimum + 1)
                for minimum in range(2, cutoff + 1)
            )
            assert direct == formula

    # Exact degree weight of the projective-infinity lattice.
    for rank in range(2, 20):
        for height in range(4, 80):
            direct = sum(
                min(first, second) - 1
                for first in range(2, height - 1)
                for second in range(2, height - first + 1)
                if first % rank == 0 and second % rank == 0
            )
            quotient = height // rank
            formula = rank * sum(
                (total * total) // 4 for total in range(2, quotient + 1)
            ) - quotient * (quotient - 1) // 2
            assert direct == formula


def run() -> None:
    if not __debug__:
        raise RuntimeError("run without -O; assertions are verification checks")
    check_combinatorial_formulas()
    primes = primes_up_to(PRIME_LIMIT)
    apery = apery_values(isqrt(PRIME_LIMIT) + 1)
    records = [prime_record(prime, apery) for prime in primes]

    payload = "\n".join(record_line(record) for record in records)
    digest = sha256(payload.encode("ascii")).hexdigest()
    assert digest == EXPECTED_DIGEST
    assert len(records) == 669

    field_sums = tuple(sum(values) for values in zip(*(astuple(r) for r in records)))
    assert field_sums[2:] == (
        43407, 677, 671, 6, 671, 3627, 676, 42731,
        12, 184, 2952, 944,
    )

    by_prime = {record.prime: record for record in records}
    expected_rows = {
        73: (73, 8, 19, 2, 0, 2, 0, 16, 11, 8, 1, 6, 15, 2),
        131: (131, 11, 17, 6, 6, 0, 6, 6, 0, 17, 0, 0, 0, 3),
        653: (653, 25, 38, 16, 16, 0, 16, 16, 0, 38, 0, 0, 0, 4),
        1049: (1049, 32, 46, 2, 2, 0, 2, 2, 0, 46, 0, 0, 0, 1),
        3331: (3331, 57, 187, 0, 0, 0, 0, 1225, 98, 89, 1, 50, 1225, 1),
        4283: (4283, 65, 203, 4, 2, 2, 2, 1084, 93, 110, 1, 47, 1081, 2),
    }
    for prime, expected in expected_rows.items():
        assert astuple(by_prime[prime]) == expected

    assert max(records, key=lambda record: record.root_mass).prime == 4283
    assert max(record.root_mass for record in records) == 203
    assert max(records, key=lambda record: record.energy).prime == 653
    assert max(record.energy for record in records) == 16
    assert max(record.column_pairs_nonboundary for record in records) == 16
    assert max(record.boundary_mass for record in records) == 98
    assert max(record.polluted_mass for record in records) == 50
    assert max(record.polluted_pairs for record in records) == 1225
    assert max(record.max_unpolluted_column for record in records) == 4
    assert max(record.window_polluted_count for record in records) == 1
    assert all(
        record.root_mass
        == record.boundary_mass + record.low_fiber_mass
        for record in records
    )

    # Compare ratios by exact cross multiplication.
    max_r_ratio = records[0]
    for record in records[1:]:
        if (
            record.root_mass * max_r_ratio.height
            > max_r_ratio.root_mass * record.height
        ):
            max_r_ratio = record
    assert max_r_ratio.prime == 3331
    max_e_ratio = records[0]
    for record in records[1:]:
        if (
            record.energy**2 * max_e_ratio.height**3
            > max_e_ratio.energy**2 * record.height**3
        ):
            max_e_ratio = record
    assert max_e_ratio.prime == 131

    assert apery[7] == 584307365 == 3331 * 175415
    assert by_prime[3331].polluted_mass == 57 - 7

    print("ENERGY COLUMN IDENTITIES PASS")
    print(f"primes<={PRIME_LIMIT}: {len(records)}")
    print(f"canonical_sha256={digest}")
    print("max R=203 at p=4283; max R/H=187/57 at p=3331")
    print("max E=16 at p=653; max E/H^(3/2)=6/11^(3/2) at p=131")
    print("max off-boundary column-pair sum=16 at p=653")
    print("max polluted k=50 and C(k,2)=1225 at p=3331")
    print("ENERGY_VERIFY PASS")


if __name__ == "__main__":
    try:
        run()
    except Exception:
        print("ENERGY_VERIFY FAIL", file=sys.stderr)
        raise
