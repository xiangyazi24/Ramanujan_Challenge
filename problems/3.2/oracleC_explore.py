#!/usr/bin/env python3
"""Finite exploration for the Oracle C marked-coordinate audit.

The computations deliberately distinguish five different notions.

* The Apéry number is the constant term of the power of one fixed Laurent
  polynomial in three variables.  Torus orthogonality therefore gives an
  exact marked coordinate whose ambient dimension, exponent box, divisor,
  coefficient rank, and tame ramification data are independent of ``p,j``.

* ``B_j = (j!)^3 b_j`` is the (1,1)-entry of a length-``j`` transfer
  product.  The state dimension is two, but the word length grows with ``j``.
* The Apéry recurrence has order two with degree-three polynomial
  coefficients.  This is a difference module, not a trace sheaf on a
  ``j``-line.
* The raw finite Mellin marker has a canonical group-algebra support.  Its
  support and its constant-coefficient cyclic recurrence order grow with
  ``p``.
* A single nontrivial Kummer sheaf on the original ``t``-line has rank one and
  geometric conductor three (rank + two tame singularities).  No sheaf on a
  ``j``-line whose trace is ``b_j`` has been specified, so its conductor is
  not a computable invariant.

All arithmetic is exact.  Numerical output is reconnaissance, never an input
to the propositions in ``oracleC_result.tex``.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from math import gcd, isqrt
from typing import Iterable, Sequence


DEFAULT_LIMIT = 2_000
SAMPLE_PRIMES = (5, 7, 11, 23, 31, 181, 827, 1999)
TORIC_SAMPLE_PRIMES = (5, 7, 11, 13, 17, 19, 23, 29, 31)


@dataclass(frozen=True)
class PrimeRecord:
    prime: int
    zero_count: int
    transfer_cases: int
    bm_b: int
    bm_B: int
    interpolation_degree_b: int
    interpolation_degree_B: int
    hasse_at_one: int
    mellin_support: int
    support_lower_bound: int
    epsilon: int


def primes_up_to(limit: int) -> list[int]:
    """Return all primes in ``[5, limit]``."""

    if limit < 5:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for divisor in range(2, isqrt(limit) + 1):
        if sieve[divisor]:
            start = divisor * divisor
            sieve[start : limit + 1 : divisor] = b"\x00" * (
                (limit - start) // divisor + 1
            )
    return [prime for prime in range(5, limit + 1) if sieve[prime]]


def prime_divisors(number: int) -> list[int]:
    divisors: list[int] = []
    candidate = 2
    while candidate * candidate <= number:
        if number % candidate == 0:
            divisors.append(candidate)
            while number % candidate == 0:
                number //= candidate
        candidate += 1
    if number > 1:
        divisors.append(number)
    return divisors


def primitive_root(prime: int) -> int:
    """Return the least primitive root of the odd prime ``prime``."""

    divisors = prime_divisors(prime - 1)
    for candidate in range(2, prime):
        if all(pow(candidate, (prime - 1) // divisor, prime) != 1 for divisor in divisors):
            return candidate
    raise ArithmeticError(f"no primitive root found modulo {prime}")


def apery_integers(max_index: int) -> list[int]:
    """Return the integer Apéry numbers ``b_0,...,b_max_index``."""

    values = [1]
    if max_index == 0:
        return values
    values.append(5)
    for index in range(1, max_index):
        numerator = (
            polynomial_P(index) * values[index]
            - index**3 * values[index - 1]
        )
        denominator = (index + 1) ** 3
        if numerator % denominator:
            raise ArithmeticError("Apéry recurrence failed to divide")
        values.append(numerator // denominator)
    return values


def polynomial_P(index: int) -> int:
    return 34 * index**3 + 51 * index**2 + 27 * index + 5


def apery_mod_prime(prime: int) -> list[int]:
    """Return ``b_0,...,b_{p-1}`` in ``F_p`` from the divided recurrence."""

    if prime < 5:
        raise ValueError("the Oracle C scan starts at p=5")
    values = [1, 5 % prime]
    for index in range(1, prime - 1):
        denominator = pow(index + 1, 3, prime)
        numerator = (
            polynomial_P(index) * values[index]
            - index**3 * values[index - 1]
        ) % prime
        values.append(numerator * pow(denominator, -1, prime) % prime)
    return values


def renormalized_values(values: Sequence[int], prime: int) -> list[int]:
    """Return ``B_j=(j!)^3 b_j`` in ``F_p``."""

    factorial = 1
    output: list[int] = []
    for index, value in enumerate(values):
        if index:
            factorial = factorial * index % prime
        output.append(pow(factorial, 3, prime) * value % prime)
    return output


Matrix = tuple[tuple[int, int], tuple[int, int]]


def matrix_multiply(left: Matrix, right: Matrix, prime: int) -> Matrix:
    return (
        (
            (left[0][0] * right[0][0] + left[0][1] * right[1][0])
            % prime,
            (left[0][0] * right[0][1] + left[0][1] * right[1][1])
            % prime,
        ),
        (
            (left[1][0] * right[0][0] + left[1][1] * right[1][0])
            % prime,
            (left[1][0] * right[0][1] + left[1][1] * right[1][1])
            % prime,
        ),
    )


def transfer_matrix(index: int, prime: int) -> Matrix:
    return (
        (polynomial_P(index) % prime, -pow(index, 6, prime) % prime),
        (1, 0),
    )


def audit_transfer(prime: int, renormalized: Sequence[int]) -> int:
    """Build every prefix product and return the number of checked markers."""

    product: Matrix = ((1, 0), (0, 1))
    for index in range(prime):
        if index:
            product = matrix_multiply(
                transfer_matrix(index - 1, prime), product, prime
            )
        if product[0][0] != renormalized[index]:
            raise ArithmeticError((prime, index, "wrong (1,1) transfer entry"))
        if index and product != (
            (renormalized[index], 0),
            (renormalized[index - 1], 0),
        ):
            raise ArithmeticError((prime, index, "wrong full prefix product"))
    return prime


def berlekamp_massey(sequence: Sequence[int], prime: int) -> int:
    """Return the finite-prefix constant linear complexity over ``F_p``."""

    connection = [1]
    previous = [1]
    complexity = 0
    shift = 1
    old_discrepancy = 1
    for index, value in enumerate(sequence):
        discrepancy = value
        for offset in range(1, complexity + 1):
            discrepancy = (
                discrepancy + connection[offset] * sequence[index - offset]
            ) % prime
        if discrepancy == 0:
            shift += 1
            continue
        saved = connection[:]
        factor = discrepancy * pow(old_discrepancy, -1, prime) % prime
        required = len(previous) + shift
        if len(connection) < required:
            connection.extend([0] * (required - len(connection)))
        for offset, coefficient in enumerate(previous):
            connection[offset + shift] = (
                connection[offset + shift] - factor * coefficient
            ) % prime
        if 2 * complexity <= index:
            complexity = index + 1 - complexity
            previous = saved
            old_discrepancy = discrepancy
            shift = 1
        else:
            shift += 1
    return complexity


def interpolation_degree(values: Sequence[int], prime: int) -> int:
    """Degree of the unique polynomial of degree < p on all of ``F_p``."""

    differences = list(values)
    degree = -1
    for order in range(len(values)):
        if differences[0] % prime:
            degree = order
        differences = [
            (differences[index + 1] - differences[index]) % prime
            for index in range(len(differences) - 1)
        ]
        if not differences:
            break
    return degree


def polynomial_evaluate(coefficients: Sequence[int], argument: int, prime: int) -> int:
    value = 0
    for coefficient in reversed(coefficients):
        value = (value * argument + coefficient) % prime
    return value


def toric_constant_term(index: int) -> int:
    """Return ``CT(Lambda**index)`` by the exact binomial expansion.

    Here

    ``Lambda=(1+x)(1+y)(1+z)((1+y)(1+z)+xyz)/(xyz)``.

    Expanding the second factor according to the number ``k`` of selected
    ``xyz`` terms gives the displayed sum.  The change of variable
    ``ell=index-k`` identifies it with the Apéry number.
    """

    from math import comb

    return sum(
        comb(index, k) ** 2 * comb(2 * index - k, index - k) ** 2
        for k in range(index + 1)
    )


def toric_lambda(x: int, y: int, z: int, prime: int) -> int:
    """Evaluate the fixed Laurent polynomial ``Lambda`` on the torus."""

    numerator = (
        (1 + x)
        * (1 + y)
        * (1 + z)
        * ((1 + y) * (1 + z) + x * y * z)
    )
    return numerator * pow(x * y * z % prime, -1, prime) % prime


def toric_markers(prime: int, last_index: int | None = None) -> list[int]:
    """Return ``-sum Lambda(x,y,z)^j`` for a range of character indices.

    The zeroth power is the constant Laurent polynomial one, including at
    points where ``Lambda`` vanishes.  This convention is the one needed for
    the endpoint ``b_0=1``.
    """

    if last_index is None:
        last_index = prime - 2
    if not 0 <= last_index <= prime - 2:
        raise ValueError("toric no-alias range is 0 <= j <= p-2")
    markers = [0] * (last_index + 1)
    for x in range(1, prime):
        for y in range(1, prime):
            for z in range(1, prime):
                value = toric_lambda(x, y, z, prime)
                power = 1
                for index in range(last_index + 1):
                    markers[index] = (markers[index] - power) % prime
                    power = power * value % prime
    return markers


def toric_complexity_profile() -> dict[str, int]:
    """Fixed input data controlling the toric Kummer realization.

    The four zero components are ``x=-1``, ``y=-1``, ``z=-1``, and
    ``(1+y)(1+z)+xyz=0``.  The six boundary components in ``(P1)^3`` are
    the zero and infinity divisors of the three coordinates.  A log
    resolution can add components, but it is a single fixed resolution and
    hence does not introduce any dependence on ``p`` or ``j``.
    """

    return {
        "ambient_dimension": 3,
        "laurent_half_width": 1,
        "coefficient_rank": 1,
        "zero_components": 4,
        "boundary_components": 6,
        "swan": 0,
    }


def toric_zero_count(point_field_size: int) -> int:
    """Exact number of torus points where ``Lambda`` vanishes."""

    q = point_field_size
    return 4 * q**2 - 14 * q + 13


def toric_complement_count(point_field_size: int) -> int:
    """Exact point count of ``U=(G_m)^3 minus V(Lambda)``."""

    q = point_field_size
    return q**3 - 7 * q**2 + 17 * q - 14


def legendre_symbol(value: int, prime: int) -> int:
    residue = pow(value % prime, (prime - 1) // 2, prime)
    if residue == prime - 1:
        return -1
    return residue


def analyze_prime(prime: int) -> PrimeRecord:
    values = apery_mod_prime(prime)
    renormalized = renormalized_values(values, prime)
    transfer_cases = audit_transfer(prime, renormalized)
    hasse_values = [
        polynomial_evaluate(values, argument, prime)
        for argument in range(1, prime)
    ]
    epsilon = (1 - legendre_symbol(-6, prime)) // 2
    support_lower_bound = (prime - 1) // 2 - epsilon
    mellin_support = sum(value != 0 for value in hasse_values)
    if mellin_support < support_lower_bound:
        raise ArithmeticError((prime, "Mellin support violates factorization bound"))
    if values[0] != 1 or values[-1] != 1:
        raise ArithmeticError((prime, "endpoint normalization failed"))
    if values != list(reversed(values)):
        raise ArithmeticError((prime, "palindromy failed"))
    return PrimeRecord(
        prime=prime,
        zero_count=sum(value == 0 for value in values),
        transfer_cases=transfer_cases,
        bm_b=berlekamp_massey(values, prime),
        bm_B=berlekamp_massey(renormalized, prime),
        interpolation_degree_b=interpolation_degree(values, prime),
        interpolation_degree_B=interpolation_degree(renormalized, prime),
        hasse_at_one=sum(values) % prime,
        mellin_support=mellin_support,
        support_lower_bound=support_lower_bound,
        epsilon=epsilon,
    )


def analyze_primes(primes: Sequence[int], workers: int) -> list[PrimeRecord]:
    if workers <= 1:
        return [analyze_prime(prime) for prime in primes]
    try:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            return sorted(
                pool.map(analyze_prime, primes), key=lambda record: record.prime
            )
    except (OSError, PermissionError):
        # Some locked-down runners disallow POSIX semaphores.  The sequential
        # path performs exactly the same deterministic computation.
        return [analyze_prime(prime) for prime in primes]


def constant_order_two_residuals() -> tuple[int, int, int]:
    """Residuals forced by a hypothetical constant order-two recurrence."""

    values = apery_integers(5)
    first = 2 * values[4] - 45 * values[3] + 79 * values[2]
    second = 2 * values[5] - 45 * values[4] + 79 * values[3]
    return first, second, gcd(first, second)


def kummer_geometric_conductor(character_index: int, prime: int) -> int:
    """Rank + number of singularities + Swan for the Kummer factor alone."""

    if not 0 <= character_index < prime - 1:
        raise ValueError("Kummer exponents are taken modulo p-1")
    return 1 if character_index == 0 else 3


def j_trace_sheaf_conductor() -> int:
    """Refuse to assign a conductor to an object that has not been constructed."""

    raise NotImplementedError(
        "NOT_DEFINED: no j-space sheaf with trace j -> b_j is specified"
    )


def smooth_locus_counterexample() -> dict[str, object]:
    """Return the exact p=31 boundary failure for the smooth K3 locus."""

    prime = 31
    values = apery_mod_prime(prime)
    roots = (14, 20)
    hasse_at_roots = tuple(
        polynomial_evaluate(values, root, prime) for root in roots
    )
    output: dict[int, tuple[int, int, int]] = {}
    for index in (8, 22):
        full_sum = sum(
            polynomial_evaluate(values, argument, prime)
            * pow(argument, -index, prime)
            for argument in range(1, prime)
        ) % prime
        boundary = sum(
            value * pow(root, -index, prime)
            for root, value in zip(roots, hasse_at_roots)
        ) % prime
        output[index] = (full_sum, boundary, (full_sum - boundary) % prime)
    return {
        "prime": prime,
        "roots": roots,
        "hasse_at_roots": hasse_at_roots,
        "zero_indices": (8, 22),
        "sums": output,
    }


def comma_join(values: Iterable[int]) -> str:
    return ",".join(str(value) for value in values)


def print_report(records: Sequence[PrimeRecord], limit: int, workers: int) -> None:
    primes = [record.prime for record in records]
    total_cases = sum(record.transfer_cases for record in records)
    bm_b_exceptions = [
        record.prime
        for record in records
        if record.bm_b != (record.prime + 1) // 2
    ]
    bm_B_exceptions = [
        record.prime
        for record in records
        if record.bm_B != (record.prime + 1) // 2
    ]
    degree_counts: dict[int, int] = {}
    for record in records:
        drop = record.prime - 1 - record.interpolation_degree_b
        degree_counts[drop] = degree_counts.get(drop, 0) + 1
    h1_zero = sum(record.hasse_at_one == 0 for record in records)
    minimum_support_record = min(
        records, key=lambda record: record.mellin_support / (record.prime - 1)
    )
    first_residual, second_residual, residual_gcd = constant_order_two_residuals()

    print("ORACLE C FINITE EXPLORATION")
    print(f"limit={limit} workers={workers}")
    print(
        f"primes={len(records)} first={primes[0]} last={primes[-1]} "
        f"transfer_cases={total_cases}"
    )
    print("transfer_marker=(M(j-1)...M(0))[1,1]=B_j (one-based matrix entry)")
    print(
        "constant_order2_residuals="
        f"{first_residual},{second_residual} gcd={residual_gcd}"
    )
    print(
        f"BM_b_default={(len(records) - len(bm_b_exceptions))} "
        f"exceptions={comma_join(bm_b_exceptions)}"
    )
    print(
        f"BM_B_default={(len(records) - len(bm_B_exceptions))} "
        f"exceptions={comma_join(bm_B_exceptions)}"
    )
    print(
        "interpolation_degree_b_drop_counts="
        + ",".join(f"{drop}:{degree_counts[drop]}" for drop in sorted(degree_counts))
    )
    print(f"H_p(1)_zero_count={h1_zero}")
    print(
        "minimum_mellin_support_ratio="
        f"{minimum_support_record.mellin_support}/"
        f"{minimum_support_record.prime - 1} at p={minimum_support_record.prime}"
    )
    print("raw_support_factorization_bound=PASS")
    toric_profile = toric_complexity_profile()
    toric_samples_pass = True
    for prime in TORIC_SAMPLE_PRIMES:
        toric_samples_pass &= toric_markers(prime) == apery_mod_prime(prime)[:-1]
    if not toric_samples_pass:
        raise ArithmeticError("toric marked-coordinate sample failed")
    print(
        "toric_marker_samples="
        f"{comma_join(TORIC_SAMPLE_PRIMES)} "
        f"{'PASS' if toric_samples_pass else 'FAIL'}"
    )
    print(
        "toric_complexity="
        f"dimension={toric_profile['ambient_dimension']} "
        f"half_width={toric_profile['laurent_half_width']} "
        f"rank={toric_profile['coefficient_rank']} "
        f"input_divisor_components="
        f"{toric_profile['zero_components'] + toric_profile['boundary_components']} "
        f"Swan={toric_profile['swan']}"
    )
    print("toric_U_point_count=q^3-7q^2+17q-14")
    print("toric_raw_endpoint=-#U=14(mod p), separate endpoint=1")
    print("kummer_factor_geometric_conductor=1(trivial),3(nontrivial)")
    try:
        j_trace_sheaf_conductor()
    except NotImplementedError as error:
        print(f"j_trace_sheaf_conductor={error}")
    boundary = smooth_locus_counterexample()
    print(
        "p31_smooth_boundary="
        f"roots={boundary['roots']} H={boundary['hasse_at_roots']} "
        f"sums={boundary['sums']}"
    )
    print("samples=p,Z,L_b,L_B,deg_I_b,support,lower_bound")
    record_by_prime = {record.prime: record for record in records}
    for prime in SAMPLE_PRIMES:
        if prime not in record_by_prime:
            continue
        record = record_by_prime[prime]
        print(
            f"{prime},{record.zero_count},{record.bm_b},{record.bm_B},"
            f"{record.interpolation_degree_b},{record.mellin_support},"
            f"{record.support_lower_bound}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=DEFAULT_LIMIT)
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="processes used for independent prime computations",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    primes = primes_up_to(args.limit)
    if not primes:
        raise SystemExit("no primes p>=5 in the requested range")
    records = analyze_primes(primes, args.workers)
    print_report(records, args.limit, args.workers)


if __name__ == "__main__":
    main()
