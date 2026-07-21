#!/usr/bin/env python3
"""Exact verification for ``oracleC_result.tex``.

The verifier checks algebraic identities and finite-field computations only.
It independently expands the fixed toric Laurent polynomial underlying the
marked coordinate.  It does not treat finite reconnaissance as a proof of a
cohomological bound, and it explicitly refuses to manufacture a conductor
for a nonexistent sheaf on the character-index line.
"""

from __future__ import annotations

from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from math import comb
from pathlib import Path
import sys
import traceback
from typing import Callable, Sequence


if not __debug__:
    print("FAIL: oracleC_verify.py refuses python -O because assertions are checks")
    raise SystemExit(2)


HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
SPEC = spec_from_file_location("oracleC_explore", HERE / "oracleC_explore.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load oracleC_explore.py")
EXPLORE = module_from_spec(SPEC)
sys.modules[SPEC.name] = EXPLORE
SPEC.loader.exec_module(EXPLORE)


Polynomial = list[int]
IntegerMatrix = list[list[int]]


def direct_apery(index: int) -> int:
    return sum(
        comb(index, summation_index) ** 2
        * comb(index + summation_index, summation_index) ** 2
        for summation_index in range(index + 1)
    )


LaurentPolynomial = dict[tuple[int, int, int], int]


def laurent_multiply(
    left: LaurentPolynomial, right: LaurentPolynomial
) -> LaurentPolynomial:
    output: LaurentPolynomial = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                left_exponent[coordinate] + right_exponent[coordinate]
                for coordinate in range(3)
            )
            output[exponent] = (
                output.get(exponent, 0) + left_coefficient * right_coefficient
            )
    return {exponent: value for exponent, value in output.items() if value}


def laurent_power(base: LaurentPolynomial, exponent: int) -> LaurentPolynomial:
    output: LaurentPolynomial = {(0, 0, 0): 1}
    while exponent:
        if exponent & 1:
            output = laurent_multiply(output, base)
        exponent >>= 1
        if exponent:
            base = laurent_multiply(base, base)
    return output


def toric_lambda_polynomial() -> LaurentPolynomial:
    output: LaurentPolynomial = {(0, 0, 0): 1}
    for factor in (
        {(0, 0, 0): 1, (1, 0, 0): 1},
        {(0, 0, 0): 1, (0, 1, 0): 1},
        {(0, 0, 0): 1, (0, 0, 1): 1},
        {
            (0, 0, 0): 1,
            (0, 1, 0): 1,
            (0, 0, 1): 1,
            (0, 1, 1): 1,
            (1, 1, 1): 1,
        },
    ):
        output = laurent_multiply(output, factor)
    return {
        (exponent[0] - 1, exponent[1] - 1, exponent[2] - 1): coefficient
        for exponent, coefficient in output.items()
    }


def check_toric_marked_coordinate() -> None:
    laurent = toric_lambda_polynomial()
    for coordinate in range(3):
        exponents = [exponent[coordinate] for exponent in laurent]
        assert min(exponents) == -1
        assert max(exponents) == 1

    for index in range(9):
        expanded_constant_term = laurent_power(laurent, index).get((0, 0, 0), 0)
        assert expanded_constant_term == EXPLORE.toric_constant_term(index)
        assert expanded_constant_term == direct_apery(index)

    for prime in EXPLORE.TORIC_SAMPLE_PRIMES:
        markers = EXPLORE.toric_markers(prime)
        values = EXPLORE.apery_mod_prime(prime)
        assert markers == values[:-1]
        assert values[-1] == values[0] == 1
        zero_count = sum(
            EXPLORE.toric_lambda(x, y, z, prime) == 0
            for x in range(1, prime)
            for y in range(1, prime)
            for z in range(1, prime)
        )
        assert zero_count == EXPLORE.toric_zero_count(prime)
        assert (prime - 1) ** 3 - zero_count == EXPLORE.toric_complement_count(prime)

    prime = 7
    raw_last_coordinate = -sum(
        pow(EXPLORE.toric_lambda(x, y, z, prime), prime - 1, prime)
        for x in range(1, prime)
        for y in range(1, prime)
        for z in range(1, prime)
    ) % prime
    assert raw_last_coordinate == 0
    assert EXPLORE.apery_mod_prime(prime)[-1] == 1

    prime = 11
    quadratic_sum = sum(
        EXPLORE.legendre_symbol(EXPLORE.toric_lambda(x, y, z, prime), prime)
        for x in range(1, prime)
        for y in range(1, prime)
        for z in range(1, prime)
    )
    assert quadratic_sum == 33
    assert EXPLORE.apery_mod_prime(prime)[5] == 0

    assert EXPLORE.toric_complexity_profile() == {
        "ambient_dimension": 3,
        "laurent_half_width": 1,
        "coefficient_rank": 1,
        "zero_components": 4,
        "boundary_components": 6,
        "swan": 0,
    }


def matrix_multiply_integer(
    left: tuple[tuple[int, int], tuple[int, int]],
    right: tuple[tuple[int, int], tuple[int, int]],
) -> tuple[tuple[int, int], tuple[int, int]]:
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def check_integer_recurrences() -> None:
    expected = [1, 5, 73, 1445, 33001, 819005]
    values = EXPLORE.apery_integers(15)
    assert values[: len(expected)] == expected
    assert values == [direct_apery(index) for index in range(16)]
    renormalized = [
        direct_apery(index)
        * product(range(1, index + 1)) ** 3
        for index in range(16)
    ]
    for index in range(1, 15):
        assert (
            renormalized[index + 1]
            == EXPLORE.polynomial_P(index) * renormalized[index]
            - index**6 * renormalized[index - 1]
        )


def product(values: Sequence[int] | range) -> int:
    result = 1
    for value in values:
        result *= value
    return result


def check_transfer_all_primes() -> None:
    cases = 0
    for prime in EXPLORE.primes_up_to(2_000):
        values = EXPLORE.apery_mod_prime(prime)
        renormalized = EXPLORE.renormalized_values(values, prime)
        cases += EXPLORE.audit_transfer(prime, renormalized)
        assert all(
            (renormalized[index] == 0) == (values[index] == 0)
            for index in range(prime)
        )
    assert len(EXPLORE.primes_up_to(2_000)) == 301
    assert cases == 277_045


def check_fundamental_matrix() -> None:
    values = EXPLORE.apery_integers(14)
    factorial = 1
    B = [1]
    for index in range(1, 15):
        factorial *= index
        B.append(factorial**3 * values[index])
    D = [0, 1]
    for index in range(1, 14):
        D.append(EXPLORE.polynomial_P(index) * D[index] - index**6 * D[index - 1])
    fundamental = ((5, 1), (1, 0))
    assert fundamental == ((B[1], D[1]), (B[0], D[0]))
    for index in range(1, 14):
        assert fundamental == (
            (B[index], D[index]),
            (B[index - 1], D[index - 1]),
        )
        assert (
            fundamental[0][0] * fundamental[1][1]
            - fundamental[0][1] * fundamental[1][0]
            == -(product(range(1, index)) ** 6)
        )
        matrix = (
            (EXPLORE.polynomial_P(index), -(index**6)),
            (1, 0),
        )
        fundamental = matrix_multiply_integer(matrix, fundamental)


def bareiss_determinant(matrix: IntegerMatrix) -> int:
    if not matrix:
        return 1
    work = [row[:] for row in matrix]
    size = len(work)
    assert all(len(row) == size for row in work)
    sign = 1
    denominator = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (
                    row
                    for row in range(pivot_index + 1, size)
                    if work[row][pivot_index] != 0
                ),
                None,
            )
            if swap is None:
                return 0
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                assert numerator % denominator == 0
                work[row][column] = numerator // denominator
        denominator = pivot
        for row in range(pivot_index + 1, size):
            work[row][pivot_index] = 0
    return sign * work[-1][-1]


def recurrence_ansatz_matrix(values: Sequence[int], degree: int, rows: int) -> IntegerMatrix:
    return [
        [
            values[index + shift] * index**power
            for shift in range(3)
            for power in range(degree + 1)
        ]
        for index in range(rows)
    ]


def rank_mod_prime(matrix: IntegerMatrix, prime: int) -> int:
    work = [[entry % prime for entry in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0]) if work else 0
    rank = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, row_count) if work[row][column]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [value * inverse % prime for value in work[rank]]
        for row in range(row_count):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                (left - factor * right) % prime
                for left, right in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


def check_recurrence_complexity() -> None:
    values = EXPLORE.apery_integers(30)
    first, second, common = EXPLORE.constant_order_two_residuals()
    assert (first, second, common) == (6744, 267120, 24)
    expected_determinant = -(
        2**26
        * 3**12
        * 5**6
        * 7**2
        * 26309
        * 50077
        * 171131
    )
    degree_two = recurrence_ansatz_matrix(values, degree=2, rows=9)
    assert bareiss_determinant(degree_two) == expected_determinant
    coefficient_vector = (
        1,
        3,
        3,
        1,
        -117,
        -231,
        -153,
        -34,
        8,
        12,
        6,
        1,
    )
    degree_three = recurrence_ansatz_matrix(values, degree=3, rows=25)
    assert all(
        sum(coefficient * entry for coefficient, entry in zip(coefficient_vector, row))
        == 0
        for row in degree_three
    )
    for prime in EXPLORE.primes_up_to(2_000):
        if prime >= 11:
            assert rank_mod_prime(degree_two, prime) == 9
        if prime >= 13:
            assert rank_mod_prime(degree_three[: min(25, prime - 2)], prime) == 11


def poly_trim(polynomial: Polynomial) -> Polynomial:
    while len(polynomial) > 1 and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def poly_add(left: Sequence[int], right: Sequence[int]) -> Polynomial:
    output = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        output[index] += value
    for index, value in enumerate(right):
        output[index] += value
    return poly_trim(output)


def poly_scale(polynomial: Sequence[int], scalar: int) -> Polynomial:
    return poly_trim([scalar * value for value in polynomial])


def poly_multiply(left: Sequence[int], right: Sequence[int]) -> Polynomial:
    output = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            output[left_degree + right_degree] += left_value * right_value
    return poly_trim(output)


def poly_power(polynomial: Sequence[int], exponent: int) -> Polynomial:
    output = [1]
    base = list(polynomial)
    while exponent:
        if exponent & 1:
            output = poly_multiply(output, base)
        exponent >>= 1
        if exponent:
            base = poly_multiply(base, base)
    return output


def poly_evaluate_integer(polynomial: Sequence[int], argument: int) -> int:
    value = 0
    for coefficient in reversed(polynomial):
        value = value * argument + coefficient
    return value


def shifted_P(shift: int) -> Polynomial:
    linear = [shift, 1]
    return poly_add(
        poly_add(
            poly_scale(poly_power(linear, 3), 34),
            poly_scale(poly_power(linear, 2), 51),
        ),
        poly_add(poly_scale(linear, 27), [5]),
    )


PolyMatrix = tuple[
    tuple[Polynomial, Polynomial], tuple[Polynomial, Polynomial]
]


def poly_matrix_multiply(left: PolyMatrix, right: PolyMatrix) -> PolyMatrix:
    def entry(row: int, column: int) -> Polynomial:
        return poly_add(
            poly_multiply(left[row][0], right[0][column]),
            poly_multiply(left[row][1], right[1][column]),
        )

    return ((entry(0, 0), entry(0, 1)), (entry(1, 0), entry(1, 1)))


def check_shifted_products_and_poles() -> None:
    identity: PolyMatrix = (([1], [0]), ([0], [1]))
    product_matrix = identity
    N = [[0], [1]]
    integer_values = EXPLORE.apery_integers(12)
    factorials = [product(range(1, index + 1)) for index in range(13)]
    B = [
        factorials[index] ** 3 * integer_values[index]
        for index in range(13)
    ]
    for length in range(1, 10):
        shift = length - 1
        matrix: PolyMatrix = (
            (shifted_P(shift), poly_scale(poly_power([shift, 1], 6), -1)),
            ([1], [0]),
        )
        product_matrix = poly_matrix_multiply(matrix, product_matrix)
        expected_degrees = (
            (3 * length, 3 * length + 3),
            (3 * length - 3, 3 * length if length >= 2 else -1),
        )
        actual_degrees = tuple(
            tuple(-1 if entry == [0] else len(entry) - 1 for entry in row)
            for row in product_matrix
        )
        assert actual_degrees == expected_degrees

        next_N = poly_add(
            poly_multiply(shifted_P(length), N[-1]),
            poly_scale(poly_multiply(poly_power([length, 1], 6), N[-2]), -1),
        )
        N.append(next_N)
        h = length + 1
        assert len(N[h]) - 1 == 3 * (h - 1)
        for endpoint in range(1, h + 1):
            assert poly_evaluate_integer(N[h], -endpoint) == (
                (-1) ** (endpoint - 1) * B[endpoint - 1] * B[h - endpoint]
            )


def raw_mellin_sum(values: Sequence[int], prime: int, index: int) -> int:
    return -sum(
        EXPLORE.polynomial_evaluate(values, argument, prime)
        * pow(argument, -index, prime)
        for argument in range(1, prime)
    ) % prime


def check_mellin_and_group_algebra() -> None:
    for prime in (5, 7, 11, 17, 31, 101):
        values = EXPLORE.apery_mod_prime(prime)
        assert raw_mellin_sum(values, prime, 0) == 2
        for index in range(1, prime - 1):
            assert raw_mellin_sum(values, prime, index) == values[index]
        generator = EXPLORE.primitive_root(prime)
        hasse_on_powers = [
            EXPLORE.polynomial_evaluate(values, pow(generator, exponent, prime), prime)
            for exponent in range(prime - 1)
        ]
        support = sum(value != 0 for value in hasse_on_powers)
        epsilon = (1 - EXPLORE.legendre_symbol(-6, prime)) // 2
        assert support >= (prime - 1) // 2 - epsilon
        for index in range(prime - 1):
            group_value = -sum(
                value * pow(generator, -exponent * index, prime)
                for exponent, value in enumerate(hasse_on_powers)
            ) % prime
            assert group_value == raw_mellin_sum(values, prime, index)

    prime = 11
    values = EXPLORE.apery_mod_prime(prime)
    generator = EXPLORE.primitive_root(prime)
    assert values[5] == 0
    assert EXPLORE.polynomial_evaluate(values, 1, prime) == 3
    omitted_t_one = -sum(
        EXPLORE.polynomial_evaluate(values, pow(generator, exponent, prime), prime)
        * pow(generator, -5 * exponent, prime)
        for exponent in range(1, prime - 1)
    ) % prime
    assert omitted_t_one == 3


def check_smooth_locus_counterexample() -> None:
    record = EXPLORE.smooth_locus_counterexample()
    assert record == {
        "prime": 31,
        "roots": (14, 20),
        "hasse_at_roots": (7, 7),
        "zero_indices": (8, 22),
        "sums": {8: (0, 11, 20), 22: (0, 11, 20)},
    }


def check_finite_complexity_samples() -> None:
    expected = {
        5: (3, 3),
        7: (4, 4),
        11: (6, 6),
        19: (10, 9),
        23: (11, 12),
        47: (23, 24),
        67: (33, 34),
        181: (91, 91),
        827: (413, 414),
        1999: (1000, 1000),
    }
    for prime, target in expected.items():
        values = EXPLORE.apery_mod_prime(prime)
        renormalized = EXPLORE.renormalized_values(values, prime)
        assert (
            EXPLORE.berlekamp_massey(values, prime),
            EXPLORE.berlekamp_massey(renormalized, prime),
        ) == target
    for prime in (5, 7, 11, 31, 181, 827, 1999):
        values = EXPLORE.apery_mod_prime(prime)
        degree = EXPLORE.interpolation_degree(values, prime)
        assert degree in (prime - 1, prime - 3)
        assert (degree == prime - 3) == (sum(values) % prime == 0)


def check_kummer_conductor_scope() -> None:
    for prime in EXPLORE.primes_up_to(2_000):
        conductors = [
            EXPLORE.kummer_geometric_conductor(index, prime)
            for index in range(prime - 1)
        ]
        assert conductors.count(1) == 1
        assert conductors.count(3) == prime - 2
    try:
        EXPLORE.j_trace_sheaf_conductor()
    except NotImplementedError as error:
        assert str(error).startswith("NOT_DEFINED:")
    else:
        raise AssertionError("a conductor was assigned to an unspecified j-sheaf")

    # Grothendieck--Ogg--Shafarevich with s=4 and zero Swan conductor.
    for rank, generic_degree, total_bound in ((3, 6, 12), (22, 44, 88)):
        assert rank * (4 - 2) == generic_degree
        assert rank * 4 == total_bound


CHECKS: tuple[tuple[str, Callable[[], None]], ...] = (
    ("integer Apéry recurrences", check_integer_recurrences),
    ("fixed toric marked coordinate", check_toric_marked_coordinate),
    ("all p<=2000 transfer markers", check_transfer_all_primes),
    ("rank-two fundamental matrix", check_fundamental_matrix),
    ("recurrence coefficient complexity", check_recurrence_complexity),
    ("shifted product degrees and uncancelled poles", check_shifted_products_and_poles),
    ("endpoint-correct Mellin/group-algebra identity", check_mellin_and_group_algebra),
    ("p=31 smooth-locus failure", check_smooth_locus_counterexample),
    ("finite linear/interpolation complexity samples", check_finite_complexity_samples),
    ("Kummer conductor scope", check_kummer_conductor_scope),
)


def main() -> None:
    failures = 0
    for name, check in CHECKS:
        try:
            check()
        except Exception:  # noqa: BLE001 - verifier must report every failed check
            failures += 1
            print(f"FAIL: {name}")
            traceback.print_exc()
        else:
            print(f"PASS: {name}")
    if failures:
        print(f"FAIL: {failures} Oracle C verification group(s) failed")
        raise SystemExit(1)
    print(f"PASS: all {len(CHECKS)} Oracle C verification groups")


if __name__ == "__main__":
    main()
