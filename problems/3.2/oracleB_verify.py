#!/usr/bin/env python3
"""Exact verification for ``oracleB_result.tex``.

The full ``p <= 20000`` zero-set scan is finite evidence only.  The proof
checks below verify the algebraic identities used in the rigorous fixed-anchor
proposition and the counterexample separating coefficient zeros from roots.
Every check prints PASS/FAIL; any failure gives a nonzero exit status.
"""

from __future__ import annotations

from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from math import comb
from pathlib import Path
import sys
import traceback
from typing import Callable, Sequence


if not __debug__:
    print("FAIL: oracleB_verify.py refuses python -O because assertions are checks")
    raise SystemExit(2)


HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
SPEC = spec_from_file_location("oracleB_explore", HERE / "oracleB_explore.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load oracleB_explore.py")
EXPLORE = module_from_spec(SPEC)
sys.modules[SPEC.name] = EXPLORE
SPEC.loader.exec_module(EXPLORE)


Polynomial = list[int]


def poly_trim(polynomial: Polynomial) -> Polynomial:
    while len(polynomial) > 1 and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def poly_mul(left: Sequence[int], right: Sequence[int], prime: int) -> Polynomial:
    product = [0] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            product[left_degree + right_degree] = (
                product[left_degree + right_degree]
                + left_coefficient * right_coefficient
            ) % prime
    return poly_trim(product)


def poly_divmod(
    dividend: Sequence[int], divisor: Sequence[int], prime: int
) -> tuple[Polynomial, Polynomial]:
    numerator = poly_trim(list(dividend))
    denominator = poly_trim(list(divisor))
    if denominator == [0]:
        raise ZeroDivisionError
    if len(numerator) < len(denominator):
        return [0], numerator
    quotient = [0] * (len(numerator) - len(denominator) + 1)
    inverse_lead = pow(denominator[-1], -1, prime)
    for shift in range(len(quotient) - 1, -1, -1):
        coefficient = numerator[shift + len(denominator) - 1] * inverse_lead % prime
        quotient[shift] = coefficient
        for degree, value in enumerate(denominator):
            numerator[shift + degree] = (
                numerator[shift + degree] - coefficient * value
            ) % prime
    return poly_trim(quotient), poly_trim(numerator[: len(denominator) - 1] or [0])


def poly_monic(polynomial: Sequence[int], prime: int) -> Polynomial:
    trimmed = poly_trim(list(polynomial))
    inverse = pow(trimmed[-1], -1, prime)
    return [(coefficient * inverse) % prime for coefficient in trimmed]


def poly_gcd(left: Sequence[int], right: Sequence[int], prime: int) -> Polynomial:
    first = poly_trim(list(left))
    second = poly_trim(list(right))
    while second != [0]:
        _, remainder = poly_divmod(first, second, prime)
        first, second = second, remainder
    return poly_monic(first, prime)


def poly_mod(
    polynomial: Sequence[int], modulus: Sequence[int], prime: int
) -> Polynomial:
    return poly_divmod(polynomial, modulus, prime)[1]


def poly_powmod(
    base: Sequence[int], exponent: int, modulus: Sequence[int], prime: int
) -> Polynomial:
    result = [1]
    power = poly_mod(base, modulus, prime)
    while exponent:
        if exponent & 1:
            result = poly_mod(poly_mul(result, power, prime), modulus, prime)
        exponent >>= 1
        if exponent:
            power = poly_mod(poly_mul(power, power, prime), modulus, prime)
    return result


def poly_derivative(polynomial: Sequence[int], prime: int) -> Polynomial:
    if len(polynomial) <= 1:
        return [0]
    return poly_trim(
        [degree * polynomial[degree] % prime for degree in range(1, len(polynomial))]
    )


def poly_evaluate(polynomial: Sequence[int], argument: int, prime: int) -> int:
    value = 0
    for coefficient in reversed(polynomial):
        value = (value * argument + coefficient) % prime
    return value


def formal_square_root(square: Sequence[int], degree: int, prime: int) -> Polynomial:
    """Return the constant-one square root, if its full square agrees."""

    if square[0] % prime != 1:
        raise ArithmeticError("formal square root requires constant coefficient one")
    inverse_two = pow(2, -1, prime)
    root = [1]
    padded = list(square) + [0] * max(0, 2 * degree + 1 - len(square))
    for index in range(1, degree + 1):
        convolution = sum(
            root[left] * root[index - left] for left in range(1, index)
        )
        root.append((padded[index] - convolution) * inverse_two % prime)
    if poly_mul(root, root, prime) != poly_trim(padded[: 2 * degree + 1]):
        raise ArithmeticError("candidate formal square root fails in high degree")
    return root


def direct_apery(index: int) -> int:
    return sum(
        comb(index, summation_index) ** 2
        * comb(index + summation_index, summation_index) ** 2
        for summation_index in range(index + 1)
    )


def check_apery_recurrence() -> None:
    expected_prefix = [1, 5, 73, 1445, 33001, 819005]
    values = EXPLORE.apery_integers(12)
    assert values[: len(expected_prefix)] == expected_prefix
    assert values == [direct_apery(index) for index in range(13)]
    assert all(value > 0 for value in values)


def check_full_zero_bank() -> None:
    records = [
        EXPLORE.ZeroRecord(prime, EXPLORE.apery_zero_set(prime))
        for prime in EXPLORE.primes_up_to(20_000)
        if prime >= 5
    ]
    EXPLORE.validate_zero_records(records)
    histogram = Counter(len(record.zeros) for record in records)
    assert len(records) == 2260
    assert histogram == Counter({0: 1356, 2: 695, 4: 176, 6: 27, 8: 4, 1: 2})
    assert [record.prime for record in records if len(record.zeros) % 2] == [11, 3137]
    assert next(record.zeros for record in records if record.prime == 181) == (
        19,
        47,
        133,
        161,
    )
    binary, digest = EXPLORE.load_binary_zero_sets(
        HERE / "data_zp_pairs.bin", 20_000
    )
    assert digest == "8746d0b400c1b669b001eae955c602908a10c9ee4cb3cac62c6676ea2ddd874d"
    for record in records:
        expected = binary.get(record.prime, ())
        if record.prime == 5:
            expected = (1, 3)
        assert record.zeros == expected


def check_fixed_anchor_lucas() -> None:
    primes = [prime for prime in EXPLORE.primes_up_to(1000) if prime >= 5]
    zero_sets = {prime: set(EXPLORE.apery_zero_set(prime)) for prime in primes}
    values = EXPLORE.apery_integers(100)
    for anchor in range(-50, 101):
        flattened = anchor if anchor >= 0 else -anchor - 1
        for prime in primes:
            hit = anchor % prime in zero_sets[prime]
            assert not hit or values[flattened] % prime == 0
            if prime > abs(anchor):
                assert hit == (values[flattened] % prime == 0)


def check_root_coefficient_separation() -> None:
    prime = 7
    hasse = [value % prime for value in EXPLORE.apery_integers(prime - 1)]
    assert hasse == [1, 5, 3, 3, 3, 5, 1]
    factor = poly_mul([(-1) % prime, 1], [1, 0, 1], prime)
    assert poly_mul(factor, factor, prime) == hasse
    assert EXPLORE.apery_zero_set(prime) == ()
    assert [x for x in range(prime) if poly_evaluate(hasse, x, prime) == 0] == [1]


def check_small_sym2_identities() -> None:
    primes = [prime for prime in EXPLORE.primes_up_to(101) if prime >= 5]
    integer_values = EXPLORE.apery_integers(100)
    plus_classes = {1, 11, 17, 19}
    minus_classes = {5, 7, 13, 23}
    for prime in primes:
        hasse = [value % prime for value in integer_values[:prime]]
        epsilon = int(prime % 24 in EXPLORE.CORRECTED_SQUARE_CLASSES)
        delta = [1, (-34) % prime, 1]
        if epsilon:
            quotient, remainder = poly_divmod(hasse, delta, prime)
            assert remainder == [0]
        else:
            quotient = hasse
        degree = (prime - 1 - 2 * epsilon) // 2
        square_root = formal_square_root(quotient, degree, prime)
        reconstructed = poly_mul(square_root, square_root, prime)
        if epsilon:
            reconstructed = poly_mul(delta, reconstructed, prime)
        assert reconstructed == hasse
        assert len(poly_gcd(square_root, poly_derivative(square_root, prime), prime)) == 1
        assert len(poly_gcd(square_root, delta, prime)) == 1
        coordinate = [0, 1]
        assert poly_powmod(coordinate, prime * prime, square_root, prime) == coordinate
        frobenius_once = poly_powmod(coordinate, prime, square_root, prime)
        linear_degree = len(
            poly_gcd(
                square_root,
                [
                    (
                        (frobenius_once[index] if index < len(frobenius_once) else 0)
                        - (coordinate[index] if index < len(coordinate) else 0)
                    )
                    % prime
                    for index in range(max(len(frobenius_once), len(coordinate)))
                ],
                prime,
            )
        ) - 1
        assert (degree - linear_degree) % 2 == 0
        reciprocal = list(reversed(square_root))
        if prime % 24 in plus_classes:
            assert reciprocal == square_root
        elif prime % 24 in minus_classes:
            assert reciprocal == [(-coefficient) % prime for coefficient in square_root]
        else:
            raise AssertionError("unclassified prime residue")


def check_proposed_legendre_pencil_is_not_identified() -> None:
    """Check the raw p=7 Hasse invariant of the equation named in the spec.

    This does not rule out an undisplayed base change/gauge.  It verifies the
    narrower statement made in the result: in the printed t-coordinate the
    equation's raw Deuring polynomial is not the Apéry H_7.
    """

    prime = 7
    t = [0, 1]
    t_squared = poly_mul(t, t, prime)
    size = max(len(t), len(t_squared))
    lam = [0] * size
    for degree in range(size):
        lam[degree] = (
            (t[degree] if degree < len(t) else 0)
            - (t_squared[degree] if degree < len(t_squared) else 0)
        ) % prime
    legendre_hasse = [1]
    power = [1]
    deuring_coefficients = [1, 2, 2, 1]
    for exponent, coefficient in enumerate(deuring_coefficients):
        if exponent:
            power = poly_mul(power, lam, prime)
        if len(legendre_hasse) < len(power):
            legendre_hasse.extend([0] * (len(power) - len(legendre_hasse)))
        for degree, value in enumerate(power):
            legendre_hasse[degree] = (
                legendre_hasse[degree] + coefficient * value
            ) % prime
    # The loop added the exponent-zero term twice; remove the extra 1.
    legendre_hasse[0] = (legendre_hasse[0] - 1) % prime
    apery_hasse = [value % prime for value in EXPLORE.apery_integers(6)]
    assert poly_trim(legendre_hasse) != apery_hasse
    assert [
        value
        for value in range(prime)
        if poly_evaluate(legendre_hasse, value, prime) == 0
    ] == [4]


def check_result_quantifiers() -> None:
    result = (HERE / "oracleB_result.tex").read_text(encoding="utf-8")
    required_fragments = (
        r"\forall c\in\mathbb Z",
        r"\sup_{m\in I_N}",
        r"[t^j]H_p(t)",
        r"\label{eq:oracleB-mh2}",
        "STALL REPORT",
        "Remark~\\ref{rem:orbit}",
        "Remark~\\ref{rem:open}",
        "Remark~\\ref{rem:squareness}",
    )
    for fragment in required_fragments:
        assert fragment in result, fragment


CHECKS: tuple[tuple[str, Callable[[], None]], ...] = (
    ("Apéry definition and integral recurrence", check_apery_recurrence),
    ("p<=20000 zero bank and reflection", check_full_zero_bank),
    ("fixed-anchor Lucas implication", check_fixed_anchor_lucas),
    ("coefficient zeros differ from Hasse roots", check_root_coefficient_separation),
    ("small-prime Sym^2 identities", check_small_sym2_identities),
    ("raw proposed pencil mismatch", check_proposed_legendre_pencil_is_not_identified),
    ("result quantifiers and citations", check_result_quantifiers),
)


def main() -> None:
    failures = 0
    for name, check in CHECKS:
        try:
            check()
        except Exception:  # noqa: BLE001 - a verifier must report every failure
            failures += 1
            print(f"FAIL: {name}")
            traceback.print_exc()
        else:
            print(f"PASS: {name}")
    if failures:
        print(f"FAIL: {failures}/{len(CHECKS)} checks failed")
        raise SystemExit(1)
    print(f"PASS: all {len(CHECKS)} oracle-B checks")


if __name__ == "__main__":
    main()
