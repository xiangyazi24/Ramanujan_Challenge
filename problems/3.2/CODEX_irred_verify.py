#!/usr/bin/env python3
"""Exact certificates for the Apery gap-cover monodromy report.

The expensive critical-value calculations are done modulo good primes.  This
is a proof over Q, not a heuristic: reduction preserves the expected degree,
so squarefreeness and pairwise coprimality modulo one prime rule out a repeated
factor or a common factor in Q[T].  A second prime is used as an independent
implementation check.
"""

from __future__ import annotations

import argparse
import hashlib
import warnings
from math import factorial, prod
from typing import Dict, List, Sequence, Tuple

from sympy import Poly, ZZ, diff, expand, factor_list, primerange, symbols
from sympy.utilities.exceptions import SymPyDeprecationWarning


x, y, T = symbols("x y T")
warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)


def apery_P(z):
    return 34 * z**3 + 51 * z**2 + 27 * z + 5


def integer_numerators(height: int) -> List[Poly]:
    numerators = [Poly(0, x, domain=ZZ), Poly(1, x, domain=ZZ)]
    for h in range(1, height):
        numerators.append(
            Poly(apery_P(x + h), x, domain=ZZ) * numerators[h]
            - Poly((x + h) ** 6, x, domain=ZZ) * numerators[h - 1]
        )
    return numerators


def mod_numerators(height: int, prime: int) -> List[Poly]:
    numerators = [
        Poly(0, x, modulus=prime),
        Poly(1, x, modulus=prime),
    ]
    for h in range(1, height):
        numerators.append(
            Poly(apery_P(x + h), x, modulus=prime) * numerators[h]
            - Poly((x + h) ** 6, x, modulus=prime) * numerators[h - 1]
        )
    return numerators


def q_poly(h: int, *, modulus: int | None = None) -> Poly:
    domain_args = {"modulus": modulus} if modulus is not None else {"domain": ZZ}
    q = Poly(1, x, **domain_args)
    for j in range(1, h + 1):
        q *= Poly(x + j, x, **domain_args)
    return q


def apery_numbers(height: int) -> List[int]:
    values = [1]
    if height == 0:
        return values
    values.append(5)
    for n in range(1, height):
        numerator = int(apery_P(n)) * values[n] - n**3 * values[n - 1]
        denominator = (n + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values


def check_exact_structure(numerators: Sequence[Poly], height: int) -> None:
    b = apery_numbers(height)
    polar_checks = 0
    for h in range(1, height + 1):
        assert numerators[h].degree() == 3 * (h - 1)
        reflected = Poly(
            expand(numerators[h].as_expr().subs(x, -x - h - 1)),
            x,
            domain=ZZ,
        )
        assert reflected == ((-1) ** (h - 1)) * numerators[h]
        for j in range(1, h + 1):
            q_without_j = prod(m - j for m in range(1, h + 1) if m != j)
            pole_coefficient = int(numerators[h].eval(-j)) // q_without_j**3
            assert pole_coefficient == b[j - 1] * b[h - j]
            polar_checks += 1
    print(
        "exact_structure: "
        f"degrees/reflection h<= {height}; polar_coefficients={polar_checks}"
    )


def check_local_tangent_cones(numerators: Sequence[Poly]) -> None:
    checks = 0
    for h, k in ((2, 3), (3, 5)):
        qh = q_poly(h).as_expr()
        qk = q_poly(k).as_expr()
        F = expand(
            numerators[h].as_expr() * qk.subs(x, y) ** 3
            - numerators[k].as_expr().subs(x, y) * qh**3
        )
        derivatives = {
            (x_order, y_order): diff(F, x, x_order, y, y_order)
            for x_order in range(4)
            for y_order in range(4 - x_order)
        }
        for j in range(1, h + 1):
            qh_without = prod(m - j for m in range(1, h + 1) if m != j)
            for i in range(1, k + 1):
                qk_without = prod(m - i for m in range(1, k + 1) if m != i)
                expected_u3 = -int(numerators[k].eval(-i)) * qh_without**3
                expected_v3 = int(numerators[h].eval(-j)) * qk_without**3
                for orders, derivative in derivatives.items():
                    x_order, y_order = orders
                    value = int(derivative.subs({x: -j, y: -i}))
                    total_order = x_order + y_order
                    if total_order < 3:
                        assert value == 0
                    elif orders == (3, 0):
                        assert value == factorial(3) * expected_u3
                    elif orders == (0, 3):
                        assert value == factorial(3) * expected_v3
                    else:
                        assert value == 0
                checks += 1
    print(f"local_tangent_cones: pairs=(2,3),(3,5); pole_pairs={checks}")


def poly_add_mod(first: List[int], second: List[int], prime: int) -> List[int]:
    result = [0] * max(len(first), len(second))
    for index, value in enumerate(first):
        result[index] = (result[index] + value) % prime
    for index, value in enumerate(second):
        result[index] = (result[index] + value) % prime
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_scale_mod(values: List[int], scalar: int, prime: int) -> List[int]:
    return [(scalar * value) % prime for value in values]


def poly_mul_mod(first: List[int], second: List[int], prime: int) -> List[int]:
    result = [0] * (len(first) + len(second) - 1)
    for i, first_value in enumerate(first):
        for j, second_value in enumerate(second):
            result[i + j] = (result[i + j] + first_value * second_value) % prime
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def interpolate_consecutive(values: Sequence[int], prime: int) -> Poly:
    """Interpolate f(0),...,f(d) in F_p using the binomial basis."""

    differences = [value % prime for value in values]
    result = [0]
    basis = [1]  # binomial(T, 0)
    for order in range(len(values)):
        result = poly_add_mod(
            result,
            poly_scale_mod(basis, differences[0], prime),
            prime,
        )
        differences = [
            (differences[index + 1] - differences[index]) % prime
            for index in range(len(differences) - 1)
        ]
        if order + 1 < len(values):
            basis = poly_scale_mod(
                poly_mul_mod(basis, [(-order) % prime, 1], prime),
                pow(order + 1, -1, prime),
                prime,
            )
    expression = sum(coefficient * T**degree for degree, coefficient in enumerate(result))
    return Poly(expression, T, modulus=prime)


def exact_small_critical_values(numerators: Sequence[Poly]) -> Dict[int, Poly]:
    exact: Dict[int, Poly] = {}
    for h in range(2, 5):
        qh = q_poly(h)
        Ah = numerators[h].diff() * qh - 3 * numerators[h] * qh.diff()
        resultant = Poly(
            (numerators[h] - Poly(T * qh.as_expr() ** 3, x, domain=ZZ[T])).resultant(Ah),
            T,
            domain=ZZ,
        )
        assert resultant.degree() == 4 * h - 4
        exact[h] = resultant
    print("exact_resultant_baseline: C_h over ZZ for h=2,3,4")
    return exact


def canonical_digest(polynomials: Dict[int, Poly], prime: int) -> str:
    records = []
    for h in sorted(polynomials):
        polynomial = polynomials[h]
        coefficients = [
            str(int(polynomial.nth(degree)) % prime)
            for degree in range(polynomial.degree() + 1)
        ]
        records.append(f"{h}:" + ",".join(coefficients))
    return hashlib.sha256(";".join(records).encode("ascii")).hexdigest()


def critical_value_certificate(
    height: int,
    prime: int,
    exact_small: Dict[int, Poly],
) -> Dict[int, Poly]:
    assert prime > 4 * height + 10
    numerators = mod_numerators(height, prime)
    critical_values: Dict[int, Poly] = {}
    qh = Poly(1, x, modulus=prime)
    for h in range(1, height + 1):
        qh *= Poly(x + h, x, modulus=prime)
        if h == 1:
            continue
        Dh = qh**3
        Ah = numerators[h].diff() * qh - 3 * numerators[h] * qh.diff()
        expected_degree = 4 * h - 4
        # At T=0 the degree in x drops by three.  Calling resultant after
        # specializing would then omit the corresponding power of lc(A_h),
        # so interpolate at the nonzero points 1,...,d+1 and shift back.
        values = [
            int((numerators[h] - value * Dh).resultant(Ah)) % prime
            for value in range(1, expected_degree + 2)
        ]
        shifted_Ch = interpolate_consecutive(values, prime)
        Ch = Poly(expand(shifted_Ch.as_expr().subs(T, T - 1)), T, modulus=prime)
        assert Ch.degree() == expected_degree
        assert all(int(Ch.nth(degree)) % prime == 0 for degree in range(1, expected_degree, 2))
        assert int(Ch.eval(0)) % prime != 0
        assert Ch.gcd(Ch.diff()).degree() == 0
        assert numerators[h].gcd(numerators[h].diff()).degree() == 0
        repeated_part = Ch.gcd(Ch.diff())
        radical = Ch.exquo(repeated_part)
        repeated_roots = radical.gcd(repeated_part)
        simple_part = radical.exquo(repeated_roots)
        if int(simple_part.eval(0)) % prime == 0:
            simple_part = simple_part.exquo(Poly(T, T, modulus=prime))
        simple_nonzero_degree = simple_part.degree()
        assert simple_nonzero_degree >= 2 * h - 1
        for value in (expected_degree + 2, expected_degree + 3, expected_degree + 7):
            direct = int((numerators[h] - value * Dh).resultant(Ah)) % prime
            assert int(Ch.eval(value)) % prime == direct
        if h in exact_small:
            assert Poly(exact_small[h].as_expr(), T, modulus=prime) == Ch
        critical_values[h] = Ch
        short_hash = canonical_digest({h: Ch}, prime)[:16]
        print(
            f"mod_certificate p={prime} h={h}: "
            f"deg={Ch.degree()} C0={int(Ch.eval(0)) % prime} "
            f"squarefree=yes simple_nonzero={simple_nonzero_degree} "
            f"hash={short_hash}"
        )

    bad_pairs: List[Tuple[int, int, int]] = []
    for h in range(2, height + 1):
        for k in range(h + 1, height + 1):
            common_degree = critical_values[h].gcd(critical_values[k]).degree()
            if common_degree != 0:
                bad_pairs.append((h, k, common_degree))
    assert not bad_pairs, bad_pairs
    pair_count = (height - 1) * (height - 2) // 2
    print(
        f"mod_certificate p={prime}: pairwise_coprime_pairs={pair_count}; "
        f"global_sha256={canonical_digest(critical_values, prime)}"
    )
    return critical_values


def factor_degrees_mod(polynomial: Poly, prime: int) -> List[int] | None:
    if int(polynomial.LC()) % prime == 0:
        return None
    _, factors = factor_list(polynomial.as_expr(), modulus=prime)
    if any(exponent != 1 for _, exponent in factors):
        return None
    return sorted(Poly(factor, x, modulus=prime).degree() for factor, _ in factors)


def check_frobenius_specializations(numerators: Sequence[Poly]) -> None:
    """Independent S_(3h) certificates for h<=6 at the specialization T=1."""

    for h in range(2, 7):
        degree = 3 * h
        specialized = numerators[h] - q_poly(h) ** 3
        _, rational_factors = factor_list(specialized.as_expr())
        assert len(rational_factors) == 1
        assert rational_factors[0][1] == 1
        assert Poly(rational_factors[0][0], x, domain=ZZ).degree() == degree

        almost_full = None
        transposition_power = None
        for prime in primerange(19, 500):
            cycle_type = factor_degrees_mod(specialized, prime)
            if cycle_type is None:
                continue
            if almost_full is None and cycle_type == [1, degree - 1]:
                almost_full = (prime, cycle_type)
            if (
                transposition_power is None
                and cycle_type.count(2) == 1
                and all(length % 2 == 1 for length in cycle_type if length != 2)
            ):
                transposition_power = (prime, cycle_type)
            if almost_full is not None and transposition_power is not None:
                break
        assert almost_full is not None
        assert transposition_power is not None
        print(
            f"frobenius h={h} T=1: irreducible_Q=yes; "
            f"p={almost_full[0]} type={almost_full[1]}; "
            f"p={transposition_power[0]} type={transposition_power[1]}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=16)
    parser.add_argument("--primes", type=int, nargs="+", default=[1009, 65537])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    assert args.height >= 6
    numerators = integer_numerators(args.height)
    check_exact_structure(numerators, args.height)
    check_local_tangent_cones(numerators)
    exact_small = exact_small_critical_values(numerators)
    check_frobenius_specializations(numerators)
    for prime in args.primes:
        critical_value_certificate(args.height, prime, exact_small)
    print("ALL_IRRED_CERTIFICATES_OK")


if __name__ == "__main__":
    main()
