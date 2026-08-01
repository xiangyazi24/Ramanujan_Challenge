#!/usr/bin/env python3
"""Finite-field and formal-series checks for the Franel--Mellin object.

This script deliberately checks only algebraic identities.  It does not claim
that a mod-p Hasse--Witt value is, by itself, an l-adic trace function.
"""

from math import comb

import sympy as sp


def primes_up_to(bound: int) -> list[int]:
    answer = []
    for candidate in range(2, bound + 1):
        if all(candidate % divisor for divisor in range(2, int(candidate**0.5) + 1)):
            answer.append(candidate)
    return answer


def franel(index: int) -> int:
    return sum(comb(index, k) ** 3 for k in range(index + 1))


def apery(index: int) -> int:
    return sum(
        comb(index, k) ** 2 * comb(index + k, k) ** 2
        for k in range(index + 1)
    )


def evaluate(coefficients: list[int], argument: int, prime: int) -> int:
    value = 0
    for coefficient in reversed(coefficients):
        value = (value * argument + coefficient) % prime
    return value


def legendre(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1


def check_characteristic_zero_pullback() -> None:
    x = sp.symbols("x")
    cutoff = 24
    h = sum(sp.Integer(franel(n)) * x**n for n in range(cutoff))
    f = sum(sp.Integer(apery(n)) * sp.symbols("t") ** n for n in range(cutoff))
    phi = x * (1 - 8 * x) / (1 + x)
    difference = sp.series(f.subs(sp.symbols("t"), phi) - (1 + x) * h**2, x, 0, cutoff)
    assert difference.removeO().expand() == 0
    print(f"VERIFIED characteristic-zero pullback through O(x^{cutoff})")


def check_cover_discriminant() -> None:
    x, t = sp.symbols("x t")
    cover_polynomial = 8 * x**2 + (t - 1) * x + t
    discriminant = sp.discriminant(cover_polynomial, x)
    assert sp.expand(discriminant) == t**2 - 34 * t + 1
    phi = x * (1 - 8 * x) / (1 + x)
    critical_numerator = sp.factor(sp.together(sp.diff(phi, x))).as_numer_denom()[0]
    assert sp.expand(critical_numerator) == -8 * x**2 - 16 * x + 1
    q_pullback = sp.factor((phi**2 - 34 * phi + 1))
    assert sp.simplify(q_pullback - ((1 - 16 * x - 8 * x**2) / (1 + x)) ** 2) == 0
    print("VERIFIED cover discriminant q(t)=t^2-34t+1 and its square pullback")


def check_lucas_dwork() -> None:
    for prime in [5, 7, 11, 13, 17, 19, 29, 37]:
        coefficients = [franel(n) % prime for n in range(4 * prime)]
        truncation = coefficients[:prime]
        product = [0] * (4 * prime)
        for left, left_value in enumerate(truncation):
            for right in range(4):
                degree = left + prime * right
                if degree < len(product):
                    product[degree] = left_value * coefficients[right] % prime
        assert product == coefficients
    print("VERIFIED h(x)=H_p(x)h(x)^p mod p through degree 4p-1 for p=5,7,11,13,17,19,29,37")


def check_toric_hasse_point_count() -> None:
    # The reflexive-hexagon compactification contributes its six rational
    # boundary points.  The resulting smooth curve is elliptic away from the
    # four Picard--Fuchs singular parameters.
    checked = 0
    for prime in [5, 7, 11, 13, 17, 19, 23]:
        hp = [franel(n) % prime for n in range(prime)]
        for parameter in range(1, prime):
            if (parameter + 1) * (8 * parameter - 1) % prime == 0:
                continue
            torus_points = 0
            for u in range(1, prime):
                for v in range(1, prime):
                    lam = (1 + u) * (1 + v) * (1 + pow(u * v, -1, prime))
                    if (1 - parameter * lam) % prime == 0:
                        torus_points += 1
            compact_points = torus_points + 6
            frobenius_trace = prime + 1 - compact_points
            assert frobenius_trace % prime == evaluate(hp, parameter, prime)
            checked += 1
    print(f"VERIFIED Franel toric Hasse/point-count congruence at {checked} smooth fibers")


def check_pushforward_and_mellin() -> None:
    tested_pairs = 0
    for prime in [5, 7, 11, 13, 17, 19, 23, 29, 37]:
        hp = [franel(n) % prime for n in range(prime)]
        ap = [apery(n) % prime for n in range(prime)]
        pushforward = [0] * prime
        for x_value in range(prime):
            if (1 + x_value) % prime == 0:
                continue
            phi_value = x_value * (1 - 8 * x_value) * pow(1 + x_value, -1, prime) % prime
            h_value = evaluate(hp, x_value, prime)
            pushforward[phi_value] = (pushforward[phi_value] + h_value**2) % prime

        for t_value in range(prime):
            q_value = (t_value * t_value - 34 * t_value + 1) % prime
            a_value = evaluate(ap, t_value, prime)
            assert pushforward[t_value] == (1 + legendre(q_value, prime)) * a_value % prime
            virtual_trace = (-pushforward[t_value] + legendre(q_value, prime) * a_value) % prime
            assert virtual_trace == -a_value % prime

        for exponent in range(1, prime - 1):
            mellin = sum(
                (-pushforward[t_value]
                 + legendre(t_value * t_value - 34 * t_value + 1, prime)
                 * evaluate(ap, t_value, prime))
                * pow(t_value, prime - 1 - exponent, prime)
                for t_value in range(1, prime)
            ) % prime
            assert mellin == ap[exponent]
            tested_pairs += 1
    print(
        "VERIFIED pushforward=(1+chi(q))A_p, virtual cancellation=-A_p, "
        f"and Mellin=b_r for {tested_pairs} (p,r) pairs"
    )


def main() -> None:
    check_characteristic_zero_pullback()
    check_cover_discriminant()
    check_lucas_dwork()
    check_toric_hasse_point_count()
    check_pushforward_and_mellin()


if __name__ == "__main__":
    main()
