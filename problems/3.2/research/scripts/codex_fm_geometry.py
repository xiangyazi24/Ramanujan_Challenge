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


def check_apery_laurent_model() -> None:
    u, v, w = sp.symbols("u v w")
    lam = (u + v) * (w + 1) * (u + v + w) * (v + w + 1) / (u * v * w)
    for exponent in range(7):
        expanded = sp.expand(lam**exponent)
        constant_term = expanded.coeff(u, 0).coeff(v, 0).coeff(w, 0)
        assert constant_term == apery(exponent)
    print("VERIFIED CT Lambda_A^n=A_n for 0<=n<=6")


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


def check_cfvz_rational_pullback() -> None:
    """Check the rational identity after clearing its fixed denominator.

    Lucas gives A_p(phi(x))=(1+x)^(1-p) H_p(x)^2.  Clearing the
    denominator turns this into a polynomial identity of degree at most
    2p-2, so the check is not restricted to F_p-rational evaluations.
    """
    for prime in [5, 7, 11, 13, 17, 19, 23, 29, 37]:
        hp = [franel(n) % prime for n in range(prime)]
        ap = [apery(n) % prime for n in range(prime)]
        left = [0] * (2 * prime - 1)
        for n, a_value in enumerate(ap):
            for j in range(n + 1):
                first = comb(n, j) * pow(-8, j, prime)
                for k in range(prime - n):
                    degree = n + j + k
                    left[degree] = (
                        left[degree]
                        + a_value * first * comb(prime - 1 - n, k)
                    ) % prime
        right = [0] * (2 * prime - 1)
        for i, left_value in enumerate(hp):
            for j, right_value in enumerate(hp):
                right[i + j] = (right[i + j] + left_value * right_value) % prime
        assert left == right
    print(
        "VERIFIED A_p(phi)=(1+x)^(1-p)H_p^2 in F_p(x) "
        "for p=5,7,11,13,17,19,23,29,37"
    )


def check_explicit_elliptic_family() -> None:
    """Check the Weierstrass Franel family and its Hasse congruence.

    The family is E_u: y^2 + (1-2u)xy + u^2 y = x^3.  It is the
    Tate-normal-form family with parameter z=27u^2/(1-2u)^3 after scaling.
    Unlike that rational presentation, this model remains valid at u=1/2.
    """
    u, v = sp.symbols("u v")
    a1 = 1 - 2 * u
    a3 = u**2
    b2 = a1**2
    b4 = a1 * a3
    b6 = a3**2
    c4 = sp.expand(b2**2 - 24 * b4)
    discriminant = sp.factor(-8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6)
    assert discriminant == u**6 * (u + 1) ** 2 * (1 - 8 * u)
    assert c4.subs(u, 0) != 0
    assert c4.subs(u, -1) != 0
    assert c4.subs(u, sp.Rational(1, 8)) != 0

    # For x=v^-2 X and y=v^-3 Y at infinity, Delta and c4 acquire
    # factors v^12 and v^4.  Their orders are 3 and 0, respectively.
    delta_at_infinity = sp.factor(v**12 * discriminant.subs(u, 1 / v))
    c4_at_infinity = sp.factor(v**4 * c4.subs(u, 1 / v))
    assert sp.limit(delta_at_infinity / v**3, v, 0) != 0
    assert sp.limit(c4_at_infinity, v, 0) != 0
    print("VERIFIED Delta(E_u)=u^6(1+u)^2(1-8u) and fiber types I6,I2,I1,I3")

    checked = 0
    for prime in primes_up_to(101):
        if prime < 5:
            continue
        hp = [franel(n) % prime for n in range(prime)]
        singular = {0, prime - 1, pow(8, -1, prime)}
        for parameter in range(prime):
            if parameter in singular:
                continue
            local_a1 = (1 - 2 * parameter) % prime
            local_a3 = parameter * parameter % prime
            points = 1  # the section at infinity
            for x_value in range(prime):
                linear_y = (local_a1 * x_value + local_a3) % prime
                quadratic_discriminant = (
                    linear_y * linear_y + 4 * x_value**3
                ) % prime
                points += 1 + legendre(quadratic_discriminant, prime)
            frobenius_trace = prime + 1 - points
            assert frobenius_trace % prime == evaluate(hp, parameter, prime)
            checked += 1
    print(
        "VERIFIED a_p(E_u)=H_p(u) mod p at every smooth u over all primes "
        f"5<=p<=101 ({checked} fibers, including u=1/2)"
    )


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
    check_apery_laurent_model()
    check_cover_discriminant()
    check_cfvz_rational_pullback()
    check_explicit_elliptic_family()
    check_lucas_dwork()
    check_toric_hasse_point_count()
    check_pushforward_and_mellin()


if __name__ == "__main__":
    main()
