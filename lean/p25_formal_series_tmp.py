#!/usr/bin/env python3
"""Temporary exact formal-series calculations for the P2.5 error cone."""

from fractions import Fraction as F
from math import comb
import sys

import sympy as s

hS, AS, BS = s.symbols("h A B")
nS = hS - 2
Q = s.Rational
P0 = s.Matrix([
    [(2*nS+5)*(nS+3)**2*(136*nS**4+1424*nS**3+5548*nS**2+9551*nS+6141),
     384*nS**6+6384*nS**5+44168*nS**4+162698*nS**3+336377*nS**2+369933*nS+169011,
     480*nS**4+4980*nS**3+19210*nS**2+32690*nS+20730],
    [(nS+2)**2*(nS+3)**2*(4*nS+10)*(48*nS**3+386*nS**2+1017*nS+879),
     (nS+2)**2*(272*nS**5+3848*nS**4+21732*nS**3+61184*nS**2+85761*nS+47808),
     (nS+2)**2*(320*nS**3+2540*nS**2+6610*nS+5640)],
    [(4*nS+10)*(nS+2)**2*(nS+3)**2*(32*nS**4+302*nS**3+1037*nS**2+1530*nS+813),
     (nS+2)**2*(192*nS**6+2984*nS**5+19116*nS**4+64452*nS**3+120256*nS**2+117279*nS+46476),
     (nS+2)**2*(16*nS**5+408*nS**4+2912*nS**3+8884*nS**2+12254*nS+6240)],
])
ZS = BS + 4*AS
xnS = Q(5, 4)*hS + AS
ynS = 2*hS*hS + Q(8, 3)*AS*hS + ZS
hpS = hS + 1
dS = s.expand(hS**4*P0[0, 0] - hS*xnS*P0[1, 0] - ynS*P0[2, 0])
e1S = s.expand(-hS**4*P0[0, 1] + hS*xnS*P0[1, 1] + ynS*P0[2, 1])
e2S = s.expand(-hS**4*P0[0, 2] + hS*xnS*P0[1, 2] + ynS*P0[2, 2])
naS = s.expand(hpS**3*e1S - Q(5, 4)*hpS*dS)
nzS = s.expand(hpS**2*(hpS**2*(e2S-Q(8, 3)*e1S)+Q(4, 3)*dS))
nbS = s.expand(nzS - 4*naS)


def affine_coefficients(expr):
    polynomial = s.Poly(expr, AS, BS)
    result = []
    for monomial in (1, AS, BS):
        coefficient = s.Poly(polynomial.coeff_monomial(monomial), hS)
        result.append({int(k[0]): F(int(v.p), int(v.q))
                       for k, v in coefficient.terms()})
    return result


D, NA, NB = map(affine_coefficients, (dS, naS, nbS))
LO = -20
HI = 40


def add(*serieses):
    result = {}
    for series in serieses:
        for exponent, coefficient in series.items():
            result[exponent] = result.get(exponent, F(0)) + coefficient
    return {k: v for k, v in result.items() if v and LO <= k <= HI}


def scale(series, coefficient):
    return {k: v*coefficient for k, v in series.items() if v*coefficient}


def multiply(left, right):
    result = {}
    for i, a in left.items():
        for j, b in right.items():
            if LO <= i+j <= HI:
                result[i+j] = result.get(i+j, F(0)) + a*b
    return {k: v for k, v in result.items() if v}


def power(series, exponent):
    result = {0: F(1)}
    for _ in range(exponent):
        result = multiply(result, series)
    return result


def shift(series, exponent):
    return {i+exponent: v for i, v in series.items()
            if LO <= i+exponent <= HI}


def subtract(left, right):
    return add(left, scale(right, -1))


def ordinary(coefficients):
    return {i: value if isinstance(value, s.Basic) else F(value)
            for i, value in enumerate(coefficients) if value != 0}


def polynomial_at_inverse(poly):
    return {-i: value for i, value in poly.items()}


z_over_one_plus_z = {i: F((-1)**(i-1)) for i in range(1, HI+1)}


def shifted_series(coefficients):
    result = {}
    for i, value in enumerate(coefficients):
        coefficient = value if isinstance(value, s.Basic) else F(value)
        result = add(result, scale(power(z_over_one_plus_z, i), coefficient))
    return result


def affine_series(coefficients, a_series, b_series):
    return add(polynomial_at_inverse(coefficients[0]),
               multiply(polynomial_at_inverse(coefficients[1]), a_series),
               multiply(polynomial_at_inverse(coefficients[2]), b_series))


def residual(a_coefficients, b_coefficients):
    a_series = ordinary(a_coefficients)
    b_series = ordinary(b_coefficients)
    d = affine_series(D, a_series, b_series)
    na = affine_series(NA, a_series, b_series)
    nb = affine_series(NB, a_series, b_series)
    return (shift(subtract(multiply(d, shifted_series(a_coefficients)), na), 12),
            shift(subtract(multiply(d, shifted_series(b_coefficients)), nb), 13))


def coefficient_pair(a_coefficients, b_coefficients, exponent):
    ea, eb = residual(a_coefficients, b_coefficients)
    return ea.get(exponent, F(0)), eb.get(exponent, F(0))


def polynomial_add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            result[exponent] = result.get(exponent, F(0)) + coefficient
    return {k: v for k, v in result.items() if v}


def polynomial_scale(polynomial, coefficient):
    return {k: v*coefficient for k, v in polynomial.items() if v*coefficient}


def polynomial_multiply(left, right):
    result = {}
    for i, a in left.items():
        for j, b in right.items():
            result[i+j] = result.get(i+j, F(0)) + a*b
    return {k: v for k, v in result.items() if v}


def corner_objects(b_value, v_value):
    a_series = {0: F(3, 16),
                -1: -F(3, 8)*b_value-F(39, 128), -2: v_value}

    def evaluate(coefficients):
        return polynomial_add(
            coefficients[0], polynomial_multiply(coefficients[1], a_series),
            polynomial_scale(coefficients[2], b_value))

    denominator = evaluate(D)
    numerator_a = evaluate(NA)
    numerator_b = evaluate(NB)
    h_plus_one = {1: F(1), 0: F(1)}
    numerator_v = polynomial_add(
        polynomial_multiply(polynomial_multiply(h_plus_one, h_plus_one),
            polynomial_add(numerator_a,
                polynomial_scale(denominator, -F(3, 16)))),
        polynomial_multiply(h_plus_one,
            polynomial_add(polynomial_scale(numerator_b, F(3, 8)),
                polynomial_scale(denominator, F(39, 128)))))
    return denominator, numerator_b, numerator_v


def corner_objects_w(b_value, w_value):
    a_series = {0: F(3, 16),
                -1: -F(3, 8)*b_value-F(39, 128),
                -2: F(3, 4)*b_value+F(191, 256), -3: w_value}

    def evaluate(coefficients):
        return polynomial_add(
            coefficients[0], polynomial_multiply(coefficients[1], a_series),
            polynomial_scale(coefficients[2], b_value))

    denominator = evaluate(D)
    numerator_a = evaluate(NA)
    numerator_b = evaluate(NB)
    h_plus_one = {1: F(1), 0: F(1)}
    numerator_v = polynomial_add(
        polynomial_multiply(polynomial_multiply(h_plus_one, h_plus_one),
            polynomial_add(numerator_a,
                polynomial_scale(denominator, -F(3, 16)))),
        polynomial_multiply(h_plus_one,
            polynomial_add(polynomial_scale(numerator_b, F(3, 8)),
                polynomial_scale(denominator, F(39, 128)))))
    numerator_w = polynomial_multiply(h_plus_one,
        polynomial_add(numerator_v,
            polynomial_scale(numerator_b, -F(3, 4)),
            polynomial_scale(denominator, -F(191, 256))))
    return denominator, numerator_b, numerator_w


def shifted_polynomial_coefficients(polynomial, start):
    minimum = min(polynomial)
    if minimum < 0:
        polynomial = {k-minimum: v for k, v in polynomial.items()}
    result = {}
    for exponent, coefficient in polynomial.items():
        for j in range(exponent+1):
            result[j] = result.get(j, F(0)) + coefficient*comb(exponent, j)*F(start)**(exponent-j)
    return result


def eventually_nonnegative(polynomial, start):
    return all(value >= 0
               for value in shifted_polynomial_coefficients(polynomial, start).values())


def test_cone(bounds, start_n, verbose=False):
    b_lower, b_upper, v_lower, v_upper = bounds
    success = True
    for b_value in (b_lower, b_upper):
        for v_value in (v_lower, v_upper):
            denominator, numerator_b, numerator_v = corner_objects(b_value, v_value)
            tests = [denominator,
                polynomial_add(numerator_b, polynomial_scale(denominator, -b_lower)),
                polynomial_add(polynomial_scale(denominator, b_upper), polynomial_scale(numerator_b, -1)),
                polynomial_add(numerator_v, polynomial_scale(denominator, -v_lower)),
                polynomial_add(polynomial_scale(denominator, v_upper), polynomial_scale(numerator_v, -1))]
            flags = [eventually_nonnegative(test, start_n+2) for test in tests]
            if verbose:
                print(b_value, v_value, flags,
                      [min(shifted_polynomial_coefficients(test, start_n+2).values())
                       for test in tests])
            success = success and all(flags)
    return success


def test_cone_w(bounds, start_n, verbose=False):
    b_lower, b_upper, w_lower, w_upper = bounds
    success = True
    for b_value in (b_lower, b_upper):
        for w_value in (w_lower, w_upper):
            denominator, numerator_b, numerator_w = corner_objects_w(b_value, w_value)
            tests = [denominator,
                polynomial_add(numerator_b, polynomial_scale(denominator, -b_lower)),
                polynomial_add(polynomial_scale(denominator, b_upper), polynomial_scale(numerator_b, -1)),
                polynomial_add(numerator_w, polynomial_scale(denominator, -w_lower)),
                polynomial_add(polynomial_scale(denominator, w_upper), polynomial_scale(numerator_w, -1))]
            flags = [eventually_nonnegative(test, start_n+2) for test in tests]
            if verbose:
                print(b_value, w_value, flags,
                      [min(shifted_polynomial_coefficients(test, start_n+2).values())
                       for test in tests])
            success = success and all(flags)
    return success


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "relations":
    t, a2, b1 = s.symbols("t a2 b1")
    a1 = (-39-48*t)/128
    equations = coefficient_pair([F(3, 16), a1, a2], [t, b1], 2)
    print("equations", *(s.factor(equation) for equation in equations))
    print("solution", s.solve(equations, (a2, b1), dict=True))


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "cone2":
    b_intervals = [
        (F(-6, 5), F(-9, 8)),
        (F(-19, 16), F(-9, 8)),
        (F(-7, 6)-F(1, 32), F(-7, 6)+F(1, 32)),
        (F(-7, 6)-F(1, 64), F(-7, 6)+F(1, 64)),
        (F(-7, 6)-F(1, 128), F(-7, 6)+F(1, 128)),
    ]
    w_intervals = [
        (F(0), F(1, 12)),
        (F(0), F(1, 16)),
        (F(0), F(1, 20)),
        (F(1, 100), F(1, 12)),
        (F(1, 100), F(1, 20)),
        (F(1, 128), F(3, 64)),
    ]
    for start in (3, 4, 5, 8, 10, 20, 50, 100, 200):
        for b_interval in b_intervals:
            for w_interval in w_intervals:
                bounds = b_interval + w_interval
                if test_cone_w(bounds, start):
                    print("FOUND", start, bounds)
                    raise SystemExit
    print("diagnostic")
    test_cone_w((F(-6, 5), F(-9, 8), F(0), F(1, 12)), 3, True)


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "cone":
    b_intervals = [
        (F(-6, 5), F(-9, 8)),
        (F(-19, 16), F(-9, 8)),
        (F(-7, 6)-F(1, 32), F(-7, 6)+F(1, 32)),
        (F(-7, 6)-F(1, 64), F(-7, 6)+F(1, 64)),
        (F(-7, 6)-F(1, 128), F(-7, 6)+F(1, 128)),
    ]
    v_intervals = [
        (F(-1, 6), F(-1, 16)),
        (F(-3, 20), F(-1, 10)),
        (F(-9, 64), F(-7, 64)),
        (F(-17, 128), F(-15, 128)),
        (F(-1, 7), F(-1, 10)),
        (F(-9, 64), F(-15, 128)),
    ]
    for start in (3, 4, 5, 8, 10, 20, 50, 100, 200):
        for b_interval in b_intervals:
            for v_interval in v_intervals:
                bounds = b_interval + v_interval
                if test_cone(bounds, start):
                    print("FOUND", start, bounds)
                    raise SystemExit
    print("diagnostic")
    test_cone((F(-6, 5), F(-9, 8), F(-1, 6), F(-1, 16)), 3, True)


if __name__ == "__main__" and len(sys.argv) == 1:
    a = [F(3, 16)]
    b = []
    print("a0", a[0], "residual", coefficient_pair(a, [0], 0))
    for exponent in range(1, 10):
        base_a = a + [F(0)]
        base_b = b + [F(0)]
        constant = coefficient_pair(base_a, base_b, exponent)
        at_a = coefficient_pair(a + [F(1)], base_b, exponent)
        at_b = coefficient_pair(base_a, b + [F(1)], exponent)
        column_a = (at_a[0]-constant[0], at_a[1]-constant[1])
        column_b = (at_b[0]-constant[0], at_b[1]-constant[1])
        determinant = column_a[0]*column_b[1] - column_a[1]*column_b[0]
        print("step", exponent, "constant", constant,
              "columns", column_a, column_b, "det", determinant)
        if determinant == 0:
            break
        next_a = (-constant[0]*column_b[1]+constant[1]*column_b[0])/determinant
        next_b = (-column_a[0]*constant[1]+column_a[1]*constant[0])/determinant
        a.append(next_a)
        b.append(next_b)
        print("  a", exponent, next_a, "b", exponent-1, next_b,
              "residual", coefficient_pair(a, b, exponent))
