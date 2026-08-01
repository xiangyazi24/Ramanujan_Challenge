#!/usr/bin/env python3
"""Check the proposed moving cross-product cone exactly."""

import sympy as s
import sys

n = s.symbols("n", nonnegative=True, integer=True)
P = s.Matrix([
    [(2*n+5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
     384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
     480*n**4+4980*n**3+19210*n**2+32690*n+20730],
    [(n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
     (n+2)**2*(272*n**5+3848*n**4+21732*n**3+61184*n**2+85761*n+47808),
     (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)],
    [(4*n+10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
     (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
     (n+2)**2*(16*n**5+408*n**4+2912*n**3+8884*n**2+12254*n+6240)],
])

if len(sys.argv) > 1 and sys.argv[1] == "limit":
    alpha, gamma = s.symbols("alpha gamma")
    x1 = (s.Rational(5, 4) + alpha / (n + 1)) / (n + 2)**2
    x2 = (2 + s.Rational(1, 2) / (n + 1) + gamma / (n + 1)**2) / (n + 2)**2
    y0 = P[0, 0] - x1 * P[1, 0] - x2 * P[2, 0]
    y1 = -P[0, 1] + x1 * P[1, 1] + x2 * P[2, 1]
    y2 = -P[0, 2] + x1 * P[1, 2] + x2 * P[2, 2]
    alpha_next = (n + 2) * ((n + 3)**2 * y1 / y0 - s.Rational(5, 4))
    gamma_next = (n + 2) * (
        (n + 2) * ((n + 3)**2 * y2 / y0 - 2) - s.Rational(1, 2)
    )
    print("alpha map", s.factor(s.limit(alpha_next, n, s.oo)))
    print("gamma map", s.factor(s.limit(gamma_next, n, s.oo)))
    sys.exit(0)

corners = [(s.Rational(a), s.Rational(b))
           for a in (1, 2) for b in (s.Rational(3, 2), s.Rational(5, 2))]

for a, b in corners:
    x1 = a / (n + 2)**2
    x2 = b / (n + 2)**2
    y0 = P[0, 0] - x1 * P[1, 0] - x2 * P[2, 0]
    y1 = -P[0, 1] + x1 * P[1, 1] + x2 * P[2, 1]
    y2 = -P[0, 2] + x1 * P[1, 2] + x2 * P[2, 2]
    tests = {
        "y0": y0,
        "y1 lower": y1 - y0 / (n + 3)**2,
        "y1 upper": 2 * y0 / (n + 3)**2 - y1,
        "y2 lower": y2 - s.Rational(3, 2) * y0 / (n + 3)**2,
        "y2 upper": s.Rational(5, 2) * y0 / (n + 3)**2 - y2,
    }
    print("corner", a, b)
    for name, value in tests.items():
        numerator, denominator = s.cancel(value).as_numer_denom()
        coefficients = s.Poly(numerator, n).all_coeffs()
        print(name, "positive coefficients", all(c >= 0 for c in coefficients),
              "factor", s.factor(numerator), "denominator", s.factor(denominator))

print("correlated cone")
a_lower = s.Rational(5, 4)
a_upper = s.Rational(5, 4) + 1 / (3 * (n + 1))
next_a_lower = s.Rational(5, 4)
next_a_upper = s.Rational(5, 4) + 1 / (3 * (n + 2))
for label, a in (("lower", a_lower), ("upper", a_upper)):
    for r in (s.Rational(4, 3), s.Rational(5, 3)):
        x1 = a / (n + 2)**2
        x2 = r * x1
        y0 = P[0, 0] - x1 * P[1, 0] - x2 * P[2, 0]
        y1 = -P[0, 1] + x1 * P[1, 1] + x2 * P[2, 1]
        y2 = -P[0, 2] + x1 * P[1, 2] + x2 * P[2, 2]
        tests = {
            "y0": y0,
            "a lower": (n + 3)**2 * y1 - next_a_lower * y0,
            "a upper": next_a_upper * y0 - (n + 3)**2 * y1,
            "r lower": y2 - s.Rational(4, 3) * y1,
            "r upper": s.Rational(5, 3) * y1 - y2,
        }
        print("corner", label, r)
        for name, value in tests.items():
            numerator, denominator = s.cancel(value).as_numer_denom()
            coefficients = s.Poly(numerator, n).all_coeffs()
            print(name, "positive coefficients", all(c >= 0 for c in coefficients),
                  "factor", s.factor(numerator), "denominator", s.factor(denominator))

print("asymptotic rectangle")
for alpha in (s.Rational(1, 6), s.Rational(1, 3)):
    for beta in (s.Rational(0), s.Rational(1, 2)):
        x1 = (s.Rational(5, 4) + alpha / (n + 1)) / (n + 2)**2
        x2 = (2 + beta / (n + 1)) / (n + 2)**2
        y0 = P[0, 0] - x1 * P[1, 0] - x2 * P[2, 0]
        y1 = -P[0, 1] + x1 * P[1, 1] + x2 * P[2, 1]
        y2 = -P[0, 2] + x1 * P[1, 2] + x2 * P[2, 2]
        next_a_lower = s.Rational(5, 4) + s.Rational(1, 6) / (n + 2)
        next_a_upper = s.Rational(5, 4) + s.Rational(1, 3) / (n + 2)
        next_b_lower = 2
        next_b_upper = 2 + s.Rational(1, 2) / (n + 2)
        tests = {
            "y0": y0,
            "alpha lower": (n + 3)**2 * y1 - next_a_lower * y0,
            "alpha upper": next_a_upper * y0 - (n + 3)**2 * y1,
            "beta lower": (n + 3)**2 * y2 - next_b_lower * y0,
            "beta upper": next_b_upper * y0 - (n + 3)**2 * y2,
        }
        print("corner", alpha, beta)
        for name, value in tests.items():
            numerator, denominator = s.cancel(value).as_numer_denom()
            coefficients = s.Poly(numerator, n).all_coeffs()
            print(name, "positive coefficients", all(c >= 0 for c in coefficients),
                  "factor", s.factor(numerator), "denominator", s.factor(denominator))
