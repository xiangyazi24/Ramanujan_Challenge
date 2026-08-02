#!/usr/bin/env python3
"""Interpolate a 4F3-to-P2.5 gauge from its exact initial projection."""

import sympy as sp


x = sp.symbols("x")


def challenge(n):
    n = sp.sympify(n)
    entries = [
        (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
        384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
        -(480*n**4+4980*n**3+19210*n**2+32690*n+20730),
        (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
        (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
        (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
        (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
        (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240),
    ]
    delta = -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
    return sp.Matrix(3, 3, entries) / delta


def source(n):
    n = sp.sympify(n)
    d0 = (n+1)**2*(2*n+3)**3
    d1 = (n+1)**3*(2*n+3)**4
    entries = [
        (136*n**5+1864*n**4+10080*n**3+26561*n**2+33582*n+15877)/d0,
        -(384*n**6+5904*n**5+38064*n**4+130814*n**3+250809*n**2+251505*n+101045)/d1,
        -(288*n**4+2868*n**3+10458*n**2+16218*n+8718)/d1,
        -(96*n**6+1788*n**5+13634*n**4+54220*n**3+117726*n**2+130736*n+56800)/d0,
        (272*n**7+5528*n**6+47756*n**5+227184*n**4+641473*n**3+1070450*n**2+969568*n+362144)/d1,
        (192*n**5+2892*n**4+16890*n**3+47388*n**2+62976*n+30912)/d1,
        (64*n**7+1508*n**6+15082*n**5+82562*n**4+265482*n**3+497398*n**2+497104*n+200800)/d0,
        -(192*n**8+4760*n**7+51588*n**6+318316*n**5+1218928*n**4+2954291*n**3+4402618*n**2+3659840*n+1282592)/d1,
        (16*n**7+136*n**6-792*n**5-14756*n**4-76886*n**3-192012*n**2-232512*n-106944)/d1,
    ]
    twist = sp.Rational(1, 2) * (2*n+3)/(n+1)
    return twist * sp.Matrix(3, 3, entries)


lower = [
    [sp.Rational(391, 28), -sp.Rational(6667, 448), sp.Rational(1665, 448)],
    [-sp.Rational(2801, 56), sp.Rational(47737, 896), -sp.Rational(11915, 896)],
]


def fit(values, max_degree=12, holdout=5):
    data = [(sp.Integer(index), value) for index, value in enumerate(values)]
    for numerator_degree in range(max_degree + 1):
        for denominator_degree in range(max_degree + 1):
            used = numerator_degree + denominator_degree + 1
            if used + holdout > len(data):
                continue
            candidate = sp.cancel(sp.rational_interpolate(
                data[:used], numerator_degree, X=x))
            if all(sp.cancel(candidate.subs(x, point) - value) == 0
                   for point, value in data[used:]):
                return numerator_degree, denominator_degree, sp.factor(candidate)
    return None


if __name__ == "__main__":
    transitions = [(challenge(index).inv(), source(index)) for index in range(28)]
    for first_row in ([1, 0, 0], [0, 1, 0], [0, 0, 1]):
        gauge = sp.Matrix([first_row, *lower])
        values = [[[] for _ in range(3)] for _ in range(3)]
        for index in range(28):
            for row in range(3):
                for column in range(3):
                    values[row][column].append(gauge[row, column])
            gauge = transitions[index][0] * gauge * transitions[index][1]
        fits = [[fit(values[row][column]) for column in range(3)]
                for row in range(3)]
        print("first row", first_row)
        print(fits)
