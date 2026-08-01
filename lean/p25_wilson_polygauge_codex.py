#!/usr/bin/env python3
"""Search polynomial intertwiners between the P2.5 and Wilson cocycles."""

import sympy as s

n = s.symbols("n")

M = s.Matrix([
    [(-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
     384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
     -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)],
    [(n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
     (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
     (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)],
    [(-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
     (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
     (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)],
])
dA = -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

d0 = 8*(n+3)**2*(4*n+11)**2
dT = d0*(4*n+13)**2
Tnum = s.Matrix([
    [(4*n+13)**2*(4*n+9)**2*(40*n**2+228*n+325),
     (4*n+13)**2*(1536*n**4+16512*n**3+66496*n**2+118896*n+79641)],
    [(4*n+9)**2*(1536*n**4+18432*n**3+82880*n**2+165504*n+123841),
     59392*n**6+1012736*n**5+7184384*n**4+27140352*n**3+
     57583336*n**2+65059404*n+30580677],
])


def coefficient_matrix(expressions, unknowns):
    equations = []
    for expression in expressions:
        polynomial = s.Poly(s.expand(expression), n)
        equations.extend(polynomial.all_coeffs())
    return s.linear_eq_to_matrix(equations, unknowns)[0]


def quotient_search(max_degree=14):
    for degree in range(max_degree + 1):
        unknowns = s.symbols(f"c0:{6*(degree+1)}")
        H = s.Matrix(3, 2, lambda i, j: sum(
            unknowns[(2*i+j)*(degree+1)+k]*n**k for k in range(degree+1)
        ))
        Hnext = H.subs(n, n+1)
        residual = dT*M*Hnext-dA*H*Tnum
        linear = coefficient_matrix(list(residual), unknowns)
        nullity = len(unknowns)-linear.rank()
        print("quotient degree", degree, "shape", linear.shape,
              "nullity", nullity, flush=True)
        if nullity:
            basis = linear.nullspace()
            for solution in basis:
                gauge = H.subs(dict(zip(unknowns, solution)))
                print("rank", gauge.rank(), flush=True)
                for row in gauge.tolist():
                    print([s.factor(value) for value in row], flush=True)
                if gauge.rank() == 2:
                    return gauge
    return None


print("found", quotient_search() is not None, flush=True)
