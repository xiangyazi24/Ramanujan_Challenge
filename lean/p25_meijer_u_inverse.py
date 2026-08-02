#!/usr/bin/env python3
"""Factor the Meijer companion Krylov inverse and sign functionals."""

import sympy as s

n = s.symbols("n")


def transition(m):
    return s.Matrix([
        [4*(m+2)*(17*m**3+111*m**2+240*m+171)/((m+1)*(m+3)*(2*m+3)*(2*m+5)),
         (m+2)*(24*m**2+101*m+102)/((m+1)*(2*m+3)),
         (m+2)*(2*m+5)*(16*m**2+81*m+90)/(2*(m+1)*(2*m+3))],
        [(96*m**4+780*m**3+2384*m**2+3273*m+1723)/((m+1)*(m+2)*(m+3)*(2*m+3)*(2*m+5)),
         (68*m**3+398*m**2+778*m+523)/(2*(m+1)*(m+2)*(2*m+3)),
         (96*m**4+884*m**3+2970*m**2+4360*m+2403)/(4*(m+1)*(m+2)*(2*m+3))],
        [-5*(24*m**2+117*m+143)/((m+1)*(m+2)*(m+3)*(2*m+3)*(2*m+5)),
         -5*(16*m+41)/(2*(m+1)*(m+2)*(2*m+3)),
         (8*m**3-44*m**2-478*m-801)/(4*(m+1)*(m+2)*(2*m+3))],
    ])


e = s.Matrix([1, 0, 0])
T = transition(n)
U = s.Matrix.hstack(e, T*e, T*transition(n+1)*e)
print("det", s.factor(U.det()))
V = s.simplify(U.inv())
for i in range(3):
    print("coordinate", i)
    for j in range(3):
        print(j, s.factor(V[j, i]))

ratios = [1, 1/((n+1)*(n+2)),
          1/((n+1)*(n+2)**2*(n+3))]
z = s.symbols("z", positive=True)
for j in range(3):
    polynomial = s.factor(sum(V[r, j]*ratios[r]*z**r for r in range(3)))
    numerator, denominator = s.cancel(polynomial).as_numer_denom()
    print("moment polynomial", j, s.factor(numerator), "/", s.factor(denominator))
    if s.degree(numerator, z) == 2:
        print("discriminant", j, s.factor(s.discriminant(numerator, z)))
