#!/usr/bin/env python3
"""Hermite reduction of the two initial one-dimensional integrands."""
import sympy as s

t=s.symbols('t')

data = {
 'den': (
   -s.Rational(375,32)*t**3/(1+t)**3,
   -420*t**6-1050*t**5-410*t**4+795*t**3+681*t**2+83*t-15,
   420*t**7+1260*t**6+900*t**5-660*t**4-1080*t**3-360*t**2),
 'num': (
   -s.Rational(1,192)*t**3/(1+t)**3,
   -866460*t**6-2164710*t**5-842230*t**4+1642725*t**3+1405263*t**2+171157*t-30921,
   866460*t**7+2597940*t**6+1852380*t**5-1365900*t**4-2229480*t**3-742680*t**2),
}

def rational_solution(rhs, extra=0):
    """Solve (1-t^2)y'-ty=rhs for rational y and optional constant rhs shift."""
    for d in range(0,6):
        for deg in range(0,18):
            cs=s.symbols(f'b0:{deg+1}')
            c=s.symbols('c') if extra else s.Integer(0)
            y=sum(cs[i]*t**i for i in range(deg+1))/(1+t)**d
            expr=s.together((1-t**2)*s.diff(y,t)-t*y-rhs+c)
            poly=s.Poly(expr.as_numer_denom()[0],t)
            sol=s.linsolve([co for co in poly.all_coeffs()], (*cs, *((c,) if extra else ())))
            if sol is not s.EmptySet and sol != s.EmptySet:
                for tup in sol:
                    if not any(x.free_symbols & set(cs) for x in tup):
                        sub=dict(zip((*cs, *((c,) if extra else ())),tup))
                        return s.factor(y.subs(sub)), (s.factor(c.subs(sub)) if extra else 0)
    raise RuntimeError('no solution')

for name,(pref,P,Q) in data.items():
    R=s.factor(pref*P)
    S=s.factor(pref*Q)
    B,c=rational_solution(S,extra=1)
    rhsA=s.factor(R+(1-t)*B/t)
    A,_=rational_solution(rhsA,extra=0)
    print(name)
    print('c =',c)
    print('B =',B)
    print('A =',A)
    print('checkL =',s.factor(c+(1-t**2)*s.diff(B,t)-t*B-S))
    print('checkR =',s.factor((1-t**2)*s.diff(A,t)-t*A-(1-t)*B/t-R))
    print('boundaries A0,A1,B0,B1 =',A.subs(t,0),A.subs(t,1),B.subs(t,0),B.subs(t,1))
