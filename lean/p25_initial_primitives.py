#!/usr/bin/env python3
"""Exact a-primitives after a=pq and t=2v/(1+v^2)."""
import sympy as s

a,t=s.symbols('a t', positive=True)
M=[a**5*(1-a)**2*t**(3+k)/(s.Integer(2)**(4+k)*(a+t)**(4+k))
   for k in range(3)]
for name,cs in [('den',(33750,126000,180000)),
                ('num',(30921,115408,164800))]:
    f=s.factor(sum(c*m for c,m in zip(cs,M)))
    F=s.integrate(f,a,risch=True)
    print(name,'integrand=',f)
    print(name,'primitive=',s.collect(s.factor(F),s.log(a+t)))
    print('check=',s.factor(s.diff(F,a)-f))
