#!/usr/bin/env python3
import sympy as s

a,t=s.symbols('a t', positive=True)
js=[]
for k in range(3):
    f=a**5*(1-a)**2*t**(3+k)/(s.Integer(2)**(4+k)*(a+t)**(4+k))
    j=s.integrate(f,(a,0,1))
    js.append(j)
    print('k',k)
    print(s.collect(s.factor(j),s.log((t+1)/t)))

for name, coeff in [('den',(33750,126000,180000)),
                    ('num',(30921,115408,164800))]:
    e=s.factor(sum(c*j for c,j in zip(coeff,js)))
    print(name)
    print(s.collect(e,[s.log(t),s.log(t+1)]))
