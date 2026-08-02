#!/usr/bin/env python3
import sympy as s

G, L = s.symbols("G L")
U = s.Matrix([
    [1, s.Rational(152,5), s.Rational(195477,175)],
    [0, s.Rational(1723,90), s.Rational(1963751,2800)],
    [0, -s.Rational(143,18), -s.Rational(165201,560)],
])
F = s.Matrix([
    150*G-128*L-s.Rational(146,3),
    s.Rational(24745,4)*G-s.Rational(14624,3)*L-s.Rational(823511,360),
    s.Rational(7225281,32)*G-s.Rational(886784,5)*L-s.Rational(2818419551,33600),
])
Y = U.T.inv()*F
Z = s.diag(1,-1,1)*Y / 96
for name, vec in (("Y/sqrtpi",Y),("Z=triple",Z)):
    print(name)
    for x in vec:
        print(s.collect(s.factor(x),[G,L]))
Apos=s.Matrix([[30921,32972,8240],[33750,36000,9000]])
print("Apos Z", [s.factor(x) for x in Apos*Z])
