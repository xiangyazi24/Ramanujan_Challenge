#!/usr/bin/env python3
"""P2.5: Factor the scalar recurrence and find its minimal form.

The recurrence c_3(n) Q̂(n+3) + c_2(n) Q̂(n+2) + c_1(n) Q̂(n+1) + c_0(n) Q̂(n) = 0
has all coefficients of degree 13.
"""
from sympy import Symbol, Poly, factor, gcd, ZZ, QQ, Rational, prod
from sympy import factorint, isprime
from math import gcd as igcd
from functools import reduce

n = Symbol('n')

# Coefficient arrays (index j = power of n, so c_i(n) = sum a[j] n^j)
c0_coeffs = [-42743162700, -206623731375, -448112471583, -579493151986, -500074412234,
             -304838513875, -135313819947, -44354837964, -10750665744, -1905022784,
             -240100240, -20397440, -1047552, -24576]
c1_coeffs = [8781630505200, 38850314624124, 78557994908508, 96136040496551, 79442239242197,
             46814452218572, 20241514501104, 6502490145168, 1552168938336, 271943188864,
             33995217088, 2871763456, 146952192, 3440640]
c2_coeffs = [-10566229124340, -43764612822972, -82725628159809, -94536939564882, -72904809920709,
             -40082159230086, -16169158004002, -4847446446296, -1080358338832, -176841798272,
             -20670362464, -1634185472, -78342144, -1720320]
c3_coeffs = [146862156672, 610678861056, 1158857071416, 1329423744980, 1029037642166,
             567735994679, 229759169143, 69074560420, 15430450432, 2530117664,
             296032016, 23408000, 1121280, 24576]

def make_poly(coeffs):
    """coeffs[j] is the coefficient of n^j."""
    return sum(c * n**j for j, c in enumerate(coeffs))

c0 = make_poly(c0_coeffs)
c1 = make_poly(c1_coeffs)
c2 = make_poly(c2_coeffs)
c3 = make_poly(c3_coeffs)

print("="*60)
print("1. Factoring each coefficient polynomial over Q[n]")
print("="*60)

for label, poly in [("c0", c0), ("c1", c1), ("c2", c2), ("c3", c3)]:
    f = factor(poly)
    print(f"\n{label}(n) = {f}")

print("\n" + "="*60)
print("2. GCD of all four polynomials")
print("="*60)

g01 = gcd(Poly(c0, n, domain='ZZ'), Poly(c1, n, domain='ZZ'))
g012 = gcd(g01, Poly(c2, n, domain='ZZ'))
g0123 = gcd(g012, Poly(c3, n, domain='ZZ'))
print(f"\ngcd(c0, c1, c2, c3) = {g0123.as_expr()}")

# After dividing by GCD, the reduced recurrence
if g0123.degree() > 0:
    c0r = Poly(c0, n, domain='ZZ') // g0123
    c1r = Poly(c1, n, domain='ZZ') // g0123
    c2r = Poly(c2, n, domain='ZZ') // g0123
    c3r = Poly(c3, n, domain='ZZ') // g0123
    print(f"\nReduced degrees: c0={c0r.degree()}, c1={c1r.degree()}, c2={c2r.degree()}, c3={c3r.degree()}")
    print("\nReduced c0 factored:", factor(c0r.as_expr()))
    print("Reduced c1 factored:", factor(c1r.as_expr()))
    print("Reduced c2 factored:", factor(c2r.as_expr()))
    print("Reduced c3 factored:", factor(c3r.as_expr()))
else:
    print("GCD is a constant — no polynomial common factor.")

print("\n" + "="*60)
print("3. Pairwise GCDs")
print("="*60)
for i, (li, pi) in enumerate([("c0", c0), ("c1", c1), ("c2", c2), ("c3", c3)]):
    for j, (lj, pj) in enumerate([("c0", c0), ("c1", c1), ("c2", c2), ("c3", c3)]):
        if j > i:
            g = gcd(Poly(pi, n, domain='ZZ'), Poly(pj, n, domain='ZZ'))
            if g.degree() > 0:
                print(f"  gcd({li}, {lj}) = {factor(g.as_expr())} (degree {g.degree()})")

print("\n" + "="*60)
print("4. Rational roots of each polynomial")
print("="*60)
from sympy import solve, Rational as R
for label, poly in [("c0", c0), ("c1", c1), ("c2", c2), ("c3", c3)]:
    p = Poly(poly, n, domain='QQ')
    roots = []
    # Check integer roots from -20 to 5
    for r in range(-20, 6):
        if p.eval(r) == 0:
            roots.append(r)
    # Check half-integer roots
    for num in range(-41, 11):
        r = Rational(num, 2)
        if p.eval(r) == 0:
            roots.append(r)
    if roots:
        print(f"  {label}: rational roots at n = {roots}")
    else:
        print(f"  {label}: no small rational roots in [-20, 5] or half-integers")

print("\n" + "="*60)
print("5. Ore factorization test (right factor order 1)")
print("="*60)
print("If L = L_2 ∘ L_1 with L_1 = a_1(n)σ + a_0(n) of degree d,")
print("then c_3(n) = b_2(n)*a_1(n+2) and c_0(n) = b_0(n)*a_0(n).")
print("\nFactoring c_3(n) and c_0(n) to find candidate splits...")

c3_factored = factor(c3)
c0_factored = factor(c0)
print(f"\nc_3(n) factored: {c3_factored}")
print(f"c_0(n) factored: {c0_factored}")

print("\n" + "="*60)
print("6. Checking c_3(n) vs c_3(n+2) factorization")
print("="*60)
# For an Ore right-factor L_1 = a_1(n)σ + a_0(n):
# c_3(n) = b_2(n) * a_1(n+2)
# c_0(n) = b_0(n) * a_0(n)
# So a_1(n+2) | c_3(n), and a_0(n) | c_0(n)
# Let's check: does any shifted factor of c_0 divide c_3?
c0_poly = Poly(c0, n, domain='ZZ')
c3_poly = Poly(c3, n, domain='ZZ')

# Factor c_0 and c_3 completely
from sympy import factor_list
print("\nc_0 factor list:")
fl0 = factor_list(c0, n)
print(f"  Content: {fl0[0]}")
for f, m in fl0[1]:
    print(f"  ({f})^{m}, degree {Poly(f, n).degree()}")

print("\nc_3 factor list:")
fl3 = factor_list(c3, n)
print(f"  Content: {fl3[0]}")
for f, m in fl3[1]:
    print(f"  ({f})^{m}, degree {Poly(f, n).degree()}")

print("\n" + "="*60)
print("7. Check if recurrence is self-adjoint or adjoint-related")
print("="*60)
# Self-adjoint: c_0(n) * c_3(n-3) = c_3(n) * c_0(n-3)? (modulo sign/normalization)
# Or more precisely: the adjoint operator L* has c_k*(n) = c_{order-k}(n-order+k)
# For order 3: c_0*(n) = c_3(n-3), c_1*(n) = c_2(n-2), etc.
c3_shifted = c3.subs(n, n-3)
ratio = Poly(c0, n, domain='QQ') / Poly(c3_shifted, n, domain='QQ')
# Simplify
from sympy import simplify, cancel
ratio_simplified = cancel(c0 / c3.subs(n, n-3))
print(f"\nc_0(n) / c_3(n-3) = {ratio_simplified}")
print("If this is a constant, the operator is essentially self-adjoint.")
