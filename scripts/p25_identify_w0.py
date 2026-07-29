#!/usr/bin/env python3
"""P2.5: Compute w_+(0) to high precision and attempt period identification.

If we can express w_+(0) in terms of known constants (G, sqrt(2), log 2, pi, ...),
the Catalan assertion (p0 - G*q0) . w_+(0) = 0 becomes provable.
"""
from mpmath import mp, mpf, matrix, catalan, sqrt, pi, log, identify, pslq

mp.dps = 200

G = catalan
sqrt2 = sqrt(2)
log2 = log(2)

def M_exact(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])

def delta(n):
    return mpf(-2) * (n+2)**2 * (n+3)**2 * (2*n+5) * (2*n+7)**2

def A_matrix(n):
    M = M_exact(n)
    d = delta(n)
    D_inv = matrix([
        [1, 0, 0],
        [0, mpf(1)/(n+1) if n >= 0 else 1, 0],
        [0, 0, mpf(1)/(n+1)**2 if n >= 0 else 1]
    ])
    D_next = matrix([
        [1, 0, 0],
        [0, mpf(n+2), 0],
        [0, 0, mpf(n+2)**2]
    ])
    return (mpf(1)/d) * D_inv * M * D_next

p0 = matrix([[30921, -32972, 8240]])
q0 = matrix([[33750, -36000, 9000]])

lam_plus = 17 + 12*sqrt2
v_plus = matrix([[2], [-sqrt2], [1]])

N_max = 200

# Backward iteration to compute w_+(0)
w = v_plus.copy()
for n in range(N_max - 1, -1, -1):
    An = A_matrix(n)
    w = An * w
    norm = max(abs(w[0,0]), abs(w[1,0]), abs(w[2,0]))
    if norm > 0:
        w = w / norm

w1, w2, w3 = w[0,0], w[1,0], w[2,0]

print(f"Precision: {mp.dps} digits, N_max = {N_max}")
print(f"w_+(0) = ({mp.nstr(w1, 60)},")
print(f"          {mp.nstr(w2, 60)},")
print(f"          {mp.nstr(w3, 60)})")
print()

# Verify L = G
num = p0[0,0]*w1 + p0[0,1]*w2 + p0[0,2]*w3
den = q0[0,0]*w1 + q0[0,1]*w2 + q0[0,2]*w3
L = num / den
res = (p0[0,0] - G*q0[0,0])*w1 + (p0[0,1] - G*q0[0,1])*w2 + (p0[0,2] - G*q0[0,2])*w3

print(f"L = {mp.nstr(L, 60)}")
print(f"G = {mp.nstr(G, 60)}")
print(f"|L - G| = {mp.nstr(abs(L-G), 15)}")
print(f"residual = {mp.nstr(res, 15)}")
print()

# Normalize w so that w3 = 1 (or w1 = 1)
if w1 != 0:
    u1 = w1/w1  # = 1
    u2 = w2/w1
    u3 = w3/w1
    print(f"Normalized w/w1: (1, {mp.nstr(u2, 50)}, {mp.nstr(u3, 50)})")
if w3 != 0:
    v1 = w1/w3
    v2 = w2/w3
    v3 = w3/w3  # = 1
    print(f"Normalized w/w3: ({mp.nstr(v1, 50)}, {mp.nstr(v2, 50)}, 1)")
print()

# Try PSLQ to identify w2/w1 in terms of known constants
print("="*60)
print("PSLQ identification attempts:")
print("="*60)

# w2/w1 in terms of {1, sqrt(2), G, pi, log(2), G*sqrt(2), pi*sqrt(2)}
r = u2  # = w2/w1
constants = [1, sqrt2, G, pi, log2, G*sqrt2, pi*sqrt2, G**2, pi**2]
names = ["1", "√2", "G", "π", "log2", "G√2", "π√2", "G²", "π²"]

# PSLQ: find integer relation among {r, 1, sqrt2, G, pi, log2, ...}
mp.dps = 180
r_high = w2/w1
r2_high = w3/w1

# Try r = w2/w1 against small basis
print(f"\nTrying r = w2/w1 = {mp.nstr(r_high, 30)}")
for basis, bnames in [
    ([r_high, mpf(1), sqrt2], ["r", "1", "√2"]),
    ([r_high, mpf(1), sqrt2, G], ["r", "1", "√2", "G"]),
    ([r_high, mpf(1), sqrt2, G, pi], ["r", "1", "√2", "G", "π"]),
    ([r_high, mpf(1), sqrt2, G, log2], ["r", "1", "√2", "G", "log2"]),
    ([r_high, mpf(1), sqrt2, G, G*sqrt2], ["r", "1", "√2", "G", "G√2"]),
    ([r_high, mpf(1), sqrt2, G, pi, log2], ["r", "1", "√2", "G", "π", "log2"]),
]:
    try:
        rel = pslq(basis)
        if rel is not None:
            terms = [f"{c}·{n}" for c, n in zip(rel, bnames) if c != 0]
            print(f"  FOUND: {' + '.join(terms)} = 0")
            if rel[0] != 0:
                val = -sum(c*v for c, v in zip(rel[1:], basis[1:])) / rel[0]
                print(f"  => r = {mp.nstr(val, 30)} (check: {mp.nstr(abs(r_high - val), 5)})")
            break
    except:
        pass

# Try r2 = w3/w1
print(f"\nTrying r2 = w3/w1 = {mp.nstr(r2_high, 30)}")
for basis, bnames in [
    ([r2_high, mpf(1), sqrt2], ["r2", "1", "√2"]),
    ([r2_high, mpf(1), sqrt2, G], ["r2", "1", "√2", "G"]),
    ([r2_high, mpf(1), sqrt2, G, pi], ["r2", "1", "√2", "G", "π"]),
    ([r2_high, mpf(1), sqrt2, G, log2], ["r2", "1", "√2", "G", "log2"]),
    ([r2_high, mpf(1), sqrt2, G, G*sqrt2], ["r2", "1", "√2", "G", "G√2"]),
    ([r2_high, mpf(1), sqrt2, G, pi, log2], ["r2", "1", "√2", "G", "π", "log2"]),
]:
    try:
        rel = pslq(basis)
        if rel is not None:
            terms = [f"{c}·{n}" for c, n in zip(rel, bnames) if c != 0]
            print(f"  FOUND: {' + '.join(terms)} = 0")
            if rel[0] != 0:
                val = -sum(c*v for c, v in zip(rel[1:], basis[1:])) / rel[0]
                print(f"  => r2 = {mp.nstr(val, 30)} (check: {mp.nstr(abs(r2_high - val), 5)})")
            break
    except:
        pass

# Also try to identify the Catalan connection directly
# (p0 - G*q0) . w = 0 can be rewritten as:
# a1*w1 + a2*w2 + a3*w3 = 0 where a_i = p0_i - G*q0_i
a1 = 30921 - G*33750
a2 = -32972 - G*(-36000)
a3 = 8240 - G*9000
print(f"\na1 = p0_1 - G*q0_1 = {mp.nstr(a1, 30)}")
print(f"a2 = p0_2 - G*q0_2 = {mp.nstr(a2, 30)}")
print(f"a3 = p0_3 - G*q0_3 = {mp.nstr(a3, 30)}")
print(f"a1*w1 + a2*w2 + a3*w3 = {mp.nstr(a1*w1 + a2*w2 + a3*w3, 15)}")

# The ratio a2/a1 and a3/a1 might be recognizable
if a1 != 0:
    print(f"\na2/a1 = {mp.nstr(a2/a1, 30)}")
    print(f"a3/a1 = {mp.nstr(a3/a1, 30)}")

# Try to identify w_+(0) components individually
print(f"\n{'='*60}")
print("Individual component identification:")
print(f"{'='*60}")
for i, wi in enumerate([w1, w2, w3]):
    result = identify(wi, tol=1e-25)
    if result:
        print(f"  w{i+1} ≈ {result}")
    else:
        print(f"  w{i+1}: no identification found")

# Also try ratios of components
print()
for label, val in [("w2/w1", w2/w1), ("w3/w1", w3/w1), ("w3/w2", w3/w2),
                    ("w1/w2", w1/w2), ("w1/w3", w1/w3)]:
    result = identify(val, tol=1e-25)
    if result:
        print(f"  {label} ≈ {result}")
