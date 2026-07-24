#!/usr/bin/env python3
"""Deep PSLQ search for w_+(0) components of P2.5.

Known: w_2 = -1 exactly. Try to identify w_1 and w_3.
"""
from mpmath import mp, mpf, matrix, catalan, sqrt, pi, log, pslq, euler, zeta, rf

mp.dps = 200

G = catalan
sqrt2 = sqrt(2)
log2 = log(2)
pi_val = pi

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

v_plus = matrix([[2], [-sqrt2], [1]])
N_max = 250

w = v_plus.copy()
for n in range(N_max - 1, -1, -1):
    An = A_matrix(n)
    w = An * w
    norm = max(abs(w[0,0]), abs(w[1,0]), abs(w[2,0]))
    if norm > 0:
        w = w / norm

w1, w2, w3 = w[0,0], w[1,0], w[2,0]
print(f"w_+(0) = ({mp.nstr(w1, 60)}, {mp.nstr(w2, 10)}, {mp.nstr(w3, 60)})")

# More constants for PSLQ
Q0_3 = log((mpf(3)+1)/(mpf(3)-1))/2  # Legendre Q_0(3) = log(2)/2
P1_3 = mpf(3)  # P_1(3) = 3
Q1_3 = mpf(3)*Q0_3 - 1  # Q_1(3) = 3*log(2)/2 - 1

K_half = mp.ellipk(mpf(1)/2)  # K(1/sqrt(2))
E_half = mp.ellipe(mpf(1)/2)
gamma4 = mp.gamma(mpf(1)/4)

# Additional: Legendre-related constants
# Q_0(3) = log(2)/2 = 0.34657...

print(f"\nConstants:")
print(f"  G = {mp.nstr(G, 30)}")
print(f"  √2 = {mp.nstr(sqrt2, 30)}")
print(f"  log2 = {mp.nstr(log2, 30)}")
print(f"  Q₀(3) = log2/2 = {mp.nstr(Q0_3, 30)}")
print(f"  Q₁(3) = 3log2/2-1 = {mp.nstr(Q1_3, 30)}")
print(f"  K(1/2) = {mp.nstr(K_half, 30)}")
print(f"  E(1/2) = {mp.nstr(E_half, 30)}")
print(f"  Γ(1/4) = {mp.nstr(gamma4, 30)}")
print()

mp.dps = 180

# System approach: instead of identifying w1, w3 separately,
# try to find the relation that makes L = G.
#
# The Catalan assertion is: (p0 - G*q0) . w = 0
# With w2 = -1: a1*w1 - a2 + a3*w3 = 0
# where a1 = 30921 - 33750G, a2 = 32972 - 36000G (note sign), a3 = 8240 - 9000G
#
# So: w3 = (a2 - a1*w1) / a3
# We just need to identify w1 (or w3).

# Try w1 against extended basis
bases = [
    # Rational + √2
    ([w1, mpf(1), sqrt2], "1, √2"),
    # Rational + √2 + log2
    ([w1, mpf(1), sqrt2, log2], "1, √2, log2"),
    # Rational + √2 + G
    ([w1, mpf(1), sqrt2, G], "1, √2, G"),
    # Rational + √2 + G + log2
    ([w1, mpf(1), sqrt2, G, log2], "1, √2, G, log2"),
    # Rational + √2 + G√2
    ([w1, mpf(1), sqrt2, G, G*sqrt2], "1, √2, G, G√2"),
    # Q_0(3) = log2/2
    ([w1, mpf(1), sqrt2, Q0_3], "1, √2, Q₀(3)"),
    # Include Q₁(3)
    ([w1, mpf(1), sqrt2, Q0_3, Q1_3], "1, √2, Q₀(3), Q₁(3)"),
    # K and E
    ([w1, mpf(1), sqrt2, K_half, E_half], "1, √2, K(1/2), E(1/2)"),
    # G and K
    ([w1, mpf(1), sqrt2, G, K_half], "1, √2, G, K(1/2)"),
    # All available
    ([w1, mpf(1), sqrt2, G, log2, K_half], "1, √2, G, log2, K(1/2)"),
    # π and √2
    ([w1, mpf(1), sqrt2, pi_val], "1, √2, π"),
    # π² and G
    ([w1, mpf(1), sqrt2, G, pi_val**2], "1, √2, G, π²"),
    # Gamma(1/4) related
    ([w1, mpf(1), sqrt2, gamma4**2/pi_val], "1, √2, Γ(1/4)²/π"),
    # Direct elliptic AGM constant
    ([w1, mpf(1), sqrt2, G, pi_val, log2, K_half, E_half], "1, √2, G, π, log2, K, E"),
]

print("PSLQ for w1:")
print("="*60)
for basis, name in bases:
    try:
        rel = pslq(basis, maxcoeff=10000)
        if rel is not None and rel[0] != 0:
            terms = []
            names = ["w1"] + name.split(", ")
            for c, n in zip(rel, names):
                if c != 0:
                    terms.append(f"({c}){n}")
            print(f"  [{name}]: {' + '.join(terms)} = 0")
            val = -sum(c*v for c, v in zip(rel[1:], basis[1:])) / rel[0]
            print(f"    => w1 = {mp.nstr(val, 30)}, check: |diff| = {mp.nstr(abs(w1-val), 5)}")
    except Exception as e:
        pass

# Same for w3
print(f"\nPSLQ for w3:")
print("="*60)
for basis, name in bases:
    basis2 = [w3] + basis[1:]
    try:
        rel = pslq(basis2, maxcoeff=10000)
        if rel is not None and rel[0] != 0:
            terms = []
            names = ["w3"] + name.split(", ")
            for c, n in zip(rel, names):
                if c != 0:
                    terms.append(f"({c}){n}")
            print(f"  [{name}]: {' + '.join(terms)} = 0")
            val = -sum(c*v for c, v in zip(rel[1:], basis2[1:])) / rel[0]
            print(f"    => w3 = {mp.nstr(val, 30)}, check: |diff| = {mp.nstr(abs(w3-val), 5)}")
    except Exception as e:
        pass

# Also try the RATIO w3/w1 (might be simpler)
print(f"\nPSLQ for w3/w1:")
print("="*60)
r31 = w3/w1
for basis, name in bases:
    basis3 = [r31] + basis[1:]
    try:
        rel = pslq(basis3, maxcoeff=10000)
        if rel is not None and rel[0] != 0:
            terms = []
            names = ["w3/w1"] + name.split(", ")
            for c, n in zip(rel, names):
                if c != 0:
                    terms.append(f"({c}){n}")
            print(f"  [{name}]: {' + '.join(terms)} = 0")
    except Exception as e:
        pass

# Try the Catalan constraint directly
# a1*w1 + a3*w3 = a2 (with known a1, a2, a3)
# This means w1 and w3 satisfy a linear relation with G-dependent coefficients.
# We need ONE more relation to determine both.
#
# The OTHER relation comes from the RECESSIVE mode: w_-(0).
# If we compute w_-(0) similarly, we get another functional.
# But the recessive mode is harder to compute (it's the smallest eigenvalue).

# Instead, let's try: does the SECOND Birkhoff functional (neutral or recessive)
# give another identifying equation?

# Actually, for the proof we don't need to identify w1 and w3 separately.
# We just need to prove (p0 - G*q0) . w_+(0) = 0.
# This is equivalent to proving that the initial conditions lie in the
# neutral+recessive subspace of the adjoint.

print(f"\n{'='*60}")
print("Looking for the proof structure:")
print(f"{'='*60}")
print(f"The equation (p0 - G*q0) . w_+(0) = 0 means:")
print(f"  {mp.nstr(mpf(30921)-G*33750, 20)} * w1 + {mp.nstr(mpf(-32972)+G*36000, 20)} * (-1) + {mp.nstr(mpf(8240)-G*9000, 20)} * w3 = 0")
print()

# Check: is there a FINITE closed-form for the inner product?
# The quantity v = (p0 - G*q0) is a vector with G-dependent entries.
# w_+(0) is determined by the recurrence.
# Their inner product being zero is the DEFINITION of G being the limit.
# Can we express this as a PERIOD integral?

# The Brafman identity + integrated K module suggests:
# The fundamental matrix at z=0 of the integrated K ODE gives the recurrence.
# The monodromy at z=rho (k=1) gives the connection matrix.
# The connection constant for the integrated solution IS G.

# So the proof strategy is:
# 1. Show the CMF recurrence matches the integrated K ODE under Brafman substitution
# 2. The monodromy theory then gives L = G automatically

# Let me verify step 1 computationally:
# The ODE is k(1-k²)Y''' + (1-3k²)Y'' - kY' = 0
# Under k = 4√(2z)/(1-z), this becomes an ODE in z.
# The power series at z=0 gives a recurrence.
# Compare with the CMF scalar recurrence.

print("This is the key: verify the CMF recurrence matches the integrated K ODE")
print("under the Brafman substitution k = 4√(2z)/(1-z).")
