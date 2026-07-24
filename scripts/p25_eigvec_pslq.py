#!/usr/bin/env python3
"""
P2.5: High-precision PSLQ on Birkhoff eigenvector components.

Key insight: k(ρ) = 1, so the singularity of the generating function
at z = ρ is controlled by K(1) = ∞ (log singularity).

Goal: identify r21 = lim u2/u1 and r31 = lim u3/u1 in terms of
known constants (G, π, √2, log(1+√2), K-values, etc.)
"""
from mpmath import mp, mpf, matrix, catalan, pi, log, sqrt, euler
from mpmath import ellipk, ellipe, pslq, identify, atanh, asinh
from mpmath import quad, inf

mp.dps = 500

def M_entries(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

G = catalan
rho = 17 - 12*sqrt(2)
xi_plus = 17 + 12*sqrt(2)

print(f"Working at {mp.dps}-digit precision")
print(f"G = {mp.nstr(G, 50)}")
print(f"ρ = {mp.nstr(rho, 50)}")

NMAX = 350  # enough for 500-digit precision

# Compute three fundamental solutions
e1 = [mpf(1), mpf(0), mpf(0)]
e2 = [mpf(0), mpf(1), mpf(0)]
e3 = [mpf(0), mpf(0), mpf(1)]

u1_vals = [e1[0]]
u2_vals = [e2[0]]
u3_vals = [e3[0]]

print(f"Computing {NMAX} CMF steps...", flush=True)
for N in range(NMAX):
    M = M_entries(mpf(N))
    d = mpf(delta_H(N))
    MH = [[M[i][j]/d for j in range(3)] for i in range(3)]

    e1_new = [sum(e1[i]*MH[i][j] for i in range(3)) for j in range(3)]
    e2_new = [sum(e2[i]*MH[i][j] for i in range(3)) for j in range(3)]
    e3_new = [sum(e3[i]*MH[i][j] for i in range(3)) for j in range(3)]

    e1, e2, e3 = e1_new, e2_new, e3_new
    u1_vals.append(e1[0])
    u2_vals.append(e2[0])
    u3_vals.append(e3[0])

    if N % 50 == 49:
        r21 = u2_vals[-1] / u1_vals[-1]
        r31 = u3_vals[-1] / u1_vals[-1]
        print(f"  N={N+1}: r21 = {mp.nstr(r21, 30)}", flush=True)

r21 = u2_vals[NMAX] / u1_vals[NMAX]
r31 = u3_vals[NMAX] / u1_vals[NMAX]

print(f"\n=== Eigenvector components ===")
print(f"r21 = {mp.nstr(r21, 80)}")
print(f"r31 = {mp.nstr(r31, 80)}")

# Verify L = G
q_dot_v = 33750 + (-36000)*r21 + 9000*r31
p_dot_v = 30921 + (-32972)*r21 + 8240*r31
L = p_dot_v / q_dot_v
err = L - G
print(f"\nL = {mp.nstr(L, 50)}")
print(f"G = {mp.nstr(G, 50)}")
print(f"|L-G| = {mp.nstr(abs(err), 10)}")

# Known constants for PSLQ basis
s2 = sqrt(2)
log2 = log(2)
log_1p_s2 = log(1 + s2)  # = arcsinh(1)
pi_val = pi
G_val = G
K_half = ellipk(mpf(1)/2)  # K(1/2)
K_s2_2 = ellipk(s2/2)  # K(√2/2)
E_half = ellipe(mpf(1)/2)

# Also compute integral G = (1/2)∫₀¹ K(k) dk
# And related integrals

print(f"\n=== PSLQ analysis ===")
print(f"Known constants:")
print(f"  G = {mp.nstr(G_val, 30)}")
print(f"  π = {mp.nstr(pi_val, 30)}")
print(f"  √2 = {mp.nstr(s2, 30)}")
print(f"  log2 = {mp.nstr(log2, 30)}")
print(f"  log(1+√2) = {mp.nstr(log_1p_s2, 30)}")
print(f"  K(1/2) = {mp.nstr(K_half, 30)}")
print(f"  K(√2/2) = {mp.nstr(K_s2_2, 30)}")
print(f"  E(1/2) = {mp.nstr(E_half, 30)}")

# Test various PSLQ bases
def try_pslq(name, target, basis, labels):
    vec = [target] + list(basis)
    result = pslq(vec, maxcoeff=10**8, maxsteps=5000)
    if result:
        terms = []
        for i, (c, l) in enumerate(zip(result[1:], labels)):
            if c != 0:
                terms.append(f"({c})·{l}")
        print(f"  {name}: {result[0]}·target + {' + '.join(terms)} = 0")
        if result[0] != 0:
            expr = " + ".join(f"({-c}/{result[0]})·{l}" for c, l in zip(result[1:], labels) if c != 0)
            print(f"    → target = {expr}")
    else:
        print(f"  {name}: no relation found")

# r21 tests
print(f"\n--- r21 = {mp.nstr(r21, 20)} ---")

try_pslq("r21 vs {G,1}", r21, [G_val, mpf(1)], ["G", "1"])
try_pslq("r21 vs {G,√2,1}", r21, [G_val, s2, mpf(1)], ["G", "√2", "1"])
try_pslq("r21 vs {G,π,1}", r21, [G_val, pi_val, mpf(1)], ["G", "π", "1"])
try_pslq("r21 vs {G,π,√2,1}", r21, [G_val, pi_val, s2, mpf(1)], ["G", "π", "√2", "1"])
try_pslq("r21 vs {G,G√2,√2,1}", r21, [G_val, G_val*s2, s2, mpf(1)], ["G", "G√2", "√2", "1"])
try_pslq("r21 vs {G,π,log2,√2,1}", r21, [G_val, pi_val, log2, s2, mpf(1)], ["G", "π", "log2", "√2", "1"])
try_pslq("r21 vs {G,K(1/2),K(√2/2),√2,1}", r21, [G_val, K_half, K_s2_2, s2, mpf(1)], ["G", "K½", "K(√2/2)", "√2", "1"])
try_pslq("r21 vs {G,log(1+√2),π,√2,1}", r21, [G_val, log_1p_s2, pi_val, s2, mpf(1)], ["G", "log(1+√2)", "π", "√2", "1"])

# r21 quadratic
try_pslq("r21² vs {r21,G,√2,1}", r21**2, [r21, G_val, s2, mpf(1)], ["r21", "G", "√2", "1"])

# r31 tests
print(f"\n--- r31 = {mp.nstr(r31, 20)} ---")
try_pslq("r31 vs {G,1}", r31, [G_val, mpf(1)], ["G", "1"])
try_pslq("r31 vs {G,√2,1}", r31, [G_val, s2, mpf(1)], ["G", "√2", "1"])
try_pslq("r31 vs {G,π,1}", r31, [G_val, pi_val, mpf(1)], ["G", "π", "1"])
try_pslq("r31 vs {G,π,√2,1}", r31, [G_val, pi_val, s2, mpf(1)], ["G", "π", "√2", "1"])
try_pslq("r31 vs {G,G√2,√2,1}", r31, [G_val, G_val*s2, s2, mpf(1)], ["G", "G√2", "√2", "1"])
try_pslq("r31 vs {G,K(1/2),K(√2/2),√2,1}", r31, [G_val, K_half, K_s2_2, s2, mpf(1)], ["G", "K½", "K(√2/2)", "√2", "1"])

# Cross-relations r21, r31
print(f"\n--- Cross-relations ---")
try_pslq("r21 vs {r31,G,√2,1}", r21, [r31, G_val, s2, mpf(1)], ["r31", "G", "√2", "1"])
try_pslq("r21·r31 vs {G,√2,1}", r21*r31, [G_val, s2, mpf(1)], ["G", "√2", "1"])
try_pslq("r21+r31 vs {G,√2,1}", r21+r31, [G_val, s2, mpf(1)], ["G", "√2", "1"])

# Direct eigenvector components
print(f"\n--- Projections ---")
print(f"q·v₊ = {mp.nstr(q_dot_v, 30)}")
print(f"p·v₊ = {mp.nstr(p_dot_v, 30)}")

try_pslq("q·v₊ vs {G,π,√2,1}", q_dot_v, [G_val, pi_val, s2, mpf(1)], ["G", "π", "√2", "1"])
try_pslq("p·v₊ vs {G,π,√2,1}", p_dot_v, [G_val, pi_val, s2, mpf(1)], ["G", "π", "√2", "1"])
try_pslq("q·v₊ vs {G²,G,√2,1}", q_dot_v, [G_val**2, G_val, s2, mpf(1)], ["G²", "G", "√2", "1"])
try_pslq("q·v₊ vs {K(1/2),G,π,√2,1}", q_dot_v, [K_half, G_val, pi_val, s2, mpf(1)], ["K½", "G", "π", "√2", "1"])

# Check individual u-ratios at specific N
print(f"\n--- Convergence check ---")
for N in [100, 200, 300, NMAX]:
    r21_N = u2_vals[N] / u1_vals[N]
    r31_N = u3_vals[N] / u1_vals[N]
    err21 = abs(r21_N - r21)
    err31 = abs(r31_N - r31)
    d21 = -mp.log10(err21) if err21 > 0 else mp.dps
    d31 = -mp.log10(err31) if err31 > 0 else mp.dps
    print(f"  N={N}: r21 stable to {mp.nstr(d21,5)} digits, r31 to {mp.nstr(d31,5)} digits")

print("\nDone.")
