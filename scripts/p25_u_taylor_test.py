#!/usr/bin/env python3
"""P2.5: Compute Taylor coefficients of U(z) = Y(k(z))/k(z) and test
whether they satisfy Q̂_n's order-3 recurrence.

U(z) = (π/4) Σ_{n≥0} C(2n,n)²·2^n/(2n+1) · z^n/(1-z)^{2n}

So [z^N] U(z) = (π/4) · r_N where:
r_N = δ_{N,0} + Σ_{n=1}^{N} C(2n,n)²·2^n/(2n+1) · C(N+n-1, N-n)

If r_N satisfies Q̂_n's recurrence, then U is a solution of L_rec,
proving L_{3,U} shares a common right factor with L_rec.
"""
from fractions import Fraction
from math import comb

# Q̂_n recurrence coefficients
c_coeffs = [
    [-170972650800, -826494925500, -1792449886332, -2317972607944, -2000297648936,
     -1219354055500, -541255279788, -177419351856, -43002662976, -7620091136,
     -960400960, -81589760, -4190208, -98304],
    [8781630505200, 38850314624124, 78557994908508, 96136040496551, 79442239242197,
     46814452218572, 20241514501104, 6502490145168, 1552168938336, 271943188864,
     33995217088, 2871763456, 146952192, 3440640],
    [-21132458248680, -87529225645944, -165451256319618, -189073879129764, -145809619841418,
     -80164318460172, -32338316008004, -9694892892592, -2160716677664, -353683596544,
     -41340724928, -3268370944, -156684288, -3440640],
    [587448626688, 2442715444224, 4635428285664, 5317694979920, 4116150568664,
     2270943978716, 919036676572, 276298241680, 61721801728, 10120470656,
     1184128064, 93632000, 4485120, 98304],
]

def eval_c(i, n):
    val = Fraction(0)
    nk = Fraction(1)
    for coeff in c_coeffs[i]:
        val += coeff * nk
        nk *= n
    return val

# Compute r_N = [z^N] U(z) / (π/4)
NMAX = 60
print("Computing U(z) Taylor coefficients r_N...", flush=True)

r = [Fraction(0)] * (NMAX + 10)
r[0] = Fraction(1)

for N in range(1, NMAX + 10):
    val = Fraction(0)
    for n in range(1, N + 1):
        c2n_n = comb(2*n, n)
        coeff = Fraction(c2n_n**2) * Fraction(2)**n / Fraction(2*n + 1)
        binom = comb(N + n - 1, N - n)
        val += coeff * Fraction(binom)
    r[N] = val
    if N <= 5 or N % 10 == 0:
        print(f"  r_{N} = {float(r[N]):.10e}", flush=True)

print(f"\nFirst few r_N values:")
for N in range(8):
    print(f"  r_{N} = {r[N]}")

# Test: does r_N satisfy Q̂_n's recurrence?
print(f"\n=== Testing r_N against Q̂_n's recurrence ===", flush=True)
all_ok = True
for n in range(NMAX):
    if n + 3 >= len(r):
        break
    res = sum(eval_c(j, Fraction(n)) * r[n + j] for j in range(4))
    if res != 0:
        if n < 5 or all_ok:
            print(f"  n={n}: residual = {res} (FAIL)")
            all_ok = False
    else:
        if n < 5:
            print(f"  n={n}: OK")

if all_ok:
    print(f"  ALL PASS for n = 0..{min(NMAX-1, len(r)-4)}")
    print("  => U(z) satisfies Q̂_n's recurrence!")
    print("  => L_{3,U} and L_rec share a common right factor")
else:
    print(f"  r_N does NOT satisfy Q̂_n's recurrence")

# Also check: does r_N satisfy the Delannoy-square recurrence?
# D_n^2 satisfies: (n+1)^2 a_{n+1} = 2(2n+1)(5*D_n - ... hmm
# Actually the Apery-like recurrence for D_n^2:
# (n+1)^2 D_{n+1}^2 = (something) — no, D_n^2 satisfies an order-3 recurrence
# Let me use the known ODE coefficients instead

# Delannoy-square recurrence: L_{3,D} in Euler form
# S_0(T) = 2T^3 - 3T^2
# S_1(T) = -70T^3 - 35T^2 + 18T + 9
# S_2(T) = 70T^3 + 175T^2 + 122T + 26
# S_3(T) = -2T^3 - 9T^2 - 12T - 5
def S0(T): return 2*T**3 - 3*T**2
def S1(T): return -70*T**3 - 35*T**2 + 18*T + 9
def S2(T): return 70*T**3 + 175*T**2 + 122*T + 26
def S3(T): return -2*T**3 - 9*T**2 - 12*T - 5

print(f"\n=== Testing r_N against Delannoy-square Euler recurrence ===", flush=True)
all_ok2 = True
for n in range(3, min(NMAX, len(r) - 1)):
    res = S0(Fraction(n)) * r[n] + S1(Fraction(n-1)) * r[n-1] + \
          S2(Fraction(n-2)) * r[n-2] + S3(Fraction(n-3)) * r[n-3]
    if res != 0:
        if n < 8 or all_ok2:
            print(f"  n={n}: residual = {float(res):.6e} (FAIL)")
            all_ok2 = False
    else:
        if n < 8:
            print(f"  n={n}: OK")

if all_ok2:
    print(f"  ALL PASS — r_N also satisfies the Delannoy Euler recurrence!")
else:
    print(f"  r_N does NOT satisfy the Delannoy Euler recurrence")

# Cross-check: compute Q̂_n and see if r_N is proportional
print(f"\n=== Comparing r_N with Q̂_n ===", flush=True)
def M_entries(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def H(nn):
    if nn == 0: return Fraction(1)
    val = Fraction(1)
    for k in range(nn):
        val *= Fraction(-16)*(k+2)**2*(k+3)**2*Fraction(2*k+5,2)*Fraction(2*k+7,2)**2
    return val

row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
qhat = [Fraction(33750)]
for n in range(20):
    M = M_entries(n)
    new_row = [Fraction(0)]*3
    for j in range(3):
        for k in range(3):
            new_row[j] += row[k]*Fraction(M[k][j])
    row = new_row
    qhat.append(Fraction(row[0], H(n+1)))

print("Q̂_n / r_N ratios (should be constant if proportional):")
for n in range(min(10, len(qhat))):
    if r[n] != 0:
        ratio = qhat[n] / r[n]
        print(f"  n={n}: Q̂_{n}/r_{n} = {ratio}")

print("\nDone.")
