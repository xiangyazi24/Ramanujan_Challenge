#!/usr/bin/env python3
"""P2.5: Test if Q̂_n ∈ span{D_n², D_n·E_n, E_n²} where D_n, E_n are
the two fundamental Delannoy solutions.

If this holds, the CMF is gauge-equivalent to Sym²(Delannoy), and the
P2.5 gap closes via the known Delannoy irrationality measure for G.
"""
from fractions import Fraction

# --- CMF data ---
def M_entries(n):
    n = Fraction(n)
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

def mat_mul(A, B):
    """3x3 matrix multiply"""
    n = len(A)
    m = len(B[0])
    k = len(B)
    C = [[Fraction(0)]*m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for l in range(k):
                C[i][j] += A[i][l] * B[l][j]
    return C

def mat_vec_mul(A, v):
    """Matrix × column vector"""
    return [sum(A[i][j]*v[j] for j in range(len(v))) for i in range(len(A))]

# A = initial matrix
A = [[Fraction(30921), Fraction(-32972), Fraction(8240)],
     [Fraction(33750), Fraction(-36000), Fraction(9000)]]

# H_n = (-16)^n * (2)_n^2 * (3)_n^2 * (5/2)_n * (7/2)_n^2
def H(n):
    val = Fraction(1)
    for k in range(n):
        val *= Fraction(-16)
        val *= Fraction(k+2)**2       # (2)_n
        val *= Fraction(k+3)**2       # (3)_n
        val *= Fraction(2*k+5, 2)     # (5/2)_n
        val *= Fraction(2*k+7, 2)**2  # (7/2)_n
    return val

# Compute Q_{N,0} = (A · M(0) · M(1) · ... · M(N-1))[1][0] and Q̂_N = Q_{N,0}/H_N
# Do this incrementally: row × M(n)
print("=== Computing Q̂_n from CMF ===", flush=True)
# Start with row 1 of A = [33750, -36000, 9000]
q_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]

Qhat = []
# N=0: Q_{0,0} = A[1][0] = 33750, H_0 = 1, Q̂_0 = 33750
Qhat.append(q_row[0] / H(0))
print(f"  Q̂[0] = {Qhat[0]}")

NMAX = 25
for N in range(NMAX):
    M = M_entries(N)
    # q_row = q_row · M(N)
    new_row = [Fraction(0)]*3
    for j in range(3):
        for k in range(3):
            new_row[j] += q_row[k] * M[k][j]
    q_row = new_row
    # Q_{N+1,0} = q_row[0], Q̂_{N+1} = Q_{N+1,0} / H_{N+1}
    Qhat.append(q_row[0] / H(N+1))

for n in range(min(6, len(Qhat))):
    print(f"  Q̂[{n}] = {float(Qhat[n]):.10e}")

# --- Delannoy solutions ---
# D_n = P_n(3) (Legendre at x=3), recurrence: (n+1)D_{n+1} = 3(2n+1)D_n - nD_{n-1}
# Solution 1: D_0=1, D_1=3 (central Delannoy numbers)
# Solution 2: E_0=0, E_1=1

def delannoy_pair(N):
    D = [Fraction(1), Fraction(3)]
    E = [Fraction(0), Fraction(1)]
    for n in range(1, N):
        coeff_main = Fraction(3*(2*n+1))
        coeff_prev = Fraction(n)
        denom = Fraction(n+1)
        D.append((coeff_main * D[n] - coeff_prev * D[n-1]) / denom)
        E.append((coeff_main * E[n] - coeff_prev * E[n-1]) / denom)
    return D, E

D, E = delannoy_pair(NMAX + 2)
print(f"\n=== Delannoy solutions ===")
for n in range(6):
    print(f"  D[{n}] = {D[n]}, E[{n}] = {float(E[n]):.6f}")

# Verify D_n are Delannoy: D_0=1, D_1=3, D_2=13, D_3=63, D_4=321
assert D[0] == 1 and D[1] == 3 and D[2] == 13 and D[3] == 63

# --- Sym² basis: g₁=D², g₂=D·E, g₃=E² ---
g1 = [D[n]**2 for n in range(len(D))]
g2 = [D[n]*E[n] for n in range(len(D))]
g3 = [E[n]**2 for n in range(len(E))]

# --- Test: Q̂_n = α·D_n² + β·D_n·E_n + γ·E_n² ---
print(f"\n=== Test: Q̂_n ∈ span{{D², D·E, E²}} ===")
# Use n=0,1,2 to solve for α,β,γ, verify at n=3,...
# System: [[g1[0], g2[0], g3[0]], [g1[1], g2[1], g3[1]], [g1[2], g2[2], g3[2]]] · [α,β,γ]ᵀ = [Q̂[0], Q̂[1], Q̂[2]]

M_sys = [[g1[0], g2[0], g3[0]],
         [g1[1], g2[1], g3[1]],
         [g1[2], g2[2], g3[2]]]
rhs = [Qhat[0], Qhat[1], Qhat[2]]

# Solve 3x3 with Fraction
def solve_3x3_exact(M, rhs):
    A = [[M[i][j] for j in range(3)] for i in range(3)]
    b = list(rhs)
    for col in range(3):
        pivot = None
        for row in range(col, 3):
            if A[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            return None
        A[col], A[pivot] = A[pivot], A[col]
        b[col], b[pivot] = b[pivot], b[col]
        for row in range(col+1, 3):
            factor = A[row][col] / A[col][col]
            for j in range(3):
                A[row][j] -= factor * A[col][j]
            b[row] -= factor * b[col]
    x = [Fraction(0)]*3
    for col in range(2, -1, -1):
        x[col] = b[col]
        for j in range(col+1, 3):
            x[col] -= A[col][j] * x[j]
        x[col] /= A[col][col]
    return x

sol = solve_3x3_exact(M_sys, rhs)
if sol is None:
    print("  System is singular!")
else:
    alpha, beta, gamma = sol
    print(f"  α = {float(alpha):.10e}")
    print(f"  β = {float(beta):.10e}")
    print(f"  γ = {float(gamma):.10e}")

    # Verify at n=3,...,NMAX
    all_match = True
    for n in range(3, len(Qhat)):
        pred = alpha * g1[n] + beta * g2[n] + gamma * g3[n]
        if pred == Qhat[n]:
            if n < 8:
                print(f"  n={n}: EXACT match ✓")
        else:
            all_match = False
            rel_err = abs(float((pred - Qhat[n]) / Qhat[n])) if Qhat[n] != 0 else float('inf')
            print(f"  n={n}: MISMATCH, rel error = {rel_err:.3e}")
            if n >= 5:
                break

    if all_match:
        print(f"\n  *** Q̂_n ∈ span{{D², D·E, E²}} CONFIRMED for n=0..{len(Qhat)-1} ***")
        print(f"  α = {alpha}")
        print(f"  β = {beta}")
        print(f"  γ = {gamma}")
        print(f"\n  This means the CMF is gauge-equivalent to Sym²(Delannoy)!")
        print(f"  The P2.5 gap closes via the known Delannoy irrationality measure for G.")
    else:
        print(f"\n  Q̂_n is NOT in span{{D², D·E, E²}}.")
        print(f"  The CMF is NOT gauge-equivalent to Sym²(Delannoy) in this simple sense.")

        # Try with shifted indices or other normalizations
        print(f"\n=== Alt test: Q̂_n = α(n)·D² + β(n)·D·E + γ(n)·E² with linear α,β,γ ===")
        # 6 unknowns: α = a₀+a₁n, β = b₀+b₁n, γ = c₀+c₁n
        # Use n=0,...,5 (6 equations)
        if len(Qhat) >= 6:
            A_sys = []
            b_sys = []
            for n in range(6):
                row = [g1[n], Fraction(n)*g1[n], g2[n], Fraction(n)*g2[n], g3[n], Fraction(n)*g3[n]]
                A_sys.append(row)
                b_sys.append(Qhat[n])

            # Solve 6x6
            import numpy as np
            A_np = np.array([[float(x) for x in row] for row in A_sys])
            b_np = np.array([float(x) for x in b_sys])
            try:
                sol6 = np.linalg.solve(A_np, b_np)
                print(f"  Linear coeffs: a₀={sol6[0]:.6e}, a₁={sol6[1]:.6e}")
                print(f"                 b₀={sol6[2]:.6e}, b₁={sol6[3]:.6e}")
                print(f"                 c₀={sol6[4]:.6e}, c₁={sol6[5]:.6e}")
                # Verify at n=6,...
                max_err = 0
                for n in range(6, min(len(Qhat), 15)):
                    pred = (sol6[0]+sol6[1]*n)*float(g1[n]) + (sol6[2]+sol6[3]*n)*float(g2[n]) + (sol6[4]+sol6[5]*n)*float(g3[n])
                    actual = float(Qhat[n])
                    if actual != 0:
                        rel = abs((pred-actual)/actual)
                        max_err = max(max_err, rel)
                        if n < 10:
                            print(f"  n={n}: rel error = {rel:.3e}")
                print(f"  Max holdout rel error: {max_err:.3e}")
            except:
                print("  6x6 solve failed")

print("\nDone.")
