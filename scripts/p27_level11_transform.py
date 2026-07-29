#!/usr/bin/env python3
"""P2.7: Verify Q4878's level-11 → P2.7 transform.

Cooper's level-11 sequence T_k satisfies:
  (k+1)³T_{k+1} = 2(2k+1)(5k²+5k+2)T_k - 8k(7k²+1)T_{k-1} + 22k(2k-1)(k-1)T_{k-2}

Poincaré polynomial: H₁₁(t) = t³ - 20t² + 56t - 44

Transform: W_n = (1/256^n) · Σ_{j=0}^{2n} C(2n,j)(-2)^{2n-j} T_j

Claim: W_n (or an Ore-transformed version) satisfies the P2.7 recurrence.
"""
from fractions import Fraction
import math

# P2.7 recurrence coefficients
def A27(n):
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

def B27(n):
    n_ = n
    return 128*(2*n_+7)**3*(2*n_+9)**3*(
        104060*n_**6 + 1077948*n_**5 + 4656738*n_**4 +
        10724703*n_**3 + 13897425*n_**2 + 9627369*n_ + 2790720)

def C27(n):
    n_ = n
    return 16*(n_+3)**4*(2*n_+9)**3*(
        3784*n_**5 + 34600*n_**4 + 125824*n_**3 +
        227612*n_**2 + 205120*n_ + 73664)

def D27(n):
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Initial values for P2.7
q = [Fraction(-215040420000),
     Fraction(-167282265043404, 905),
     Fraction(-964185327658080, 6071)]

# Extend q using recurrence: A(n)q_{n+3} + B(n)q_{n+2} + C(n)q_{n+1} + D(n)q_n = 0
# q_{n+3} = -(B(n)q_{n+2} + C(n)q_{n+1} + D(n)q_n) / A(n)
NMAX = 100
for n in range(NMAX):
    qnew = -(Fraction(B27(n))*q[n+2] + Fraction(C27(n))*q[n+1] + Fraction(D27(n))*q[n]) / Fraction(A27(n))
    q.append(qnew)

print(f"=== P2.7 q_n: computed {len(q)} values ===")
for n in range(4):
    print(f"  q[{n}] = {float(q[n]):.10e}")

# --- Cooper's level-11 sequence ---
T = [Fraction(1), Fraction(4), Fraction(28)]
for k in range(2, 300):
    # (k+1)³T_{k+1} = 2(2k+1)(5k²+5k+2)T_k - 8k(7k²+1)T_{k-1} + 22k(2k-1)(k-1)T_{k-2}
    num = (2*(2*k+1)*(5*k**2+5*k+2)*T[k]
           - 8*k*(7*k**2+1)*T[k-1]
           + 22*k*(2*k-1)*(k-1)*T[k-2])
    T.append(num / Fraction((k+1)**3))

print(f"\n=== Cooper level-11 T_k: computed {len(T)} values ===")
for k in range(6):
    print(f"  T[{k}] = {T[k]}")

# Verify T_3 = 268
assert T[3] == 268, f"T[3] = {T[3]}, expected 268"

# Verify Poincaré: T_{k+1}/T_k → dominant root of H₁₁(t) = t³-20t²+56t-44
print(f"\n  T[50]/T[49] = {float(T[50]/T[49]):.10f}")
print(f"  T[100]/T[99] = {float(T[100]/T[99]):.10f}")
print(f"  T[200]/T[199] = {float(T[200]/T[199]):.10f}")

# --- Binomial transform ---
# b_m = Σ_{j=0}^m C(m,j)(-2)^{m-j} T_j
def compute_bm(m):
    val = Fraction(0)
    for j in range(m+1):
        binom_mj = math.comb(m, j)
        val += binom_mj * Fraction(-2)**(m-j) * T[j]
    return val

# W_n = b_{2n} / 256^n
print(f"\n=== Transform W_n = b_{{2n}} / 256^n ===")
W = []
for n in range(min(NMAX + 4, 103)):
    b2n = compute_bm(2*n)
    Wn = b2n / Fraction(256)**n
    W.append(Wn)
    if n < 6 or n % 20 == 0:
        print(f"  W[{n}] = {float(Wn):.10e}")

# --- Test: does W_n satisfy P2.7 recurrence directly? ---
print(f"\n=== Test: W_n vs P2.7 recurrence ===")
for n in range(min(len(W)-3, 20)):
    res = A27(n)*W[n+3] + B27(n)*W[n+2] + C27(n)*W[n+1] + D27(n)*W[n]
    if res == 0:
        print(f"  n={n}: residual = 0 (EXACT)")
    else:
        print(f"  n={n}: residual = {float(res):.3e} (nonzero)")
    if n >= 4 and res != 0:
        print(f"  W_n does NOT satisfy P2.7 recurrence directly.")
        break

# --- Compare ratios q_n/W_n to look for pattern ---
print(f"\n=== Ratios q_n/W_n ===")
for n in range(min(len(W), 10)):
    if W[n] != 0:
        r = q[n] / W[n]
        print(f"  q[{n}]/W[{n}] = {float(r):.10e}")

# --- Search for order-2 Ore relation: q_n = h_n(a(n)W_n + b(n)W_{n+1} + c(n)W_{n+2}) ---
# If h_n = product of Pochhammer factors, start with h_n = 1 and search for polynomial a,b,c
print(f"\n=== Search: q_n = a(n)W_n + b(n)W_{n+1} + c(n)W_{n+2} (constant coefficients) ===")
# 3 unknowns a, b, c. Use n=0,1,2 to solve, verify at n=3.
from fractions import Fraction

# q[0] = a*W[0] + b*W[1] + c*W[2]
# q[1] = a*W[1] + b*W[2] + c*W[3]
# q[2] = a*W[2] + b*W[3] + c*W[4]
# Verify: q[3] = a*W[3] + b*W[4] + c*W[5]

M = [[W[0], W[1], W[2]],
     [W[1], W[2], W[3]],
     [W[2], W[3], W[4]]]
rhs = [q[0], q[1], q[2]]

# Gaussian elimination (3x3)
def solve_3x3(M, rhs):
    A = [list(row) for row in M]
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
            for j in range(col, 3):
                A[row][j] -= factor * A[col][j]
            b[row] -= factor * b[col]
    # Back substitution
    x = [Fraction(0)]*3
    for col in range(2, -1, -1):
        x[col] = b[col]
        for j in range(col+1, 3):
            x[col] -= A[col][j] * x[j]
        x[col] /= A[col][col]
    return x

sol = solve_3x3(M, rhs)
if sol:
    a, b_coeff, c = sol
    print(f"  a = {float(a):.6e}")
    print(f"  b = {float(b_coeff):.6e}")
    print(f"  c = {float(c):.6e}")
    # Verify
    for n in range(3, min(len(W)-2, 10)):
        pred = a*W[n] + b_coeff*W[n+1] + c*W[n+2]
        if q[n] == pred:
            print(f"  n={n}: EXACT match")
        else:
            print(f"  n={n}: mismatch, rel error = {float(abs((q[n]-pred)/q[n])):.3e}")

# --- Search for q_n = f(n) * W_n (proportional up to rational function) ---
print(f"\n=== Ratios q_n/W_n (looking for polynomial pattern) ===")
ratios = []
for n in range(min(len(W), 20)):
    if W[n] != 0:
        r = q[n] / W[n]
        ratios.append((n, r))
        if n < 6:
            print(f"  q[{n}]/W[{n}] = {r}")

# Check if ratios form a polynomial pattern
if len(ratios) >= 5:
    print(f"\n  Ratio differences (looking for polynomial):")
    for i in range(1, min(5, len(ratios))):
        print(f"  r[{ratios[i][0]}]/r[{ratios[i-1][0]}] = {float(ratios[i][1]/ratios[i-1][1]):.10f}")

print("\nDone.")
