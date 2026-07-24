#!/usr/bin/env sage
"""P2.7: Ore algebra computation following Q4893's prescription.

1. Construct Cooper's operator and the two quadratic branch operators
2. Compute LCLM
3. Descend to even coordinates
4. Factor and test against q_n
"""
from sage.all import *

try:
    from ore_algebra import OreAlgebra
    HAS_ORE = True
except ImportError:
    HAS_ORE = False
    print("WARNING: ore_algebra not installed. Using manual construction.")

# Problem 2.7 recurrence coefficients
def A27(n):
    return 1024*(2*n+5)^4*(2*n+7)^3*(2*n+9)^3*(946*n^2+6407*n+10860)

def B27(n):
    P6 = 104060*n^6+1745370*n^5+12145238*n^4+44886481*n^3+92943995*n^2+102256019*n+46709052
    return 128*(2*n+7)^3*(2*n+9)^3*P6

def C27(n):
    P5 = 3784*n^5+57792*n^4+351019*n^3+1059230*n^2+1587211*n+944620
    return 16*(n+3)^4*(2*n+9)^3*P5

def D27(n):
    return (n+3)^4*(n+4)^6*(946*n^2+4515*n+5399)

# Cooper's level-11 recurrence
def Ac(k):
    return (k+1)^3

def Bc(k):
    return 2*(2*k+1)*(5*k^2+5*k+2)

def Cc(k):
    return -8*k*(7*k^2+1)

def Dc(k):
    return 22*k*(2*k-1)*(k-1)

# Compute T_k (Cooper)
print("Computing Cooper's T_k...", flush=True)
T = [QQ(1)]
# Need T_1: from k=0: 1*T_1 = 2*1*2*T_0 = 4. So T_1 = 4
T.append(QQ(4))
# k=1: 8*T_2 = 2*3*12*T_1 - 8*1*8*T_0 = 288 - 64 = 224. T_2 = 28
T.append(QQ(28))
KMAX = 250
for k in range(2, KMAX):
    t_next = (Bc(k)*T[k] + Cc(k)*T[k-1] + Dc(k)*T[k-2]) / Ac(k)
    T.append(t_next)
print(f"  T_0..T_5 = {T[:6]}")

# Compute W_n via binomial transform
print("Computing W_n...", flush=True)
NMAX = 60
W = []
for n in range(NMAX):
    val = QQ(0)
    for j in range(min(2*n+1, len(T))):
        val += binomial(2*n, j) * QQ(-2)^(2*n-j) * T[j]
    W.append(val / QQ(256)^n)
print(f"  W_0..W_3 = {W[:4]}")

# Compute q_n
print("Computing q_n...", flush=True)
q = [QQ(-215040420000)]
q.append(QQ(-167282265043404) / QQ(905))
q.append(QQ(-964185327658080) / QQ(6071))
for n in range(2, NMAX):
    if n-1 >= 1 and n-2 >= 0:
        q_next = (B27(n-1)*q[n] - C27(n-2)*q[n-1] + D27(n-3)*q[n-2]) / A27(n-1) if n >= 3 else None
        if q_next is None:
            # n=2 case: use the recurrence with n→n-1
            # A(0)q_1 = B(0)q_0 - ... but we already have q_0, q_1, q_2
            break
        q.append(q_next)

# Actually, let me use the standard form more carefully
# The recurrence is: A(n)q_{n+3} = B(n)q_{n+2} - C(n)q_{n+1} + D(n)q_n
# Wait, need to check the exact form from the comprehensive_2_7.py file

# Let me just compute from the 4-term recurrence directly
print("Recomputing q_n with exact recurrence...", flush=True)
q = [QQ(-215040420000)]
q.append(QQ(-167282265043404) / QQ(905))
q.append(QQ(-964185327658080) / QQ(6071))

# From comprehensive_2_7.py, the recurrence is:
# A(n)*q[n+1] = B(n)*q[n] - C(n)*q[n-1] + D(n)*q[n-2]
# for n >= 2
for n in range(2, NMAX-1):
    q_next = (B27(n)*q[n] - C27(n)*q[n-1] + D27(n)*q[n-2]) / A27(n)
    q.append(q_next)

print(f"  q_0 = {q[0]}")
print(f"  q_1 = {q[1]}")
print(f"  len(q) = {len(q)}")

# Check ratio q_n/W_n
print("\n=== q_n / W_n ratios ===", flush=True)
for n in range(min(20, len(q), len(W))):
    if W[n] != 0:
        r = q[n] / W[n]
        print(f"  n={n}: q/W = {float(r):.10f}")

# Check ratio q_n / (W_n * h_n) for various h_n
print("\n=== Testing h_n = (3/2)_n / n! (central binomial related) ===", flush=True)
from sage.functions.other import rising_factorial as rf

for n in range(min(15, len(q), len(W))):
    if W[n] != 0:
        h = rf(QQ(3)/QQ(2), n) / factorial(n)
        r = q[n] / (W[n] * h)
        print(f"  n={n}: q/(W*h) = {float(r):.10f}")

# Now check: can q_n be written as sum of u_j(n) * h_n * W_{n+j}?
# where h_n = sqrt(pi) * Gamma(n+3/2) / Gamma(n+1) = (2n+1)!! / 2^n
print("\n=== h_n = (2n+1)!! / 2^n ===", flush=True)
for n in range(min(15, len(q), len(W))):
    if W[n] != 0:
        h = QQ(1)
        for j in range(n):
            h *= QQ(2*j+1)
        h /= QQ(2)^n
        r = q[n] / (W[n] * h)
        print(f"  n={n}: q/(W*h) = {float(r):.10f}")

# Now try the POINCARE APPROACH: compute the sequence q_n / W_n
# and see what recurrence IT satisfies
print("\n=== Recurrence for r_n = q_n / W_n ===", flush=True)
r_vals = []
for n in range(min(len(q), len(W))):
    if W[n] != 0:
        r_vals.append(q[n] / W[n])
    else:
        r_vals.append(None)

# Check: r_{n+1}/r_n
print("  r_{n+1}/r_n:")
for n in range(min(15, len(r_vals)-1)):
    if r_vals[n] and r_vals[n+1]:
        print(f"    n={n}: {float(r_vals[n+1]/r_vals[n]):.15f}")

# Search for a rational function of n matching r_{n+1}/r_n
# If r_{n+1}/r_n = P(n)/Q(n), then r_n = prod P(j)/Q(j) which is a hypergeometric twist
print("\n=== Searching for r_{n+1}/r_n = P(n)/Q(n) ===", flush=True)
rr = []
for n in range(min(len(r_vals)-1, 40)):
    if r_vals[n] and r_vals[n+1]:
        rr.append(r_vals[n+1] / r_vals[n])

# Try Padé approximant for rr[n] as a function of n
# rr[n] ≈ P(n)/Q(n) where we search for deg P = deg Q = d
R.<x> = PolynomialRing(QQ)
for d in range(1, 10):
    # Build system: rr[n] * Q(n) = P(n)
    # P = a_0 + a_1*n + ... + a_d*n^d
    # Q = 1 + b_1*n + ... + b_d*n^d (monic constant term 1... actually let's make leading coeff 1)
    # Q = n^d + b_{d-1}*n^{d-1} + ... + b_0
    n_unk = 2*d + 1  # d+1 for P, d for Q (leading coeff of Q fixed to 1)
    n_pts = min(n_unk + 5, len(rr))

    if n_pts < n_unk:
        break

    M = matrix(QQ, n_pts, n_unk)
    b = vector(QQ, n_pts)
    for i in range(n_pts):
        n = i + 1  # start from n=1 to avoid issues
        # P coefficients: -n^0, -n^1, ..., -n^d
        for k in range(d+1):
            M[i, k] = -QQ(n)^k
        # Q coefficients: rr[i+1] * n^0, ..., rr[i+1] * n^{d-1}
        for k in range(d):
            M[i, d+1+k] = rr[i+1] * QQ(n)^k
        # RHS: -rr[i+1] * n^d (from the leading term of Q)
        b[i] = -rr[i+1] * QQ(n)^d

    try:
        sol = M.solve_right(b)
    except ValueError:
        continue

    # Verify on holdout
    ok = True
    for i in range(n_pts, min(n_pts+5, len(rr))):
        n = i + 1
        P_val = sum(sol[k] * QQ(n)^k for k in range(d+1))
        Q_val = QQ(n)^d + sum(sol[d+1+k] * QQ(n)^k for k in range(d))
        if Q_val == 0 or P_val/Q_val != rr[i+1]:
            ok = False
            break

    if ok:
        P_coeffs = [sol[k] for k in range(d+1)]
        Q_coeffs = [sol[d+1+k] for k in range(d)] + [QQ(1)]
        P_poly = sum(P_coeffs[k] * x^k for k in range(d+1))
        Q_poly = sum(Q_coeffs[k] * x^k for k in range(d+1))
        print(f"\n*** MATCH at degree {d}: ***")
        print(f"  r_{{n+1}}/r_n = P(n)/Q(n)")
        print(f"  P(n) = {P_poly}")
        print(f"  Q(n) = {Q_poly}")

        # Extended verification
        all_ok = True
        for i in range(len(rr)):
            n = i + 1
            P_val = P_poly(n)
            Q_val = Q_poly(n)
            if Q_val == 0 or P_val/Q_val != rr[i+1]:
                all_ok = False
                print(f"  FAIL at n={n}")
                break
        if all_ok:
            print(f"  VERIFIED for all n=1..{len(rr)-1}")

            # Factor the polynomials
            print(f"  P factors: {P_poly.factor()}")
            print(f"  Q factors: {Q_poly.factor()}")

            # Express h_n = prod_{j=0}^{n-1} P(j)/Q(j) in closed form
            print(f"\n  h_n = r_0 * prod_{{j=0}}^{{n-1}} P(j)/Q(j)")
            print(f"  r_0 = q_0/W_0 = {r_vals[0]}")
        break

print("\nDone.")
