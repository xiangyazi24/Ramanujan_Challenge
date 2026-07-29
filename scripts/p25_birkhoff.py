#!/usr/bin/env python3
"""P2.5: Compute Birkhoff structure — asymptotic ratios and error analysis.
Goal: understand the convergence mechanism rigorously."""
import mpmath as mp
mp.mp.dps = 120

# Full CMF matrix M(n) from the challenge
def M_matrix(n):
    """Return 3x3 matrix M(n) as mpmath matrix."""
    n = mp.mpf(n)
    # Row 1
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    # Row 2: factor (n+2)^2
    m21 = (n+2)**2*(-2*n-5)*(8*n**4+80*n**3+296*n**2+482*n+291)
    m22 = (n+2)**2*(24*n**4+264*n**3+1020*n**2+1599*n+813)
    m23 = (n+2)**2*(-32*n**2-204*n-310)
    # Row 3: factor (n+2)^2
    m31 = (n+2)**2*(-2*n-5)*(56*n**4+584*n**3+2240*n**2+3740*n+2295)
    m32 = (n+2)**2*(152*n**4+1924*n**3+8776*n**2+16609*n+10884)
    m33 = (n+2)**2*(-2*n-5)*(96*n**2+556*n+810)
    return mp.matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])

# Initial matrix A (2x3)
A = mp.matrix([[30921, -32972, 8240],
               [33750, -36000, 9000]])

G = mp.catalan  # Catalan's constant

# Compute matrix product M_N = M(0)*M(1)*...*M(N-1)
def compute_product(N):
    prod = mp.eye(3)
    for k in range(N):
        prod = M_matrix(k) * prod
    return prod

print("=== Convergence verification ===")
print(f"G = {mp.nstr(G, 50)}")

errors = []
for N in [5, 10, 20, 40, 60, 80]:
    MN = compute_product(N)
    AMN = A * MN  # 2x3 matrix: rows are (P_1,P_2,P_3) and (Q_1,Q_2,Q_3)
    ratios = []
    for j in range(3):
        P = AMN[0,j]
        Q = AMN[1,j]
        ratio = P/Q
        err = abs(ratio - G)
        ratios.append(ratio)
        if N <= 20:
            print(f"  N={N:3d}, j={j+1}: P/Q - G = {mp.nstr(err, 6)}")
    # Check error decay
    if N >= 10:
        err_j1 = abs(G * AMN[1,0] - AMN[0,0])
        errors.append((N, err_j1))

print("\n=== Error decay analysis ===")
rho_theory = 17 - 12*mp.sqrt(2)
c0 = mp.mpf(-16)
c_plus = -16*(17 + 12*mp.sqrt(2))

for i in range(1, len(errors)):
    N1, e1 = errors[i-1]
    N2, e2 = errors[i]
    if e1 > 0 and e2 > 0:
        ratio = mp.log(e2/e1) / (N2 - N1)
        eff_rho = mp.exp(ratio)
        print(f"  N={N1}→{N2}: |error ratio|^(1/(N2-N1)) = {mp.nstr(abs(eff_rho), 10)}")
        print(f"    theory: |c0/c+| = {mp.nstr(abs(c0/c_plus), 10)}")

# Extract the Birkhoff coefficients
print("\n=== Birkhoff analysis: R_0(X) components ===")
p0 = mp.matrix([[30921, -32972, 8240]])  # row vector
q0 = mp.matrix([[33750, -36000, 9000]])
# R_0(X) = q0 - (1+X)*p0 = (q0 - p0) - X*p0
A_coeff = q0 - p0  # constant part of R_0(X)
B_coeff = -p0  # coefficient of X

print(f"  A_0 (const part): {[int(A_coeff[0,j]) for j in range(3)]}")
print(f"  B_0 (X coeff):    {[int(B_coeff[0,j]) for j in range(3)]}")

# Compute R_N(X) = R_0 * M_N: A_N = A_0 * M_N, B_N = B_0 * M_N
for N in [1, 2, 5, 10]:
    MN = compute_product(N)
    AN = A_coeff * MN  # row vector
    BN = B_coeff * MN  # row vector
    for j in range(3):
        a, b = AN[0,j], BN[0,j]
        # Error = G*a + (1-G)*b = integral
        err_from_integral = G * (a - b) + b  # = G*Q - P since Q=A-B, P=-B
        err_direct = G * (a-b) - (-b)
        print(f"  N={N:2d}, j={j+1}: A={mp.nstr(a,6):>15s}, B={mp.nstr(b,6):>15s}, "
              f"R_N(0)={mp.nstr(a,6):>12s}, R_N(1)={mp.nstr(a+b,6):>12s}")
