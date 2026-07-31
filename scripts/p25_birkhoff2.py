#!/usr/bin/env python3
"""P2.5: Birkhoff analysis with CORRECT matrix product order.
M_N = M(0) * M(1) * ... * M(N-1) [left to right]."""
import mpmath as mp
mp.mp.dps = 120

def M_matrix(n):
    n = mp.mpf(n)
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(-2*n-5)*(8*n**4+80*n**3+296*n**2+482*n+291)
    m22 = (n+2)**2*(24*n**4+264*n**3+1020*n**2+1599*n+813)
    m23 = (n+2)**2*(-32*n**2-204*n-310)
    m31 = (n+2)**2*(-2*n-5)*(56*n**4+584*n**3+2240*n**2+3740*n+2295)
    m32 = (n+2)**2*(152*n**4+1924*n**3+8776*n**2+16609*n+10884)
    m33 = (n+2)**2*(-2*n-5)*(96*n**2+556*n+810)
    return mp.matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])

A = mp.matrix([[30921, -32972, 8240],
               [33750, -36000, 9000]])
G = mp.catalan

def compute_product(N):
    """M_N = M(0) * M(1) * ... * M(N-1)"""
    prod = mp.eye(3)
    for k in range(N):
        prod = prod * M_matrix(k)  # CORRECT order: left to right
    return prod

print("=== Convergence verification (corrected order) ===")
print(f"G = {mp.nstr(G, 50)}")

prev_err = None
for N in [1, 2, 3, 5, 10, 15, 20, 30, 40, 50]:
    MN = compute_product(N)
    AMN = A * MN
    P1, Q1 = AMN[0,0], AMN[1,0]
    err = abs(P1/Q1 - G)
    if err > 0:
        digits = -float(mp.log10(err))
    else:
        digits = 120
    
    rate = ""
    if prev_err is not None and err > 0 and prev_err[1] > 0:
        dN = N - prev_err[0]
        log_ratio = mp.log(err / prev_err[1]) / dN
        eff_rho = mp.exp(log_ratio)
        rate = f"  rho_eff = {mp.nstr(abs(eff_rho), 6)}"
    
    print(f"  N={N:3d}: |P/Q - G| = {mp.nstr(err, 8):>15s}  ({digits:.1f} digits){rate}")
    prev_err = (N, err)

# Also check the three columns match
print("\n=== Column consistency at N=30 ===")
MN = compute_product(30)
AMN = A * MN
for j in range(3):
    P, Q = AMN[0,j], AMN[1,j]
    print(f"  j={j+1}: P/Q = {mp.nstr(P/Q, 50)}")

print(f"\n  G   = {mp.nstr(G, 50)}")
