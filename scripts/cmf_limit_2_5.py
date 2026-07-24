#!/usr/bin/env python3
"""
Problem 2.5: Compute the CMF limit correctly.

From the challenge PDF:
  M_N = M(0)·M(1)·...·M(N-1)
  A = ((30921, -32972, 8240), (33750, -36000, 9000))
  A·M_N = ((P_{N,1}, P_{N,2}, P_{N,3}), (Q_{N,1}, Q_{N,2}, Q_{N,3}))
  Prove: lim P_{N,j}/Q_{N,j} = G for j=1,2,3.

Also: verify the Brafman integral and identify the proof architecture.
"""
from mpmath import mp, mpf, matrix, catalan, log, sqrt, pi, ellipk, quad, power, fac

mp.dps = 80

def M_entry(n):
    """Return the 3x3 CMF matrix M(n) as mpmath matrix."""
    n = mpf(n)
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)
    return matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])

A = matrix([[30921, -32972, 8240], [33750, -36000, 9000]])

# Compute M_N = M(0)·M(1)·...·M(N-1)
N_max = 80
MN = matrix([[1,0,0],[0,1,0],[0,0,1]])

G_exact = catalan
print(f"Catalan G = {G_exact}")
print(f"ln(2)/2  = {log(2)/2}")
print()

for N in range(N_max):
    MN = MN * M_entry(N)
    if N in [10, 20, 30, 40, 50, 60, 70, 79]:
        AMN = A * MN
        print(f"N = {N+1}:")
        for j in range(3):
            P = AMN[0, j]
            Q = AMN[1, j]
            ratio = P / Q
            err = ratio - G_exact
            digits = -mp.log10(abs(err)) if err != 0 else mp.inf
            print(f"  P/Q col {j+1}: {mp.nstr(ratio, 30)}  ({mp.nstr(digits, 4)} digits)")
        print()

# Now verify the Brafman integral
print("=" * 60)
print("BRAFMAN INTEGRAL VERIFICATION")
print("=" * 60)

rho = 17 - 12*sqrt(2)
print(f"rho = 17-12*sqrt(2) = {rho}")

# Brafman's identity: F(z) = sum D_n^2 z^n = (2/(pi*(1-z))) * K(4*sqrt(2z)/(1-z))
# where K(k) = elliptic integral with modulus k

# Central Delannoy numbers
def delannoy(n_max):
    D = [mpf(1), mpf(3)]
    for n in range(1, n_max):
        D.append((3*(2*n+1)*D[n] - n*D[n-1]) / (n+1))
    return D

D = delannoy(200)

# Verify Brafman at z = 0.01
for z_val in [mpf('0.01'), mpf('0.1'), rho/2, rho*mpf('0.99')]:
    F_series = sum(D[n]**2 * z_val**n for n in range(200))
    k_mod = 4*sqrt(2*z_val) / (1 - z_val)
    if abs(k_mod) < 1:
        F_brafman = 2 / (pi * (1 - z_val)) * ellipk(k_mod**2)
        err = abs(F_series - F_brafman)
        print(f"z = {mp.nstr(z_val, 6)}: diff = {mp.nstr(err, 5)}")
    else:
        print(f"z = {mp.nstr(z_val, 6)}: k = {mp.nstr(k_mod, 6)} >= 1, boundary")

# The Catalan integral: G = (pi*sqrt(2)/2) * integral_0^rho [(1+z)/(sqrt(z)*(1-z))] * F(z) dz
# where F(z) = sum D_n^2 z^n

# Using Brafman substitution:
# G = (pi*sqrt(2)/2) * integral_0^rho [(1+z)/(sqrt(z)*(1-z))] * [2/(pi*(1-z))] * K(4*sqrt(2z)/(1-z)) dz
#   = sqrt(2) * integral_0^rho [(1+z)/(sqrt(z)*(1-z)^2)] * K(4*sqrt(2z)/(1-z)) dz

print()
print("Catalan integral verification:")
print(f"  G (exact) = {G_exact}")

# Direct numerical integration
def integrand(z):
    if z <= 0:
        return mpf(0)
    k_sq = 32*z / (1-z)**2
    if k_sq >= 1:
        return mpf(0)
    return (1+z) / (sqrt(z) * (1-z)**2) * ellipk(k_sq)

G_int_val = sqrt(2) * quad(integrand, [mpf('1e-15'), rho * mpf('0.9999')])
print(f"  G (integral) = {G_int_val}")
print(f"  diff         = {mp.nstr(abs(G_int_val - G_exact), 5)}")

# Alternative: use the series form
# G = (pi*sqrt(2)/2) * sum_{n>=0} D_n^2 * integral_0^rho z^{n-1/2} * (1+z)/(1-z) dz
# This separates the D_n^2 series from the integral.

# Let I_n = integral_0^rho z^{n-1/2} * (1+z)/(1-z) dz
# = integral_0^rho z^{n-1/2}/(1-z) dz + integral_0^rho z^{n+1/2}/(1-z) dz
# The function (1+z)/(1-z) = 1 + 2z/(1-z)

# Actually let me compute term by term
print()
print("Term-by-term Catalan sum:")
G_sum = mpf(0)
for n in range(150):
    # I_n = integral_0^rho z^n * (1+z)/(sqrt(z)*(1-z)) dz
    def term_integrand(z, n=n):
        return z**n * (1+z) / (sqrt(z) * (1-z))
    I_n = quad(term_integrand, [mpf('1e-15'), rho])
    contrib = D[n]**2 * I_n
    G_sum += contrib
    if n % 30 == 0 or n < 5:
        partial = pi * sqrt(2) / 2 * G_sum
        err_digits = -mp.log10(abs(partial - G_exact)) if partial != G_exact else 999
        print(f"  n={n}: partial = {mp.nstr(partial, 30)}, digits = {mp.nstr(err_digits, 4)}")

G_brafman = pi * sqrt(2) / 2 * G_sum
print(f"\n  G (Brafman sum) = {G_brafman}")
print(f"  G (exact)       = {G_exact}")
print(f"  diff            = {mp.nstr(abs(G_brafman - G_exact), 5)}")
