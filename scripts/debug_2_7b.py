#!/usr/bin/env python3
"""Debug 2.7: Check sequence ratio behavior and coefficient values."""
from mpmath import mp, mpf, zeta, log10, fabs

mp.dps = 150

def A(n):
    n = mpf(n)
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

def B(n):
    n = mpf(n)
    P6 = 104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052
    return 128*(2*n+7)**3*(2*n+9)**3*P6

def C(n):
    n = mpf(n)
    P5 = 3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620
    return 16*(n+3)**4*(2*n+9)**3*P5

def D(n):
    n = mpf(n)
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Version 1 (PDF literal): u_{n+1} = B_n/A_n u_n - C_{n-1}/A_{n-1} u_{n-1} + D_{n-2}/A_{n-2} u_{n-2}
q = [mpf(-215040420000), mpf(-1672822650043404)/mpf(905), mpf(-964185327658080)/mpf(6071)]

for n in range(2, 100):
    q_next = B(n)/A(n)*q[n] - C(n-1)/A(n-1)*q[n-1] + D(n-2)/A(n-2)*q[n-2]
    q.append(q_next)

print("Ratio q_{n+1}/q_n (should approach Poincaré root 0.8588):")
for n in [5, 10, 20, 30, 40, 50, 60, 70, 80]:
    r = q[n]/q[n-1]
    print(f"  n={n}: q[n]/q[n-1] = {mp.nstr(r, 15)}")

print(f"\nPoincaré dominant root: 0.858807735912...")
print(f"q[5] = {mp.nstr(q[5], 15)}")
print(f"q[10] = {mp.nstr(q[10], 15)}")
print(f"q[50] = {mp.nstr(q[50], 15)}")

# Check: coefficient ratios at large n
print("\nCoefficient ratios B(n)/A(n) at various n:")
for n in [5, 10, 50]:
    print(f"  n={n}: B/A = {mp.nstr(B(n)/A(n), 10)}, C/A = {mp.nstr(C(n)/A(n), 10)}, D/A = {mp.nstr(D(n)/A(n), 10)}")

# At large n: B/A → B_lead/A_lead, C/A → C_lead/A_lead, D/A → D_lead/A_lead
# B_lead/A_lead = 852459520/991952896 = 0.85949...
# Wait, this doesn't match Poincaré root 0.8588!
print(f"\nLeading ratios:")
print(f"  B_lead/A_lead = {128*104060 / (1024*946)}")
# Need to include FULL leading terms:
# A ~ 1024 * 2^4 * 2^3 * 2^3 * 946 * n^{12}  → leading coeff in n^12
# NO: A(n) has terms (2n+5)^4(2n+7)^3(2n+9)^3(946n^2+...)
# Leading in n: 2^4*n^4 * 2^3*n^3 * 2^3*n^3 * 946*n^2 = 2^{10} * 946 * n^{12}
# = 1024 * 946 * n^{12} (factor of 1024 from overall constant already included)
# So: 1024 * 1024 * 946 = but wait A already has the 1024 in front
# A ~ 1024 * (2n)^4 * (2n)^3 * (2n)^3 * 946n^2 = 1024 * 16n^4 * 8n^3 * 8n^3 * 946n^2
# = 1024 * 1024 * 946 * n^{12}
A_lead = 1024 * 16 * 8 * 8 * 946
B_lead = 128 * 8 * 8 * 104060
C_lead = 16 * 1 * 8 * 3784  # (n+3)^4~n^4=1*n^4, (2n+9)^3~8n^3
D_lead = 1 * 1 * 946  # (n+3)^4~n^4, (n+4)^6~n^6, 946n^2

print(f"  A_lead = {A_lead}")
print(f"  B_lead = {B_lead}")
print(f"  C_lead = {C_lead}")
print(f"  D_lead = {D_lead}")
print(f"  B/A = {B_lead/A_lead:.10f}")
print(f"  C/A = {C_lead/A_lead:.10e}")
print(f"  D/A = {D_lead/A_lead:.10e}")

# Poincaré polynomial: c^3 - (B/A)c^2 + ... hmm
# Actually the Poincaré polynomial comes from the leading behavior of the recurrence
# For u_{n+1} = (B/A)u_n - (C_{n-1}/A_{n-1})u_{n-1} + (D_{n-2}/A_{n-2})u_{n-2}
# We need the Poincaré polynomial of the STANDARD form:
# A_n u_{n+1} - B_n u_n + C_{n-1} u_{n-1} - D_{n-2} u_{n-2} = 0

# But the shifted indices mean the Poincaré polynomial is:
# A∞ c^3 - B∞ c^2 + C∞ c - D∞ = 0  (with the ∞ meaning limiting leading coefficients)
# BUT C_{n-1} has DIFFERENT leading coefficient than C_n shifted by 1.
# Actually (n-1+3)^4 = (n+2)^4 ~ n^4, same leading. And (2(n-1)+9)^3 = (2n+7)^3 ~ 8n^3.
# And 3784(n-1)^5 ~ 3784 n^5. So C_{n-1} ~ C_lead * n^{12}. Same.
# Similarly D_{n-2} ~ D_lead * n^{12}.

# So the Poincaré polynomial is:
# A_lead c^3 - B_lead c^2 + C_lead c - D_lead = 0
import numpy as np
roots = np.roots([A_lead, -B_lead, C_lead, -D_lead])
print(f"  Poincaré roots: {roots}")

# Now verify the B/A ratio against Poincaré root
print(f"\n  B_lead/A_lead = {B_lead/A_lead:.10f}")
print(f"  Sum of roots = {sum(roots):.10f}")
print(f"  Should equal B_lead/A_lead: {B_lead/A_lead:.10f}")

# Let me also print B(n)/A(n) for very large n
for n in [100, 200, 500]:
    print(f"  B({n})/A({n}) = {mp.nstr(B(n)/A(n), 12)}")
