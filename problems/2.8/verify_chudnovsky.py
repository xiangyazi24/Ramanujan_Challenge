#!/usr/bin/env python3
"""Problem 2.8: Verify the 4×4 CMF encodes Chudnovsky series for √10005/π.

Strategy: compute the matrix product A·M(0)·M(1)·...·M(N-1) and compare
each column's P/Q ratio with √10005/π.

The matrix M(n) entries are transcribed from the challenge PDF (pages 5-6)."""
from mpmath import mp, mpf, sqrt, pi, matrix

mp.dps = 100

R = mpf(151931373056001)

def make_M(n):
    """Construct the 4×4 matrix M(n) from the challenge PDF."""
    u = mpf(2*n + 3)
    w = u * (3*u - 2) * (3*u + 2)

    a1 = (144*R - 99)*u**5 - (288*R - 333)*u**4 + (144*R - 229)*u**3 - 114*u**2 + 40*u + 64
    a2 = (432*R - 243)*u**4 - (864*R - 909)*u**3 + (432*R - 868)*u**2 - 80*u + 272
    a3 = (432*R - 153)*u**3 - (864*R - 648)*u**2 + (432*R - 860)*u + 360

    b1 = 9*u**4 - (144*R - 63)*u**3 + 158*u**2 + 168*u + 64
    b2 = 36*u**3 + (216*R - 189)*u**2 - 316*u - 168
    b3 = 54*u**2 + (108*R - 189)*u - 158

    c1 = 18*u**5 + (54*R + 45)*u**4 - (288*R**2 - 378*R + 251)*u**3 \
         + (948*R - 1086)*u**2 + (1008*R - 1384)*u + (384*R - 576)
    c2 = (153*R - 72)*u**4 - (657*R - 702)*u**3 - (432*R**2 - 1292*R + 1069)*u**2 \
         + (2064*R - 2508)*u + (1072*R - 1512)
    c3 = (180*R - 108)*u**3 - (891*R - 864)*u**2 - (216*R**2 - 1450*R + 1385)*u \
         + (1116*R - 1422)
    c4 = (6*R - 4)*u**2 - (33*R - 32)*u - (4*R**2 - 58*R - 236337691420383)

    M = matrix(4, 4)
    M[0,0] = a1 / w;      M[0,1] = a2 / w;       M[0,2] = a3 / w;       M[0,3] = 144*R*(u-1)**2 / w
    M[1,0] = -u**3;        M[1,1] = -3*u**2;       M[1,2] = -3*u;         M[1,3] = -1
    M[2,0] = b1/(144*R);   M[2,1] = -b2/(72*R);    M[2,2] = b3/(36*R);    M[2,3] = -(2*u + 2*R - 7)/(2*R)
    M[3,0] = c1/(288*R**2); M[3,1] = c2/(144*R**2); M[3,2] = c3/(72*R**2); M[3,3] = c4/(4*R**2)

    return M

# Initial matrix A (2×4)
# Transcribed from PDF page 6 — very large integers (CORRECTED)
A1 = 37169305760442252761441
A2 = 111507917281327441564208
A3 = 111507917281327599720129
A4 = 371693057604442410917362

B1 = 1167416361542639692320
B2 = 3502249084627896132160
B3 = 3502249084627879697280
B4 = 1167416361542622723840

# Note: the paper has these as very large integers.
# Let me verify with actual values from the PDF...
# (A1,A2,A3,A4) and (B1,B2,B3,B4) as given

A_mat = matrix(2, 4)
A_mat[0,0] = mpf(A1); A_mat[0,1] = mpf(A2); A_mat[0,2] = mpf(A3); A_mat[0,3] = mpf(A4)
A_mat[1,0] = mpf(B1); A_mat[1,1] = mpf(B2); A_mat[1,2] = mpf(B3); A_mat[1,3] = mpf(B4)

target = sqrt(mpf(10005)) / pi
print(f"Target √10005/π = {target}")
print(f"R = {R}")
print(f"R - 1 = 640320^3/1728 = {mpf(640320)**3/1728}")
print(f"Actual R-1 = {R - 1}")
print(f"Match: {abs((R-1) - mpf(640320)**3/1728) < 1}")

# Compute matrix product
N = 20
prod = matrix(4, 4)
for i in range(4):
    prod[i,i] = 1  # identity

for n in range(N):
    Mn = make_M(n)
    prod = prod * Mn

result = A_mat * prod

print(f"\nA·M(0)·...·M({N-1}):")
for j in range(4):
    P = result[0,j]
    Q = result[1,j]
    if abs(Q) > 1e-10:
        ratio = P / Q
        print(f"  j={j+1}: P/Q = {ratio}")
        print(f"         diff = {ratio - target}")
