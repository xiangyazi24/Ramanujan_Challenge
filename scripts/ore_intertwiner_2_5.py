#!/usr/bin/env python3
"""Try to find the Ore intertwiner connecting CMF denominators to D_n^2.

The CMF denominator q_n should be expressible as:
  q_n = u_0(n)*h_n*D_n^2 + u_1(n)*h_{n+1}*D_{n+1}^2 + u_2(n)*h_{n+2}*D_{n+2}^2

where h_n = (-16)^n * (n!)^4 * (3/2)_n^3 is the hypergeometric twist,
and u_0, u_1, u_2 are rational functions of n.

Step 1: Compute the CMF denominator sequence from the matrix product.
Step 2: Compute D_n^2 and the twist h_n.
Step 3: Try to find rational u_0, u_1, u_2 by solving a linear system.
"""
from fractions import Fraction
from functools import reduce

def M(n):
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

def matmul(A, B):
    rows, cols = len(A), len(B[0])
    inner = len(B)
    return [[sum(A[i][k]*B[k][j] for k in range(inner)) for j in range(cols)] for i in range(rows)]

A_init = [[30921, -32972, 8240],
          [33750, -36000, 9000]]

# Compute CMF products and extract Q_{N,1} (second row, first column of A*M_N)
N_max = 25
prod = [[1,0,0],[0,1,0],[0,0,1]]
q_cmf = []
for n in range(N_max):
    AM = matmul(A_init, prod)
    q_cmf.append(AM[1][0])  # Q_{N,1}
    prod = matmul(prod, M(n))
# One more
AM = matmul(A_init, prod)
q_cmf.append(AM[1][0])

print("CMF denominators Q_{N,1} for N=0..5:")
for i in range(6):
    print(f"  Q_{i} = {q_cmf[i]}")

# Central Delannoy numbers
def delannoy(n_max):
    D = [1, 3]
    for n in range(1, n_max):
        D.append((3*(2*n+1)*D[-1] - n*D[-2]) // (n+1))
    return D

D = delannoy(N_max + 5)
D2 = [d*d for d in D]

print("\nD_n^2 for n=0..5:")
for i in range(6):
    print(f"  D_{i}^2 = {D2[i]}")

# Compute ratios q_cmf[n] / D_n^2
print("\nRatios Q_N / D_N^2:")
for i in range(min(15, len(q_cmf))):
    if D2[i] != 0:
        r = Fraction(q_cmf[i], D2[i])
        print(f"  N={i}: {float(r):.6e}  ({r})")

# The twist h_n = (-16)^n * (n!)^4 * (3/2)_n^3
# (3/2)_n = (3/2)(5/2)...(2n+1)/2 = prod_{k=0}^{n-1} (2k+3)/2
def pochhammer_3_2(n):
    """Compute (3/2)_n as a Fraction."""
    result = Fraction(1)
    for k in range(n):
        result *= Fraction(2*k+3, 2)
    return result

def factorial(n):
    result = 1
    for k in range(1, n+1):
        result *= k
    return result

def h(n):
    return (-16)**n * factorial(n)**4 * pochhammer_3_2(n)**3

print("\nTwisted ratios Q_N / (h_n * D_n^2):")
for i in range(min(10, len(q_cmf))):
    if D2[i] != 0:
        hi = h(i)
        r = Fraction(q_cmf[i]) / (hi * D2[i])
        print(f"  N={i}: {float(r):.10e}")

# Try without the twist: look for pattern in Q_N / D_N^2
print("\n--- Looking for polynomial pattern in Q_N / D_N^2 ---")
ratios = []
for i in range(min(20, len(q_cmf))):
    if D2[i] != 0:
        ratios.append(Fraction(q_cmf[i], D2[i]))

# Check if ratios are polynomial in N
# Compute successive differences
if len(ratios) > 5:
    diffs = [ratios[i+1] - ratios[i] for i in range(len(ratios)-1)]
    print("First differences Q_{N+1}/D_{N+1}^2 - Q_N/D_N^2:")
    for i in range(min(8, len(diffs))):
        print(f"  delta_{i} = {float(diffs[i]):.6e}")

    diffs2 = [diffs[i+1] - diffs[i] for i in range(len(diffs)-1)]
    print("Second differences:")
    for i in range(min(8, len(diffs2))):
        print(f"  delta2_{i} = {float(diffs2[i]):.6e}")
