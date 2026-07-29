#!/usr/bin/env python3
"""Verify P2.5 Catalan assertion via adjoint dominant functional.

Computes w_+(0) by backward iteration of A_n w_+(n+1) = w_+(n),
then checks L = p_0 . w_+(0) / q_0 . w_+(0) = G.

This is the approach from Q4841: the Catalan assertion reduces to
the single scalar identity (p_0 - G q_0) . w_+(0) = 0.
"""
import mpmath as mp

mp.mp.dps = 150

def M_exact(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return mp.matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])

def delta(n):
    return mp.mpf(-2) * (n+2)**2 * (n+3)**2 * (2*n+5) * (2*n+7)**2

def A_matrix(n):
    """Normalized one-step: A_n = (H_n/H_{n+1}) D_n^{-1} M(n) D_{n+1}."""
    M = M_exact(n)
    d = delta(n)
    D_inv = mp.matrix([
        [1, 0, 0],
        [0, mp.mpf(1)/(n+1) if n >= 0 else 1, 0],
        [0, 0, mp.mpf(1)/(n+1)**2 if n >= 0 else 1]
    ])
    D_next = mp.matrix([
        [1, 0, 0],
        [0, mp.mpf(n+2), 0],
        [0, 0, mp.mpf(n+2)**2]
    ])
    return (mp.mpf(1)/d) * D_inv * M * D_next

p0 = mp.matrix([[30921, -32972, 8240]])
q0 = mp.matrix([[33750, -36000, 9000]])
G = mp.catalan

sqrt2 = mp.sqrt(2)
lambda_plus = 17 + 12*sqrt2
v_plus = mp.matrix([[2], [-sqrt2], [1]])

N_max = 120

# Forward: compute w_+(0) via projective power iteration
# [w_+(0)] = lim_{N->inf} [A_0 A_1 ... A_{N-1} v_+]
w = v_plus.copy()
for n in range(N_max - 1, -1, -1):
    An = A_matrix(n)
    w = An * w
    # Normalize to prevent overflow/underflow
    norm = max(abs(w[0,0]), abs(w[1,0]), abs(w[2,0]))
    if norm > 0:
        w = w / norm

w1, w2, w3 = w[0,0], w[1,0], w[2,0]

# Compute the limit
num = p0[0,0]*w1 + p0[0,1]*w2 + p0[0,2]*w3
den = q0[0,0]*w1 + q0[0,1]*w2 + q0[0,2]*w3
L = num / den

# Catalan assertion check
residual = (p0[0,0] - G*q0[0,0])*w1 + (p0[0,1] - G*q0[0,1])*w2 + (p0[0,2] - G*q0[0,2])*w3

print(f"Adjoint dominant functional w_+(0) (projective):")
print(f"  w1 = {mp.nstr(w1, 40)}")
print(f"  w2 = {mp.nstr(w2, 40)}")
print(f"  w3 = {mp.nstr(w3, 40)}")
print()
print(f"L = p0.w / q0.w = {mp.nstr(L, 50)}")
print(f"G =               {mp.nstr(G, 50)}")
print(f"|L - G| = {mp.nstr(abs(L - G), 15)}")
print()
print(f"Catalan residual (p0 - G*q0).w = {mp.nstr(residual, 15)}")
print(f"Relative residual / den = {mp.nstr(abs(residual/den), 15)}")
print()

# Also verify by direct matrix product
print("--- Direct matrix product verification ---")
prod = mp.eye(3)
for n in range(80):
    prod = prod * M_exact(n)
P_row = p0 * prod
Q_row = q0 * prod
for j in range(3):
    ratio = P_row[0,j] / Q_row[0,j]
    print(f"P/Q column {j}: {mp.nstr(ratio, 50)}")
print(f"G =              {mp.nstr(G, 50)}")
