#!/usr/bin/env python3
"""Check if (S-1) is a factor of the order-3 recurrence for q_N.

The recurrence is: q_{N+3} = c₂(N)q_{N+2} + c₁(N)q_{N+1} + c₀(N)q_N
where c_i are rational functions of N.

Constants-killing: 1 - c₂(N) - c₁(N) - c₀(N) = 0 for all N.
"""
from mpmath import mp, mpf, nstr, matrix

mp.dps = 200

def M_mat(n):
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

A_mat = matrix([[mpf(30921), mpf(-32972), mpf(8240)],
                [mpf(33750), mpf(-36000), mpf(9000)]])

N_max = 40
s_vals = []
T = matrix([[1,0,0],[0,1,0],[0,0,1]])
for N in range(N_max + 1):
    AT = A_mat * T
    s_vals.append((AT[1,0], AT[1,1], AT[1,2]))
    T = T * M_mat(N)

print("=== Check 1 - c₂ - c₁ - c₀ = 0 for each N ===")
for N in range(3, 30):
    M_N = M_mat(N)
    M_N1 = M_mat(N+1)
    M_N2 = M_mat(N+2)

    m00 = M_N[0,0]; m01 = M_N[0,1]; m02 = M_N[0,2]
    m10 = M_N[1,0]; m11 = M_N[1,1]; m12 = M_N[1,2]
    m20 = M_N[2,0]; m21 = M_N[2,1]; m22 = M_N[2,2]
    m00p = M_N1[0,0]; m10p = M_N1[1,0]; m20p = M_N1[2,0]
    m01p = M_N1[0,1]; m11p = M_N1[1,1]; m21p = M_N1[2,1]
    m02p = M_N1[0,2]; m12p = M_N1[1,2]; m22p = M_N1[2,2]

    A11 = m10; A12 = m20
    A21 = m11*m10p + m12*m20p
    A22 = m21*m10p + m22*m20p
    det_A = A11*A22 - A12*A21

    coef_B1_q0 = -m00; coef_B1_q1 = mpf(1); coef_B1_q2 = mpf(0)
    coef_B2_q0 = -(m01*m10p + m02*m20p); coef_B2_q1 = -m00p; coef_B2_q2 = mpf(1)

    pN_q0 = (coef_B1_q0*A22 - coef_B2_q0*A12) / det_A
    pN_q1 = (coef_B1_q1*A22 - coef_B2_q1*A12) / det_A
    pN_q2 = (coef_B1_q2*A22 - coef_B2_q2*A12) / det_A
    rN_q0 = (A11*coef_B2_q0 - A21*coef_B1_q0) / det_A
    rN_q1 = (A11*coef_B2_q1 - A21*coef_B1_q1) / det_A
    rN_q2 = (A11*coef_B1_q2*0 + A11*coef_B2_q2 - A21*coef_B1_q2) / det_A

    # Fix rN_q2: = (A11*1 - A21*0) / det_A = A11/det_A
    rN_q2 = A11 / det_A

    pN1_q0 = m01 + pN_q0*m11 + rN_q0*m21
    pN1_q1 = pN_q1*m11 + rN_q1*m21
    pN1_q2 = pN_q2*m11 + rN_q2*m21
    rN1_q0 = m02 + pN_q0*m12 + rN_q0*m22
    rN1_q1 = pN_q1*m12 + rN_q1*m22
    rN1_q2 = pN_q2*m12 + rN_q2*m22

    pN2_q0 = pN1_q0*m11p + rN1_q0*m21p
    pN2_q1 = m01p + pN1_q1*m11p + rN1_q1*m21p
    pN2_q2 = pN1_q2*m11p + rN1_q2*m21p
    rN2_q0 = pN1_q0*m12p + rN1_q0*m22p
    rN2_q1 = m02p + pN1_q1*m12p + rN1_q1*m22p
    rN2_q2 = pN1_q2*m12p + rN1_q2*m22p

    m00pp = M_N2[0,0]; m10pp = M_N2[1,0]; m20pp = M_N2[2,0]

    c0 = m10pp*pN2_q0 + m20pp*rN2_q0
    c1 = m10pp*pN2_q1 + m20pp*rN2_q1
    c2 = m00pp + m10pp*pN2_q2 + m20pp*rN2_q2

    residual = 1 - c2 - c1 - c0

    # Also compute 1 - c₂ - c₁ - c₀ relative to max(|c_i|)
    scale = max(abs(c0), abs(c1), abs(c2))
    rel = abs(residual) / scale if scale != 0 else abs(residual)

    print(f"  N={N:2d}: 1-c₂-c₁-c₀ = {nstr(residual, 12)},  |rel| = {float(rel):.3e}")

# Also check: does v_k = q_{k+1} - q_k still grow as fast as q_k?
# If (S-1) were a factor, v_k would satisfy a lower-order recurrence
# with SLOWER growth.
print("\n=== Ratio |v_k/q_k| (should → 0 if (S-1) mode is neutral) ===")
q = [s[0] for s in s_vals]
for k in [5, 10, 15, 20, 25, 30]:
    if k+1 < len(q) and q[k] != 0:
        ratio = abs(q[k+1] - q[k]) / abs(q[k])
        print(f"  k={k}: |v_k/q_k| = {nstr(ratio, 15)}")
