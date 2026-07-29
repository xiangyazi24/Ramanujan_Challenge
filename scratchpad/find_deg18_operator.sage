#!/usr/bin/env sage
"""
Find the degree-18 polynomial-coefficient operator annihilating q_n (P2.7).

The recurrence eq:rec has ratio form with different A denominators.
The standard polynomial form should have order 3, degree 18.

Approach: compute many q_n values, then use guess() or matrix kernel.
"""
from ore_algebra import OreAlgebra, guess
import mpmath
mpmath.mp.dps = 100

# ============================================================
# Compute q_n with CORRECT recurrence
# ============================================================
N = 120  # Plenty for guessing degree 18

def A_p27(n):
    return QQ(1024) * (2*n+5)^4 * (2*n+7)^3 * (2*n+9)^3 * (946*n^2+6407*n+10860)

def B_p27(n):
    return QQ(128) * (2*n+7)^3 * (2*n+9)^3 * (104060*n^6 + 1745370*n^5 +
        12145238*n^4 + 44886481*n^3 + 92943995*n^2 + 102256019*n + 46709052)

def C_p27(n):
    return QQ(16) * (n+3)^4 * (2*n+9)^3 * (3784*n^5 + 57792*n^4 +
        351019*n^3 + 1059230*n^2 + 1587211*n + 944620)

def D_p27(n):
    return QQ(1) * (n+3)^4 * (n+4)^6 * (946*n^2 + 4515*n + 5399)

q = [QQ(0)] * N
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)

for n in range(2, N-1):
    q[n+1] = B_p27(n)/A_p27(n) * q[n] - C_p27(n-1)/A_p27(n-1) * q[n-1] + D_p27(n-2)/A_p27(n-2) * q[n-2]

print("q_0 = %s" % q[0])
print("q_3 = %s" % q[3])
print("q_10 = %s (approx %e)" % (q[10], float(q[10])))
print("q_20 (approx) = %e" % float(q[20]))

# ============================================================
# Method 1: guess() with ore_algebra
# ============================================================
Rn = PolynomialRing(QQ, 'nn')
nn_var = Rn.gen()
OS = OreAlgebra(Rn, 'Snn')
Snn = OS.gen()

print("\n=== Method 1: guess() ===")
q_for_guess = q[:80]

# Try order 3
for max_deg in [18, 20, 25, 30]:
    try:
        rec = guess(q_for_guess, OS, order=3, degree=max_deg)
        print("Found at order=3, max_deg=%d!" % max_deg)
        for j in range(4):
            print("  P_%d: degree %d" % (j, rec[j].degree()))
        # Print factored form of leading and trailing
        print("  P_3 = %s" % factor(rec[3]))
        print("  P_0 = %s" % factor(rec[0]))
        break
    except Exception as e:
        print("  order=3, max_deg=%d: %s" % (max_deg, e))

# ============================================================
# Method 2: Matrix kernel for specific degree
# ============================================================
print("\n=== Method 2: matrix kernel for degree 18 ===")
target_deg = 18
num_coeffs = 4 * (target_deg + 1)  # 4 shift terms × 19 polynomial coeffs = 76

rows = []
for n_val in range(80):
    row = []
    for j_shift in range(4):
        idx = n_val + j_shift
        if idx < N:
            for d in range(target_deg + 1):
                row.append(QQ(n_val**d) * q[idx])
        else:
            for d in range(target_deg + 1):
                row.append(QQ(0))
    rows.append(row)

M = matrix(QQ, rows)
K = M.right_kernel()
print("Kernel dimension: %d (want 1)" % K.dimension())

if K.dimension() >= 1:
    v = K.basis()[0]
    P_poly = []
    for j_shift in range(4):
        poly = sum(v[(target_deg+1)*j_shift + d] * nn_var^d for d in range(target_deg+1))
        P_poly.append(poly)
        print("  P_%d: degree %d" % (j_shift, poly.degree()))

    # Verify
    print("\n  Verification:")
    for n_val in range(5):
        val = sum(P_poly[j](nn=n_val) * q[n_val+j] for j in range(4))
        print("    n=%d: residual = %s" % (n_val, val))

    # Factor leading and trailing coefficients
    print("\n  P_3 factored: %s" % factor(P_poly[3]))
    print("  P_0 factored: %s" % factor(P_poly[0]))

    # Check if it also annihilates p_n
    p = [QQ(0)] * N
    p[0] = QQ(-612218384750)
    p[1] = QQ(-9525021973931919) / QQ(18100)
    p[2] = QQ(-29561828382772029) / QQ(65380)
    for n in range(2, N-1):
        p[n+1] = B_p27(n)/A_p27(n) * p[n] - C_p27(n-1)/A_p27(n-1) * p[n-1] + D_p27(n-2)/A_p27(n-2) * p[n-2]

    print("\n  P operator on p_n:")
    for n_val in range(5):
        val = sum(P_poly[j](nn=n_val) * p[n_val+j] for j in range(4))
        print("    n=%d: residual = %s" % (n_val, val))

    # ============================================================
    # GCRD with rec_a
    # ============================================================
    print("\n=== GCRD with rec_a ===")
    a_list = []
    for i in range(N):
        val = QQ(0)
        for k in range(i+1):
            val += binomial(i,k)^2 * binomial(i+k,i) * binomial(i+2*k,i)
        a_list.append(val)

    rec_a_op = guess(a_list[:60], OS, order=3)
    print("rec_a:")
    for j in range(4):
        print("  P_%d: degree %d" % (j, rec_a_op[j].degree()))

    # Build L_P27 as ore_algebra operator
    L_P27 = sum(P_poly[j] * Snn^j for j in range(4))
    print("\nL_P27 = %s" % L_P27)

    try:
        gcrd = L_P27.gcrd(rec_a_op)
        print("GCRD order: %d" % gcrd.order())
        if gcrd.order() > 0:
            print("GCRD = %s" % gcrd)
    except Exception as e:
        print("GCRD failed: %s" % e)

    try:
        lclm = L_P27.lclm(rec_a_op)
        print("LCLM order: %d" % lclm.order())
    except Exception as e:
        print("LCLM failed: %s" % e)

    # ============================================================
    # Convert to theta/D operators and find GF ODE
    # ============================================================
    print("\n=== GF ODE from P2.7 ===")
    try:
        theta_op_p27 = L_P27.to_T('Tz')
        print("Theta operator order: %d" % theta_op_p27.order())

        diff_op_p27 = L_P27.to_D('Dz')
        print("Diff operator order: %d" % diff_op_p27.order())

        # Indicial polynomial
        ind = diff_op_p27.indicial_polynomial(0)
        print("Indicial at z=0: %s" % ind)
        print("Roots: %s" % ind.roots(QQbar))
    except Exception as e:
        import traceback
        print("Conversion error:")
        traceback.print_exc()

    # ============================================================
    # GCRD in GF ODE form
    # ============================================================
    print("\n=== GCRD of GF ODEs ===")
    try:
        theta_a = rec_a_op.to_T('Tz')
        gcrd_gf = theta_op_p27.gcrd(theta_a)
        print("GCRD of GF theta operators: order %d" % gcrd_gf.order())
    except Exception as e:
        print("GF GCRD failed: %s" % e)
