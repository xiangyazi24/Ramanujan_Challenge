#!/usr/bin/env python3
"""P2.5: Extract the minimal ODE (θ-order 5, z-degree 11) for F̂(z).

Found by p25_euler_ode.py. Extract the null vector, clear denominators,
and analyze the ODE (singular points, local exponents).
"""

def M_entries(n):
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

def delta_int(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

from fractions import Fraction
from math import gcd as igcd

NMAX = 300
print(f"Computing Q̂_n (exact rational) for n = 0..80, mod P for n = 0..{NMAX}...")

# Exact computation for first 81 terms
q_exact = [33750, -36000, 9000]
Q0_int = [33750]
for N in range(80):
    M = M_entries(N)
    new_q = [0, 0, 0]
    for j in range(3):
        for k in range(3):
            new_q[j] += q_exact[k] * M[k][j]
    q_exact = new_q
    Q0_int.append(q_exact[0])

def H_frac(n):
    result = Fraction(1)
    for k in range(n):
        result *= Fraction(-16)
        result *= Fraction(k+2)**2
        result *= Fraction(k+3)**2
        result *= Fraction(2*k+5, 2)
        result *= Fraction(2*k+7, 2)**2
    return result

Q_hat_exact = []
for n in range(81):
    Hn = H_frac(n)
    Q_hat_exact.append(Fraction(Q0_int[n]) / Hn if Hn != 0 else Fraction(0))

# Modular computation for all terms
P = (1 << 61) - 1
def modinv(a, p=P):
    return pow(a % p, p - 2, p)

q_mod = [33750 % P, (-36000) % P, 9000 % P]
Q0_mod = [33750 % P]
H_mod = [1]
h_val = 1
for N in range(NMAX):
    M = M_entries(N)
    new_q = [0, 0, 0]
    for j in range(3):
        for k in range(3):
            new_q[j] = (new_q[j] + q_mod[k] * M[k][j]) % P
    q_mod = new_q
    Q0_mod.append(q_mod[0])
    d = delta_int(N) % P
    h_val = (h_val * d) % P
    H_mod.append(h_val)

Q_hat_mod = [(Q0_mod[n] * modinv(H_mod[n])) % P if H_mod[n] != 0 else 0
             for n in range(NMAX + 1)]
print(f"Done: {len(Q_hat_exact)} exact, {len(Q_hat_mod)} modular terms")

# Build and solve the Euler ODE system exactly
theta_order = 5
z_deg = 11
D = z_deg
order = theta_order
num_unknowns = (order + 1) * (D + 1)
print(f"\nExtracting ODE: θ-order={order}, z-deg={D}, {num_unknowns} unknowns")

# Use modular arithmetic to find the null vector structure first
m_min = D
m_max = NMAX
use_eqs = num_unknowns + 5

mat_mod = []
for eq_idx in range(use_eqs):
    m = m_min + eq_idx
    row = [0] * num_unknowns
    for k in range(order + 1):
        for j in range(D + 1):
            col = k * (D + 1) + j
            idx = m - j
            if 0 <= idx < len(Q_hat_mod):
                coeff = pow(m - j, k, P) * Q_hat_mod[idx] % P
                row[col] = coeff
    mat_mod.append(row)

# Gaussian elimination mod P, track pivots
mat_work = [row[:] for row in mat_mod]
r = 0
pivot_cols = []
for c in range(num_unknowns):
    found = -1
    for i in range(r, len(mat_work)):
        if mat_work[i][c] % P != 0:
            found = i
            break
    if found < 0:
        continue
    mat_work[r], mat_work[found] = mat_work[found], mat_work[r]
    pivot_cols.append(c)
    inv = pow(mat_work[r][c], P - 2, P)
    for j_col in range(num_unknowns):
        mat_work[r][j_col] = (mat_work[r][j_col] * inv) % P
    for i in range(len(mat_work)):
        if i != r and mat_work[i][c] % P != 0:
            factor = mat_work[i][c]
            for j_col in range(num_unknowns):
                mat_work[i][j_col] = (mat_work[i][j_col] - factor * mat_work[r][j_col]) % P
    r += 1

rank = r
free_cols = [c for c in range(num_unknowns) if c not in pivot_cols]
print(f"Modular rank={rank}, free cols={free_cols}")

# Extract modular null vector
x_mod = [0] * num_unknowns
fc = free_cols[0]
x_mod[fc] = 1
for pr_idx in range(len(pivot_cols) - 1, -1, -1):
    pc = pivot_cols[pr_idx]
    s = 0
    for j_col in range(pc + 1, num_unknowns):
        s = (s + mat_work[pr_idx][j_col] * x_mod[j_col]) % P
    x_mod[pc] = (-s) % P

# Verify modular null vector on extra equations
print("\nVerifying modular null vector...")
for m in [50, 100, 200, 290]:
    val = 0
    for k in range(order + 1):
        for j in range(D + 1):
            col = k * (D + 1) + j
            idx = m - j
            if 0 <= idx < len(Q_hat_mod):
                coeff = pow(m - j, k, P) * Q_hat_mod[idx] % P
                val = (val + x_mod[col] * coeff) % P
    print(f"  m={m}: residual = {val}")

# Now extract exact ODE using rational reconstruction
# Build exact system with first 72 equations (using exact Q̂_n for n <= 80)
print("\nBuilding exact system...")
exact_unknowns = num_unknowns
exact_eqs = min(72 - D, exact_unknowns + 5)
m_min_exact = D

mat_exact = []
for eq_idx in range(exact_eqs):
    m = m_min_exact + eq_idx
    row = [Fraction(0)] * exact_unknowns
    for k in range(order + 1):
        for j in range(D + 1):
            col = k * (D + 1) + j
            idx = m - j
            if 0 <= idx < len(Q_hat_exact):
                coeff = Fraction(m - j) ** k * Q_hat_exact[idx]
                row[col] = coeff
    mat_exact.append(row)

# Exact Gaussian elimination
mat_e = [row[:] for row in mat_exact]
r_e = 0
pivot_cols_e = []
for c in range(exact_unknowns):
    found = -1
    for i in range(r_e, len(mat_e)):
        if mat_e[i][c] != Fraction(0):
            found = i
            break
    if found < 0:
        continue
    mat_e[r_e], mat_e[found] = mat_e[found], mat_e[r_e]
    pivot_cols_e.append(c)
    pv = mat_e[r_e][c]
    for j_col in range(exact_unknowns):
        mat_e[r_e][j_col] /= pv
    for i in range(len(mat_e)):
        if i != r_e and mat_e[i][c] != Fraction(0):
            f = mat_e[i][c]
            for j_col in range(exact_unknowns):
                mat_e[i][j_col] -= f * mat_e[r_e][j_col]
    r_e += 1

exact_rank = r_e
free_cols_e = [c for c in range(exact_unknowns) if c not in pivot_cols_e]
print(f"Exact rank={exact_rank}, free cols={free_cols_e}")

if len(free_cols_e) == 0:
    print("NO exact null vector found — need more terms or different approach")
else:
    # Extract null vector
    fc_e = free_cols_e[0]
    x_exact = [Fraction(0)] * exact_unknowns
    x_exact[fc_e] = Fraction(1)
    for pr_idx in range(len(pivot_cols_e) - 1, -1, -1):
        pc = pivot_cols_e[pr_idx]
        s = sum(mat_e[pr_idx][j_col] * x_exact[j_col] for j_col in range(pc + 1, exact_unknowns))
        x_exact[pc] = -s

    # Clear denominators
    lcm_den = 1
    for c in x_exact:
        if c != 0:
            lcm_den = lcm_den * c.denominator // igcd(lcm_den, c.denominator)
    x_int = [int(c * lcm_den) for c in x_exact]
    g = 0
    for c in x_int:
        g = igcd(g, abs(c))
    if g > 0:
        x_int = [c // g for c in x_int]

    # Verify on modular data
    print("\nVerifying exact null vector on modular data...")
    ok = True
    for m in range(D, min(len(Q_hat_mod), 295)):
        val = 0
        for k in range(order + 1):
            for j in range(D + 1):
                col = k * (D + 1) + j
                idx = m - j
                if 0 <= idx < len(Q_hat_mod):
                    coeff = pow(m - j, k, P) * Q_hat_mod[idx] % P
                    val = (val + x_int[col] * coeff) % P
        if val != 0:
            print(f"  FAIL at m={m}: residual = {val}")
            ok = False
            break
    if ok:
        print(f"  VERIFIED on m={D}..{min(len(Q_hat_mod)-1, 294)}")

    # Display the ODE coefficients
    print("\n" + "="*60)
    print("MINIMAL ODE: Σ_{k=0}^5 q_k(z) θ^k F(z) = 0")
    print("="*60)
    for k in range(order + 1):
        coeffs = x_int[k * (D + 1):(k + 1) * (D + 1)]
        deg = 0
        for j in range(len(coeffs) - 1, -1, -1):
            if coeffs[j] != 0:
                deg = j
                break
        print(f"\nq_{k}(z) [z-degree {deg}]:")
        for j in range(deg + 1):
            if coeffs[j] != 0:
                print(f"  z^{j}: {coeffs[j]}")

    # Find singular points: roots of leading coefficient q_5(z)
    print("\n" + "="*60)
    print("Singular points (roots of q_5(z))")
    print("="*60)
    q5_coeffs = x_int[5 * (D + 1):6 * (D + 1)]
    print("q_5(z) coefficients:", q5_coeffs)

    try:
        from sympy import Symbol, Poly, factor, roots as sym_roots, Rational, sqrt, simplify
        z = Symbol('z')
        q5_poly = sum(c * z**j for j, c in enumerate(q5_coeffs) if c != 0)
        print(f"\nq_5(z) = {q5_poly}")
        q5_factored = factor(q5_poly)
        print(f"factored: {q5_factored}")
        rts = sym_roots(q5_poly, z)
        print(f"\nRoots of q_5:")
        for rt, mult in rts.items():
            print(f"  z = {rt} (mult {mult}), numerical = {complex(rt):.10f}")
    except Exception as e:
        print(f"Sympy factoring failed: {e}")
