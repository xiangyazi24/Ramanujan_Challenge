#!/usr/bin/env python3
"""P2.5: Numerical decomposition of the CMF error E_{N,j} = G*Q_{N,j} - P_{N,j}.

Tests whether E_{N,j} matches the Legendre neutral model V_N^0 = P_N(3)*Q_N(3)/(N+1)^2.
If E_{N,j} / V_N^0 -> d_0 (constant), then the convergence rate follows from the
explicit integral representation of V_N^0.

Also tests the direct decay: E_{N,j} * N^3 -> constant (neutral formal exponent = -3).
"""
from mpmath import mp, mpf, matrix, catalan, log, sqrt, pi

mp.dps = 80
G = catalan

def M_exact(n):
    """The 3x3 CMF matrix M(n)."""
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return matrix([
        [mpf(m11), mpf(m12), mpf(m13)],
        [mpf(m21), mpf(m22), mpf(m23)],
        [mpf(m31), mpf(m32), mpf(m33)],
    ])

def legendre_P(n, x):
    """Legendre P_n(x) via forward recurrence."""
    if n == 0: return mpf(1)
    if n == 1: return mpf(x)
    p0, p1 = mpf(1), mpf(x)
    for k in range(1, n):
        p2 = ((2*k+1)*x*p1 - k*p0) / (k+1)
        p0, p1 = p1, p2
    return p1

def legendre_Q(n, x):
    """Legendre Q_n(x) for x > 1 via forward recurrence."""
    x = mpf(x)
    q0 = log((x+1)/(x-1)) / 2
    if n == 0: return q0
    q1 = x * q0 - 1
    if n == 1: return q1
    for k in range(1, n):
        q2 = ((2*k+1)*x*q1 - k*q0) / (k+1)
        q0, q1 = q1, q2
    return q1

# Initial matrix A
p0 = [mpf(30921), mpf(-32972), mpf(8240)]
q0 = [mpf(33750), mpf(-36000), mpf(9000)]

# Poincaré roots
sqrt2 = sqrt(2)
lam_plus = 17 + 12*sqrt2
rho = 17 - 12*sqrt2

print(f"Precision: {mp.dps} digits")
print(f"G = {mp.nstr(G, 40)}")
print(f"λ₊ = {mp.nstr(lam_plus, 20)}")
print(f"ρ = {mp.nstr(rho, 20)}")
print()

# Compute CMF matrix products and extract P, Q, E sequences
Nmax = 50

# M_N = M(0) * M(1) * ... * M(N-1)
# A * M_N = [[P_{N,1}, P_{N,2}, P_{N,3}],
#            [Q_{N,1}, Q_{N,2}, Q_{N,3}]]
# P_{N,j} = p0 * M_N[:,j], Q_{N,j} = q0 * M_N[:,j]

# At N=0: A * I = A, so P_{0,j} = p0[j], Q_{0,j} = q0[j]

# Forward iteration: maintain p_row and q_row as 1x3 row vectors
# p_row(N) = p0 * M(0) * M(1) * ... * M(N-1)
p_row = list(p0)
q_row = list(q0)

print("="*80)
print(f"{'N':>3} {'E_{N,1}':>20} {'E_{N,1}·N³':>20} {'E_{N,1}/V⁰_N':>20} {'P/Q ratio':>20}")
print("="*80)

results = []

for N in range(Nmax+1):
    # Current state: p_row = p0 * M(0) * ... * M(N-1) at this N
    P_vals = list(p_row)
    Q_vals = list(q_row)

    # Error for column j=0 (first column)
    j = 0
    E_j = G * Q_vals[j] - P_vals[j]

    # Legendre model
    Pleg = legendre_P(N, 3)
    Qleg = legendre_Q(N, 3)
    V0 = Pleg * Qleg / (N+1)**2  # neutral model

    # Ratios
    ratio_N3 = E_j * (N+1)**3 if N > 0 else mpf(0)
    ratio_V0 = E_j / V0 if V0 != 0 else mpf(0)
    PQ_ratio = P_vals[j] / Q_vals[j] if Q_vals[j] != 0 else mpf(0)

    results.append((N, E_j, ratio_N3, ratio_V0, PQ_ratio))

    if N <= 5 or N % 5 == 0:
        print(f"{N:3d} {mp.nstr(E_j, 12):>20s} {mp.nstr(ratio_N3, 12):>20s} "
              f"{mp.nstr(ratio_V0, 12):>20s} {mp.nstr(PQ_ratio, 12):>20s}")

    # Advance: multiply by M(N) on the right
    if N < Nmax:
        M = M_exact(N)
        new_p = [mpf(0)]*3
        new_q = [mpf(0)]*3
        for col in range(3):
            for k in range(3):
                new_p[col] += p_row[k] * M[k, col]
                new_q[col] += q_row[k] * M[k, col]
        p_row = new_p
        q_row = new_q

print()
print("="*80)
print("Analysis of E_{N,j} / V_N^0 convergence (should stabilize if decomposition exists):")
print("="*80)
for N, E_j, ratio_N3, ratio_V0, PQ_ratio in results:
    if N >= 5 and N <= 40 and N % 5 == 0:
        # Also compute dominant normalization
        E_norm = E_j / lam_plus**N
        print(f"  N={N:3d}: E/V⁰ = {mp.nstr(ratio_V0, 25)},  E·N³ = {mp.nstr(ratio_N3, 25)},  "
              f"E/λ₊^N = {mp.nstr(E_norm, 15)}")

# Check all three columns
print()
print("="*80)
print("All three columns at N=30:")
print("="*80)

# Recompute for all columns at specific N
p_row2 = list(p0)
q_row2 = list(q0)
target_N = 30

for N in range(target_N):
    M = M_exact(N)
    new_p = [mpf(0)]*3
    new_q = [mpf(0)]*3
    for col in range(3):
        for k in range(3):
            new_p[col] += p_row2[k] * M[k, col]
            new_q[col] += q_row2[k] * M[k, col]
    p_row2 = new_p
    q_row2 = new_q

Pleg30 = legendre_P(target_N, 3)
Qleg30 = legendre_Q(target_N, 3)
V0_30 = Pleg30 * Qleg30 / (target_N+1)**2

for j in range(3):
    E_j = G * q_row2[j] - p_row2[j]
    ratio_V0 = E_j / V0_30
    ratio_N3 = E_j * (target_N+1)**3
    print(f"  j={j}: E/V⁰ = {mp.nstr(ratio_V0, 30)}, E·N³ = {mp.nstr(ratio_N3, 30)}")

# Check: do the three columns have proportional errors?
print()
E_cols = [G * q_row2[j] - p_row2[j] for j in range(3)]
if E_cols[0] != 0:
    print(f"Column ratios at N={target_N}:")
    print(f"  E_1/E_0 = {mp.nstr(E_cols[1]/E_cols[0], 30)}")
    print(f"  E_2/E_0 = {mp.nstr(E_cols[2]/E_cols[0], 30)}")
