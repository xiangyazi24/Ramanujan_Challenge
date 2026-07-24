#!/usr/bin/env python3
"""P2.5: Analyze convergence rate |E_{N,j}/Q_{N,j}| and compare with ρ^N.

The key question: is |G*Q_N - P_N| / |Q_N| = O(ρ^N) or O(N^{-3} ρ^N)?

Also check if E/Q ratios match the Legendre neutral model prediction.
"""
from mpmath import mp, mpf, matrix, catalan, log, sqrt

mp.dps = 80
G = catalan

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
    return matrix([
        [mpf(m11), mpf(m12), mpf(m13)],
        [mpf(m21), mpf(m22), mpf(m23)],
        [mpf(m31), mpf(m32), mpf(m33)],
    ])

def delta(n):
    """Pochhammer normalization factor."""
    return mpf(-2) * (n+2)**2 * (n+3)**2 * (2*n+5) * (2*n+7)**2

def M_normalized(n):
    """M_H(n) = M(n) / delta(n)."""
    return M_exact(n) / delta(n)

p0 = [mpf(30921), mpf(-32972), mpf(8240)]
q0 = [mpf(33750), mpf(-36000), mpf(9000)]

sqrt2 = sqrt(2)
lam_plus = 17 + 12*sqrt2
rho = 17 - 12*sqrt2

print(f"ρ = {mp.nstr(rho, 20)}")
print(f"log₁₀(ρ) = {mp.nstr(mp.log10(rho), 10)}")
print(f"Expected digits/step: {mp.nstr(-mp.log10(rho), 5)}")
print()

Nmax = 50

# Use NORMALIZED matrix for cleaner analysis
# With normalized M_H(n), the Poincaré roots are the actual eigenvalues at infinity.

# Forward iteration with raw matrix (simpler, just track ratios)
p_row = list(p0)
q_row = list(q0)

print(f"{'N':>3} {'|E/Q| col 0':>18} {'|E/Q|/ρ^N':>18} {'|E/Q|·N³/ρ^N':>18} {'digits':>8}")
print("="*75)

for N in range(Nmax+1):
    P_vals = list(p_row)
    Q_vals = list(q_row)

    for j in [0]:
        E_j = G * Q_vals[j] - P_vals[j]
        if Q_vals[j] != 0 and N > 0:
            ratio = abs(E_j / Q_vals[j])
            rho_N = rho**N
            normalized = ratio / rho_N
            norm_N3 = normalized * (N+1)**3
            digits = -mp.log10(ratio) if ratio > 0 else mpf(0)

            if N <= 10 or N % 5 == 0:
                print(f"{N:3d} {mp.nstr(ratio, 10):>18s} {mp.nstr(normalized, 10):>18s} "
                      f"{mp.nstr(norm_N3, 10):>18s} {mp.nstr(digits, 5):>8s}")

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

# Now do the normalized analysis with M_H
print()
print("="*75)
print("NORMALIZED MATRIX ANALYSIS (M_H = M/delta)")
print("="*75)

# With normalized matrix, forward product gives "Poincaré-scaled" sequences
# But we need to be careful about the cumulative scaling
# prod(delta(k), k=0..N-1) cancels the Pochhammer factor H_N

# Let me compute with the normalized matrix directly
p_norm = list(p0)
q_norm = list(q0)
cum_delta = mpf(1)  # cumulative delta product

print(f"\n{'N':>3} {'|Ê/Q̂| (norm)':>18} {'|Ê/Q̂|/ρ^N':>18} {'|Ê/Q̂|·N³/ρ^N':>18}")
print("-"*65)

for N in range(Nmax+1):
    P_vals = list(p_norm)
    Q_vals = list(q_norm)

    for j in [0]:
        E_j = G * Q_vals[j] - P_vals[j]
        if Q_vals[j] != 0 and N > 0:
            ratio = abs(E_j / Q_vals[j])
            rho_N = rho**N
            normalized = ratio / rho_N
            norm_N3 = normalized * (N+1)**3

            if N <= 10 or N % 5 == 0:
                print(f"{N:3d} {mp.nstr(ratio, 10):>18s} {mp.nstr(normalized, 10):>18s} "
                      f"{mp.nstr(norm_N3, 10):>18s}")

    if N < Nmax:
        MH = M_normalized(N)
        new_p = [mpf(0)]*3
        new_q = [mpf(0)]*3
        for col in range(3):
            for k in range(3):
                new_p[col] += p_norm[k] * MH[k, col]
                new_q[col] += q_norm[k] * MH[k, col]
        p_norm = new_p
        q_norm = new_q

# Cross-column analysis: are the column errors proportional?
print()
print("="*75)
print("CROSS-COLUMN ERROR RATIOS (should be N-independent if neutral mode dominates)")
print("="*75)

p_row = list(p0)
q_row = list(q0)

for N in range(41):
    if N > 0 and (N <= 5 or N % 5 == 0):
        E_cols = [G * q_row[j] - p_row[j] for j in range(3)]
        if E_cols[0] != 0:
            r1 = E_cols[1] / E_cols[0]
            r2 = E_cols[2] / E_cols[0]
            print(f"  N={N:3d}: E₁/E₀ = {mp.nstr(r1, 20)}, E₂/E₀ = {mp.nstr(r2, 20)}")

    if N < 40:
        M = M_exact(N)
        new_p = [mpf(0)]*3
        new_q = [mpf(0)]*3
        for col in range(3):
            for k in range(3):
                new_p[col] += p_row[k] * M[k, col]
                new_q[col] += q_row[k] * M[k, col]
        p_row = new_p
        q_row = new_q
