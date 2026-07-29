#!/usr/bin/env python3
"""
Search for non-scalar Ore intertwiner: cyclic-vector map from Zudilin to P2.7.

Ansatz: q̂_n = 64^n * q_n = u₀(n)*b_n + u₁(n)*b_{n+1} + u₂(n)*b_{n+2}
where u_j(n) ∈ Q(n) are rational functions.

Similarly for the numerator:
p̂_n = 64^n * p_n = v₀(n)*(b̃_n + b̃̃_n) + v₁(n)*(b̃_{n+1} + b̃̃_{n+1}) + v₂(n)*(b̃_{n+2} + b̃̃_{n+2})
                  + w₀(n)*b_n + w₁(n)*b_{n+1} + w₂(n)*b_{n+2}

If found, this transfers Zudilin's Barnes error bounds to P2.7.
"""
from fractions import Fraction as F
from mpmath import mp, mpf, matrix, lu_solve, log10, fabs

mp.dps = 100

def frac_to_mpf(x):
    """Convert Fraction to mpf."""
    if isinstance(x, F):
        return mpf(x.numerator) / mpf(x.denominator)
    return mpf(x)

# ===== Zudilin recurrence for b_n (eq 6.4 from Zudilin's paper) =====
def zudilin_terms(initial, N):
    """Compute terms of Zudilin's recurrence."""
    u = list(initial)
    for k in range(2, N):
        n = k
        QZ = 946*n*n - 731*n + 153
        MZ = 104060*n**6 + 127710*n**5 + 12788*n**4 - 34525*n**3 - 8482*n**2 + 3298*n + 1071
        NZ = 3784*n**5 - 1032*n**4 - 1925*n**3 + 853*n**2 + 328*n - 184
        RZ = 946*n*n + 1161*n + 368

        nxt = (2*MZ*u[k] - 2*n*NZ*u[k-1] + RZ*n*(n-1)**3*u[k-2])
        denom = 2*QZ*(2*n+1)*(n+1)**3
        u.append(F(nxt, denom))
    return u

# ===== P2.7 recurrence =====
def A_c(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860)
def B_c(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_c(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_c(n): return (n+3)**4*(n+4)**6*(946*n*n+4515*n+5399)

def p27_terms(initial, N):
    u = list(initial)
    for k in range(2, N):
        n = k
        nxt = F(B_c(n)*u[n] - C_c(n-1)*u[n-1]*A_c(n)//A_c(n-1) + D_c(n-2)*u[n-2]*A_c(n)//A_c(n-2), A_c(n))
        # More carefully:
        nxt = F(B_c(n), A_c(n)) * u[n] - F(C_c(n-1), A_c(n-1)) * u[n-1] + F(D_c(n-2), A_c(n-2)) * u[n-2]
        u.append(nxt)
    return u

N = 30

# Zudilin solutions
b = zudilin_terms([F(1), F(7), F(163)], N)
bt = zudilin_terms([F(0), F(23,2), F(2145,8)], N)   # b̃ (ζ(2) companion)
btt = zudilin_terms([F(0), F(17,2), F(3135,16)], N)  # b̃̃ (ζ(3) companion)

# P2.7 solutions
q = p27_terms([F(-215040420000), F(-167282265043404, 905), F(-964185327658080, 6071)], N)
p = p27_terms([F(-612218384750), F(-9525021973931919, 18100), F(-29561828382772029, 65380)], N)

# Scaled P2.7: q̂_n = 64^n * q_n
qhat = [F(64)**n * q[n] for n in range(N)]
phat = [F(64)**n * p[n] for n in range(N)]

print("=== Verify sequences ===")
print(f"b[0:5] = {[float(x) for x in b[:5]]}")
print(f"q[0:3] = {[float(x) for x in q[:3]]}")
print(f"qhat[0:3] = {[float(x) for x in qhat[:3]]}")

# ===== Search for cyclic-vector map =====
# Ansatz: qhat_n = u0(n)*b_n + u1(n)*b_{n+1} + u2(n)*b_{n+2}
# If u_j are polynomial of degree d, we need 3*(d+1) unknowns, so 3*(d+1) equations.

def search_poly_intertwiner(target, sources_list, max_deg=6):
    """
    Search for polynomial coefficients u_j(n) such that
    target[n] = sum_j u_j(n) * sources_list[j][n]

    u_j(n) = sum_{k=0}^d c_{j,k} * n^k
    """
    for deg in range(max_deg + 1):
        dim = len(sources_list) * (deg + 1)
        if dim > N - 2:
            print(f"  deg {deg}: not enough data points ({N-2} < {dim})")
            continue

        # Build matrix: row n gives the coefficients of c_{j,k}
        # target[n] = sum_j sum_k c_{j,k} * n^k * sources_list[j][n]
        rows = min(N - 2, dim + 5)  # overdetermined for verification

        A_mat = matrix(rows, dim)
        b_vec = matrix(rows, 1)

        for i in range(rows):
            n = i
            b_vec[i, 0] = frac_to_mpf(target[n])
            col = 0
            for j in range(len(sources_list)):
                for k in range(deg + 1):
                    A_mat[i, col] = mpf(n**k) * frac_to_mpf(sources_list[j][n])
                    col += 1

        # Solve least-squares (use first dim equations, check rest)
        A_sq = A_mat[:dim, :]
        b_sq = b_vec[:dim, :]

        try:
            x = lu_solve(A_sq, b_sq)
        except:
            print(f"  deg {deg}: singular system")
            continue

        # Check residuals on ALL rows
        residuals = A_mat * x - b_vec
        max_res = max(fabs(residuals[i, 0]) for i in range(rows))

        if max_res < mpf(10)**(-50):
            print(f"  deg {deg}: EXACT MATCH! max residual = {float(log10(max_res + 1e-200)):.1f} digits")
            # Print the coefficients
            col = 0
            for j in range(len(sources_list)):
                coeffs = []
                for k in range(deg + 1):
                    coeffs.append(x[col, 0])
                    col += 1
                print(f"    u_{j}(n) = {' + '.join(f'{float(c):.6g}*n^{k}' for k, c in enumerate(coeffs))}")
            return x
        else:
            rel = float(log10(max_res + 1e-200))
            print(f"  deg {deg}: max residual ~ 10^{rel:.0f}")

    return None

# Search 1: qhat_n = u0(n)*b_n + u1(n)*b_{n+1} + u2(n)*b_{n+2}
print("\n=== Search: qhat_n = u0(n)*b_n + u1(n)*b_{n+1} + u2(n)*b_{n+2} ===")
sources_shift = [[b[n] for n in range(N)],
                 [b[n+1] if n+1 < N else F(0) for n in range(N)],
                 [b[n+2] if n+2 < N else F(0) for n in range(N)]]
search_poly_intertwiner(qhat, sources_shift, max_deg=8)

# Search 2: qhat_n = u0(n)*b_n + u1(n)*b_{n-1} + u2(n)*b_{n-2}  (backward shifts)
print("\n=== Search: qhat_n = u0(n)*b_n + u1(n)*b_{n-1} + u2(n)*b_{n-2} ===")
sources_back = [[b[n] for n in range(N)],
                [b[n-1] if n >= 1 else F(0) for n in range(N)],
                [b[n-2] if n >= 2 else F(0) for n in range(N)]]
search_poly_intertwiner([qhat[n] for n in range(N)], sources_back, max_deg=8)

# Search 3: wider window — qhat_n from b_{n-2},...,b_{n+2}
print("\n=== Search: qhat_n from b_{n-2},...,b_{n+2} (5 sources) ===")
sources_wide = []
for shift in range(-2, 3):
    sources_wide.append([b[n+shift] if 0 <= n+shift < N else F(0) for n in range(N)])
search_poly_intertwiner([qhat[n] for n in range(N)], sources_wide, max_deg=4)

# Search 4: rational coefficients — try u_j(n) = P_j(n)/Q(n) with common denominator
# Let's try Q(n) = A_c(n)/1024... or simpler denominators
print("\n=== Search: qhat_n = u0(n)*b_n + u1(n)*b_{n+1} with RATIONAL u_j ===")
# Try q_n * A_c(n) = w0(n)*b_n + w1(n)*b_{n+1} + w2(n)*b_{n+2}  (polynomial w_j)
# i.e., multiply both sides by A_c(n) to clear denominators
qhat_times_A = [qhat[n] * A_c(n) for n in range(N)]
sources_times_nothing = [[b[n] for n in range(N)],
                          [b[n+1] if n+1 < N else F(0) for n in range(N)],
                          [b[n+2] if n+2 < N else F(0) for n in range(N)]]
print("  (clearing denominator by A_c(n))")
search_poly_intertwiner(qhat_times_A, sources_times_nothing, max_deg=6)

print("\n=== P2.6 series connection ===")
# P2.6: 2077/720 + sum u_j = zeta(2) + zeta(3)
# u_n satisfies: -2(n+3)^3(2n+5)(3n+5) u_n + (n+2)^2(15n^3+85n^2+155n+93) u_{n-1} - (n+1)^3(n+2)(3n+8) u_{n-2} = 0

def p26_terms(N):
    u = {1: F(-93, 4480), 2: F(-117, 14000)}
    for n in range(3, N):
        coeff_n = -2*(n+3)**3*(2*n+5)*(3*n+5)
        coeff_n1 = (n+2)**2*(15*n**3+85*n**2+155*n+93)
        coeff_n2 = -(n+1)**3*(n+2)*(3*n+8)
        # coeff_n * u_n + coeff_n1 * u_{n-1} + coeff_n2 * u_{n-2} = 0
        u[n] = F(-(coeff_n1 * u[n-1] + coeff_n2 * u[n-2]), coeff_n)
    return u

u26 = p26_terms(50)
partial_sums = [F(2077, 720)]
for j in range(1, 50):
    partial_sums.append(partial_sums[-1] + u26[j])

from mpmath import zeta
mp.dps = 50
L = zeta(2) + zeta(3)
print(f"ζ(2)+ζ(3) = {L}")
print(f"P2.6 partial sum S_10 = {float(partial_sums[10])}")
print(f"P2.6 partial sum S_20 = {float(partial_sums[20])}")
print(f"P2.6 partial sum S_30 = {float(partial_sums[30])}")
print(f"P2.6 S_30 - L = {float(frac_to_mpf(partial_sums[30]) - L)}")

# Check: is there a direct relationship between P2.6 terms and P2.7 sequences?
# E.g., q_n = some transform of partial sums?
print("\n=== Checking P2.6-P2.7 numerical relationships ===")
# Ratio q_n/S_n?
for n in range(1, 10):
    r = float(frac_to_mpf(q[n]) / frac_to_mpf(partial_sums[n]))
    print(f"  q_{n}/S_{n} = {r:.10e}")

# Ratio q_n/u26_n?
print()
for n in range(1, 10):
    r = float(frac_to_mpf(q[n]) / frac_to_mpf(u26[n]))
    print(f"  q_{n}/u_{n}(P2.6) = {r:.10e}")

# Check if the 946 in P2.7 Poincare also appears in P2.6
# P2.6 Poincare: 4λ² - 5λ + 1 = 0 → (4λ-1)(λ-1) = 0
# Discriminant: 25 - 16 = 9. No 946.
# But 946 = 2·11·43. Check if 43 appears in P2.6:
print("\n=== Factor analysis ===")
print(f"946 = 2 * 11 * 43")
for n in range(1, 10):
    val = int(u26[n].numerator)
    if val % 43 == 0:
        print(f"  u_{n} numerator {val} divisible by 43")
