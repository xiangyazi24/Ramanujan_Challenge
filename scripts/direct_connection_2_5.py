#!/usr/bin/env python3
"""
Problem 2.5: Direct connection between CMF and Sym²(Delannoy).

Test: does q(N) = r₀(N)·D_N² + r₁(N)·D_N·E_N + r₂(N)·E_N²
where r₀, r₁, r₂ are rational functions of N?

More precisely, compute the 3×3 connection matrix G(N) = Π(N)·Φ(N)⁻¹
where Π(N) is the CMF product matrix and Φ(N) is the Sym² Delannoy
fundamental matrix. If G(N) has entries that are rational functions of N,
that's the certificate.

After factorial gauge: t(N) = q(N) / [(-16)^N · (N!)^7],
check if G̃(N) = gauged_Π(N) · Φ(N)⁻¹ has rational entries.
"""
from fractions import Fraction as F
from math import factorial

def M_int(n):
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)
    return [[F(m11), F(m12), F(m13)], [F(m21), F(m22), F(m23)], [F(m31), F(m32), F(m33)]]

def mat_mul_F(A, B):
    n = len(A)
    return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

def mat_inv_F(M):
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, i = M[2]
    det = a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)
    if det == 0:
        return None
    return [[(e*i-f*h)/det, -(b*i-c*h)/det, (b*f-c*e)/det],
           [-(d*i-f*g)/det, (a*i-c*g)/det, -(a*f-c*d)/det],
           [(d*h-e*g)/det, -(a*h-b*g)/det, (a*e-b*d)/det]]

def sym2_companion(n):
    a = F(3*(2*n+1), n+1)
    b = F(-n, n+1)
    return [[a*a, 2*a*b, b*b],
            [a, b, F(0)],
            [F(1), F(0), F(0)]]

# Compute CMF product Π(N) and Sym² product Φ(N)
N_MAX = 25
I3 = [[F(1), F(0), F(0)], [F(0), F(1), F(0)], [F(0), F(0), F(1)]]

Pi = [row[:] for row in I3]  # CMF product
Phi = [row[:] for row in I3]  # Sym² product

print("=== Connection matrix G(N) = Π(N) · Φ(N)⁻¹ ===")
print("If G(N) entries are rational functions of N, we have the certificate.\n")

# Also compute Delannoy sequences for reference
D_seq = [F(1), F(3)]
E_seq = [F(0), F(1)]
for n in range(1, N_MAX + 5):
    D_seq.append((F(3*(2*n+1)) * D_seq[n] - F(n) * D_seq[n-1]) / F(n+1))
    E_seq.append((F(3*(2*n+1)) * E_seq[n] - F(n) * E_seq[n-1]) / F(n+1))

gauge_ratios = []

for N in range(N_MAX):
    if N > 0:
        Pi = mat_mul_F(Pi, M_int(N-1))
        Phi = mat_mul_F(Phi, sym2_companion(N-1))

    Phi_inv = mat_inv_F(Phi)
    if Phi_inv is None:
        print(f"N={N}: Φ singular!")
        continue

    G = mat_mul_F(Pi, Phi_inv)

    if N <= 8:
        print(f"G({N}):")
        for i in range(3):
            entries = []
            for j in range(3):
                v = G[i][j]
                if v == 0:
                    entries.append("0")
                else:
                    # Show numerator/denominator digit counts
                    nd = len(str(abs(v.numerator)))
                    dd = len(str(abs(v.denominator)))
                    entries.append(f"[{nd}d/{dd}d]")
            print(f"  {entries}")

    # Check digit growth rate
    max_num_digits = max(len(str(abs(G[i][j].numerator))) for i in range(3) for j in range(3) if G[i][j] != 0)
    max_den_digits = max(len(str(abs(G[i][j].denominator))) for i in range(3) for j in range(3) if G[i][j] != 0)

    if N >= 2:
        gauge_ratios.append((N, max_num_digits, max_den_digits))

    if N <= 10 or N % 5 == 0:
        print(f"  N={N}: max digits = {max_num_digits}(num)/{max_den_digits}(den)")

# Analyze growth
print("\n=== Growth rate of G(N) entries ===")
for i in range(1, len(gauge_ratios)):
    N, nd, dd = gauge_ratios[i]
    N_prev, nd_prev, dd_prev = gauge_ratios[i-1]
    if nd_prev > 0:
        print(f"  N={N}: num {nd} digits (growth {nd-nd_prev}/step), den {dd} digits (growth {dd-dd_prev}/step)")

# NOW try factorial gauge
print("\n\n=== FACTORIAL GAUGE: G̃(N) = gauged_Π(N) · Φ(N)⁻¹ ===")
print("Gauge: multiply Π row by 1/[(-16)^N · (N!)^7]")

Pi = [row[:] for row in I3]
Phi = [row[:] for row in I3]

for N in range(N_MAX):
    if N > 0:
        Pi = mat_mul_F(Pi, M_int(N-1))
        Phi = mat_mul_F(Phi, sym2_companion(N-1))

    # Apply gauge to ALL of Pi
    if N > 0:
        g = F((-16)**N) * F(factorial(N))**7
    else:
        g = F(1)

    gauged_Pi = [[Pi[i][j] / g for j in range(3)] for i in range(3)]

    Phi_inv = mat_inv_F(Phi)
    if Phi_inv is None:
        continue

    G_tilde = mat_mul_F(gauged_Pi, Phi_inv)

    if N <= 8:
        print(f"\nG̃({N}):")
        for i in range(3):
            entries = []
            for j in range(3):
                v = G_tilde[i][j]
                if v == 0:
                    entries.append("0")
                else:
                    nd = len(str(abs(v.numerator)))
                    dd = len(str(abs(v.denominator)))
                    entries.append(f"[{nd}d/{dd}d]")
            print(f"  {entries}")

    max_num_digits = max(len(str(abs(G_tilde[i][j].numerator))) for i in range(3) for j in range(3) if G_tilde[i][j] != 0)
    max_den_digits = max(len(str(abs(G_tilde[i][j].denominator))) for i in range(3) for j in range(3) if G_tilde[i][j] != 0)

    if N <= 10 or N % 5 == 0:
        print(f"  N={N}: max digits = {max_num_digits}(num)/{max_den_digits}(den)")

# Try the CORRECT gauge: extract ALL Pochhammer factors
print("\n\n=== POCHHAMMER GAUGE ===")
print("Using the factored form, gauge = product of Pochhammer factors")

def pochhammer_gauge(N):
    """
    From factored c₀(N): the Pochhammer factors are
    (N+1)(N+2)(N+3)^5(N+4)^3(N+5)^2(2N+3)(2N+5)^2(2N+7)^4(2N+9)^3

    The gauge ratio g(N+1)/g(N) should absorb these.
    For a degree-7 ratio, we need to select 7 factors.

    With P₆(N) = P₆'(N+1), the natural choice is:
    r(N) = -16 * (N+1)(N+3/2)(N+2)(N+5/2)(N+3)(N+7/2)(N+4) / [(N+?)^7]

    But we need to figure out the denominator from c₃'s structure.

    Actually, let's just try the ratio c₀(N)/c₃(N) approach.
    """
    if N == 0:
        return F(1)
    # g(N) = ∏_{k=0}^{N-1} r(k) where r(k) = -c₀(k)/c₃(k+2) or similar
    # For now, just use (-16)^N * (N!)^7
    return F((-16)**N) * F(factorial(N))**7

# Instead, let's try to identify r_i(N) directly.
# Express q(N) = r₀(N)*D_N² + r₁(N)*D_N*E_N + r₂(N)*E_N²
# by using all 9 entries of the CMF product matrix.
print("\n\n=== DIRECT R(N) IDENTIFICATION ===")
print("For each N, solve Π(N) = R(N) · [D², DE, E²; ...] where R(N) is 3×3")

Pi = [row[:] for row in I3]

for N in range(min(N_MAX, 15)):
    if N > 0:
        Pi = mat_mul_F(Pi, M_int(N-1))

    D_N = D_seq[N]
    E_N = E_seq[N]

    # The Sym² solution matrix at step N:
    # Column 1 (from initial [1,0,0]): [D_N², D_N*something, ...]
    # Actually, the fundamental matrix of the Sym² system starting from I
    # is just the product of Sym² companion matrices.
    # But I want to express CMF entries in terms of D_N², D_N*E_N, E_N².

    # For each column j of Π(N), solve:
    # Π(N)[:,j] = a·D_N² + b·D_N·E_N + c·E_N²
    # But this is underdetermined (3 unknowns, 3 equations from 3 rows)
    # Wait, that's exactly determined!

    # Set up: [D_N², D_N*E_N, E_N²]  <- this is just a vector, not a matrix
    # We need THREE linearly independent Sym² vectors to invert.

    # The Sym² fundamental matrix Φ(N) gives us this:
    # Φ(N)[:,0] = first solution starting from [1,0,0]
    # Φ(N)[:,1] = second solution starting from [0,1,0]
    # etc.

    # So Φ(N) is already the matrix of {D², DE, E²} in some basis.
    # And G(N) = Π(N) · Φ(N)⁻¹ is what we computed above.

    # Let me instead look at specific RATIOS of G entries.
    # If G(N) = diag(f₁(N), f₂(N), f₃(N)) · C for some constant C,
    # then f_i(N) should be identifiable.

    pass

# Let me try a completely different approach.
# Instead of the matrix gauge, try to find rational r₀, r₁, r₂ such that:
# q[0][N] = r₀(N) * D_N² + r₁(N) * D_N * E_N + r₂(N) * E_N²
# using the FACT that q satisfies L₂₅ and {D², DE, E²} satisfies L_Sym².
# This means: L₂₅(r₀ D² + r₁ DE + r₂ E²) = 0
# which gives a COUPLED system for r₀, r₁, r₂.

# But a simpler approach: just compute r₀, r₁, r₂ at EACH N by using
# three independent initial conditions for the CMF.

print("\nUsing three independent CMF starting vectors:")
print("Column 0, 1, 2 of CMF product → three independent solutions\n")

Pi = [row[:] for row in I3]
Phi = [row[:] for row in I3]

prev_G = None
for N in range(min(N_MAX, 20)):
    if N > 0:
        Pi = mat_mul_F(Pi, M_int(N-1))
        Phi = mat_mul_F(Phi, sym2_companion(N-1))

    Phi_inv = mat_inv_F(Phi)
    if Phi_inv is None:
        continue

    G = mat_mul_F(Pi, Phi_inv)

    # Check if G(N)/G(N-1) has a simple structure
    if prev_G is not None and N <= 12:
        G_inv_prev = mat_inv_F(prev_G)
        if G_inv_prev is not None:
            ratio = mat_mul_F(G, G_inv_prev)
            # If ratio is diagonal, that's great
            diag = [ratio[i][i] for i in range(3)]
            off_diag = [ratio[i][j] for i in range(3) for j in range(3) if i != j]
            is_diag = all(x == 0 for x in off_diag)
            if is_diag:
                print(f"  N={N}: G(N)/G(N-1) is DIAGONAL: {[str(d) for d in diag]}")
            else:
                print(f"  N={N}: G(N)/G(N-1) is NOT diagonal")
                # Show the matrix
                for i in range(3):
                    row_str = [str(ratio[i][j]) if abs(ratio[i][j]) < 10**10 else f"[{len(str(abs(ratio[i][j].numerator)))}d]" for j in range(3)]
                    print(f"    {row_str}")

    prev_G = G

print("\n\nDone.")
