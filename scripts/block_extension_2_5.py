#!/usr/bin/env python3
"""
Problem 2.5: Block-extension certificate.

The twisted operator L₂₅♯ has Poincaré roots {1, 17+12√2, 17-12√2}.
Root 1 = "neutral" direction. The CMF is a nontrivial extension of Sym²(Delannoy)
by this neutral direction. G enters as the accumulated cocycle.

Strategy:
1. Extract the neutral (Poincaré root 1) solution of L₂₅♯
2. Express CMF solutions in terms of Sym² solutions + neutral solution
3. Identify the extension cocycle
4. Show the accumulated cocycle gives G

STEP 1: Find the neutral solution of L₂₅♯.
The recurrence ℓ₀(n)f(n) + ℓ₁(n)f(n+1) + ℓ₂(n)f(n+2) + ℓ₃(n)f(n+3) = 0
has a solution ~O(1) as n→∞. Find it.
"""
from fractions import Fraction as F
from math import sqrt, log
import sys

def pochhammer_int(a_num, a_den, n):
    result = F(1)
    for k in range(n):
        result *= F(a_num + k * a_den, a_den)
    return result

def H(n):
    if n == 0: return F(1)
    neg16_n = F(-16)**n
    poch_2 = pochhammer_int(2, 1, n)
    poch_3 = pochhammer_int(3, 1, n)
    poch_5_2 = pochhammer_int(5, 2, n)
    poch_7_2 = pochhammer_int(7, 2, n)
    return neg16_n * poch_2**2 * poch_3**2 * poch_5_2 * poch_7_2**2

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
    return [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]

def mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

N_MAX = 80

# Compute Delannoy solutions
D = [F(1), F(3)]
E_del = [F(0), F(1)]
for n in range(1, N_MAX + 8):
    D.append((F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1))
    E_del.append((F(3*(2*n+1)) * E_del[n] - F(n) * E_del[n-1]) / F(n+1))

DD = [D[n]**2 for n in range(N_MAX+8)]
DE = [D[n]*E_del[n] for n in range(N_MAX+8)]
EE = [E_del[n]**2 for n in range(N_MAX+8)]

# CMF products and twisted solutions
I3 = [[1,0,0],[0,1,0],[0,0,1]]
prod_mat = [row[:] for row in I3]
q_raw = [[], [], []]
for N in range(N_MAX + 5):
    for k in range(3):
        q_raw[k].append(prod_mat[k][0])
    if N < N_MAX + 4:
        prod_mat = mat_mul(prod_mat, M_int(N))

H_vals = [H(n) for n in range(N_MAX + 5)]
q_tw = [[], [], []]
for k in range(3):
    for n in range(N_MAX + 5):
        q_tw[k].append(F(q_raw[k][n]) / H_vals[n])

# ---- Find the neutral solution ----
# q♯_k(n) are linear combinations of the three Poincaré modes:
# q♯_k(n) = A_k · (17+12√2)^n + B_k · (17-12√2)^n + C_k · 1^n
#
# The dominant mode is (17+12√2)^n ≈ 33.97^n.
# q♯_k(n) / D_n² → A_k / (leading coeff of D_n² in dominant mode)
# because D_n² ~ C · (17+12√2)^n (from D_n ~ C' · (3+2√2)^n).

# Compute connection coefficients: q♯_k(n) / D_n²
print("=== Connection coefficients q♯_k(n) / D(n)² ===")
for k in range(3):
    print(f"\n  q♯_{k}(n) / D(n)²:")
    for n in [10, 20, 30, 40, 50]:
        if n < N_MAX and DD[n] != 0:
            r = q_tw[k][n] / DD[n]
            print(f"    n={n}: {float(r):.15f}")

# The ratio should converge to A_k as n→∞.
# Then q♯_k(n) - A_k · D_n² has the subdominant mode eliminated.
# The next mode is (17-12√2)^n ≈ 0.0294^n (recessive).
# So (q♯_k(n) - A_k · D_n²) / D_n² → 0 exponentially.

# But actually, D_n² is one solution and the "neutral" solution (Poincaré root 1)
# is another. Let me think about what the 3 solutions of L₂₅♯ really are.

# The Sym²(Delannoy) solutions D², DE, E² have growth rates:
# D² ~ C₁ · (17+12√2)^n  (dominant)
# DE ~ C₂ · (17+12√2)^n  (also dominant, since D ~ (3+2√2)^n, E ~ constant × (3+2√2)^n)
# E² ~ C₃ · (17+12√2)^n  (also dominant!)
# Wait, E_n/D_n → log(2)/2, so E_n ~ (log(2)/2) · D_n. Thus all three are dominant.

# Actually: D_n ~ c₁ · α^n + c₂ · β^n where α = 3+2√2, β = 3-2√2.
# E_n = (log(2)/2)·D_n - Q_n(3), where Q_n(3) ~ c · β^n.
# So E_n ~ (log(2)/2)·c₁·α^n + [stuff of order β^n].

# Therefore:
# D² ~ c₁² · α^{2n}  [dominant, growth rate α² = 17+12√2]
# DE ~ (log(2)/2)·c₁² · α^{2n}
# E² ~ (log(2)/2)² · c₁² · α^{2n}

# So all three Sym² solutions grow like α^{2n} = (17+12√2)^n!
# The Sym² recurrence's three Poincaré roots are 17+12√2, 17-12√2, 1.
# The root 1 corresponds to... what?

# Actually, the Sym² solutions are NOT just D², DE, E².
# The Sym²(Delannoy) recurrence has 3 solutions with growth rates:
# s₁ ~ (17+12√2)^n (from D²)
# s₂ ~ 1^n (the NEUTRAL direction)
# s₃ ~ (17-12√2)^n (the recessive direction)

# But D², DE, E² are ALL dominated by (17+12√2)^n.
# The neutral solution (Poincaré root 1) must be a specific linear combination
# of D², DE, E² that cancels the dominant mode.

# Let me find it. If D_n ~ a·α^n + b·β^n, E_n ~ q·a·α^n + ... (where q = log(2)/2)
# then D² ~ a²·α^{2n}, DE ~ a²q·α^{2n}, E² ~ a²q²·α^{2n}.
# So E² - 2q·DE + q²·D² ~ 0 · α^{2n} (dominant cancels).
# But this is (E - qD)² ~ [Q_n(3)]² which decays like β^{2n} = (17-12√2)^n.
# That's the RECESSIVE direction, not neutral.

# For the NEUTRAL direction (root 1), I need a combination that:
# - cancels the dominant (17+12√2)^n mode
# - cancels the recessive (17-12√2)^n mode
# - leaves just the O(1) part

# Hmm, but if D² ~ a₀·ρ^{-n} + a₁ + a₂·ρ^n where ρ = 17-12√2,
# then D², DE, E² are three independent solutions, and the neutral solution
# is a combination that keeps only the a₁ part.

# Actually, more carefully:
# D_n = A·α^n + B·β^n where α = 3+2√2, β = 3-2√2 (= 1/α)
# Then:
# D² = A²α^{2n} + 2AB(αβ)^n + B²β^{2n} = A²α^{2n} + 2AB + B²β^{2n}
#                                             since αβ = 9-8 = 1
# Similarly DE = (Aq·α^n + ...)(A·α^n + B·β^n) = ...
# E² = ...

# KEY: D_n² has the term 2AB · (αβ)^n = 2AB · 1^n = constant!
# So D² itself contains the neutral (Poincaré root 1) component!

# The decomposition into Poincaré modes is:
# D² = A² · (17+12√2)^n + 2AB · 1^n + B² · (17-12√2)^n
# where A and B are the connection coefficients of D_n.

# So the "neutral" solution of the Sym² recurrence is:
# u_neutral(n) = D(n)^2 - (A²/B²)·β^{2n}·something... no this isn't right.

# Let me just compute the "Poincaré mode 1" component of D², DE, E².

# From Legendre: D_n = P_n(3), E_n = (log(2)/2)·P_n(3) - Q_n(3)
# P_n(x) ~ (x+√(x²-1))^{n+1/2} / (2√π·n^{1/2}·(x²-1)^{1/4})  as n→∞
# Q_n(x) ~ (x-√(x²-1))^{n+1/2} · √π / (2·n^{1/2}·(x²-1)^{1/4})  as n→∞

# At x=3: α = 3+2√2, β = 3-2√2
# P_n(3) ~ α^{n+1/2} / (2√π·n^{1/2}·(8)^{1/4}) [with α^{1/2} and 8^{1/4} = 2^{3/4}]
# Q_n(3) ~ β^{n+1/2} · √π / (2·n^{1/2}·8^{1/4})

# So:
# D_n² = P_n(3)² ~ α^{2n+1} / (4π·n·8^{1/2}) = α^{2n+1} / (4π·n·2√2)
# But also has a cross term 2·P_n·(something recessive) and a double-recessive term.

# The cross term 2·C_dom·C_rec · (αβ)^n = 2·C_dom·C_rec · 1^n is the NEUTRAL part.

# Let's compute this numerically.
# D_n = A_dom · α^n + A_rec · β^n  (approximately, ignoring n-dependent factors)
# Then D_n² = A_dom² α^{2n} + 2·A_dom·A_rec·(αβ)^n + A_rec²·β^{2n}
#           = A_dom²·(17+12√2)^n + 2·A_dom·A_rec + A_rec²·(17-12√2)^n

# So the Poincaré-1 component of D² is exactly 2·A_dom·A_rec (up to n-dependent corrections).

# For D·E: E_n ~ q·D_n - Q_n(3) where Q_n ~ B_rec · β^n
# So D·E ~ (q·A_dom·α^n + ...)(A_dom·α^n + A_rec·β^n)
# = q·A_dom²·α^{2n} + (q·A_dom·A_rec + A_dom·(...)·1^n + ...

# This is getting complicated. Let me just compute numerically.

# The idea: at large n, the ratio q♯_k(n)/D_n² converges to the
# "connection coefficient" with the dominant Poincaré mode.
# The RESIDUAL q♯_k(n) - (conn_coeff)·D_n² should behave like O(1)
# (the neutral mode) plus O(ρ^n) (the recessive mode).

# Compute the residuals
print("\n=== Residual q♯_0(n) - c₀·D(n)² ===")
# Estimate c₀ from large n
c0_est = q_tw[0][60] / DD[60]
print(f"c₀ estimate (n=60): {float(c0_est):.15f}")

for n in [20, 30, 40, 50, 60]:
    resid = q_tw[0][n] - c0_est * DD[n]
    # This should be ~ const + O(ρ^n)
    print(f"  n={n}: residual = {float(resid):.8e}, residual/DD[n] = {float(resid/DD[n]):.6e}")

# The residual should converge to a constant as n→∞
print("\n=== Neutral-mode extraction ===")
# q♯_0(n) = c₀·D² + neutral + recessive
# q♯_0(n) - c₀·D² ≈ neutral_coeff + recessive(n)
# Use two large n values to estimate the neutral coefficient:
# At large n, recessive is negligible.

for K in range(3):
    c_est = q_tw[K][70] / DD[70]
    residuals = []
    for n in [30, 40, 50, 60, 70]:
        resid = q_tw[K][n] - c_est * DD[n]
        residuals.append((n, float(resid)))

    print(f"\n  q♯_{K}: c_dom = {float(c_est):.12f}")
    for n, r in residuals:
        print(f"    n={n}: residual = {r:.10e}")

    # The residuals should converge to a constant
    if len(residuals) >= 2:
        diff = abs(residuals[-1][1] - residuals[-2][1])
        print(f"    Last two residuals differ by: {diff:.6e}")

# Now the key question: what linear combination of q♯_0, q♯_1, q♯_2
# has the dominant mode EXACTLY zero, leaving only neutral + recessive?
# From the CMF matrix: q♯_k(n) = P_n[0][k] / H_n
# The "sum" s_n = Σ_k P_n[0][k] = sum of row 0 entries of the product matrix
# Divided by H_n.

# Actually, what specific combination gives the CMF limit?
# Let me check: what does the PROBLEM say the CMF converges to?
# The CMF is typically set up so that the ratio of specific entries → G.

# Check: row 0 of M(0)·M(1)·...·M(n-1)
print("\n=== CMF row 0 entries ===")
prod_mat2 = [row[:] for row in I3]
for n in range(51):
    prod_mat2 = mat_mul(prod_mat2, M_int(n))

r0 = prod_mat2[0]
print(f"P(51)[0] = [{r0[0]}, {r0[1]}, {r0[2]}]")
print(f"Ratios:")
if r0[0] != 0:
    print(f"  P₀₁/P₀₀ = {float(F(r0[1])/F(r0[0])):.15f}")
    print(f"  P₀₂/P₀₀ = {float(F(r0[2])/F(r0[0])):.15f}")
    # Also check other combinations
    print(f"  (P₀₁+P₀₂)/P₀₀ = {float(F(r0[1]+r0[2])/F(r0[0])):.15f}")

print(f"\n  Catalan G = {0.915965594177219015054603514932384110774:.15f}")
print(f"  ln(2)/2  = {log(2)/2:.15f}")
print(f"  2G       = {2*0.915965594177219015054603514932384110774:.15f}")

# Check if P₀₁/P₀₀ → G (or some other combination)
# From the problem: the CMF approaches G. Let's check all plausible ratios.

# Also check the second and third rows
for row_idx in range(3):
    r = prod_mat2[row_idx]
    print(f"\nRow {row_idx}: [{float(F(r[0])):.6e}, {float(F(r[1])):.6e}, {float(F(r[2])):.6e}]")
    if r[0] != 0:
        for j in range(1, 3):
            ratio = float(F(r[j]) / F(r[0]))
            print(f"  P[{row_idx}][{j}]/P[{row_idx}][0] = {ratio:.15f}")

print("\nDone.")
