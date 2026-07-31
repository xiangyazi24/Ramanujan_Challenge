#!/usr/bin/env sage
"""
Numerical test: does p_n/q_n → ζ(2)+ζ(3) as n → ∞?

If yes, then c_*(p) = (ζ(2)+ζ(3)) · c_*(q), which implies c₀(e) = 0.

The P2.7 recurrence:
  u_{n+3} = (B(n+2)/A(n+2)) u_{n+2} - (C(n+1)/A(n+1)) u_{n+1} + (D(n)/A(n)) u_n

where A, B, C, D are the polynomial coefficients from proof.tex.
"""
from sage.all import *

# Work at 600-digit precision
PREC = 2000  # bits
RR = RealField(PREC)

def A(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def B(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D(n): return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Compute q_n (denominator sequence) using exact rationals
# Initial conditions for q: q_0, q_1, q_2 from proof.tex
# The inner sum: a_n = sum_{k=0}^{n} C(n,k)^2 C(n+k,n) C(n+2k,n)
# q_n = C(2n,n) * a_n

def inner_sum(n):
    s = QQ(0)
    for k in range(n+1):
        s += binomial(n,k)**2 * binomial(n+k,n) * binomial(n+2*k,n)
    return s

# Compute first few q_n
print("Computing q_n sequence...", flush=True)
N_MAX = 200
q = [QQ(0)] * (N_MAX + 1)
for nn in range(min(N_MAX+1, 50)):
    q[nn] = QQ(binomial(2*nn, nn)) * inner_sum(nn)

# Extend using recurrence
for nn in range(3, N_MAX+1):
    if q[nn] == 0:  # not computed by direct formula
        q[nn] = (B(nn-1)/A(nn-1)) * q[nn-1] - (C(nn-2)/A(nn-2)) * q[nn-2] + (D(nn-3)/A(nn-3)) * q[nn-3]

print(f"q_0 = {q[0]}")
print(f"q_1 = {q[1]}")
print(f"q_2 = {q[2]}")
print(f"q_3 = {q[3]}")

# Verify recurrence on computed terms
print("\nVerifying recurrence on q_n...", flush=True)
for nn in range(3, 49):
    lhs = q[nn]
    rhs = (B(nn-1)/A(nn-1)) * q[nn-1] - (C(nn-2)/A(nn-2)) * q[nn-2] + (D(nn-3)/A(nn-3)) * q[nn-3]
    if lhs != rhs:
        print(f"  MISMATCH at n={nn}")
        break
else:
    print("  Recurrence verified on q_0...q_48")

# Now compute p_n. This requires the "numerator" sequence.
# p_n = sum_{k=0}^{n} C(n,k)^2 C(n+k,n) C(n+2k,n) * [H(n,k)]
# where H(n,k) involves harmonic numbers.
# For the P2.7 problem, p_n is defined so that p_n/q_n → ζ(2)+ζ(3).
#
# Actually, we need to be more careful. Let me read the definition from proof.tex.
# The p_n sequence is NOT directly from a hypergeometric sum with harmonic numbers.
# Instead, we should compute it from the recurrence with specific initial conditions.
#
# From the proof: e_n = p_n - L*q_n satisfies the same recurrence.
# And |e_n|^{1/n} → |μ±| << μ₀.
# So p_n/q_n → L = ζ(2)+ζ(3) if c₀(e) = 0.
#
# To compute p_n, we need its initial conditions from proof.tex.
# Let me check what p_0, p_1, p_2 are.

# Actually, from the standard Apéry-like construction, we should read
# the initial values from the proof file.
# For now, let's use a different approach: compute the adjoint bracket
# to extract c_*(p)/c_*(q) directly.

# Alternative approach: just compute q_n to high index and look at ratios.
# The ratio q_{n+1}/q_n → μ₀ (the dominant root).

print("\n--- Ratio q_{n+1}/q_n convergence ---")
for nn in [10, 20, 50, 100, 150, 199]:
    if q[nn] != 0 and q[nn+1] != 0:
        ratio = RR(q[nn+1]) / RR(q[nn])
        print(f"  n={nn}: q_{{n+1}}/q_n = {ratio:.15f}")

# Compute μ₀ numerically
P_cubic = 4*x**3 - 220*x**2 + 8*x - 1
CC = ComplexField(PREC)
R_poly.<x> = PolynomialRing(QQ)
roots = (4*x**3 - 220*x**2 + 8*x - 1).roots(CC, multiplicities=False)
mu0 = max(roots, key=lambda z: z.real())
print(f"\nμ₀ = {RR(mu0.real()):.15f}")
print(f"|μ±| = {abs(roots[1]):.15f}")

# Compute σ₀
sigma0 = 24*(4*mu0 - 1)/(220*mu0**2 - 16*mu0 + 3)
print(f"σ₀ = {RR(sigma0.real()):.15f}")

# Extract dominant coefficient: c_*(q) = lim q_n / (μ₀^n * n^{-σ₀})
print("\n--- Dominant coefficient extraction ---")
for nn in [50, 100, 150, 199]:
    if q[nn] != 0:
        c_q = RR(q[nn]) / (RR(mu0.real())**nn * RR(nn)**(-RR(sigma0.real())))
        print(f"  n={nn}: c_*(q) ≈ {c_q:.15e}")

# Now we need p_n. Let me try to read initial conditions from proof.tex.
print("\n--- Need p_n initial conditions from proof.tex ---")
print("(Will read the file to find p_0, p_1, p_2)")

print("\nDone with q-sequence analysis.")
