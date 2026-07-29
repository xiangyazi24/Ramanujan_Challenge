#!/usr/bin/env python3
"""
Test harmonic-enriched summands for P2.7.

Compute:
  a_n = Σ_k C(n,k)² C(n+k,n) C(n+2k,n)           (AESZ #209)
  S₂(n) = Σ_k C(n,k)² C(n+k,n) C(n+2k,n) · H_k^{(2)}
  S₃(n) = Σ_k C(n,k)² C(n+k,n) C(n+2k,n) · H_k^{(3)}
  D(n) = Σ_k T₀(n,k) · [digamma derivatives]

Check if α·S₂(n) + β·S₃(n) + γ·a_n satisfies P2.7.
"""
import mpmath
mpmath.mp.dps = 100

def binom(n, k):
    if k < 0 or k > n: return 0
    r = 1
    for i in range(k): r = r * (n - i) // (i + 1)
    return r

def harmonic(k, j):
    """H_k^{(j)} = Σ_{m=1}^k 1/m^j"""
    return sum(mpmath.mpf(1)/mpmath.mpf(m)**j for m in range(1, k+1))

N = 25

# AESZ #209 and harmonic-weighted sums
print("Computing AESZ #209 and harmonic sums...")
a = [mpmath.mpf(0)] * N
S2 = [mpmath.mpf(0)] * N
S3 = [mpmath.mpf(0)] * N

for n in range(N):
    for k in range(n+1):
        T = mpmath.mpf(binom(n,k)**2 * binom(n+k,n) * binom(n+2*k,n))
        a[n] += T
        S2[n] += T * harmonic(k, 2)
        S3[n] += T * harmonic(k, 3)
    if n % 5 == 0:
        print(f"  n={n}: a={float(a[n]):.6e}, S2={float(S2[n]):.6e}, S3={float(S3[n]):.6e}")

# P2.7 recurrence
def A_p(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def B_p(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_p(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_p(n): return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

q = [mpmath.mpf(0)] * N
q[0] = mpmath.mpf('-215040420000')
q[1] = mpmath.mpf('-167282265043404') / mpmath.mpf('905')
q[2] = mpmath.mpf('-964185327658080') / mpmath.mpf('6071')
for n in range(2, N-1):
    q[n+1] = (mpmath.mpf(B_p(n))/A_p(n)*q[n]
              - mpmath.mpf(C_p(n-1))/A_p(n-1)*q[n-1]
              + mpmath.mpf(D_p(n-2))/A_p(n-2)*q[n-2])

# Check if α S₂ + β S₃ + γ a satisfies P2.7
# For this, we need to find α,β,γ such that the sequence u_n = α S₂(n) + β S₃(n) + γ a(n)
# satisfies u_{n+1} = (B/A) u_n - (C/A) u_{n-1} + (D/A) u_{n-2}
# This is linear in (α,β,γ), so for each n we get a constraint.

print("\n" + "="*70)
print("Checking if linear combo α·S₂ + β·S₃ + γ·a satisfies P2.7")
print("="*70)

# Define residual: R(u, n) = u_{n+1} - (B/A)u_n + (C/A)u_{n-1} - (D/A)u_{n-2}
# For u = α S₂ + β S₃ + γ a:
# R = α R(S₂,n) + β R(S₃,n) + γ R(a,n)
# We want R=0 for all n, so R(S₂,n), R(S₃,n), R(a,n) must be proportional.

def recurrence_residual(seq, n):
    """R(u,n) = u_{n+1} - (B/A)u_n + (C/A)u_{n-1} - (D/A)u_{n-2}"""
    return (seq[n+1]
            - mpmath.mpf(B_p(n))/A_p(n) * seq[n]
            + mpmath.mpf(C_p(n-1))/A_p(n-1) * seq[n-1]
            - mpmath.mpf(D_p(n-2))/A_p(n-2) * seq[n-2])

print("\nP2.7 residuals for individual sequences:")
print(f"{'n':>4s} | {'R(a,n)':>14s} | {'R(S₂,n)':>14s} | {'R(S₃,n)':>14s} | {'R(S₂)/R(a)':>14s} | {'R(S₃)/R(a)':>14s}")
print("-"*90)
for n in range(2, min(20, N-2)):
    Ra = recurrence_residual(a, n)
    RS2 = recurrence_residual(S2, n)
    RS3 = recurrence_residual(S3, n)
    r2a = float(RS2/Ra) if abs(Ra) > 1e-50 else float('inf')
    r3a = float(RS3/Ra) if abs(Ra) > 1e-50 else float('inf')
    print(f"{n:4d} | {float(Ra):+14.6e} | {float(RS2):+14.6e} | {float(RS3):+14.6e} | {r2a:+14.10f} | {r3a:+14.10f}")

# If R(S₂)/R(a) is constant = -β/γ and R(S₃)/R(a) = -γ'/γ, then we have a solution.
# Let's check constancy more carefully.
print("\n--- Constancy check of R(S₂)/R(a) ---")
ratios_2a = []
ratios_3a = []
for n in range(2, min(18, N-2)):
    Ra = recurrence_residual(a, n)
    RS2 = recurrence_residual(S2, n)
    RS3 = recurrence_residual(S3, n)
    if abs(Ra) > 1e-50:
        ratios_2a.append(float(RS2/Ra))
        ratios_3a.append(float(RS3/Ra))

if len(ratios_2a) > 1:
    print(f"  R(S₂)/R(a) range: [{min(ratios_2a):.10f}, {max(ratios_2a):.10f}]")
    print(f"  R(S₃)/R(a) range: [{min(ratios_3a):.10f}, {max(ratios_3a):.10f}]")
    print(f"  Variation in R(S₂)/R(a): {max(ratios_2a)-min(ratios_2a):.2e}")
    print(f"  Variation in R(S₃)/R(a): {max(ratios_3a)-min(ratios_3a):.2e}")

# Also try: does q_n itself lie in span(a, S₂, S₃)?
print("\n" + "="*70)
print("Does q_n ∈ span(a, S₂, S₃)?")
print("="*70)
# q_n = α a_n + β S₂(n) + γ S₃(n)
# Use n=0,1,2 to solve for α,β,γ, then check n=3,...
M = mpmath.matrix(3, 3)
rhs = mpmath.matrix(3, 1)
for i in range(3):
    M[i,0] = a[i]
    M[i,1] = S2[i]
    M[i,2] = S3[i]
    rhs[i,0] = q[i]

try:
    sol = mpmath.lu_solve(M, rhs)
    alpha, beta, gamma = sol[0,0], sol[1,0], sol[2,0]
    print(f"  α = {float(alpha):.15e}")
    print(f"  β = {float(beta):.15e}")
    print(f"  γ = {float(gamma):.15e}")

    print("\n  Verification:")
    for n in range(min(15, N)):
        pred = alpha * a[n] + beta * S2[n] + gamma * S3[n]
        diff = abs(pred - q[n])
        rel = float(diff / abs(q[n])) if abs(q[n]) > 1e-50 else float(diff)
        tag = " *** MATCH ***" if rel < 1e-20 else ""
        print(f"    n={n:2d}: rel err = {rel:.6e}{tag}")
except Exception as e:
    print(f"  Solve failed: {e}")

# Also try the digamma derivative sum
print("\n" + "="*70)
print("Digamma derivative sum")
print("="*70)
# D(n) = Σ_k T₀(n,k) · [2(ψ(n+k+1)-ψ(k+1)) + (ψ(n+2k+1)-ψ(2k+1))]
# where ψ = digamma = Γ'/Γ

D = [mpmath.mpf(0)] * N
for n in range(N):
    for k in range(n+1):
        T = mpmath.mpf(binom(n,k)**2 * binom(n+k,n) * binom(n+2*k,n))
        dig_term = (2 * (mpmath.digamma(n+k+1) - mpmath.digamma(k+1))
                    + (mpmath.digamma(n+2*k+1) - mpmath.digamma(2*k+1)))
        D[n] += T * dig_term
    if n % 5 == 0:
        print(f"  n={n}: D={float(D[n]):.6e}")

# Check if q ∈ span(a, S₂, S₃, D)
print("\n  q_n ∈ span(a, S₂, S₃, D)?")
M4 = mpmath.matrix(4, 4)
rhs4 = mpmath.matrix(4, 1)
for i in range(4):
    M4[i,0] = a[i]
    M4[i,1] = S2[i]
    M4[i,2] = S3[i]
    M4[i,3] = D[i]
    rhs4[i,0] = q[i]

try:
    sol4 = mpmath.lu_solve(M4, rhs4)
    print(f"  α={float(sol4[0,0]):.10e}, β={float(sol4[1,0]):.10e}, γ={float(sol4[2,0]):.10e}, δ={float(sol4[3,0]):.10e}")

    for n in range(min(12, N)):
        pred = sol4[0,0]*a[n] + sol4[1,0]*S2[n] + sol4[2,0]*S3[n] + sol4[3,0]*D[n]
        diff = abs(pred - q[n])
        rel = float(diff / abs(q[n])) if abs(q[n]) > 1e-50 else float(diff)
        tag = " *** MATCH ***" if rel < 1e-20 else ""
        print(f"    n={n:2d}: rel err = {rel:.6e}{tag}")
except Exception as e:
    print(f"  Solve failed: {e}")

# Also try: H_n sums (harmonic numbers of n, not k)
print("\n" + "="*70)
print("Harmonic numbers of n (not k)")
print("="*70)
# S₂ₙ(n) = a_n · H_n^{(2)}
# S₃ₙ(n) = a_n · H_n^{(3)}
S2n = [a[n] * harmonic(n, 2) for n in range(N)]
S3n = [a[n] * harmonic(n, 3) for n in range(N)]

M6 = mpmath.matrix(6, 6)
rhs6 = mpmath.matrix(6, 1)
seqs = [a, S2, S3, D, S2n, S3n]
for i in range(6):
    for j in range(6):
        M6[i,j] = seqs[j][i]
    rhs6[i,0] = q[i]

try:
    sol6 = mpmath.lu_solve(M6, rhs6)
    print("  Coefficients:", [float(sol6[j,0]) for j in range(6)])
    print("  Verification:")
    for n in range(min(15, N)):
        pred = sum(sol6[j,0] * seqs[j][n] for j in range(6))
        diff = abs(pred - q[n])
        rel = float(diff / abs(q[n])) if abs(q[n]) > 1e-50 else float(diff)
        tag = " *** MATCH ***" if rel < 1e-15 else ""
        print(f"    n={n:2d}: rel err = {rel:.6e}{tag}")
except Exception as e:
    print(f"  Solve failed: {e}")
