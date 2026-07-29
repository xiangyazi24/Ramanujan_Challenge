"""
Compute Zudilin companion sequences p^{(2)}_n, p^{(3)}_n using
Miller backward iteration on the Zudilin recurrence.

Strategy: Since p^{(j)}_n/b_n → ζ(j), the sequence p^{(j)}_n has
the same dominant mode as b_n. By seeding backward iteration with
b_N·ζ(j) at large N, we recover p^{(j)} in the dominant direction.
"""
import mpmath
mpmath.mp.dps = 300  # 300 digits

from mpmath import mpf, mp, zeta, pi, fac, binomial
from fractions import Fraction as Q
from math import comb

# Zudilin recurrence coefficients (order 3):
# Z3(n)*b_{n+1} + Z2(n)*b_n + Z1(n)*b_{n-1} + Z0(n)*b_{n-2} = 0
# Equivalently: b_{n+1} = -Z2(n)/Z3(n)*b_n - Z1(n)/Z3(n)*b_{n-1} - Z0(n)/Z3(n)*b_{n-2}

def Z3(n):
    return 2*(946*n**2-731*n+153)*(2*n+1)*(n+1)**3

def Z2(n):
    return -2*(104060*n**6+127710*n**5+12788*n**4-34525*n**3-8482*n**2+3298*n+1071)

def Z1(n):
    return 2*n*(3784*n**5-1032*n**4-1925*n**3+853*n**2+328*n-184)

def Z0(n):
    return -(946*n**2+1161*n+368)*n*(n-1)**3

def compute_b(n):
    return sum(comb(n,k)**2 * comb(n+k,n) * comb(n+2*k,n) for k in range(n+1))

# Compute b_n exactly
N_exact = 50
b_exact = [Q(compute_b(n)) for n in range(N_exact)]

# Verify Zudilin recurrence on b_n
for n in range(2, 10):
    res = Z3(n)*int(b_exact[n+1]) + Z2(n)*int(b_exact[n]) + Z1(n)*int(b_exact[n-1]) + Z0(n)*int(b_exact[n-2])
    assert res == 0, f"Zudilin recurrence failed at n={n}"
print("Zudilin recurrence verified on b_n.")

# High-precision constants
z2 = float(zeta(2))
z3 = float(zeta(3))
z2_mp = zeta(2)
z3_mp = zeta(3)

# Forward compute b_n at high precision
b_mp = [mpf(int(b_exact[n])) for n in range(N_exact)]
# Extend b_mp using recurrence
for n in range(N_exact - 1, 200):
    # b_{n+1} = -Z2(n)/Z3(n)*b_n - Z1(n)/Z3(n)*b_{n-1} - Z0(n)/Z3(n)*b_{n-2}
    bn1 = (-mpf(Z2(n))/mpf(Z3(n)) * b_mp[n]
           - mpf(Z1(n))/mpf(Z3(n)) * b_mp[n-1]
           - mpf(Z0(n))/mpf(Z3(n)) * b_mp[n-2])
    b_mp.append(bn1)

N = len(b_mp) - 1

# Strategy: p^{(2)}_n satisfies the same recurrence as b_n.
# For large n: p^{(2)}_n ≈ ζ(2) · b_n (since r^{(2)}_n → 0).
#
# Forward iteration from initial p^{(2)}_0, p^{(2)}_1, p^{(2)}_2 amplifies the
# dominant mode correctly, so:
# If we set p^{(2)}_n = ζ(2)·b_n for n = 0,1,2 and propagate forward,
# we get p^{(2)}_n ≈ ζ(2)·b_n for all n (since linearity and same recurrence).
# But this is trivially ζ(2)·b_n for ALL n, not the correct p^{(2)}_n!
#
# The actual p^{(2)}_n differs from ζ(2)·b_n by a recessive term.
# To find the correct p^{(2)}_n, we need the initial conditions.
#
# Alternative: use the EXPLICIT FORMULA from the hypergeometric summand.
# p^{(2)}_n = Σ_{k=0}^n F(n,k) · σ^{(2)}_k where σ^{(2)}_k involves harmonic numbers.

# Let me try the parameter-derivative approach:
# F(n,k) = C(n,k)^2 C(n+k,n) C(n+2k,n)
#
# Take d/da of (-a)_k at a = n:
# d/da [(-a)_k] = (-a)_k * Σ_{j=0}^{k-1} 1/(j-a)
# At a = n: (-n)_k = (-1)^k n!/(n-k)! and Σ_{j=0}^{k-1} 1/(j-n) = -Σ_{j=1}^{n} 1/(n-j+1)... complicated.
#
# Actually, there's a simpler formula. In the Apéry context:
# For b_n = Σ_k C(n,k)^2 C(n+k,k)^2 (Apéry for ζ(3)):
# p_n = Σ_k C(n,k)^2 C(n+k,k)^2 · [Σ_{m=1}^n 1/m^3 + Σ_{m=1}^k (-1)^{m-1}/(2m^3 C(n,m)C(n+m,m))]
#
# For our ₅F₄, the companion sequences involve:
# p^{(2)}_n = Σ_k F(n,k) · [some combination of H_k^{(2)}]
# p^{(3)}_n = Σ_k F(n,k) · [some combination of H_k^{(3)}]
#
# But the EXACT formulas depend on which parameters we differentiate.

# Let me try: d/d(lower parameter 1) of ₅F₄ at (lower)_1 = 1
# F(n,k) with (1)_k in denominator → d/(1)_k = -(1)_k · ψ(k+1)
# So the derivative w.r.t. lower parameter 1 at value 1 gives:
# Σ_k F(n,k) · [-ψ(k+1) + ψ(1)] = -Σ_k F(n,k) · H_k
# where H_k = harmonic number = ψ(k+1) - ψ(1) = Σ_{j=1}^k 1/j
#
# But this gives ζ(1)-related, not ζ(2)-related.
# For ζ(2): need second derivative, or derivative of upper parameter...

# The clearest approach: differentiate the UPPER parameters -n, -n at a point
# that introduces 1/(k+1)^2 type sums.
#
# Parameter derivative of (-n)_k w.r.t. n at integer n:
# (-n)_k = (-n)(-n+1)...(-n+k-1)
# d/dn [(-n)_k] = (-n)_k · [-1/(-n) + (-1)/(-n+1) + ... + (-1)/(-n+k-1)]
#               = (-n)_k · Σ_{j=0}^{k-1} (-1)/(-n+j)
#               = (-n)_k · Σ_{j=0}^{k-1} 1/(n-j)
#               = (-n)_k · [1/n + 1/(n-1) + ... + 1/(n-k+1)]
# For k ≤ n, this is = (-n)_k · [H_n - H_{n-k}]
#
# So d/dn of F(n,k) involves [H_n - H_{n-k}] type sums from the two (-n)_k factors,
# and d/dn of (n+1)_k and ((n+1)/2)_k and ((n+2)/2)_k from the other factors.

# This is getting complex. Let me just COMPUTE the companions numerically
# using the RESIDUE approach.

# R_n(t) from Zudilin's construction: let me compute using the Barnes integral form.
# Instead, let me use the fact that the error r^{(j)}_n is small and RATIONAL.

# r^{(2)}_n = b_n·ζ(2) - p^{(2)}_n → 0
# r^{(3)}_n = b_n·ζ(3) - p^{(3)}_n → 0
#
# Since r^{(j)}_n → 0 exponentially, for large n:
# p^{(j)}_n ≈ b_n·ζ(j)
# with the error being exponentially small.
#
# The KEY: r^{(j)}_n = b_n·ζ(j) - p^{(j)}_n is a SPECIFIC recessive solution.
# In the 3D solution space, the dominant mode is b_n, and there are 2 recessive modes.
# r^{(j)}_n is a specific linear combination of the 2 recessive modes.
#
# To compute r^{(j)}_n exactly, I can use MILLER BACKWARD ITERATION for the
# Zudilin recurrence, starting from large N with seed values, and normalizing.

# Miller backward for the Zudilin recurrence:
# The Zudilin recurrence for forward: Z3(n)*y_{n+1} + Z2(n)*y_n + Z1(n)*y_{n-1} + Z0(n)*y_{n-2} = 0
# Rewrite for backward: y_{n-2} = -(Z3(n)*y_{n+1} + Z2(n)*y_n + Z1(n)*y_{n-1}) / Z0(n)
#
# But Z0(n) = -(946n²+1161n+368)·n·(n-1)³ has zeros at n=0 and n=1.
# So backward iteration works for n ≥ 2.
#
# Starting from large N, seed with arbitrary values for y_N, y_{N-1}, y_{N-2},
# and iterate backward. The result converges to the recessive direction.

# Actually, for the Zudilin recurrence, backward iteration suppresses the dominant
# mode (ν₀^n) and amplifies the recessive modes. But there are TWO recessive modes
# (ν₁^n and ν₂^n with |ν₁| ≈ |ν₂|). So backward iteration gives a SPECIFIC
# linear combination of the two recessive modes, depending on the seed.

# What I actually want: the ADJOINT slow solution w^{(0)} such that the bracket
# J(w^{(0)}, y) = c₀(y) for any solution y.
#
# For the Zudilin recurrence, the adjoint recurrence is obtained by applying
# the formal adjoint.

# Let me try a simpler approach: compute p^{(2)}_n + p^{(3)}_n directly.
#
# E_n^Z = (p^{(2)}_n + p^{(3)}_n) - (ζ(2)+ζ(3))·b_n
# This is a specific recessive solution of the Zudilin recurrence.
#
# I can compute E_n^Z via Miller backward iteration for the Zudilin recurrence.
# Then p^{(2)}_n + p^{(3)}_n = (ζ(2)+ζ(3))·b_n + E_n^Z.
#
# Miller backward for E_n^Z:
# Seed with (1, 0, 0) at n=N, N-1, N-2, iterate backward, then normalize.
# Normalization: for n=0, the "true" E_0^Z = (p^{(2)}_0 + p^{(3)}_0) - (ζ(2)+ζ(3))·b_0
#               = (p^{(2)}_0 + p^{(3)}_0) - (ζ(2)+ζ(3))
# But I don't know p^{(2)}_0 + p^{(3)}_0!
#
# Hmm, circular. Let me think again...

# Actually, r^{(j)}_n are defined by specific integral formulas.
# For the ₅F₄ series, there should be explicit formulas like:
# r^{(2)}_n = Σ_{k=0}^n F(n,k)·(-1)^k·... / something
#
# Let me try computing using PARTIAL FRACTIONS at the singularities.

# Different approach: just compute the parameter derivatives of F(n,k) directly.
# d/d(lower_1) ₅F₄ at lower_1 = 1:
# The lower Pochhammer (1)_k contributes -ψ(k+1) + ψ(1) = -H_k per term.
# So: d/d(c₁) [₅F₄(...; c₁, 1, 1, 1/2; 1)]|_{c₁=1}
#   = -Σ_k F(n,k) · H_k
# where H_k = 1 + 1/2 + ... + 1/k.

# Similarly, d²/dc₁² gives Σ F(n,k) · (H_k² + H_k^{(2)}) where H_k^{(2)} = Σ 1/j².

# For ζ(2), we typically need Σ F(n,k) · H_k^{(2)} type sums.
# For ζ(3), we need Σ F(n,k) · H_k^{(3)} type sums.

# Let me just compute these directly for small n.

print("\n=== Computing parameter-derivative companions ===")

def F_nk(n, k):
    return comb(n,k)**2 * comb(n+k,n) * comb(n+2*k,n)

def H(k, s):
    """Generalized harmonic number H_k^{(s)} = Σ_{j=1}^k 1/j^s"""
    return sum(Q(1, j**s) for j in range(1, k+1))

# Companion 1: c1_n = Σ_k F(n,k) · H_k^{(1)} (first harmonic)
# Companion 2: c2_n = Σ_k F(n,k) · H_k^{(2)} (generalized harmonic for ζ(2))
# Companion 3: c3_n = Σ_k F(n,k) · H_k^{(3)} (generalized harmonic for ζ(3))

# Also try: upper parameter derivatives
# c4_n = Σ_k F(n,k) · [2(H_n - H_{n-k})] from the two (-n)_k upper params
# c5_n = Σ_k F(n,k) · [H_{n+k} - H_n] from (n+1)_k upper param

N_comp = 20
b_vals = [Q(compute_b(n)) for n in range(N_comp)]

for name, hfunc in [
    ("H_k^(1)", lambda k: H(k, 1)),
    ("H_k^(2)", lambda k: H(k, 2)),
    ("H_k^(3)", lambda k: H(k, 3)),
    ("2(H_n-H_{n-k})", lambda k, n=None: Q(2)*(H(n,1) - H(n-k,1))),
    ("H_{n+k}-H_n", lambda k, n=None: H(n+k,1) - H(n,1)),
    ("H_{n+2k}-H_n", lambda k, n=None: H(n+2*k,1) - H(n,1)),
]:
    needs_n = name.find('n') >= 0
    comp = []
    for n in range(N_comp):
        if n == 0 and needs_n:
            # Avoid H(0-k,...) issues
            s = Q(0)
            for k in range(n+1):
                try:
                    if needs_n:
                        hval = hfunc(k, n=n)
                    else:
                        hval = hfunc(k)
                    s += Q(F_nk(n, k)) * hval
                except:
                    pass
            comp.append(s)
            continue
        s = Q(0)
        for k in range(n+1):
            try:
                if needs_n:
                    hval = hfunc(k, n=n)
                else:
                    hval = hfunc(k)
                s += Q(F_nk(n, k)) * hval
            except:
                s += Q(0)
        comp.append(s)

    # Check ratio comp_n / b_n → ?
    print(f"\n  Companion: Σ F(n,k) · {name}")
    for n in range(min(15, N_comp)):
        if b_vals[n] != 0:
            ratio = float(comp[n]) / float(b_vals[n])
            print(f"    n={n}: ratio = {ratio:.15f}")
    # Check if ratio → ζ(2) or ζ(3) or ζ(2)+ζ(3)
    if N_comp > 10 and b_vals[10] != 0:
        r10 = float(comp[10]) / float(b_vals[10])
        print(f"    → Compare: ζ(2)={z2:.15f}, ζ(3)={z3:.15f}, ζ(2)+ζ(3)={z2+z3:.15f}")
