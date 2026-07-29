#!/usr/bin/env python3
"""P2.5: Test the Jacobi parameter-derivative kernel from Q4873.

The idea: Jacobi polynomials P_n^{(0,ε)}(1-2X) at ε=0 are Legendre.
Differentiating in ε produces Catalan's constant G.
The derivative kernel R_n^J has an exact Catalan-moment identity.
Test whether this matches the CMF error ê_n = G·Q̂_n - P̂_n.
"""
import mpmath as mp
mp.mp.dps = 150

G = mp.catalan

# === Part 1: CMF sequences ===

def pochhammer(a, n):
    r = mp.mpf(1)
    for k in range(n):
        r *= (a + k)
    return r

def H_n(n):
    return (mp.mpf(-16))**n * pochhammer(2, n)**2 * pochhammer(3, n)**2 * \
           pochhammer(mp.mpf(5)/2, n) * pochhammer(mp.mpf(7)/2, n)**2

def M_entries(n):
    n = mp.mpf(n)
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

A_init = [[mp.mpf(30921), mp.mpf(-32972), mp.mpf(8240)],
          [mp.mpf(33750), mp.mpf(-36000), mp.mpf(9000)]]

# Compute Q̂_n, P̂_n (H_n-normalized)
# v_n = (1, P_n, Q_n) row vector, A·M(0)·M(1)·...·M(n-1) gives v_n
# Actually: v_0 = A[0], v_1 = A[1], then v_{n+1} = M(n)·v_n... no.
# The CMF structure: [1, P_n, Q_n] = [1, P_0, Q_0] · M(0)·M(1)·...·M(n-1)
# Wait, let me re-derive from the recurrence structure.
#
# The 3×3 CMF gives a recurrence:
# c₃(n) u_{n+1} + c₂(n) u_n + c₁(n) u_{n-1} + c₀(n) u_{n-2} = 0
# for both P and Q sequences.
#
# Let me just use the shift recurrence directly.

# From the recurrence coefficients (degree 13 in n):
c_coeffs = [
    [-170972650800, -826494925500, -1792449886332, -2317972607944, -2000297648936,
     -1219354055500, -541255279788, -177419351856, -43002662976, -7620091136,
     -960400960, -81589760, -4190208, -98304],
    [8781630505200, 38850314624124, 78557994908508, 96136040496551, 79442239242197,
     46814452218572, 20241514501104, 6502490145168, 1552168938336, 271943188864,
     33995217088, 2871763456, 146952192, 3440640],
    [-21132458248680, -87529225645944, -165451256319618, -189073879129764, -145809619841418,
     -80164318460172, -32338316008004, -9694892892592, -2160716677664, -353683596544,
     -41340724928, -3268370944, -156684288, -3440640],
    [587448626688, 2442715444224, 4635428285664, 5317694979920, 4116150568664,
     2270943978716, 919036676572, 276298241680, 61721801728, 10120470656,
     1184128064, 93632000, 4485120, 98304],
]

def c_poly(j, n):
    val = mp.mpf(0)
    nn = mp.mpf(n)
    for k, coeff in enumerate(c_coeffs[j]):
        val += mp.mpf(coeff) * nn**k
    return val

# Compute from initial values.
# We need P_n, Q_n for n = 0, 1, 2, ...
# The recurrence is: c₀(n) u_{n-2} + c₁(n) u_{n-1} + c₂(n) u_n + c₃(n) u_{n+1} = 0
# So u_{n+1} = -(c₂(n) u_n + c₁(n) u_{n-1} + c₀(n) u_{n-2}) / c₃(n)

# Initial: from the two initial rows A and the matrix M(0)
# v_0 = A[0] = [30921, -32972, 8240]
# v_1 = A[1] = [33750, -36000, 9000]
# The columns are [1, P, Q]: so P_0 = -32972, Q_0 = 8240
# But actually the normalization might be different.
# Let me think... The CMF structure:
# The 3×3 matrix M(n) acts on the vector (1, P_n, Q_n).
# So if A[0] = (a₀₀, a₀₁, a₀₂) gives (1, P_0, Q_0):
# P_0/1 = a₀₁/a₀₀ = -32972/30921
# Q_0/1 = a₀₂/a₀₀ = 8240/30921
# Similarly for row 1:
# P_1/1 = -36000/33750 = -16/15
# Q_1/1 = 9000/33750 = 4/15

# Actually, in the CMF framework, the rows of A give the
# period vectors at n=0,1. Let me interpret:
# Row 0: (1, P_0, Q_0) means the "constant" component is 1
# But actually these might be: (Q_0^{(0)}, Q_0^{(1)}, Q_0^{(2)})
# three different sequences.

# Let me just use the shift recurrence to compute sequences.
# I need P_n and Q_n. Let me define them from the initial row data.

# From the CMF, the vector v_n = (1, P_n, Q_n) satisfies v_n = v_0 · M(0)·...·M(n-1)
# Check: v_0 = A[0], v_1 = A[0] · M(0) should equal A[1].

v0 = A_init[0][:]
M0 = M_entries(0)
# v1_computed = v0 · M(0)
v1_comp = [sum(v0[j] * M0[j][i] for j in range(3)) for i in range(3)]
print("v0 · M(0) =", [mp.nstr(x, 15) for x in v1_comp])
print("A[1]      =", [mp.nstr(x, 15) for x in A_init[1]])

# Check if they match (up to scaling)
for i in range(3):
    if v1_comp[i] != 0:
        ratio = A_init[1][i] / v1_comp[i]
        print(f"  ratio[{i}] = {mp.nstr(ratio, 15)}")

# Let me just extract P and Q from the initial rows and the recurrence.
# If the first column is 1, then:
# Seq a: 30921, 33750, ...
# Seq b: -32972, -36000, ...
# Seq c: 8240, 9000, ...
# All three satisfy the same recurrence.
# Q̂_n = a_n / H_n, P̂_n = b_n / H_n (or some assignment of columns)

# But for checking the gap ê_n = G·Q̂_n - P̂_n, I need to know which
# columns correspond to the approximation p_n/q_n → G.

# Let me compute the sequences from the recurrence and check which
# ratio → G.

# Three sequences from the three columns of A:
seqs = [[A_init[0][j], A_init[1][j]] for j in range(3)]

# Need one more initial value. From v_2 = v_1 · M(1):
M1 = M_entries(1)
v1 = A_init[1][:]
v2 = [sum(v1[j] * M1[j][i] for j in range(3)) for i in range(3)]
for j in range(3):
    seqs[j].append(v2[j])

print("\nInitial values:")
for j in range(3):
    print(f"  seq[{j}]: {[mp.nstr(x, 15) for x in seqs[j]]}")

# Now extend using the shift recurrence
NMAX = 20
for n in range(2, NMAX):
    # v_{n+1} = v_n · M(n)
    vn = [seqs[j][-1] for j in range(3)]
    Mn = M_entries(n)
    vn1 = [sum(vn[j] * Mn[j][i] for j in range(3)) for i in range(3)]
    for j in range(3):
        seqs[j].append(vn1[j])

# Check which pair of sequences gives ratio → G
print("\nRatios seq[j]/seq[k] for large n:")
for j in range(3):
    for k in range(3):
        if j != k and seqs[k][-1] != 0:
            ratio = seqs[j][-1] / seqs[k][-1]
            print(f"  seq[{j}]/seq[{k}] at n={NMAX}: {mp.nstr(ratio, 20)}", end="")
            if abs(ratio - G) < mp.mpf(10)**(-5):
                print(" ← close to G!", end="")
            print()

# Identify: which ratio → G?
print(f"\nG = {mp.nstr(G, 20)}")

# === Part 2: Jacobi derivative kernel ===
# P_n^{(0,ε)}(z) satisfies the 3-term Jacobi recurrence
# At ε=0: P_n^{(0,0)}(z) = P_n(z) (Legendre)
# P_0(z) = 1, P_1(z) = z
# (n+1)P_{n+1}(z) = (2n+1)z P_n(z) - n P_{n-1}(z)

# ∂_ε P_n^{(0,ε)}(z)|_{ε=0} = ?
# From the Jacobi recurrence:
# 2(n+1)(n+ε+1)(2n+ε+1) P_{n+1} = [(2n+ε+1)(2n+ε+2)(2n+ε) z - ε² ...]  P_n - ...
# This is complicated. Let me use a different approach.

# The Jacobi polynomial has the hypergeometric representation:
# P_n^{(α,β)}(z) = C(n,α,β) · ₂F₁(-n, n+α+β+1; α+1; (1-z)/2)
# where C = (α+1)_n / n!

# At α=0, β=ε:
# P_n^{(0,ε)}(z) = ₂F₁(-n, n+ε+1; 1; (1-z)/2)

# So ∂_ε P_n^{(0,ε)}(z)|_{ε=0} = ∂_ε ₂F₁(-n, n+ε+1; 1; (1-z)/2)|_{ε=0}

# ₂F₁(-n, b; 1; x) = Σ_{k=0}^n (-n)_k (b)_k / (1)_k / k! · x^k

# ∂/∂b [₂F₁(-n, b; 1; x)] = Σ_{k=0}^n (-n)_k (b)_k / k!² · x^k · [ψ(b+k) - ψ(b)]
# where ψ is the digamma function.

# At b = n+1:
# ∂_ε|₀ = Σ_{k=0}^n (-n)_k (n+1)_k / k!² · x^k · [ψ(n+1+k) - ψ(n+1)]

# But x = (1-z)/2. At z = 1-2X, x = X. At X = -1: x = -1. At z=3: x = -1.

# Let me compute these at z = 3 (i.e. x = (1-3)/2 = -1):

print("\n=== Jacobi derivative kernel ===")

# Legendre at z=3: P_n(3) = ₂F₁(-n, n+1; 1; -1) = Σ_{k=0}^n (-n)_k(n+1)_k/k!² · (-1)^k
def legendre_at_3(n):
    val = mp.mpf(0)
    term = mp.mpf(1)
    for k in range(n + 1):
        val += term
        if k < n:
            term *= mp.mpf(-n + k) * mp.mpf(n + 1 + k) / mp.mpf(k + 1)**2 * mp.mpf(-1)
    return val

# ∂_ε P_n^{(0,ε)}(3)|_{ε=0}
def jacobi_deriv_at_3(n):
    val = mp.mpf(0)
    for k in range(n + 1):
        # (-n)_k (n+1)_k / k!² · (-1)^k · [ψ(n+1+k) - ψ(n+1)]
        poch_neg_n = mp.mpf(1)
        for j in range(k):
            poch_neg_n *= (-n + j)
        poch_n1 = mp.mpf(1)
        for j in range(k):
            poch_n1 *= (n + 1 + j)
        fac_k_sq = mp.factorial(k)**2
        psi_diff = mp.digamma(n + 1 + k) - mp.digamma(n + 1)
        val += poch_neg_n * poch_n1 / fac_k_sq * mp.mpf(-1)**k * psi_diff
    return val

print("\nLegendre P_n(3) and ∂_ε P_n^{(0,ε)}(3)|_{ε=0}:")
for n in range(8):
    Pn = legendre_at_3(n)
    dPn = jacobi_deriv_at_3(n)
    print(f"  n={n}: P_n(3) = {mp.nstr(Pn, 15)}, ∂_ε = {mp.nstr(dPn, 15)}")

# Delannoy numbers D_n = P_n(3) (Legendre polynomials at z=3)
# So D_n = Σ_{k=0}^n C(n,k)² 2^k (central Delannoy)
# This is well-known: D_n = P_n(3)

# The Catalan integral: G = ∫₀¹ (-log t)/(1+t²) dt
# The moment functional: C[t^k] = (-1)^k [G - Σ_{j<k} (-1)^j/(2j+1)²]
# So C[P_n(1-2t)] = ?

# Actually, the relevant integral for the Jacobi kernel involves
# ∫₀¹ (-log t)/(1+t²) R_n^J(t²) dt
# where R_n^J is derived from the ε-derivative.

# From Q4873: R_n^J(X) is the polynomial obtained from the derivative
# of the Padé error. The exact formula involves:
# R_n^J(X) comes from ∂_ε [Q_n(ε)·P_n^{(0,ε)}(1-2X) - P_n(ε)] |_{ε=0}
# where Q_n(ε), P_n(ε) are the Padé approximants.

# Let me take a simpler approach: directly compute the CMF error
# and the Catalan moment of Legendre-type polynomials.

# === Part 3: Compute CMF error and test candidates ===

# Identify which sequences are P and Q
# seq[0] starts: 30921, 33750, ...
# seq[1] starts: -32972, -36000, ...
# seq[2] starts: 8240, 9000, ...

# Compute H_n-normalized versions
print("\n=== H_n-normalized sequences ===")
for j in range(3):
    print(f"seq[{j}]_hat (first 5):", end=" ")
    for n in range(5):
        Hn = H_n(n)
        print(f"{mp.nstr(seqs[j][n]/Hn, 10)}", end=" ")
    print()

# Check: for which (j,k), seq[j]/seq[k] → G as n → ∞?
# G ≈ 0.9159655941772190...
# From the ratios above, find the pair.

# Let me just compute ê_n = G · seq[k] - seq[j] for each candidate pair
# and check if it decays like n^{-3} (neutral mode)
print("\n=== Testing error sequences ===")
for j in range(3):
    for k in range(3):
        if j == k:
            continue
        # e_n = G * seq[k][n] - seq[j][n]
        errors = [G * seqs[k][n] - seqs[j][n] for n in range(NMAX + 1)]
        # Check decay
        if errors[-1] != 0 and errors[-2] != 0:
            ratio = errors[-1] / errors[-2]
            # Check if errors are O(ρ^n) or O(n^{-p})
            if abs(ratio) < 0.5:
                print(f"  G*seq[{k}] - seq[{j}]: ratio = {mp.nstr(ratio, 10)} (exponential decay)")
            elif abs(ratio) > 0.9 and abs(ratio) < 1.1:
                # Polynomial decay: check power
                r = mp.log(abs(errors[-1]) / abs(errors[-5])) / mp.log(mp.mpf(NMAX) / mp.mpf(NMAX - 4))
                print(f"  G*seq[{k}] - seq[{j}]: power law ~ n^{mp.nstr(r, 6)}")
            else:
                print(f"  G*seq[{k}] - seq[{j}]: ratio = {mp.nstr(ratio, 10)}")

# === Part 4: Catalan moment of Legendre polynomials ===
# C[f] = ∫₀¹ (-log t)/(1+t²) f(t²) dt
# For f(x) = x^k: C[x^k] = ∫₀¹ (-log t)/(1+t²) t^{2k} dt
#   = ∫₀¹ t^{2k} (-log t) Σ_{j≥0} (-t²)^j dt
#   = Σ_{j≥0} (-1)^j ∫₀¹ t^{2k+2j} (-log t) dt
#   = Σ_{j≥0} (-1)^j / (2k+2j+1)²

# C[x^k] = Σ_{j≥0} (-1)^j / (2k+2j+1)²
# This is related to Catalan's constant via partial sums.

# Let me compute C[P_n(1-2x)] numerically (quadrature)
def catalan_moment(poly_coeffs):
    """Compute ∫₀¹ (-log t)/(1+t²) · poly(t²) dt where poly has given coefficients."""
    def integrand(t):
        x = t**2
        poly_val = sum(c * x**k for k, c in enumerate(poly_coeffs))
        return (-mp.log(t)) / (1 + t**2) * poly_val
    return mp.quad(integrand, [0, 1])

# Legendre P_n(1-2X) as polynomial in X
def legendre_poly_shifted(n):
    """Return coefficients [c_0, ..., c_n] of P_n(1-2X) as polynomial in X."""
    # P_n(z) where z = 1-2X → polynomial in X
    # Use recurrence: P_0 = 1, P_1 = 1-2X
    # (n+1)P_{n+1} = (2n+1)(1-2X)P_n - n P_{n-1}
    if n == 0:
        return [mp.mpf(1)]
    if n == 1:
        return [mp.mpf(1), mp.mpf(-2)]

    prev2 = [mp.mpf(1)]  # P_0
    prev1 = [mp.mpf(1), mp.mpf(-2)]  # P_1

    for m in range(1, n):
        # P_{m+1} = ((2m+1)(1-2X) P_m - m P_{m-1}) / (m+1)
        # (1-2X) * P_m: multiply poly by (1-2X)
        shifted = [mp.mpf(0)] * (len(prev1) + 1)
        for k, c in enumerate(prev1):
            shifted[k] += c  # * 1
            shifted[k+1] += -2 * c  # * (-2X)

        new = [mp.mpf(0)] * (len(shifted))
        for k in range(len(shifted)):
            new[k] += mp.mpf(2*m+1) * shifted[k] / mp.mpf(m+1)
        for k in range(len(prev2)):
            new[k] -= mp.mpf(m) * prev2[k] / mp.mpf(m+1)

        prev2 = prev1
        prev1 = new

    return prev1

# Compute Catalan moments of P_n(1-2X) for X = t²
print("\n=== Catalan moments of P_n(1-2t²) ===")
cat_moments = []
for n in range(12):
    coeffs = legendre_poly_shifted(n)
    # We need ∫₀¹ (-log t)/(1+t²) P_n(1-2t²) dt
    cm = catalan_moment(coeffs)
    cat_moments.append(cm)
    Pn3 = legendre_at_3(n)
    print(f"  n={n}: C[P_n] = {mp.nstr(cm, 15)}, P_n(3) = {mp.nstr(Pn3, 15)}, ratio C/P = {mp.nstr(cm/Pn3 if Pn3 != 0 else 0, 15)}")

# Check if C[P_n] = G·P_n(3) - something rational
print("\n=== Checking C[P_n(1-2t²)] = G·D_n - r_n ===")
for n in range(12):
    Dn = legendre_at_3(n)
    r_n = G * Dn - cat_moments[n]
    print(f"  n={n}: G·D_n - C[P_n] = {mp.nstr(r_n, 25)}")

# === Part 5: Compare CMF error with Catalan moments ===
# The CMF error ê_n should be related to some Catalan moment.
# Find which sequence combination matches.
print("\n=== Comparing CMF error to Catalan moments ===")

# First, find the correct (P, Q) assignment by checking which ratio → G
best_j, best_k = None, None
best_diff = mp.mpf(1)
for j in range(3):
    for k in range(3):
        if j == k:
            continue
        r = seqs[j][NMAX] / seqs[k][NMAX]
        d = abs(r - G)
        if d < best_diff:
            best_diff = d
            best_j, best_k = j, k

print(f"Best: seq[{best_j}]/seq[{best_k}] → G (diff = {mp.nstr(best_diff, 6)})")
print(f"  → P_n = seq[{best_j}], Q_n = seq[{best_k}]")

# Compute ê_n = G·Q_n - P_n (unnormalized)
e_n = [G * seqs[best_k][n] - seqs[best_j][n] for n in range(min(12, NMAX+1))]
print(f"\nCMF error ê_n (unnormalized):")
for n in range(min(12, len(e_n))):
    print(f"  n={n}: ê_n = {mp.nstr(e_n[n], 20)}")

# Normalize by H_n
e_hat = [e_n[n] / H_n(n) for n in range(min(12, len(e_n)))]
print(f"\nH_n-normalized error ê̂_n:")
for n in range(min(12, len(e_hat))):
    print(f"  n={n}: ê̂_n = {mp.nstr(e_hat[n], 20)}")

# Compare ê̂_n with Catalan moments C[P_n]
print(f"\nRatio ê̂_n / C[P_n(1-2t²)]:")
for n in range(min(12, len(e_hat), len(cat_moments))):
    if cat_moments[n] != 0:
        ratio = e_hat[n] / cat_moments[n]
        print(f"  n={n}: ê̂_n / C[P_n] = {mp.nstr(ratio, 15)}")

# Try: ê̂_n / D_n² (Delannoy squared)
print(f"\nRatio ê̂_n / D_n²:")
for n in range(min(12, len(e_hat))):
    Dn = legendre_at_3(n)
    if Dn != 0:
        ratio = e_hat[n] / Dn**2
        print(f"  n={n}: ê̂_n / D_n² = {mp.nstr(ratio, 15)}")

# Try: ê_n (unnormalized) vs D_n²
print(f"\nRatio ê_n / D_n²:")
for n in range(min(12, len(e_n))):
    Dn = legendre_at_3(n)
    if Dn != 0:
        ratio = e_n[n] / Dn**2
        print(f"  n={n}: ê_n / D_n² = {mp.nstr(ratio, 15)}")

print("\nDone.")
