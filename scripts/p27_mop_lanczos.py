#!/usr/bin/env python3
"""P2.7: Compute step-line MOP polynomials using the Lanczos/Gauss-Borel approach.
This avoids the ill-conditioned determinant formula.

For the measures μ₂ = (-log t)dt, μ₃ = ½(log²t)dt on [0,1],
compute the type-II step-line MOP polynomials Q_N(z) using
numerical integration for the inner products.

The recurrence coefficients are computed by Gram-Schmidt orthogonalization.
"""
import mpmath as mp
mp.mp.dps = 250

def inner_product_2(f_coeffs, g_coeffs):
    """<f, g>₂ = ∫₀¹ f(t)·g(t)·(-log t) dt"""
    # f(t) = Σ f_k t^k, g(t) = Σ g_k t^k
    # <f,g>₂ = Σ f_j g_k ∫₀¹ t^{j+k}(-log t) dt = Σ f_j g_k / (j+k+1)²
    val = mp.mpf(0)
    for j, fj in enumerate(f_coeffs):
        for k, gk in enumerate(g_coeffs):
            val += fj * gk / mp.mpf(j+k+1)**2
    return val

def inner_product_3(f_coeffs, g_coeffs):
    """<f, g>₃ = ∫₀¹ f(t)·g(t)·½(log²t) dt"""
    val = mp.mpf(0)
    for j, fj in enumerate(f_coeffs):
        for k, gk in enumerate(g_coeffs):
            val += fj * gk / mp.mpf(j+k+1)**3
    return val

def poly_mult_z(coeffs):
    """Multiply polynomial by z: [c₀, c₁, ...] → [0, c₀, c₁, ...]"""
    return [mp.mpf(0)] + list(coeffs)

def poly_add(a, b):
    """Add two polynomials"""
    n = max(len(a), len(b))
    result = [mp.mpf(0)] * n
    for i in range(len(a)):
        result[i] += a[i]
    for i in range(len(b)):
        result[i] += b[i]
    return result

def poly_scale(coeffs, s):
    """Scale polynomial by scalar s"""
    return [c * s for c in coeffs]

def poly_eval(coeffs, z):
    """Evaluate polynomial at z"""
    val = mp.mpf(0)
    zp = mp.mpf(1)
    for c in coeffs:
        val += c * zp
        zp *= z
    return val

# Step line: alternating (n₂, n₃)
# N=0: (0,0), N=1: (1,0), N=2: (1,1), N=3: (2,1), N=4: (2,2), ...
def step_multi_index(N):
    n2 = (N + 1) // 2
    n3 = N // 2
    return n2, n3

# The step-line MOP Q_N has degree N and satisfies:
# ∫₀¹ t^k Q_N(t) (-log t) dt = 0 for k = 0, ..., n₂-1
# ∫₀¹ t^k Q_N(t) ½(log²t) dt = 0 for k = 0, ..., n₃-1

# Use Gram-Schmidt to build Q_N incrementally.
# The nearest-neighbor recurrence on the step line:
# z·Q_N(z) = Q_{N+1}(z) + β_N·Q_N(z) + α_N·Q_{N-1}(z) + γ_N·Q_{N-2}(z)
#
# But the exact recurrence structure depends on the step pattern.
# Let me use a direct Gram-Schmidt approach instead.

# Gram-Schmidt: start with monomials 1, z, z², ...
# Project each onto the orthogonal complement of the previous Q_j's
# under ALL the inner products <·,·>₂ and <·,·>₃ according to the step line

NMAX = 80

# Store Q polynomials and their inner products
Q_polys = []  # Q_polys[N] = list of coefficients
Q_at_1 = []

# To build Q_N: start with z^N, subtract projections
# The orthogonality conditions for Q_N on the step line (n₂, n₃):
# <t^k, Q_N>₂ = 0 for k = 0,...,n₂-1
# <t^k, Q_N>₃ = 0 for k = 0,...,n₃-1
# These are N linear conditions on N+1 coefficients (degree-N polynomial).
# With the monic normalization (leading coeff = 1), this determines Q_N uniquely.

print("=== Computing step-line MOP via constrained linear system ===", flush=True)

for N in range(NMAX + 1):
    n2, n3 = step_multi_index(N)

    # Q_N(t) = t^N + c_{N-1} t^{N-1} + ... + c_0 (monic)
    # Conditions:
    # ∫₀¹ t^k Q_N(t) (-log t) dt = 0 for k = 0,...,n₂-1
    # ∫₀¹ t^k Q_N(t) ½(log²t) dt = 0 for k = 0,...,n₃-1
    # Total: n₂ + n₃ = N conditions, N unknowns c₀,...,c_{N-1}

    if N == 0:
        Q_polys.append([mp.mpf(1)])
        Q_at_1.append(mp.mpf(1))
        print(f"  N={N}: Q(1) = {mp.nstr(Q_at_1[-1], 15)}", flush=True)
        continue

    # Build the linear system
    # Unknown: c = [c_0, c_1, ..., c_{N-1}]
    # Q_N(t) = c_0 + c_1 t + ... + c_{N-1} t^{N-1} + t^N

    # Condition from μ₂: ∫₀¹ t^k · (c_0 + ... + c_{N-1} t^{N-1} + t^N) · (-log t) dt = 0
    # = Σ_j c_j ∫₀¹ t^{k+j} (-log t) dt + ∫₀¹ t^{k+N} (-log t) dt = 0
    # = Σ_j c_j / (k+j+1)² + 1/(k+N+1)² = 0

    # Similarly for μ₃: Σ_j c_j / (k+j+1)³ + 1/(k+N+1)³ = 0

    A_mat = mp.matrix(N, N)
    b_vec = mp.matrix(N, 1)

    row = 0
    for k in range(n2):
        for j in range(N):
            A_mat[row, j] = mp.mpf(1) / mp.mpf(k+j+1)**2
        b_vec[row] = -mp.mpf(1) / mp.mpf(k+N+1)**2
        row += 1

    for k in range(n3):
        for j in range(N):
            A_mat[row, j] = mp.mpf(1) / mp.mpf(k+j+1)**3
        b_vec[row] = -mp.mpf(1) / mp.mpf(k+N+1)**3
        row += 1

    assert row == N, f"Expected {N} conditions, got {row}"

    try:
        c = mp.lu_solve(A_mat, b_vec)
        coeffs = [c[j] for j in range(N)] + [mp.mpf(1)]
        Q_polys.append(coeffs)
        Q_at_1.append(poly_eval(coeffs, mp.mpf(1)))

        if N <= 20 or N % 5 == 0:
            print(f"  N={N}: Q(1) = {mp.nstr(Q_at_1[-1], 15)}, ({n2},{n3})", flush=True)
    except Exception as e:
        print(f"  N={N}: FAILED ({e})", flush=True)
        break

if len(Q_at_1) < 7:
    print("Not enough values for recurrence search.")
    import sys; sys.exit(0)

# Compute ratios
print(f"\n=== Q_N(1) ratios ===")
for N in range(1, len(Q_at_1)):
    if Q_at_1[N-1] != 0:
        r = Q_at_1[N] / Q_at_1[N-1]
        if N < 10 or N % 5 == 0:
            print(f"  Q_{N}/Q_{N-1} = {mp.nstr(r, 15)}")

# 4-term recurrence search
print(f"\n=== Searching for 4-term recurrence ===")
nvals = len(Q_at_1)
for deg in range(1, 15):
    nparams = 4 * (deg + 1) - 1
    neq = nvals - 3
    if neq < nparams + 3:
        continue

    A_sys = []
    b_sys = []
    for n in range(neq):
        row = []
        for j in range(4):
            for k in range(deg + 1):
                if j == 3 and k == deg:
                    continue
                row.append(mp.mpf(n)**k * Q_at_1[n + j])
        A_sys.append(row)
        b_sys.append(-mp.mpf(n)**deg * Q_at_1[n + 3])

    A_solve = mp.matrix([r[:nparams] for r in A_sys[:nparams]])
    b_solve = mp.matrix([b_sys[i] for i in range(nparams)])

    try:
        sol = mp.lu_solve(A_solve, b_solve)
    except:
        print(f"  degree {deg}: singular")
        continue

    # Verify on held-out equations
    max_res = mp.mpf(0)
    for n in range(nparams, min(neq, nparams + 5)):
        res = b_sys[n]
        for i in range(nparams):
            res -= A_sys[n][i] * sol[i]
        max_res = max(max_res, abs(res))

    print(f"  degree {deg}: max holdout residual = {mp.nstr(max_res, 6)}", flush=True)

    if max_res < mp.mpf(10)**(-100):
        print(f"  *** FOUND polynomial degree {deg} ***")

        coeffs = {}
        idx = 0
        for j in range(4):
            for k in range(deg + 1):
                if j == 3 and k == deg:
                    coeffs[(j, k)] = mp.mpf(1)
                else:
                    coeffs[(j, k)] = sol[idx]
                    idx += 1

        lc = [coeffs[(j, deg)] for j in range(4)]
        print(f"  Leading coefficients: {[mp.nstr(c, 10) for c in lc]}")

        char_poly = [lc[j]/lc[3] for j in range(4)]
        print(f"  Poincaré: λ³ + {mp.nstr(char_poly[2], 10)}λ² + {mp.nstr(char_poly[1], 10)}λ + {mp.nstr(char_poly[0], 10)} = 0")

        roots = mp.polyroots([char_poly[3], char_poly[2], char_poly[1], char_poly[0]])
        print(f"  Roots: {[mp.nstr(r, 10) for r in roots]}")
        print(f"  |Roots|: {[mp.nstr(abs(r), 10) for r in roots]}")
        break

print("\nDone.")
