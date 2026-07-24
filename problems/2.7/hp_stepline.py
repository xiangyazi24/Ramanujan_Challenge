"""
Two-measure step-line Hermite-Padé for simultaneous ζ(2) + ζ(3).

Measures: dμ₂ = -log(t) dt, dμ₃ = (1/2)log²(t) dt on [0,1]
Moments: μ₂(k) = 1/(k+1)², μ₃(k) = 1/(k+1)³

Step-line path in (d₂, d₃):
(0,0) → (1,0) → (1,1) → (2,1) → (2,2) → (3,2) → ...

At step m with (d₂, d₃), find Q_m(t) = Σ_{j=0}^m a_j t^j satisfying:
  ∫₀¹ Q_m(t) · t^k · dμ₂ = 0  for k = 0,...,d₂-1
  ∫₀¹ Q_m(t) · t^k · dμ₃ = 0  for k = 0,...,d₃-1
"""
from fractions import Fraction as Q

def mu2(k):
    """Moment of dμ₂: ∫₀¹ t^k · (-log t) dt = 1/(k+1)²"""
    return Q(1, (k+1)**2)

def mu3(k):
    """Moment of dμ₃: ∫₀¹ t^k · (1/2)log²(t) dt = 1/(k+1)³"""
    return Q(1, (k+1)**3)

def stepline_degrees(m):
    """(d₂, d₃) at step m on the step-line path."""
    # Path: (0,0), (1,0), (1,1), (2,1), (2,2), (3,2), (3,3), ...
    if m == 0:
        return (0, 0)
    d2 = (m + 1) // 2
    d3 = m // 2
    return (d2, d3)

def compute_hp_poly(m):
    """
    Compute the HP polynomial Q_m(t) = Σ a_j t^j of degree m.
    Normalized so that the leading coefficient is 1.
    """
    d2, d3 = stepline_degrees(m)
    # Total orthogonality conditions: d2 + d3 = m
    # Number of coefficients: m + 1
    # So we have m equations and m+1 unknowns → 1D solution space

    # Set up matrix: A[i, :] · coeffs = 0
    # For dμ₂: ∫ (Σ a_j t^j) · t^k · dμ₂ = Σ a_j · μ₂(j+k) = 0 for k=0,...,d₂-1
    # For dμ₃: ∫ (Σ a_j t^j) · t^k · dμ₃ = Σ a_j · μ₃(j+k) = 0 for k=0,...,d₃-1

    n_vars = m + 1  # coefficients a_0, ..., a_m
    n_eqs = d2 + d3  # = m

    mat = []
    for k in range(d2):
        row = [mu2(j + k) for j in range(n_vars)]
        mat.append(row)
    for k in range(d3):
        row = [mu3(j + k) for j in range(n_vars)]
        mat.append(row)

    # Solve with a_m = 1 (monic normalization)
    # Move a_m to RHS: Σ_{j=0}^{m-1} A[i,j] · a_j = -A[i,m]
    if m == 0:
        return [Q(1)]

    lhs = [[mat[i][j] for j in range(m)] for i in range(n_eqs)]
    rhs = [-mat[i][m] for i in range(n_eqs)]

    # Gaussian elimination over Q
    aug = [lhs[i][:] + [rhs[i]] for i in range(n_eqs)]
    n_cols = m

    for col in range(min(n_eqs, n_cols)):
        # Find pivot
        pivot = None
        for row in range(col, n_eqs):
            if aug[row][col] != Q(0):
                pivot = row
                break
        if pivot is None:
            print(f"  m={m}: no pivot at col {col}")
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]
        # Eliminate
        for row in range(n_eqs):
            if row != col and aug[row][col] != Q(0):
                factor = aug[row][col] / aug[col][col]
                for j in range(n_cols + 1):
                    aug[row][j] -= factor * aug[col][j]

    # Extract solution
    coeffs = [Q(0)] * (m + 1)
    coeffs[m] = Q(1)  # monic
    for col in range(min(n_eqs, n_cols)):
        coeffs[col] = aug[col][n_cols] / aug[col][col]

    return coeffs

# Compute HP polynomials for several steps
print("=== Step-line HP polynomials ===")
b_vals = []  # Q_m(0) = a_0
for m in range(25):
    d2, d3 = stepline_degrees(m)
    coeffs = compute_hp_poly(m)
    if coeffs is None:
        print(f"  m={m}: FAILED")
        b_vals.append(None)
        continue
    b_m = coeffs[0]  # Q_m(0)
    b_vals.append(b_m)
    print(f"  m={m}: (d₂,d₃)=({d2},{d3}), Q_m(0) = {float(b_m):.10e}")

# Now compute the "error": b_m · (ζ(2)+ζ(3)) - p_m where p_m = Σ a_j (μ₂(j) + μ₃(j))
from math import pi
zeta2 = pi**2 / 6
zeta3_approx = 1.2020569031595942  # ζ(3)
L = zeta2 + zeta3_approx

print("\n=== Convergence to ζ(2)+ζ(3) ===")
for m in range(25):
    if b_vals[m] is None or b_vals[m] == 0:
        continue
    coeffs = compute_hp_poly(m)
    p2_m = sum(coeffs[j] * mu2(j) for j in range(m+1))
    p3_m = sum(coeffs[j] * mu3(j) for j in range(m+1))
    p_sum = p2_m + p3_m
    ratio = float(p_sum / b_vals[m])
    err = ratio - L
    print(f"  m={m}: (p₂+p₃)/b = {ratio:.15f}, err = {err:.6e}")

# Check if b_m satisfies a 4-term recurrence
print("\n=== Searching for 4-term recurrence (order 3) ===")
# b_{m+3} = α(m) b_{m+2} + β(m) b_{m+1} + γ(m) b_m
# Try polynomial coefficients of degree d
valid_b = [(i, b_vals[i]) for i in range(len(b_vals)) if b_vals[i] is not None]
if len(valid_b) >= 10:
    indices = [v[0] for v in valid_b]
    vals = [v[1] for v in valid_b]
    print(f"  Have {len(vals)} values")

    # Check ratios
    print("  Successive ratios b_{m+1}/b_m:")
    for i in range(min(15, len(vals)-1)):
        if vals[i] != 0:
            print(f"    m={indices[i]}: {float(vals[i+1]/vals[i]):.10f}")

    # Try: coefficients are polynomials in m of degree d
    # Equations: c₃(m)·b_{m+3} + c₂(m)·b_{m+2} + c₁(m)·b_{m+1} + c₀(m)·b_m = 0
    for deg in range(8):
        n_coeffs = 4 * (deg + 1)  # 4 polynomial coefficients, each degree d
        n_available = len(vals) - 3
        if n_available < n_coeffs + 2:
            break

        # Build system
        mat = []
        for idx in range(n_available):
            m = indices[idx]
            row = []
            for shift in range(4):  # b_{m+shift} for shift=0,1,2,3
                for power in range(deg + 1):
                    row.append(Q(m) ** power * vals[idx + shift])
            mat.append(row)

        # Find nullspace
        aug = [row[:] for row in mat]
        n_rows = len(aug)
        n_cols_mat = n_coeffs

        pivot_cols = []
        r = 0
        for c in range(n_cols_mat):
            pivot = None
            for i in range(r, n_rows):
                if aug[i][c] != Q(0):
                    pivot = i
                    break
            if pivot is None:
                continue
            aug[r], aug[pivot] = aug[pivot], aug[r]
            pivot_cols.append(c)
            for i in range(n_rows):
                if i != r and aug[i][c] != Q(0):
                    factor = aug[i][c] / aug[r][c]
                    for j in range(n_cols_mat):
                        aug[i][j] -= factor * aug[r][j]
            r += 1

        rank = r
        nullity = n_cols_mat - rank
        if nullity > 0:
            print(f"  deg={deg}: rank={rank}, nullity={nullity} — RECURRENCE EXISTS")
            # Extract the null vector
            free_cols = [c for c in range(n_cols_mat) if c not in pivot_cols]
            if free_cols:
                # Set first free variable to 1, rest to 0
                sol = [Q(0)] * n_cols_mat
                fc = free_cols[0]
                sol[fc] = Q(1)
                # Back-substitute
                for idx_r in range(rank - 1, -1, -1):
                    pc = pivot_cols[idx_r]
                    val_s = Q(0)
                    for c in range(n_cols_mat):
                        if c != pc:
                            val_s += aug[idx_r][c] * sol[c]
                    sol[pc] = -val_s / aug[idx_r][pc]

                # Print coefficients
                print(f"  Recurrence coefficients (deg={deg}):")
                for shift in range(4):
                    poly_coeffs = [sol[shift*(deg+1)+p] for p in range(deg+1)]
                    if any(c != 0 for c in poly_coeffs):
                        terms = []
                        for p, c in enumerate(poly_coeffs):
                            if c != 0:
                                terms.append(f"({c})*m^{p}")
                        print(f"    c_{shift}(m) = {' + '.join(terms)}")

                # Verify
                ok = True
                for idx in range(n_available):
                    m = indices[idx]
                    total = Q(0)
                    for shift in range(4):
                        c_val = Q(0)
                        for power in range(deg+1):
                            c_val += sol[shift*(deg+1)+power] * Q(m)**power
                        total += c_val * vals[idx+shift]
                    if total != Q(0):
                        ok = False
                        break
                print(f"    Verified: {'YES' if ok else 'NO'}")
            break
        else:
            print(f"  deg={deg}: full rank (no recurrence)")
