#!/usr/bin/env python3
"""
Problem 2.5: Algebraic Ore intertwiner via the E0-E5 equations from Q4792.

Set up L₂₅♯ ∘ T = R ∘ L_Sym² as an operator identity in Q(n)[S].
Eliminate r₀, r₁, r₂ using E0, E1, E5.
Remaining E2, E3, E4 are 3 functional equations for t₀, t₁, t₂.

Test order 1 first (t₂ = 0), then order 2.
Use polynomial ansatz with increasing degree.
"""
from sympy import (Symbol, Rational, factor, expand, Poly, cancel,
                   simplify, gcd as sgcd, lcm as slcm, collect,
                   degree, numer, denom, div, symbols, solve)

n = Symbol('n')

# ---- L_Sym² coefficients ----
Q35 = 35*n**2 + 140*n + 131
a0 = -(2*n+5)*(n+1)**2
a1 = (2*n+5)*Q35
a2 = -(2*n+3)*Q35
a3 = (2*n+3)*(n+3)**2

# ---- L₂₅♯ coefficients (from twisted_L25_extract.py) ----
# Extracted with normalization: leading coeff of ℓ₃ = 1 (in the extraction)
# Let me use the UN-scaled versions (multiply by the denom factor later)

# From the output:
# ℓ₀ = -(N+1)(N+3)²(N+4)(2N+3)(2N+5)(2N+7) · P₆(N) / 24576
# ℓ₁ = (N+2)²(N+4)(2N+7) · P₉(N) / 98304
# ℓ₂ = -(N+2)(N+3)(2N+9) · P₁₀(N) / 49152
# ℓ₃ = (N+2)(N+3)(N+4)²(2N+7)(2N+11)² · P₆'(N) / 24576

P6 = 3072*n**6 + 74112*n**5 + 738544*n**4 + 3890106*n**3 + 11417947*n**2 + 17696904*n + 11307715
P9 = (1720320*n**9 + 53692416*n**8 + 735843584*n**7 + 5810223840*n**6
      + 29119642544*n**5 + 96028512072*n**4 + 208287850700*n**3
      + 286416594222*n**2 + 226466421477*n + 78407415225)
P10 = (860160*n**10 + 30999552*n**9 + 498082432*n**8 + 4696686576*n**7
       + 28770039448*n**6 + 119560001580*n**5 + 341146546318*n**4
       + 659435701854*n**3 + 825643834707*n**2 + 603914277213*n + 195670909710)
P6p = 3072*n**6 + 55680*n**5 + 414064*n**4 + 1615610*n**3 + 3483853*n**2 + 3929280*n + 1806156

# Scale to clear the 1/24576 etc factors. Multiply all by 98304 (= LCM of 24576, 98304, 49152)
# 98304 / 24576 = 4
# 98304 / 98304 = 1
# 98304 / 49152 = 2
# 98304 / 24576 = 4
ell0 = -4*(n+1)*(n+3)**2*(n+4)*(2*n+3)*(2*n+5)*(2*n+7)*P6
ell1 = (n+2)**2*(n+4)*(2*n+7)*P9
ell2 = -2*(n+2)*(n+3)*(2*n+9)*P10
ell3 = 4*(n+2)*(n+3)*(n+4)**2*(2*n+7)*(2*n+11)**2*P6p

print("=== Checking P₆(n) = P₆'(n+1) ===")
diff = expand(P6 - P6p.subs(n, n+1))
print(f"P₆(n) - P₆'(n+1) = {diff}")
print(f"{'✓ Confirmed!' if diff == 0 else '✗ MISMATCH!'}")

# ---- Shift function ----
def sh(expr, k):
    """Shift: f(n) -> f(n+k)"""
    return expand(expr.subs(n, n+k))

# ---- Ore identity: L₂₅♯ ∘ T = R ∘ L_Sym² ----
# Equations E0-E5:
# E0: ℓ₀t₀ = r₀a₀
# E1: ℓ₀t₁ + ℓ₁·t₀[1] = r₀a₁ + r₁·a₀[1]
# E2: ℓ₀t₂ + ℓ₁·t₁[1] + ℓ₂·t₀[2] = r₀a₂ + r₁·a₁[1] + r₂·a₀[2]
# E3: ℓ₁·t₂[1] + ℓ₂·t₁[2] + ℓ₃·t₀[3] = r₀a₃ + r₁·a₂[1] + r₂·a₁[2]
# E4: ℓ₂·t₂[2] + ℓ₃·t₁[3] = r₁·a₃[1] + r₂·a₂[2]
# E5: ℓ₃·t₂[3] = r₂·a₃[2]

# From E0: r₀ = ℓ₀t₀/a₀
# From E5: r₂ = ℓ₃·t₂[3]/a₃[2]
# From E1: r₁ = (ℓ₀t₁ + ℓ₁·t₀[1] - r₀·a₁)/a₀[1]

print("\n=== Testing ORDER 1 (t₂ = 0) ===")

# With t₂ = 0, r₂ = 0
# r₀ = ℓ₀t₀/a₀
# r₁ = (ℓ₀t₁ + ℓ₁·t₀[1] - ℓ₀t₀a₁/a₀) / a₀[1]
# E2: ℓ₁·t₁[1] + ℓ₂·t₀[2] = r₀·a₂ + r₁·a₁[1]
# E3: ℓ₂·t₁[2] + ℓ₃·t₀[3] = r₀·a₃ + r₁·a₂[1]
# E4: ℓ₃·t₁[3] = r₁·a₃[1]

# Substitute r₀ and r₁, multiply through to clear denominators.
# The key denominator from r₀ is a₀ = -(2n+5)(n+1)²
# and from r₁ is a₀(n+1) = -(2n+7)(n+2)²

# Try polynomial ansatz for t₀, t₁ with increasing degree
for d in range(12):
    t0_coeffs = symbols(f't0_0:{d+1}')
    t1_coeffs = symbols(f't1_0:{d+1}')

    t0_poly = sum(t0_coeffs[m] * n**m for m in range(d+1))
    t1_poly = sum(t1_coeffs[m] * n**m for m in range(d+1))

    # r₀ = ℓ₀ · t₀ / a₀
    r0_num = ell0 * t0_poly
    r0_den = a0

    # r₁ = (ℓ₀·t₁ + ℓ₁·t₀[1] - r₀·a₁) / a₀[1]
    # = (ℓ₀·t₁ + ℓ₁·t₀[1] - ℓ₀·t₀·a₁/a₀) / a₀[1]
    # = (ℓ₀·t₁·a₀ + ℓ₁·t₀[1]·a₀ - ℓ₀·t₀·a₁) / (a₀·a₀[1])
    r1_num = expand(ell0 * t1_poly * a0 + ell1 * sh(t0_poly, 1) * a0 - ell0 * t0_poly * a1)
    r1_den = expand(a0 * sh(a0, 1))

    # E4: ℓ₃·t₁[3] = r₁·a₃[1]
    # ℓ₃·t₁[3] · a₀·a₀[1] = r1_num · a₃[1]
    E4_lhs = expand(ell3 * sh(t1_poly, 3) * a0 * sh(a0, 1))
    E4_rhs = expand(r1_num * sh(a3, 1))
    E4_eq = expand(E4_lhs - E4_rhs)

    # Collect as polynomial in n
    E4_poly = Poly(E4_eq, n)

    # Each coefficient of n^k must be zero
    E4_system = E4_poly.all_coeffs()

    all_vars = list(t0_coeffs) + list(t1_coeffs)

    # Filter to linear equations
    linear_eqs = []
    for eq in E4_system:
        if eq == 0:
            continue
        # Check it's linear in the unknowns
        p = Poly(eq, *all_vars)
        if p.total_degree() <= 1:
            linear_eqs.append(eq)
        else:
            break
    else:
        # All equations are linear
        sol = solve(linear_eqs, all_vars, dict=True)
        if sol:
            print(f"  d={d}: E4 has solution with {len(sol)} parameter families")
            # Check E2 and E3 too

            # E2: ℓ₁·t₁[1] + ℓ₂·t₀[2] = r₀·a₂ + r₁·a₁[1]
            # Multiply through by a₀·a₀[1]:
            E2_lhs = expand((ell1 * sh(t1_poly, 1) + ell2 * sh(t0_poly, 2)) * a0 * sh(a0, 1))
            E2_rhs_r0 = expand(ell0 * t0_poly * a2 * sh(a0, 1) / a0 * a0 * sh(a0, 1))  # hmm this is wrong

            # Let me redo more carefully.
            # E2: ℓ₁t₁[1] + ℓ₂t₀[2] - (ℓ₀t₀/a₀)a₂ - r₁·a₁[1] = 0
            # Multiply by a₀·a₀[1]:
            # (ℓ₁t₁[1] + ℓ₂t₀[2])·a₀·a₀[1] - ℓ₀t₀·a₂·a₀[1] - r1_num·a₁[1] = 0
            E2_eq = expand(
                (ell1 * sh(t1_poly, 1) + ell2 * sh(t0_poly, 2)) * a0 * sh(a0, 1)
                - ell0 * t0_poly * a2 * sh(a0, 1)
                - r1_num * sh(a1, 1)
            )

            E2_poly = Poly(E2_eq, n)
            E2_system = E2_poly.all_coeffs()

            # E3: ℓ₂t₁[2] + ℓ₃t₀[3] - (ℓ₀t₀/a₀)a₃ - r₁·a₂[1] = 0
            E3_eq = expand(
                (ell2 * sh(t1_poly, 2) + ell3 * sh(t0_poly, 3)) * a0 * sh(a0, 1)
                - ell0 * t0_poly * a3 * sh(a0, 1)
                - r1_num * sh(a2, 1)
            )

            E3_poly = Poly(E3_eq, n)
            E3_system = E3_poly.all_coeffs()

            # Combine all equations
            all_eqs = list(E4_system) + list(E2_system) + list(E3_system)
            all_eqs = [eq for eq in all_eqs if eq != 0]

            print(f"  d={d}: total equations from E2+E3+E4: {len(all_eqs)}, unknowns: {len(all_vars)}")

            full_sol = solve(all_eqs, all_vars, dict=True)
            if full_sol:
                print(f"  d={d}: *** SOLUTION FOUND! ***")
                for s in full_sol[:1]:
                    print(f"    {s}")
                break
            else:
                print(f"  d={d}: E4 has solutions but full system (E2+E3+E4) is inconsistent")
        else:
            if d <= 5:
                print(f"  d={d}: E4 alone has no solution")
        continue

    if d <= 5:
        print(f"  d={d}: E4 has nonlinear terms (unexpected)")

print("\n=== Testing ORDER 2 (t₀, t₁, t₂ all nonzero) ===")
for d in range(6):
    t0_c = symbols(f'p0_0:{d+1}')
    t1_c = symbols(f'p1_0:{d+1}')
    t2_c = symbols(f'p2_0:{d+1}')

    t0_p = sum(t0_c[m] * n**m for m in range(d+1))
    t1_p = sum(t1_c[m] * n**m for m in range(d+1))
    t2_p = sum(t2_c[m] * n**m for m in range(d+1))

    # r₀ = ℓ₀t₀/a₀
    # r₂ = ℓ₃t₂[3]/a₃[2]
    # r₁ = (ℓ₀t₁ + ℓ₁t₀[1] - ℓ₀t₀a₁/a₀) / a₀[1]
    # = (ℓ₀t₁a₀ + ℓ₁t₀[1]a₀ - ℓ₀t₀a₁) / (a₀·a₀[1])

    r1_num_expr = expand(ell0 * t1_p * a0 + ell1 * sh(t0_p, 1) * a0 - ell0 * t0_p * a1)
    r2_num_expr = expand(ell3 * sh(t2_p, 3))
    r2_den_expr = sh(a3, 2)

    # E2: (ℓ₀t₂ + ℓ₁t₁[1] + ℓ₂t₀[2])·a₀·a₀[1]·a₃[2]
    #   - ℓ₀t₀·a₂·a₀[1]·a₃[2]
    #   - r1_num·a₁[1]·a₃[2]
    #   - r2_num·a₀[2]·a₀·a₀[1] = 0

    common = a0 * sh(a0, 1) * r2_den_expr

    E2_eq = expand(
        (ell0 * t2_p + ell1 * sh(t1_p, 1) + ell2 * sh(t0_p, 2)) * common
        - ell0 * t0_p * a2 * sh(a0, 1) * r2_den_expr
        - r1_num_expr * sh(a1, 1) * r2_den_expr
        - r2_num_expr * sh(a0, 2) * a0 * sh(a0, 1)
    )

    E3_eq = expand(
        (ell1 * sh(t2_p, 1) + ell2 * sh(t1_p, 2) + ell3 * sh(t0_p, 3)) * common
        - ell0 * t0_p * a3 * sh(a0, 1) * r2_den_expr
        - r1_num_expr * sh(a2, 1) * r2_den_expr
        - r2_num_expr * sh(a1, 2) * a0 * sh(a0, 1)
    )

    E4_eq = expand(
        (ell2 * sh(t2_p, 2) + ell3 * sh(t1_p, 3)) * common
        - r1_num_expr * sh(a3, 1) * r2_den_expr
        - r2_num_expr * sh(a2, 2) * a0 * sh(a0, 1)
    )

    all_vars = list(t0_c) + list(t1_c) + list(t2_c)

    print(f"\n  d={d}: setting up equations...")

    E2_poly = Poly(E2_eq, n)
    E3_poly = Poly(E3_eq, n)
    E4_poly = Poly(E4_eq, n)

    all_eqs = list(E2_poly.all_coeffs()) + list(E3_poly.all_coeffs()) + list(E4_poly.all_coeffs())
    all_eqs = [eq for eq in all_eqs if eq != 0]

    print(f"  d={d}: {len(all_eqs)} equations, {len(all_vars)} unknowns")

    sol = solve(all_eqs, all_vars, dict=True)
    if sol:
        print(f"  d={d}: *** SOLUTION FOUND! ***")
        for s in sol[:1]:
            for var, val in s.items():
                print(f"    {var} = {val}")
        break
    else:
        print(f"  d={d}: no solution")

print("\nDone.")
