#!/usr/bin/env python3
"""Compute deformed tetrahedron shapes for the 7_2 knot along the
real A-polynomial branch, and evaluate the Rogers dilogarithm at endpoints.

Gluing equations from SnapPy (4 tetrahedra):
  edge 0: [1,1,0, 2,0,0, 1,0,1, 0,0,0]
  edge 1: [0,1,0, 0,0,1, 1,0,0, 0,0,2]
  edge 2: [1,0,1, 0,2,1, 0,1,1, 1,0,0]
  edge 3: [0,0,1, 0,0,0, 0,1,0, 1,2,0]
  merid:  [0,0,0,-1,0,0, 0,1,0, 0,0,0]
  long:   [0,-2,-1,-2,2,2, 0,3,-2, 0,0,0]

Column format: for tet j, (a_j, b_j, c_j) are exponents of
  z_j^a * (1/(1-z_j))^b * ((z_j-1)/z_j)^c

Meridian: z1^{-1} * (1/(1-z2))^1 = M^2
  => M^2 = 1/((1-z2)*z1)

Longitude: z0'^{-2} z0''^{-1} z1^{-2} z1'^2 z1''^2 z2'^3 z2''^{-2} = L
"""
import mpmath as mp

mp.mp.dps = 50

# Gluing equations matrix (each row = edge/cusp equation)
# Columns: (a0,b0,c0, a1,b1,c1, a2,b2,c2, a3,b3,c3)
GLUING = [
    [1,1,0, 2,0,0, 1,0,1, 0,0,0],  # edge 0
    [0,1,0, 0,0,1, 1,0,0, 0,0,2],  # edge 1
    [1,0,1, 0,2,1, 0,1,1, 1,0,0],  # edge 2
    [0,0,1, 0,0,0, 0,1,0, 1,2,0],  # edge 3
    [0,0,0,-1,0,0, 0,1,0, 0,0,0],  # meridian
    [0,-2,-1,-2,2,2, 0,3,-2, 0,0,0],  # longitude
]


def edge_product(z, row):
    """Compute product of z^a * z'^b * z''^c for all tetrahedra."""
    prod = mp.mpf(1)
    for j in range(4):
        a, b, c = row[3*j], row[3*j+1], row[3*j+2]
        zj = z[j]
        if a != 0:
            prod *= zj ** a
        if b != 0:
            prod *= (1 / (1 - zj)) ** b
        if c != 0:
            prod *= ((zj - 1) / zj) ** c
    return prod


def edge_equations(z, M2):
    """Return residuals for 3 independent edge equations + meridian."""
    # Edge 0: z0*z0'*z1^2*z2*z2'' = 1
    eq0 = edge_product(z, GLUING[0]) - 1
    # Edge 1: z0'*z1''*z2*z3''^2 = 1
    eq1 = edge_product(z, GLUING[1]) - 1
    # Edge 3: z0''*z2'*z3*z3'^2 = 1
    eq3 = edge_product(z, GLUING[3]) - 1
    # Meridian: z1^{-1} * z2' = M^2
    merid = edge_product(z, GLUING[4]) - M2
    return [eq0, eq1, eq3, merid]


def solve_shapes(M2_val, z_init):
    """Solve for shapes given M^2 value."""
    def sys(z0, z1, z2, z3):
        z = [z0, z1, z2, z3]
        return edge_equations(z, M2_val)

    result = mp.findroot(sys, z_init, tol=mp.power(10, -40))
    return list(result)


def log_holonomy(z, row):
    """Compute sum of a*log(z) + b*log(z') + c*log(z'') for all tets."""
    total = mp.mpf(0)
    for j in range(4):
        a, b, c = row[3*j], row[3*j+1], row[3*j+2]
        zj = z[j]
        if a != 0:
            total += a * mp.log(zj)
        if b != 0:
            total += b * mp.log(mp.mpf(1) / (1 - zj))
        if c != 0:
            total += c * mp.log((zj - 1) / zj)
    return total


def rogers_dilog(z):
    """Rogers dilogarithm R(z) = Li2(z) + (1/2)*log(z)*log(1-z)."""
    return mp.polylog(2, z) + mp.log(z) * mp.log(1 - z) / 2


def main():
    # First find shapes at the geometric point (complex)
    # Geometric shapes from SnapPy:
    z_geom = [
        mp.mpc('0.97968392714', '0.59056955984'),
        mp.mpc('0.25132270106', '0.45131497073'),
        mp.mpc('0.0581813774', '1.6912791495'),
        mp.mpc('1.1636911715', '0.5641856323'),
    ]

    # Verify edge equations at geometric point (M=1 for complete structure)
    print("=== Verification at geometric point ===")
    for i in range(4):
        prod = edge_product(z_geom, GLUING[i])
        print(f"Edge {i}: product = {mp.nstr(prod, 15)} (should be 1)")

    M2_geom = edge_product(z_geom, GLUING[4])
    L_geom = edge_product(z_geom, GLUING[5])
    print(f"M^2 = {mp.nstr(M2_geom, 15)} (should be 1)")
    print(f"L = {mp.nstr(L_geom, 15)} (should be 1)")

    # Endpoints from the A-polynomial
    s_alpha = mp.findroot(
        lambda s: mp.polyval(
            [1, -2, -3, 2, 2, 8, 6, 1, 5, -4, 0, 0, 0,
             -4, 5, 1, 6, 8, 2, 2, -3, -2, 1, 0, 0],
            s  # This isn't right... use the verify script's approach
        ),
        mp.mpf('0.591'),
    )

    # Actually let me just use hardcoded high-precision values from the verify script
    from verify_3_1_apoly import A, F_alpha_exact, F_beta_exact

    s_alpha = mp.findroot(F_alpha_exact, (mp.mpf('0.58'), mp.mpf('0.60')),
                          tol=mp.power(10, -40))
    s_beta = mp.findroot(F_beta_exact, (mp.mpf('0.40'), mp.mpf('0.42')),
                         tol=mp.power(10, -40))

    M_alpha = s_alpha**2
    L_alpha = s_alpha
    M_beta = s_beta
    L_beta = s_beta

    print(f"\n=== Endpoints ===")
    print(f"s_alpha = {mp.nstr(s_alpha, 30)}")
    print(f"M_alpha = {mp.nstr(M_alpha, 30)}, L_alpha = {mp.nstr(L_alpha, 30)}")
    print(f"s_beta = {mp.nstr(s_beta, 30)}")
    print(f"M_beta = L_beta = {mp.nstr(M_beta, 30)}")

    M2_alpha = M_alpha**2
    M2_beta = M_beta**2

    # Now try to find REAL shapes at the endpoints
    # At the beta endpoint, shapes should be real
    # Start with initial guesses near 0.5
    print(f"\n=== Solving for shapes at beta endpoint ===")
    print(f"M^2 = {mp.nstr(M2_beta, 20)}")

    # Try various initial guesses for real shapes
    for trial, z_init in enumerate([
        (0.3, 0.5, 0.7, 0.4),
        (0.1, 0.8, 0.9, 0.2),
        (0.9, 0.2, 0.1, 0.8),
        (0.5, 0.3, 0.6, 0.5),
        (2.0, 0.1, 0.3, -0.5),
        (-0.5, 0.8, 0.95, 3.0),
        (0.4, 0.6, 0.8, 0.3),
        (1.5, 0.4, 0.2, -1.0),
    ]):
        try:
            z_init_mp = [mp.mpf(x) for x in z_init]
            z_sol = solve_shapes(M2_beta, z_init_mp)
            # Check if solution is real
            is_real = all(abs(mp.im(z)) < 1e-30 for z in z_sol)
            # Verify edge equations
            resid = max(abs(r) for r in edge_equations(z_sol, M2_beta))
            if resid < 1e-30:
                print(f"  Trial {trial}: CONVERGED (resid={float(resid):.1e}, real={is_real})")
                for j in range(4):
                    print(f"    z_{j} = {mp.nstr(z_sol[j], 25)}")
                # Check longitude
                L_computed = edge_product(z_sol, GLUING[5])
                print(f"    L computed = {mp.nstr(L_computed, 20)}")
                print(f"    L target  = {mp.nstr(L_beta, 20)}")
                print(f"    L diff    = {mp.nstr(abs(L_computed - L_beta), 10)}")

                if is_real:
                    # Compute Rogers dilogarithm sum
                    R_sum = sum(rogers_dilog(mp.re(z)) for z in z_sol)
                    print(f"    Sum R(z_j) = {mp.nstr(R_sum, 25)}")
                    print(f"    pi^2/6     = {mp.nstr(mp.pi**2/6, 25)}")
                print()
        except Exception as e:
            pass

    print(f"\n=== Solving for shapes at alpha endpoint ===")
    print(f"M^2 = {mp.nstr(M2_alpha, 20)}")

    for trial, z_init in enumerate([
        (0.3, 0.5, 0.7, 0.4),
        (0.1, 0.8, 0.9, 0.2),
        (0.9, 0.2, 0.1, 0.8),
        (0.5, 0.3, 0.6, 0.5),
        (2.0, 0.1, 0.3, -0.5),
        (-0.5, 0.8, 0.95, 3.0),
        (0.4, 0.6, 0.8, 0.3),
        (1.5, 0.4, 0.2, -1.0),
    ]):
        try:
            z_init_mp = [mp.mpf(x) for x in z_init]
            z_sol = solve_shapes(M2_alpha, z_init_mp)
            is_real = all(abs(mp.im(z)) < 1e-30 for z in z_sol)
            resid = max(abs(r) for r in edge_equations(z_sol, M2_alpha))
            if resid < 1e-30:
                print(f"  Trial {trial}: CONVERGED (resid={float(resid):.1e}, real={is_real})")
                for j in range(4):
                    print(f"    z_{j} = {mp.nstr(z_sol[j], 25)}")
                L_computed = edge_product(z_sol, GLUING[5])
                print(f"    L computed = {mp.nstr(L_computed, 20)}")
                print(f"    L target  = {mp.nstr(L_alpha, 20)}")
                print(f"    L diff    = {mp.nstr(abs(L_computed - L_alpha), 10)}")
                print()
        except Exception as e:
            pass


if __name__ == "__main__":
    main()
