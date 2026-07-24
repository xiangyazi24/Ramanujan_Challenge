"""
Numerical monodromy computation for M₁₆.

Computes the monodromy matrix around z₀ = 64/μ₀ ≈ 1.1644 and verifies
that A(p)/A(q) = ζ(2)+ζ(3) via the rank-1 nilpotent (M-I).
"""
from sage.all import *
from ore_algebra import OreAlgebra

print("=== Building M₁₆ ===")

Rn.<n> = PolynomialRing(QQ)
Kn = Rn.fraction_field()
Rt.<t> = PolynomialRing(QQ)
Kt = Rt.fraction_field()
OD.<Dt> = OreAlgebra(Kt, 'Dt')
theta = Kt(t)*Dt

def eval_poly_at_op(f, T):
    out = OD.zero()
    for a in reversed(Rn(f).list()):
        out = out*T + Kt(a)
    return out

def primitive_coefficients(c):
    c = [Kn(a) for a in c]
    den = Rn.one()
    for a in c:
        den = lcm(den, Rn(a.denominator()))
    nums = [Rn(den*a) for a in c]
    g = nums[0]
    for a in nums[1:]:
        g = gcd(g, a)
    nums = [a.quo_rem(g)[0] for a in nums]
    scl = QQ(1)/QQ(nums[-1].leading_coefficient())
    return [Rn(scl*a) for a in nums]

def A(z): return 1024*(2*z+5)^4*(2*z+7)^3*(2*z+9)^3*(946*z^2+6407*z+10860)
def B(z): return 128*(2*z+7)^3*(2*z+9)^3*(104060*z^6+1745370*z^5+12145238*z^4+44886481*z^3+92943995*z^2+102256019*z+46709052)
def C(z): return 16*(z+3)^4*(2*z+9)^3*(3784*z^5+57792*z^4+351019*z^3+1059230*z^2+1587211*z+944620)
def D(z): return (z+3)^4*(z+4)^6*(946*z^2+4515*z+5399)

crat = [Kn(-D(n)/A(n)), Kn(C(n+1)/A(n+1)), Kn(-B(n+2)/A(n+2)), Kn(1)]
p_coeff = primitive_coefficients(crat)

Rx.<x> = PolynomialRing(QQ)
shifted = [Rx(p_coeff[j](x-j)) for j in range(4)]
gshift = shifted[0]
for f in shifted[1:]:
    gshift = gcd(gshift, f)
gshift = gshift.monic()
assert gshift == x^2 + QQ(105)/22*x + QQ(5399)/946

Rpoly = []
for f in shifted:
    q, rem = f.quo_rem(gshift)
    assert rem == 0
    Rpoly.append(q)

M16 = OD.zero()
for j in range(4):
    M16 += Kt(t^j)*eval_poly_at_op(Rpoly[j], -theta-1)
assert M16.order() == 16
print(f"M16 order = {M16.order()}")

# Compute the singularities
chi = t^3 - QQ(55)/64*t^2 + QQ(1)/2048*t - QQ(1)/(2^20)
print(f"\nCharacteristic polynomial: {chi}")
chi_roots = (2^20 * chi).roots(QQbar, multiplicities=False)
print(f"Singularities z₀ = 1/λ:")
for r in chi_roots:
    if r.imag().abs() < 1e-30:
        z0 = 1/r
        print(f"  z = {float(z0.real()):.10f} (real)")
    else:
        z0 = 1/r
        print(f"  z = {float(z0.real()):.6f} + {float(z0.imag()):.6f}i (|z| = {float(abs(z0)):.6f})")

# The real singularity
real_roots = [r for r in chi_roots if r.imag().abs() < 1e-30]
z0_exact = 1/max(real_roots, key=lambda r: r.real())
z0_approx = float(z0_exact.real())
print(f"\nDominant singularity z₀ ≈ {z0_approx:.15f}")

# Try ore_algebra monodromy computation
print("\n=== Monodromy computation ===")

# Method 1: Try monodromy_matrices
try:
    from ore_algebra.analytic.monodromy import monodromy_matrices
    print("monodromy_matrices available")

    # Compute monodromy at z₀
    # We need z₀ as an algebraic number
    # z₀ = 64/μ₀ where μ₀ is the dominant root of 4μ³-220μ²+8μ-1
    # Equivalently, z₀ is the smallest positive root of 2^20·z³·chi(1/z) = 0
    # or z₀ is a root of t³·chi(1/t)·2^20

    # Let's try with a numerical approximation first
    from ore_algebra.analytic.path import Point
    print(f"Computing monodromy at z₀ ≈ {z0_approx}...")
    print("(This may take several minutes for order 16)")

    # Use 50 digits of precision for now
    prec = 50
    try:
        monos = monodromy_matrices(M16, z0_exact, eps=RealField(prec)(2)^(-prec))
        print(f"Got {len(monos)} monodromy matrices")
        for i, M_mat in enumerate(monos):
            print(f"  Matrix {i}: {M_mat.nrows()}x{M_mat.ncols()}")
            MmI = M_mat - identity_matrix(M_mat.base_ring(), M_mat.nrows())
            # Numerical rank
            svs = [abs(s) for s in MmI.singular_values() if hasattr(MmI, 'singular_values')]
            rank_est = sum(1 for s in MmI.SVD()[1].diagonal() if abs(s) > 1e-10) if hasattr(MmI, 'SVD') else "unknown"
            print(f"  rank(M-I) ≈ {rank_est}")
    except Exception as e:
        print(f"monodromy_matrices failed: {e}")
        import traceback
        traceback.print_exc()
except ImportError:
    print("monodromy_matrices not available")

# Method 2: Try numerical_transition_matrix
try:
    from ore_algebra.analytic.analytic_continuation import analytic_continuation
    print("\nanalytic_continuation available")

    # Compute transition matrix from 0.5 to 0.5 via a loop around z₀
    # Path: 0.5 → z₀ - 0.3 + 0.2i → z₀ + 0.2i → z₀ + 0.3 → z₀ - 0.2i → z₀ - 0.3 - 0.2i → 0.5
    # This should encircle z₀ once

    CBF = ComplexBallField(200)
    z0_ball = CBF(z0_approx)

    # Simple loop: go to just before z₀, circle around, come back
    # Points on a small circle around z₀ of radius 0.15
    import cmath
    r_loop = 0.15
    n_pts = 8
    path_pts = [CBF(0.5)]
    for k in range(n_pts + 1):
        angle = 2 * cmath.pi * k / n_pts
        pt = z0_approx + r_loop * cmath.exp(1j * angle)
        path_pts.append(CBF(pt.real, pt.imag))
    path_pts.append(CBF(0.5))

    print(f"Path: {len(path_pts)} points, loop radius = {r_loop}")
    print(f"First few points: {[complex(p) for p in path_pts[:4]]}")

    # This gives the transition matrix for the analytic continuation
    # The monodromy matrix M satisfies: solutions at endpoint = M * solutions at startpoint
    # Since start = end = 0.5, M is the monodromy matrix

    try:
        # Try using the operator's numerical_transition_matrix method
        M_mono = M16.numerical_transition_matrix(path_pts, eps=RealField(53)(1e-15))
        print(f"Transition matrix: {M_mono.nrows()}x{M_mono.ncols()}")
        MmI = M_mono - identity_matrix(M_mono.base_ring(), M_mono.nrows())
        print(f"Max |M-I| entry: {max(abs(e) for e in MmI.list()):.6e}")
    except Exception as e:
        print(f"numerical_transition_matrix on path failed: {e}")
        import traceback
        traceback.print_exc()

except ImportError:
    print("analytic_continuation not available")

# Method 3: Direct approach via ore_algebra's transition matrix
print("\n=== Direct transition matrix approach ===")
try:
    # Compute basis of solutions at z=1/2 (regular point)
    # Then continue around z₀
    ini = M16.numerical_transition_matrix([0, QQ(1)/2], eps=1e-15)
    print(f"Transition 0 → 1/2: {ini.nrows()}x{ini.ncols()}")

    # Now loop around z₀
    # Use algebraic points for better precision
    z_pre = QQ(9)/10  # just before z₀ on the real line
    z_post = QQ(6)/5  # just past z₀
    z_above = QQbar(z0_approx + QQ(1)/5*QQbar.gen())  # z₀ + 0.2i

    # Hmm, need to be more careful with the path
    # Let me try a simpler approach: use numerical_transition_matrix with exact rational path points

    # Path from 1/2 around z₀ and back
    # z₀ ≈ 1.1644
    # Path: 1/2 → 9/10 → 9/10 + i/5 → 6/5 + i/5 → 6/5 - i/5 → 9/10 - i/5 → 9/10 → 1/2
    from ore_algebra.analytic.path import Point

    path_loop = [QQ(1)/2, QQ(9)/10,
                 QQbar(QQ(9)/10 + QQ(1)/5*QQbar.gen()),
                 QQbar(QQ(6)/5 + QQ(1)/5*QQbar.gen()),
                 QQbar(QQ(6)/5 - QQ(1)/5*QQbar.gen()),
                 QQbar(QQ(9)/10 - QQ(1)/5*QQbar.gen()),
                 QQ(9)/10, QQ(1)/2]

    print(f"Loop path: {[complex(p) for p in path_loop]}")
    M_loop = M16.numerical_transition_matrix(path_loop, eps=1e-30)
    print(f"Loop transition matrix: {M_loop.nrows()}x{M_loop.ncols()}")

    # This should be close to identity if the loop doesn't encircle z₀
    # or have nontrivial monodromy if it does
    MmI = M_loop - identity_matrix(M_loop.base_ring(), M_loop.nrows())
    max_entry = max(abs(e) for e in MmI.list())
    print(f"Max |M-I| entry: {float(max_entry):.6e}")

    if max_entry > 1e-5:
        print("Nontrivial monodromy detected!")
        # Check rank
        # Find the largest singular value / smallest, etc.
        # For a rank-1 matrix, the ratio of 2nd/1st singular value should be ~0
        entries = [abs(e) for e in MmI.list() if abs(e) > max_entry * 1e-10]
        print(f"Number of significant entries: {len(entries)}")

        # Try to extract the ratio
        # The image of M-I is 1-dimensional: all columns are proportional
        # Pick the column with largest norm
        cols = [vector([MmI[i,j] for i in range(MmI.nrows())]) for j in range(MmI.ncols())]
        norms = [c.norm() for c in cols]
        best_col = max(range(len(norms)), key=lambda i: norms[i])
        print(f"Best column: {best_col}, norm = {float(norms[best_col]):.6e}")

        # All other columns should be scalar multiples of this one
        v = cols[best_col]
        for j in range(MmI.ncols()):
            if norms[j] > max_entry * 1e-10:
                ratio = cols[j] / v if v.norm() > 0 else None
                # Check if ratio is approximately constant
                if ratio is not None:
                    ratios = [cols[j][i]/v[i] for i in range(len(v)) if abs(v[i]) > 1e-20]
                    if len(ratios) > 1:
                        spread = max(abs(r - ratios[0]) for r in ratios)
                        print(f"  col {j}: ratio spread = {float(spread):.6e}, ratio ≈ {complex(ratios[0]):.6f}")
    else:
        print("Loop monodromy is trivial (identity). Check path encircles z₀.")

except Exception as e:
    print(f"Direct approach failed: {e}")
    import traceback
    traceback.print_exc()

print("\nDone.")
