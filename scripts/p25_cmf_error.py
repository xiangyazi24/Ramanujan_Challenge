#!/usr/bin/env python3
"""P2.5: Correctly compute CMF error sequences and test kernel candidates.

Structure: A·M_N gives 2×3 matrix with rows P_{N,j} and Q_{N,j}.
Error: ê_{N,j} = G·Q_{N,j} - P_{N,j}
Test: does ê_{N,j} = ∫₀¹ (-log t)/(1+t²) · K_N(t²) dt for some known kernel K_N?
"""
import mpmath as mp
mp.mp.dps = 100

G = mp.catalan

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

def mat_vec_right(row, M):
    """row · M (1×3 × 3×3 = 1×3)"""
    return [sum(row[j] * M[j][i] for j in range(3)) for i in range(3)]

# Initial rows
p0 = [mp.mpf(30921), mp.mpf(-32972), mp.mpf(8240)]
q0 = [mp.mpf(33750), mp.mpf(-36000), mp.mpf(9000)]

NMAX = 25
p_rows = [p0[:]]
q_rows = [q0[:]]

for n in range(NMAX):
    Mn = M_entries(n)
    p_rows.append(mat_vec_right(p_rows[-1], Mn))
    q_rows.append(mat_vec_right(q_rows[-1], Mn))

# Verify convergence P_{N,j}/Q_{N,j} → G
print("=== Convergence check: P_{N,j}/Q_{N,j} → G ===")
print(f"G = {mp.nstr(G, 30)}")
for N in [0, 1, 2, 3, 5, 10, 15, 20, 25]:
    ratios = [p_rows[N][j] / q_rows[N][j] if q_rows[N][j] != 0 else 0 for j in range(3)]
    diffs = [abs(ratios[j] - G) for j in range(3)]
    digits = [-mp.log10(d) if d > 0 else 999 for d in diffs]
    print(f"  N={N:3d}: digits = [{mp.nstr(digits[0],5)}, {mp.nstr(digits[1],5)}, {mp.nstr(digits[2],5)}]")

# Compute errors ê_{N,j} = G·Q_{N,j} - P_{N,j}
print("\n=== Error sequences ê_{N,j} = G·Q_{N,j} - P_{N,j} ===")
errors = []
for N in range(NMAX + 1):
    e = [G * q_rows[N][j] - p_rows[N][j] for j in range(3)]
    errors.append(e)

for N in range(min(16, NMAX + 1)):
    print(f"  N={N:3d}: [{mp.nstr(errors[N][0], 12)}, {mp.nstr(errors[N][1], 12)}, {mp.nstr(errors[N][2], 12)}]")

# Error ratios (should approach c₀/c₊ = ρ ≈ 0.0294 for exponential part)
rho = mp.mpf(17) - 12*mp.sqrt(2)
print(f"\nρ = 17 - 12√2 = {mp.nstr(rho, 15)}")
print(f"\nError ratios ê_{'{N+1}'}/ê_N (should → ρ for dominant exponential mode):")
for N in range(1, min(16, NMAX)):
    if errors[N][0] != 0:
        r = errors[N+1][0] / errors[N][0]
        print(f"  N={N:3d}: ratio = {mp.nstr(r, 15)}")

# Normalize by H_N to get ê̂
def pochhammer(a, n):
    r = mp.mpf(1)
    for k in range(n):
        r *= (a + k)
    return r

def H_n(n):
    return (mp.mpf(-16))**n * pochhammer(2, n)**2 * pochhammer(3, n)**2 * \
           pochhammer(mp.mpf(5)/2, n) * pochhammer(mp.mpf(7)/2, n)**2

print(f"\n=== H_n-normalized errors ê̂_N = ê_N / H_N ===")
e_hat = []
for N in range(min(16, NMAX + 1)):
    Hn = H_n(N)
    eh = [errors[N][j] / Hn for j in range(3)]
    e_hat.append(eh)
    if N < 12:
        print(f"  N={N:3d}: [{mp.nstr(eh[0], 12)}, {mp.nstr(eh[1], 12)}, {mp.nstr(eh[2], 12)}]")

# Check: ê̂_N should decay like n^{-3} × ρ^N (from formal index -3)
print(f"\n=== Decay analysis of ê̂_N ===")
print("N^3 × ê̂_N (should stabilize for the neutral mode component):")
for N in range(3, min(16, len(e_hat))):
    scaled = [N**3 * e_hat[N][j] for j in range(3)]
    print(f"  N={N:3d}: [{mp.nstr(scaled[0], 12)}, {mp.nstr(scaled[1], 12)}, {mp.nstr(scaled[2], 12)}]")

# R₀(X) polynomial
# R_0(X) = q0 - (1+X)·p0 = (q0 - p0) - X·p0
# R_0(X) = [33750-30921, -36000+32972, 9000-8240] + X·[-30921, 32972, -8240]
# = [2829, -3028, 760] + X·[-30921, 32972, -8240]
r0_const = [q0[j] - p0[j] for j in range(3)]
r0_coeff = [-p0[j] for j in range(3)]
print(f"\nR_0(X) constant term: {[mp.nstr(c, 8) for c in r0_const]}")
print(f"R_0(X) X-coefficient: {[mp.nstr(c, 8) for c in r0_coeff]}")

# R_{N,j}(X) = (R_0(X) · M_N)_j
# At X = t², the integral ∫₀¹ (-log t)/(1+t²) R_{N,j}(t²) dt should equal ê_{N,j}
# This is the moment lift identity (already proved as Theorem 2, CIRCULAR)

# Let's verify this identity numerically
print(f"\n=== Verify moment lift identity (should be exact) ===")
for N in range(5):
    # R_N(X) = [r0_const + X·r0_coeff] · M(0)·...·M(N-1)
    # = (R_N_const + X · R_N_linear)
    r_const = r0_const[:]
    r_linear = r0_coeff[:]
    for n in range(N):
        Mn = M_entries(n)
        r_const = mat_vec_right(r_const, Mn)
        r_linear = mat_vec_right(r_linear, Mn)

    for j in range(3):
        A_Nj = r_const[j]
        B_Nj = r_linear[j]
        # ∫₀¹ (-log t)/(1+t²) (A + B t²) dt = A·G + B·(1-G)
        integral_val = A_Nj * G + B_Nj * (1 - G)
        diff = integral_val - errors[N][j]
        print(f"  N={N}, j={j}: integral = {mp.nstr(integral_val, 15)}, error = {mp.nstr(errors[N][j], 15)}, diff = {mp.nstr(diff, 6)}")

# Now test: does ê_{N,j} match the Catalan moment of ANY independently defined polynomial?
# Key idea from Q4873: use the Jacobi derivative P_n^{(0,ε)} differentiated at ε=0

# But the CMF kernel R_{N,j} is LINEAR in X (only degree 1).
# A Legendre polynomial of degree n has degree n.
# So the kernel can't literally be a Legendre polynomial for n > 1.

# However, the CMF might FACTOR through a Legendre system:
# R_{N,j}(X) might be a WIDTH-2 linear combination of Legendre evaluations.

# Actually, the problem is more subtle. The error ê_{N,j} = G·Q_{N,j} - P_{N,j}
# is just a NUMBER, not a function. The integral representation is:
# ê_{N,j} = ∫₀¹ (-log t)/(1+t²) · [A_{N,j} + B_{N,j}·t²] dt
# where A_{N,j}, B_{N,j} are the coefficients of R_{N,j}(X).

# The "kernel" is just a linear function! So the question is:
# Can we find an INDEPENDENT integral ∫₀¹ f(t) K_N(t) dt = ê_{N,j}
# where K_N is NOT defined via the CMF?

# For the Apéry case: ê_n = ∫₀¹∫₀¹ x^n(1-x)^n y^n(1-y)^n / (1-xy)^{n+1} dx dy
# This integral is INDEPENDENTLY defined (no CMF matrices), and
# the integrand's N-dependent part is explicit.

# For P2.5: the analogous construction would be a Beukers-type integral
# with explicit N-dependence. From Q4873, the Jacobi derivative gives one.

# Let me check: does D_n² (Delannoy squared) satisfy the SAME recurrence as Q_{N,j}?
print(f"\n=== Test: do Delannoy-squared values satisfy the CMF recurrence? ===")

def delannoy(n):
    """Central Delannoy number D_n = P_n(3)"""
    if n == 0: return mp.mpf(1)
    if n == 1: return mp.mpf(3)
    d_prev2, d_prev1 = mp.mpf(1), mp.mpf(3)
    for m in range(1, n):
        d_new = (3*(2*m+1)*d_prev1 - m*d_prev2) / (m+1)
        d_prev2 = d_prev1
        d_prev1 = d_new
    return d_prev1

# Compute D_n² for n=0..15
D_sq = [delannoy(n)**2 for n in range(16)]
print(f"D_n²: {[mp.nstr(d, 8) for d in D_sq[:8]]}")

# The Q_{N,j} (column j) should satisfy the ORDER-3 scalar recurrence
# c₃(N)·u_{N+1} + c₂(N)·u_N + c₁(N)·u_{N-1} + c₀(N)·u_{N-2} = 0

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

# Test if Q_{N,j} satisfy the recurrence
print(f"\n=== Verify Q_{'{N,j}'} satisfy the scalar recurrence ===")
for j in range(3):
    print(f"  Column j={j}:")
    for N in range(3, min(10, NMAX - 1)):
        res = (c_poly(3, N) * q_rows[N+1][j] + c_poly(2, N) * q_rows[N][j]
             + c_poly(1, N) * q_rows[N-1][j] + c_poly(0, N) * q_rows[N-2][j])
        # Relative residual
        scale = abs(c_poly(3, N) * q_rows[N+1][j])
        if scale > 0:
            rel = abs(res) / scale
            if N < 6:
                print(f"    N={N}: rel residual = {mp.nstr(rel, 6)}")

# Test if D_n² satisfy the SAME recurrence
print(f"\n=== Does D_n² satisfy the CMF recurrence? ===")
for N in range(3, 12):
    res = (c_poly(3, N) * D_sq[N+1] + c_poly(2, N) * D_sq[N]
         + c_poly(1, N) * D_sq[N-1] + c_poly(0, N) * D_sq[N-2])
    print(f"  N={N}: residual = {mp.nstr(res, 15)}")

# Q4873 approach: compute the JACOBI derivative's Catalan moment
# and compare with ê_{N,j}
print(f"\n=== Jacobi derivative kernel vs CMF error ===")

# The Jacobi derivative at z=3, ε=0:
# ∂_ε P_n^{(0,ε)}(3)|_{ε=0}
def jacobi_deriv_at_3(n):
    val = mp.mpf(0)
    for k in range(n + 1):
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

# The "Jacobi error" = G·D_n - ∂_ε P_n^{(0,ε)}(3)|_{ε=0} ?
# Actually, the Jacobi derivative gives ∂_ε [P_n^{(0,ε)}(3) / n!] or similar.
# The Catalan moment functional C[t^k] = ∫₀¹ (-log t)/(1+t²) t^{2k} dt
# = Σ_{j≥0} (-1)^j / (2k+2j+1)²

# G·P_n(3) - S_n where S_n = Σ_k F(n,k) C_k (the partial Catalan sums)
# But this is the Delannoy-squared carrier from the proof.

# For the CMF error, ê_{N,j} is defined by the matrix product.
# Let me check: is ê_{N,1} proportional to ∫ (-log t)/(1+t²) · [some simple function of t] dt?

# Since R_{N,j}(t²) = A_{N,j} + B_{N,j} t², the integral is A·G + B·(1-G).
# The error ê_{N,j} = (A_{N,j} - B_{N,j})G + B_{N,j} = Q_{N,j}·G - (-B_{N,j}) = Q·G - P.
# This is trivially exact but circular.

# The NON-CIRCULAR approach needs an integral where the N-dependence
# is in a SEPARATE, independently defined integrand.

# From the Delannoy-squared carrier:
# ê_n^{Del} = Σ_k F(n,k) C_k  where F(n,k) = 2^k C(2k,k) C(n,k) C(n+k,k)
# and C_k = Σ_{j<k} (-1)^j/(2j+1)²

# This gives P_n^{Del} - G·D_n² = ∫₀¹ (-log t)/(1+t²) · f_n(t) dt
# with an explicit f_n. BUT the Delannoy squared doesn't satisfy
# the CMF recurrence (as we just tested), so this is a DIFFERENT
# system with slower convergence (σ ≈ 0.882 vs ρ ≈ 0.029).

# The REAL question for P2.5 is: can we find an independent integral
# for the FAST CMF sequences (convergence rate ρ ≈ 0.029)?

# From Q4872: the natural parent is a ₃F₂ CMF. The trajectory matrix
# of the ₃F₂ parent should match M(n) after gauge.

# Let me compute: what is ê_{N,1} / (D_N² · ρ^N)?
# If the error has the form ê_{N,j} ~ C · D_N² · ρ^N · something,
# this ratio should stabilize.

print(f"\n=== Error normalized by D_N² × ρ^N ===")
for N in range(1, min(16, NMAX + 1)):
    Dn_sq = delannoy(N)**2
    rho_N = rho**N
    for j in range(3):
        if Dn_sq * rho_N != 0:
            ratio = errors[N][j] / (Dn_sq * rho_N)
            if j == 0:
                print(f"  N={N:3d}: ratio_j0 = {mp.nstr(ratio, 15)}")

print("\nDone.")
