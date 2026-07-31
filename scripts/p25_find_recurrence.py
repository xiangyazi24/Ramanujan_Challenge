#!/usr/bin/env python3
"""
P2.5: Find the scalar recurrence for Q_{N,1} with correct degree allocation.
The Newton polygon with constant slope -7 gives:
  deg(c₃)=7, deg(c₂)=14, deg(c₁)=21, deg(c₀)=28

Then check if f(z) = Σ Q_n z^n satisfies the pullback of the integrated-K ODE.
"""
from fractions import Fraction
import sys

def M_entries(n):
    n = Fraction(n)
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

NTERMS = 200
print(f"Computing {NTERMS} terms of Q_{{N,1}}...", flush=True)
q_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
Q = [q_row[0]]
for N in range(NTERMS - 1):
    M = M_entries(N)
    q_new = [sum(q_row[i]*M[i][j] for i in range(3)) for j in range(3)]
    q_row = q_new
    Q.append(q_row[0])
print(f"  Done. |Q[{NTERMS-1}]| has {len(str(abs(Q[NTERMS-1].numerator)))} digits", flush=True)

# Try finding recurrence with increasing degree ranges
# Newton polygon slope -7: deg(c_k) for shift k should be 28 - 7k
# But the UNNORMALIZED Q_n might have different pattern.
# Try: for each max total degree D, search for order-3 recurrence
# with c₃ deg d₃, c₂ deg d₂, c₁ deg d₁, c₀ deg d₀

def find_recurrence(Q, degs, verbose=True):
    """Find recurrence with given degree pattern"""
    d3, d2, d1, d0 = degs
    total_unk = (d3+1) + (d2+1) + (d1+1) + (d0+1)
    neqs = total_unk + 5
    
    if len(Q) < neqs + 3:
        return None
    
    rows = []
    for n in range(neqs):
        row = []
        for j in range(d3+1):
            row.append(Fraction(n)**j * Q[n+3])
        for j in range(d2+1):
            row.append(Fraction(n)**j * Q[n+2])
        for j in range(d1+1):
            row.append(Fraction(n)**j * Q[n+1])
        for j in range(d0+1):
            row.append(Fraction(n)**j * Q[n])
        rows.append(row)
    
    m = len(rows)
    nc = total_unk
    mat = [list(row) for row in rows]
    
    pivot_cols = []
    for col in range(nc):
        found = -1
        for r in range(len(pivot_cols), m):
            if mat[r][col] != 0:
                found = r
                break
        if found == -1:
            continue
        mat[len(pivot_cols)], mat[found] = mat[found], mat[len(pivot_cols)]
        pivot_row = len(pivot_cols)
        pivot_cols.append(col)
        pivot_val = mat[pivot_row][col]
        for r in range(m):
            if r != pivot_row and mat[r][col] != 0:
                factor = mat[r][col] / pivot_val
                for c2 in range(nc):
                    mat[r][c2] -= factor * mat[pivot_row][c2]
    
    rank = len(pivot_cols)
    nullity = nc - rank
    
    if verbose:
        print(f"  degs=({d3},{d2},{d1},{d0}), unk={total_unk}, rank={rank}, null={nullity}")
    
    if nullity != 1:
        return None
    
    # Extract kernel
    free_col = list(set(range(nc)) - set(pivot_cols))[0]
    kernel = [Fraction(0)] * nc
    kernel[free_col] = Fraction(1)
    for i in range(rank-1, -1, -1):
        pc = pivot_cols[i]
        val = sum(mat[i][c2] * kernel[c2] for c2 in range(nc) if c2 != pc)
        kernel[pc] = -val / mat[i][pc]
    
    # Parse coefficients
    c3 = kernel[:d3+1]
    c2 = kernel[d3+1:d3+1+d2+1]
    c1 = kernel[d3+1+d2+1:d3+1+d2+1+d1+1]
    c0 = kernel[d3+1+d2+1+d1+1:]
    
    # Verify
    def peval(coeffs, n):
        return sum(c * Fraction(n)**k for k, c in enumerate(coeffs))
    
    for n in range(neqs, min(neqs+20, len(Q)-3)):
        val = peval(c3,n)*Q[n+3] + peval(c2,n)*Q[n+2] + peval(c1,n)*Q[n+1] + peval(c0,n)*Q[n]
        if val != 0:
            if verbose:
                print(f"    VERIFY FAIL at n={n}")
            return None
    
    if verbose:
        print(f"    VERIFIED for 20 extra equations!")
    return (c3, c2, c1, c0)

print("\nSearching for scalar recurrence...", flush=True)

# First try the correct Newton polygon: (7,14,21,28)
result = find_recurrence(Q, (7, 14, 21, 28))
if result:
    print("Found with Newton polygon (7,14,21,28)!")
else:
    # Try equal degrees
    for d in range(5, 30):
        result = find_recurrence(Q, (d, d, d, d), verbose=False)
        if result:
            print(f"Found with equal degree d={d}!")
            break
    if not result:
        # Try systematically: for total <= 100
        for total in range(20, 100, 4):
            # Allocate uniformly
            d = total // 4
            result = find_recurrence(Q, (d, d, d, d), verbose=False)
            if result:
                print(f"Found with equal degree d={d}, total={4*(d+1)}")
                break
            # Try Newton-like: d₃=d, d₂=d+7, d₁=d+14, d₀=d+21
            for base in range(0, 10):
                degs = (base, base+7, base+14, base+21)
                if sum(d+1 for d in degs) > len(Q) - 8:
                    continue
                result = find_recurrence(Q, degs, verbose=False)
                if result:
                    print(f"Found with Newton-shifted base={base}: degs={degs}")
                    break
            if result:
                break

if result is None:
    print("No recurrence found!")
    sys.exit(1)

c3, c2, c1, c0 = result

# Poincaré analysis
def actual_deg(coeffs):
    d = len(coeffs) - 1
    while d > 0 and coeffs[d] == 0:
        d -= 1
    return d

d3 = actual_deg(c3)
d2 = actual_deg(c2)
d1 = actual_deg(c1)
d0 = actual_deg(c0)
print(f"\nActual degrees: c₃={d3}, c₂={d2}, c₁={d1}, c₀={d0}")
print(f"Degree pattern: ({d3},{d2},{d1},{d0})")

# Poincaré polynomial: leading coeff of each
lc3 = c3[d3]
lc2 = c2[d2]
lc1 = c1[d1]
lc0 = c0[d0]

print(f"\nLeading coefficients:")
print(f"  lc₃ (shift 3, deg {d3}) = {float(lc3):.10e}")
print(f"  lc₂ (shift 2, deg {d2}) = {float(lc2):.10e}")
print(f"  lc₁ (shift 1, deg {d1}) = {float(lc1):.10e}")
print(f"  lc₀ (shift 0, deg {d0}) = {float(lc0):.10e}")

# Poincaré polynomial: lc₃·ξ³ + lc₂·ξ² + lc₁·ξ + lc₀ = 0
a = lc2/lc3
b = lc1/lc3
c = lc0/lc3
print(f"\nPoincaré poly: ξ³ + ({float(a):.6f})ξ² + ({float(b):.6f})ξ + ({float(c):.6f})")

# Expected: (c+16)(c²+544c+256) = c³ + 560c² + 8960c + 4096
# With shifts based on which convention...
print(f"Expected (from proof): ξ³ + 560ξ² + 8960ξ + 4096")
print(f"Ratios to expected: {float(a)/560:.10f}, {float(b)/8960:.10f}, {float(c)/4096:.10f}")

# Now check ODE: does Σ Q_n z^n satisfy the pullback of k(1-k²)Y'''+...=0?
print("\n=== ODE check ===")
from mpmath import mp, mpf, sqrt as msqrt, catalan
mp.dps = 100

# Compute f, f', f'', f''' at several z-values by power series
def eval_gf(Q, z, deriv=0):
    """Evaluate Σ n^(deriv-factor) Q_n z^n or its derivatives"""
    val = mpf(0)
    for n in range(len(Q)):
        qn = mpf(int(Q[n].numerator)) / mpf(int(Q[n].denominator))
        if deriv == 0:
            val += qn * z**n
        elif deriv == 1:
            if n >= 1:
                val += qn * n * z**(n-1)
        elif deriv == 2:
            if n >= 2:
                val += qn * n*(n-1) * z**(n-2)
        elif deriv == 3:
            if n >= 3:
                val += qn * n*(n-1)*(n-2) * z**(n-3)
    return val

z_test = mpf('0.005')
f0 = eval_gf(Q, z_test, 0)
f1 = eval_gf(Q, z_test, 1)
f2 = eval_gf(Q, z_test, 2)
f3 = eval_gf(Q, z_test, 3)

# k(z) = 4√(2z)/(1-z)
k = 4*msqrt(2*z_test)/(1-z_test)
# dk/dz = 2√2(1+z)/(√z·(1-z)²)  
dk = 2*msqrt(2)*(1+z_test)/(msqrt(z_test)*(1-z_test)**2)

# Y = f(z), dY/dk = f'(z) / dk/dz
Y1 = f1/dk

# d²Y/dk² using chain rule properly
# dk/dz = k'(z), d²k/dz² = k''(z)
# Y' = f'/k', Y'' = (f'' - Y'·k'')/k'^2 = (f''·k' - f'·k'')/k'^3
h = z_test * mpf('1e-15')
dk_p = 2*msqrt(2)*(1+z_test+h)/(msqrt(z_test+h)*(1-z_test-h)**2)
dk_m = 2*msqrt(2)*(1+z_test-h)/(msqrt(z_test-h)*(1-z_test+h)**2)
d2k = (dk_p - dk_m)/(2*h)

Y2 = (f2 - Y1*d2k)/dk**2

# Y''' by finite difference of Y'' at nearby z
delta = z_test * mpf('1e-8')
# Compute Y2 at z+delta and z-delta
for label, z_shifted in [("z+δ", z_test+delta), ("z-δ", z_test-delta)]:
    f1s = eval_gf(Q, z_shifted, 1)
    f2s = eval_gf(Q, z_shifted, 2)
    ks = 4*msqrt(2*z_shifted)/(1-z_shifted)
    dks = 2*msqrt(2)*(1+z_shifted)/(msqrt(z_shifted)*(1-z_shifted)**2)
    Y1s = f1s/dks
    hs = z_shifted * mpf('1e-15')
    dksp = 2*msqrt(2)*(1+z_shifted+hs)/(msqrt(z_shifted+hs)*(1-z_shifted-hs)**2)
    dksm = 2*msqrt(2)*(1+z_shifted-hs)/(msqrt(z_shifted-hs)*(1-z_shifted+hs)**2)
    d2ks = (dksp - dksm)/(2*hs)
    Y2s = (f2s - Y1s*d2ks)/dks**2
    if label == "z+δ":
        Y2_plus = Y2s
        dk_at_plus = dks
    else:
        Y2_minus = Y2s
        dk_at_minus = dks

Y3 = (Y2_plus - Y2_minus)/(2*delta) / dk

# Check integrated-K ODE: k(1-k²)Y''' + (1-3k²)Y'' - kY' = 0
ode_val = k*(1-k**2)*Y3 + (1-3*k**2)*Y2 - k*Y1
print(f"At z={mp.nstr(z_test,4)}:")
print(f"  k = {mp.nstr(k,15)}")
print(f"  |ODE residual| = {mp.nstr(abs(ode_val), 10)}")
print(f"  |k·Y'| = {mp.nstr(abs(k*Y1), 10)}")
print(f"  relative = {mp.nstr(abs(ode_val/(k*Y1+1)), 10)}")

# Also check the Gauss/K ODE: k(1-k²)Y'' + (1-3k²)Y' - kY = 0
# (this would be for f' = Σ n Q_n z^{n-1}, not f itself)
# Actually, the integrated-K factors as L_K ∘ D. So Y' should satisfy L_K.
# Y' = f'(z)/k'(z) is the derivative of Y w.r.t. k.
# f' satisfies a certain ODE in z.
# Y' = f'/k' should satisfy k(1-k²)W'' + (1-3k²)W' - kW = 0
ode_gauss = k*(1-k**2)*Y2 + (1-3*k**2)*Y1 - k*f0  # wait, need W = Y', W' = Y'', W itself is Y'
# W = Y', W' = Y'', and L_K(W) = k(1-k²)W'' + (1-3k²)W' - kW
# But W = Y' = dY/dk, W' = d²Y/dk², W'' = d³Y/dk³ (w.r.t. k)
# Hmm, the notation is confusing. Let me be careful.
# L_K = k(1-k²)D² + (1-3k²)D - k where D = d/dk
# L_K(W) = k(1-k²)W'' + (1-3k²)W' - kW where primes are d/dk
# With W = Y' (= dY/dk = Y1), W' = Y'' (= Y2), W'' = Y''' (= Y3):
ode_gauss = k*(1-k**2)*Y3 + (1-3*k**2)*Y2 - k*Y1  # this is L_K(Y') = L_K ∘ D (Y)

# But this is the SAME as the integrated-K ODE!
# L_int(Y) = L_K(D(Y)) = L_K(Y') = k(1-k²)Y''' + (1-3k²)Y'' - kY' = 0

print(f"\nNote: the integrated-K ODE IS L_K applied to Y', which is the same expression.")
print(f"If the ODE residual is small, f(z) satisfies the pullback, proving the module is M_int.")

# The relative residual tells us how close we are
# Finite difference errors limit precision. Let me also try at a different z.
for z_test2 in [mpf('0.001'), mpf('0.01'), mpf('0.02')]:
    f0 = eval_gf(Q, z_test2, 0)
    f1 = eval_gf(Q, z_test2, 1)
    f2 = eval_gf(Q, z_test2, 2)
    f3 = eval_gf(Q, z_test2, 3)
    k = 4*msqrt(2*z_test2)/(1-z_test2)
    dk = 2*msqrt(2)*(1+z_test2)/(msqrt(z_test2)*(1-z_test2)**2)
    h = z_test2 * mpf('1e-15')
    dk_p = 2*msqrt(2)*(1+z_test2+h)/(msqrt(z_test2+h)*(1-z_test2-h)**2)
    dk_m = 2*msqrt(2)*(1+z_test2-h)/(msqrt(z_test2-h)*(1-z_test2+h)**2)
    d2k = (dk_p - dk_m)/(2*h)
    Y1 = f1/dk
    Y2 = (f2 - Y1*d2k)/dk**2
    delta = z_test2 * mpf('1e-8')
    f1p = eval_gf(Q, z_test2+delta, 1)
    f2p = eval_gf(Q, z_test2+delta, 2)
    dkp = 2*msqrt(2)*(1+z_test2+delta)/(msqrt(z_test2+delta)*(1-z_test2-delta)**2)
    hp = (z_test2+delta)*mpf('1e-15')
    dkpp = 2*msqrt(2)*(1+z_test2+delta+hp)/(msqrt(z_test2+delta+hp)*(1-z_test2-delta-hp)**2)
    dkpm = 2*msqrt(2)*(1+z_test2+delta-hp)/(msqrt(z_test2+delta-hp)*(1-z_test2-delta+hp)**2)
    d2kp = (dkpp-dkpm)/(2*hp)
    Y1p = f1p/dkp
    Y2p = (f2p - Y1p*d2kp)/dkp**2
    
    f1m = eval_gf(Q, z_test2-delta, 1)
    f2m = eval_gf(Q, z_test2-delta, 2)
    dkm = 2*msqrt(2)*(1+z_test2-delta)/(msqrt(z_test2-delta)*(1-z_test2+delta)**2)
    hm = (z_test2-delta)*mpf('1e-15')
    dkmp = 2*msqrt(2)*(1+z_test2-delta+hm)/(msqrt(z_test2-delta+hm)*(1-z_test2+delta-hm)**2)
    dkmm = 2*msqrt(2)*(1+z_test2-delta-hm)/(msqrt(z_test2-delta-hm)*(1-z_test2+delta+hm)**2)
    d2km = (dkmp-dkmm)/(2*hm)
    Y1m = f1m/dkm
    Y2m = (f2m - Y1m*d2km)/dkm**2
    
    Y3 = (Y2p - Y2m)/(2*delta)/dk
    
    ode_val = k*(1-k**2)*Y3 + (1-3*k**2)*Y2 - k*Y1
    scale = abs(k*Y1) + abs(k*(1-k**2)*Y3) + 1
    print(f"  z={mp.nstr(z_test2,4)}: |ODE|/scale = {mp.nstr(abs(ode_val/scale), 6)}")

print("\nDone.")
