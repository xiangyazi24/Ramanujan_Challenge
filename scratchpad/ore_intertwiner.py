#!/usr/bin/env python3
"""
Numerical search for sequence-level Ore intertwiner.

Given:
  b_n = Cauchy convolution C*a (AESZ #209 twisted by formal power correction)
  q_n = P2.7 sequence

Search for T such that (T·b)_n = const · q_n, where:
  (T·b)_n = Σ_{j=0}^r Σ_{k=0}^d  t_{j,k} · n^k · b_{n-j}

Also try direct sequence transforms:
  q_n = Σ_{j=0}^r  p_j(n) · a_{n-j}    (direct from AESZ, no Cauchy twist)
"""
import mpmath
mpmath.mp.dps = 150

def binom(n, k):
    if k < 0 or k > n: return 0
    r = 1
    for i in range(k): r = r * (n - i) // (i + 1)
    return r

N = 40

# AESZ #209
print("Computing AESZ #209...")
a = [mpmath.mpf(0)] * N
for n in range(N):
    a[n] = sum(mpmath.mpf(binom(n,k)**2 * binom(n+k,n) * binom(n+2*k,n)) for k in range(n+1))

# P2.7
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

# ============================================================
# APPROACH 1: Direct search for q_n = Σ p_j(n) · a_{n-j}
# This is the "direct Ore intertwiner" without Cauchy twist
# ============================================================
print("\n" + "="*70)
print("APPROACH 1: q_n = Σ_{j=0}^r p_j(n) · a_{n-j}")
print("="*70)

for r in range(4):  # shift order 0,1,2,3
    for d in range(6):  # polynomial degree 0,1,...,5
        num_params = (r+1) * (d+1)
        # Need at least num_params equations
        n_start = r  # n must be >= r for a_{n-j} to exist
        n_end = min(n_start + num_params + 5, N-1)
        num_eqs = n_end - n_start

        if num_eqs < num_params:
            continue

        # Build overdetermined system
        M = mpmath.matrix(num_eqs, num_params)
        rhs = mpmath.matrix(num_eqs, 1)

        for idx, n in enumerate(range(n_start, n_end)):
            col = 0
            for j in range(r+1):
                if n-j < 0:
                    for k in range(d+1):
                        M[idx, col] = 0
                        col += 1
                else:
                    for k in range(d+1):
                        M[idx, col] = mpmath.mpf(n)**k * a[n-j]
                        col += 1
            rhs[idx, 0] = q[n]

        # Solve least squares: M^T M x = M^T rhs
        try:
            MTM = M.T * M
            MTr = M.T * rhs
            x = mpmath.lu_solve(MTM, MTr)

            # Check residual
            pred = M * x
            max_rel = 0
            for idx in range(num_eqs):
                if abs(rhs[idx,0]) > 1e-50:
                    rel = float(abs((pred[idx,0] - rhs[idx,0])/rhs[idx,0]))
                    max_rel = max(max_rel, rel)

            if max_rel < 1e-10:
                print(f"\n  *** HIT: r={r}, d={d} (max rel err = {max_rel:.2e}) ***")
                for j in range(r+1):
                    coeffs = [float(x[(d+1)*j+k, 0]) for k in range(d+1)]
                    print(f"    p_{j}(n) = {coeffs}")

                # Extended verification
                ok = True
                for n in range(n_end, min(N-1, n_end+10)):
                    val = mpmath.mpf(0)
                    for j in range(r+1):
                        if n-j >= 0:
                            for k in range(d+1):
                                val += x[(d+1)*j+k, 0] * mpmath.mpf(n)**k * a[n-j]
                    rel = float(abs((val - q[n])/q[n])) if abs(q[n]) > 1e-50 else float(abs(val))
                    if rel > 1e-10:
                        ok = False
                        break
                if ok:
                    print(f"    Extended check PASSED (n up to {min(N-1, n_end+10)-1})")
                else:
                    print(f"    Extended check FAILED")
        except:
            pass

# ============================================================
# APPROACH 2: q_n / h_n = Σ p_j(n) · (a_{n-j} / h_{n-j})
# where h_n is the Pochhammer gauge
# ============================================================
print("\n" + "="*70)
print("APPROACH 2: Pochhammer-gauged search")
print("="*70)

# h_n = 2^{-20n} · (3)_n^4 · (4)_n^6 / [(5/2)_n^4 · (7/2)_n^3 · (9/2)_n^3]
def pochhammer(a_val, n_val):
    """(a)_n = a(a+1)...(a+n-1)"""
    result = mpmath.mpf(1)
    for i in range(n_val):
        result *= (a_val + i)
    return result

h = [mpmath.mpf(0)] * N
for n in range(N):
    h[n] = (mpmath.mpf(2)**(-20*n)
            * pochhammer(mpmath.mpf(3), n)**4
            * pochhammer(mpmath.mpf(4), n)**6
            / (pochhammer(mpmath.mpf('5/2'), n)**4
               * pochhammer(mpmath.mpf('7/2'), n)**3
               * pochhammer(mpmath.mpf('9/2'), n)**3))

# Gauged sequences
a_g = [a[n] / h[n] if abs(h[n]) > 1e-200 else mpmath.mpf(0) for n in range(N)]
q_g = [q[n] / h[n] if abs(h[n]) > 1e-200 else mpmath.mpf(0) for n in range(N)]

print("  a_gauged ratios:")
for n in range(1, min(10, N)):
    if abs(a_g[n-1]) > 1e-50:
        r = float(a_g[n] / a_g[n-1])
        print(f"    a_g[{n}]/a_g[{n-1}] = {r:.10f}")

print("  q_gauged ratios:")
for n in range(1, min(10, N)):
    if abs(q_g[n-1]) > 1e-50:
        r = float(q_g[n] / q_g[n-1])
        print(f"    q_g[{n}]/q_g[{n-1}] = {r:.10f}")

# Ratio q_g / a_g
print("\n  Ratio q_gauged / a_gauged:")
for n in range(min(15, N)):
    if abs(a_g[n]) > 1e-50:
        r = float(q_g[n] / a_g[n])
        print(f"    n={n:2d}: {r:+.15e}")

# ============================================================
# APPROACH 3: Search for q_n in terms of SHIFTED a_n
# q_n = Σ α_j · a_{n+j} with RATIONAL α_j
# ============================================================
print("\n" + "="*70)
print("APPROACH 3: q_n = Σ α_j · a_{n+j} (shift search)")
print("="*70)

for r in range(1, 8):  # shifts -r..r
    shifts = list(range(-r, r+1))
    num_params = len(shifts)
    n_start = r
    n_end = min(n_start + num_params + 5, N - r)
    num_eqs = n_end - n_start

    if num_eqs < num_params:
        continue

    M = mpmath.matrix(num_eqs, num_params)
    rhs = mpmath.matrix(num_eqs, 1)

    for idx, n in enumerate(range(n_start, n_end)):
        for col, j in enumerate(shifts):
            M[idx, col] = a[n+j]
        rhs[idx, 0] = q[n]

    try:
        MTM = M.T * M
        MTr = M.T * rhs
        x = mpmath.lu_solve(MTM, MTr)

        pred = M * x
        max_rel = 0
        for idx in range(num_eqs):
            if abs(rhs[idx,0]) > 1e-50:
                rel = float(abs((pred[idx,0] - rhs[idx,0])/rhs[idx,0]))
                max_rel = max(max_rel, rel)

        if max_rel < 1e-10:
            print(f"  *** HIT: shifts {shifts} (max rel = {max_rel:.2e}) ***")
            for col, j in enumerate(shifts):
                print(f"    α[{j:+d}] = {float(x[col,0]):.15e}")
    except:
        pass

print("\n  No constant-coefficient shift formula found (expected).")

# ============================================================
# APPROACH 4: q_n = R(n) · a_n + S(n) · a_{n-1} + T(n) · a_{n-2}
# where R,S,T are RATIONAL functions of n
# ============================================================
print("\n" + "="*70)
print("APPROACH 4: polynomial-coefficient shift (r=2, d up to 8)")
print("="*70)

for d in range(9):
    num_params = 3 * (d+1)
    n_start = 2
    n_end = min(n_start + num_params + 3, N-1)
    num_eqs = n_end - n_start

    if num_eqs < num_params:
        continue

    M = mpmath.matrix(num_eqs, num_params)
    rhs = mpmath.matrix(num_eqs, 1)

    for idx, n in enumerate(range(n_start, n_end)):
        col = 0
        for j in range(3):
            for k in range(d+1):
                M[idx, col] = mpmath.mpf(n)**k * a[n-j]
                col += 1
        rhs[idx, 0] = q[n]

    try:
        MTM = M.T * M
        MTr = M.T * rhs
        x = mpmath.lu_solve(MTM, MTr)

        pred = M * x
        max_rel = 0
        for idx in range(num_eqs):
            if abs(rhs[idx,0]) > 1e-50:
                rel = float(abs((pred[idx,0] - rhs[idx,0])/rhs[idx,0]))
                max_rel = max(max_rel, rel)

        if max_rel < 1e-10:
            print(f"  *** HIT: d={d} (max rel = {max_rel:.2e}) ***")
            for j in range(3):
                coeffs = [float(x[(d+1)*j+k, 0]) for k in range(d+1)]
                print(f"    p_{j}(n) = {coeffs}")
        elif d <= 3:
            print(f"  d={d}: max rel = {max_rel:.2e}")
    except Exception as e:
        if d <= 3:
            print(f"  d={d}: {str(e)[:60]}")

print("\n  (If no hits, P2.7 is NOT a finite polynomial-coefficient")
print("   shift of AESZ #209 of order ≤2 and degree ≤8.)")
