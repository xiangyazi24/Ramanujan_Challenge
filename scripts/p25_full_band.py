#!/usr/bin/env python3
"""P2.5: Full contiguous band search with R^J, K^J, P^J families.

From Q4886: search for R_n^cand = Σ a_r(n) R_{n+r}^J + Σ b_r(n) K_{n+r}^J + Σ c_r(n) P_{n+r}^J
that satisfies the CMF recurrence.
"""
from fractions import Fraction
import math

def rising_frac(a, k):
    result = Fraction(1)
    for i in range(k):
        result *= (a + i)
    return result

def odd_harmonic(m):
    return sum(Fraction(1, 2*r+1) for r in range(m))

def J_and_dJ(N):
    pref = rising_frac(Fraction(1,2), N) / Fraction(math.factorial(N))
    J = [Fraction(0)] * (N + 1)
    dJ = [Fraction(0)] * (N + 1)
    for k in range(N + 1):
        c = (Fraction(-1)**k * Fraction(math.comb(N, k)) * pref
             * rising_frac(Fraction(N) + Fraction(1,2), k)
             / rising_frac(Fraction(1,2), k))
        dl = (2 * odd_harmonic(N)
              + 2 * sum(Fraction(1, 2*N+2*r+1) for r in range(k))
              - 2 * odd_harmonic(k))
        J[k] = c
        dJ[k] = c * dl
    return J, dJ

def poly_eval(c, x):
    r = Fraction(0)
    xp = Fraction(1)
    for v in c:
        r += v * xp
        xp *= x
    return r

def poly_deriv(c):
    return [Fraction(k+1)*c[k+1] for k in range(len(c)-1)] if len(c)>1 else [Fraction(0)]

def poly_scale(c, s):
    return [v*s for v in c]

def poly_add(a, b):
    n = max(len(a), len(b))
    r = [Fraction(0)]*n
    for i in range(len(a)): r[i] += a[i]
    for i in range(len(b)): r[i] += b[i]
    return r

def poly_sub(a, b):
    return poly_add(a, poly_scale(b, Fraction(-1)))

def poly_mul(a, b):
    if not a or not b: return [Fraction(0)]
    r = [Fraction(0)]*(len(a)+len(b)-1)
    for i in range(len(a)):
        for j in range(len(b)):
            r[i+j] += a[i]*b[j]
    return r

def poly_divmod(num, den):
    num = list(num)
    dd = len(den)-1
    q = [Fraction(0)]*max(0, len(num)-dd)
    for i in range(len(q)-1,-1,-1):
        q[i] = num[i+dd]/den[dd]
        for j in range(dd+1):
            num[i+j] -= q[i]*den[j]
    return q, num[:dd] if dd>0 else []

def catalan_pair(k):
    q = Fraction((-1)**k)
    p = q * sum(Fraction((-1)**j,(2*j+1)**2) for j in range(k))
    return q, p

def moment_pair(coeffs):
    q = p = Fraction(0)
    for k, a in enumerate(coeffs):
        if a == 0: continue
        qk, pk = catalan_pair(k)
        q += a*qk; p += a*pk
    return q, p

NMAX = 40
print("Computing families...", flush=True)

# Precompute J, B, dJ for all N
J_cache = {}
B_cache = {}
kappa_cache = {}
for N in range(NMAX + 5):
    J_c, dJ_c = J_and_dJ(N)
    B_val = poly_eval(J_c, Fraction(-1))
    J_cache[N] = (J_c, dJ_c)
    B_cache[N] = B_val
    kappa_cache[N] = Fraction(4*N+1, 2)
    if N % 10 == 0:
        print(f"  J_{N} done", flush=True)

# Family R: R_n^J
def R_pair(N):
    J_c, dJ_c = J_cache[N]
    B = B_cache[N]; kappa = kappa_cache[N]
    dB = poly_eval(dJ_c, Fraction(-1))
    num = poly_sub(poly_scale(dJ_c, B), poly_scale(J_c, dB))
    C, _ = poly_divmod(num, [Fraction(1), Fraction(1)])
    BJ = poly_scale(J_c, B)
    Cp = poly_deriv(C)
    twoXCp = [Fraction(0)] + poly_scale(Cp, Fraction(2))
    inner = poly_add(twoXCp, C)
    hi = poly_scale(inner, Fraction(1,2))
    xp1h = poly_mul(hi, [Fraction(1), Fraction(1)])
    R = poly_scale(poly_sub(BJ, xp1h), kappa)
    return moment_pair(R)

# Family P: P_n^J = κ_n B_n J_n
def P_pair(N):
    J_c, _ = J_cache[N]
    B = B_cache[N]; kappa = kappa_cache[N]
    P_poly = poly_scale(J_c, kappa * B)
    return moment_pair(P_poly)

# Family K: K_n^J = Σ_{k=0}^n κ_k B_k J_k
# Moment pair computed incrementally
K_q_cum = Fraction(0)
K_p_cum = Fraction(0)
K_pairs = {}
for N in range(NMAX + 5):
    J_c, _ = J_cache[N]
    kB = kappa_cache[N] * B_cache[N]
    # moment of κ_N B_N J_N(X)
    q_raw, p_raw = moment_pair(poly_scale(J_c, kB))
    K_q_cum += q_raw
    K_p_cum += p_raw
    K_pairs[N] = (K_q_cum, K_p_cum)

print("Computing family pairs...", flush=True)
R_pairs = {}
P_pairs = {}
for N in range(NMAX + 5):
    R_pairs[N] = R_pair(N)
    P_pairs[N] = P_pair(N)

print("  Done.", flush=True)

# CMF pairs
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

def delta_H(n):
    n = Fraction(n)
    return Fraction(-2)*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

def MH_at(n):
    M = M_entries(n)
    d = delta_H(n)
    return [[M[i][j]/d for j in range(3)] for i in range(3)]

print("Computing CMF pairs...", flush=True)
p_row = [Fraction(30921), Fraction(-32972), Fraction(8240)]
q_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
cmf_q = []; cmf_p = []
for N in range(NMAX):
    cmf_q.append(q_row[0]); cmf_p.append(p_row[0])
    MH = MH_at(N)
    p_row = [sum(p_row[i]*MH[i][j] for i in range(3)) for j in range(3)]
    q_row = [sum(q_row[i]*MH[i][j] for i in range(3)) for j in range(3)]

# Band search with all three families
families = {
    'R': lambda n: R_pairs[n],
    'P': lambda n: P_pairs[n],
    'K': lambda n: K_pairs[n],
}

print("\n=== Full band search ===", flush=True)

for width in range(3):
    for deg in range(4):
        labels = []
        for fname in ['R', 'K', 'P']:
            for r in range(width + 1):
                for d in range(deg + 1):
                    labels.append((fname, r, d))

        n_unknowns = len(labels)
        n_train = n_unknowns + 5
        n_holdout = 3

        if n_train + n_holdout + width >= NMAX:
            continue

        # Build system (both q and p)
        A_rows = []; b_vec = []
        for n in range(n_train + n_holdout):
            row_q = []; row_p = []
            for fname, r, d in labels:
                pair = families[fname](n + r)
                row_q.append(Fraction(n)**d * pair[0])
                row_p.append(Fraction(n)**d * pair[1])
            A_rows.append(row_q); b_vec.append(cmf_q[n])
            A_rows.append(row_p); b_vec.append(cmf_p[n])

        # Solve
        n_eq = 2 * n_train
        aug = [list(A_rows[i]) + [b_vec[i]] for i in range(n_eq)]

        pivot_cols = []; row_idx = 0
        for col in range(n_unknowns):
            found = -1
            for rr in range(row_idx, n_eq):
                if aug[rr][col] != 0:
                    found = rr; break
            if found == -1: continue
            aug[row_idx], aug[found] = aug[found], aug[row_idx]
            piv = aug[row_idx][col]
            for j2 in range(n_unknowns+1): aug[row_idx][j2] /= piv
            for rr in range(n_eq):
                if rr == row_idx: continue
                if aug[rr][col] == 0: continue
                f = aug[rr][col]
                for j2 in range(n_unknowns+1): aug[rr][j2] -= f*aug[row_idx][j2]
            pivot_cols.append(col); row_idx += 1

        rank = len(pivot_cols)
        if rank < n_unknowns:
            if deg == 0:
                print(f"  w={width}, d={deg}: rank={rank}<{n_unknowns}", flush=True)
            continue

        consistent = True
        for i in range(rank, n_eq):
            if aug[i][n_unknowns] != 0:
                consistent = False; break

        if not consistent:
            if deg <= 1 or deg % 2 == 0:
                print(f"  w={width}, d={deg}: INCONSISTENT", flush=True)
            continue

        x = [Fraction(0)] * n_unknowns
        for pi, pc in enumerate(pivot_cols):
            x[pc] = aug[pi][n_unknowns]

        # Holdout
        ok = True
        for n in range(n_train, n_train + n_holdout):
            pq = pp = Fraction(0)
            for idx, (fname, r, d) in enumerate(labels):
                pair = families[fname](n + r)
                pq += x[idx] * Fraction(n)**d * pair[0]
                pp += x[idx] * Fraction(n)**d * pair[1]
            if pq != cmf_q[n] or pp != cmf_p[n]:
                ok = False; break

        if ok:
            print(f"\n*** MATCH: width={width}, deg={deg} ***", flush=True)
            for idx, (fname, r, d) in enumerate(labels):
                if x[idx] != 0:
                    print(f"  {fname}_{r}, n^{d}: {x[idx]}")

            # Full verify
            all_ok = True
            for n in range(n_train+n_holdout, NMAX-width):
                pq = Fraction(0)
                for idx, (fname, r, d) in enumerate(labels):
                    pair = families[fname](n+r)
                    pq += x[idx] * Fraction(n)**d * pair[0]
                if pq != cmf_q[n]:
                    all_ok = False
                    print(f"  FAIL at n={n}")
                    break
            if all_ok:
                print(f"  Verified ALL n=0..{NMAX-width-1}")
            break
        else:
            if deg <= 1:
                print(f"  w={width}, d={deg}: holdout fails", flush=True)
    else:
        continue
    break

print("\nDone.")
