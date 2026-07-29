#!/usr/bin/env python3
"""P2.5: Fixed band search — handle underdetermined systems.

Check consistency of underdetermined systems and extract solutions.
"""
from fractions import Fraction
import math

def rising_frac(a, k):
    result = Fraction(1)
    for i in range(k): result *= (a + i)
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
        J[k] = c; dJ[k] = c * dl
    return J, dJ

def poly_eval(c, x):
    r = Fraction(0); xp = Fraction(1)
    for v in c: r += v * xp; xp *= x
    return r

def poly_deriv(c):
    return [Fraction(k+1)*c[k+1] for k in range(len(c)-1)] if len(c)>1 else [Fraction(0)]

def poly_scale(c, s): return [v*s for v in c]
def poly_add(a, b):
    n = max(len(a), len(b)); r = [Fraction(0)]*n
    for i in range(len(a)): r[i] += a[i]
    for i in range(len(b)): r[i] += b[i]
    return r
def poly_sub(a, b): return poly_add(a, poly_scale(b, Fraction(-1)))
def poly_mul(a, b):
    if not a or not b: return [Fraction(0)]
    r = [Fraction(0)]*(len(a)+len(b)-1)
    for i in range(len(a)):
        for j in range(len(b)): r[i+j] += a[i]*b[j]
    return r
def poly_divmod(num, den):
    num = list(num); dd = len(den)-1
    q = [Fraction(0)]*max(0, len(num)-dd)
    for i in range(len(q)-1,-1,-1):
        q[i] = num[i+dd]/den[dd]
        for j in range(dd+1): num[i+j] -= q[i]*den[j]
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

NMAX = 45

print("Precomputing...", flush=True)
J_cache = {}; B_cache = {}; kappa_cache = {}
for N in range(NMAX + 5):
    J_c, dJ_c = J_and_dJ(N)
    B_cache[N] = poly_eval(J_c, Fraction(-1))
    kappa_cache[N] = Fraction(4*N+1, 2)
    J_cache[N] = (J_c, dJ_c)
    if N % 15 == 0: print(f"  J_{N}", flush=True)

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

def P_pair(N):
    J_c, _ = J_cache[N]
    P_poly = poly_scale(J_c, kappa_cache[N] * B_cache[N])
    return moment_pair(P_poly)

R_p = {}; P_p = {}; K_p = {}
Kq = Kp = Fraction(0)
for N in range(NMAX + 5):
    R_p[N] = R_pair(N)
    P_p[N] = P_pair(N)
    J_c, _ = J_cache[N]
    kB = kappa_cache[N] * B_cache[N]
    qr, pr = moment_pair(poly_scale(J_c, kB))
    Kq += qr; Kp += pr
    K_p[N] = (Kq, Kp)

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
    M = M_entries(n); d = delta_H(n)
    return [[M[i][j]/d for j in range(3)] for i in range(3)]

print("CMF pairs...", flush=True)
pr = [Fraction(30921), Fraction(-32972), Fraction(8240)]
qr = [Fraction(33750), Fraction(-36000), Fraction(9000)]
cq = []; cp = []
for N in range(NMAX):
    cq.append(qr[0]); cp.append(pr[0])
    MH = MH_at(N)
    pr = [sum(pr[i]*MH[i][j] for i in range(3)) for j in range(3)]
    qr = [sum(qr[i]*MH[i][j] for i in range(3)) for j in range(3)]

fams = {'R': lambda n: R_p[n], 'P': lambda n: P_p[n], 'K': lambda n: K_p[n]}

def solve_augmented(A_rows, b_vec, n_unknowns):
    """Return (rank, aug_rank, solution_or_None)."""
    m = len(A_rows)
    aug = [list(A_rows[i]) + [b_vec[i]] for i in range(m)]

    pivot_cols = []; row_idx = 0
    for col in range(n_unknowns + 1):  # include RHS column
        found = -1
        for rr in range(row_idx, m):
            if aug[rr][col] != 0: found = rr; break
        if found == -1: continue
        aug[row_idx], aug[found] = aug[found], aug[row_idx]
        piv = aug[row_idx][col]
        for j2 in range(n_unknowns + 1): aug[row_idx][j2] /= piv
        for rr in range(m):
            if rr == row_idx: continue
            if aug[rr][col] == 0: continue
            f = aug[rr][col]
            for j2 in range(n_unknowns + 1): aug[rr][j2] -= f * aug[row_idx][j2]
        pivot_cols.append(col); row_idx += 1

    rank_a = sum(1 for c in pivot_cols if c < n_unknowns)
    rank_aug = len(pivot_cols)

    if rank_aug > rank_a:
        return rank_a, rank_aug, None  # inconsistent

    # Extract solution (set free vars to 0)
    x = [Fraction(0)] * n_unknowns
    for pi, pc in enumerate(pivot_cols):
        if pc < n_unknowns:
            x[pc] = aug[pi][n_unknowns]
    return rank_a, rank_aug, x

print("\n=== Full band search (with underdetermined handling) ===", flush=True)

for width in range(5):
    for deg in range(6):
        labels = []
        for fname in ['R', 'K', 'P']:
            for r in range(width + 1):
                for d in range(deg + 1):
                    labels.append((fname, r, d))

        n_unk = len(labels)
        n_train = max(n_unk + 3, 15)
        n_hold = 5

        if n_train + n_hold + width >= NMAX: continue

        # Build system
        rows = []; bv = []
        for n in range(n_train):
            rq = []; rp = []
            for fname, r, d in labels:
                pair = fams[fname](n + r)
                rq.append(Fraction(n)**d * pair[0])
                rp.append(Fraction(n)**d * pair[1])
            rows.append(rq); bv.append(cq[n])
            rows.append(rp); bv.append(cp[n])

        rank_a, rank_aug, x = solve_augmented(rows, bv, n_unk)

        if x is None:
            if deg <= 1 or deg % 3 == 0:
                print(f"  w={width}, d={deg}: inconsistent (rank {rank_a}, aug {rank_aug})", flush=True)
            continue

        # Verify holdout
        ok = True
        for n in range(n_train, n_train + n_hold):
            pq = pp = Fraction(0)
            for idx, (fname, r, d) in enumerate(labels):
                pair = fams[fname](n + r)
                pq += x[idx] * Fraction(n)**d * pair[0]
                pp += x[idx] * Fraction(n)**d * pair[1]
            if pq != cq[n] or pp != cp[n]:
                ok = False; break

        if ok:
            print(f"\n*** MATCH: w={width}, d={deg}, unknowns={n_unk}, rank={rank_a} ***", flush=True)
            for idx, (fname, r, d) in enumerate(labels):
                if x[idx] != 0:
                    print(f"  {fname}_{r} n^{d}: {x[idx]}")

            # Full verify
            all_ok = True
            for n in range(n_train + n_hold, NMAX - width):
                pq = Fraction(0)
                for idx, (fname, r, d) in enumerate(labels):
                    pair = fams[fname](n+r)
                    pq += x[idx] * Fraction(n)**d * pair[0]
                if pq != cq[n]:
                    all_ok = False
                    print(f"  FAIL at n={n}")
                    break
            if all_ok:
                print(f"  Verified ALL n=0..{NMAX-width-1}")
            break
        else:
            if deg <= 1:
                print(f"  w={width}, d={deg}: holdout fails (rank={rank_a})", flush=True)

    else:
        continue
    break

print("\nDone.")
