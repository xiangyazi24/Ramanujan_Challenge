#!/usr/bin/env python3
"""Creative telescoping attempt for Problem 2.5 (Catalan's constant G).

Goal: find a rational/closed-form kernel R(n,k) with
    Qhat_n = sum_{k=0}^n F_D(n,k) * R(n,k),
where F_D(n,k) = 2^k C(2k,k) C(n,k) C(n+k,k)   (so sum_k F_D(n,k) = D_n^2,
D_n = central Delannoy), and Qhat_n is the CMF denominator sequence
(pure powers of 2 in the denominator, Poincare roots {1, (3+-2sqrt2)^2}).

Strategy (all exact Fraction arithmetic):
  [1] Build Qhat_n for n = 0..NMAX from the CMF matrix M_H(n) = M(n)/delta_H(n),
      initial Q-row [33750, -36000, 9000], and sanity-check everything we
      "know" (2-power denominators, Delannoy identity, AZ constant term,
      Poincare ratio -> (3+2sqrt2)^2 = 33.97...).
  [2] Diagonal probe: R(n,n) = Qhat_n / F_D(n,n) if all other R(n,k) = 0.
  [3] Single-sequence probes: for each basis kernel h_j(k,n), form
      G_j(n) = sum_k F_D(n,k) h_j(k,n) and test whether Qhat_n / G_j(n)
      is a fixed rational function of n (exact linear rational fit).
  [4] Tiered ansatz: R(n,k) = sum_j a_j(n) h_j(k,n) with a_j(n) drawn from
      increasing spaces of weight functions of n (poly, 1/(n+1), 1/(2n+1),
      harmonic H_n, ...).  Each tier is an exact linear system over the
      training range n = 0..N_TRAIN; consistent solutions are re-verified
      on the held-out range n = N_TRAIN+1..NMAX.
  [5] Numeric least squares (if numpy present) to measure how much of Qhat
      each tier can capture even when the exact system is inconsistent.

Basis kernels h_j(k,n):
  1, H_k, H_{n-k}, H_{n+k}, 1/(2k+1), S_k = sum_{j<k} (-1)^j/(2j+1)^2,
  k/(n+1), C(2k,k)/4^k
(S_k is the natural G-partial-sum candidate; extended tiers add a few
alternating variants.)
"""

from fractions import Fraction
from math import comb, isqrt
import sys
import time

F = Fraction

NMAX = 50        # compute Qhat_n, F_D(n,k) for n = 0..NMAX
N_TRAIN = 44     # exact linear systems use n = 0..N_TRAIN; rest is hold-out

T0 = time.time()


def log(msg=""):
    print(msg, flush=True)


# ----------------------------------------------------------------------
# [1] CMF recurrence  ->  Qhat_n
# ----------------------------------------------------------------------

def M_entries(n):
    n = F(n)
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]


def delta_H(n):
    n = F(n)
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2


def compute_qhat(nmax):
    v = [F(33750), F(-36000), F(9000)]
    q = [v[0]]
    for n in range(nmax):
        M = M_entries(n)
        d = delta_H(n)
        v = [sum(v[k]*M[k][j] for k in range(3))/d for j in range(3)]
        q.append(v[0])
    return q


log("[1] Computing Qhat_n from the CMF (exact) ...")
Qhat = compute_qhat(NMAX)

# sanity: denominators are pure powers of 2
bad = []
for n, q in enumerate(Qhat):
    d = q.denominator
    if d & (d-1) != 0:  # not a power of two
        bad.append(n)
if bad:
    log(f"  !! FAIL: Qhat_n has non-2-power denominators at n={bad[:5]} ...")
    log("  (convention mismatch -- aborting)")
    sys.exit(1)
log(f"  OK: all Qhat_n (n<= {NMAX}) have pure 2-power denominators.")
log(f"  Qhat_0..4 = {[str(q) for q in Qhat[:5]]}")
r = Qhat[NMAX] / Qhat[NMAX-1]
log(f"  Qhat_{NMAX}/Qhat_{NMAX-1} = {float(r):.10f}   (target (3+2sqrt2)^2 = 33.9705627485)")

# Delannoy numbers and F_D
def FD(n, k):
    return (1 << k) * comb(2*k, k) * comb(n, k) * comb(n+k, k)


Dn = [sum(comb(n, k)**2 * (1 << k) for k in range(n+1)) for n in range(NMAX+1)]
assert Dn[:6] == [1, 3, 13, 63, 321, 1683], "Delannoy sanity failed"
for n in range(NMAX+1):
    assert Dn[n]**2 == sum(FD(n, k) for k in range(n+1)), f"D^2 identity fails n={n}"
log("  OK: D_n central Delannoy, and D_n^2 = sum_k F_D(n,k) verified.")

# AZ phase check: D_n = CT_u[ phi(u)^n ], phi(u) = (1+u)(u+2)/u
for n in range(11):
    # ((1+u)(u+2))^n coefficient of u^n
    poly = [1]
    base = [2, 3, 1]  # (1+u)(u+2) = 2 + 3u + u^2
    for _ in range(n):
        new = [0]*(len(poly)+2)
        for i, c in enumerate(poly):
            for j, b in enumerate(base):
                new[i+j] += c*b
        poly = new
    assert poly[n] == Dn[n], f"AZ phase CT check fails n={n}"
log("  OK: D_n = CT_u[phi(u)^n], phi(u) = (1+u)(u+2)/u  (checked n<=10).")


# ----------------------------------------------------------------------
# [2] Diagonal probe: R(n,n) = Qhat_n / F_D(n,n)
# ----------------------------------------------------------------------

log("\n[2] Diagonal probe R(n,n) = Qhat_n / F_D(n,n):")
diag = []
for n in range(min(16, NMAX)+1):
    rv = Qhat[n] / FD(n, n)
    diag.append(rv)
    log(f"  n={n:2d}: R = {rv}  ~ {float(rv):+.6e}")
log("  ratios R(n+1,n+1)/R(n,n):")
for n in range(1, min(12, len(diag)-1)+1):
    if diag[n-1] != 0:
        log(f"    n={n:2d}: {float(diag[n]/diag[n-1]):+.6f}")


# ----------------------------------------------------------------------
# basis kernels h_j(k,n) and the weighted sums G_j(n)
# ----------------------------------------------------------------------

H = [F(0)]
for i in range(1, 2*NMAX+2):
    H.append(H[-1] + F(1, i))

S = [F(0)]  # S_k = sum_{j=0}^{k-1} (-1)^j/(2j+1)^2  -> G
for j in range(0, NMAX+1):
    S.append(S[-1] + F((-1)**j, (2*j+1)**2))

CB = [F(comb(2*k, k), 4**k) for k in range(NMAX+1)]  # C(2k,k)/4^k

T2 = [F(0)]  # T_k = sum_{j<k} 1/(2j+1)^2   (non-alternating)
for j in range(0, NMAX+1):
    T2.append(T2[-1] + F(1, (2*j+1)**2))
LZ = [F(0)]  # Leibniz partial sum: sum_{j<k} (-1)^j/(2j+1)  -> pi/4
for j in range(0, NMAX+1):
    LZ.append(LZ[-1] + F((-1)**j, 2*j+1))

BASIS = [
    ("1",        lambda n, k: F(1)),
    ("H_k",      lambda n, k: H[k]),
    ("H_{n-k}",  lambda n, k: H[n-k]),
    ("H_{n+k}",  lambda n, k: H[n+k]),
    ("1/(2k+1)", lambda n, k: F(1, 2*k+1)),
    ("S_k",      lambda n, k: S[k]),
    ("k/(n+1)",  lambda n, k: F(k, n+1)),
    ("C(2k,k)/4^k", lambda n, k: CB[k]),
]

EXTRA_BASIS = [
    ("(-1)^k",       lambda n, k: F((-1)**k)),
    ("(-1)^k S_k",   lambda n, k: F((-1)**k) * S[k]),
    ("(-1)^k/(2k+1)", lambda n, k: F((-1)**k, 2*k+1)),
    ("k",            lambda n, k: F(k)),
    ("k^2",          lambda n, k: F(k*k)),
    ("1/(2k-1)",     lambda n, k: F(1, 2*k-1)),
    ("H_{n+k}-H_{n-k}", lambda n, k: H[n+k]-H[n-k]),
    ("4^k/C(2k,k)",  lambda n, k: 1/CB[k]),
    ("T_k",          lambda n, k: T2[k]),          # sum_{j<k} 1/(2j+1)^2
    ("(-1)^k T_k",   lambda n, k: F((-1)**k)*T2[k]),
    ("Leib_k",       lambda n, k: LZ[k]),          # sum_{j<k} (-1)^j/(2j+1)
]

log("\n[3] Building weighted sums G_j(n) = sum_k F_D(n,k) h_j(k,n) ...")
FDtab = [[FD(n, k) for k in range(n+1)] for n in range(NMAX+1)]


def weighted_sum(basis_fn, n):
    return sum(FDtab[n][k] * basis_fn(n, k) for k in range(n+1))


ALL_BASIS = BASIS + EXTRA_BASIS
Gtab = {}
for name, fn in ALL_BASIS:
    Gtab[name] = [weighted_sum(fn, n) for n in range(NMAX+1)]
log(f"  done ({time.time()-T0:.1f}s).  G_0(n) = D_n^2 check: "
    f"{'OK' if Gtab['1'] == [F(d*d) for d in Dn] else 'FAIL'}")


# ----------------------------------------------------------------------
# exact linear algebra helpers
# ----------------------------------------------------------------------

def rref_solve(A, b):
    """Exact RREF on [A|b]. Returns (consistent, rank, particular_solution)."""
    m = len(A)
    ncols = len(A[0]) if m else 0
    rows = [list(A[i]) + [b[i]] for i in range(m)]
    piv_cols = []
    r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, m) if rows[i][c] != 0), None)
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(m):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [x - f*y for x, y in zip(rows[i], rows[r])]
        piv_cols.append(c)
        r += 1
        if r == m:
            break
    consistent = all(rows[i][ncols] == 0 for i in range(r, m))
    sol = None
    if consistent:
        sol = [F(0)]*ncols
        for i, c in enumerate(piv_cols):
            sol[c] = rows[i][ncols]
    return consistent, r, sol


def rational_fit(pairs, dp, dq, verify_pairs):
    """Try s(n) = p(n)/q(n), deg p<=dp, deg q<=dq, exact.
    pairs: [(n, s_n)] training; verify on verify_pairs. Returns (p,q) or None."""
    ncols = (dp+1) + (dq+1)
    A = []
    for n, s in pairs:
        row = [F(n)**i for i in range(dp+1)] + [-s * F(n)**i for i in range(dq+1)]
        A.append(row)
    # homogeneous nullspace via RREF
    m = len(A)
    rows = [list(r) for r in A]
    piv_cols = []
    r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, m) if rows[i][c] != 0), None)
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        pv = rows[r][c]
        rows[r] = [x/pv for x in rows[r]]
        for i in range(m):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [x - f*y for x, y in zip(rows[i], rows[r])]
        piv_cols.append(c)
        r += 1
        if r == m:
            break
    free = [c for c in range(ncols) if c not in piv_cols]
    if not free:
        return None
    v = [F(0)]*ncols
    v[free[0]] = F(1)
    for i in range(len(piv_cols)-1, -1, -1):
        pc = piv_cols[i]
        v[pc] = -sum(rows[i][c]*v[c] for c in free)
    p = v[:dp+1]
    q = v[dp+1:]
    if all(x == 0 for x in q):
        return None
    # verify
    for n, s in verify_pairs:
        qn = sum(q[i]*F(n)**i for i in range(dq+1))
        pn = sum(p[i]*F(n)**i for i in range(dp+1))
        if qn == 0 or pn/qn != s:
            return None
    return p, q


# ----------------------------------------------------------------------
# [3] single-sequence probes: Qhat_n / G_j(n) rational in n?
# ----------------------------------------------------------------------

log("\n[3] Single-kernel probes: is Qhat_n = a(n) * G_j(n) with a(n) rational?")
DEG = 8
for name, _ in ALL_BASIS:
    Gj = Gtab[name]
    pairs = [(n, Qhat[n]/Gj[n]) for n in range(NMAX+1) if Gj[n] != 0]
    if len(pairs) < 2*DEG + 6:
        log(f"  {name:16s}: skipped (too many zeros)")
        continue
    train = pairs[:2*DEG+4]
    verify = pairs[2*DEG+4:]
    res = rational_fit(train, DEG, DEG, verify)
    if res:
        p, q = res
        log(f"  {name:16s}: *** RATIONAL FIT FOUND ***  p={p}  q={q}")
    else:
        v0 = float(pairs[-1][1]) if abs(pairs[-1][1]) < F(10)**300 else None
        log(f"  {name:16s}: no rational a(n) up to deg {DEG}."
            f"  (Qhat/G at n={pairs[-1][0]}: {v0:+.4e})" if v0 is not None else
            f"  {name:16s}: no rational a(n) up to deg {DEG}.")

# special look: u_n = Qhat_n / D_n^2 (same Poincare growth 33.97^n)
log("\n  u_n = Qhat_n / D_n^2 (both grow like (3+2sqrt2)^{2n}):")
for n in [0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50]:
    if n <= NMAX:
        u = Qhat[n] / F(Dn[n]**2)
        log(f"    n={n:2d}: u_n = {float(u):+.10f}")


# ----------------------------------------------------------------------
# [4] tiered creative-telescoping ansatz (exact)
# ----------------------------------------------------------------------

WEIGHTS = {
    "1":        lambda n: F(1),
    "n":        lambda n: F(n),
    "n^2":      lambda n: F(n*n),
    "n^3":      lambda n: F(n**3),
    "1/(n+1)":  lambda n: F(1, n+1),
    "1/(2n+1)": lambda n: F(1, 2*n+1),
    "1/(2n+3)": lambda n: F(1, 2*n+3),
    "H_n":      lambda n: H[n],
    "H_{2n+1}": lambda n: H[2*n+1],
}

# kernels suggested by the window-fit structure in section [6]:
# the free fit locks a[H_k] ~ -a[H_{n+k}] and zeroes S_k / H_{n-k}
REDUCED_BASIS = [
    ("1",            lambda n, k: F(1)),
    ("H_{n+k}-H_k",  lambda n, k: H[n+k] - H[k]),
    ("1/(2k+1)",     lambda n, k: F(1, 2*k+1)),
    ("k/(n+1)",      lambda n, k: F(k, n+1)),
    ("C(2k,k)/4^k",  lambda n, k: CB[k]),
]
for _name, _fn in REDUCED_BASIS:
    if _name not in Gtab:
        Gtab[_name] = [weighted_sum(_fn, n) for n in range(NMAX+1)]

TIERS = [
    ("A: a_j const",              BASIS,             ["1"]),
    ("B: a_j linear",             BASIS,             ["1", "n"]),
    ("C: a_j quadratic",          BASIS,             ["1", "n", "n^2"]),
    ("D: a_j cubic",              BASIS,             ["1", "n", "n^2", "n^3"]),
    ("E: quad + 1/(n+1),1/(2n+1)", BASIS,            ["1", "n", "n^2", "1/(n+1)", "1/(2n+1)"]),
    ("F: {1,n} x harmonic-n",     BASIS,             ["1", "n", "H_n", "H_{2n+1}"]),
    ("G: extended kernels, {1,n}", BASIS+EXTRA_BASIS, ["1", "n"]),
    ("H: reduced kernels, rich n-weights", REDUCED_BASIS,
     ["1", "n", "n^2", "1/(n+1)", "1/(2n+1)", "1/(2n+3)", "H_n", "H_{2n+1}"]),
]

log("\n[4] Tiered exact ansatz  Qhat_n = sum_j a_j(n) G_j(n),  "
    f"train n=0..{N_TRAIN}, verify n={N_TRAIN+1}..{NMAX}:")

found_any = False
for tier_name, basis_list, weight_names in TIERS:
    cols = []       # (label, sequence over n)
    for bname, _ in basis_list:
        for wname in weight_names:
            wfn = WEIGHTS[wname]
            seq = [Gtab[bname][n] * wfn(n) for n in range(NMAX+1)]
            cols.append((f"{wname} * {bname}", seq))
    nunk = len(cols)
    if nunk > N_TRAIN - 2:  # need >= 3 surplus training equations
        log(f"  Tier {tier_name}: skipped ({nunk} unknowns vs {N_TRAIN+1} equations)")
        continue
    A = [[col[1][n] for col in cols] for n in range(N_TRAIN+1)]
    b = [Qhat[n] for n in range(N_TRAIN+1)]
    t1 = time.time()
    consistent, rank, sol = rref_solve(A, b)
    status = ""
    if consistent:
        ok = all(
            sum(sol[j]*cols[j][1][n] for j in range(nunk)) == Qhat[n]
            for n in range(N_TRAIN+1, NMAX+1)
        )
        if ok:
            found_any = True
            log(f"  Tier {tier_name}: *** CONSISTENT + VERIFIED on hold-out *** "
                f"({nunk} unknowns, rank {rank}, {time.time()-t1:.1f}s)")
            for j in range(nunk):
                if sol[j] != 0:
                    log(f"      a[{cols[j][0]}] = {sol[j]}")
        else:
            log(f"  Tier {tier_name}: consistent on training but FAILS hold-out "
                f"(overfit; {nunk} unknowns, rank {rank}, {time.time()-t1:.1f}s)")
    else:
        log(f"  Tier {tier_name}: INCONSISTENT "
            f"({nunk} unknowns, rank {rank}, {time.time()-t1:.1f}s)")


# ----------------------------------------------------------------------
# [5] numeric least squares: how much can each tier capture?
# ----------------------------------------------------------------------

log("\n[5] Numeric least-squares residuals (rows normalized by Qhat_n):")
try:
    import numpy as np
    for tier_name, basis_list, weight_names in TIERS:
        cols = []
        for bname, _ in basis_list:
            for wname in weight_names:
                wfn = WEIGHTS[wname]
                cols.append([Gtab[bname][n] * wfn(n) for n in range(NMAX+1)])
        nunk = len(cols)
        rows_n = [n for n in range(NMAX+1) if Qhat[n] != 0]
        Anp = np.array([[float(cols[j][n]/Qhat[n]) for j in range(nunk)]
                        for n in rows_n])
        bnp = np.ones(len(rows_n))
        x, *_ = np.linalg.lstsq(Anp, bnp, rcond=None)
        resid = Anp @ x - bnp
        log(f"  Tier {tier_name}: max |rel resid| = {np.max(np.abs(resid)):.3e}, "
            f"rms = {np.sqrt(np.mean(resid**2)):.3e}")
except ImportError:
    log("  (numpy not available; skipped)")

# ----------------------------------------------------------------------
# [6] window-fit probe: is Qhat = (Tier-B combination) + subdominant mode?
# ----------------------------------------------------------------------
# The recurrence has Poincare roots {1, (3-2sqrt2)^2, (3+2sqrt2)^2}. If
# Qhat_n = (clean combination of G_j with linear-in-n coefficients)
#          + (solution component with subdominant growth ~1^n or 0.029^n),
# then an exact square fit on a window of large n gives coefficients that
# converge (window to window) and the residual at small n exposes the
# growth rate of the unmodeled component.

log("\n[6] Window-fit probe (Tier B: {1,n} x 8 kernels, 16 unknowns):")
tierB_cols = []
for bname, _ in BASIS:
    for wname in ["1", "n"]:
        wfn = WEIGHTS[wname]
        tierB_cols.append((f"{wname} * {bname}",
                           [Gtab[bname][n] * wfn(n) for n in range(NMAX+1)]))
nunkB = len(tierB_cols)

def window_fit(n_lo):
    ns = list(range(n_lo, n_lo + nunkB))
    A = [[c[1][n] for c in tierB_cols] for n in ns]
    b = [Qhat[n] for n in ns]
    consistent, rank, sol = rref_solve(A, b)
    return sol if consistent else None

for W in (20, 35):
    if W + nunkB - 1 > NMAX:
        continue
    sol = window_fit(W)
    if sol is None:
        log(f"  window n={W}..{W+nunkB-1}: singular/inconsistent square system")
        continue
    log(f"  window n={W}..{W+nunkB-1}: exact square fit found. Coefficients "
        f"(float, and limit_denominator guess):")
    for j in range(nunkB):
        fv = float(sol[j])
        guess = F(fv).limit_denominator(10**6)
        log(f"      a[{tierB_cols[j][0]:16s}] = {fv:+.12e}   ~ {guess}")
    # exact residual outside the window
    log(f"    residual e_n = Qhat_n - model_n (exact, from window n>={W}):")
    prev = None
    for n in range(0, min(W, 20)):
        e = Qhat[n] - sum(sol[j]*tierB_cols[j][1][n] for j in range(nunkB))
        fe = float(e)
        rel = float(e / Qhat[n]) if Qhat[n] != 0 else float('nan')
        ratio = (fe/prev) if (prev not in (None, 0.0)) else None
        log(f"      n={n:2d}: e = {fe:+.6e}  rel = {rel:+.3e}"
            + (f"  e_n/e_(n-1) = {ratio:+.4f}" if ratio is not None else ""))
        prev = fe

log(f"\nTotal time: {time.time()-T0:.1f}s")
if not found_any:
    log("RESULT: no exact creative-telescoping representation found in the "
        "tested ansatz spaces (see tier statuses above).")
