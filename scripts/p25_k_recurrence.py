#!/usr/bin/env python3
"""
P2.5: Find a linear recurrence with polynomial coefficients in k for the
Delannoy-basis coefficients f(k) (from Q-hat) and g(k) (from P-hat).

Q_hat_N = sum_{k=0}^N f(k) * B(N,k),   B(N,k) = 2^k C(2k,k) C(N,k) C(N+k,k)
P_hat_N = sum_{k=0}^N g(k) * B(N,k)

Sequences considered (all exact Fractions):
  f_e1 : Delannoy coefficients of the e1-trajectory Q_hat^{e1}_N = (Prod M_H)_{1,1}
         (f(0)=1, f(1)=5749/3136, f(2)=16811771/4572288)
  f_q  : coefficients of the combined q-row Q_hat_N,  q=(33750,-36000,9000)
  g_p  : coefficients of the combined p-row P_hat_N,  p=(30921,-32972,8240)
         (g_p(k)/f_q(k) -> Catalan G geometrically)

Search: c_r(k) s(k+r) + ... + c_0(k) s(k) = 0 with deg c_j <= d.
Grid: the task grid r=2,3,4 / d=1..10 first (exact), then an extended grid
prescreened mod large primes with rational reconstruction, verified exactly.

RESOLUTION (2026-07-14): this grid is too small -- no recurrence exists with
r <= 8, d <= 24 (also ruled out rigorously: boxes (16,16), (20,12), (6,160)
have exact nullity 0 mod p). The true minimal recurrences, found by
ore_algebra guessing on 401 exact terms and verified exactly on all of them
(see p25_k_recurrence_guess2.sage / _guess3.sage):
  f_e1: order 8, degree 42     f_q: order 8, degree 41
  g_p : order 8, degree 42
  all three have Poincare polynomial (xi-1)^2 (xi+1/8)^6.
  A common operator LCLM(Lf_q, Lg_p) of order 9, degree 43 annihilates both
  f_q and g_p (verified exactly); GCRD has order 7.
  eps_k = g_p/f_q - G decays with ratio -> -1/8 (the Poincare root), NOT
  17-12*sqrt(2). f_e1, f_q, g_p > 0 for all k <= 400.
Also proved along the way: the inverse Delannoy triangle is proper
hypergeometric:
  B^{-1}(k,N) = (-1)^{k+N} (2N+1) C(k,N)
                / (2^k C(2k,k) (k+N+1) C(k+N,N)),
verified exactly as two-sided inverse up to k=20 (p25_k_recurrence_deep2.py),
so f(k) = sum_N B^{-1}(k,N) Qhat_N is holonomic by Zeilberger theory.
"""
from fractions import Fraction as F
from math import comb, gcd, lcm
from functools import reduce
import sys, time

KMAX = 140  # compute coefficients for k = 0..KMAX

# ----------------------------------------------------------------------
# CMF M_H(n) (copied from p25_partial_sum_check.py)
# ----------------------------------------------------------------------
def M_entries(n):
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
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

def delannoy_summand(N, k):
    if k < 0 or k > N:
        return F(0)
    return F(2**k * comb(2*k, k) * comb(N, k) * comb(N+k, k))

# ----------------------------------------------------------------------
# Step 1: CMF trajectories and triangular inversion
# ----------------------------------------------------------------------
t0 = time.time()
print(f"Step 1: computing trajectories for N = 0..{KMAX} (exact) ...")
sys.stdout.flush()

rows = {
    'e1': [F(1), F(0), F(0)],
    'e2': [F(0), F(1), F(0)],
    'e3': [F(0), F(0), F(1)],
}
history = {key: [list(v)] for key, v in rows.items()}

for N in range(KMAX):
    M = M_entries(N)
    d = F(delta_H(N))
    MH = [[F(M[i][j]) / d for j in range(3)] for i in range(3)]
    for key in rows:
        r = rows[key]
        new_r = [sum(r[i]*MH[i][k] for i in range(3)) for k in range(3)]
        rows[key] = new_r
        history[key].append(list(new_r))

q = [F(33750), F(-36000), F(9000)]
p = [F(30921), F(-32972), F(8240)]

Qe1_vals = [history['e1'][N][0] for N in range(KMAX+1)]
Q_vals = [sum(q[j] * history['e'+str(j+1)][N][0] for j in range(3))
          for N in range(KMAX+1)]
P_vals = [sum(p[j] * history['e'+str(j+1)][N][0] for j in range(3))
          for N in range(KMAX+1)]

def decompose_in_delannoy(vals, NMAX):
    coeffs = []
    for K in range(NMAX + 1):
        rhs = vals[K]
        for k in range(K):
            rhs -= coeffs[k] * delannoy_summand(K, k)
        coeffs.append(rhs / delannoy_summand(K, K))
    return coeffs

f_e1 = decompose_in_delannoy(Qe1_vals, KMAX)
f_q  = decompose_in_delannoy(Q_vals, KMAX)
g_p  = decompose_in_delannoy(P_vals, KMAX)

assert f_e1[0] == F(1), f_e1[0]
assert f_e1[1] == F(5749, 3136), f_e1[1]
assert f_e1[2] == F(16811771, 4572288), f_e1[2]
print("  sanity check on f(0), f(1), f(2) [e1-trajectory]: OK")
print(f"  coefficients computed for k = 0..{KMAX}  ({time.time()-t0:.1f}s)")
for nm, seq in (("f_e1", f_e1), ("f_q", f_q), ("g_p", g_p)):
    neg = [k for k in range(KMAX+1) if seq[k] <= 0]
    print(f"  {nm}(k) > 0 for all computed k:",
          "YES" if not neg else f"NO, fails at {neg}")
print()
sys.stdout.flush()

# ----------------------------------------------------------------------
# Linear algebra: exact and mod-p null spaces
# ----------------------------------------------------------------------
def nullspace_exact(A, ncols):
    A = [row[:] for row in A]
    m = len(A)
    pivots = []
    prow = 0
    for col in range(ncols):
        piv = next((r for r in range(prow, m) if A[r][col] != 0), None)
        if piv is None:
            continue
        A[prow], A[piv] = A[piv], A[prow]
        pv = A[prow][col]
        A[prow] = [x / pv for x in A[prow]]
        for r in range(m):
            if r != prow and A[r][col] != 0:
                fac = A[r][col]
                A[r] = [A[r][j] - fac * A[prow][j] for j in range(ncols)]
        pivots.append((prow, col))
        prow += 1
        if prow == m:
            break
    pivot_cols = {c for _, c in pivots}
    basis = []
    for fc in [c for c in range(ncols) if c not in pivot_cols]:
        v = [F(0)] * ncols
        v[fc] = F(1)
        for r, c in pivots:
            v[c] = -A[r][fc]
        basis.append(v)
    return basis

def nullspace_modp(A, ncols, p):
    """A: rows of ints mod p. Returns basis of null space mod p."""
    A = [row[:] for row in A]
    m = len(A)
    pivots = []
    prow = 0
    for col in range(ncols):
        piv = next((r for r in range(prow, m) if A[r][col] % p != 0), None)
        if piv is None:
            continue
        A[prow], A[piv] = A[piv], A[prow]
        inv = pow(A[prow][col], p - 2, p)
        A[prow] = [(x * inv) % p for x in A[prow]]
        for r in range(m):
            if r != prow and A[r][col] % p != 0:
                fac = A[r][col]
                A[r] = [(A[r][j] - fac * A[prow][j]) % p for j in range(ncols)]
        pivots.append((prow, col))
        prow += 1
        if prow == m:
            break
    pivot_cols = {c for _, c in pivots}
    basis = []
    for fc in [c for c in range(ncols) if c not in pivot_cols]:
        v = [0] * ncols
        v[fc] = 1
        for r, c in pivots:
            v[c] = (-A[r][fc]) % p
        basis.append(v)
    return basis

def normalize_int_vector(v):
    dens = [x.denominator for x in v if x != 0]
    L = reduce(lcm, dens, 1)
    ints = [int(x * L) for x in v]
    g = reduce(gcd, (abs(x) for x in ints if x != 0), 0)
    if g:
        ints = [x // g for x in ints]
    for x in ints:
        if x != 0:
            if x < 0:
                ints = [-y for y in ints]
            break
    return ints

def rational_reconstruct(a, m):
    """Find n/d = a mod m with |n|, d <= sqrt(m/2). Returns Fraction or None."""
    a %= m
    bound = int((m // 2) ** 0.5)
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1 > bound:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    if abs(s1) > bound or s1 == 0:
        return None
    if gcd(r1, abs(s1)) != 1:
        return None
    return F(r1, s1)

# ----------------------------------------------------------------------
# System construction
# ----------------------------------------------------------------------
def build_rows(seq, r, d, kmin, kmax_row, to_val):
    """Rows for k = kmin..kmax_row. Column order: j (shift) outer, k-power inner.
    to_val maps a Fraction to the field element."""
    A = []
    for k in range(kmin, kmax_row + 1):
        row = []
        for j in range(r + 1):
            fk = to_val(seq[k + j])
            kp = 1
            for _ in range(d + 1):
                row.append(fk * kp if not isinstance(fk, F) else fk * kp)
                kp *= k
        A.append(row)
    return A

def frac_modp(x, p):
    return (x.numerator % p) * pow(x.denominator % p, p - 2, p) % p

# ----------------------------------------------------------------------
# Search
# ----------------------------------------------------------------------
PRIMES = [(1 << 61) - 1, (1 << 61) - 31, (1 << 61) - 45]  # 2^61-1 etc. (primes)

def poly_eval(coeffs, k):
    acc = 0
    kp = 1
    for c in coeffs:
        acc += c * kp
        kp *= k
    return acc

def verify_recurrence(seq, r, polys, kmin=0):
    """Check sum_j c_j(k) seq(k+j) = 0 for all k in [kmin, len-1-r]."""
    for k in range(kmin, len(seq) - r):
        if sum(F(poly_eval(polys[j], k)) * seq[k + j] for j in range(r + 1)) != 0:
            return False, k
    return True, None

def earliest_valid_k(seq, r, polys):
    """Largest suffix on which the recurrence holds; returns smallest kmin."""
    ok_from = None
    for k in range(len(seq) - r - 1, -1, -1):
        if sum(F(poly_eval(polys[j], k)) * seq[k + j] for j in range(r + 1)) != 0:
            return k + 1
        ok_from = k
    return ok_from

def search_modp(seq, orders, degrees, kmin, name):
    """Prescreen (r,d) grid mod PRIMES[0]; on hit, reconstruct rationally from
    all PRIMES via CRT, return integer polynomial coefficients."""
    K = len(seq) - 1
    p0 = PRIMES[0]
    seq_mod = {p: [frac_modp(x, p) for x in seq] for p in PRIMES}
    # order the grid by number of unknowns (Occam first)
    grid = sorted(((r, d) for r in orders for d in degrees),
                  key=lambda t: ((t[0]+1)*(t[1]+1), t[0]))
    for r, d in grid:
        nunk = (r + 1) * (d + 1)
        kmax_row = K - r
        nrows = kmax_row - kmin + 1
        if nrows < nunk + 8:
            continue
        A0 = []
        sm = seq_mod[p0]
        for k in range(kmin, kmax_row + 1):
            row = []
            for j in range(r + 1):
                fk = sm[k + j]
                kp = 1
                for _ in range(d + 1):
                    row.append(fk * kp % p0)
                    kp *= k
            A0.append(row)
        basis0 = nullspace_modp(A0, nunk, p0)
        if not basis0:
            continue
        print(f"  [{name}] mod-p hit: r={r}, d={d}, nullity={len(basis0)} "
              f"(rows={nrows}, unknowns={nunk})")
        sys.stdout.flush()
        # reconstruct one null vector from the CRT of all primes
        # normalize each prime's vector so its first nonzero coord (same index) = 1
        idx = next(i for i, x in enumerate(basis0[0]) if x % p0 != 0)
        vecs = {}
        good = True
        for pp in PRIMES:
            sm2 = seq_mod[pp]
            A = []
            for k in range(kmin, kmax_row + 1):
                row = []
                for j in range(r + 1):
                    fk = sm2[k + j]
                    kp = 1
                    for _ in range(d + 1):
                        row.append(fk * kp % pp)
                        kp *= k
                A.append(row)
            b = nullspace_modp(A, nunk, pp)
            if len(b) != len(basis0):
                good = False
                break
            v = b[0]
            if v[idx] % pp == 0:
                good = False
                break
            inv = pow(v[idx], pp - 2, pp)
            vecs[pp] = [(x * inv) % pp for x in v]
        if not good:
            print(f"  [{name}]   (inconsistent across primes; skipping)")
            continue
        # CRT
        Mprod = 1
        for pp in PRIMES:
            Mprod *= pp
        crt = []
        for i in range(nunk):
            acc = 0
            for pp in PRIMES:
                Mi = Mprod // pp
                acc = (acc + vecs[pp][i] * Mi * pow(Mi % pp, pp - 2, pp)) % Mprod
            crt.append(acc)
        rec = [rational_reconstruct(a, Mprod) for a in crt]
        if any(x is None for x in rec):
            print(f"  [{name}]   (rational reconstruction failed; skipping)")
            continue
        ints = normalize_int_vector(rec)
        polys = [ints[j*(d+1):(j+1)*(d+1)] for j in range(r + 1)]
        return r, d, polys, len(basis0)
    return None

def poly_str(coeffs):
    terms = []
    for m, c in enumerate(coeffs):
        if c == 0:
            continue
        if m == 0:
            terms.append(f"{c}")
        elif m == 1:
            terms.append(f"{c}*k")
        else:
            terms.append(f"{c}*k^{m}")
    return (" + ".join(terms)).replace("+ -", "- ") if terms else "0"

def poincare_analysis(polys):
    def deg(pj):
        dd = -1
        for m, c in enumerate(pj):
            if c != 0:
                dd = m
        return dd
    D = max(deg(pj) for pj in polys)
    lead = [(pj[D] if len(pj) > D else 0) for pj in polys]
    g = reduce(gcd, (abs(x) for x in lead if x != 0), 0)
    if g:
        lead = [x // g for x in lead]
    if lead[-1] < 0:
        lead = [-x for x in lead]
    terms = []
    for j, c in enumerate(lead):
        if c == 0:
            continue
        terms.append(f"{c}" if j == 0 else (f"{c}*xi" if j == 1 else f"{c}*xi^{j}"))
    pstr = (" + ".join(terms)).replace("+ -", "- ")
    roots = None
    try:
        import numpy as np
        roots = sorted(np.roots(list(reversed([float(c) for c in lead]))),
                       key=lambda z: abs(z))
    except ImportError:
        pass
    return lead, pstr, roots

def report(seq, name):
    print(f"--- Searching recurrence for {name}(k), k = 0..{len(seq)-1} ---")
    sys.stdout.flush()
    # Stage A: the task grid, recurrence required from k=0 (exact prescreen mod p)
    res = search_modp(seq, orders=(2, 3, 4), degrees=range(1, 11), kmin=0,
                      name=name)
    if res is None:
        print(f"  [{name}] task grid (r=2..4, d=1..10, from k=0): no recurrence")
        # Stage B: extended grid, allow validity only from k >= 5
        res = search_modp(seq, orders=(2, 3, 4, 5, 6, 7, 8),
                          degrees=range(1, 25), kmin=5, name=name)
    if res is None:
        print(f"  [{name}] extended grid (r=2..8, d=1..24, from k=5): "
              f"no recurrence either")
        print()
        return None
    r, d, polys, nsdim = res
    print(f"  [{name}] FOUND recurrence: order r={r}, coefficient degree d={d}, "
          f"null-space dim={nsdim}")
    for j in range(r + 1):
        print(f"    c_{j}(k) = {poly_str(polys[j])}")
    k0 = earliest_valid_k(seq, r, polys)
    ok, badk = verify_recurrence(seq, r, polys, kmin=k0 if k0 else 0)
    print(f"  [{name}] EXACT verification: holds for all k = {k0}.."
          f"{len(seq)-1-r} " + ("(PASS)" if ok else f"FAIL at {badk}"))
    if k0 and k0 > 0:
        print(f"  [{name}]   note: recurrence starts holding at k = {k0}")
    lead, pstr, roots = poincare_analysis(polys)
    print(f"  [{name}] Poincare (leading-coeff) polynomial: {pstr}")
    if roots is not None:
        print(f"  [{name}] Poincare roots:")
        for z in roots:
            if abs(z.imag) < 1e-10:
                print(f"      {z.real:.12f}")
            else:
                print(f"      {z.real:.12f} + {z.imag:.12f}i")
    print()
    sys.stdout.flush()
    return r, d, polys

print("Step 2: recurrence search (mod-p prescreen + CRT rational "
      "reconstruction + exact verification)")
res_f = report(f_e1, "f_e1")
res_g = report(g_p, "g_p")
res_fq = report(f_q, "f_q")

# ----------------------------------------------------------------------
# Step 3: compare
# ----------------------------------------------------------------------
print("=" * 70)
if res_f and res_g:
    rf, df, pf = res_f
    rg, dg, pg = res_g
    if rf == rg and pf == pg:
        print("SAME RECURRENCE: f(k) and g(k) satisfy the IDENTICAL recurrence")
        print(f"  order r = {rf}, coefficient degree d = {df}")
    else:
        print("The minimal recurrences found differ in form:")
        print(f"  f_e1: r={rf}, d={df}; g_p: r={rg}, d={dg}")
        okg, badk = verify_recurrence(g_p, rf, pf,
                                      kmin=earliest_valid_k(g_p, rf, pf) or 0)
        k0g = earliest_valid_k(g_p, rf, pf)
        print(f"  Does g_p satisfy f_e1's recurrence (for k >= {k0g})? "
              f"{'YES' if okg and k0g is not None else 'NO'}")
        k0f = earliest_valid_k(f_e1, rg, pg)
        okf, badk = verify_recurrence(f_e1, rg, pg, kmin=k0f or 0)
        print(f"  Does f_e1 satisfy g_p's recurrence (for k >= {k0f})? "
              f"{'YES' if okf and k0f is not None else 'NO'}")
    if res_fq:
        rq, dq, pq = res_fq
        print(f"  f_q recurrence identical to f_e1's: "
              f"{'YES' if (rq, pq) == (rf, pf) else 'NO'}")
    elif res_f:
        k0q = earliest_valid_k(f_q, rf, pf)
        okq, _ = verify_recurrence(f_q, rf, pf, kmin=k0q or 0)
        print(f"  Does f_q satisfy f_e1's recurrence (k >= {k0q})? "
              f"{'YES' if okq and k0q is not None else 'NO'}")
else:
    print("At least one sequence had no recurrence in the searched range.")

print()
print("g_p(k)/f_q(k) tail (float):")
for k in range(KMAX - 4, KMAX + 1):
    if f_q[k] != 0:
        print(f"  k={k}: {float(g_p[k]/f_q[k]):.18f}")
G = 0.915965594177219015054603514932
print(f"  Catalan G:  {G:.18f}")
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done.")
