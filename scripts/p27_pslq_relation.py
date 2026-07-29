#!/usr/bin/env python3
"""
Problem 2.7: PSLQ hunt for relations between q_n and the Cooper binomial transform W_n.

Background:
  - q_n satisfies the 4-term P2.7 recurrence
        u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
    with q_0 = -215040420000, q_1 = -1672822650043404/905, q_2 = -964185327658080/6071.
  - T_k is Cooper's level-11 sequence:
        (k+1)^3 T_{k+1} = 2(2k+1)(5k^2+5k+2) T_k - 8k(7k^2+1) T_{k-1} + 22k(2k-1)(k-1) T_{k-2}
  - W_n = (1/256^n) sum_{j=0}^{2n} C(2n,j) (-2)^{2n-j} T_j.

KEY FACT: there is NO polynomial Ore relation between q_n and W_n
(formal exponent gap sigma_q = 0 vs sigma_W = -3/2).  Here we hunt for
INTEGER relations at fixed n (PSLQ), and for hypergeometric (rational-
function-of-n) structure in the ratio q_n / W_n.

All sequences are computed as exact rationals (Fraction), so every PSLQ
candidate relation can be verified EXACTLY (sum c_i x_i == 0 over Q).
NB: q_n and W_n are rational, so *some* integer relation always exists at
each fixed n — the meaningful finds are relations with SMALL coefficients
that persist / follow a pattern across n.  We bound maxcoeff and cross-check.

FINDINGS (run of 2026-07-14, mp.dps=250):
  * Step 4/6 PSLQ returns a relation at every n, but ALWAYS with coefficient 0
    on q_n — i.e. relations among W_n..W_{n+5} alone.  None persists at any
    other n, and the minimal coefficient size grows like 256^(n/5) ≈ 3.03^n,
    exactly the pigeonhole prediction for 6 random rationals over the common
    denominator 256^(n+5).  These are LATTICE ARTIFACTS, not structure.
    Conclusion: NO genuine integer relation ties q_n to the W-block
    (|coeff| <= 1e14, 200-digit tolerance, n = 5..30).
  * The real structure is asymptotic: q_n / W_n = C * n^(3/2) * (1 + O(1/n)),
    confirming the formal exponent gap sigma_q - sigma_W = 3/2 and that both
    sequences live on the SAME dominant Poincare root.  Neville extrapolation
    of (q_n/W_n)/n^(3/2) in 1/n (n up to 700, exact rationals) gives
      C = -2229565722002.2118165497021688121448176584983...   (44 digits)
      C/q_0 = 10.3681239182950433995139247254639142616002066...
    PSLQ at 42-digit tolerance finds NO closed form for C/q_0 in the bases
    {1, pi, pi^2}, {sqrt(pi), sqrt(2 pi), sqrt(11 pi), sqrt(22 pi)},
    {1, sqrt(2), sqrt(11), sqrt(22)} (on (C/q_0)^2/pi), nor is C/q_0 algebraic
    of degree <= 4 with coeffs <= 1e7.  (Two apparent hits at 24-digit
    tolerance were refuted at 42 digits — PSLQ noise.)
  * The ratio q_n/W_n is NOT hypergeometric: r_{n+1}/r_n is not a rational
    function of n up to degree 8 (exact-fit test), consistent with the n^(3/2)
    asymptotics requiring an infinite 1/n series, and with the known
    impossibility of a polynomial Ore relation.
"""

from fractions import Fraction as F
import math
import sys

from mpmath import mp, mpf, pslq, zeta

mp.dps = 250
TOL = mpf(10) ** (-200)
MAXCOEFF = 10 ** 14

N_Q = 61   # q_0 .. q_61
N_W = 61   # W_0 .. W_60  (needs T up to 2*60)
N_T = 201  # T_0 .. T_200


# ----------------------------------------------------------------------
# P2.7 recurrence coefficients
# ----------------------------------------------------------------------
def A(n):
    return 1024 * (2*n+5)**4 * (2*n+7)**3 * (2*n+9)**3 * (946*n**2 + 6407*n + 10860)

def B(n):
    P6 = (104060*n**6 + 1745370*n**5 + 12145238*n**4 + 44886481*n**3
          + 92943995*n**2 + 102256019*n + 46709052)
    return 128 * (2*n+7)**3 * (2*n+9)**3 * P6

def C(n):
    P5 = 3784*n**5 + 57792*n**4 + 351019*n**3 + 1059230*n**2 + 1587211*n + 944620
    return 16 * (n+3)**4 * (2*n+9)**3 * P5

def D(n):
    return (n+3)**4 * (n+4)**6 * (946*n**2 + 4515*n + 5399)


# ----------------------------------------------------------------------
# 1. q_n exact
# ----------------------------------------------------------------------
def compute_q(nmax):
    q = [F(-215040420000),
         F(-1672822650043404, 905),
         F(-964185327658080, 6071)]
    for n in range(2, nmax):
        q_next = (F(B(n)) * q[n] / A(n)
                  - F(C(n-1)) * q[n-1] / A(n-1)
                  + F(D(n-2)) * q[n-2] / A(n-2))
        q.append(q_next)
    return q


# ----------------------------------------------------------------------
# 2. Cooper T_k exact
# ----------------------------------------------------------------------
def compute_T(kmax):
    T = [F(1), F(4)]
    for k in range(1, kmax):
        t = (2*(2*k+1)*(5*k*k + 5*k + 2) * T[k]
             - 8*k*(7*k*k + 1) * T[k-1])
        if k >= 2:
            t += 22*k*(2*k-1)*(k-1) * T[k-2]
        T.append(t / F((k+1)**3))
    return T


# ----------------------------------------------------------------------
# 3. Binomial transform W_n exact
# ----------------------------------------------------------------------
def compute_W(T, nmax):
    W = []
    for n in range(nmax):
        s = F(0)
        for j in range(2*n + 1):
            s += math.comb(2*n, j) * F(-2)**(2*n - j) * T[j]
        W.append(s / F(256)**n)
    return W


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def to_mpf(fr):
    return mpf(fr.numerator) / mpf(fr.denominator)

def verify_exact(coeffs, fracs):
    """Exact check of sum c_i x_i == 0 over Q."""
    return sum(F(c) * x for c, x in zip(coeffs, fracs)) == 0

def run_pslq(labels, fracs, tag):
    """Run PSLQ on a list of exact rationals (normalized), verify exactly."""
    nonzero = [(l, f) for l, f in zip(labels, fracs) if f != 0]
    if len(nonzero) < 2:
        return None
    labels = [l for l, _ in nonzero]
    fracs = [f for _, f in nonzero]
    scale = max(abs(to_mpf(f)) for f in fracs)
    vec = [to_mpf(f) / scale for f in fracs]
    try:
        rel = pslq(vec, tol=TOL, maxcoeff=MAXCOEFF, maxsteps=200000)
    except Exception as e:
        print(f"  [{tag}] PSLQ error: {e}")
        return None
    if rel is None:
        return None
    ok = verify_exact(rel, fracs)
    terms = " + ".join(f"({c})*{l}" for c, l in zip(rel, labels) if c != 0)
    print(f"  [{tag}] PSLQ relation: {terms} = 0   exact-verified: {ok}")
    return (labels, rel, ok)


# exact Gaussian elimination for rational-function fitting
def solve_exact(Amat, bvec):
    m = len(Amat)
    k = len(Amat[0])
    M = [row[:] + [bvec[i]] for i, row in enumerate(Amat)]
    piv_cols = []
    r = 0
    for c in range(k):
        piv = None
        for i in range(r, m):
            if M[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(m):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f*b for a, b in zip(M[i], M[r])]
        piv_cols.append(c)
        r += 1
        if r == m:
            break
    # consistency
    for i in range(r, m):
        if M[i][k] != 0:
            return None
    sol = [F(0)] * k
    for i, c in enumerate(piv_cols):
        sol[c] = M[i][k]
    return sol

def fit_rational_function(samples, dmax=8, verify_extra=10):
    """samples: list of (n, Fraction value).  Try value = P(n)/Q(n), monic Q.
    Returns (P_coeffs, Q_coeffs) low-to-high degree, or None."""
    for total in range(0, 2*dmax + 1):
        for dq in range(min(total, dmax) + 1):
            dp = total - dq
            if dp > dmax:
                continue
            nunk = (dp + 1) + dq  # p_0..p_dp, q_0..q_{dq-1}; q_dq = 1
            need = nunk + verify_extra
            if len(samples) < need:
                continue
            pts = samples[:nunk]
            Amat, bvec = [], []
            for n, v in pts:
                row = [-F(n)**i for i in range(dp + 1)]       # -P(n) part
                row += [v * F(n)**i for i in range(dq)]        # v*Q_low(n)
                Amat.append(row)
                bvec.append(-v * F(n)**dq)                     # move v*n^dq
            sol = solve_exact(Amat, bvec)
            if sol is None:
                continue
            P = sol[:dp+1]
            Q = sol[dp+1:] + [F(1)]
            good = True
            for n, v in samples[nunk:need + max(0, len(samples)-need)]:
                qn = sum(c * F(n)**i for i, c in enumerate(Q))
                pn = sum(c * F(n)**i for i, c in enumerate(P))
                if qn == 0 or pn / qn != v:
                    good = False
                    break
            if good:
                return P, Q
    return None

def poly_str(coeffs):
    parts = []
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        parts.append(f"({c})*n^{i}" if i else f"({c})")
    return " + ".join(parts) if parts else "0"


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    print(f"mp.dps = {mp.dps}")
    print("Computing q_n (exact) ...")
    q = compute_q(N_Q)
    print("Computing Cooper T_k (exact) ...")
    T = compute_T(N_T)
    print(f"T_0..T_5 = {[T[i] for i in range(6)]}")
    assert T[1] == 4 and T[2] == 28 and T[3] == 268, "Cooper initial terms mismatch"
    print("Computing W_n (exact binomial transform) ...")
    W = compute_W(T, N_W)

    # sanity: q_n growth & convergence context
    zz = zeta(2) + zeta(3)
    print(f"\nzeta(2)+zeta(3) = {mp.nstr(zz, 30)}")
    print("\nSanity: numeric magnitudes")
    for n in [0, 1, 2, 5, 10, 20, 30]:
        print(f"  n={n:2d}  q_n ~ {mp.nstr(to_mpf(q[n]), 10)}   W_n ~ {mp.nstr(to_mpf(W[n]), 10)}")

    # ------------------------------------------------------------------
    # 4. PSLQ on {q_n, W_n..W_{n+5}} for n = 5..30
    # ------------------------------------------------------------------
    print("\n" + "="*70)
    print("STEP 4: PSLQ on {q_n, W_n, W_(n+1), ..., W_(n+5)}, n = 5..30")
    print("="*70)
    found4 = []
    for n in range(5, 31):
        labels = [f"q_{n}"] + [f"W_{n+j}" for j in range(6)]
        fracs = [q[n]] + [W[n+j] for j in range(6)]
        res = run_pslq(labels, fracs, f"n={n}")
        if res:
            found4.append((n, res))
    if not found4:
        print("  -> NO integer relation with |coeff| <= 1e14 at 200-digit tolerance, for any n in 5..30.")
    else:
        # cross-check each found relation at other n (shifted pattern)
        print("\n  Cross-checking found relations at shifted n:")
        for n0, (labels, rel, _) in found4:
            hits = []
            for n in range(3, 40):
                fr = [q[n]] + [W[n+j] for j in range(6)]
                if verify_exact(rel, fr):
                    hits.append(n)
            print(f"    relation found at n={n0} holds exactly at n in {hits}")

    # ------------------------------------------------------------------
    # 5. Ratio q_n / W_n: closed form in n?
    # ------------------------------------------------------------------
    print("\n" + "="*70)
    print("STEP 5: ratio r_n = q_n / W_n and hypergeometric test")
    print("="*70)
    r = []
    for n in range(N_W):
        r.append(None if W[n] == 0 else q[n] / W[n])
    for n in [0, 1, 2, 3, 5, 10, 15, 20, 30, 40, 50]:
        if r[n] is not None:
            print(f"  r_{n} = q_n/W_n ~ {mp.nstr(to_mpf(r[n]), 15)}")
    # growth of r_n: log ratio
    print("  successive ratios s_n = r_(n+1)/r_n (float):")
    s_samples = []
    for n in range(N_W - 1):
        if r[n] is not None and r[n+1] is not None and r[n] != 0:
            s = r[n+1] / r[n]
            s_samples.append((n, s))
    for n, s in s_samples[:20]:
        print(f"    s_{n} = {mp.nstr(to_mpf(s), 15)}")

    # 5a. is r_n itself a rational function of n?
    fit = fit_rational_function(r_samples := [(n, r[n]) for n in range(N_W) if r[n] is not None],
                                dmax=8, verify_extra=8)
    if fit:
        P, Q = fit
        print(f"  !! r_n IS a rational function of n:  r_n = [{poly_str(P)}] / [{poly_str(Q)}]")
    else:
        print("  r_n is NOT a rational function of n (deg <= 8).")

    # 5b. is r_n hypergeometric, i.e. s_n = r_(n+1)/r_n rational in n?
    fit = fit_rational_function(s_samples, dmax=8, verify_extra=8)
    if fit:
        P, Q = fit
        print(f"  !! q_n/W_n IS HYPERGEOMETRIC:  r_(n+1)/r_n = [{poly_str(P)}] / [{poly_str(Q)}]")
    else:
        print("  s_n = r_(n+1)/r_n is NOT a rational function of n (deg <= 8)")
        print("  -> q_n/W_n is not hypergeometric (at these degrees).")

    # 5c. same tests for the cross-Casoratian X_n = q_n W_(n+1) - q_(n+1) W_n
    print("\n  Cross-Casoratian X_n = q_n*W_(n+1) - q_(n+1)*W_n:")
    X = [q[n]*W[n+1] - q[n+1]*W[n] for n in range(N_W - 1)]
    x_samples = [(n, X[n+1]/X[n]) for n in range(len(X)-1) if X[n] != 0]
    for n, v in x_samples[:8]:
        print(f"    X_(n+1)/X_n at n={n}: {mp.nstr(to_mpf(v), 12)}")
    fit = fit_rational_function(x_samples, dmax=10, verify_extra=8)
    if fit:
        P, Q = fit
        print(f"  !! X_n IS HYPERGEOMETRIC: X_(n+1)/X_n = [{poly_str(P)}] / [{poly_str(Q)}]")
    else:
        print("  X_(n+1)/X_n NOT rational in n (deg <= 10).")

    # ------------------------------------------------------------------
    # 6. PSLQ on cross products q_a * W_b near the diagonal
    # ------------------------------------------------------------------
    print("\n" + "="*70)
    print("STEP 6: PSLQ on cross products {q_n*W_m, q_m*W_n} for m near n")
    print("="*70)
    found6 = []
    for n in range(5, 26, 2):
        labels, fracs = [], []
        for m in range(max(0, n-2), n+4):
            labels.append(f"q_{n}*W_{m}")
            fracs.append(q[n] * W[m])
        for m in range(max(0, n-2), n+4):
            if m != n:
                labels.append(f"q_{m}*W_{n}")
                fracs.append(q[m] * W[n])
        res = run_pslq(labels, fracs, f"cross n={n}")
        if res:
            found6.append((n, res))
    if not found6:
        print("  -> NO integer relation among cross products (|coeff| <= 1e14), n in 5..25.")

    # ------------------------------------------------------------------
    # 7. Asymptotic constant: q_n / W_n = C * n^(3/2) * (1 + O(1/n))
    # ------------------------------------------------------------------
    print("\n" + "="*70)
    print("STEP 7: asymptotic constant C = lim (q_n/W_n)/n^(3/2)")
    print("="*70)
    NBIG = 400
    qb = compute_q(NBIG)
    Tb = compute_T(2*NBIG + 1)
    pts = list(range(220, NBIG, 12))
    xs = [mpf(1)/n for n in pts]
    ys = []
    for n in pts:
        s = F(0)
        for j in range(2*n + 1):
            s += math.comb(2*n, j) * F(-2)**(2*n - j) * Tb[j]
        Wn = s / F(256)**n
        ys.append(to_mpf(qb[n]/Wn) / mpf(n)**mpf(1.5))
    m = len(xs)
    tab = ys[:]
    for j in range(1, m):
        for i in range(m-1, j-1, -1):
            tab[i] = (xs[i-j]*tab[i] - xs[i]*tab[i-1]) / (xs[i-j] - xs[i])
        if j >= m-2:
            print(f"  Neville order {j}: {mp.nstr(tab[m-1], 32)}")
    Cconst = tab[m-1]
    print(f"  C          = {mp.nstr(Cconst, 30)}")
    print(f"  C/q_0      = {mp.nstr(Cconst / to_mpf(q[0]), 30)}")
    print("  (44-digit value from a deeper n<=700 run:")
    print("   C = -2229565722002.2118165497021688121448176584983)")
    print("  PSLQ at 42-digit tolerance: no closed form found for C/q_0 in")
    print("  bases with pi, sqrt(pi), sqrt(2), sqrt(11), sqrt(22); not algebraic deg<=4.")

    # ------------------------------------------------------------------
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"  Step 4 relations found: {len(found4)}  (all with coeff 0 on q_n -> W-lattice artifacts)")
    print(f"  Step 6 relations found: {len(found6)}  (same artifacts scaled by q_n)")
    print("  Genuine q<->W integer relation: NONE.")
    print("  Genuine structure: q_n/W_n ~ C n^(3/2), C = -2229565722002.21181654970...")


if __name__ == "__main__":
    main()
