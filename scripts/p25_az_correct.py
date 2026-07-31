#!/usr/bin/env python3
"""
Problem 2.5 (Catalan G): AZ certificate search with the CORRECT carrier phase.

    Phi(u,v) = phi(u)*phi(v),   phi(u) = (1+u)(u+2)/u = 3 + u + 2/u

AZ equation for a row-valued numerator P_n in Mat_{1x3}:

    Phi(u,v) * P_{n+1}(u,v,t) - P_n(u,v,t) * M_H(n) = D_u(A_n) + D_v(B_n) + D_t(C_n)

Covariant derivatives (connection from log-derivative of the AZ measure):
    D_u = d_u + n(u^2-2)/(u(u+1)(u+2)) - 1/u
    D_v = d_v + n(v^2-2)/(v(v+1)(v+2)) - 1/v
    D_t = d_t + eps/t - 2t/(1+t^2)

M_H(n) = M(n)/delta_H(n), delta_H(n) = -2(n+2)^2(n+3)^2(2n+5)(2n+7)^2.

KEY FIX (v2): certificate supports must be WIDER than P to absorb the extra
degrees introduced by Phi and the geometric factor G.

After clearing by delta_H * G:
  - Phi*G has u-degree 0..4  -->  Phi*G*P at u contributes [a_min, a_max+4]
  - G*M*P at u contributes [a_min+1, a_max+3]
  - D_u(A) after clearing contributes [a_A_min, a_A_max+2]
  => need a_A_max >= a_P_max + 2, a_A_min <= a_P_min
  Similarly for v/B.  For t: need c_C_max >= c_P_max + 1.
"""

import argparse
import random
import sys
import time
from fractions import Fraction
from math import gcd, isqrt

import numpy as np

PRIME1 = 2147483647          # 2^31 - 1
PRIME2 = 2147483629          # also prime

# --------------------------------------------------------------------------
# Laurent polynomial dictionaries: {(a,b,c): coeff}.  p=0 means exact ints.
# --------------------------------------------------------------------------

def lmul(d1, d2, p=0):
    out = {}
    for k1, v1 in d1.items():
        for k2, v2 in d2.items():
            k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])
            w = out.get(k, 0) + v1 * v2
            out[k] = w % p if p else w
    return {k: v for k, v in out.items() if v}


def ladd(*ds, p=0):
    out = {}
    for d in ds:
        for k, v in d.items():
            w = out.get(k, 0) + v
            out[k] = w % p if p else w
    return {k: v for k, v in out.items() if v}


def lscale(d, s, p=0):
    if p:
        s %= p
    out = {}
    for k, v in d.items():
        w = v * s % p if p else v * s
        if w:
            out[k] = w
    return out


def lshift(d, m):
    a, b, c = m
    return {(k[0] + a, k[1] + b, k[2] + c): v for k, v in d.items()}


def tag(d, comp):
    """Attach component index: {(comp,a,b,c): v}"""
    return {(comp,) + k: v for k, v in d.items()}


# fixed geometric polynomials -----------------------------------------------
GU   = {(3,0,0):1, (2,0,0):3, (1,0,0):2}            # u(u+1)(u+2)
GV   = {(0,3,0):1, (0,2,0):3, (0,1,0):2}            # v(v+1)(v+2)
GT   = {(0,0,3):1, (0,0,1):1}                       # t(1+t^2)
U1U2 = {(2,0,0):1, (1,0,0):3, (0,0,0):2}            # (u+1)(u+2)
V1V2 = {(0,2,0):1, (0,1,0):3, (0,0,0):2}            # (v+1)(v+2)
U2M2 = {(2,0,0):1, (0,0,0):-2}                      # u^2 - 2
V2M2 = {(0,2,0):1, (0,0,0):-2}                      # v^2 - 2
ONE_T2 = {(0,0,0):1, (0,0,2):1}                     # 1 + t^2

# Phi(u,v) * G  -- after clearing, this is a polynomial (no negative powers)
# Phi = phi(u)*phi(v) = [(1+u)(u+2)/u]*[(1+v)(v+2)/v]
# Phi*G = (1+u)^2*(u+2)^2 * (1+v)^2*(v+2)^2 * t(1+t^2)
def _compute_phig():
    # (1+u)^2 = u^2+2u+1, (u+2)^2 = u^2+4u+4
    u1sq = {(2,0,0):1, (1,0,0):2, (0,0,0):1}
    u2sq = {(2,0,0):1, (1,0,0):4, (0,0,0):4}
    phiu = lmul(u1sq, u2sq)  # u^4 + 6u^3 + 13u^2 + 12u + 4
    v1sq = {(0,2,0):1, (0,1,0):2, (0,0,0):1}
    v2sq = {(0,2,0):1, (0,1,0):4, (0,0,0):4}
    phiv = lmul(v1sq, v2sq)
    return lmul(lmul(phiu, phiv), GT)

PHIG = _compute_phig()
G    = lmul(lmul(GU, GV), GT)
GVGT = lmul(GV, GT)
GUGT = lmul(GU, GT)
GUGV = lmul(GU, GV)


# --------------------------------------------------------------------------
# The CMF matrix M(n) and normalizer delta_H(n)
# --------------------------------------------------------------------------

def Mmat(n):
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


def deltaH(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2


# --------------------------------------------------------------------------
# Build per-unknown "terms": lists of (scalar_id, laurent-vector dict).
# Now with SEPARATE supports for P, A, B, C.
# --------------------------------------------------------------------------

def make_sup(umin, umax, vmin, vmax, cmin, cmax):
    return [(a, b, c) for a in range(umin, umax+1)
                      for b in range(vmin, vmax+1)
                      for c in range(cmin, cmax+1)]


def build_terms(dP, dC, eps_mode, eps_val, supP, supA, supB, supC,
                uv_exp=-1):
    """Return (unknowns, terms) with SEPARATE supports for each block.

    supP, supA, supB, supC: lists of (a,b,c) tuples.
    """
    blocks = [('P', dP, supP), ('A', dC, supA), ('B', dC, supB), ('C', dC, supC)]
    if eps_mode == 'free':
        blocks.append(('C2', dC, supC))
    connT = {(0,0,0): eps_val, (0,0,2): eps_val - 2}  # eps(1+t^2) - 2t^2
    unknowns, terms = [], []
    for blk, deg, sup in blocks:
        for k in range(deg + 1):
            for i in range(3):
                for m in sup:
                    a, b, c = m
                    tl = []
                    if blk == 'P':
                        tl.append((('phi', k), tag(lshift(PHIG, m), i)))
                        Gs = lshift(G, m)
                        for j in range(3):
                            tl.append((('m', k, i, j), tag(Gs, j)))
                    elif blk == 'A':
                        # D_u cleared:
                        # GVGT * [GU*d_u(A) + (n*U2M2 + uv_exp*U1U2)*A].
                        # The historical carrier has uv_exp=-1.  The value -2
                        # is the ordinary-form carrier whose P-period extracts
                        # the Delannoy constant term when P contains a factor uv.
                        # GU * d_u(u^a) = a * u(u+1)(u+2) * u^{a-1}
                        base = ladd(lscale(lshift(GU, (a-1, b, c)), a),
                                    lscale(lshift(U1U2, m), uv_exp))
                        tl.append((('c0', k), tag(lmul(GVGT, base), i)))
                        tl.append((('c1', k), tag(lmul(GVGT, lshift(U2M2, m)), i)))
                    elif blk == 'B':
                        base = ladd(lscale(lshift(GV, (a, b-1, c)), b),
                                    lscale(lshift(V1V2, m), uv_exp))
                        tl.append((('c0', k), tag(lmul(GUGT, base), i)))
                        tl.append((('c1', k), tag(lmul(GUGT, lshift(V2M2, m)), i)))
                    elif blk == 'C':
                        base = ladd(lscale(lshift(GT, (a, b, c-1)), c),
                                    lshift(connT, m))
                        tl.append((('c0', k), tag(lmul(GUGV, base), i)))
                    elif blk == 'C2':
                        # eps * (1+t^2) * C  contribution, with eps treated as free
                        tl.append((('c0', k), tag(lmul(GUGV, lshift(ONE_T2, m)), i)))
                    unknowns.append((blk, k, i, m))
                    terms.append([(sid, d) for sid, d in tl if d])
    return unknowns, terms


def scalars(n0, maxk, p):
    dH = deltaH(n0) % p
    Mn = [[x % p for x in row] for row in Mmat(n0)]
    sc = {}
    for k in range(maxk + 1):
        nk = pow(n0, k, p)
        sc[('phi', k)] = dH * pow(n0 + 1, k, p) % p
        sc[('c0', k)] = (-dH * nk) % p
        sc[('c1', k)] = (-dH * nk % p) * (n0 % p) % p
        for i in range(3):
            for j in range(3):
                sc[('m', k, i, j)] = (-nk * Mn[i][j]) % p
    return sc


def gen_sample_triplets(terms, sc, rowidx, p):
    rows, cols, vals = [], [], []
    for col, tl in enumerate(terms):
        acc = {}
        for sid, d in tl:
            s = sc[sid]
            if s == 0:
                continue
            for key, v in d.items():
                acc[key] = (acc.get(key, 0) + s * v) % p
        for key, v in acc.items():
            if v:
                rows.append(rowidx[key])
                cols.append(col)
                vals.append(v)
    return (np.array(rows, dtype=np.int64),
            np.array(cols, dtype=np.int64),
            np.array(vals, dtype=np.int64))


# --------------------------------------------------------------------------
# Mod-p linear algebra (numpy int64; p < 2^31 so products fit in int64)
# --------------------------------------------------------------------------

def rref_mod(A, p):
    A = A % p
    m, n = A.shape
    piv, r = [], 0
    for c in range(n):
        if r >= m:
            break
        nz = np.flatnonzero(A[r:, c])
        if nz.size == 0:
            continue
        i = r + int(nz[0])
        if i != r:
            A[[r, i]] = A[[i, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        f = A[:, c].copy()
        f[r] = 0
        nzr = np.flatnonzero(f)
        if nzr.size:
            A[nzr] = (A[nzr] - np.outer(f[nzr], A[r])) % p
        piv.append(c)
        r += 1
    return A[:r], piv


def nullspace_mod(A, p):
    R, piv = rref_mod(A, p)
    n = A.shape[1]
    pivset = set(piv)
    free = [c for c in range(n) if c not in pivset]
    if not free:
        return np.zeros((n, 0), dtype=np.int64)
    V = np.zeros((n, len(free)), dtype=np.int64)
    for jf, f in enumerate(free):
        V[f, jf] = 1
        for ri, c in enumerate(piv):
            V[c, jf] = (-int(R[ri, f])) % p
    return V


def matvec_mod(rows, cols, vals, v, nrows, p):
    """Sparse (rows,cols,vals) @ v mod p, overflow-safe via 16-bit limb split."""
    vlo = v & 0xFFFF
    vhi = v >> 16
    rlo = np.zeros(nrows, dtype=np.int64)
    rhi = np.zeros(nrows, dtype=np.int64)
    np.add.at(rlo, rows, vals * vlo[cols])
    np.add.at(rhi, rows, vals * vhi[cols])
    return (rlo % p + ((rhi % p) << 16)) % p


def dense_matvec_mod(A, v, p):
    """A @ v mod p for int64 A,v < p, overflow-safe (limb split on v)."""
    vlo = v & 0xFFFF
    vhi = v >> 16
    return ((A @ vlo) % p + (((A @ vhi) % p) << 16)) % p


# --------------------------------------------------------------------------
# Rational reconstruction
# --------------------------------------------------------------------------

def ratrec(a, p, bound=None):
    a %= p
    if bound is None:
        bound = isqrt(p // 2)
    g0, g1 = p, a
    s0, s1 = 0, 1
    while g1 > bound:
        q = g0 // g1
        g0, g1 = g1, g0 - q * g1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound or gcd(g1, abs(s1)) != 1:
        return None
    num, den = (g1, s1) if s1 > 0 else (-g1, -s1)
    if num % p != a * den % p:
        return None
    return Fraction(num, den)


def crt(a1, a2, p1, p2):
    m = pow(p1, -1, p2)
    return (a1 + p1 * ((a2 - a1) * m % p2)) % (p1 * p2)


# --------------------------------------------------------------------------
# One full solve for a given (dP, dC, eps_mode) and prime p
# --------------------------------------------------------------------------

def solve_config(dP, dC, eps_mode, p, seed=0,
                 Pu=(-2,2), Pv=(-2,2), Pt=(0,2),
                 Au=None, Av=None, At=None,
                 Bu=None, Bv=None, Bt=None,
                 Cu=None, Cv=None, Ct=None,
                 uv_exp=-1,
                 verbose=True):
    """
    Solve the AZ linear system mod p.

    Support ranges for each block default to the necessary widths:
      A_u = [Pu[0], Pu[1]+2], etc.
    """
    t0 = time.time()
    rng = np.random.default_rng(seed + p)
    if eps_mode == 'generic':
        eps_val = random.Random(seed + p).randrange(10**6, p - 10**6)
    elif eps_mode == 'free':
        eps_val = 0
    else:
        eps_val = int(eps_mode)

    # Default certificate supports: widen to absorb Phi*G*P terms
    if Au is None: Au = (Pu[0], Pu[1] + 2)
    if Av is None: Av = Pv
    if At is None: At = Pt
    if Bu is None: Bu = Pu
    if Bv is None: Bv = (Pv[0], Pv[1] + 2)
    if Bt is None: Bt = Pt
    if Cu is None: Cu = Pu
    if Cv is None: Cv = Pv
    if Ct is None: Ct = (Pt[0], Pt[1] + 1)

    supP = make_sup(*Pu, *Pv, *Pt)
    supA = make_sup(*Au, *Av, *At)
    supB = make_sup(*Bu, *Bv, *Bt)
    supC = make_sup(*Cu, *Cv, *Ct)

    unknowns, terms = build_terms(
        dP, dC, eps_mode, eps_val, supP, supA, supB, supC, uv_exp=uv_exp)
    N = len(unknowns)
    maxk = max(dP, dC)

    # global row map (sample-independent monomial support)
    keys = set()
    for tl in terms:
        for _, d in tl:
            keys.update(d.keys())
    rowkeys = sorted(keys)
    rowidx = {k: idx for idx, k in enumerate(rowkeys)}
    nrows = len(rowkeys)

    nsamples = 8 + maxk + 3
    samples = list(range(1, nsamples + 1))

    # compressed accumulator: nc must be >= N so the random projection
    # preserves the full system's rank w.h.p.  (When nrows < N, the old
    # nc = min(N+80, nrows+80) was too small and produced spurious nullity.)
    nc = N + 80
    acc = np.zeros((nc, N), dtype=np.int64)
    trips = []
    for n0 in samples:
        sc = scalars(n0, maxk, p)
        rows, cols, vals = gen_sample_triplets(terms, sc, rowidx, p)
        trips.append((n0, rows, cols, vals))
        Z = np.zeros((nrows, N), dtype=np.int64)
        Z[rows, cols] = vals
        R = rng.integers(1, 1 << 16, size=(nc, nrows)).astype(np.float64)
        Zlo = (Z & 0xFFFF).astype(np.float64)
        Zhi = (Z >> 16).astype(np.float64)
        lo = np.rint(R @ Zlo).astype(np.int64) % p
        hi = np.rint(R @ Zhi).astype(np.int64) % p
        acc = (acc + lo + ((hi << 16) % p)) % p

    V = nullspace_mod(acc, p)
    nullity = V.shape[1]

    # verify a random nullspace combination against the FULL system
    ok = True
    if nullity:
        coeffs = rng.integers(1, p, size=nullity).astype(np.int64)
        vtest = dense_matvec_mod(V, coeffs, p)
        for n0, rows, cols, vals in trips:
            r = matvec_mod(rows, cols, vals, vtest, nrows, p)
            if np.any(r):
                ok = False
                break

    # dimension of realizable P-space
    Pidx = np.array([idx for idx, u in enumerate(unknowns) if u[0] == 'P'],
                    dtype=np.int64)
    Pdim, Pbasis = 0, None
    if nullity:
        Pmat = V[Pidx, :] % p
        Rp, pivp = rref_mod(Pmat.T.copy(), p)
        Pdim = len(pivp)
        Pbasis = Rp

    nP = len(supP) * 3 * (dP + 1)
    nA = len(supA) * 3 * (dC + 1)
    nB = len(supB) * 3 * (dC + 1)
    nC = len(supC) * 3 * (dC + 1)
    if verbose:
        print(f"  [p={p}] dP={dP} dC={dC} eps={eps_mode}: "
              f"N={N} (P:{nP} A:{nA} B:{nB} C:{nC}), "
              f"rows/sample={nrows}, samples={nsamples}, "
              f"nullity={nullity}, P-dim={Pdim}, "
              f"verify={'OK' if ok else 'FAIL'}, {time.time()-t0:.1f}s")
        sys.stdout.flush()

    return {
        'unknowns': unknowns, 'terms': terms,
        'supP': supP, 'supA': supA, 'supB': supB, 'supC': supC,
        'p': p,
        'rowidx': rowidx, 'nrows': nrows, 'trips': trips,
        'V': V, 'nullity': nullity, 'verify_ok': ok,
        'Pidx': Pidx, 'Pdim': Pdim, 'Pbasis': Pbasis,
        'eps_val': eps_val, 'dP': dP, 'dC': dC, 'eps_mode': eps_mode,
    }


def full_vector_for_P(res, ptarget):
    """Find nullspace combination whose P-part equals ptarget (mod p)."""
    p = res['p']
    V, Pidx = res['V'], res['Pidx']
    Pmat = V[Pidx, :] % p
    aug = np.concatenate([Pmat, ptarget.reshape(-1, 1)], axis=1)
    R, piv = rref_mod(aug, p)
    ncol = Pmat.shape[1]
    if ncol in piv:
        return None
    x = np.zeros(ncol, dtype=np.int64)
    for ri, c in enumerate(piv):
        x[c] = int(R[ri, ncol])
    return dense_matvec_mod(V, x, p)


def describe_P(res, prow):
    """Human-readable description of a P-space basis vector."""
    p = res['p']
    unknowns = res['unknowns']
    Punk = [u for u in unknowns if u[0] == 'P']
    out = [[], [], []]
    for idx, val in enumerate(prow):
        if val == 0:
            continue
        _, k, i, (a, b, c) = Punk[idx]
        r = ratrec(int(val), p)
        cs = str(r) if r is not None else f"{int(val)} (mod p)"
        mono = ""
        if k:
            mono += f"n^{k}*" if k > 1 else "n*"
        for s, e in (("u", a), ("v", b), ("t", c)):
            if e:
                mono += f"{s}^{e}" if abs(e) != 1 else (f"{s}" if e == 1 else f"{s}^-1")
        out[i].append(f"{cs}" + (f"*{mono}" if mono else ""))
    return ["  comp %d:  %s" % (i, "  +  ".join(out[i]) if out[i] else "0")
            for i in range(3)]


# --------------------------------------------------------------------------
# Exact verification over Z
# --------------------------------------------------------------------------

def exact_residual(unknowns_terms, xexact, eps_exact, n0, dP, dC):
    """Evaluate cleared equation with exact integer/Fraction unknown vector at n."""
    unknowns, terms = unknowns_terms
    dH = deltaH(n0)
    Mn = Mmat(n0)
    maxk = max(dP, dC)
    sc = {}
    for k in range(maxk + 1):
        sc[('phi', k)] = dH * (n0 + 1) ** k
        sc[('c0', k)] = -dH * n0 ** k
        sc[('c1', k)] = -dH * n0 ** (k + 1)
        for i in range(3):
            for j in range(3):
                sc[('m', k, i, j)] = -(n0 ** k) * Mn[i][j]
    acc = {}
    for col, tl in enumerate(terms):
        x = xexact[col]
        if x == 0:
            continue
        for sid, d in tl:
            s = sc[sid] * x
            if s == 0:
                continue
            for key, v in d.items():
                acc[key] = acc.get(key, 0) + s * v
    return {k: v for k, v in acc.items() if v}


# --------------------------------------------------------------------------
# Self-test of the equation builder against sympy
# --------------------------------------------------------------------------

def selftest():
    import sympy as sp
    print("Self-test: comparing cleared-equation builder against sympy ...")
    u, v, t = sp.symbols('u v t')
    supP = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (0, 1)]
    supA = [(a, b, c) for a in (-1, 0, 1, 2, 3) for b in (-1, 0, 1) for c in (0, 1)]
    supB = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1, 2, 3) for c in (0, 1)]
    supC = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (0, 1, 2)]
    dP = dC = 1
    eps = 1
    rnd = random.Random(42)
    unknowns, terms = build_terms(dP, dC, str(eps), eps, supP, supA, supB, supC)
    x = [rnd.randrange(-5, 6) for _ in unknowns]
    n0 = 3

    R = exact_residual((unknowns, terms), x, eps, n0, dP, dC)

    # sympy reference
    def blk_poly(blk, nval):
        rows = [sp.Integer(0)] * 3
        for idx, (b, k, i, (a, bb, c)) in enumerate(unknowns):
            if b != blk or x[idx] == 0:
                continue
            rows[i] += x[idx] * nval**k * u**a * v**bb * t**c
        return sp.Matrix(1, 3, rows)

    n = sp.Integer(n0)
    Pn, Pn1 = blk_poly('P', n), blk_poly('P', n + 1)
    A, B, C = blk_poly('A', n), blk_poly('B', n), blk_poly('C', n)
    MH = sp.Matrix(Mmat(n0)) / deltaH(n0)
    Phi = (1 + u) * (u + 2) / u * (1 + v) * (v + 2) / v

    def Du(F): return F.diff(u) + (n * (u**2 - 2) / (u * (u + 1) * (u + 2)) - 1 / u) * F
    def Dv(F): return F.diff(v) + (n * (v**2 - 2) / (v * (v + 1) * (v + 2)) - 1 / v) * F
    def Dt(F): return F.diff(t) + (sp.Integer(eps) / t - 2 * t / (1 + t**2)) * F

    E = Phi * Pn1 - Pn * MH - Du(A) - Dv(B) - Dt(C)
    Gsym = u * (u + 1) * (u + 2) * v * (v + 1) * (v + 2) * t * (1 + t**2)
    for comp in range(3):
        lhs = sp.expand(sp.cancel(deltaH(n0) * Gsym * E[comp]))
        ours = sum(sp.Integer(val) * u**k[1] * v**k[2] * t**k[3]
                   for k, val in R.items() if k[0] == comp)
        diff = sp.expand(lhs - ours)
        if diff != 0:
            print(f"  SELF-TEST FAILED at component {comp}:")
            print("  diff =", diff)
            return False
    print("  Self-test PASSED (all 3 components match sympy).")
    return True


# --------------------------------------------------------------------------
# Main search
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dmax', type=int, default=2)
    ap.add_argument('--eps', type=str, default='free,generic,0,1,-1,2')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--dcplus', type=int, default=0,
                    help='certificate degree = dP + dcplus')
    ap.add_argument('--pu', type=str, default='-2,2', help='P u-range')
    ap.add_argument('--pv', type=str, default='-2,2', help='P v-range')
    ap.add_argument('--pt', type=str, default='0,2', help='P t-range')
    ap.add_argument('--uv-exp', type=int, default=-1,
                    help='u and v monomial exponents in the ordinary-form carrier')
    args = ap.parse_args()

    if args.selftest:
        if not selftest():
            sys.exit(1)
        return

    Pu = tuple(int(x) for x in args.pu.split(','))
    Pv = tuple(int(x) for x in args.pv.split(','))
    Pt = tuple(int(x) for x in args.pt.split(','))

    eps_modes = args.eps.split(',')
    print("=" * 78)
    print("P2.5 AZ certificate search (v2: correct certificate support widths)")
    print(f"  Phase: Phi(u,v) = phi(u)*phi(v),  phi(u) = (1+u)(u+2)/u")
    print(f"  Ordinary-form carrier exponents: u^{args.uv_exp} v^{args.uv_exp}")
    print(f"  P support: u in {Pu}, v in {Pv}, t in {Pt}")
    print(f"  A support: u in ({Pu[0]},{Pu[1]+2}), v in {Pv}, t in {Pt}")
    print(f"  B support: u in {Pu}, v in ({Pv[0]},{Pv[1]+2}), t in {Pt}")
    print(f"  C support: u in {Pu}, v in {Pv}, t in ({Pt[0]},{Pt[1]+1})")
    print("=" * 78)
    sys.stdout.flush()

    hits = []
    for dP in range(0, args.dmax + 1):
        dC = dP + args.dcplus
        print(f"\n--- n-degree dP={dP}, dC={dC} ---")
        sys.stdout.flush()
        for em in eps_modes:
            res = solve_config(dP, dC, em, PRIME1, seed=args.seed,
                               Pu=Pu, Pv=Pv, Pt=Pt, uv_exp=args.uv_exp)
            if not res['verify_ok']:
                print("    WARNING: compressed nullspace failed full verification; "
                      "retrying with new seed")
                sys.stdout.flush()
                res = solve_config(dP, dC, em, PRIME1, seed=args.seed + 1000,
                                   Pu=Pu, Pv=Pv, Pt=Pt, uv_exp=args.uv_exp)
            if res['Pdim'] > 0:
                # confirm with second prime
                res2 = solve_config(dP, dC, em, PRIME2, seed=args.seed,
                                    Pu=Pu, Pv=Pv, Pt=Pt, uv_exp=args.uv_exp)
                agree = (res2['Pdim'] == res['Pdim']
                         and res2['nullity'] == res['nullity'])
                print(f"    >>> HIT: P-dim={res['Pdim']} (eps mode {em}); "
                      f"second prime agreement: {agree} "
                      f"(p2 nullity={res2['nullity']}, P-dim={res2['Pdim']})")
                sys.stdout.flush()
                print(f"\n  === P-space basis (canonical RREF, {res['Pdim']} vector(s)) ===")
                for bi in range(min(res['Pdim'], 5)):
                    print(f"  P basis vector #{bi}:")
                    for line in describe_P(res, res['Pbasis'][bi]):
                        print(line)
                sys.stdout.flush()
                hits.append((dP, em, res, res2))
        if hits:
            print(f"\nFound solution(s) at dP={dP}; stopping degree escalation.")
            sys.stdout.flush()
            break

    print("\n" + "=" * 78)
    if hits:
        print("RESULT: non-trivial AZ certificate(s) with P != 0 FOUND:")
        for dP, em, res, _ in hits:
            print(f"  dP={dP}, eps mode={em}, P-dim={res['Pdim']}, "
                  f"nullity={res['nullity']}")
    else:
        print("RESULT: NO non-trivial P != 0 certificate in the searched space")
        print(f"  (P support: u in {Pu}, v in {Pv}, t in {Pt}; "
              f"n-degree <= {args.dmax};")
        print(f"   eps modes {eps_modes}; incl. 'free' relaxation covering ALL eps).")
    print("=" * 78)
    sys.stdout.flush()


if __name__ == '__main__':
    main()
