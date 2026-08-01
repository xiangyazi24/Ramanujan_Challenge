#!/usr/bin/env python3
"""CRON_q6429_verify.py — machine audit of ANSWER Q6429 (q23 symbolic factorization, verdict DEAD).

Verifies, in pure Python over F_p for p in {13,17,29}, every numerical/identity claim:
  1. Apery row mod p, palindromy b_{N-r}=b_r (N=p-1), Z_p sets.
  2. Endpoint-defect identity: X^N Arow(1/X) = Arow + X^N - 1 (truncated row NOT reciprocal).
  3. Coefficient/value transposition: C(g^{-r}) = N * b_r, hence
     Rcorrect = prod_r C(g^{-r}) = N^N * prod b_r  (the TRUE Mellin/diagonal product),
     vs Rrow = prod_{t!=0} Arow(t) (the WRITTEN row resultant). Fourier-dual, NOT equal.
  4. Decisive p=17 witness: Z_17={3,13} nonempty => Rcorrect=0, yet Rrow=2 (a unit).
  5. Chebyshev/Dickson factorization of Rrow:   Rrow = (-1)^m Res(D_m-2,Q-1) Res(D_m+2,Q+1).
  6. Chebyshev factorization of Rcorrect = even*odd = Fourier alias (parity split of r).
  7. CFVZ square: B/q^eps is a perfect square in F_p[X], eps by p mod 24, q=X^2-34X+1.
Expected scalars (from Q6429 §7):
  p=13: Z=[],     Rrow=3,  Rcorrect=4, row_cheb=(1,3),  parity=(1,4)
  p=17: Z=[3,13], Rrow=2,  Rcorrect=0, row_cheb=(1,2),  parity=(2,0)
  p=29: Z=[],     Rrow=16, Rcorrect=5, row_cheb=(9,5),  parity=(4,23)
"""
import sys

def inv(a, p): return pow(a, p - 2, p)

# ---------- polynomial helpers over F_p (coefficient lists, index = degree) ----------
def ptrim(f):
    while f and f[-1] == 0: f.pop()
    return f

def padd(f, g, p):
    n = max(len(f), len(g)); r = [0]*n
    for i in range(n):
        r[i] = ((f[i] if i < len(f) else 0) + (g[i] if i < len(g) else 0)) % p
    return ptrim(r)

def psub(f, g, p):
    n = max(len(f), len(g)); r = [0]*n
    for i in range(n):
        r[i] = ((f[i] if i < len(f) else 0) - (g[i] if i < len(g) else 0)) % p
    return ptrim(r)

def pmul(f, g, p):
    if not f or not g: return []
    r = [0]*(len(f)+len(g)-1)
    for i, a in enumerate(f):
        if a:
            for j, b in enumerate(g):
                r[i+j] = (r[i+j] + a*b) % p
    return ptrim(r)

def pscal(c, f, p): return ptrim([(c*a) % p for a in f])

def peval(f, x, p):
    r = 0
    for a in reversed(f): r = (r*x + a) % p
    return r

def pdivmod(f, g, p):
    f = f[:]; q = [0]*(max(len(f)-len(g)+1, 0))
    ginv = inv(g[-1], p)
    while len(f) >= len(g) and f:
        c = (f[-1]*ginv) % p; d = len(f)-len(g)
        q[d] = c
        for i in range(len(g)):
            f[d+i] = (f[d+i] - c*g[i]) % p
        ptrim(f)
    return ptrim(q), f

def resultant(f, g, p):
    """Res(f,g) over F_p via Euclidean algorithm."""
    f = f[:]; g = g[:]
    if not f or not g: return 0
    res = 1
    while True:
        if len(g) == 1:  # g constant, nonzero
            return (res * pow(g[0], len(f)-1, p)) % p
        _, r = pdivmod(f, g, p)
        if not r: return 0
        # Res(f,g) = (-1)^{deg f * deg g} Res(g,f);  Res(g, r) relation:
        # Res(f,g) = lc(g)^{deg f - deg r} * (-1)^{deg f * deg g} * Res(g, r)
        df, dg, dr = len(f)-1, len(g)-1, len(r)-1
        res = (res * pow(g[-1], df - dr, p) * pow(-1, df*dg, p)) % p
        f, g = g, r

def poly_sqrt(f, p):
    """Return g with g^2 = f over F_p, or None. f assumed nonzero."""
    if (len(f)-1) % 2: return None
    k = (len(f)-1)//2
    lc = f[-1]
    # sqrt of lc mod p (p odd): try s = lc^((p+1)/4) if p%4==3 else Tonelli (small p: brute force)
    s = None
    for c in range(1, p):
        if (c*c) % p == lc: s = c; break
    if s is None: return None
    g = [0]*(k+1); g[k] = s
    for d in range(k-1, -1, -1):
        # coefficient of X^{k+d} in g^2 must equal f[k+d]
        acc = 0
        for i in range(d+1, k+1):
            j = k+d-i
            if 0 <= j <= k and j > d: acc = (acc + g[i]*g[j]) % p
            elif 0 <= j <= k and j == d and i != j: pass
        # g^2 coeff at k+d = sum_{i+j=k+d} g_i g_j = 2 g_d g_k + sum_{i,j>d}
        acc2 = 0
        for i in range(k+1):
            j = k+d-i
            if 0 <= j <= k and i > d and j > d: acc2 = (acc2 + g[i]*g[j]) % p
        g[d] = ((f[k+d] - acc2) * inv(2*g[k] % p, p)) % p
    return g if pmul(g, g, p) == ptrim(f[:]) else None

# ---------- Apery ----------
def apery_row(p):
    b = [0]*p
    b[0] = 1 % p; b[1] = 5 % p
    for n in range(1, p-1):
        Pn = (34*n**3 + 51*n**2 + 27*n + 5) % p
        b[n+1] = ((Pn*b[n] - pow(n, 3, p)*b[n-1]) * inv(pow(n+1, 3, p), p)) % p
    return b

EXPECT = {
    13: dict(Z=[],      Rrow=3,  Rcorrect=4, row_cheb=(1,3), parity=(1,4)),
    17: dict(Z=[3,13],  Rrow=2,  Rcorrect=0, row_cheb=(1,2), parity=(2,0)),
    29: dict(Z=[],      Rrow=16, Rcorrect=5, row_cheb=(9,5), parity=(4,23)),
}

def gen_of(p):
    for g in range(2, p):
        seen, x = set(), 1
        for _ in range(p-1):
            x = (x*g) % p
            if x in seen: break
            seen.add(x)
        if len(seen) == p-1: return g
    raise RuntimeError

fails = 0
def chk(cond, msg, pp):
    global fails
    tag = "OK " if cond else "FAIL"
    if not cond: fails += 1
    print(f"  [{tag}] p={pp}: {msg}")

for p in (13, 17, 29):
    print(f"== p={p} ==")
    N = p-1; m = N//2; g = gen_of(p)
    b = apery_row(p)
    # 1. palindromy incl endpoint b_N (b has indices 0..p-1 = 0..N)
    chk(all(b[N-r] == b[r] for r in range(N+1)), "palindromy b_{N-r}=b_r (full row incl endpoint)", p)
    Z = [r for r in range(N) if b[r] == 0]
    chk(Z == EXPECT[p]['Z'], f"Z_p = {Z} (expected {EXPECT[p]['Z']})", p)
    # 2. endpoint defect: X^N Arow(1/X) = Arow + X^N - 1
    Arow = ptrim([b[r] for r in range(N)])
    Brow = ptrim([b[r] for r in range(N+1)])
    rev_A = ptrim([ (Arow[N-i] if 0 <= N-i < len(Arow) else 0) for i in range(N+1) ])  # X^N Arow(1/X)
    rhs = padd(Arow[:], [p-1] + [0]*(N-1) + [1], p)  # Arow + X^N - 1
    chk(rev_A == rhs, "endpoint defect X^N*Arow(1/X) = Arow + X^N - 1", p)
    # full B reciprocal
    rev_B = ptrim([Brow[N-i] if 0 <= N-i < len(Brow) else 0 for i in range(N+1)])
    chk(rev_B == ptrim(Brow[:]), "full B_p reciprocal", p)
    # 3. row resultant vs corrected
    Rrow = 1
    for t in range(1, p): Rrow = (Rrow * peval(Arow, t, p)) % p
    chk(Rrow == EXPECT[p]['Rrow'], f"Rrow = {Rrow} (expected {EXPECT[p]['Rrow']})", p)
    c = [peval(Arow, pow(g, j, p), p) for j in range(N)]
    C = ptrim(c[:])
    ok_transp = all(peval(C, pow(g, (-r) % (p-1) if False else p-1-(r % (p-1)), p), p) == (N*b[r]) % p for r in range(N))
    # careful: g^{-r} = g^{(p-1-r) mod (p-1)}
    ok_transp = all(peval(C, pow(g, (p-1-r) % (p-1), p), p) == (N*b[r]) % p for r in range(N))
    chk(ok_transp, "transposition C(g^{-r}) = N*b_r for all r", p)
    Rcorrect = 1
    for r in range(N): Rcorrect = (Rcorrect * peval(C, pow(g, (p-1-r) % (p-1), p), p)) % p
    chk(Rcorrect == EXPECT[p]['Rcorrect'], f"Rcorrect = {Rcorrect} (expected {EXPECT[p]['Rcorrect']})", p)
    prod_b = 1
    for r in range(N): prod_b = (prod_b * b[r]) % p
    chk(Rcorrect == (pow(N, N, p) * prod_b) % p, "Rcorrect = N^N * prod b_r", p)
    # decisive witness
    if p == 17:
        chk(Rrow != 0 and Rcorrect == 0, "p=17 witness: Rrow unit but Rcorrect=0 (|Z_17|=2)", p)
        chk(all(x != 0 for x in c), "p=17: A_row has NO zero on F_p^* (value list nonzero)", p)
    # 5. Dickson/Chebyshev
    D = [[2 % p], [0, 1]]
    for j in range(1, m+1):
        D.append(psub(pmul([0,1], D[-1], p), D[-2], p))
    Q = [b[m] % p]
    for j in range(1, m+1):
        Q = padd(Q, pscal(b[m-j], D[j], p), p)
    row_plus  = resultant(psub(D[m], [2], p), psub(Q, [1], p), p)
    row_minus = resultant(padd(D[m], [2], p), padd(Q, [1], p), p)
    chk((row_plus, row_minus) == EXPECT[p]['row_cheb'],
        f"row Chebyshev factors ({row_plus},{row_minus}) (expected {EXPECT[p]['row_cheb']})", p)
    chk(Rrow == (pow(-1, m, p) * row_plus * row_minus) % p,
        "Rrow = (-1)^m * Res(D_m-2,Q-1) * Res(D_m+2,Q+1)", p)
    # 6. corrected-side parity factorization
    chk(all(c[j] == c[(N-j) % N] for j in range(N)), "value vector c_j symmetric under j -> N-j", p)
    S = [c[0] % p]
    for j in range(1, m):
        S = padd(S, pscal(c[j], D[j], p), p)
    d = c[m]
    even_f = resultant(psub(D[m], [2], p), padd(S, [d], p), p)
    odd_f  = resultant(padd(D[m], [2], p), psub(S, [d], p), p)
    chk((even_f, odd_f) == EXPECT[p]['parity'],
        f"parity factors ({even_f},{odd_f}) (expected {EXPECT[p]['parity']})", p)
    ev = 1; od = 1
    for r in range(0, N, 2): ev = (ev * N * b[r]) % p
    for r in range(1, N, 2): od = (od * N * b[r]) % p
    chk(even_f == ev and odd_f == od, "even/odd factors = parity-split products of N*b_r", p)
    chk(Rcorrect == (even_f * odd_f) % p, "Rcorrect = even*odd", p)
    # 7. CFVZ square
    q = [1, (-34) % p, 1]
    eps = 0 if p % 24 in (1, 5, 7, 11) else 1
    quot = Brow[:]
    for _ in range(eps):
        quot, rem = pdivmod(quot, q, p)
        chk(rem == [], f"q^{eps} divides B_p exactly", p)
    sq = poly_sqrt(quot, p)
    chk(sq is not None, f"B_p/q^{eps} is a perfect square in F_p[X] (eps={eps})", p)

print()
print("ALL PASS" if fails == 0 else f"{fails} FAILURES")
sys.exit(0 if fails == 0 else 1)
