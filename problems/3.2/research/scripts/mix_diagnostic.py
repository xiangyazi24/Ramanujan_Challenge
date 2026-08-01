#!/usr/bin/env python3
"""(MIX) diagnostic: average root count of the mixed determinant
C_{h,u;t}(x) = det[ T_t(x) w_h(x) | w_u(x+t) ], w_d = (B_d, 1-A_d)^T.
If avg roots per triple stays O(1) as H grows => (MIX) plausible => 3/2 alive.
If it grows ~ H => two-base filter saturates degree bound => need third base/Fourier."""
import numpy as np

def run(p, HMAX, G):
    # tables A[d][x], B[d][x] for d <= HMAX+G
    D = HMAX + G
    A = np.zeros((D+1, p), dtype=np.int64); B = np.zeros((D+1, p), dtype=np.int64)
    A[0] = 1; B[0] = 0
    inv = np.zeros(p, dtype=np.int64)
    inv[1:] = np.array([pow(int(i), p-2, p) for i in range(1, p)], dtype=np.int64)
    xs = np.arange(p, dtype=np.int64)
    aa = np.zeros(p, dtype=np.int64); bb = np.zeros(p, dtype=np.int64)
    den = (xs+1) % p
    ok = den != 0
    d3inv = np.zeros(p, dtype=np.int64)
    d3inv[ok] = inv[den[ok]]**3 % p
    aa[ok] = (2*xs[ok]+1)*(17*xs[ok]**2+17*xs[ok]+5) % p * d3inv[ok] % p
    bb[ok] = (-(xs[ok]**3)) % p * d3inv[ok] % p
    # A_d(x) recursion: A_d(x) = a(x+d-1) A_{d-1}(x) + beta(x+d-1) A_{d-2}(x)
    Am1 = np.ones(p, dtype=np.int64); Bm1 = np.zeros(p, dtype=np.int64)  # d=0
    sh = lambda arr, k: np.roll(arr, -k)
    A[1] = sh(aa,0); B[1] = sh(bb,0)
    for d in range(2, D+1):
        av = sh(aa, d-1); bv = sh(bb, d-1)
        A[d] = (av*A[d-1] + bv*A[d-2]) % p
        B[d] = (av*B[d-1] + bv*B[d-2]) % p
    # T_t entries: [[A_t, B_t],[A_{t-1}, B_{t-1}]]
    stats = {}
    for H in (8, 12, 16, HMAX):
        tot_roots = 0; n_triples = 0
        for t in range(1, G+1):
            for h in range(1, H):
                wh1 = B[h]; wh2 = (1 - A[h]) % p
                v1 = (A[t]*wh1 + B[t]*wh2) % p
                v2 = (A[t-1]*wh1 + B[t-1]*wh2) % p if t >= 1 else wh1
                for u in range(1, H):
                    wu1 = sh(B[u], t); wu2 = sh((1 - A[u]) % p, t)
                    C = (v1*wu2 - v2*wu1) % p
                    # exclude wrapping x region and identically-zero (partial matching)
                    body = C[: p - (max(h, t+u) + 2)]
                    z = int(np.count_nonzero(body == 0))
                    if z > len(body) // 2:   # identically zero (bad alignment)
                        continue
                    tot_roots += z; n_triples += 1
        stats[H] = (tot_roots, n_triples, tot_roots/max(1,n_triples))
    return stats

for p in (1009, 5003):
    st = run(p, 20, 20)
    print(f"p={p}: " + "  ".join(f"[H={H}: avg roots/triple={a:.2f} ({r} roots / {n} triples)]"
          for H,(r,n,a) in st.items()))
