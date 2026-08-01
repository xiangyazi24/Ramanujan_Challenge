#!/usr/bin/env python3
"""(SG) diagnostic: Lambda(p,H) = sigma_1(A_H)^2 / H for the lag-incidence matrix
A_H[r,h] = 1[b_{r+h} = b_r], h < H. Bounded/polylog => (SG) plausible => 3/2 alive.
Also Q = |A 1|^2/H^2 (all-ones form), D = nnz/H, and top-vector lag concentration."""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import svds

def apery_mod(p):
    b = [1, 5]
    for n in range(1, p-1):
        b.append(((2*n+1)*(17*n*n+17*n+5) % p * b[n] - pow(n,3,p)*b[n-1]) * pow(pow(n+1,3,p), p-2, p) % p)
    return b

for p in (1009, 4001, 16001, 40009, 100003):
    b = apery_mod(p)
    H = int(p ** 0.5)
    rows, cols = [], []
    for h in range(1, H):
        bh = b[h:]
        for r in range(len(b) - h):
            if b[r+h] == b[r]:
                rows.append(r); cols.append(h-1)
    A = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(b), H-1))
    if A.nnz < 3:
        print(f"p={p}: too sparse ({A.nnz})"); continue
    k = 1
    try:
        u, s, vt = svds(A, k=1, which='LM')
    except Exception as e:
        print(f"p={p}: svds fail {e}"); continue
    sigma = s[0]
    Lam = sigma**2 / H
    ones = np.ones(H-1)
    Q = np.linalg.norm(A @ ones)**2 / H**2
    D = A.nnz / H
    v = np.abs(vt[0]); v = v / v.sum()
    top3 = np.argsort(-v)[:3]
    print(f"p={p} H={H}: nnz={A.nnz} sigma1^2={sigma**2:.1f} Lambda={Lam:.3f} Q={Q:.3f} D={D:.2f} "
          f"top lags h={[int(t)+1 for t in top3]} weights={[round(float(v[t]),3) for t in top3]}")
