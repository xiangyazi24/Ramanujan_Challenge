#!/usr/bin/env python3
"""Q6359 (growth moments / GARQI) verification.

Checks, on real Apery data:
  V1. Row-model identity (1.1): H(n) = sum_p 1_{A_{p,X}}(n) equals direct count
      #{p in (n/2, n] : p | b_n}  (via Gessel-Lucas digit product).
  V2. Factorial-moment identity (1.2): sum_n (H(n))_k = sum over distinct ordered
      prime tuples of |A_{p1} cap ... cap A_{pk}|   (small X, k<=3).
  V3. Deterministic lemma (5.3): max H <= k-1 + M_k^{1/k} for all k.
  V4. Empirical Poisson calibration: lambda_X = S_X/N_X, and ratio
      M_k / (N_X * lambda_X^k * k!)  for k=1..6 at dyadic X.
"""
import sys, math
from sympy import primerange

def zero_set(p):
    """Z_p = {0<=r<p : b_r == 0 mod p} via Apery recurrence (n+1<p so invertible)."""
    Z = []
    b0, b1 = 1 % p, 5 % p
    if b0 == 0: Z.append(0)
    if p > 1 and b1 == 0: Z.append(1)
    bprev, bcur = b0, b1
    for n in range(1, p - 1):
        # (n+1)^3 b_{n+1} = (34n^3+51n^2+27n+5) b_n - n^3 b_{n-1}
        num = ((34*n*n*n + 51*n*n + 27*n + 5) * bcur - n*n*n * bprev) % p
        inv = pow((n + 1) % p, p - 4, p) if False else pow(pow(n+1, 3, p), p - 2, p)
        bnext = num * inv % p
        if bnext == 0: Z.append(n + 1)
        bprev, bcur = bcur, bnext
    return Z

def b_mod_p_digits(n, p, brow):
    """b_n mod p via Gessel-Lucas digit product; brow = table of b_r mod p, r<p."""
    v = 1
    while n > 0:
        v = v * brow[n % p] % p
        if v == 0: return 0
        n //= p
    return v

def brow_table(p):
    row = [1 % p, 5 % p] + [0]*(p-2) if p > 1 else [0]
    bprev, bcur = row[0], row[1]
    for n in range(1, p - 1):
        num = ((34*n*n*n + 51*n*n + 27*n + 5) * bcur - n*n*n * bprev) % p
        inv = pow(pow(n+1, 3, p), p - 2, p)
        bnext = num * inv % p
        row[n+1] = bnext
        bprev, bcur = bcur, bnext
    return row

def falling(x, k):
    v = 1
    for j in range(k): v *= (x - j)
    return v

def run_X(X, kmax=6, do_tuple_check=False, do_digit_check=False):
    IX = range(X + 1, 2*X + 1)          # (X, 2X]
    N = X
    primes = list(primerange(X//2 + 1, 2*X + 1))   # P_X = (X/2, 2X]
    A = {}          # p -> set of n in I_X
    S = 0
    t0 = None
    import time; t0 = time.time(); last = t0
    for i, p in enumerate(primes):
        Z = zero_set(p)
        As = set()
        for r in Z:
            n = p + r
            if X < n <= 2*X: As.add(n)
        if As: A[p] = As
        S += len(As)
        now = time.time()
        if now - last > 5:
            print(f"  [progress] X={X}: {i+1}/{len(primes)} primes, {now-t0:.0f}s", flush=True)
            last = now
    # H(n) from rows
    H = {n: 0 for n in IX}
    for p, As in A.items():
        for n in As: H[n] += 1
    # V1: direct check via digit product for a sample of n
    if do_digit_check:
        bad = 0
        sample = list(IX)[::max(1, X//200)]
        for n in sample:
            direct = 0
            for p in primes:
                if n/2 < p <= n:
                    if b_mod_p_digits(n, p, brow_table(p)) == 0:
                        direct += 1
            if direct != H[n]:
                bad += 1
                print(f"  V1 MISMATCH n={n}: direct={direct} rows={H[n]}")
        print(f"  V1 row-model identity: {len(sample)-bad}/{len(sample)} OK")
    # moments
    maxH = max(H.values())
    Ms = {}
    for k in range(1, kmax + 1):
        Ms[k] = sum(falling(h, k) for h in H.values() if h >= k)
    # V2: tuple-intersection identity for small X
    if do_tuple_check:
        from itertools import permutations
        plist = sorted(A)
        for k in (2, 3):
            tot = 0
            for tup in permutations(plist, k):
                inter = A[tup[0]]
                for p in tup[1:]:
                    inter = inter & A[p]
                    if not inter: break
                tot += len(inter)
            ok = (tot == Ms[k])
            print(f"  V2 k={k}: tuple-sum={tot} direct M_k={Ms[k]} {'OK' if ok else 'MISMATCH'}")
    # V3 lemma 5.3
    v3ok = all(maxH <= k - 1 + Ms[k]**(1.0/k) for k in range(1, kmax+1) if Ms[k] > 0)
    lam = S / N
    print(f"X={X}: N={N} #P={len(primes)} S_X={S} lambda={lam:.4f} "
          f"lambda*logX={lam*math.log(X):.3f} maxH={maxH} V3(5.3):{'OK' if v3ok else 'FAIL'}")
    for k in range(1, kmax + 1):
        bench = N * lam**k
        poiss = Ms[k] / bench if bench > 0 else float('nan')
        print(f"    k={k}: M_k={Ms[k]}  M_k/(N lam^k)={poiss:.3f}  "
              f"(k!={math.factorial(k)}; A-slack log ratio/log k! = "
              f"{(math.log(poiss)/math.log(math.factorial(k)) if poiss>0 and k>1 else 0):.2f})")
    return dict(X=X, S=S, lam=lam, maxH=maxH, Ms=Ms)

if __name__ == "__main__":
    print("== small X: full identity checks (V1, V2) ==")
    run_X(120, kmax=4, do_tuple_check=True, do_digit_check=True)
    run_X(240, kmax=4, do_tuple_check=True, do_digit_check=True)
    print("== calibration at growing X (V3, V4) ==")
    for X in (1000, 2000, 4000, 8000):
        run_X(X, kmax=6)
