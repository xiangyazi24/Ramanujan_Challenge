#!/usr/bin/env python3
"""CRON_stratify_t34.py — 异常分层实验 (Q6511 §9.3 第5条 + AO.7 盲点) + T3/T4 模型对表.

对每个素数:
  T3/p, T4/p (模型预言: T3/p->4, T4/p->5 [Poisson(1/2)+镜像, E[(2X)_3]=8*1/8+12*1/4=4, E[(2X)_4]=10 每半格点 => /p 取半 => 5])
  T3 分层:
    - mirror 层: 三元组含镜像对 (n_i+n_j == p-1 mod p-1) vs 纯生日
    - 奇偶层: 三元组 gap (h,k) = (n2-n1, n3-n2) 的奇偶类 (ee/eo/oe/oo)
    - 小gap层: 三元组含 gap<=14 (结构根区) 的份额, 按 (-51|p) 分类
  全部与镜像随机基准对照 (3 draws).
判据: 若某薄层 data/base 比值显著偏离 1, 即 Q6511 所言"总均值密度为零但高影响的异常族"。
"""
import sys, time
from collections import defaultdict
import numpy as np

def primes_in(a, b):
    sieve = np.ones(b + 1, dtype=bool); sieve[:2] = False
    for i in range(2, int(b ** 0.5) + 1):
        if sieve[i]: sieve[i * i::i] = False
    return [int(x) for x in np.nonzero(sieve)[0] if x >= a]

def orbit_keys(p):
    N = p - 1
    b = [0]*N; c = [0]*N
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    for n in range(1, N - 1):
        Pn = (34*n*n*n + 51*n*n + 27*n + 5) % p
        n3 = (n*n*n) % p
        inv = pow((n+1)**3 % p, p-2, p)
        b[n+1] = ((Pn*b[n] - n3*b[n-1]) % p) * inv % p
        c[n+1] = ((Pn*c[n] - n3*c[n-1]) % p) * inv % p
    key = [0]*N
    for n in range(N):
        key[n] = (b[n] * pow(c[n], p-2, p)) % p if c[n] != 0 else p
    return key

def mirror_random_key(p, rng):
    N = p - 1
    key = [-1]*N
    for n in range(N):
        if key[n] >= 0: continue
        m = (p - 1 - n) % N
        v = int(rng.integers(0, p + 1))
        key[n] = v; key[m] = v
    return key

def stats(key, p):
    N = p - 1
    cls = defaultdict(list)
    for n, v in enumerate(key): cls[v].append(n)
    T3 = 0; T4 = 0
    strat = dict(mir3=0, pure3=0, ee=0, eo=0, oe=0, oo=0, smallgap=0, total3=0)
    for v, lst in cls.items():
        m = len(lst)
        if m >= 3: T3 += m*(m-1)*(m-2)
        if m >= 4: T4 += m*(m-1)*(m-2)*(m-3)
        if m < 3: continue
        lst.sort()
        L = len(lst)
        for i in range(L):
            for j in range(i+1, L):
                for k in range(j+1, L):
                    n1, n2, n3 = lst[i], lst[j], lst[k]
                    strat['total3'] += 1
                    has_mirror = any(((a + b) % N) == (p - 1) % N for a, b in ((n1,n2),(n1,n3),(n2,n3)))
                    if has_mirror: strat['mir3'] += 1
                    else: strat['pure3'] += 1
                    h, kk = n2-n1, n3-n2
                    strat['ee' if h%2==0 and kk%2==0 else 'eo' if h%2==0 else 'oe' if kk%2==0 else 'oo'] += 1
                    if h <= 14 or kk <= 14 or (n3-n1) <= 14: strat['smallgap'] += 1
    return T3, T4, strat

def leg51(p):
    return pow(-51 % p, (p-1)//2, p) == 1

def main():
    calib = len(sys.argv) > 1 and sys.argv[1] == 'calib'
    plist = [3001, 3011] if calib else primes_in(3000, 4200)
    nrand = 2 if calib else 3
    rng = np.random.default_rng(20260801)
    rows = []
    t0 = time.time()
    for i, p in enumerate(plist):
        key = orbit_keys(p)
        T3, T4, st = stats(key, p)
        bT3 = []; bT4 = []; bst = defaultdict(list)
        for _ in range(nrand):
            kb = mirror_random_key(p, rng)
            t3b, t4b, stb = stats(kb, p)
            bT3.append(t3b); bT4.append(t4b)
            for k2, v2 in stb.items(): bst[k2].append(v2)
        rows.append((p, T3, T4, st, np.mean(bT3), np.mean(bT4), {k2: np.mean(v2) for k2, v2 in bst.items()}, leg51(p)))
        if (i+1) % 20 == 0 or calib:
            print(f"[{i+1}/{len(plist)}] p={p} T3/p={T3/p:.3f} (base {np.mean(bT3)/p:.3f}) T4/p={T4/p:.3f} (base {np.mean(bT4)/p:.3f})", flush=True)
    # aggregate
    P = np.array([r[0] for r in rows], dtype=float)
    T3s = np.array([r[1] for r in rows]); T4s = np.array([r[2] for r in rows])
    bT3s = np.array([r[4] for r in rows]); bT4s = np.array([r[5] for r in rows])
    print("\n=== T3/T4 model check ===", flush=True)
    print(f"T3/p: data {np.mean(T3s/P):.4f} base {np.mean(bT3s/P):.4f} (Poisson预言 4.0 = (p+1)E[(2X)_3]/p)")
    print(f"T4/p: data {np.mean(T4s/P):.4f} base {np.mean(bT4s/P):.4f} (预言 10.0 = (p+1)E[(2X)_4]/p; bin数=p+1 非 N/2, Q6532 的 5 与 Q6534 的 24 均判错)")
    print("\n=== T3 stratification (mean fraction of triples, data vs base) ===")
    for k2 in ['mir3','pure3','ee','eo','oe','oo','smallgap']:
        dfrac = np.mean([r[3][k2]/max(1,r[3]['total3']) for r in rows])
        bfrac = np.mean([r[6][k2]/max(1,r[6]['total3']) for r in rows])
        print(f"{k2}: data {dfrac:.4f} base {bfrac:.4f} ratio {dfrac/max(bfrac,1e-12):.3f}")
    # split by (-51|p)
    for tag, sel in [("(-51|p)=+1", [r for r in rows if r[7]]), ("(-51|p)=-1", [r for r in rows if not r[7]])]:
        if not sel: continue
        Ps = np.array([r[0] for r in sel], dtype=float)
        print(f"{tag}: n={len(sel)} T3/p={np.mean([r[1] for r in sel]/Ps):.4f} smallgap_frac={np.mean([r[3]['smallgap']/max(1,r[3]['total3']) for r in sel]):.4f}")
    print(f"done {time.time()-t0:.1f}s", flush=True)

if __name__ == '__main__':
    main()
