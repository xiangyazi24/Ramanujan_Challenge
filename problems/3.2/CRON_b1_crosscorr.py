#!/usr/bin/env python3
"""CRON_b1_crosscorr.py — avenue (b) B1: 跨 (h1,h2) 根集相关性经验测量.

对象: root set R_h = {r: p | N_h(r)} = {r: π(r)=π(r+h)} (轨道碰撞口径, 已证等价).
问题: 方差路线的非对角项 ⟺ 不同 h 的根集是否(近)独立.
注意: 相关性有一个**被迫下限** — 传递性: π(r)=π(r+h1) ∧ π(r)=π(r+h2) ⟹ r+h1 是
N_{h2−h1} 的根. 所以 J(h1,h2)=|R_{h1}∩R_{h2}| 数的是三重碰撞, 随机基准必须用
"同样的轨道模型" (镜像 2-对-1 + 均匀半轨道) 而不是独立伯努利.

测量:
 (1) 每素数: 对 h ∈ H 计算 R_h, 全部 pair 的 J(h1,h2); 归一 ρ = J·L/(|R_h1|·|R_h2|),
     L = 有效区间长. 独立 ⟹ ρ≈1.
 (2) 三阶矩 T3 = Σ_v m(m−1)(m−2) (= 有序三重碰撞总数 = Σ_{h1≠h2} J 的对角化身).
 (3) 随机基准: 均匀随机半轨道 (N=(p−1)/2 点均匀入 p+1 槽, 镜像翻倍), 同样统计.
输出逐素数进度 + 聚合表. 单写者: cron.
"""
import random, sys
from statistics import mean
from sympy import primerange

def orbit(p):
    pts = []
    b0, b1, c0, c1 = 1 % p, 5 % p, 0, 1
    def key(x, y):
        return (1, y * pow(x, -1, p) % p) if x != 0 else (0, 1)
    pts.append(key(b0, c0)); pts.append(key(b1, c1))
    bm2, bm1, cm2, cm1 = b0, b1, c0, c1
    for n in range(2, p - 1):
        A = (34*n*n*n - 51*n*n + 27*n - 5) % p
        B = ((n-1)**3) % p
        inv = pow((n*n*n) % p, -1, p)
        bn = (A*bm1 - B*bm2) * inv % p
        cn = (A*cm1 - B*cm2) * inv % p
        pts.append((1, cn * pow(bn, -1, p) % p) if bn != 0 else (0, 1))
        bm2, bm1, cm2, cm1 = bm1, bn, cm1, cn
    return pts

def random_orbit(p, rng):
    half = [(1, rng.randrange(p)) if rng.randrange(p+1) else (0, 1) for _ in range((p-1)//2)]
    # 镜像: pts[p-1-n] = pts[n] (γ=id 逐点恒等, AK.2)
    pts = [None]*(p-1)
    for n, v in enumerate(half):
        pts[n] = v
        pts[p-2-n] = v          # n ↔ p-2-n? 镜像是 n ↔ p-1-n, 索引 0..p-2: n ↔ p-1-n 越界1位
    # 修正: 真镜像 c_{p-1-n}=c_n 作用在 n=0..p-1, 轨道只到 p-2; 用 n ↔ p-2-n 近似(差一格不影响统计)
    if (p-1) % 2 == 1:
        pts[(p-1)//2] = half[-1] if half else (0, 1)
    for i in range(p-1):
        if pts[i] is None:
            pts[i] = (1, rng.randrange(p))
    return pts

def stats_for(pts, p, H):
    from collections import defaultdict
    pos = defaultdict(list)
    for n, v in enumerate(pts):
        pos[v].append(n)
    T3 = sum(len(l)*(len(l)-1)*(len(l)-2) for l in pos.values())
    roots = {}
    for h in H:
        if h >= p - 2: continue
        roots[h] = set(r for r in range(0, p-1-h) if pts[r] == pts[r+h])
    rows = []
    hs = sorted(roots)
    for i in range(len(hs)):
        for j in range(i+1, len(hs)):
            h1, h2 = hs[i], hs[j]
            L = p - 1 - h2
            R1 = set(r for r in roots[h1] if r < L)
            R2 = set(r for r in roots[h2] if r < L)
            if not R1 or not R2: continue
            J = len(R1 & R2)
            rows.append((h1, h2, J, len(R1)*len(R2)/L))
    return T3, rows

BANDS = [(3000, 4200)]
H = [1, 2, 3, 5, 8, 13, 21, 55]
agg_J, agg_E, agg_T3, agg_T3r, agg_Jr, agg_Er = 0, 0.0, [], [], 0, 0.0
pair_acc = {}
rng = random.Random(20260801)
ps = [q for lo, hi in BANDS for q in primerange(lo, hi)]
print(f"primes: {len(ps)} in {BANDS}, H={H}", flush=True)
for i, p in enumerate(ps):
    pts = orbit(p)
    T3, rows = stats_for(pts, p, H + [p//3, (p-1)//2])
    agg_T3.append(T3/p)
    for h1, h2, J, E in rows:
        agg_J += J; agg_E += E
        k = (h1, h2) if h2 in H else ('large', 'large')
        a = pair_acc.setdefault(k, [0, 0.0]); a[0] += J; a[1] += E
    rpts = random_orbit(p, rng)
    T3r, rrows = stats_for(rpts, p, H + [p//3, (p-1)//2])
    agg_T3r.append(T3r/p)
    for h1, h2, J, E in rrows:
        agg_Jr += J; agg_Er += E
    if i % 10 == 0:
        print(f"  {i+1}/{len(ps)} p={p}: T3/p={T3/p:.3f} (rand {T3r/p:.3f}) "
              f"cumJ={agg_J} cumE={agg_E:.1f} (rand J={agg_Jr} E={agg_Er:.1f})", flush=True)

print("\n== AGGREGATE ==", flush=True)
print(f"  data:  ΣJ={agg_J}  ΣE[J]_indep={agg_E:.1f}  ratio={agg_J/agg_E if agg_E else float('nan'):.3f}", flush=True)
print(f"  rand:  ΣJ={agg_Jr}  ΣE[J]_indep={agg_Er:.1f}  ratio={agg_Jr/agg_Er if agg_Er else float('nan'):.3f}", flush=True)
print(f"  T3/p:  data mean={mean(agg_T3):.4f}   rand mean={mean(agg_T3r):.4f}", flush=True)
print("\n== per-pair (small h) ratio J/E ==", flush=True)
for k in sorted(pair_acc, key=str):
    J, E = pair_acc[k]
    if E > 5:
        print(f"  {k}: J={J} E={E:.1f} ratio={J/E:.3f}", flush=True)
