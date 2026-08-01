#!/usr/bin/env python3
"""CRON_rh_variance.py — attack #3 前置: (p,h) 族 R_h 分布/方差机器实验.

R_h = #{r: π(r)=π(r+h)} (0<=r<r+h<=p-2). 轨道语言的结构预言:
 - 镜像逐点恒等 π(p−1−n)=π(n) ⟹ 每个 n∈[1,(p−3)/2] 给一个强制碰撞, gap h=p−1−2n **恒偶**
   (p 奇), 且每个偶 h∈[2,p−3] 恰好一次 ⟹ R_h^{even} = 1 + X_h, R_h^{odd} = X_h,
   X_h ≈ 真碰撞贡献, 均值 ~ (Σ真碰撞)/(#gaps) ≈ (p/8)/(p−2) ≈ 1/8.
 - 攻击 #3 的问题 = X_h 是否跨 h Poisson (方差=均值, 无算术超额).
统计: 逐素数把 {R_h} 按奇偶分类, 聚合分布 histogram + 均值/方差; 镜像随机基准同算.
单写者: cron.
"""
import random
from collections import Counter, defaultdict
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
    pts = [None] * (p - 1)
    for n in range((p - 1) // 2):
        v = (1, rng.randrange(p)) if rng.randrange(p + 1) else (0, 1)
        pts[n] = v
        pts[p - 2 - n] = v
    for i in range(p - 1):
        if pts[i] is None:
            pts[i] = (1, rng.randrange(p))
    return pts

def rh_profile(pts, p, src):
    pos = defaultdict(list)
    for n, v in enumerate(pts):
        pos[v].append(n)
    Rh = Counter()
    first_gap = Counter()   # P 部分 (life [GAP-QRLL]): 每个 r 的首返回 gap 分布
    for l in pos.values():
        for i in range(len(l)):
            for j in range(i + 1, len(l)):
                Rh[l[j] - l[i]] += 1
        for i in range(len(l) - 1):
            first_gap[l[i+1] - l[i]] += 1   # 相邻出现 = 首返回
    ev, od = [], []
    for h in range(1, p - 2):
        (ev if h % 2 == 0 else od).append(Rh.get(h, 0))
    return ev, od, first_gap

BAND = (3000, 4200)
ps = list(primerange(*BAND))
rng = random.Random(31415)
acc = {('data','even'): Counter(), ('data','odd'): Counter(),
       ('rand','even'): Counter(), ('rand','odd'): Counter()}
mv = {k: [] for k in acc}
print(f"primes: {len(ps)} in {BAND}", flush=True)
fg_acc = {'data': Counter(), 'rand': Counter()}
fg_H = {'data': [], 'rand': []}   # 归一化首返回累积: P_p(H)/H @ H=√p
for i, p in enumerate(ps):
    for src, pts in (('data', orbit(p)), ('rand', random_orbit(p, rng))):
        ev, od, fg = rh_profile(pts, p, src)
        for par, lst in (('even', ev), ('odd', od)):
            acc[(src, par)].update(lst)
            m = mean(lst); v = mean(x*x for x in lst) - m*m
            mv[(src, par)].append((m, v))
        H = int(p ** 0.5)
        fg_acc[src].update({min(g, 10**9): c for g, c in fg.items()})
        fg_H[src].append(sum(c for g, c in fg.items() if g <= H) / H)
    if i % 20 == 0:
        print(f"  {i+1}/{len(ps)} p={p}", flush=True)

print("\n== R_h 分布 (聚合计数: R_h值 -> #gaps) ==", flush=True)
for k in acc:
    tot = sum(acc[k].values())
    hist = {r: acc[k][r] for r in sorted(acc[k])}
    print(f"  {k}: {hist}  (N={tot})", flush=True)
print("\n== 均值/方差 (跨素数平均) ==", flush=True)
for k in mv:
    ms = mean(x[0] for x in mv[k]); vs = mean(x[1] for x in mv[k])
    print(f"  {k}: mean(R_h)={ms:.4f}  var(R_h)={vs:.4f}  (Poisson 检验: 偶类 var 应≈mean−1? "
          f"精确: R=1+X ⟹ var(R)=var(X); 奇类 var≈mean)", flush=True)
print("\n== P 部分: 首返回 gap 分布 (life [GAP-QRLL] P_p(H)) ==", flush=True)
for src in ('data', 'rand'):
    small = {g: fg_acc[src][g] for g in sorted(fg_acc[src]) if g <= 20}
    tot = sum(fg_acc[src].values())
    print(f"  {src}: 小 gap(≤20) {small} | 总首返回数={tot} | "
          f"mean[P_p(√p)/√p]={mean(fg_H[src]):.4f}", flush=True)
print("DONE_RH", flush=True)
