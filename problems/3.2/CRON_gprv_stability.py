#!/usr/bin/env python3
"""CRON_gprv_stability.py — Q6516 dispatch-priority #1: GPRV 统计量 G/p 跨尺度稳定性.

G(p) = Σ_{h=1}^{p-3} (R_h − μ_h)²,  R_h = #{r: π(r)=π(r+h)},
μ_h = 奇偶类模型期望: 偶 h → 1 + 2λ_h, 奇 h → 2λ_h, λ_h = 真碰撞对(镜像对成双)均率.
这里用最诚实的口径: μ_h 取奇偶类经验均值(数据自定, 不引入模型参数),
G/p 稳定 ⟺ GPRV θ=0 经验成立. 波段: 3k / 10k / 30k / 100k 各 8 素数.
单写者: cron.
"""
from collections import Counter, defaultdict
from statistics import mean
from sympy import primerange, nextprime

def orbit(p):
    pts = [None]*(p-1)
    b0, b1, c0, c1 = 1 % p, 5 % p, 0, 1
    def key(x, y):
        return (1, y * pow(x, -1, p) % p) if x != 0 else (0, 1)
    pts[0] = key(b0, c0); pts[1] = key(b1, c1)
    bm2, bm1, cm2, cm1 = b0, b1, c0, c1
    for n in range(2, p - 1):
        A = (34*n*n*n - 51*n*n + 27*n - 5) % p
        B = ((n-1)**3) % p
        inv = pow((n*n*n) % p, -1, p)
        bn = (A*bm1 - B*bm2) * inv % p
        cn = (A*cm1 - B*cm2) * inv % p
        pts[n] = (1, cn * pow(bn, -1, p) % p) if bn != 0 else (0, 1)
        bm2, bm1, cm2, cm1 = bm1, bn, cm1, cn
    return pts

def G_over_p(p):
    pos = defaultdict(list)
    for n, v in enumerate(orbit(p)):
        pos[v].append(n)
    Rh = Counter()
    for l in pos.values():
        for i in range(len(l)):
            for j in range(i+1, len(l)):
                Rh[l[j]-l[i]] += 1
    ev = [Rh.get(h, 0) for h in range(2, p-2, 2)]
    od = [Rh.get(h, 0) for h in range(1, p-2, 2)]
    me, mo = mean(ev), mean(od)
    G = sum((x-me)**2 for x in ev) + sum((x-mo)**2 for x in od)
    return G / p, me, mo

for base in (3000, 10000, 30000, 100000):
    vals = []
    p = base
    for _ in range(8):
        p = int(nextprime(p))
        g, me, mo = G_over_p(p)
        vals.append(g)
        print(f"  p={p}: G/p={g:.4f} (μ_even={me:.3f} μ_odd={mo:.3f})", flush=True)
    print(f"== band {base}: mean G/p = {mean(vals):.4f} ==", flush=True)
print("DONE_GPRV", flush=True)
