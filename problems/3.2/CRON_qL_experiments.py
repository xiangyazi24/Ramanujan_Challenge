#!/usr/bin/env python3
"""CRON_qL_experiments.py — Q6511 §9.3 实验批 (1)(3)(4)(6):

 (1) 口径锁: 3·Σ_{h≠k}J_nonwrap = T3, p=101/211 全滞后字面枚举.
 (3) dyadic 尺度曲线: Q_p(H)=Σ_r C(d_H(r),2), H=2,4,...,~p/2; 报 Q_p(H)/(H+H²/p)
     — (PAIR-QRLL) 的直接经验形状 (family 定理测试, 固定十个 h 测不了).
 (4) A_p(H)=Σ_{h≤H}R_h²; 报 A_p(H)/H 平坦性 — (SAME-LAG-BDH) 证据.
 (6) 径向 Fourier 恒等式: p²E^π=(p−1)N²+Σ_ℓ(pM_ℓ−N)² — p=5 玩具(答案给的向量)+真轨道核对.
单写者: cron.
"""
from collections import Counter, defaultdict
from sympy import isprime

def orbit_pos(p):
    pos = defaultdict(list)
    b0, b1, c0, c1 = 1 % p, 5 % p, 0, 1
    def key(x, y):
        return (1, y * pow(x, -1, p) % p) if x != 0 else (0, 1)
    pos[key(b0, c0)].append(0); pos[key(b1, c1)].append(1)
    bm2, bm1, cm2, cm1 = b0, b1, c0, c1
    for n in range(2, p - 1):
        A = (34*n*n*n - 51*n*n + 27*n - 5) % p
        B = ((n-1)**3) % p
        inv = pow((n*n*n) % p, -1, p)
        bn = (A*bm1 - B*bm2) * inv % p
        cn = (A*cm1 - B*cm2) * inv % p
        pos[(1, cn * pow(bn, -1, p) % p) if bn != 0 else (0, 1)].append(n)
        bm2, bm1, cm2, cm1 = bm1, bn, cm1, cn
    return pos

fails = 0
def chk(c, m):
    global fails
    print(("  [OK ] " if c else "  [FAIL] ") + m, flush=True)
    if not c: fails += 1

# (1) 口径锁
print("== (1) 口径锁: 3·ΣJ_nonwrap(有序 h≠k) = T3 ==", flush=True)
for p in (101, 211):
    pos = orbit_pos(p)
    T3 = sum(len(l)*(len(l)-1)*(len(l)-2) for l in pos.values())
    # 字面枚举: 同 base r=l[i], 两不同滞后 (h,k) 有序对
    J_sum = 0
    for l in pos.values():
        m = len(l)
        for i in range(m):
            t = m - 1 - i           # base l[i] 之后的出现数
            J_sum += t * (t - 1)    # 有序 (h,k), h≠k
    chk(3 * J_sum == T3, f"p={p}: 3·{J_sum} == T3={T3}")

# (3)+(4) dyadic 曲线
print("\n== (3) Q_p(H)/(H+H²/p) 与 (4) A_p(H)/H, dyadic H ==", flush=True)
for p in (3001, 10007, 30011, 100003):
    assert isprime(p)
    pos = orbit_pos(p)
    base_gaps = defaultdict(list)
    Rh = Counter()
    for l in pos.values():
        for i in range(len(l)):
            for j in range(i+1, len(l)):
                g = l[j] - l[i]
                base_gaps[l[i]].append(g)
                Rh[g] += 1
    Hs, row_q, row_a = [], [], []
    H = 2
    while H <= p // 2:
        Q = 0
        for gaps in base_gaps.values():
            d = sum(1 for g in gaps if g <= H)
            Q += d * (d - 1) // 2
        A = sum(Rh[h]**2 for h in range(1, H + 1))
        Hs.append(H); row_q.append(Q / (H + H * H / p)); row_a.append(A / H)
        H *= 4
    print(f"  p={p}:", flush=True)
    print(f"    H      : {Hs}", flush=True)
    print(f"    Q-ratio: {[f'{x:.3f}' for x in row_q]}", flush=True)
    print(f"    A/H    : {[f'{x:.3f}' for x in row_a]}", flush=True)

# (6) 径向恒等式
print("\n== (6) 径向 Fourier 恒等式 p²E^π=(p−1)N²+Σ_ℓ(pM_ℓ−N)² ==", flush=True)
# p=5 玩具: 向量 (1,0),(2,0),(0,1) — 方向重数 2,1; E=5 (答案 §7.1)
p, N = 5, 3
mult = {(1, 0): 2, (0, 1): 1}
E = sum(m*m for m in mult.values())
tot = 0
for l in [(1, t) for t in range(p)] + [(0, 1)]:
    M = mult.get(l, 0)
    tot += (p*M - N)**2
chk(p*p*E == (p-1)*N*N + tot, f"p=5 玩具: {p*p*E} == {(p-1)*N*N}+{tot}")
for p in (101, 211):
    pos = orbit_pos(p)
    N = p - 1
    E = sum(len(l)**2 for l in pos.values())
    tot = 0
    for l in [(1, t) for t in range(p)] + [(0, 1)]:
        M = len(pos.get(l, []))
        tot += (p*M - N)**2
    chk(p*p*E == (p-1)*N*N + tot, f"p={p}: 径向恒等式 (E^π={E})")

print(f"\nFAILS={fails}", flush=True)
