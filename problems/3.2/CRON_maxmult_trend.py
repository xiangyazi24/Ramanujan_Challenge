#!/usr/bin/env python3
"""CRON_maxmult_trend.py — 复核 life §95/§96 (2026-08-01):

 (i)  自指恒等式 (§96): 轨道重数 mult(v) = 该解的零点数.
      对 v=[α:β], 解 y_n = β·b_n − α·c_n 有 (y_0,y_1)=(β, 5β−α);
      mult(v) := #{0<=n<=p−2: π(n)=v} 应等于 #{n: y_n≡0}.  (线性代数上恒真; 机器逐点验证 p=101, 211.)
      特例: v=[0:1]=L ⟹ mult(L)=|Z_p| (Apéry 零点数).  max_v mult(v) = 全解族 max|Z|.
 (ii) §95 双口径: E^off = Σ_v m(m−1) ≈ 2p,  Σ_v m² = E^π ≈ 3p,  差 = p+1 个轨道点 (对角).
 (iii) max mult 的 p-趋势: life 测 p≈10⁴/10⁵/10⁶ 单素数 max=12 恒定, §96 修正预期为
      极值律 ~log p/loglog p (增长, 非有界).  这里按波段采样多素数, 看 per-p max 的
      带内分布与跨带漂移.  (Poisson(λ) 极值 M_N ~ 逆 Gamma 尾: 对 N≈p 个解,
      P[某解零点数≥k] ≈ p·λ^k/k! → max ≈ 解 λ^k/k! = 1/p, 即 k ~ log p/loglog p.)

输出: 每素数一行 (进度), 每波段小结.  单写者: cron.
"""
import sys
from sympy import primerange

def orbit_stats(p):
    """返回 (maxmult, E_pi, E_off, Zp, spec_even_share, mults_dict_or_None)."""
    # 轨道 π(n)=[b_n:c_n], n=0..p-2;  n³y_n=(34n³−51n²+27n−5)y_{n−1}−(n−1)³y_{n−2}
    mult = {}
    b0, b1 = 1 % p, 5 % p
    c0, c1 = 0, 1
    def key(x, y):
        if x != 0:
            return (1, y * pow(x, -1, p) % p)   # [1 : y/x]
        return (0, 1)                            # [0:1] = L
    mult[key(b0, c0)] = 1
    k1 = key(b1, c1)
    mult[k1] = mult.get(k1, 0) + 1
    Zp = (1 if b0 == 0 else 0) + (1 if b1 == 0 else 0)
    bm2, bm1, cm2, cm1 = b0, b1, c0, c1
    for n in range(2, p - 1):
        A = (34 * n * n * n - 51 * n * n + 27 * n - 5) % p
        B = ((n - 1) ** 3) % p
        inv = pow((n * n * n) % p, -1, p)
        bn = (A * bm1 - B * bm2) * inv % p
        cn = (A * cm1 - B * cm2) * inv % p
        k = key(bn, cn)
        mult[k] = mult.get(k, 0) + 1
        if bn == 0:
            Zp += 1
        bm2, bm1, cm2, cm1 = bm1, bn, cm1, cn
    E_pi = sum(m * m for m in mult.values())
    E_off = sum(m * (m - 1) for m in mult.values())
    even_pts = sum(1 for m in mult.values() if m % 2 == 0)
    return max(mult.values()), E_pi, E_off, Zp, even_pts / len(mult), mult

def zeros_of_solution(p, y0, y1):
    z = (1 if y0 % p == 0 else 0) + (1 if y1 % p == 0 else 0)
    ym2, ym1 = y0 % p, y1 % p
    for n in range(2, p - 1):
        A = (34 * n * n * n - 51 * n * n + 27 * n - 5) % p
        B = ((n - 1) ** 3) % p
        yn = (A * ym1 - B * ym2) * pow((n * n * n) % p, -1, p) % p
        if yn == 0:
            z += 1
        ym2, ym1 = ym1, yn
    return z

fails = 0
def chk(c, m):
    global fails
    print(("  [OK ] " if c else "  [FAIL] ") + m, flush=True)
    if not c:
        fails += 1

# ---------- (i) 自指恒等式逐点验证 ----------
print("== (i) §96 自指恒等式: mult(v) = 零点数(解 y0=β, y1=5β−α), p=101/211 全 v ==", flush=True)
for p in (101, 211):
    _, _, _, Zp, _, mult = orbit_stats(p)
    ok = True
    # 全部 P^1 点: [1:t] ↦ α=1? 注意 key 约定 [x:y]: x=1 ⟹ v=[1:t] 即 α=1, β=t... 统一用 (α,β)=v=[α:β]
    # key(1,t) 代表 [1:t] ⟹ α=1, β=t;  key(0,1)=[0:1] ⟹ α=0, β=1
    for t in range(p):
        v_mult = mult.get((1, t), 0)
        z = zeros_of_solution(p, t, (5 * t - 1) % p)     # α=1, β=t
        if v_mult != z:
            ok = False
            print(f"    MISMATCH p={p} v=[1:{t}]: mult={v_mult} zeros={z}", flush=True)
    zL = zeros_of_solution(p, 1, 5)                       # v=[0:1]: α=0, β=1
    if mult.get((0, 1), 0) != zL or zL != Zp:
        ok = False
    chk(ok, f"p={p}: mult(v)=zeros(v) 全 {p+1} 点; mult(L)={mult.get((0,1),0)}=|Z_p|={Zp}")

# ---------- (ii)+(iii) 波段采样 ----------
BANDS = [(1000, 1400), (10000, 10600), (100000, 100700), (1000000, 1002200)]
print("\n== (ii)+(iii) 波段采样: per-p max mult / E^π/p / E^off/p / 偶份额 ==", flush=True)
summary = []
for lo, hi in BANDS:
    ps = list(primerange(lo, hi))
    maxes, epis, eoffs = [], [], []
    for i, p in enumerate(ps):
        mm, E_pi, E_off, Zp, evsh, _ = orbit_stats(p)
        maxes.append(mm)
        epis.append(E_pi / p)
        eoffs.append(E_off / p)
        print(f"  band[{lo}] {i+1}/{len(ps)} p={p}: max={mm} E^π/p={E_pi/p:.3f} "
              f"E^off/p={E_off/p:.3f} |Z_p|={Zp} even={evsh:.4f}", flush=True)
    from statistics import mean
    summary.append((lo, len(ps), max(maxes), mean(maxes), mean(epis), mean(eoffs)))
    print(f"  == band [{lo},{hi}) n={len(ps)}: MAX(max)={max(maxes)} mean(max)={mean(maxes):.2f} "
          f"E^π/p={mean(epis):.4f} E^off/p={mean(eoffs):.4f} ==", flush=True)

print("\n== SUMMARY (band_lo, #p, MAX, mean_max, E^π/p, E^off/p) ==", flush=True)
import math
for lo, n, mx, mn, ep, eo in summary:
    pred = math.log(lo) / math.log(math.log(lo))
    print(f"  {lo:>8} n={n:<3} MAX={mx:<3} mean_max={mn:6.2f} E^π/p={ep:.4f} E^off/p={eo:.4f} "
          f"[log p/loglog p = {pred:.2f}]", flush=True)
print(f"\nFAILS={fails}", flush=True)
sys.exit(1 if fails else 0)
