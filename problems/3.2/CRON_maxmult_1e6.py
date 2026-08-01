#!/usr/bin/env python3
"""CRON_maxmult_1e6.py — 10⁶ 波段 max-mult 极值律采样（优化版, 接替 CRON_maxmult_trend.py 的
band[1000000]; 前三波段结果以原 log 为准）.

优化: 整数化身递推 B_n=A(n)B_{n−1}−(n−1)⁶B_{n−2}, C 同 (C_0=0, C_1=1·1³=1) — 无逆元;
π(n)=[B_n:C_n]. 归一化 key 用 Montgomery 批量求逆 (一次 pow + 3(N−1) 乘).
验证锚: 与原脚本在 p=1000003 交叉核对 (max/E^π/E^off/|Z_p| 必须逐项相等).
"""
import sys, math
from statistics import mean
from sympy import primerange

def orbit_stats_fast(p):
    N = p - 1
    B = [0] * N
    C = [0] * N
    B[0], B[1] = 1 % p, 5 % p          # (0!)³·1, (1!)³·5
    C[0], C[1] = 0, 1
    for n in range(2, N):
        A = (34*n*n*n - 51*n*n + 27*n - 5) % p
        D = ((n-1)**3 % p)
        D = D * D % p                   # (n-1)^6
        B[n] = (A * B[n-1] - D * B[n-2]) % p
        C[n] = (A * C[n-1] - D * C[n-2]) % p
    # 批量求逆全部非零 B[n]
    idx = [n for n in range(N) if B[n]]
    pref = [1] * (len(idx) + 1)
    for i, n in enumerate(idx):
        pref[i+1] = pref[i] * B[n] % p
    inv_all = pow(pref[-1], p - 2, p)
    invs = [0] * len(idx)
    for i in range(len(idx) - 1, -1, -1):
        invs[i] = pref[i] * inv_all % p
        inv_all = inv_all * B[idx[i]] % p
    mult = {}
    L = p  # sentinel key for [0:1]
    Zp = 0
    j = 0
    for n in range(N):
        if B[n]:
            k = C[n] * invs[j] % p
            j += 1
        else:
            k = L
            Zp += 1
        mult[k] = mult.get(k, 0) + 1
    E_pi = sum(m*m for m in mult.values())
    E_off = E_pi - N
    return max(mult.values()), E_pi, E_off, Zp

# 交叉核对锚: p=1000003 原脚本已出 max=12 E^π/p=3.002 E^off/p=2.002 |Z_p|=2
mm, E_pi, E_off, Zp = orbit_stats_fast(1000003)
print(f"ANCHOR p=1000003: max={mm} E^π/p={E_pi/1000003:.3f} E^off/p={E_off/1000003:.3f} |Z_p|={Zp} "
      f"(expect 12 / 3.002 / 2.002 / 2)", flush=True)
assert mm == 12 and Zp == 2, "ANCHOR MISMATCH — abort"

ps = list(primerange(1000000, 1002200))
maxes, epis = [], []
for i, p in enumerate(ps):
    mm, E_pi, E_off, Zp = orbit_stats_fast(p)
    maxes.append(mm); epis.append(E_pi / p)
    print(f"  1e6 {i+1}/{len(ps)} p={p}: max={mm} E^π/p={E_pi/p:.3f} |Z_p|={Zp}", flush=True)
print(f"\n== BAND [1000000,1002200) n={len(ps)}: MAX(max)={max(maxes)} mean_max={mean(maxes):.2f} "
      f"E^π/p={mean(epis):.4f} [log p/loglog p = {math.log(1e6)/math.log(math.log(1e6)):.2f}] ==", flush=True)
print("DONE_1E6", flush=True)
