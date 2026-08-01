#!/usr/bin/env python3
"""CRON_weyl2d.py — attack #2 前置: 二维联合 Weyl 和 T(α,β)=Σ_{n=0}^{p-2} e_p(α b_n + β c_n) 全谱.

E^π 与 T 的关系: 碰撞检测在 (b,c) 平面的 Fourier 侧 ⟹ 平方根相消 max|T|≪p^{1/2+ε}
将严格强于单坐标 [GAP-RES-WEYL]. 本实验: 全 (α,β) 谱 = 点集 {(b_n,c_n)} 计数测度的
2D DFT (numpy FFT2, O(p²log p)). 统计 max_{(α,β)≠(0,0)}|T|/√p 与 |T|² 分布尾部,
对照: 随机点集 (p−1 个均匀 iid 点) + 理论 Parseval Σ|T|²=p²·(p−1)/p... (直接机器对账).
素数样本: 211, 499, 1009, 2003 (p² 网格内存 ~ p²·16B, 2003²≈64MB ok).
单写者: cron.
"""
import numpy as np
from sympy import isprime

def orbit_bc(p):
    b = np.zeros(p - 1, dtype=np.int64)
    c = np.zeros(p - 1, dtype=np.int64)
    b[0], b[1] = 1, 5
    c[0], c[1] = 0, 1
    for n in range(2, p - 1):
        A = (34 * n * n * n - 51 * n * n + 27 * n - 5) % p
        D = pow(n, 3, p)
        Dinv = pow(D, p - 2, p)
        B = pow(n - 1, 3, p)
        b[n] = (A * b[n - 1] - B * b[n - 2]) * Dinv % p
        c[n] = (A * c[n - 1] - B * c[n - 2]) * Dinv % p
    return b % p, c % p

def spectrum(px, py, p):
    grid = np.zeros((p, p))
    np.add.at(grid, (px, py), 1.0)
    T = np.fft.fft2(grid)
    absT = np.abs(T)
    absT[0, 0] = 0.0          # 去掉主项
    return absT

rng = np.random.default_rng(20260801)
print("p | max|T|/√p (orbit) | max|T|/√p (rand) | √(2 log p²)≈GUE基准 | Σ|T|²/p² 对账(=N−N²/p²·p²...直接Parseval)", flush=True)
for p in (211, 499, 1009, 2003):
    assert isprime(p)
    b, c = orbit_bc(p)
    A = spectrum(b.astype(int), c.astype(int), p)
    mx = A.max() / np.sqrt(p)
    # Parseval: Σ_{(α,β)}|T|² = p²·Σ_x m(x)² (m=点重数);  去主项后 = p²·(Σm²) − N²
    N = p - 1
    par = (A ** 2).sum()
    m2 = (par + N * N) / (p * p)   # 应 = Σ_x m(x)² (二维点重数, 轨道单射时 = N)
    rx = rng.integers(0, p, N); ry = rng.integers(0, p, N)
    Ar = spectrum(rx, ry, p)
    mxr = Ar.max() / np.sqrt(p)
    bench = np.sqrt(2 * np.log(p * p))
    print(f"  p={p}: orbit {mx:.3f} | rand {mxr:.3f} | extreme-bench {bench:.3f} | Σm²(2D)={m2:.1f} (N={N})", flush=True)
print("DONE_W2D", flush=True)
