#!/usr/bin/env python3
"""CRON_radon_spectrum.py — Q6527 §8.3 二变量 Artin-Schreier/Radon 谱探针 (优化 O(p^2 log p)).

D(r,s) = b_r c_s - b_s c_r mod p (索引按 Q6527 口径 mod p 回绕, b,c 定义在 0..p-1).
A(t,h) = sum_r e_p(t D(r, r+h)),  F(t,xi) = sum_h A(t,h) e_p(xi h).
统计: M_p = max_{t!=0, xi} |F(t,xi)|/p 随 p 的增长曲线.
判读: M_p ~ sqrt(2 log p) 级(随机水平) vs p^eta 增长(单值性路线判死) vs O(1) 有界(有界导子有戏).
另输出 max_t |A(t,h)|/sqrt(p) @ h=1,2,3 (Weil 对照: 小 h 应为 O(1))."""
import time
import numpy as np

def apery_pair(p):
    b = np.zeros(p, dtype=np.int64); c = np.zeros(p, dtype=np.int64)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    for n in range(1, p-1):
        Pn = (34*n**3 + 51*n*n + 27*n + 5) % p
        n3 = n**3 % p
        inv = pow((n+1)**3 % p, p-2, p)
        b[n+1] = (Pn*b[n] - n3*b[n-1]) % p * inv % p
        c[n+1] = (Pn*c[n] - n3*c[n-1]) % p * inv % p
    return b, c

def radon(p):
    b, c = apery_pair(p)
    # hist[h, d] = #{r: D(r, (r+h)%p) = d}
    hist = np.zeros((p, p), dtype=np.float64)
    for h in range(p):
        s = np.roll(np.arange(p), -h)  # s = (r+h) % p
        D = (b * c[s] - b[s] * c) % p
        np.add.at(hist[h], D, 1.0)
    # A[h, t] = sum_d hist[h,d] e_p(t d)  -> FFT along d
    A = np.fft.fft(hist, axis=1)  # A[h, t], t=0..p-1 (e^{-2pi i td/p} 约定, 模长统计不受影响)
    # F[xi, t] = sum_h A[h,t] e_p(xi h) -> FFT along h
    F = np.fft.fft(A, axis=0)     # F[xi, t]
    Fa = np.abs(F) / p
    Mp = Fa[:, 1:].max()          # t != 0, 全部 xi
    # Weil 对照
    weil = [np.abs(A[h, 1:]).max() / np.sqrt(p) for h in (1, 2, 3)]
    # 剔除 xi=0 的表面项单独看
    surf = Fa[0, 1:].max()
    return Mp, surf, weil

if __name__ == '__main__':
    for p in [101, 211, 401, 809, 1601, 3001]:
        t0 = time.time()
        Mp, surf, weil = radon(p)
        import math
        print(f"p={p}: M_p={Mp:.3f}  (sqrt(2 ln p)={math.sqrt(2*math.log(p)):.3f})  surface(xi=0)={surf:.3f}  maxA/sqrt(p)@h=1,2,3={weil[0]:.2f},{weil[1]:.2f},{weil[2]:.2f}  [{time.time()-t0:.1f}s]", flush=True)
