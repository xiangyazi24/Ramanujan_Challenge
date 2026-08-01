#!/usr/bin/env python3
"""CRON_bt_fullscan.py — [GAP-2D-SQRT] 全频率扫描 (Q6526 必补实验#1).

B_t(H) = sum_{h<=H} sum_{0<=r<=p-2-h} e_p(t * D_h(r)), D_h(r)=b_r c_{r+h}-b_{r+h} c_r (真轨道口径).
对全部 t 同时算: 直方图 cnt_H[d] 增量累积 + 一次 FFT.
输出: 每 (p, H): max_{t!=0} |B_t|/sqrt(pH), 中位数, 与 t=1,2,5 (AQ.3 旧口径) 对照.
判据: 猜想要求 max_t 一致 (pH)^{1/2} p^eps —— 若 max/sqrt(pH) 随 H 或 p 幂增长, [GAP-2D-SQRT] 假.
H 扫到 p^0.6 (超临界线之上, Q6526 要求)."""
import time
import numpy as np

def apery_pair(p):
    b = np.zeros(p - 1, dtype=np.int64); c = np.zeros(p - 1, dtype=np.int64)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    for n in range(1, p - 2):
        Pn = (34*n**3 + 51*n*n + 27*n + 5) % p
        n3 = n**3 % p
        inv = pow((n+1)**3 % p, p - 2, p)
        b[n+1] = (Pn*b[n] - n3*b[n-1]) % p * inv % p
        c[n+1] = (Pn*c[n] - n3*c[n-1]) % p * inv % p
    return b, c

def scan(p):
    b, c = apery_pair(p)
    N = p - 1
    Hmax = int(p ** 0.6)
    Hs = []
    H = 16
    while H <= Hmax: Hs.append(H); H *= 2
    if Hs[-1] != Hmax: Hs.append(Hmax)
    cnt = np.zeros(p, dtype=np.float64)
    results = []
    h = 1
    for H in Hs:
        while h <= H:
            D = (b[:N-h] * c[h:] - b[h:] * c[:N-h]) % p
            np.add.at(cnt, D, 1.0)
            h += 1
        B = np.fft.fft(cnt)
        Ba = np.abs(B[1:]) / np.sqrt(p * H)
        r135 = [float(Ba[t-1]) for t in (1, 2, 5)]
        results.append((H, float(Ba.max()), float(np.median(Ba)), r135))
    return results

if __name__ == '__main__':
    for p in [10007, 30011, 100003]:
        t0 = time.time()
        res = scan(p)
        print(f"p={p} (p^0.5={int(p**0.5)}):", flush=True)
        for H, mx, med, r135 in res:
            tag = 'SUPER' if H > p ** 0.5 else '     '
            print(f"  H={H:6d} {tag} max|B|/sqrt(pH)={mx:7.3f} median={med:.3f} t=1,2,5: {r135[0]:.2f},{r135[1]:.2f},{r135[2]:.2f}", flush=True)
        print(f"  [{time.time()-t0:.0f}s]", flush=True)
