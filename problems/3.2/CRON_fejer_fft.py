#!/usr/bin/env python3
"""CRON_fejer_fft.py — Q6515 §8.2 时间低频谱账 (第一优先实验).

对每个素数 p:
  轨道 pi(n)=[b_n:c_n], n=0..p-2 (N=p-1 点).
  类内有序对 mod-p 差分自相关 C_p(h);  A_p(k) = FFT_h(C_p) = sum_v |hat 1_v(k)|^2.
  口径锁: A_p(0)=E^pi;  sum_k A_p(k) = p*N;  Fejer 恒等式
    (1/p) sum_k (A_p(k)-N) F_H(k) = 2 sum_{h=1..H} (1-h/(H+1)) C_p(h)   [C_p 口径]
  输出: L_p(H) (data vs 镜像随机基准), 低频带 |k|<=p/H 残差 (max/L1/L2).

判据 (Q6515 §8.2): 若 P/H~1.444 的轻微排斥来自真正时间谱隙, 低频残差应低于镜像
随机基准, 而不仅是边缘分布略低.
"""
import sys, math, time
import numpy as np

def primes_in(a, b):
    sieve = np.ones(b + 1, dtype=bool); sieve[:2] = False
    for i in range(2, int(b ** 0.5) + 1):
        if sieve[i]: sieve[i * i::i] = False
    return [int(x) for x in np.nonzero(sieve)[0] if x >= a]

def orbit_keys(p):
    """返回长度 N=p-1 的数组 key[n]: pi(n) 的射影键 (0..p-1 = b*c^{-1}; p = 无穷远 c=0)."""
    N = p - 1
    b = np.zeros(N, dtype=np.int64); c = np.zeros(N, dtype=np.int64)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    for n in range(1, N - 1):
        Pn = (34 * n * n * n + 51 * n * n + 27 * n + 5) % p
        n3 = (n * n * n) % p
        inv = pow((n + 1) ** 3 % p, p - 2, p)
        b[n + 1] = ((Pn * b[n] - n3 * b[n - 1]) % p) * inv % p
        c[n + 1] = ((Pn * c[n] - n3 * c[n - 1]) % p) * inv % p
    key = np.empty(N, dtype=np.int64)
    for n in range(N):
        key[n] = (b[n] * pow(int(c[n]), p - 2, p)) % p if c[n] != 0 else p
    return key

def autocorr_Cp(key, p):
    """C_p(h), h=0..p-1: 类内有序对 (n,m) 的 n-m mod p 差分计数. 成本 O(E^pi)."""
    N = p - 1
    from collections import defaultdict
    cls = defaultdict(list)
    for n in range(N): cls[int(key[n])].append(n)
    C = np.zeros(p, dtype=np.float64)
    for v, lst in cls.items():
        arr = np.array(lst)
        if len(arr) == 1: C[0] += 1; continue
        d = (arr[:, None] - arr[None, :]) % p
        np.add.at(C, d.ravel(), 1.0)
    return C

def mirror_random_key(p, rng):
    """γ-对称随机模型: 镜像类代表 n 与 p-1-n (mod p-1) 同值, 值均匀取 P^1 (p+1 个)."""
    N = p - 1
    key = np.full(N, -1, dtype=np.int64)
    for n in range(N):
        if key[n] >= 0: continue
        m = (p - 1 - n) % N
        v = rng.integers(0, p + 1)
        key[n] = v; key[m] = v
    return key

def analyze(C, p, H_list):
    N = p - 1
    A = np.fft.fft(C).real  # A_p(k); C 对称 => 实
    Epi = C[0] if False else A.mean() * 0 + C.sum() * 0 + A[0]  # A[0]=sum C = E^pi
    out = {'Epi': float(A[0])}
    # 口径锁: sum_k A_p(k) = p * C_p(0) = p*N
    parseval = A.sum()
    out['parseval_ok'] = abs(parseval - p * N) < 1e-3 * p * N
    res = A - N  # A_p(k) - N
    ks = np.fft.fftfreq(p, d=1.0 / p)  # 频率标号 (整数, 含负)
    for H in H_list:
        # Fejer 核 F_H(k) = sum_{|h|<=H} (1-|h|/(H+1)) e_p(kh) 的闭式:
        # F_H(k) = (1/(H+1)) * (sin(pi (H+1) k /p) / sin(pi k/p))^2, k!=0; F_H(0)=H+1
        kk = np.arange(p)
        with np.errstate(divide='ignore', invalid='ignore'):
            num = np.sin(np.pi * (H + 1) * kk / p) ** 2
            den = np.sin(np.pi * kk / p) ** 2
            F = num / den / (H + 1)
        F[0] = H + 1
        L = (res * F).sum() / p
        # 直接口径 (RHS): 2 sum_{h=1..H} (1-h/(H+1)) C(h)
        hh = np.arange(1, H + 1)
        Ldir = 2.0 * ((1 - hh / (H + 1)) * C[1:H + 1]).sum() + (C[0] - N)  # C[0]=N => 0
        band = int(max(1, p // H))
        lowk = np.r_[1:band + 1]  # 正频率 1..p/H (对称)
        r = res[lowk]
        out[f'H{H}'] = dict(L=float(L), Ldir=float(Ldir), band_max=float(np.abs(r).max()),
                            band_L1=float(np.abs(r).mean()), band_L2=float(np.sqrt((r * r).mean())))
    return out

def main():
    t0 = time.time()
    args = sys.argv[1:]
    if args and args[0] == 'calib':
        plist = [3001]
        nrand = 2
    else:
        plist = primes_in(3000, 4200)
        plist += [10007, 30011]
        nrand = 3
    rng = np.random.default_rng(20260801)
    agg = {}
    for i, p in enumerate(plist):
        H_list = [8, 32, 128, int(math.isqrt(p))]
        key = orbit_keys(p)
        C = autocorr_Cp(key, p)
        d = analyze(C, p, H_list)
        # 镜像随机基准
        base = []
        for j in range(nrand):
            kb = mirror_random_key(p, rng)
            Cb = autocorr_Cp(kb, p)
            base.append(analyze(Cb, p, H_list))
        row = {'p': p, 'data': d, 'base': base}
        agg[p] = row
        # 进度 + 单素数摘要
        b0 = base[0]
        msg = f"[{i+1}/{len(plist)}] p={p} E^pi/p={d['Epi']/p:.3f} parseval={'OK' if d['parseval_ok'] else 'FAIL'}"
        for H in H_list:
            dd, bb = d[f'H{H}'], np.mean([b[f'H{H}']['L'] for b in base])
            chk = 'OK' if abs(dd['L'] - dd['Ldir']) < 1e-6 * max(1, abs(dd['L'])) + 1e-6 else 'MISMATCH'
            msg += f" | H={H}: L={dd['L']:.1f} (base {bb:.1f}) id:{chk}"
        print(msg, flush=True)
    # 聚合
    print("\n=== AGGREGATE (data/base ratios) ===", flush=True)
    for H_ix, Hname in enumerate(['H8', 'H32', 'H128']):
        rL, rM, rL2 = [], [], []
        for p, row in agg.items():
            if Hname not in row['data']: continue
            d = row['data'][Hname]
            bL = np.mean([b[Hname]['L'] for b in row['base']])
            bM = np.mean([b[Hname]['band_max'] for b in row['base']])
            b2 = np.mean([b[Hname]['band_L2'] for b in row['base']])
            if bL: rL.append(d['L'] / bL)
            if bM: rM.append(d['band_max'] / bM)
            if b2: rL2.append(d['band_L2'] / b2)
        print(f"{Hname}: L ratio mean={np.mean(rL):.4f}+-{np.std(rL):.4f} | band_max ratio={np.mean(rM):.4f} | band_L2 ratio={np.mean(rL2):.4f}", flush=True)
    print(f"done in {time.time()-t0:.1f}s", flush=True)

if __name__ == '__main__':
    main()
