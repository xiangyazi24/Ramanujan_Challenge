#!/usr/bin/env python3
"""CRON_qM_audit24.py — Q6509/Q6517 AUDIT-2 与 AUDIT-4.

AUDIT-2 (定理审计, 应零失败): p∈{13,17,29}, 0≤v≤p−2:
  p³c[p+v] ≡ (p³c[p])·K_v(p)  (mod p^7),
  K_0=1, K_1(T)=P(T)/(T+1)³, (T+v+1)³K_{v+1}=P(T+v)K_v−(T+v)³K_{v−1}  [精确有理数求值 T=p]
  附: mod p³ jet 律 p³c[p+v] ≡ b_v+2pD_v+p²E_v; mod p⁴ 加 p³(F_v+κ_p b_v) — 用 K_v(T) 系数 (v≤2 已知值对表).
AUDIT-4 (分叉实验): p∈{101,211}: jet 塔 D,E,F mod p (K_v 递推于 F_p[T]/(T⁴));
  (a) 碰撞对 (r,s) 上 jet_D/E/F(r,s)=b_rD_s−b_sD_r 等的消失率 vs 随机基线 ~1/p;
  (b) c 向量是否落在 span{b,D,E,F} (秩测试, v 限安全区 v<p/2).
"""
from fractions import Fraction
from math import factorial
from collections import defaultdict

def solve_pair_frac(N):
    b = [Fraction(1), Fraction(5)]; c = [Fraction(0), Fraction(1)]
    for n in range(2, N + 1):
        A = 34*n**3 - 51*n**2 + 27*n - 5
        B = (n - 1)**3
        b.append((A*b[n-1] - B*b[n-2]) / n**3)
        c.append((A*c[n-1] - B*c[n-2]) / n**3)
    return b, c

def vp_frac(fr, p):
    if fr == 0: return 10**9
    n, d = fr.numerator, fr.denominator
    v = 0
    while n % p == 0: n //= p; v += 1
    while d % p == 0: d //= p; v -= 1
    return v

fails = 0
def chk(c, m):
    global fails
    if not c:
        fails += 1
        print(f"  [FAIL] {m}", flush=True)

def P_(x):  # works for Fraction / int / poly-coeff lists via caller
    return 34*x**3 + 51*x**2 + 27*x + 5

print("== AUDIT-2: p³c[p+v] ≡ (p³c[p])·K_v(p) mod p^7 ==", flush=True)
for p in (13, 17, 29):
    b, c = solve_pair_frac(2*p + 2)
    Km1, K0v = None, Fraction(1)                    # K_0=1
    K = [Fraction(1), P_(Fraction(p)) / (Fraction(p) + 1)**3]
    for v in range(1, p - 1):
        T = Fraction(p)
        Kn = (P_(T + v)*K[v] - (T + v)**3*K[v-1]) / (T + v + 1)**3
        K.append(Kn)
    lead = c[p] * p**3
    nok = 0
    for v in range(0, p - 1):
        diff = c[p + v]*p**3 - lead*K[v]
        if vp_frac(diff, p) < 7:
            chk(False, f"p={p} v={v}: vp={vp_frac(diff,p)} < 7")
        else:
            nok += 1
    print(f"  p={p}: {nok}/{p-1} 个 v 全过 mod p^7", flush=True)

print("\n== AUDIT-2 附: jet 律 mod p³/p⁴ (v=0,1,2 已知 jet 值) ==", flush=True)
JET = {0: (0, 0, 0), 1: (6, 0, Fraction(-7)), 2: (105, 72, Fraction(-1011, 8))}
from sympy import bernoulli, Rational
for p in (13, 17, 29):
    b, c = solve_pair_frac(p + 3)
    Bp3 = Rational(bernoulli(p - 3))
    kappa = int(Bp3.p) * pow(int(Bp3.q), -1, p) % p
    kappa = (p - 2) * pow(3, -1, p) * kappa % p     # −(2/3)B_{p−3}
    for v in (0, 1, 2):
        D, E, F = JET[v]
        lhs = c[p + v] * p**3
        rhs3 = b[v] + 2*p*D + p*p*E
        chk(vp_frac(lhs - rhs3, p) >= 3, f"p={p} v={v}: mod p³ jet 律")
        rhs4 = rhs3 + p**3*(F + kappa*b[v])
        chk(vp_frac(lhs - rhs4, p) >= 4, f"p={p} v={v}: mod p⁴ jet+κ 律")
    print(f"  p={p}: v=0,1,2 mod p³/p⁴ 全过 (κ_p={kappa})", flush=True)

print("\n== AUDIT-4: jet 消失率与 span 测试 ==", flush=True)
def jets_mod_p(p, vmax):
    """K_v(T) in F_p[T]/(T^4): 返回 (b, D, E, F) 数组 mod p (b=K_v(0))."""
    def pmul(a, bb):
        out = [0]*4
        for i in range(4):
            for j in range(4 - i):
                out[i+j] = (out[i+j] + a[i]*bb[j]) % p
        return out
    def cube_lin(s):   # (T+s)^3 as coeff list mod p
        return [pow(s, 3, p), 3*s*s % p, 3*s % p, 1]
    def inv_series(a):  # invert a[0]+a[1]T+... mod T^4, a[0]≠0
        i0 = pow(a[0], -1, p)
        out = [i0, 0, 0, 0]
        for k in range(1, 4):
            s = 0
            for j in range(1, k + 1):
                s = (s + a[j]*out[k-j]) % p
            out[k] = (-i0 * s) % p
        return out
    def P_series(s):    # P(T+s) mod T^4
        return [P_(s) % p, (102*s*s + 102*s + 27) % p, (102*s + 51) % p, 34]
    K = [[1, 0, 0, 0], pmul(P_series(0), inv_series(cube_lin(1)))]
    for v in range(1, vmax):
        num = [(x - y) % p for x, y in zip(pmul(P_series(v), K[v]), pmul(cube_lin(v), K[v-1]))]
        K.append(pmul(num, inv_series(cube_lin(v + 1))))
    b_ = [k[0] for k in K]
    D_ = [k[1] * pow(2, -1, p) % p for k in K]      # K'(0)=2D ⟹ D=[T¹]/2
    E_ = [k[2] for k in K]
    F_ = [k[3] for k in K]
    return b_, D_, E_, F_

def orbit_bc_mod(p):
    b = [1, 5 % p]; c = [0, 1]
    for n in range(2, p - 1):
        A = (34*n**3 - 51*n**2 + 27*n - 5) % p
        B = (n - 1)**3 % p
        inv = pow(n**3 % p, -1, p)
        b.append((A*b[n-1] - B*b[n-2]) * inv % p)
        c.append((A*c[n-1] - B*c[n-2]) * inv % p)
    return b, c

for p in (101, 211):
    vmax = (p - 1)//2 - 1                            # 安全区
    bK, D, E, F = jets_mod_p(p, vmax)
    b, c = orbit_bc_mod(p)
    chk(all(bK[v] == b[v] % p for v in range(vmax)), f"p={p}: K_v(0)=b_v 自检")
    # 碰撞对 (安全区内)
    pos = defaultdict(list)
    for n in range(vmax):
        key = (1, c[n] * pow(b[n], -1, p) % p) if b[n] % p else (0, 1)
        pos[key].append(n)
    pairs = [(l[i], l[j]) for l in pos.values() for i in range(len(l)) for j in range(i+1, len(l))]
    nD = sum(1 for r, s in pairs if (b[r]*D[s] - b[s]*D[r]) % p == 0)
    nE = sum(1 for r, s in pairs if (b[r]*E[s] - b[s]*E[r]) % p == 0)
    nF = sum(1 for r, s in pairs if (b[r]*F[s] - b[s]*F[r]) % p == 0)
    # span 测试: c ∈ span{b,D,E,F}? 增广秩比较 (模 p 高斯消元)
    def rank(rows):
        M = [r[:] for r in rows]; R = len(M); C = len(M[0]); rk = 0
        for col in range(C):
            piv = next((i for i in range(rk, R) if M[i][col] % p), None)
            if piv is None: continue
            M[rk], M[piv] = M[piv], M[rk]
            inv = pow(M[rk][col], -1, p)
            M[rk] = [x*inv % p for x in M[rk]]
            for i in range(R):
                if i != rk and M[i][col] % p:
                    f = M[i][col]
                    M[i] = [(x - f*y) % p for x, y in zip(M[i], M[rk])]
            rk += 1
            if rk == R: break
        return rk
    base = [bK[:vmax], D[:vmax], E[:vmax], F[:vmax]]
    r4 = rank(base); r5 = rank(base + [c[:vmax]])
    print(f"  p={p}: 碰撞对={len(pairs)}; jet 消失率 D={nD}/{len(pairs)} E={nE}/{len(pairs)} "
          f"F={nF}/{len(pairs)} (随机基线≈{len(pairs)/p:.2f}); span: rank(b,D,E,F)={r4}, "
          f"加 c 后={r5} ⟹ c{'∈' if r5==r4 else '∉'}span", flush=True)

print(f"\nFAILS={fails}", flush=True)
