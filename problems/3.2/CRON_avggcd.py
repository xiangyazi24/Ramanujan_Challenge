#!/usr/bin/env python3
"""CRON_avggcd.py — [AVG-GCD_1/2] dyadic 机器 (Q6521 §6.2).

gap 多项式 (AJ.1, 已符号验证 h<=12):
  N_1 = 1, N_2(r) = P(r+1), N_h(r) = P(r+1)*N_{h-1}(r+1) - (r+2)^6*N_{h-2}(r+2)
  p | N_h(r) <=> pi(r)=pi(r+h);  deg N_h = 3(h-1)
  偶 h: M_h = N_h/(2r+h+1); 奇 h: M_h = N_h
输出: dyadic H: weighted=sum deg gcd(M_h,M_k), support=#{gcd 非平凡}
  归一 weighted/H^{5/2}, support/H^{3/2} ([AVG-GCD_1/2] 目标口径)
先做轨道交叉验证: N_h 的 F_p 根数(distinct) vs 轨道 R_h, h<=10 必须全等。
"""
import sys, time

def poly_trim(a):
    while a and a[-1] == 0: a.pop()
    return a

def poly_add(a, b, p):
    n = max(len(a), len(b)); r = [0]*n
    for i, x in enumerate(a): r[i] = x
    for i, x in enumerate(b): r[i] = (r[i] + x) % p
    return poly_trim(r)

def poly_scale(a, s, p):
    return poly_trim([x*s % p for x in a])

def poly_mul(a, b, p):
    if not a or not b: return []
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i+j] = (r[i+j] + x*y) % p
    return poly_trim(r)

def poly_shift(a, s, p):
    """a(r+s)"""
    # Horner in (r+s): evaluate coefficients via binomial — do repeated synthetic
    r = list(a)
    n = len(r)
    if n == 0: return r
    # Taylor shift by s: classic algorithm
    for i in range(n-1):
        for j in range(n-2, i-1, -1):
            r[j] = (r[j] + s*r[j+1]) % p
    return poly_trim(r)

def poly_mod(a, b, p):
    a = list(a); db = len(b)-1; inv = pow(b[-1], p-2, p)
    while len(a)-1 >= db and a:
        c = a[-1]*inv % p; da = len(a)-1
        for i in range(db+1):
            a[da-db+i] = (a[da-db+i] - c*b[i]) % p
        poly_trim(a)
    return a

def poly_gcd(a, b, p):
    a, b = list(a), list(b)
    while b:
        a, b = b, poly_mod(a, b, p)
    return a

def poly_divexact(a, b, p):
    """exact division a/b"""
    a = list(a); q = [0]*(len(a)-len(b)+1); db = len(b)-1; inv = pow(b[-1], p-2, p)
    while len(a)-1 >= db and a:
        c = a[-1]*inv % p; da = len(a)-1
        q[da-db] = c
        for i in range(db+1):
            a[da-db+i] = (a[da-db+i] - c*b[i]) % p
        poly_trim(a)
    assert not a, "not exact"
    return poly_trim(q)

def gap_polys(p, hmax):
    """N_h for h=1..hmax over F_p, h-方向递推 (Casoratian 第一性推导):
       g(s)=b_r c_s - b_s c_r 满足 s-递推, g(r)=0 => N_h := 清分母连分式解,
       N_1=1, N_2=P(r+1), N_{h+1}(r) = P(r+h)*N_h(r) - (r+h)^6*N_{h-1}(r).
       P(x)=34x^3+51x^2+27x+5, coeffs low->high: [5,27,51,34]"""
    Pp = [5 % p, 27 % p, 51 % p, 34 % p]
    N = {1: [1]}
    N[2] = poly_shift(Pp, 1, p)  # P(r+1)
    for h in range(2, hmax):
        Ph = poly_shift(Pp, h, p)  # P(r+h)
        rh = [h % p, 1]            # (r+h)
        rh6 = [1]
        for _ in range(6): rh6 = poly_mul(rh6, rh, p)
        t1 = poly_mul(Ph, N[h], p)
        t2 = poly_mul(rh6, N[h-1], p)
        N[h+1] = poly_add(t1, poly_scale(t2, p-1, p), p)
    return N

def orbit_Rh(p, hmax):
    b = [0]*(p-1); c = [0]*(p-1)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    for n in range(1, p-2):
        Pn = (34*n**3 + 51*n*n + 27*n + 5) % p
        n3 = n**3 % p
        inv = pow((n+1)**3 % p, p-2, p)
        b[n+1] = (Pn*b[n] - n3*b[n-1]) % p * inv % p
        c[n+1] = (Pn*c[n] - n3*c[n-1]) % p * inv % p
    key = [(b[n]*pow(c[n], p-2, p)) % p if c[n] else p for n in range(p-1)]
    R = {}
    for h in range(1, hmax+1):
        R[h] = sum(1 for r in range(p-1-h) if key[r] == key[r+h])
    return R

def count_roots(poly, p):
    """distinct roots in F_p: deg gcd(x^p - x, poly). x^p mod poly via repeated squaring."""
    if len(poly) <= 1: return 0
    # compute x^p mod poly
    def mulmod(a, b): return poly_mod(poly_mul(a, b, p), poly, p)
    result = [1]; base = [0, 1]; e = p
    while e:
        if e & 1: result = mulmod(result, base)
        base = mulmod(base, base); e >>= 1
    xp_minus_x = poly_add(result, poly_scale([0, 1], p-1, p), p)
    g = poly_gcd(poly, xp_minus_x, p)
    return len(g)-1

def main():
    t0 = time.time()
    primes = [3001, 4001, 10007]
    Hs = [8, 16, 32]
    hmax = 64
    for p in primes:
        N = gap_polys(p, hmax)
        # 交叉验证 h<=10
        R = orbit_Rh(p, 10)
        ok = True
        for h in range(2, 11):
            nr = count_roots(N[h], p)
            # 有效区间: 根 r 须在 [0, p-2-h]; 多项式根在全 F_p — R_h 只数区间内。近似核对: nr >= R[h]
            if nr < R[h]:
                print(f"  VALIDATION FAIL p={p} h={h}: poly roots {nr} < orbit R_h {R[h]}", flush=True)
                ok = False
        print(f"p={p} validation (poly roots >= orbit R_h, h=2..10): {'OK' if ok else 'FAIL'}", flush=True)
        # M_h
        M = {}
        for h in range(2, hmax+1):
            if h % 2 == 0:
                lin = [(h+1) % p, 2]  # 2r+h+1
                try:
                    M[h] = poly_divexact(N[h], lin, p)
                except AssertionError:
                    M[h] = N[h]
                    print(f"  NOTE p={p} h={h}: (2r+h+1) does not divide N_h — using N_h", flush=True)
            else:
                M[h] = N[h]
        for H in Hs:
            weighted = 0; support = 0; npairs = 0
            for h in range(H+1, 2*H+1):
                for k in range(h+1, 2*H+1):
                    g = poly_gcd(M[h], M[k], p)
                    d = len(g)-1
                    npairs += 1
                    if d > 0:
                        weighted += d; support += 1
            print(f"  H={H}: pairs={npairs} weighted={weighted} support={support} | weighted/H^2.5={weighted/H**2.5:.4f} support/H^1.5={support/H**1.5:.4f}", flush=True)
    print(f"done {time.time()-t0:.1f}s", flush=True)

if __name__ == '__main__':
    main()
