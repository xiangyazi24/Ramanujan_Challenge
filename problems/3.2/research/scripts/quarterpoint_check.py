#!/usr/bin/env python3
"""Independent check of cron's quarter-point law:
p = 5 mod 24 => tau_{(p-1)/4} = 0 mod p (tau = sqrt(F));
p = 23 mod 24 => sigma_{(p-3)/4} = 0 mod p (sigma = sqrt(F/q), q = 1-34x+x^2);
p = 1, 19 mod 24 => no such vanishing. Also verify the half-integer recurrence for tau over Q."""
from fractions import Fraction as F

def apery(N):
    b = [1, 5]
    for n in range(1, N):
        num = (2*n+1)*(17*n*n+17*n+5)*b[n] - n**3*b[n-1]
        q, r = divmod(num, (n+1)**3); assert r == 0
        b.append(q)
    return b

def sqrt_series_mod(coeffs, N, p):
    s = [1] + [0]*(N-1); inv2 = pow(2, p-2, p)
    for n in range(1, N):
        acc = coeffs[n] if n < len(coeffs) else 0
        t = sum(s[i]*s[n-i] for i in range(1, n)) % p
        s[n] = (acc - t) % p * inv2 % p
    return s

def div_series_mod(a, den, N, p):
    c = [0]*N; d0inv = pow(den[0], p-2, p)
    for n in range(N):
        t = a[n] if n < len(a) else 0
        t = (t - sum(den[i]*c[n-i] for i in range(1, min(n, len(den)-1)+1))) % p
        c[n] = t * d0inv % p
    return c

def sieve(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2, int(n**.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]

b = apery(600)
stats = {5: [0,0], 23: [0,0], 1: [0,0], 19: [0,0]}
for p in sieve(1000):
    m = p % 24
    if m not in stats or p < 7: continue
    bm = [v % p for v in b[:p]]
    if m == 5:
        tau = sqrt_series_mod(bm, (p-1)//4 + 1, p)
        z = tau[(p-1)//4] == 0
    elif m == 23:
        q = [1, (-34) % p, 1]
        Fq = div_series_mod(bm, q, (p-3)//4 + 1, p)
        sig = sqrt_series_mod(Fq, (p-3)//4 + 1, p)
        z = sig[(p-3)//4] == 0
    else:  # 1, 19: check BOTH candidate positions do NOT vanish
        tau = sqrt_series_mod(bm, (p-1)//4 + 2, p)
        q = [1, (-34) % p, 1]
        Fq = div_series_mod(bm, q, (p-1)//4 + 2, p)
        sig = sqrt_series_mod(Fq, (p-1)//4 + 2, p)
        z = (tau[(p-1)//4] == 0) or (sig[(p-3)//4] == 0)
    stats[m][0 if z else 1] += 1
print("class p%24=5  (tau quarter zero): zero/nonzero =", stats[5])
print("class p%24=23 (sigma quarter zero): zero/nonzero =", stats[23])
print("class p%24=1  (should NOT vanish): zero/nonzero =", stats[1])
print("class p%24=19 (should NOT vanish): zero/nonzero =", stats[19])

# half-integer recurrence for tau over Q: 4(j+2)^2 tau_{j+2} = 2(68j^2+170j+107) tau_{j+1} - (2j+1)^2 tau_j
# build sqrt(F) over Q
N = 40
tauQ = [F(1)] + [F(0)]*(N-1)
for n in range(1, N):
    acc = F(b[n]) - sum(tauQ[i]*tauQ[n-i] for i in range(1, n))
    tauQ[n] = acc / 2
ok = all(4*(j+2)**2*tauQ[j+2] == 2*(68*j*j+170*j+107)*tauQ[j+1] - (2*j+1)**2*tauQ[j] for j in range(N-2))
print("half-integer recurrence for tau = sqrt(F) over Q (j <", N-2, "):", "VERIFIED" if ok else "FAIL")
