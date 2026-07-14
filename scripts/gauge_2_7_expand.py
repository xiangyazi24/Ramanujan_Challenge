#!/usr/bin/env python3
"""Problem 2.7: Extract asymptotic expansion of gauge ratio."""
from mpmath import mp, mpf, binomial, nstr

mp.dps = 150

# ---- Cooper level-11 sequence ----
def cooper_seq(N):
    T = [mpf(0)] * (N+1)
    T[0] = 1; T[1] = 4; T[2] = 28
    for n in range(2, N):
        a = (2*n+1)*(10*n**2 + 10*n + 4)
        b = n*(-56*n**2 - 8)
        c = 22*n*(2*n-1)*(n-1)
        T[n+1] = (a*T[n] + b*T[n-1] + c*T[n-2]) / (n+1)**3
    return T

def centered_binom_transform(T, N_out):
    KT = []
    for n in range(N_out):
        s = mpf(0)
        for k in range(2*n+1):
            s += binomial(2*n, k) * mpf(-2)**(2*n-k) * T[k]
        s /= mpf(256)**n
        KT.append(s)
    return KT

def A27(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def B27(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C27(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D27(n): return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

def q27_seq(N):
    q = [mpf(0)] * (N+1)
    q[0] = mpf('-215040420000')
    q[1] = mpf('-167282265043404') / 905
    q[2] = mpf('-964185327658080') / 6071
    for n in range(2, N):
        q[n+1] = (B27(n)/A27(n))*q[n] - (C27(n-1)/A27(n-1))*q[n-1] + (D27(n-2)/A27(n-2))*q[n-2]
    return q

N_cooper = 200
N_out = 60
T = cooper_seq(N_cooper)
KT = centered_binom_transform(T, N_out)
q = q27_seq(N_out)

g = [q[n] / KT[n] if KT[n] != 0 else None for n in range(N_out)]

# Extract expansion: g_{n+1}/g_n = 1 + c_1/n + c_2/n^2 + c_3/n^3 + ...
# Richardson extrapolation: compute (r-1)*n for many n values, then extrapolate
ratios = []
for n in range(1, N_out-1):
    if g[n] and g[n+1] and g[n] != 0:
        ratios.append((n, g[n+1]/g[n]))

# c_1: limit of (r-1)*n as n→∞
print("=== c_1 extraction (limit of (r-1)*n) ===")
c1_vals = [(n, (r-1)*n) for n, r in ratios]
for n, v in c1_vals[-10:]:
    print(f"  n={n}: {nstr(v, 30)}")
c1 = mpf(3)/2
print(f"\nc_1 = 3/2 = {nstr(c1, 30)}")

# c_2: limit of ((r-1)*n - c_1) * n as n→∞
print("\n=== c_2 extraction ===")
c2_vals = [(n, ((r-1)*n - c1)*n) for n, r in ratios]
for n, v in c2_vals[-10:]:
    print(f"  n={n}: {nstr(v, 30)}")
# Richardson: extrapolate pairs
if len(c2_vals) >= 2:
    n1, v1 = c2_vals[-2]
    n2, v2 = c2_vals[-1]
    c2_rich = (n2*v2 - n1*v1) / (n2 - n1)
    print(f"\nRichardson extrapolation: c_2 ≈ {nstr(c2_rich, 20)}")

# c_3: limit of (((r-1)*n - c_1)*n - c_2) * n
print("\n=== c_3 extraction ===")
# Use the Richardson-extrapolated c_2
c2 = c2_rich
c3_vals = [(n, (((r-1)*n - c1)*n - c2)*n) for n, r in ratios]
for n, v in c3_vals[-10:]:
    print(f"  n={n}: {nstr(v, 30)}")

# Try recognizing c_2 as a simple fraction
print("\n=== Rationality test for c_2 ===")
# c_2 candidates: 1/8, 3/8, 5/8, 7/8, 1/4, 3/4, etc.
for num in range(-20, 20):
    for den in [1, 2, 4, 8, 16]:
        if den > 0:
            cand = mpf(num) / den
            if abs(cand - c2_rich) < 0.02:
                print(f"  {num}/{den} = {nstr(cand, 10)}, diff = {nstr(abs(cand - c2_rich), 10)}")

# Also check: is g_n a product of multiple Pochhammer ratios?
# g_{n+1}/g_n = prod_{j} (n+a_j) / prod_{j} (n+b_j)
# If there are 2 pairs: (n+a1)(n+a2) / ((n+b1)(n+b2))
# = 1 + (a1+a2-b1-b2)/n + (a1a2 - b1b2 + (a1+a2)(sum) - ...)/n^2 + ...
# sum(a) - sum(b) = c_1 = 3/2
# And c_2 = (a1*a2 - b1*b2)/? No, let me expand properly.
# Let f(n) = (1 + a1/n)(1 + a2/n) / ((1 + b1/n)(1 + b2/n))
# ≈ 1 + (a1+a2-b1-b2)/n + (a1a2 - b1b2 - (b1+b2)(a1+a2-b1-b2) + ...)/n^2

# Actually let's try: sum a_j = S, sum b_j = S - 3/2
# Then c_2 = (S(S-1) - sum a_j*a_k) ... this gets complicated.
# Just see if c_2 = 3/2 * something simple.

print(f"\nc_2 / c_1 = {nstr(c2_rich / c1, 20)}")
print(f"c_2 / (c_1^2) = {nstr(c2_rich / c1**2, 20)}")
print(f"c_2 - c_1^2/2 = {nstr(c2_rich - c1**2/2, 20)}")
print(f"c_2 - 9/8 = {nstr(c2_rich - mpf(9)/8, 20)}")
print(f"c_2 + 3/8 = {nstr(c2_rich + mpf(3)/8, 20)}")
print(f"4*c_2 = {nstr(4*c2_rich, 20)}")
print(f"8*c_2 = {nstr(8*c2_rich, 20)}")
