#!/usr/bin/env python3
"""Problem 2.7: Numerically identify the gauge between (KT)_n and q_n.

(KT)_n = (1/256^n) sum_{k=0}^{2n} C(2n,k)(-2)^{2n-k} T_k
where T_k is the Cooper level-11 sequence.
q_n is the challenge denominator sequence.

We compute both and look for g_n = q_n / (KT)_n.
"""
from mpmath import mp, mpf, binomial, fac, nstr, log, gamma, pi, zeta

mp.dps = 150

# ---- Cooper level-11 sequence ----
# (n+1)^3 T_{n+1} = (2n+1)(10n^2+10n+4) T_n + n(-56n^2-8) T_{n-1}
#                    + 22 n(2n-1)(n-1) T_{n-2}
# T_0=1, T_1=4, T_2=28

def cooper_seq(N):
    T = [mpf(0)] * (N+1)
    T[0] = 1
    T[1] = 4
    T[2] = 28
    for n in range(2, N):
        a = (2*n+1)*(10*n**2 + 10*n + 4)
        b = n*(-56*n**2 - 8)
        c = 22*n*(2*n-1)*(n-1)
        T[n+1] = (a*T[n] + b*T[n-1] + c*T[n-2]) / (n+1)**3
    return T

# ---- Centered binomial transform ----
def centered_binom_transform(T, N_out):
    """(KT)_n = (1/256^n) sum_{k=0}^{2n} C(2n,k)(-2)^{2n-k} T_k"""
    KT = []
    for n in range(N_out):
        s = mpf(0)
        for k in range(2*n+1):
            s += binomial(2*n, k) * mpf(-2)**(2*n-k) * T[k]
        s /= mpf(256)**n
        KT.append(s)
    return KT

# ---- Problem 2.7 recurrence ----
def A27(n):
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

def B27(n):
    return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)

def C27(n):
    return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)

def D27(n):
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

def q27_seq(N):
    """Compute Problem 2.7 denominator sequence."""
    q = [mpf(0)] * (N+1)
    # Initial values from the summary
    q[0] = mpf('-215040420000')
    q[1] = mpf('-167282265043404') / 905
    q[2] = mpf('-964185327658080') / 6071
    for n in range(2, N):
        q[n+1] = (B27(n)/A27(n))*q[n] - (C27(n-1)/A27(n-1))*q[n-1] + (D27(n-2)/A27(n-2))*q[n-2]
    return q

# ---- Compute both ----
N_cooper = 120
N_out = 40

print("Computing Cooper sequence...")
T = cooper_seq(N_cooper)
print(f"T[0..5] = {[int(T[i]) for i in range(6)]}")

print("\nComputing centered binomial transform...")
KT = centered_binom_transform(T, N_out)
print(f"(KT)[0..5]:")
for i in range(6):
    print(f"  (KT)[{i}] = {nstr(KT[i], 30)}")

print("\nComputing Problem 2.7 sequence...")
q = q27_seq(N_out)
print(f"q[0..5]:")
for i in range(6):
    print(f"  q[{i}] = {nstr(q[i], 30)}")

# ---- Gauge g_n = q_n / (KT)_n ----
print("\n=== Gauge g_n = q_n / (KT)_n ===")
g = []
for n in range(N_out):
    if KT[n] != 0:
        gn = q[n] / KT[n]
        g.append(gn)
        if n < 15:
            print(f"  g[{n}] = {nstr(gn, 30)}")
    else:
        g.append(None)

# ---- Gauge ratio g_{n+1}/g_n ----
print("\n=== Gauge ratio g_{n+1}/g_n ===")
for n in range(min(25, len(g)-1)):
    if g[n] is not None and g[n+1] is not None and g[n] != 0:
        r = g[n+1] / g[n]
        print(f"  g[{n+1}]/g[{n}] = {nstr(r, 30)}")

# ---- Try to identify the ratio as (n+a)/(n+b) ----
print("\n=== (g[n+1]/g[n] - 1) * n ===")
for n in range(1, min(25, len(g)-1)):
    if g[n] is not None and g[n+1] is not None and g[n] != 0:
        r = g[n+1] / g[n]
        val = (r - 1) * n
        print(f"  n={n}: {nstr(val, 30)}")

# ---- Try Pochhammer form: g_n = Γ(n+a)/Γ(n+b) ----
print("\n=== Try g_n = C * Γ(n+a) / Γ(n+b) with a-b=3/2 ===")
# g_{n+1}/g_n = (n+a)/(n+b)
# At n=1: (1+a)/(1+b) = g[1]/g[0]
# With a = b + 3/2, we get (1+b+3/2)/(1+b) = (b+5/2)/(b+1)
if g[0] is not None and g[1] is not None and g[0] != 0:
    r01 = g[1] / g[0]
    print(f"g[1]/g[0] = {nstr(r01, 30)}")
    # (b+5/2)/(b+1) = r01 => b+5/2 = r01*(b+1) => b(1-r01) = r01 - 5/2
    # b = (r01 - 5/2) / (1 - r01) = (5/2 - r01) / (r01 - 1)
    b_val = (mpf(5)/2 - r01) / (r01 - 1)
    a_val = b_val + mpf(3)/2
    print(f"Solved: a = {nstr(a_val, 30)}, b = {nstr(b_val, 30)}")

    # Verify at other indices
    print("\nVerification:")
    for n in range(2, min(20, len(g)-1)):
        if g[n] is not None and g[n+1] is not None and g[n] != 0:
            predicted = (n + a_val) / (n + b_val)
            actual = g[n+1] / g[n]
            rel_err = float(abs((predicted - actual) / actual))
            print(f"  n={n}: predicted={nstr(predicted,15)}, actual={nstr(actual,15)}, rel_err={rel_err:.6e}")

# ---- More general: try product of Pochhammer ratios ----
# g_n = C * prod_{j} Γ(n+a_j) / Γ(n+b_j)
# g_{n+1}/g_n = prod (n+a_j) / prod (n+b_j)
# Expand: prod (n+a_j) / prod (n+b_j) ≈ 1 + (sum a_j - sum b_j)/n + ...
# We need sum(a_j - b_j) = 3/2

# Let's compute higher-order corrections
print("\n=== Higher-order expansion of gauge ratio ===")
print("g_{n+1}/g_n = 1 + c_1/n + c_2/n^2 + ...")
for k in range(1, 6):
    print(f"\n  c_{k} estimates from various n:")
    for n in [10, 15, 20, 25, 30]:
        if n < len(g)-1 and g[n] is not None and g[n+1] is not None and g[n] != 0:
            r = g[n+1] / g[n]
            # Successively subtract known terms
            remainder = r - 1
            for j in range(1, k):
                # We'd need to know c_j; just print the raw expansion
                pass
            val = remainder * n**k
            print(f"    n={n}: (g ratio - 1) * n^{k} = {nstr(val, 20)}")
