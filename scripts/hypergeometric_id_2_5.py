#!/usr/bin/env python3
"""
Problem 2.5: Identify the hypergeometric structure by analyzing
the CMF correction terms delta_N = P_N/Q_N - P_{N-1}/Q_{N-1}.

If G = sum_{N>=0} delta_N and delta_{N+1}/delta_N is a rational function
of N, then the series is a hypergeometric function and can be identified.
"""
from mpmath import mp, mpf, catalan, nstr, matrix, log10, fabs

mp.dps = 200

def M(n):
    n = mpf(n)
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)
    return matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])

A = matrix([[mpf(30921), mpf(-32972), mpf(8240)],
            [mpf(33750), mpf(-36000), mpf(9000)]])

# Compute convergents P_N/Q_N
print("=== CMF Correction Term Analysis ===")
print(f"G = {nstr(catalan, 50)}")

N_MAX = 70
prod = matrix([[1,0,0],[0,1,0],[0,0,1]])
convergents = []

for N in range(N_MAX):
    if N > 0:
        prod = prod * M(N-1)
    AM = A * prod
    # Use first column: P = AM[0,0], Q = AM[1,0]
    P = AM[0, 0]
    Q = AM[1, 0]
    if Q != 0:
        convergents.append((N, P, Q, P/Q))

# Correction terms delta_N = conv[N] - conv[N-1]
print("\n=== Correction terms delta_N ===")
deltas = []
for i in range(1, len(convergents)):
    N, P, Q, r = convergents[i]
    _, _, _, r_prev = convergents[i-1]
    delta = r - r_prev
    deltas.append((N, delta))

# Ratios delta_{N+1}/delta_N
print("\n=== Ratio delta_{N+1}/delta_N ===")
ratios = []
for i in range(1, len(deltas)):
    N, d = deltas[i]
    N_prev, d_prev = deltas[i-1]
    if d_prev != 0:
        ratio = d / d_prev
        ratios.append((N, ratio))
        if N <= 15 or N % 10 == 0:
            print(f"  N={N}: ratio = {nstr(ratio, 30)}")

# For a hypergeometric term, ratio should be a rational function of N:
# delta_{N+1}/delta_N = P(N)/Q(N) where P, Q are polynomials
# This means: ratio(N) * N^k → constant as N → infty
# The leading behavior tells us the "convergence rate"

print("\n=== Asymptotic analysis of ratio ===")
# Expect ratio → c as N → infty (where c is the convergence rate)
# c should be related to the Poincaré root ratio
for N, r in ratios[-10:]:
    print(f"  N={N}: ratio = {nstr(r, 30)}, |ratio| = {nstr(fabs(r), 15)}")

# The ratio should approach -16 * (subdominant/dominant) Poincaré root
# = -16 * (17-12√2) / (17+12√2) ... hmm, or simply 1/(convergence rate)

# Actually, for a CMF, the correction terms decay as |c_subdominant/c_dominant|^N
# Poincaré roots: -16, -16(17+12√2), -16(17-12√2)
# Convergence rate = c_sub / c_dom

from mpmath import sqrt
alpha_p = 17 + 12*sqrt(2)
alpha_m = 17 - 12*sqrt(2)
c0, cp, cm = -16, -16*alpha_p, -16*alpha_m

print(f"\nPoincaré roots: c0={nstr(c0,10)}, c+={nstr(cp,10)}, c-={nstr(cm,10)}")
print(f"Convergence rate c-/c0 = {nstr(cm/c0,15)} = {nstr(alpha_m,15)}")
print(f"c0/c+ = {nstr(c0/cp,15)} = {nstr(1/alpha_p,15)}")

# The convergent P_N/Q_N converges at rate |c_min/c_mid|
# For 3x3 CMF: convergence involves TWO subdominant modes
# The dominant column gives the faster convergence

# Let me check: ratio * N^k for various k
print("\n=== ratio * N^k test ===")
for k in range(15):
    vals = [(N, r * N**k) for N, r in ratios[-5:]]
    spread = max(fabs(v) for _, v in vals) / min(fabs(v) for _, v in vals) if vals else 0
    if spread < 1.1:  # nearly constant
        avg = sum(v for _, v in vals) / len(vals)
        print(f"  k={k}: nearly constant! avg = {nstr(avg, 20)}, spread = {nstr(spread, 5)}")
        break
    elif k <= 3 or spread < 2:
        avg = sum(v for _, v in vals) / len(vals)
        print(f"  k={k}: avg = {nstr(avg, 15)}, spread = {nstr(spread, 5)}")

# Extract the "term" structure: what is delta_N * (-16)^N * N^sigma?
print("\n=== Normalized correction terms ===")
# Expected: delta_N ~ C * convergence_rate^N * N^sigma
# Try: delta_N / (1/16)^N * N^sigma ... hmm convergence_rate is the ratio of
# two consecutive corrections, which we already computed.

# Let me compute delta_N * (c+/c-)^N to "flatten" the exponential decay
# Actually, delta_N decays exponentially. The ratio tells us the rate.
# From the data: ratio approaches some limit r_inf.
# Then delta_N ~ const * r_inf^N.

# Let me check r_inf:
print("\n=== Convergence of ratio to limit ===")
for i in range(-10, 0):
    N, r = ratios[i]
    print(f"  N={N}: ratio = {nstr(r, 40)}")

# Now try to identify r_inf as a recognizable number
r_last = ratios[-1][1]
print(f"\nr_inf ≈ {nstr(r_last, 30)}")
print(f"1/r_inf ≈ {nstr(1/r_last, 30)}")

# Check against known values
from mpmath import pi, e, euler, ln
candidates = {
    "1/16": mpf(1)/16,
    "-1/16": mpf(-1)/16,
    "(17-12√2)/16": alpha_m/16,
    "-(17-12√2)/16": -alpha_m/16,
    "alpha_m": alpha_m,
    "-alpha_m": -alpha_m,
    "1/alpha_p": 1/alpha_p,
    "-1/alpha_p": -1/alpha_p,
    "c-/c0": cm/c0,
    "c0/c+": c0/cp,
    "-1/(4+3√2)²": -1/(4+3*sqrt(2))**2,
}

for name, val in candidates.items():
    if fabs(r_last - val) < mpf(10)**(-10):
        print(f"  MATCH: r_inf = {name} = {nstr(val, 20)}, diff = {nstr(fabs(r_last-val), 10)}")
    elif fabs(fabs(r_last) - fabs(val)) < mpf(10)**(-5):
        print(f"  CLOSE: |r_inf| ≈ |{name}| = {nstr(fabs(val), 20)}")

# Also check: is delta_N a clean series term?
# Compute delta_N / delta_0 and see if it has Pochhammer structure
print("\n=== Term structure: delta_N / delta_0 ===")
d0 = deltas[0][1]
for i in range(min(20, len(deltas))):
    N, d = deltas[i]
    ratio = d / d0
    print(f"  N={N}: delta_N/delta_0 = {nstr(ratio, 25)}")
