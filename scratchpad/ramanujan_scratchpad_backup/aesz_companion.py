"""
Search for the ζ-companion of the AESZ #209 period.

The AESZ #209 period is A_n = C(2n,n) * Σ C(n,k)^2 * C(n+k,n) * C(n+2k,n).
The P2.7 sequences q_n, p_n are related to A_n and its ζ-companion by:
  q_n = N * g(n) * A_n  (up to normalization N)
  p_n = N * g(n) * B_n  (B_n = ζ-companion)

where g(n) = n^2 + 105/22*n + 5399/946 is the Ore gauge.

Step 1: Find N by comparing q_n with g(n)*A_n
Step 2: Compute B_n = p_n / (N * g(n))
Step 3: Check if B_n can be expressed as a harmonic-decorated sum
"""
from fractions import Fraction
from math import comb
import mpmath

mpmath.mp.dps = 100

def A_period(n):
    """AESZ #209 holomorphic period"""
    s = Fraction(0)
    for k in range(n+1):
        s += Fraction(comb(n,k)**2 * comb(n+k,n) * comb(n+2*k,n))
    return Fraction(comb(2*n,n)) * s

def g(n):
    """Ore gauge: g(n) = n^2 + 105/22*n + 5399/946"""
    return Fraction(n*n) + Fraction(105, 22) * n + Fraction(5399, 946)

# P2.7 initial conditions
p = [Fraction(-612218384750), Fraction(-9525021973931919, 18100), Fraction(-29561828382772029, 65380)]
q = [Fraction(-215040420000), Fraction(-167282265043404, 905), Fraction(-964185327658080, 6071)]

# Step 1: Compute normalization N = q_n / (g(n) * A_n)
print("=== Step 1: Normalization ===")
for n in range(3):
    A_n = A_period(n)
    g_n = g(n)
    ratio = q[n] / (g_n * A_n)
    print(f"n={n}: A_{n} = {A_n}, g({n}) = {g_n}")
    print(f"  q_{n} / (g_{n} * A_{n}) = {ratio}")
    print(f"  = {float(ratio):.6f}")
    print()

N_values = [q[n] / (g(n) * A_period(n)) for n in range(3)]
print(f"N values: {[float(x) for x in N_values]}")
if N_values[0] == N_values[1] == N_values[2]:
    print("✓ All ratios match! N =", N_values[0])
    N = N_values[0]
else:
    print("✗ Ratios don't match")
    print(f"N[0] = {N_values[0]}")
    print(f"N[1] = {N_values[1]}")
    print(f"N[2] = {N_values[2]}")
    print()

    # Maybe the gauge is different. Try q_n / A_n directly
    print("Trying without gauge (q_n / A_n):")
    for n in range(3):
        r = q[n] / A_period(n)
        print(f"  n={n}: {r} = {float(r):.6f}")

    # Maybe the AESZ recurrence needs different initial conditions
    # Let's compute the AESZ recurrence from the P2.7 recurrence
    # v_n = q_n / g(n) should satisfy the AESZ recurrence
    print("\nComputing v_n = q_n / g(n):")
    for n in range(3):
        v = q[n] / g(n)
        A_n = A_period(n)
        r = v / A_n
        print(f"  n={n}: v_{n} = {v}, v/A = {r} = {float(r):.6f}")
    v_ratios = [q[n] / (g(n) * A_period(n)) for n in range(3)]
    print(f"  v/A ratios: {[float(r) for r in v_ratios]}")

# Step 2: Regardless, compute B_n = p_n / q_n * A_n or similar
print("\n=== Step 2: p/q ratio check ===")
zeta_sum = mpmath.zeta(2) + mpmath.zeta(3)
print(f"ζ(2)+ζ(3) = {float(zeta_sum):.15f}")
for n in range(3):
    print(f"n={n}: p/q = {float(p[n]/q[n]):.15f}")

# Step 3: Extend sequences using recurrence and check more terms
print("\n=== Step 3: Extended sequences ===")

def A_coeff(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def B_coeff(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_coeff(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_coeff(n): return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

def extend(u0, u1, u2, n_max):
    seq = [Fraction(u0), Fraction(u1), Fraction(u2)]
    for n in range(2, n_max):
        u_new = (Fraction(B_coeff(n), A_coeff(n)) * seq[-1]
                - Fraction(C_coeff(n-1), A_coeff(n-1)) * seq[-2]
                + Fraction(D_coeff(n-2), A_coeff(n-2)) * seq[-3])
        seq.append(u_new)
    return seq

q_ext = extend(q[0], q[1], q[2], 8)
p_ext = extend(p[0], p[1], p[2], 8)

print("Convergence p_n/q_n → ζ(2)+ζ(3):")
for n in range(min(8, len(q_ext))):
    ratio = float(p_ext[n] / q_ext[n])
    print(f"  n={n}: p/q = {ratio:.15f}, diff = {ratio - float(zeta_sum):.6e}")

# Step 4: Compute v_n = q_n/g(n) and check against A_n
print("\n=== Step 4: v_n = q_n/g(n) vs A_n ===")
for n in range(min(6, len(q_ext))):
    v_n = q_ext[n] / g(n)
    A_n = A_period(n)
    if A_n != 0:
        ratio = v_n / A_n
        print(f"  n={n}: v/A = {ratio}")
    else:
        print(f"  n={n}: A=0")

# Step 5: If v_n is NOT proportional to A_n, try other gauges
print("\n=== Step 5: Trying R(n) = 946*g(n) = 946n²+4515n+5399 ===")
def R(n): return 946*n**2 + 4515*n + 5399
for n in range(min(6, len(q_ext))):
    v_n = q_ext[n] / Fraction(R(n))
    A_n = A_period(n)
    if A_n != 0:
        ratio = v_n / A_n
        print(f"  n={n}: q/(R*A) = {ratio}")

# Step 6: Maybe q_n itself satisfies AESZ #209 recurrence directly
# (the P2.7 recurrence IS the AESZ #209 recurrence, just different normalization)
# Let's check: does A_n satisfy the P2.7 recurrence?
print("\n=== Step 6: Does A_n satisfy P2.7 recurrence? ===")
A_vals = [A_period(n) for n in range(8)]
for n in range(2, 5):
    lhs = Fraction(A_coeff(n)) * A_vals[n+1] - Fraction(B_coeff(n)) * A_vals[n] + Fraction(C_coeff(n-1)) * A_vals[n-1] - Fraction(D_coeff(n-2)) * A_vals[n-2]
    print(f"  n={n}: residual = {lhs}")

# Step 7: Check A_n * g(n) with the recurrence
print("\n=== Step 7: Does g(n)*A_n satisfy P2.7 recurrence? ===")
gA_vals = [g(n) * A_period(n) for n in range(8)]
for n in range(2, 5):
    lhs = Fraction(A_coeff(n)) * gA_vals[n+1] - Fraction(B_coeff(n)) * gA_vals[n] + Fraction(C_coeff(n-1)) * gA_vals[n-1] - Fraction(D_coeff(n-2)) * gA_vals[n-2]
    print(f"  n={n}: residual = {lhs}")

print("\n=== Step 8: Factor analysis of q_0 ===")
import sympy
q0_abs = abs(int(q[0]))
print(f"|q_0| = {q0_abs}")
print(f"Factored: {sympy.factorint(q0_abs)}")
A0 = int(A_period(0))
g0 = g(0)
print(f"A_0 = {A0}, g(0) = {g0}")
print(f"g(0) * A_0 = {g0}")
print(f"q_0 / g(0) = {q[0] / g0}")
