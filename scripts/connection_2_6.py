#!/usr/bin/env python3
"""Problem 2.6: Verify the connection formula U(1) = ζ(2)+ζ(3)-2077/720.
The generating function U(x) = Σ u_n x^n satisfies a specific ODE.
We verify by computing U(1) two ways:
1. Direct summation (series)
2. High-precision Poincaré root analysis

Also: improve the series convergence via Richardson extrapolation."""
from mpmath import mp, mpf, zeta, nsum, inf, power, log, richardson, shanks

mp.dps = 60

# Recurrence for u_n
def compute_u(N):
    """Compute u_1, ..., u_N using the recurrence."""
    u = [mpf(0), mpf(-93)/4480, mpf(-117)/14000]  # u[0] unused, u[1], u[2]

    for n in range(3, N+1):
        A = 2*(n+3)**3*(2*n+5)*(3*n+5)
        B = (n+2)**2*(15*n**3 + 85*n**2 + 155*n + 93)
        C = (n+1)**3*(n+2)*(3*n+8)
        un = (B*u[-1] - C*u[-2]) / A
        u.append(un)

    return u

print("=== Problem 2.6: Connection formula verification ===")

target = zeta(2) + zeta(3)
target_sum = target - mpf(2077)/720
print(f"Target: ζ(2)+ζ(3) = {target}")
print(f"Target sum: Σu_j = ζ(2)+ζ(3)-2077/720 = {target_sum}")

# Compute many terms
N = 500
u = compute_u(N)

# Partial sums
partial_sums = []
s = mpf(0)
for j in range(1, N+1):
    s += u[j]
    partial_sums.append(s)

print(f"\nPartial sums S_N = Σ_{{j=1}}^N u_j:")
for n in [50, 100, 200, 300, 400, 500]:
    if n <= N:
        diff = partial_sums[n-1] - target_sum
        print(f"  S_{n:3d} - target = {mp.nstr(diff, 6, strip_zeros=False)}")

# Poincaré root analysis: u_{n+1}/u_n → 1/4
print(f"\nRatio u_{{n+1}}/u_n (should → 1/4 = 0.25):")
for n in [10, 20, 50, 100, 200, 500]:
    if n < N:
        ratio = u[n+1] / u[n]
        print(f"  n={n:3d}: u_{{n+1}}/u_n = {ratio}")

# Richardson extrapolation on partial sums
# Since u_n ~ C * (1/4)^n * n^α, the tail Σ_{j>N} u_j ~ C' * (1/4)^N * N^β
# Richardson extrapolation accelerates geometric convergence
print(f"\nRichardson extrapolation (Euler-Maclaurin acceleration):")
print(f"Using Shanks transformation on partial sums...")

# Shanks transformation (epsilon algorithm)
def shanks_transform(sums, order=5):
    """Shanks transformation to accelerate convergence."""
    n = len(sums)
    if n < 2*order + 1:
        return sums[-1]

    # Use the last 2*order+1 partial sums
    e = [[mpf(0)] * (2*order + 2) for _ in range(2*order + 2)]
    for i in range(2*order + 1):
        e[i][0] = mpf(0)
        e[i][1] = sums[n - 2*order - 1 + i]

    for j in range(2, 2*order + 1):
        for i in range(0, 2*order + 1 - j):
            diff = e[i+1][j-1] - e[i][j-1]
            if abs(diff) < mpf(10)**(-mp.dps + 5):
                e[i][j] = e[i+1][j-1]
            else:
                e[i][j] = e[i][j-2] + 1/diff

    # The best estimate is e[0][2*order] (even columns are the accelerated sums)
    return e[0][2*order]

# Apply Shanks to partial sums at various N
for N_use in [100, 200, 300, 400, 500]:
    if N_use <= N:
        result = shanks_transform(partial_sums[:N_use], order=8)
        diff = result - target_sum
        print(f"  Shanks(N={N_use:3d}, order=8): diff = {diff}")

# Also try simple Richardson: assume tail ~ C * q^N with q = 1/4
# Then S_N = L - C * q^N / (1-q) approximately
# So L ≈ (S_{2N} - q^N * S_N) / (1 - q^N)
print(f"\nRichardson (geometric q=1/4):")
q = mpf(1)/4
for N_use in [100, 200]:
    if 2*N_use <= N:
        S1 = partial_sums[N_use - 1]
        S2 = partial_sums[2*N_use - 1]
        L = (S2 - q**N_use * S1) / (1 - q**N_use)
        diff = L - target_sum
        print(f"  N={N_use}: Richardson estimate diff = {diff}")

print(f"\nDirect sum S_500: {partial_sums[-1]}")
print(f"Target:           {target_sum}")
