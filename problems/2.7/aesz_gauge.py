"""
Test: 256^n · q_n vs A_n = C(2n,n) · b_n

Since AESZ roots = 256 × P2.7 roots, the sequence g_n = 256^n · q_n
should have the same dominant exponential as A_n.

Also test: is there a polynomial gauge R(n) such that
g_n = R(n) · A_n  or  g_n = Σ R_i(n) · A_{n+i} ?
"""
from fractions import Fraction as Q
from math import comb, factorial

def pochhammer(a, n):
    result = Q(1)
    for i in range(n):
        result *= Q(a) + i
    return result

def A_c(n):
    n = Q(n)
    return (1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860))

def B_c(n):
    n = Q(n)
    return (128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052))

def C_c(n):
    n = Q(n)
    return (16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620))

def D_c(n):
    n = Q(n)
    return ((n+3)**4*(n+4)**6*(946*n*n+4515*n+5399))

def compute_b(n):
    return sum(comb(n,k)**2 * comb(n+k,n) * comb(n+2*k,n) for k in range(n+1))

N = 40

# Exact b_n and A_n = C(2n,n) · b_n
b = [Q(compute_b(n)) for n in range(N)]
A = [Q(comb(2*n, n)) * b[n] for n in range(N)]

# Exact q_n
q = [Q(-215040420000), Q(-167282265043404, 905), Q(-964185327658080, 6071)]
for i in range(3, N):
    n = i - 1
    new_q = Q(B_c(n), A_c(n)) * q[-1] - Q(C_c(n-1), A_c(n-1)) * q[-2] + Q(D_c(n-2), A_c(n-2)) * q[-3]
    q.append(new_q)

# p_n
p = [Q(-612218384750), Q(-9525021973931919, 18100), Q(-29561828382772029, 65380)]
for i in range(3, N):
    n = i - 1
    new_p = Q(B_c(n), A_c(n)) * p[-1] - Q(C_c(n-1), A_c(n-1)) * p[-2] + Q(D_c(n-2), A_c(n-2)) * p[-3]
    p.append(new_p)

# g_n = 256^n · q_n
g = [Q(256**n) * q[n] for n in range(N)]

print("=== g_n = 256^n · q_n vs A_n = C(2n,n) · b_n ===")
print("\ng_n / A_n:")
for n in range(25):
    if A[n] != 0:
        ratio = g[n] / A[n]
        print(f"  n={n}: {float(ratio):.15e}")

# Successive ratios of g_n/A_n
print("\n(g_{n+1}/A_{n+1}) / (g_n/A_n):")
for n in range(1, 20):
    if A[n] != 0 and A[n+1] != 0 and g[n] != 0:
        r0 = g[n] / A[n]
        r1 = g[n+1] / A[n+1]
        print(f"  n={n}: {float(r1/r0):.15f}")

# Try: g_n = R(n) · A_n with polynomial R(n)
# f_n = g_n/A_n should be polynomial if simple gauge exists
print("\n=== Testing polynomial gauge g_n = R(n) · A_n ===")
f = [g[n] / A[n] if A[n] != 0 else Q(0) for n in range(N)]

# Check differences: Δf, Δ²f, etc.
df = [f[n+1] - f[n] for n in range(N-1)]
d2f = [df[n+1] - df[n] for n in range(N-2)]
d3f = [d2f[n+1] - d2f[n] for n in range(N-3)]

print("\nΔ(g/A):")
for n in range(10):
    print(f"  n={n}: {float(df[n]):.10e}")

print("\nΔ²(g/A):")
for n in range(10):
    print(f"  n={n}: {float(d2f[n]):.10e}")

# g_n/A_n is clearly not polynomial. Try: does g_n satisfy the AESZ recurrence?
# AESZ #209: c₃(n)A_{n+3} + c₂(n)A_{n+2} + c₁(n)A_{n+1} + c₀(n)A_n = 0
def c3_aesz(n):
    n = Q(n)
    return (n+3)**4 * (946*n**2 + 3053*n + 2475)

def c2_aesz(n):
    n = Q(n)
    return -2*(104060*n**6 + 1376430*n**5 + 7533488*n**4 + 21825779*n**3 + 35282480*n**2 + 30160606*n + 10648707)

def c1_aesz(n):
    n = Q(n)
    return 4*(2*n+3)*(3784*n**5 + 36808*n**4 + 141179*n**3 + 267255*n**2 + 250336*n + 93060)

def c0_aesz(n):
    n = Q(n)
    return -4*(n+1)**2*(2*n+1)*(2*n+3)*(946*n**2 + 4945*n + 6474)

print("\n=== AESZ recurrence verification on A_n ===")
for n in range(5):
    res = c3_aesz(n)*A[n+3] + c2_aesz(n)*A[n+2] + c1_aesz(n)*A[n+1] + c0_aesz(n)*A[n]
    print(f"  n={n}: {res}")

print("\n=== Does g_n = 256^n · q_n satisfy AESZ? ===")
for n in range(5):
    res = c3_aesz(n)*g[n+3] + c2_aesz(n)*g[n+2] + c1_aesz(n)*g[n+1] + c0_aesz(n)*g[n]
    print(f"  n={n}: {float(res):.6e}")

# Also try the GAUGED AESZ: the 256^n factor changes the recurrence
# If A_n satisfies AESZ, then G_n = 256^n · A_n satisfies a gauged version
# with roots scaled by 256.
# The gauged version: replace ρ_A by 256·ρ_A in the char poly,
# which gives roots 256·220 ≈ 56320, etc.
# But q_n satisfies P2.7 with roots ρ₀ ≈ 0.859, not 56320.
# So 256^n · q_n has dominant rate 256^n · 0.859^n = (256·0.859)^n ≈ 220^n = ρ_{A,0}^n.
# This is the SAME rate as A_n.

# Search for Ore intertwiner: g_n = u₀(n)A_n + u₁(n)A_{n+1} + u₂(n)A_{n+2}
print("\n=== Searching for g_n = u₀(n)A_n + u₁(n)A_{n+1} + u₂(n)A_{n+2} ===")
for deg in range(8):
    n_unknowns = 3*(deg+1)
    n_eqs = min(n_unknowns + 5, N - 2)

    mat = []
    rhs = []
    for n_idx in range(n_eqs):
        row = []
        for i in range(3):
            for j in range(deg+1):
                row.append(Q(n_idx)**j * A[n_idx+i])
        mat.append(row)
        rhs.append(g[n_idx])

    # Gaussian elimination
    aug = [mat[i][:] + [rhs[i]] for i in range(n_eqs)]
    m = n_eqs
    nc = n_unknowns

    pivot_cols = []
    r = 0
    for c in range(nc):
        pivot = None
        for i in range(r, m):
            if aug[i][c] != Q(0):
                pivot = i
                break
        if pivot is None:
            continue
        aug[r], aug[pivot] = aug[pivot], aug[r]
        pivot_cols.append(c)
        for i in range(m):
            if i != r and aug[i][c] != Q(0):
                factor = aug[i][c] / aug[r][c]
                for j in range(nc+1):
                    aug[i][j] -= factor * aug[r][j]
        r += 1

    consistent = True
    for i in range(r, m):
        if aug[i][nc] != Q(0):
            consistent = False
            break

    if consistent and r == nc:
        sol = [Q(0)] * nc
        for idx, c in enumerate(pivot_cols):
            sol[c] = aug[idx][nc] / aug[idx][c]

        ok = True
        for n_idx in range(n_eqs, min(N-2, n_eqs+5)):
            pred = Q(0)
            for i in range(3):
                u_i = Q(0)
                for j in range(deg+1):
                    u_i += sol[i*(deg+1)+j] * Q(n_idx)**j
                pred += u_i * A[n_idx+i]
            if pred != g[n_idx]:
                ok = False
                break

        if ok:
            print(f"\n*** FOUND solution at degree {deg}! ***")
            for i in range(3):
                coeffs = [sol[i*(deg+1)+j] for j in range(deg+1)]
                print(f"  u_{i}(n) = {coeffs}")
        else:
            print(f"  deg={deg}: fails verification")
    elif consistent:
        print(f"  deg={deg}: underdetermined")
    else:
        print(f"  deg={deg}: inconsistent")
