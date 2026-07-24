"""
Search for Ore intertwiner between P2.7 and twisted Zudilin operator.

Q5175 §7: Define ĥ_n = 64^{-n} · (5/2)_n / n! · b_n
where b_n = Σ C(n,k)²C(n+k,n)C(n+2k,n).

Search for U = u₀(n) + u₁(n)S + u₂(n)S² such that
q_n = u₀(n)ĥ_n + u₁(n)ĥ_{n+1} + u₂(n)ĥ_{n+2}
"""
from fractions import Fraction as Q
from math import comb, factorial

def pochhammer(a, n):
    """(a)_n = a(a+1)...(a+n-1)"""
    result = Q(1)
    for i in range(n):
        result *= Q(a) + i
    return result

# P2.7 coefficients
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

# Compute b_n
def compute_b(n):
    return sum(comb(n,k)**2 * comb(n+k,n) * comb(n+2*k,n) for k in range(n+1))

# Compute sequences
N = 35
b = [Q(compute_b(n)) for n in range(N)]
print("b[0..5]:", [int(b[i]) for i in range(6)])

# Twisted: ĥ_n = 64^{-n} · (5/2)_n / n! · b_n
h_hat = []
for n in range(N):
    val = Q(1, 64**n) * pochhammer(Q(5,2), n) / Q(factorial(n)) * b[n]
    h_hat.append(val)
print("ĥ[0..5]:", [float(h_hat[i]) for i in range(6)])

# P2.7 forward: q_n
q = [Q(-215040420000), Q(-167282265043404, 905), Q(-964185327658080, 6071)]
for i in range(3, N):
    n = i - 1
    new_q = Q(B_c(n), A_c(n)) * q[-1] - Q(C_c(n-1), A_c(n-1)) * q[-2] + Q(D_c(n-2), A_c(n-2)) * q[-3]
    q.append(new_q)

# P2.7 forward: p_n
p = [Q(-612218384750), Q(-9525021973931919, 18100), Q(-29561828382772029, 65380)]
for i in range(3, N):
    n = i - 1
    new_p = Q(B_c(n), A_c(n)) * p[-1] - Q(C_c(n-1), A_c(n-1)) * p[-2] + Q(D_c(n-2), A_c(n-2)) * p[-3]
    p.append(new_p)

print("\nq[0..3]:", [float(q[i]) for i in range(4)])
print("ĥ[0..3]:", [float(h_hat[i]) for i in range(4)])

# Test: q_n = u₀(n)ĥ_n + u₁(n)ĥ_{n+1} + u₂(n)ĥ_{n+2}
# For polynomial u_i(n) = Σ a_{i,j} n^j

# Also test UNTWISTED: f_n = q_n · 64^n · n! / (5/2)_n = u₀(n)b_n + u₁(n)b_{n+1} + u₂(n)b_{n+2}
f = []
for n in range(N):
    val = q[n] * Q(64**n) * Q(factorial(n)) / pochhammer(Q(5,2), n)
    f.append(val)

print("\nf_n = q_n * 64^n * n! / (5/2)_n:")
print("f[0..5]:", [float(f[i]) for i in range(6)])

# Check ratios f_{n+1}/f_n
print("\nf_{n+1}/f_n ratios:")
for n in range(min(10, N-1)):
    if f[n] != 0:
        print(f"  n={n}: {float(f[n+1]/f[n]):.10f}")

# Search for polynomial u_i(n) of degree d
# q_n = u₀(n)ĥ_n + u₁(n)ĥ_{n+1} + u₂(n)ĥ_{n+2}
# Equivalently: f_n = u₀(n)b_n + u₁(n)b_{n+1} + u₂(n)b_{n+2}

print("\n=== Searching for polynomial intertwiner ===")
for deg in range(8):
    # u_i(n) = Σ_{j=0}^{deg} c_{i,j} n^j
    # Total unknowns: 3*(deg+1)
    num_unknowns = 3*(deg+1)
    # Equations: f_n = u₀(n)b_n + u₁(n)b_{n+1} + u₂(n)b_{n+2} for n=0,...,num_eqs-1
    num_eqs = min(num_unknowns + 5, N - 2)

    # Build matrix A where A[n, :] corresponds to equation at n
    # unknowns: c_{0,0}, c_{0,1}, ..., c_{0,deg}, c_{1,0}, ..., c_{1,deg}, c_{2,0}, ..., c_{2,deg}
    mat = []
    rhs = []
    for n in range(num_eqs):
        row = []
        for i in range(3):
            for j in range(deg+1):
                row.append(Q(n)**j * b[n+i])
        # For h_hat version: row would use h_hat[n+i]
        mat.append(row)
        rhs.append(f[n])

    # Solve via Gaussian elimination over Q
    # Augmented matrix
    aug = [mat[i] + [rhs[i]] for i in range(num_eqs)]
    m = num_eqs
    nc = num_unknowns

    # Forward elimination
    pivot_cols = []
    r = 0
    for c in range(nc):
        # Find pivot
        pivot = None
        for i in range(r, m):
            if aug[i][c] != Q(0):
                pivot = i
                break
        if pivot is None:
            continue
        aug[r], aug[pivot] = aug[pivot], aug[r]
        pivot_cols.append(c)
        # Eliminate
        for i in range(m):
            if i != r and aug[i][c] != Q(0):
                factor = aug[i][c] / aug[r][c]
                for j in range(nc+1):
                    aug[i][j] -= factor * aug[r][j]
        r += 1

    # Check consistency
    consistent = True
    for i in range(r, m):
        if aug[i][nc] != Q(0):
            consistent = False
            break

    if consistent and r == nc:
        # Unique solution
        sol = [Q(0)] * nc
        for idx, c in enumerate(pivot_cols):
            sol[c] = aug[idx][nc] / aug[idx][c]

        # Verify on remaining points
        ok = True
        for n in range(num_eqs, min(N-2, num_eqs+5)):
            pred = Q(0)
            for i in range(3):
                u_i = Q(0)
                for j in range(deg+1):
                    u_i += sol[i*(deg+1)+j] * Q(n)**j
                pred += u_i * b[n+i]
            if pred != f[n]:
                ok = False
                break

        if ok:
            print(f"\n*** FOUND solution at degree {deg}! ***")
            for i in range(3):
                coeffs = [sol[i*(deg+1)+j] for j in range(deg+1)]
                print(f"  u_{i}(n) = {' + '.join(f'({c})*n^{j}' for j,c in enumerate(coeffs) if c != 0)}")
            # Verify all
            max_check = min(N-2, 30)
            all_ok = True
            for n in range(max_check):
                pred = Q(0)
                for i in range(3):
                    u_i = Q(0)
                    for j in range(deg+1):
                        u_i += sol[i*(deg+1)+j] * Q(n)**j
                    pred += u_i * b[n+i]
                if pred != f[n]:
                    all_ok = False
                    print(f"  FAILED at n={n}")
                    break
            if all_ok:
                print(f"  Verified for n=0..{max_check-1}")
        else:
            print(f"  deg={deg}: solution found but fails verification beyond training set")
    elif consistent:
        print(f"  deg={deg}: underdetermined (rank {r} < {nc})")
    else:
        print(f"  deg={deg}: inconsistent (no solution)")

# Also try: u_i(n) = P_i(n) / common_denom(n) with denominator being a product of small linear factors
print("\n=== Trying rational intertwiner with denominator (n+1)(n+2)(n+3) ===")
for deg_num in range(6):
    num_unknowns = 3*(deg_num+1)
    num_eqs = min(num_unknowns + 5, N - 5)

    mat = []
    rhs = []
    for n in range(num_eqs):
        denom = Q((n+1)*(n+2)*(n+3))
        row = []
        for i in range(3):
            for j in range(deg_num+1):
                row.append(Q(n)**j * b[n+i] / denom)
        mat.append(row)
        rhs.append(f[n])

    # Gaussian elimination (same as above)
    aug = [mat[ii] + [rhs[ii]] for ii in range(num_eqs)]
    m_rows = num_eqs
    nc2 = num_unknowns

    pivot_cols2 = []
    r2 = 0
    for c in range(nc2):
        pivot = None
        for ii in range(r2, m_rows):
            if aug[ii][c] != Q(0):
                pivot = ii
                break
        if pivot is None:
            continue
        aug[r2], aug[pivot] = aug[pivot], aug[r2]
        pivot_cols2.append(c)
        for ii in range(m_rows):
            if ii != r2 and aug[ii][c] != Q(0):
                factor = aug[ii][c] / aug[r2][c]
                for jj in range(nc2+1):
                    aug[ii][jj] -= factor * aug[r2][jj]
        r2 += 1

    consistent2 = True
    for ii in range(r2, m_rows):
        if aug[ii][nc2] != Q(0):
            consistent2 = False
            break

    if consistent2 and r2 == nc2:
        sol2 = [Q(0)] * nc2
        for idx, c in enumerate(pivot_cols2):
            sol2[c] = aug[idx][nc2] / aug[idx][c]

        ok2 = True
        for n in range(num_eqs, min(N-5, num_eqs+3)):
            denom = Q((n+1)*(n+2)*(n+3))
            pred = Q(0)
            for i in range(3):
                u_i = Q(0)
                for j in range(deg_num+1):
                    u_i += sol2[i*(deg_num+1)+j] * Q(n)**j
                pred += u_i * b[n+i] / denom
            if pred != f[n]:
                ok2 = False
                break

        if ok2:
            print(f"\n*** FOUND rational solution, numerator deg {deg_num}, denom (n+1)(n+2)(n+3)! ***")
        else:
            print(f"  num_deg={deg_num}: fails verification")
    elif consistent2:
        print(f"  num_deg={deg_num}: underdetermined")
    else:
        print(f"  num_deg={deg_num}: inconsistent")
