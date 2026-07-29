#!/usr/bin/env python3
"""P2.5: Find the minimal ODE for F̂(z) = Σ Q̂_n z^n.

Search for ODE of order 3 with polynomial z-coefficients of degree D:
Σ_{k=0}^3 p_k(z) F^{(k)}(z) = 0
where p_k(z) = Σ_{j=0}^D α_{k,j} z^j.

The coefficient of z^m gives:
Σ_{k=0}^3 Σ_{j=0}^D α_{k,j} · falling(m-j+k, k) · Q̂_{m-j+k} = 0
"""
from fractions import Fraction
from math import gcd as igcd

def M_exact_int(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def H_frac(n):
    result = Fraction(1)
    for k in range(n):
        result *= Fraction(-16)
        result *= Fraction(k+2)**2
        result *= Fraction(k+3)**2
        result *= Fraction(2*k+5, 2)
        result *= Fraction(2*k+7, 2)**2
    return result

# Compute Q̂_n for n = 0..80
print("Computing Q̂_n for n = 0..80...")
q = [33750, -36000, 9000]
Q0_int = [33750]
for N in range(80):
    M = M_exact_int(N)
    new_q = [0, 0, 0]
    for j in range(3):
        for k in range(3):
            new_q[j] += q[k] * M[k][j]
    q = new_q
    Q0_int.append(q[0])

Q_hat = []
for n in range(81):
    Hn = H_frac(n)
    Q_hat.append(Fraction(Q0_int[n]) / Hn if Hn != Fraction(0) else Fraction(Q0_int[n]))
print(f"  Done: {len(Q_hat)} terms")

def falling(n, k):
    """Falling factorial n^{(k)} = n(n-1)...(n-k+1)"""
    result = 1
    for i in range(k):
        result *= (n - i)
    return result

def search_ode(Q_hat, order, max_D, verbose=True):
    """Search for minimal ODE of given order with z-degree D."""
    N = len(Q_hat)

    for D in range(2, max_D + 1):
        num_unknowns = (order + 1) * (D + 1)

        # Equation for each m: sum_{k=0}^{order} sum_{j=0}^D alpha_{k,j} * falling(m-j+k,k) * Q_hat[m-j+k] = 0
        # Valid range: m-j+k >= 0 for all j,k => m >= D (worst: j=D, k=0)
        # and m-j+k < N for all j,k => m+order < N (worst: j=0, k=order)
        m_min = D
        m_max = N - order - 1

        num_eqs = m_max - m_min + 1
        if num_eqs < num_unknowns + 3:
            print(f"  D={D}: not enough equations ({num_eqs} for {num_unknowns} unknowns)")
            continue

        # Use modular arithmetic first for speed: work mod a large prime
        P = 2**61 - 1  # Mersenne prime

        # Convert Q_hat to mod P
        Q_mod = []
        for q in Q_hat:
            num = q.numerator % P
            den_inv = pow(q.denominator % P, P - 2, P)
            Q_mod.append((num * den_inv) % P)

        # Build matrix mod P
        use_eqs = min(num_eqs, num_unknowns + 5)
        mat = []
        for eq_idx in range(use_eqs):
            m = m_min + eq_idx
            row = [0] * num_unknowns
            for k in range(order + 1):
                for j in range(D + 1):
                    col = k * (D + 1) + j
                    idx = m - j + k
                    if 0 <= idx < N:
                        ff = falling(m - j + k, k)
                        row[col] = (ff * Q_mod[idx]) % P
            mat.append(row)

        # Gaussian elimination mod P
        r = 0
        pivot_cols = []
        for c in range(num_unknowns):
            found = -1
            for i in range(r, len(mat)):
                if mat[i][c] % P != 0:
                    found = i
                    break
            if found < 0:
                continue
            mat[r], mat[found] = mat[found], mat[r]
            pivot_cols.append(c)
            inv = pow(mat[r][c], P - 2, P)
            for j in range(num_unknowns):
                mat[r][j] = (mat[r][j] * inv) % P
            for i in range(len(mat)):
                if i != r and mat[i][c] % P != 0:
                    factor = mat[i][c]
                    for j in range(num_unknowns):
                        mat[i][j] = (mat[i][j] - factor * mat[r][j]) % P
            r += 1

        rank = r
        null_dim = num_unknowns - rank
        if verbose:
            print(f"  D={D}: {num_unknowns} unknowns, {use_eqs} eqs, rank={rank}, null_dim={null_dim}")

        if null_dim > 0:
            # Verify with exact arithmetic
            print(f"  *** Potential ODE found at D={D}! Verifying with exact arithmetic...")
            mat_exact = []
            for eq_idx in range(min(num_eqs, num_unknowns + 10)):
                m = m_min + eq_idx
                row = [Fraction(0)] * (num_unknowns + 1)
                for k in range(order + 1):
                    for j in range(D + 1):
                        col = k * (D + 1) + j
                        idx = m - j + k
                        if 0 <= idx < N:
                            ff = falling(m - j + k, k)
                            row[col] = Fraction(ff) * Q_hat[idx]
                mat_exact.append(row)

            # Gaussian elimination
            m_e = len(mat_exact)
            nc = num_unknowns
            pivot_cols_e = []
            r_e = 0
            for c in range(nc):
                found = -1
                for i in range(r_e, m_e):
                    if mat_exact[i][c] != Fraction(0):
                        found = i
                        break
                if found < 0:
                    continue
                mat_exact[r_e], mat_exact[found] = mat_exact[found], mat_exact[r_e]
                pivot_cols_e.append(c)
                pv = mat_exact[r_e][c]
                for j in range(nc + 1):
                    mat_exact[r_e][j] /= pv
                for i in range(m_e):
                    if i != r_e and mat_exact[i][c] != Fraction(0):
                        f = mat_exact[i][c]
                        for j in range(nc + 1):
                            mat_exact[i][j] -= f * mat_exact[r_e][j]
                r_e += 1

            null_dim_exact = nc - r_e
            print(f"    Exact rank={r_e}, null_dim={null_dim_exact}")

            if null_dim_exact > 0:
                # Extract null vector
                free_cols = [c for c in range(nc) if c not in pivot_cols_e]
                x = [Fraction(0)] * nc
                fc = free_cols[0]
                x[fc] = Fraction(1)
                for pr_idx in range(len(pivot_cols_e) - 1, -1, -1):
                    pc = pivot_cols_e[pr_idx]
                    s = sum(mat_exact[pr_idx][j] * x[j] for j in range(pc + 1, nc))
                    x[pc] = -s

                # Verify on ALL data
                max_res = Fraction(0)
                for eq_idx in range(num_eqs):
                    m = m_min + eq_idx
                    val = Fraction(0)
                    for k in range(order + 1):
                        for j in range(D + 1):
                            col = k * (D + 1) + j
                            idx = m - j + k
                            if 0 <= idx < N:
                                ff = falling(m - j + k, k)
                                val += x[col] * Fraction(ff) * Q_hat[idx]
                    max_res = max(max_res, abs(val))

                if max_res == Fraction(0):
                    print(f"\n*** CONFIRMED: Minimal ODE of order {order}, z-degree {D} ***")

                    # Clear denominators and extract polynomials
                    lcm_den = 1
                    for c in x:
                        if c != 0:
                            lcm_den = lcm_den * c.denominator // igcd(lcm_den, c.denominator)
                    x_int = [int(c * lcm_den) for c in x]
                    g = 0
                    for c in x_int:
                        g = igcd(g, abs(c))
                    if g > 0:
                        x_int = [c // g for c in x_int]

                    for k in range(order + 1):
                        coeffs = x_int[k * (D + 1):(k + 1) * (D + 1)]
                        # Find actual degree
                        deg = 0
                        for j in range(len(coeffs) - 1, -1, -1):
                            if coeffs[j] != 0:
                                deg = j
                                break
                        print(f"  p_{k}(z) = degree {deg}: {coeffs[:deg+1]}")

                    return D, x_int
                else:
                    print(f"    Residual nonzero: {float(max_res):.2e} — false positive from mod arithmetic")

    return None, None

print("\nSearching for minimal ODE of order 3...")
D_min, coeffs = search_ode(Q_hat, 3, 20)

if D_min is None:
    print("\nNo order-3 ODE found with z-degree ≤ 20.")
    print("Trying order 2 (in case module has rank 2)...")
    D_min2, coeffs2 = search_ode(Q_hat, 2, 25)
