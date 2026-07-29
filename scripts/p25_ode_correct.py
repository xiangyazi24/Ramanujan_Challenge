#!/usr/bin/env python3
"""P2.5: Find minimal ODE — CORRECT formulation.

Recurrence: order 3, degree 13 → ODE: order ≤ 13, z-degree ≤ 3.
Search with z-degree fixed at D (≤ 3) and vary ODE order.
"""
import sys

def M_entries(n):
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

def delta_int(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

P = (1 << 61) - 1

def modinv(a, p=P):
    return pow(a % p, p - 2, p)

NMAX = 300
print(f"Computing Q̂_n mod P for n = 0..{NMAX}...")

q = [33750 % P, (-36000) % P, 9000 % P]
Q0_mod = [33750 % P]
H_mod = [1]
h_val = 1

for N in range(NMAX):
    M = M_entries(N)
    new_q = [0, 0, 0]
    for j in range(3):
        for k in range(3):
            new_q[j] = (new_q[j] + q[k] * M[k][j]) % P
    q = new_q
    Q0_mod.append(q[0])
    d = delta_int(N) % P
    h_val = (h_val * d) % P
    H_mod.append(h_val)

Q_hat_mod = []
for n in range(NMAX + 1):
    if H_mod[n] == 0:
        Q_hat_mod.append(0)
    else:
        Q_hat_mod.append((Q0_mod[n] * modinv(H_mod[n])) % P)

print(f"Done: {len(Q_hat_mod)} terms")

def falling(n, k):
    result = 1
    for i in range(k):
        result *= (n - i)
    return result

def search_ode_fixed_zdeg(Q_mod, z_deg, order_range, prime=P):
    """Search for ODE with fixed z-degree D and varying order."""
    N = len(Q_mod)
    D = z_deg

    for order in order_range:
        # ODE: sum_{k=0}^{order} p_k(z) F^{(k)}(z) = 0
        # p_k(z) = sum_{j=0}^D alpha_{k,j} z^j
        # Unknowns: (order+1) * (D+1)
        num_unknowns = (order + 1) * (D + 1)

        # Coefficient of z^m: sum_{k,j} alpha_{k,j} * falling(m-j+k, k) * Q[m-j+k]
        m_min = D
        m_max = N - order - 1
        num_eqs = m_max - m_min + 1

        if num_eqs < num_unknowns + 3:
            print(f"  z_deg={D}, order={order:2d}: insufficient data ({num_eqs} eqs for {num_unknowns} unknowns)")
            continue

        use_eqs = min(num_eqs, num_unknowns + 10)

        mat = []
        for eq_idx in range(use_eqs):
            m = m_min + eq_idx
            row = [0] * num_unknowns
            for k in range(order + 1):
                for j in range(D + 1):
                    col = k * (D + 1) + j
                    idx = m - j + k
                    if 0 <= idx < N:
                        ff = falling(m - j + k, k) % prime
                        row[col] = (ff * Q_mod[idx]) % prime
            mat.append(row)

        # Gaussian elimination mod prime
        r = 0
        pivot_cols = []
        for c in range(num_unknowns):
            found = -1
            for i in range(r, len(mat)):
                if mat[i][c] % prime != 0:
                    found = i
                    break
            if found < 0:
                continue
            mat[r], mat[found] = mat[found], mat[r]
            pivot_cols.append(c)
            inv = pow(mat[r][c], prime - 2, prime)
            for j in range(num_unknowns):
                mat[r][j] = (mat[r][j] * inv) % prime
            for i in range(len(mat)):
                if i != r and mat[i][c] % prime != 0:
                    factor = mat[i][c]
                    for j in range(num_unknowns):
                        mat[i][j] = (mat[i][j] - factor * mat[r][j]) % prime
            r += 1

        rank = r
        null_dim = num_unknowns - rank
        print(f"  z_deg={D}, order={order:2d}: {num_unknowns:4d} unknowns, {use_eqs:4d} eqs, rank={rank:4d}, null_dim={null_dim}")

        if null_dim > 0:
            return order, null_dim
    return None, 0

# Search with z-degree = 3 (from recurrence order 3)
print("\n" + "="*60)
print("z-degree = 3, varying ODE order")
print("="*60)
for zd in [3]:
    o, nd = search_ode_fixed_zdeg(Q_hat_mod, zd, range(1, 14))
    if o:
        print(f"\n*** ODE found: z-degree={zd}, order={o}, null_dim={nd} ***")
        break

# Also try z-degree = 4 (might have fewer apparent singularities)
print("\n" + "="*60)
print("z-degree = 4, varying ODE order")
print("="*60)
for zd in [4]:
    o, nd = search_ode_fixed_zdeg(Q_hat_mod, zd, range(1, 14))
    if o:
        print(f"\n*** ODE found: z-degree={zd}, order={o}, null_dim={nd} ***")
        break

# And z-degree = 5
print("\n" + "="*60)
print("z-degree = 5, varying ODE order")
print("="*60)
for zd in [5]:
    o, nd = search_ode_fixed_zdeg(Q_hat_mod, zd, range(1, 14))
    if o:
        print(f"\n*** ODE found: z-degree={zd}, order={o}, null_dim={nd} ***")
        break
