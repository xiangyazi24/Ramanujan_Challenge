#!/usr/bin/env python3
"""P2.5: Find ODE for F̂(z) using EULER operator θ = z d/dz.

Recurrence (order 3, degree 13) → Euler ODE (θ-order 13, z-degree 3).
Equation: Σ_{k=0}^{order} Σ_{j=0}^{D} β_{k,j} (m-j)^k Q̂_{m-j} = 0
"""

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

Q_hat_mod = [(Q0_mod[n] * modinv(H_mod[n])) % P if H_mod[n] != 0 else 0
             for n in range(NMAX + 1)]
print(f"Done: {len(Q_hat_mod)} terms")

def search_euler_ode(Q_mod, theta_order, z_deg, prime=P):
    """Search for Euler ODE: Σ_{k=0}^{order} q_k(z) θ^k F = 0.
    Coefficient of z^m: Σ_{k,j} β_{k,j} (m-j)^k Q̂_{m-j} = 0
    """
    N = len(Q_mod)
    D = z_deg
    order = theta_order
    num_unknowns = (order + 1) * (D + 1)

    m_min = D
    m_max = N - 1
    num_eqs = m_max - m_min + 1

    print(f"  θ-order={order}, z-deg={D}: {num_unknowns} unknowns, {num_eqs} available eqs")

    if num_eqs < num_unknowns + 3:
        print("  Insufficient data!")
        return None

    use_eqs = min(num_eqs, num_unknowns + 10)

    mat = []
    for eq_idx in range(use_eqs):
        m = m_min + eq_idx
        row = [0] * num_unknowns
        for k in range(order + 1):
            for j in range(D + 1):
                col = k * (D + 1) + j
                idx = m - j
                if 0 <= idx < N:
                    coeff = pow(m - j, k, prime) * Q_mod[idx] % prime
                    row[col] = coeff
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
        for j_col in range(num_unknowns):
            mat[r][j_col] = (mat[r][j_col] * inv) % prime
        for i in range(len(mat)):
            if i != r and mat[i][c] % prime != 0:
                factor = mat[i][c]
                for j_col in range(num_unknowns):
                    mat[i][j_col] = (mat[i][j_col] - factor * mat[r][j_col]) % prime
        r += 1

    rank = r
    null_dim = num_unknowns - rank
    print(f"  rank={rank}, null_dim={null_dim}")
    return null_dim

# The direct translation: θ-order 13, z-degree 3
print("\n" + "="*60)
print("Direct ODE from recurrence: θ-order 13, z-degree 3")
print("="*60)
nd = search_euler_ode(Q_hat_mod, 13, 3)

# Verify the direct translation works
# The ODE from recurrence is: c_3(m) Q̂_m + c_2(m-1) Q̂_{m-1} + c_1(m-2) Q̂_{m-2} + c_0(m-3) Q̂_{m-3} = 0
# In Euler form: Σ_{k=0}^{13} [β_{k,0} m^k Q̂_m + β_{k,1} (m-1)^k Q̂_{m-1} + β_{k,2} (m-2)^k Q̂_{m-2} + β_{k,3} (m-3)^k Q̂_{m-3}] = 0

# The original recurrence c_3(n) Q̂_{n+3} + ... maps to:
# c_3(m) Q̂_m + c_2(m-1) Q̂_{m-1} + c_1(m-2) Q̂_{m-2} + c_0(m-3) Q̂_{m-3} = 0
# So β_{k,j} encodes c_{3-j}(m-j) = Σ_k a_{3-j,k} (m-j)^k

print("\nVerifying direct translation numerically...")
from fractions import Fraction
c_coeffs = [
    [-42743162700, -206623731375, -448112471583, -579493151986, -500074412234,
     -304838513875, -135313819947, -44354837964, -10750665744, -1905022784,
     -240100240, -20397440, -1047552, -24576],  # c_0
    [8781630505200, 38850314624124, 78557994908508, 96136040496551, 79442239242197,
     46814452218572, 20241514501104, 6502490145168, 1552168938336, 271943188864,
     33995217088, 2871763456, 146952192, 3440640],  # c_1
    [-10566229124340, -43764612822972, -82725628159809, -94536939564882, -72904809920709,
     -40082159230086, -16169158004002, -4847446446296, -1080358338832, -176841798272,
     -20670362464, -1634185472, -78342144, -1720320],  # c_2
    [146862156672, 610678861056, 1158857071416, 1329423744980, 1029037642166,
     567735994679, 229759169143, 69074560420, 15430450432, 2530117664,
     296032016, 23408000, 1121280, 24576],  # c_3
]

def eval_poly_mod(coeffs, x, p):
    result = 0
    xk = 1
    for c in coeffs:
        result = (result + c * xk) % p
        xk = xk * x % p
    return result

# Check: c_3(m) Q̂_m + c_2(m-1) Q̂_{m-1} + c_1(m-2) Q̂_{m-2} + c_0(m-3) Q̂_{m-3} = 0
print("ODE recurrence verification:")
for m in [10, 50, 100, 200]:
    if m >= 3 and m < len(Q_hat_mod):
        val = 0
        for j in range(4):
            c_idx = 3 - j  # c_3 for j=0, c_2 for j=1, etc.
            c_val = eval_poly_mod(c_coeffs[c_idx], m - j, P)
            val = (val + c_val * Q_hat_mod[m - j]) % P
        print(f"  m={m}: residual mod P = {val}")

# Now search for MINIMAL Euler ODE (lower θ-order)
print("\n" + "="*60)
print("Searching for MINIMAL Euler ODE (θ-order < 13)")
print("="*60)
for theta_ord in range(1, 14):
    for z_d in range(3, 15):
        nd = search_euler_ode(Q_hat_mod, theta_ord, z_d)
        if nd and nd > 0:
            print(f"\n*** FOUND: θ-order={theta_ord}, z-deg={z_d}, null_dim={nd} ***")
            break
    else:
        continue
    break
