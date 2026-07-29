#!/usr/bin/env python3
"""P2.5: Find minimal ODE using modular arithmetic with many terms.

Compute Q̂_n mod P for n = 0..N, then search for ODE of order 3.
H_n is an integer: H_{n+1}/H_n = delta(n) = -2(n+2)^2(n+3)^2(2n+5)(2n+7)^2.
"""
import sys

def M_entries(n):
    """Return M(n) entries as integers."""
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
    """delta(n) = -2*(n+2)^2*(n+3)^2*(2n+5)*(2n+7)^2"""
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

P = (1 << 61) - 1  # Mersenne prime 2^61 - 1

def modinv(a, p=P):
    return pow(a % p, p - 2, p)

# Compute Q0_int[n] mod P and H_n mod P incrementally
NMAX = 300
print(f"Computing Q̂_n mod P for n = 0..{NMAX}...")

q = [33750 % P, (-36000) % P, 9000 % P]
Q0_mod = [33750 % P]
H_mod = [1]  # H_0 = 1

h_val = 1  # H_n mod P, incrementally

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

    if N % 50 == 49:
        print(f"  N={N+1}")

# Compute Q̂_n mod P = Q0_mod[n] * modinv(H_mod[n])
Q_hat_mod = []
for n in range(NMAX + 1):
    if H_mod[n] == 0:
        Q_hat_mod.append(0)
    else:
        Q_hat_mod.append((Q0_mod[n] * modinv(H_mod[n])) % P)

print(f"  Done: {len(Q_hat_mod)} terms")

# Verify first few values match exact computation
from fractions import Fraction
def H_frac(n):
    result = Fraction(1)
    for k in range(n):
        result *= Fraction(-16)
        result *= Fraction(k+2)**2
        result *= Fraction(k+3)**2
        result *= Fraction(2*k+5, 2)
        result *= Fraction(2*k+7, 2)**2
    return result

q_exact = [33750, -36000, 9000]
print("\nVerification (first 5 terms):")
for n in range(5):
    if n == 0:
        q0 = 33750
    else:
        M = M_entries(n-1)
        new_q = [0, 0, 0]
        for j in range(3):
            for k in range(3):
                new_q[j] += q_exact[k] * M[k][j]
        q_exact = new_q
        q0 = q_exact[0]

    Hn = H_frac(n)
    q_hat_exact = Fraction(q0) / Hn if Hn != 0 else Fraction(0)
    q_hat_mod_val = (q0 * modinv(int(Hn)) if Hn.denominator == 1 else
                     (q0 * int(Hn.denominator) * modinv(int(Hn.numerator)))) % P
    q_hat_exact_mod = (q_hat_exact.numerator * modinv(q_hat_exact.denominator)) % P
    match = Q_hat_mod[n] == q_hat_exact_mod
    print(f"  n={n}: mod_computed={Q_hat_mod[n]}, mod_exact={q_hat_exact_mod}, match={match}")

def falling(n, k):
    result = 1
    for i in range(k):
        result *= (n - i)
    return result

# Search for ODE of order 3 with z-degree D
def search_ode_mod(Q_mod, order, D_range, prime=P):
    """Search for minimal ODE using modular Gaussian elimination."""
    N = len(Q_mod)

    for D in D_range:
        num_unknowns = (order + 1) * (D + 1)
        m_min = D
        m_max = N - order - 1
        num_eqs = m_max - m_min + 1

        if num_eqs < num_unknowns + 3:
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
        print(f"  D={D:3d}: {num_unknowns:4d} unknowns, {use_eqs:4d} eqs, rank={rank:4d}, null_dim={null_dim}")

        if null_dim > 0:
            return D, null_dim
    return None, 0

print("\n" + "="*60)
print("Searching for minimal ODE of order 3...")
print("="*60)
D_found, nd = search_ode_mod(Q_hat_mod, 3, range(2, 60))

if D_found:
    print(f"\n*** ODE found at z-degree D = {D_found}, null_dim = {nd} ***")
else:
    print("\nNo order-3 ODE found up to z-degree 59.")
    print("\nTrying order 1 (in case the recurrence is reducible)...")
    D1, nd1 = search_ode_mod(Q_hat_mod, 1, range(2, 80))
    if D1:
        print(f"\n*** Order-1 ODE found at z-degree D = {D1}, null_dim = {nd1} ***")
