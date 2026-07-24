#!/usr/bin/env python3
"""
P2.5: Search for k-recurrence using modular arithmetic.
Compute α₁(k) mod p for large prime p, then search for
recurrence with polynomial coefficients.

With modular arithmetic, we can compute 200+ terms fast.
"""
import sys

P = (1 << 61) - 1  # Mersenne prime

def modinv(a, p=P):
    return pow(a % p, p - 2, p)

def M_entries_mod(n, p=P):
    n = n % p
    m11 = ((-2*n-5) * pow(n+3, 2, p) % p * (136*pow(n,4,p)+1424*pow(n,3,p)+5548*pow(n,2,p)+9551*n+6141)) % p
    m12 = (384*pow(n,6,p)+6384*pow(n,5,p)+44168*pow(n,4,p)+162698*pow(n,3,p)+336377*pow(n,2,p)+369933*n+169011) % p
    m13 = (-(480*pow(n,4,p)+4980*pow(n,3,p)+19210*pow(n,2,p)+32690*n+20730)) % p
    m21 = (pow(n+2, 2, p) * pow(n+3, 2, p) % p * (4*n+10) % p * (48*pow(n,3,p)+386*pow(n,2,p)+1017*n+879)) % p
    m22 = (pow(n+2, 2, p) * (-272*pow(n,5,p)-3848*pow(n,4,p)-21732*pow(n,3,p)-61184*pow(n,2,p)-85761*n-47808)) % p
    m23 = (pow(n+2, 2, p) * (320*pow(n,3,p)+2540*pow(n,2,p)+6610*n+5640)) % p
    m31 = ((-4*n-10) * pow(n+2, 2, p) % p * pow(n+3, 2, p) % p * (32*pow(n,4,p)+302*pow(n,3,p)+1037*pow(n,2,p)+1530*n+813)) % p
    m32 = (pow(n+2, 2, p) * (192*pow(n,6,p)+2984*pow(n,5,p)+19116*pow(n,4,p)+64452*pow(n,3,p)+120256*pow(n,2,p)+117279*n+46476)) % p
    m33 = (pow(n+2, 2, p) * (-16*pow(n,5,p)-408*pow(n,4,p)-2912*pow(n,3,p)-8884*pow(n,2,p)-12254*n-6240)) % p
    return [[m11 % p, m12 % p, m13 % p],
            [m21 % p, m22 % p, m23 % p],
            [m31 % p, m32 % p, m33 % p]]

def delta_H_mod(n, p=P):
    return (-2 * pow(n+2, 2, p) % p * pow(n+3, 2, p) % p * (2*n+5) % p * pow(2*n+7, 2, p)) % p

def B_mod(N, k, p=P):
    """B(N,k) = 2^k C(2k,k) C(N,k) C(N+k,k) mod p"""
    if k < 0 or k > N:
        return 0
    # Compute using factorials mod p
    result = pow(2, k, p)
    # C(2k,k) = (2k)! / (k!)^2
    # C(N,k) = N! / (k!(N-k)!)
    # C(N+k,k) = (N+k)! / (k! N!)
    # Total = 2^k (2k)! N! (N+k)! / (k!^4 (N-k)! N!)
    #       = 2^k (2k)! (N+k)! / (k!^4 (N-k)!)
    # But easier to compute each binomial separately

    # C(2k,k) mod p
    num = 1
    for i in range(1, k+1):
        num = num * (k + i) % p * modinv(i) % p
    result = result * num % p

    # C(N,k) mod p
    num = 1
    for i in range(k):
        num = num * (N - i) % p * modinv(i + 1) % p
    result = result * num % p

    # C(N+k,k) mod p
    num = 1
    for i in range(k):
        num = num * (N + k - i) % p * modinv(i + 1) % p
    result = result * num % p

    return result

KMAX = 200
print(f"Computing α₁(k) mod {P} for k=0..{KMAX}...", flush=True)

# Compute u₁(N) = e₁ · Π M_H · e₁ mod p
row = [1, 0, 0]
u1_vals = [row[0]]

for N in range(KMAX):
    M = M_entries_mod(N)
    d = delta_H_mod(N)
    dinv = modinv(d)
    MH = [[M[i][j] * dinv % P for j in range(3)] for i in range(3)]
    new_row = [sum(row[i] * MH[i][j] % P for i in range(3)) % P for j in range(3)]
    row = new_row
    u1_vals.append(row[0])
    if (N+1) % 50 == 0:
        print(f"  N={N+1} done", flush=True)

# Triangular inversion: α₁(k) from u₁(N) mod p
print(f"\nTriangular inversion mod p...", flush=True)
alpha1 = []
for K in range(KMAX + 1):
    rhs = u1_vals[K]
    for k in range(K):
        rhs = (rhs - alpha1[k] * B_mod(K, k)) % P
    bKK = B_mod(K, K)
    alpha1.append(rhs * modinv(bKK) % P)
    if (K+1) % 50 == 0:
        print(f"  k={K+1} done", flush=True)

print(f"\nFirst 5 α₁(k) mod p: {alpha1[:5]}")

# Search for recurrence: Σ_i c_i(k) α₁(k+i) = 0
# c_i(k) = Σ_j a_{ij} k^j, polynomials of degree ≤ d
def search_rec(seq, name, max_order=6, max_deg=25):
    L = len(seq)
    print(f"\nSearching k-recurrence for {name} (L={L})...", flush=True)

    for order in range(2, max_order + 1):
        for deg in range(1, max_deg + 1):
            n_unknowns = (order + 1) * (deg + 1)
            n_eqs = L - order
            if n_eqs < n_unknowns + 2:
                continue

            # Build matrix mod P
            mat = []
            for k_val in range(min(n_eqs, n_unknowns + 5)):
                row = []
                for i in range(order + 1):
                    s = seq[k_val + i]
                    for j in range(deg + 1):
                        row.append(pow(k_val, j, P) * s % P)
                mat.append([x % P for x in row])

            # Gaussian elimination mod P
            nrows = len(mat)
            ncols = n_unknowns
            pivot_rows = 0
            for col in range(ncols):
                found = False
                for r in range(pivot_rows, nrows):
                    if mat[r][col] % P != 0:
                        mat[pivot_rows], mat[r] = mat[r], mat[pivot_rows]
                        inv = modinv(mat[pivot_rows][col])
                        for r2 in range(nrows):
                            if r2 != pivot_rows and mat[r2][col] % P != 0:
                                factor = mat[r2][col] * inv % P
                                for c2 in range(ncols):
                                    mat[r2][c2] = (mat[r2][c2] - factor * mat[pivot_rows][c2]) % P
                        pivot_rows += 1
                        found = True
                        break

            rank = pivot_rows
            nullity = ncols - rank
            if nullity > 0:
                print(f"  order={order}, deg={deg}: unknowns={n_unknowns}, "
                      f"rank={rank}, nullity={nullity}")
                if nullity == 1:
                    print(f"  *** FOUND UNIQUE: order {order}, deg {deg} ***")
                    return order, deg

    print(f"  No recurrence found")
    return None, None

search_rec(alpha1, "α₁(k)", max_order=5, max_deg=22)

print("\nDone.")
