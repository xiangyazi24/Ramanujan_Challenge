#!/usr/bin/env python3
"""P2.5: Test whether the Christoffel-Darboux polynomial K_n(X) gives
the CMF error via the Catalan moment functional.

C[R] = integral_0^1 (-log t)/(1+t^2) * R(t^2) dt
C[X^k] = (-1)^k * (G - sum_{j=0}^{k-1} (-1)^j/(2j+1)^2)
C[R] = G * R(-1) - A[R]   where A[R] is rational

K_n(X) = sum_{k=0}^n (2k+1) * D_k * P_k(1-2X)
K_n(-1) = sum_{k=0}^n (2k+1) * D_k^2

Test: is K_n(-1) proportional to Q_hat_n?
If so, C[K_n] gives G * Q_hat_n - P_hat_n.
"""
from fractions import Fraction
import sys

# Legendre polynomials evaluated at (1-2X), stored as polynomial coefficients
# P_k(1-2X) is a polynomial of degree k in X
def legendre_poly_shifted(k):
    """Return coefficients [c_0, c_1, ..., c_k] of P_k(1-2X) as polynomial in X."""
    if k == 0:
        return [Fraction(1)]
    if k == 1:
        return [Fraction(1), Fraction(-2)]

    # Use recurrence: (k+1)P_{k+1}(x) = (2k+1)x P_k(x) - k P_{k-1}(x)
    # with x = 1-2X
    prev = [Fraction(1)]  # P_0
    curr = [Fraction(1), Fraction(-2)]  # P_1 = 1-2X

    for n in range(1, k):
        # P_{n+1}(1-2X) = ((2n+1)(1-2X) P_n(1-2X) - n P_{n-1}(1-2X)) / (n+1)
        # Multiply curr by (2n+1)(1-2X) = (2n+1) - 2(2n+1)X
        new_len = len(curr) + 1
        mult = [Fraction(0)] * new_len
        for i, c in enumerate(curr):
            mult[i] += Fraction(2*n+1) * c
            mult[i+1] += Fraction(-2*(2*n+1)) * c

        # Subtract n * prev
        for i, c in enumerate(prev):
            mult[i] -= Fraction(n) * c

        # Divide by (n+1)
        for i in range(len(mult)):
            mult[i] /= Fraction(n+1)

        prev = curr
        curr = mult

    return curr

# Central Delannoy numbers D_k = P_k(3)
def delannoy(N):
    D = [Fraction(1)]
    if N == 0:
        return D
    D.append(Fraction(3))
    for k in range(1, N):
        D.append(Fraction(3*(2*k+1)) * D[k] / Fraction(k+1) - Fraction(k) * D[k-1] / Fraction(k+1))
    return D

# Catalan moment: C[X^k] = (-1)^k * (G - sum_{j=0}^{k-1} (-1)^j / (2j+1)^2)
# So C[X^k] = (-1)^k * G - (-1)^k * sum_{j<k} (-1)^j/(2j+1)^2
# = (-1)^k * G - rational
def catalan_moment_pair(k):
    """Return (q, p) such that C[X^k] = q*G - p."""
    q = Fraction((-1)**k)
    partial = sum(Fraction((-1)**j, (2*j+1)**2) for j in range(k))
    p = Fraction((-1)**k) * partial
    return q, p

# Build K_n(X) = sum_{k=0}^n (2k+1) D_k P_k(1-2X)
def christoffel_darboux(n, D_vals):
    """Return polynomial coefficients of K_n(X)."""
    max_deg = n
    coeffs = [Fraction(0)] * (max_deg + 1)
    for k in range(n + 1):
        weight = Fraction(2*k + 1) * D_vals[k]
        pk = legendre_poly_shifted(k)
        for i, c in enumerate(pk):
            coeffs[i] += weight * c
    return coeffs

# CMF M(n) matrix
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
    return [[Fraction(m11),Fraction(m12),Fraction(m13)],
            [Fraction(m21),Fraction(m22),Fraction(m23)],
            [Fraction(m31),Fraction(m32),Fraction(m33)]]

# Pochhammer normalization
def H(nn):
    if nn == 0:
        return Fraction(1)
    val = Fraction(1)
    for k in range(nn):
        val *= Fraction(-16) * Fraction(k+2)**2 * Fraction(k+3)**2 * Fraction(2*k+5,2) * Fraction(2*k+7,2)**2
    return val

# Compute Q_hat and P_hat from CMF
def compute_cmf_sequences(N):
    """Compute Q_hat_n, P_hat_n for n=0..N."""
    row_q = [Fraction(33750), Fraction(-36000), Fraction(9000)]
    row_p = [Fraction(30921), Fraction(-32972), Fraction(8240)]

    qhat = [row_q[0]]
    phat = [row_p[0]]

    for n in range(N):
        M = M_entries(n)
        new_q = [sum(row_q[k]*M[k][j] for k in range(3)) for j in range(3)]
        new_p = [sum(row_p[k]*M[k][j] for k in range(3)) for j in range(3)]
        row_q = new_q
        row_p = new_p
        h = H(n+1)
        qhat.append(Fraction(row_q[0], h))
        phat.append(Fraction(row_p[0], h))

    return qhat, phat

print("=== P2.5: Testing Christoffel-Darboux kernel ===", flush=True)

NMAX = 12

# Compute Delannoy numbers
D = delannoy(NMAX + 2)
print(f"Delannoy numbers D_0..D_{NMAX}: {[int(d) for d in D[:NMAX+1]]}")

# Compute K_n(-1) = sum (2k+1) D_k^2
print(f"\nK_n(-1) values:")
kn_at_minus1 = []
for n in range(NMAX + 1):
    val = sum(Fraction(2*k+1) * D[k]**2 for k in range(n+1))
    kn_at_minus1.append(val)
    print(f"  K_{n}(-1) = {val}")

# Compute CMF Q_hat_n
print(f"\nComputing CMF Q̂_n...", flush=True)
qhat, phat = compute_cmf_sequences(NMAX)

print(f"\nQ̂_n values:")
for n in range(min(8, len(qhat))):
    print(f"  Q̂_{n} = {qhat[n]}")

# Compare K_n(-1) / Q_hat_n
print(f"\n=== Comparing K_n(-1) with Q̂_n ===")
for n in range(min(len(kn_at_minus1), len(qhat))):
    if qhat[n] != 0:
        ratio = kn_at_minus1[n] / qhat[n]
        print(f"  n={n}: K_{n}(-1)/Q̂_{n} = {ratio} = {float(ratio):.10f}")

# Now compute Catalan moments C[K_n] = G * K_n(-1) - A[K_n]
print(f"\n=== Computing Catalan moments of K_n ===", flush=True)
for n in range(min(8, NMAX + 1)):
    print(f"  Building K_{n}...", end=" ", flush=True)
    kn_coeffs = christoffel_darboux(n, D)

    # C[K_n] = sum_k c_k * C[X^k] = G * K_n(-1) - A[K_n]
    q_total = Fraction(0)
    p_total = Fraction(0)
    for k, c in enumerate(kn_coeffs):
        qk, pk = catalan_moment_pair(k)
        q_total += c * qk
        p_total += c * pk

    # Verify q_total = K_n(-1)
    kn_eval = sum(c * Fraction(-1)**k for k, c in enumerate(kn_coeffs))
    assert q_total == kn_eval, f"Mismatch: q={q_total}, K(-1)={kn_eval}"

    print(f"C[K_{n}] = {q_total}·G - {p_total}")

    # Compare with ê_n = G*Q̂_n - P̂_n
    if n < len(qhat):
        print(f"          ê_{n}  = {qhat[n]}·G - {phat[n]}")
        if q_total == qhat[n]:
            print(f"          ✓ q matches!")
            if p_total == phat[n]:
                print(f"          ✓ p matches! EXACT MATCH!")
            else:
                print(f"          ✗ p differs by {p_total - phat[n]}")
        else:
            r = q_total / qhat[n] if qhat[n] != 0 else None
            print(f"          q ratio = {r} = {float(r):.10f}" if r else "          Q̂ = 0")

print("\nDone.")
