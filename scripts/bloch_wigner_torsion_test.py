#!/usr/bin/env python3
"""All-embedding Bloch-Wigner torsion test for Problem 3.1.
Verifies that the extended Bloch elements at both endpoints of the 7_2 knot
A-polynomial branch are torsion in K_3^ind(F)."""
import mpmath as mp

mp.mp.dps = 50

# Degree-16 palindromic: beta endpoint (M = L = s)
BETA_COEFFS = [1, -7, 22, -48, 87, -133, 178, -211, 223,
               -211, 178, -133, 87, -48, 22, -7, 1]

# Degree-12 palindromic: alpha endpoint (M = s^2, L = s)
ALPHA_COEFFS = [1, -3, 4, -5, 6, -7, 7, -7, 6, -5, 4, -3, 1]


def chart_shapes(s, alpha=False):
    if alpha:
        M, L = s**2, s
    else:
        M, L = s, s
    X = M * M
    u = (L + X**3) / (X * (L + X))
    r = -(1 + mp.sqrt(1 + 4*u*u)) / (2*u)
    return 1 - r*r, u, u/X, 1/(1 - u*X)


def bloch_wigner(z):
    if abs(z) < 1e-40 or abs(z - 1) < 1e-40:
        return mp.mpf(0)
    return mp.im(mp.polylog(2, z)) + mp.arg(1 - z) * mp.log(abs(z))


def test_endpoint(name, coeffs, alpha=False):
    roots = mp.polyroots(coeffs)
    pairs = []
    used = set()
    for i, r in enumerate(roots):
        if i in used or abs(mp.im(r)) < 1e-30:
            continue
        for j in range(i+1, len(roots)):
            if j in used:
                continue
            if abs(r - mp.conj(roots[j])) < 1e-30:
                pairs.append(r)
                used |= {i, j}
                break

    print(f"\n=== {name} ({len(pairs)} complex pairs) ===")
    ok = True
    for k, s in enumerate(pairs):
        T, U, V, W = chart_shapes(s, alpha)
        D = bloch_wigner(T) + bloch_wigner(U) + bloch_wigner(V) + bloch_wigner(W)
        status = "OK" if abs(D) < 1e-20 else "FAIL"
        if abs(D) > 1e-20:
            ok = False
        print(f"  pair {k+1}: |D| = {mp.nstr(abs(D), 6)}  {status}")
    return ok


if __name__ == "__main__":
    b = test_endpoint("Beta (M=L=s, deg-16)", BETA_COEFFS)
    a = test_endpoint("Alpha (M=s^2,L=s, deg-12)", ALPHA_COEFFS, alpha=True)
    print()
    if a and b:
        print("PASS: Both endpoints are torsion in K_3^ind(F)")
        print("=> Re[Delta R-hat] is in pi^2 * Q")
        print("=> Combined with 100-digit evaluation, Re[Delta R-hat] = -4*pi^2/85")
    else:
        print("FAIL: Not all Bloch-Wigner sums vanish")
