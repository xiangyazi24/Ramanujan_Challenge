#!/usr/bin/env python3
"""
MOP for mixed Beukers kernel — pure Python/mpmath version.
"""
import mpmath
mpmath.mp.dps = 200

def psi_tail(j, m):
    """ψ_j(m) = ζ(j) - H_{m-1}^{(j)}"""
    if j == 2:
        return mpmath.zeta(2) - sum(mpmath.mpf(1)/mpmath.mpf(k)**2 for k in range(1, m))
    elif j == 3:
        return mpmath.zeta(3) - sum(mpmath.mpf(1)/mpmath.mpf(k)**3 for k in range(1, m))

def compute_mop(n, n1, n2):
    """Compute monic MOP Q_n(t) = t^n + c_{n-1}t^{n-1} + ... + c_0."""
    assert n1 + n2 == n and n1 >= 1 and n2 >= 1

    A = mpmath.matrix(n, n)
    b = mpmath.matrix(n, 1)

    for row in range(n1):
        s = row + 1
        for k in range(n):
            A[row, k] = psi_tail(2, k + s)
        b[row, 0] = -psi_tail(2, n + s)

    for row in range(n2):
        s = row + 1
        for k in range(n):
            A[n1 + row, k] = psi_tail(3, k + s)
        b[n1 + row, 0] = -psi_tail(3, n + s)

    c_vec = mpmath.lu_solve(A, b)
    c = [c_vec[k, 0] for k in range(n)] + [mpmath.mpf(1)]
    Qn1 = sum(c)
    return c, Qn1

# P2.7 coefficients
def A_p27(n): return 1024 * (2*n+5)**4 * (2*n+7)**3 * (2*n+9)**3 * (946*n**2+6407*n+10860)
def B_p27(n): return 128 * (2*n+7)**3 * (2*n+9)**3 * (104060*n**6 + 1745370*n**5 + 12145238*n**4 + 44886481*n**3 + 92943995*n**2 + 102256019*n + 46709052)
def C_p27(n): return 16 * (n+3)**4 * (2*n+9)**3 * (3784*n**5 + 57792*n**4 + 351019*n**3 + 1059230*n**2 + 1587211*n + 944620)
def D_p27(n): return (n+3)**4 * (n+4)**6 * (946*n**2 + 4515*n + 5399)

# P2.7 reference
print("=== P2.7 reference ===")
q = [mpmath.mpf(0)] * 20
q[0] = mpmath.mpf('-215040420000')
q[1] = mpmath.mpf('-167282265043404') / mpmath.mpf('905')
q[2] = mpmath.mpf('-964185327658080') / mpmath.mpf('6071')
for n in range(2, 19):
    q[n+1] = (mpmath.mpf(B_p27(n))/A_p27(n) * q[n]
              - mpmath.mpf(C_p27(n-1))/A_p27(n-1) * q[n-1]
              + mpmath.mpf(D_p27(n-2))/A_p27(n-2) * q[n-2])
for n in range(8):
    print(f"  q[{n+1}]/q[{n}] = {float(q[n+1]/q[n]):.12f}")

# Main computation
print("\n" + "="*70)
print("MOP with polygamma moments")
print("="*70)

max_n = 14

for step_name, n1_fn in [
    ("balanced", lambda n: (n+1)//2),
    ("ζ(2)-heavy (n-1,1)", lambda n: n-1),
    ("ζ(3)-heavy (1,n-1)", lambda n: 1),
    ("2:1", lambda n: (2*n+2)//3),
]:
    print(f"\n--- {step_name} ---")
    results = []
    for n in range(2, max_n + 1):
        n1 = n1_fn(n)
        n2 = n - n1
        if n2 < 1 or n1 < 1:
            continue
        try:
            c, Qn1 = compute_mop(n, n1, n2)
            results.append((n, Qn1, c))
            print(f"  n={n:2d} (n1={n1},n2={n2}): Q_n(1) = {float(Qn1):+.12e}")
        except Exception as e:
            print(f"  n={n:2d}: FAILED: {str(e)[:80]}")

    if len(results) >= 2:
        print("  Growth ratios:")
        for i in range(len(results) - 1):
            na, va, _ = results[i]
            nb, vb, _ = results[i+1]
            if abs(va) > 1e-50 and nb == na + 1:
                r = float(vb / va)
                print(f"    n={na}→{nb}: {r:+.10f}  (target: ~0.859)")

# Error computation
print("\n" + "="*70)
print("Error decay")
print("="*70)
L = mpmath.zeta(2) + mpmath.zeta(3)
print(f"L = ζ(2)+ζ(3) = {float(L):.20f}")

for step_name, n1_fn in [
    ("balanced", lambda n: (n+1)//2),
    ("ζ(2)-heavy", lambda n: n-1),
    ("ζ(3)-heavy", lambda n: 1),
]:
    print(f"\n--- {step_name} ---")
    errors = []
    for n in range(2, min(max_n + 1, 12)):
        n1 = n1_fn(n)
        n2 = n - n1
        if n2 < 1 or n1 < 1:
            continue
        try:
            c, Qn1 = compute_mop(n, n1, n2)
            p_n = mpmath.mpf(0)
            for k in range(n+1):
                Hk2 = sum(mpmath.mpf(1)/mpmath.mpf(j)**2 for j in range(1, k+1))
                Hk3 = sum(mpmath.mpf(1)/mpmath.mpf(j)**3 for j in range(1, k+1))
                p_n += c[k] * (Hk2 + Hk3)
            e_n = Qn1 * L - p_n
            errors.append((n, e_n, Qn1))
            rel = float(abs(e_n / Qn1)) if abs(Qn1) > 1e-100 else 0
            print(f"  n={n:2d}: |e_n| = {float(abs(e_n)):.6e}, |e_n/Q_n(1)| = {rel:.6e}")
        except Exception as e:
            print(f"  n={n:2d}: {str(e)[:60]}")

    if len(errors) >= 2:
        print("  Error decay ratios:")
        for i in range(len(errors) - 1):
            na, ea, _ = errors[i]
            nb, eb, _ = errors[i+1]
            if abs(ea) > 1e-200 and nb == na + 1:
                r = float(abs(eb) / abs(ea))
                print(f"    n={na}→{nb}: {r:.10f}  (target: ~0.067)")

# P2.7 recurrence check
print("\n" + "="*70)
print("P2.7 recurrence check for balanced MOP Q_n(1)")
print("="*70)
vals = []
for n in range(2, max_n + 1):
    n1 = (n+1)//2; n2 = n - n1
    if n2 < 1: continue
    try:
        c, Qn1 = compute_mop(n, n1, n2)
        vals.append((n, Qn1))
    except:
        pass

for i in range(3, len(vals)):
    nc = vals[i][0]
    if not (vals[i-1][0] == nc-1 and vals[i-2][0] == nc-2 and vals[i-3][0] == nc-3):
        continue
    u = [vals[i-3+j][1] for j in range(4)]
    nn = nc - 3
    pred = (mpmath.mpf(B_p27(nn+2))/A_p27(nn+2) * u[2]
            - mpmath.mpf(C_p27(nn+1))/A_p27(nn+1) * u[1]
            + mpmath.mpf(D_p27(nn))/A_p27(nn) * u[0])
    res = u[3] - pred
    rel = float(abs(res / u[3])) if abs(u[3]) > 1e-100 else float(abs(res))
    tag = " *** MATCH ***" if rel < 1e-10 else ""
    print(f"  n={nc}: rel residual = {rel:.6e}{tag}")
