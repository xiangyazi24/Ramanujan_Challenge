import sympy as sp
from sympy import Poly, GF
X = sp.Symbol('X'); Y = sp.Symbol('Y'); Z = sp.Symbol('Z')
def A_(Y_): return 34*Y_**3 - 51*Y_**2 + 27*Y_ - 5

def N_mod_p(h, p):
    Pm = [Poly(0, X, modulus=p), Poly(1, X, modulus=p)]
    for k in range(2, h+1):
        Pk = Poly(A_(X+k), X, modulus=p)*Pm[-1] - Poly((X+k-1)**6, X, modulus=p)*Pm[-2]
        Pm.append(Pk)
    return Pm[h]

def analyze(p, h):
    Nh = N_mod_p(h, p)
    roots = Nh.ground_roots()
    nroots = sum(roots.values())
    # fold: c = -(h+1)/2 mod p
    inv2 = pow(2, -1, p); c = (-(h+1)*inv2) % p
    Nc = Poly(Nh.as_expr().subs(X, c + Y), Y, modulus=p)
    coeffs = Nc.all_coeffs()
    deg = Nc.degree()
    odd_part_zero = all(int(coeffs[deg-i]) % p == 0 for i in range(1, deg+1, 2)) if h%2==1 else None
    # for odd h: N(c+Y) = H(Y^2)
    if h % 2 == 1:
        Hc = [coeffs[deg-i] for i in range(deg, -1, -2)]  # coeffs of Y^deg ... down step2
        H = Poly(Hc, Z, modulus=p)
        fl = H.factor_list()[1]
        pat = sorted(f.degree() for f,_ in fl)
        zroots = [r for r,m in H.ground_roots().items() for _ in range(m)]
        qr = [z for z in zroots if z != 0 and pow(int(z) % p, (p-1)//2, p) == 1]
        return dict(R=nroots, deg=Nh.degree(), fold_ok=True, Hdeg=H.degree(), pattern=pat[:12], nz=len(zroots), nQR=len(qr))
    else:
        return dict(R=nroots, deg=Nh.degree(), even=True)

import time
for (p, hs) in [(3001, [13, 11, 15, 127]), (10007, [127, 439, 125, 129])]:
    for h in hs:
        t=time.time()
        try:
            d = analyze(p, h)
            print(p, h, d, f"{time.time()-t:.1f}s", flush=True)
        except Exception as e:
            print(p, h, "ERR", e, flush=True)
