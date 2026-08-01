# Pointwise test: A_p(phi(x)) == H_p(x)^2 for x in F_p \ {-1}, phi(x)=x(1-8x)/(1+x)
from math import comb
for p in [13, 29, 37, 101]:
    b = [1, 5]
    for n in range(1, p):
        b.append(((34*n**3+51*n**2+27*n+5)*b[n] - n**3*b[n-1]) * pow((n+1)**3, p-2, p) % p)
    A = b[:p]                       # A_p coeffs
    f = [sum(comb(n,k)**3 for k in range(n+1)) % p for n in range(p)]  # Franel
    ok, bad = 0, []
    for x in range(p):
        if (1+x) % p == 0: continue
        phi = x*(1-8*x) % p * pow(1+x, p-2, p) % p
        lhs = sum(A[n]*pow(phi, n, p) for n in range(p)) % p
        H = sum(f[n]*pow(x, n, p) for n in range(p)) % p
        if lhs == H*H % p: ok += 1
        else: bad.append((x, phi, lhs, H*H % p))
    print(f"p={p}: pointwise OK {ok}/{p-1}", ("BAD:"+str(bad[:3]) if bad else ""))
