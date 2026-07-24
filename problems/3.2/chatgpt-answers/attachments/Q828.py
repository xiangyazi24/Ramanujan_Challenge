# Run with: sage q828_exact_check.sage
from sage.all import *

HMAX = 300
BMAX = 8
SHOW_H = [20, 30, 50, 75, 100, 150, 200, 250, 300]

K.<rt2> = QuadraticField(2)
lam = K(17) + K(12)*rt2
gam = lam - 1
s = lam/gam
eta = 1/(2*gam)

# Apéry numbers, exactly.
A = [ZZ(0)]*(HMAX + 1)
A[0] = 1
A[1] = 5
for m in range(1, HMAX):
    num = (34*m^3 + 51*m^2 + 27*m + 5)*A[m] - m^3*A[m-1]
    den = (m+1)^3
    assert num % den == 0
    A[m+1] = num // den

# Forward differences c_k = Delta^k A_0, exactly.
c = []
for k in range(HMAX + 1):
    c.append(sum((-1)^(k-j)*binomial(k,j)*A[j] for j in range(k+1)))
    assert c[-1] > 0


def M(k, ell):
    if ell > k:
        return ZZ(0)
    return binomial(k, ell)*sum(
        binomial(ell,t)*c[k-ell+t] for t in range(ell+1)
    )


def exact_ratio(H, b):
    a = H-b
    n = 3*H + 1

    if b == 0:
        q = [QQ(1)]
    else:
        R = matrix(QQ, [
            [M(k,ell) for ell in range(b)]
            for k in range(a+1, H+1)
        ])
        rhs = vector(QQ, [-M(k,b) for k in range(a+1, H+1)])
        assert R.det() != 0
        q = list(R.solve_right(rhs)) + [QQ(1)]

    d = [
        sum(q[ell]*M(k,ell) for ell in range(b+1))
        for k in range(a+1)
    ]

    # Verify the annihilated high rows exactly.
    for k in range(a+1, H+1):
        assert sum(q[ell]*M(k,ell) for ell in range(b+1)) == 0

    Phat = sum(d[k]*binomial(n,k) for k in range(a+1))
    pred = ((-1)^b * s^b * K(c[a]) * K(binomial(n,a))
            * (1-eta)^(-b-1))
    rho = K(Phat)/pred
    return rho, Phat, q


print("selected exact ratios")
print("H,b,coeff_1,coeff_sqrt2,decimal,error")
for H in SHOW_H:
    for b in range(min(BMAX,H)+1):
        rho, Phat, q = exact_ratio(H,b)
        # rho = rho[0] + rho[1]*sqrt(2), both rational.
        dec = rho.n(80)
        err = abs(dec - 1)
        print(H, b, rho[0], rho[1], dec, err, sep=",")

print("\nworst deviation for each fixed strip, 8 <= H <= 300")
for b in range(BMAX+1):
    worst = RealField(100)(0)
    where = None
    for H in range(max(8,b), HMAX+1):
        rho, _, _ = exact_ratio(H,b)
        err = abs(rho.n(100) - 1)
        if err > worst:
            worst = err
            where = H
    print("b =", b, "max_error =", worst, "at H =", where)
