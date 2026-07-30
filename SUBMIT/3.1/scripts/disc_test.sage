## Decisive discreteness test for the alpha endpoint representation.
##
## Every elliptic element of a DISCRETE subgroup of PSL_2(R) has finite order,
## so its trace is zeta + zeta^{-1} for a root of unity zeta.  If the trace of
## the commutator word is elliptic but NOT of that form, the image is
## nondiscrete and no orbifold / Brooks-Goldman route can apply.

R.<a> = PolynomialRing(QQ)
falpha = (a^12 - 3*a^11 + 4*a^10 - 5*a^9 + 6*a^8 - 7*a^7
          + 7*a^6 - 7*a^5 + 6*a^4 - 5*a^3 + 4*a^2 - 3*a + 1)
print("f_alpha irreducible:", falpha.is_irreducible(), " degree:", falpha.degree())

K.<aa> = NumberField(falpha, embedding=RealField(200)("0.590989428670256"))
print("K degree:", K.degree(), " signature:", K.signature())

t  = aa^(-4)
nu = (t - 1) * (1 - aa) / (aa*t + 1)
Z  = nu^2 - (t + 1/t - 2)*nu + 2

RR200 = RealField(200)
print()
print("t  =", RR200(t))
print("nu =", RR200(nu))
print("Z  =", RR200(Z))
print("|Z| < 2 (elliptic)?", abs(RR200(Z)) < 2)

RX.<X> = PolynomialRing(QQ)
fZ = RX(Z.minpoly())
D = fZ.degree()
print()
print("minimal polynomial of Z, degree", D)
print(fZ)
print("Z is an algebraic integer:", Z.is_integral())

if not Z.is_integral():
    print()
    print(">>> Z is NOT an algebraic integer, so it cannot be zeta + zeta^{-1}.")
    print(">>> The element has infinite order  ==>  IMAGE IS NONDISCRETE.")
else:
    print()
    print("searching m with phi(m) <= 2D =", 2*D, ", m <= 8D^2 =", 8*D^2)
    hits = []
    for m in range(3, 8*D^2 + 1):
        if euler_phi(m) <= 2*D:
            C.<zeta> = CyclotomicField(m)
            psi = RX((zeta + zeta^-1).minpoly())
            g = gcd(fZ, psi)
            if g.degree() > 0:
                hits.append((m, g))
    print("hits:", hits)
    if not hits:
        print()
        print(">>> NO root-of-unity trace  ==>  infinite elliptic order")
        print(">>> IMAGE IS NONDISCRETE.")
    else:
        print()
        print(">>> Finite order candidate found; necessary but NOT sufficient")
        print(">>> for discreteness.  Needs Poincare-polygon verification.")
