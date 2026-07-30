## Avenue (a): explicit denominator bound for the torsion of K_3^ind(F).
##
## Zickert (J. reine angew. Math. 704 (2015) 21-54, Thm 1.1): the extended Bloch
## group is isomorphic to K_3^ind(F).  Merkurjev-Suslin:
##      |K_3^ind(F)_tors| = w_2(F)
## where w_2(F) = prod_p p^{nu_p},  nu_p = max{ nu : Gal(F(mu_{p^nu})/F) is
## killed by 2 }, equivalently the largest p^nu with
## [F(zeta_{p^nu} + zeta_{p^nu}^{-1}) : F] = 1, i.e. the maximal real cyclotomic
## subfield contained in F.
##
## Practical computation: for each prime power p^nu whose maximal real
## cyclotomic field Q(zeta_{p^nu})^+ has degree dividing [F:Q], test whether
## that real cyclotomic field embeds in F.

R.<z> = PolynomialRing(QQ)
falpha = (z^12 - 3*z^11 + 4*z^10 - 5*z^9 + 6*z^8 - 7*z^7
          + 7*z^6 - 7*z^5 + 6*z^4 - 5*z^3 + 4*z^2 - 3*z + 1)
fbeta  = (z^16 - 7*z^15 + 22*z^14 - 48*z^13 + 87*z^12 - 133*z^11 + 178*z^10
          - 211*z^9 + 223*z^8 - 211*z^7 + 178*z^6 - 133*z^5 + 87*z^4
          - 48*z^3 + 22*z^2 - 7*z + 1)

def w2_bound(f, name):
    K.<t> = NumberField(f)
    d = K.degree()
    disc = K.discriminant()
    print("="*64)
    print(name, ": degree", d)
    print("  discriminant =", disc)
    print("  factored     =", factor(disc))
    print("  signature    =", K.signature())
    # w_2 always contains the factor 24 for any number field (w_2(Q) = 24)
    w2 = 1
    details = []
    # candidate prime powers: the maximal real cyclotomic subfield
    # Q(zeta_m)^+ has degree phi(m)/2 (m>2); need phi(m)/2 | d
    for p in prime_range(2, 200):
        nu_max = 0
        for nu in range(1, 8):
            m = p^nu
            deg = euler_phi(m)//2 if m > 2 else 1
            if deg == 0: deg = 1
            if d % deg != 0:
                break
            # necessary: p ramifies in F unless deg == 1
            if deg > 1 and disc % p != 0:
                break
            # test embedding of Q(zeta_m)^+ into K
            if deg == 1:
                nu_max = nu
                continue
            g = (CyclotomicField(m).gen() + CyclotomicField(m).gen()^-1).minpoly()
            if len(g.roots(K)) > 0:
                nu_max = nu
            else:
                break
        if nu_max:
            w2 *= p^nu_max
            details.append("%d^%d" % (p, nu_max))
    print("  w_2(F) =", w2, "=", " * ".join(details) if details else "1")
    return w2

wa = w2_bound(falpha, "F_alpha = Q(a)")
wb = w2_bound(fbeta,  "F_beta  = Q(b)")

print()
print("="*64)
print("Denominator bound for Re[Delta R]/pi^2")
print("="*64)
Q = lcm(wa, wb)
print("  the difference of two torsion classes has order dividing lcm:")
print("  Q = lcm(%d, %d) = %d" % (wa, wb, Q))
print()
err = RealField(60)("1.6332113e-301")
print("  measured error |value - (-4/85)| =", err)
print("  reconstruction is unique when err < 1/(2 Q^2):")
print("     1/(2 Q^2) =", RealField(60)(1)/(2*Q^2))
print("     err       =", err)
ok = err < RealField(60)(1)/(2*Q^2)
print("     CONDITION SATISFIED:", ok)
print()
print("  margin: Q could be as large as", RealField(30)(sqrt(1/(2*err))))
print("  actual Q =", Q, " -> margin factor", RealField(20)(sqrt(1/(2*err))/Q))
