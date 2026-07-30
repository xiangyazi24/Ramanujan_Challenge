## A4 repair attempt (a): is 1 + 4u^2 a SQUARE in the endpoint field?
## If yes, the shapes live in F itself and the discriminants used are correct.
## Also recompute w_2 WITH the leading factor 2 (Merkurjev-Suslin), and check
## the mandatory 24 | w_2.

R.<z> = PolynomialRing(QQ)
falpha = (z^12 - 3*z^11 + 4*z^10 - 5*z^9 + 6*z^8 - 7*z^7
          + 7*z^6 - 7*z^5 + 6*z^4 - 5*z^3 + 4*z^2 - 3*z + 1)
fbeta  = (z^16 - 7*z^15 + 22*z^14 - 48*z^13 + 87*z^12 - 133*z^11 + 178*z^10
          - 211*z^9 + 223*z^8 - 211*z^7 + 178*z^6 - 133*z^5 + 87*z^4
          - 48*z^3 + 22*z^2 - 7*z + 1)

def check(f, name, alpha_chart):
    print("="*66)
    print(name)
    K.<t> = NumberField(f)
    if alpha_chart:
        M = t^2; L = t
    else:
        M = t;   L = t
    X = M^2
    u = (L + X^3)/(X*(L + X))
    d = 1 + 4*u^2
    print("  u   =", u)
    print("  1+4u^2 in K :", d)
    issq = d.is_square()
    print("  IS 1+4u^2 A SQUARE IN K ?  ->", issq)
    if issq:
        s = d.sqrt()
        print("     sqrt =", s)
        print("     check s^2 == d :", s^2 == d)
    else:
        # the shape field is the quadratic extension K(sqrt d)
        print("  => shapes live in the quadratic extension K(sqrt(1+4u^2))")
        S.<Y> = PolynomialRing(K)
        L2.<w> = K.extension(Y^2 - d)
        print("     shape field degree over Q :", K.degree()*2)
    return K, u, d

Ka, ua, da = check(falpha, "F_alpha = Q(a),  alpha chart (M=a^2, L=a)", True)
Kb, ub, db = check(fbeta,  "F_beta  = Q(b),  beta chart  (M=L=b)",     False)

print()
print("="*66)
print("w_2 WITH the Merkurjev-Suslin leading factor 2")
print("="*66)

def w2(K, name):
    d = K.degree()
    disc = K.discriminant()
    w = 2                      # the leading factor
    parts = ["2"]
    for p in prime_range(2, 300):
        numax = 0
        for nu in range(1, 8):
            m = p^nu
            deg = max(euler_phi(m)//2, 1)
            if d % deg != 0: break
            if deg > 1 and disc % p != 0: break
            if deg == 1:
                numax = nu; continue
            g = (CyclotomicField(m).gen() + CyclotomicField(m).gen()^-1).minpoly()
            if len(g.roots(K)) > 0: numax = nu
            else: break
        if numax:
            w *= p^numax; parts.append("%d^%d" % (p, numax))
    print("  %s : degree %d, disc %s" % (name, d, factor(disc)))
    print("     w_2 = %d = %s" % (w, " * ".join(parts)))
    print("     24 | w_2 ?  %s   (w_2/24 = %s)" % (w % 24 == 0, w/24 if w % 24 == 0 else "-"))
    return w

wa = w2(Ka, "F_alpha")
wb = w2(Kb, "F_beta")
Q = lcm(wa, wb)
print()
print("  lcm =", Q)
print("  85 | lcm ?", Q % 85 == 0, "  quotient", Q//85 if Q % 85 == 0 else "-")
err = RealField(60)("1.1e-301")
thr = RealField(60)(1)/(2*Q^2)
print("  threshold 1/(2Q^2) =", thr.n(20))
print("  err                =", err)
print("  RECONSTRUCTION STILL VALID:", err < thr)
