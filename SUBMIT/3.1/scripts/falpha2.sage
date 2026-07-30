R.<M,L> = PolynomialRing(QQ, 2)
A = (L^5
     + L^4*( M^14 - M^12 + 3*M^4 + 4*M^2 - 2 )
     + L^3*( -2*M^18 + 5*M^16 + M^14 - 4*M^12 + 6*M^8 + 5*M^6 + 2*M^4 - 4*M^2 + 1 )
     + L^2*( M^22 - 4*M^20 + 2*M^18 + 5*M^16 + 6*M^14 - 4*M^10 + M^8 + 5*M^6 - 2*M^4 )
     + L  *( -2*M^22 + 4*M^20 + 3*M^18 - M^10 + M^8 )
     + M^22)

T.<a> = PolynomialRing(QQ)
Aalpha = T(A.subs({M: a^2, L: a}))

target = RR("0.590989428670256")
print("=== factors of A(a^2,a), locating the endpoint by real roots ===")
falpha = None
for f, e in Aalpha.factor():
    rts = [r for r in f.roots(RR, multiplicities=False)]
    near = [r for r in rts if abs(r - target) < 1e-9]
    mark = "  <== ENDPOINT" if near else ""
    print("degree %2d  mult %d  real roots: %s%s"
          % (f.degree(), e, [RR(r).n(20) for r in rts][:6], mark))
    if near:
        falpha = f

print()
if falpha is None:
    print("NOT FOUND")
else:
    c = falpha.coefficients(sparse=False)
    print("f_alpha, degree", falpha.degree())
    print(falpha)
    print()
    print("coefficients ascending:", c)
    print("palindromic:", c == list(reversed(c)))
    print("irreducible:", falpha.is_irreducible())
    # certified root
    K.<aa> = NumberField(falpha, embedding=RealField(200)(target))
    print("number field degree:", K.degree())
    print("signature (r1, r2):", K.signature())
    print("a  =", RealField(60)(aa))
    print("a^2=", RealField(60)(aa^2), "   (alpha)")
    # save for the next stage
    with open("/tmp/falpha.txt", "w") as fh:
        fh.write(str(falpha) + "\n")
        fh.write(",".join(str(x) for x in c) + "\n")
    print()
    print("saved to /tmp/falpha.txt")
