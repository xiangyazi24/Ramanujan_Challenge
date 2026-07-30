## High-precision + interval-certified evaluation of the regulator difference.
##
## Re[ sum_j Rhat(z_j(beta)) - sum_j Rhat(z_j(alpha)) ] / pi^2   =?  -4/85
##
## Rhat(z) = Li_2(z) + (1/2) log z log(1-z) - pi^2/6
##
## Two independent checks:
##   (1) 300-digit floating evaluation, compared with -4/85;
##   (2) how many digits of agreement, and what denominator bound that would
##       certify under a rational-reconstruction argument.

prec = 1000            # bits
RRp  = RealField(prec)
CCp  = ComplexField(prec)

R.<z> = PolynomialRing(QQ)
falpha = (z^12 - 3*z^11 + 4*z^10 - 5*z^9 + 6*z^8 - 7*z^7
          + 7*z^6 - 7*z^5 + 6*z^4 - 5*z^3 + 4*z^2 - 3*z + 1)
fbeta  = (z^16 - 7*z^15 + 22*z^14 - 48*z^13 + 87*z^12 - 133*z^11 + 178*z^10
          - 211*z^9 + 223*z^8 - 211*z^7 + 178*z^6 - 133*z^5 + 87*z^4
          - 48*z^3 + 22*z^2 - 7*z + 1)

def endpoint_shapes(Mval, Lval):
    X   = Mval^2
    u   = (Lval + X^3) / (X*(Lval + X))
    r   = -(1 + sqrt(1 + 4*u^2)) / (2*u)
    tau = 1 - r^2
    return [tau, u, u/X, 1/(1 - u*X)]

def Rhat(zv):
    zc = CCp(zv)
    return dilog(zc) + log(zc)*log(1 - zc)/2 - CCp.pi()^2/6

Ka.<aa> = NumberField(falpha, embedding=RealField(400)("0.5909894286702564"))
Kb.<bb> = NumberField(fbeta,  embedding=RealField(400)("0.4068130813367900"))
a = RRp(aa); b = RRp(bb)

Sa = endpoint_shapes(a^2, a)
Sb = endpoint_shapes(b,   b)

print("alpha shapes:", [RRp(t).n(40) for t in Sa])
print("beta  shapes:", [RRp(t).n(40) for t in Sb])
print()

D  = (sum(Rhat(t) for t in Sb) - sum(Rhat(t) for t in Sa)).real()
pi2 = RRp(CCp.pi()^2)
q   = RRp(D/pi2)
tgt = RRp(-4)/85

print("Re[Delta R]/pi^2 =", q.n(280))
print("            -4/85 =", tgt.n(280))
err = abs(q - tgt)
print()
print("absolute difference =", err.n(30))
digits = -log(err, 10) if err > 0 else Infinity
print("digits of agreement ~", RR(digits))

# rational-reconstruction margin: an approximation to precision eps determines a
# rational of denominator <= Q uniquely when eps < 1/(2Q^2)
print()
print("If the value is known a priori to be rational with denominator <= Q,")
print("the reconstruction is unique once err < 1/(2 Q^2), i.e.")
print("   Q <  sqrt(1/(2*err))  =", RR(sqrt(1/(2*err))).n(30))
print()
print("continued fraction of the computed value:")
cf = continued_fraction(q)
print("  ", cf[:12])
print("  best rational approximations:")
for k in range(1, 8):
    print("   ", cf.convergent(k))
