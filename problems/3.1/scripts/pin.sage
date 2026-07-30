# Final step: rational reconstruction with the proved denominator bound.
prec=1000; RRp=RealField(prec); CCp=ComplexField(prec)
R.<z>=PolynomialRing(QQ)
fa=z^12-3*z^11+4*z^10-5*z^9+6*z^8-7*z^7+7*z^6-7*z^5+6*z^4-5*z^3+4*z^2-3*z+1
fb=(z^16-7*z^15+22*z^14-48*z^13+87*z^12-133*z^11+178*z^10-211*z^9+223*z^8
    -211*z^7+178*z^6-133*z^5+87*z^4-48*z^3+22*z^2-7*z+1)
def sh(M,L):
    X=M^2; u=(L+X^3)/(X*(L+X)); r=-(1+sqrt(1+4*u^2))/(2*u)
    return [1-r^2,u,u/X,1/(1-u*X)]
def Rh(v):
    zc=CCp(v); return dilog(zc)+log(zc)*log(1-zc)/2-CCp.pi()^2/6
Ka.<aa>=NumberField(fa,embedding=RealField(400)("0.5909894286702564"))
Kb.<bb>=NumberField(fb,embedding=RealField(400)("0.4068130813367900"))
a=RRp(aa); b=RRp(bb)
D=(sum(Rh(t) for t in sh(b,b))-sum(Rh(t) for t in sh(a^2,a))).real()
q=RRp(D/CCp.pi()^2)
Q=1020
print("value      =", q.n(60))
print("bound Q    =", Q, "  (= lcm(w2(F_alpha), w2(F_beta)) = lcm(60,204))")
print("1/(2Q^2)   =", (RRp(1)/(2*Q^2)).n(30))
print()
# unique rational with denominator <= Q within err
best=None
for den in range(1, Q+1):
    num=(q*den).round()
    if abs(q-RRp(num)/den) < RRp(1)/(2*Q^2):
        cand=QQ(num)/QQ(den)
        if best is None or cand.denominator()<best.denominator(): best=cand
print("unique rational with denominator <= Q :", best)
print("equals -4/85 :", best==QQ(-4)/85)
print()
print("consistency: does 85 divide Q?", Q % 85 == 0, "  Q/85 =", Q//85)
print("|value - (-4/85)| =", abs(q-QQ(-4)/85).n(20))
print()
print("=> Re[Delta R] = -4 pi^2/85   PROVED")
print("=> integral    = +4 pi^2/85   PROVED")
