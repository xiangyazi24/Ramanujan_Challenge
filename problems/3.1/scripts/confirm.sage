R.<z> = PolynomialRing(QQ)
fa = z^12-3*z^11+4*z^10-5*z^9+6*z^8-7*z^7+7*z^6-7*z^5+6*z^4-5*z^3+4*z^2-3*z+1
fb = (z^16-7*z^15+22*z^14-48*z^13+87*z^12-133*z^11+178*z^10-211*z^9+223*z^8
      -211*z^7+178*z^6-133*z^5+87*z^4-48*z^3+22*z^2-7*z+1)
CCp = ComplexField(300)
for f,name in [(fa,"f_alpha"),(fb,"f_beta")]:
    rs = f.roots(CCp, multiplicities=False)
    on = sum(1 for t in rs if abs(abs(t)-1) < 1e-80)
    rl = sum(1 for t in rs if abs(t.imag()) < 1e-80)
    print("%s: degree %d -> %d on |z|=1, %d real, total %d"%(name,f.degree(),on,rl,len(rs)))
    # g = the trace polynomial of a + 1/a
    S.<w> = PolynomialRing(QQ)
    K.<t> = NumberField(f)
    g = (t + 1/t).minpoly()
    print("   g(w) = %s   degree %d"%(g, g.degree()))
    rr = g.roots(RealField(300), multiplicities=False)
    print("   g totally real: %s;  roots in [-2,2]: %d;  outside: %d"
          %(len(rr)==g.degree(), sum(1 for r in rr if -2<=r<=2), sum(1 for r in rr if not(-2<=r<=2))))
    print("   => predicts %d on the circle, %d real  --  MATCHES: %s"
          %(2*sum(1 for r in rr if -2<=r<=2), 2*sum(1 for r in rr if not(-2<=r<=2)),
            on==2*sum(1 for r in rr if -2<=r<=2)))
    print()
