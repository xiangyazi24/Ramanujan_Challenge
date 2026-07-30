R.<z> = PolynomialRing(QQ)
falpha = (z^12-3*z^11+4*z^10-5*z^9+6*z^8-7*z^7+7*z^6-7*z^5+6*z^4-5*z^3+4*z^2-3*z+1)
fbeta  = (z^16-7*z^15+22*z^14-48*z^13+87*z^12-133*z^11+178*z^10-211*z^9+223*z^8
          -211*z^7+178*z^6-133*z^5+87*z^4-48*z^3+22*z^2-7*z+1)
def cert(f, name, alpha_chart):
    print("="*70); print(name)
    K.<t> = NumberField(f)
    M, L = (t^2, t) if alpha_chart else (t, t)
    X = M^2
    u = (L + X^3)/(X*(L + X)); d = 1 + 4*u^2
    assert d.is_square()
    s = d.sqrt()
    up = R(u.polynomial().list()); sp = R(s.polynomial().list())
    print("  u(z) =", up)
    print("  s(z) =", sp)
    print("  s integral coeffs:", all(c in ZZ for c in sp.list()),
          " u integral coeffs:", all(c in ZZ for c in up.list()))
    rem = sp^2 - (1 + 4*up^2)
    q, r = rem.quo_rem(f)
    print("  s^2 - (1+4u^2) = c(z)*f(z)  exactly:", r == 0)
    print("  cofactor c(z) =", q)
    # also certify u equals the chart as a rational function
    numer = L + X^3; denom = X*(L+X)
    nn = R((numer).polynomial().list()) if hasattr(numer,'polynomial') else None
    print("  denom nonzero in K:", denom != 0)
    # chart identity: up * denom - numer == 0 mod f
    dnum = R((denom).polynomial().list()); nnum = R((numer).polynomial().list())
    chk = (up*dnum - nnum)
    q2, r2 = chk.quo_rem(f)
    print("  u(z)*denom(z) - numer(z) = c2(z)*f(z):", r2 == 0)
    return up, sp, q
cert(falpha, "ALPHA", True)
cert(fbeta,  "BETA",  False)
