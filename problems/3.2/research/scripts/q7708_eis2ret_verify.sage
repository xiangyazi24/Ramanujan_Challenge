#!/usr/bin/env sage
"""Q7708 exact verifier: actual-source EIS-2RET identities.

This checks finite data and algebraic identities used in the Q7708 note:
  * actual g and Xi recurrences;
  * Franel pullback and first-block Cartier formula for g_m;
  * reflected Hasse-block coefficient formula for a two-return Green period;
  * p=41 zero-period guard and the p=2237 later-return nonzero tests;
  * low-q level-six identities H^2=E^2 Delta and M4=(5-t)Psi;
  * the modular constant-term formula g_m=CT(t^{-m}/E) and
    CT(M4*t^{-m}/E)=-delta_{m,1} for small m.

The universal theorems are proved in the accompanying note; these computations
are reproducibility checks, not evidence promoted to proof.
"""

from sage.all import *


def P(n):
    n = ZZ(n)
    return 34*n^3 + 51*n^2 + 27*n + 5


def apery_mod(p):
    Fp = GF(p)
    b = [Fp(1), Fp(5)]
    for n in range(2, p):
        b.append((Fp(P(n-1))*b[n-1] - Fp((n-1)^3)*b[n-2]) / Fp(n^3))
    return b


def qsqrt_inv_mod(p):
    Fp = GF(p)
    q = [Fp(1), Fp(17)]
    for n in range(1, p-1):
        q.append((Fp(17*(2*n+1))*q[n] - Fp(n)*q[n-1]) / Fp(n+1))
    return q


def actual_source_mod(p, b):
    Fp = GF(p)
    R = PowerSeriesRing(Fp, 't', default_prec=p)
    F = R(b).add_bigoh(p)
    Q = R(qsqrt_inv_mod(p)).add_bigoh(p)
    G = (Q/(F*F)).add_bigoh(p)
    return [Fp(G[n]) for n in range(p)], G


def xi_mod(p, b, g):
    Fp = GF(p)
    xi = [Fp(-1)]
    for n in range(1, p):
        xi.append(xi[-1] - Fp(5)*g[n]*b[n-1])
    return xi


def franel_mod(p):
    Fp = GF(p)
    f = [Fp(1), Fp(2)]
    for n in range(1, p-1):
        f.append((Fp(7*n*n+7*n+2)*f[n] + Fp(8*n*n)*f[n-1]) / Fp((n+1)^2))
    return f


def check_pullback_cartier(p, sample_rows=None):
    Fp = GF(p)
    b = apery_mod(p)
    g, _ = actual_source_mod(p, b)
    hcoef = franel_mod(p)
    R = PowerSeriesRing(Fp, 'x', default_prec=p)
    x = R.gen()
    Hp = R(hcoef).add_bigoh(p)
    invH4 = (Hp^4).inverse_of_unit()
    phi = (x*(1-Fp(8)*x)/(1+x)).add_bigoh(p)
    rho = ((1+x)^(3*p-3) * invH4).add_bigoh(p)

    if sample_rows is None:
        hz = [r for r in range(1,p) if b[r] == 0]
        sample_rows = sorted(set([0,1,2,p-2,p-1] + hz))
    for m in sample_rows:
        assert 0 <= m < p
        # Lagrange/Franel coefficient formula.
        lag = ((1+x)^(m-2) * (1-Fp(8)*x)^(-m-1) * invH4).add_bigoh(m+1)
        assert lag[m] == g[m], (p,m,'lagrange',lag[m],g[m])
        # First-block Frobenius/Cartier formula.
        car = (rho * phi^(p-m-1)).add_bigoh(p)
        assert car[p-1] == g[m], (p,m,'cartier',car[p-1],g[m])
    print('PULLBACK_CARTIER', p, tuple(sample_rows))


def source_factor_series(p):
    Fp = GF(p)
    b = apery_mod(p)
    R = PowerSeriesRing(Fp, 't', default_prec=p)
    t = R.gen()
    A = R(b).add_bigoh(p)
    D = R(1) - Fp(34)*t + t^2
    Rp = (D^((p-1)//2) / (A*A)).add_bigoh(p)
    g, G = actual_source_mod(p,b)
    assert all(Rp[n] == g[n] for n in range(p))
    return b,g,xi_mod(p,b,g),Rp


def block_poly(p,b,r,s,R):
    t=R.gen()
    C=R(0)
    for m in range(r+1,s+1):
        C += b[m-1]*t^(p-1-m)
    return C


def check_return_pair(p,r,s,expect_zero=None):
    Fp=GF(p)
    b,g,xi,Rp=source_factor_series(p)
    assert 0<r<s<p
    assert b[r]==0 and b[s]==0
    assert b[r-1]!=0 and b[r+1]!=0 and b[s-1]!=0 and b[s+1]!=0
    S=sum((b[m-1]*g[m] for m in range(r+1,s+1)),Fp(0))
    assert xi[s]-xi[r] == -Fp(5)*S

    R=Rp.parent(); t=R.gen()
    C=block_poly(p,b,r,s,R)
    coeff=(Rp*C).add_bigoh(p)[p-1]
    assert coeff==S

    # Hasse reciprocity turns C into one contiguous reflected coefficient block.
    Crefl=R(0)
    for k in range(p-s,p-r):
        Crefl += b[k]*t^(k-1)
    assert C==Crefl
    assert b[p-r-1]==b[r]==0
    assert b[p-s-1]==b[s]==0
    h=s-r
    inner=R(0)
    for j in range(h-1):
        inner += b[s-1-j]*t^j
    assert C == t^(p-s-1)*inner
    assert inner[0]==b[s-1]!=0
    assert inner[h-2]==b[r+1]!=0

    if expect_zero is not None:
        assert (S==0) == expect_zero
    print('RETURN_PAIR',p,r,s,'XI',int(xi[r]),int(xi[s]),'S',int(S),'BLOCK_DEG',h-2)
    return int(S)


# Mandatory actual-source guard: a zero Green interval need not start at Xi=0.
assert check_return_pair(41,10,30,True)==0

# Known common prime 2237 has two later Hasse returns.  EIS-2RET requires both
# corresponding actual Cartier-block coefficients to be nonzero; verify exactly.
b2237,g2237,xi2237,R2237=source_factor_series(2237)
hz2237=[r for r in range(1,2237) if b2237[r]==0]
assert hz2237 == [23,492,1744,2213]
assert xi2237[492]==0
assert [r for r in hz2237 if xi2237[r]==0] == [492]
S1=check_return_pair(2237,492,1744,False)
S2=check_return_pair(2237,492,2213,False)
print('P2237_LATER_SHOOTING',S1,S2)

# A few exact pullback/Cartier checks, including all Hasse rows in these primes.
for p in (7,17,41,181):
    check_pullback_cartier(p)


# ---------- small characteristic-zero q-series audit ----------

def balanced_product(R, factors, prec):
    if not factors:
        return R(1).add_bigoh(prec)
    layer=[R(f).add_bigoh(prec) for f in factors]
    while len(layer)>1:
        nxt=[]
        for i in range(0,len(layer),2):
            nxt.append((layer[i]*layer[i+1]).add_bigoh(prec) if i+1<len(layer) else layer[i])
        layer=nxt
    return layer[0]


def sparse_factor(R,step,e,sign,prec):
    q=R.gen(); z=R(0)
    for j in range(e+1):
        if step*j>=prec: break
        z += QQ(binomial(e,j)*sign^j)*q^(step*j)
    return z


def product_quotient(R,specs,prec):
    num=[]; den=[]
    for d,e,sign in specs:
        for n in range(1,(prec-1)//d+1):
            (num if e>0 else den).append(sparse_factor(R,d*n,abs(e),sign,prec))
    return (balanced_product(R,num,prec)/balanced_product(R,den,prec)).add_bigoh(prec)


def e4_series(R,d,prec):
    q=R.gen(); out=R(1)
    for n in range(1,(prec-1)//d+1):
        out += 240*sigma(n,3)*q^(d*n)
    return out.add_bigoh(prec)


def check_level6_qseries(prec=14):
    R=PowerSeriesRing(QQ,'q',default_prec=prec); q=R.gen()
    T=product_quotient(R,((3,12,1),(1,-12,1)),prec)
    t=(q*T).add_bigoh(prec)
    E=product_quotient(R,((2,7,-1),(3,7,-1),(1,-5,-1),(6,-5,-1)),prec)
    H=(1+q*T.derivative()/T).add_bigoh(prec)
    Psi=(E*H).add_bigoh(prec)
    Delta=(1-34*t+t*t).add_bigoh(prec)
    assert (H*H-E*E*Delta).add_bigoh(prec)==0

    M4=(-3*e4_series(R,1,prec)+4*e4_series(R,2,prec)-9*e4_series(R,3,prec)+108*e4_series(R,6,prec))/20
    M4=M4.add_bigoh(prec)
    assert (M4-(5-t)*Psi).add_bigoh(prec)==0

    # E=F(t), to the available precision.
    B=[ZZ(1),ZZ(5)]
    for n in range(2,prec):
        B.append((P(n-1)*B[n-1]-(n-1)^3*B[n-2])//(n^3))
    Ft=R(0)
    for n in range(prec):
        Ft += B[n]*t^n
    assert (E-Ft).add_bigoh(prec)==0

    # g(t(q))=1/Psi and g_m=CT_q(t^{-m}/E).
    Rt=PowerSeriesRing(QQ,'z',default_prec=8); z=Rt.gen()
    Fsmall=Rt(B[:8]); Dsmall=Rt(1)-34*z+z^2
    Gsmall=(Dsmall^QQ(-1)/2)/(Fsmall*Fsmall)
    # Sage's fractional power above can be version-sensitive; independent
    # recurrence/inversion path for exact coefficients:
    qs=[QQ(1),QQ(17)]
    for n in range(1,7):
        qs.append((QQ(17*(2*n+1))*qs[n]-QQ(n)*qs[n-1])/QQ(n+1))
    Qsmall=Rt(qs[:8]); Gsmall=(Qsmall/(Fsmall*Fsmall)).add_bigoh(8)
    assert ((Gsmall(t))*Psi-1).add_bigoh(8)==0
    for m in range(1,6):
        A=(T^(-m)/E).add_bigoh(m+1)
        assert A[m]==Gsmall[m]
        MA=(M4*T^(-m)/E).add_bigoh(m+1)
        assert MA[m] == (-1 if m==1 else 0)
    print('LEVEL6_QSERIES',prec,'M4_OVER_PSI=5-t','CT_G_OK')


check_level6_qseries()
print('Q7708_EIS2RET_VERIFY=PASS')
