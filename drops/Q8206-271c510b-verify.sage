#!/usr/bin/env sage
"""Q8206 exact verifier. No harmonic_number() and no floating arithmetic.

Run: sage drops/Q8206-271c510b-verify.sage [--prime-limit 150]
Finite loops are regression/counterexample searches, not proofs.
"""
from argparse import ArgumentParser
from sage.all import ZZ, QQ, GF, PolynomialRing, binomial, factorial, prime_range


def Htab(n):
    H = [QQ(0)]
    for j in range(1, int(n) + 1):
        H.append(H[-1] + QQ(1) / ZZ(j))
    return H


def b(n):
    return ZZ(sum(binomial(n,k)^2 * binomial(n+k,k)^2 for k in range(n+1)))


def D(n):
    H = Htab(2*n); ans = QQ(0)
    for k in range(n+1):
        t = ZZ(binomial(n,k)^2 * binomial(n+k,k)^2)
        ans += QQ(t) * (H[n+k] - H[n-k])
    return ans


def P(n): return ZZ(34*n^3 + 51*n^2 + 27*n + 5)
def Pd(n): return ZZ(102*n^2 + 102*n + 27)


def qmod(x, p):
    x = QQ(x); p = ZZ(p)
    a = ZZ(x.numerator()) % p; d = ZZ(x.denominator()) % p
    assert d != 0, (x, p)
    return a * d.inverse_mod(p) % p


def compose(f, u, R):
    out = R.zero(); power = R.one()
    for c in f.list():
        out += R(c) * power; power *= u
    return out


def build(p):
    p = ZZ(p); N = (p-1)//2; K = GF(p); R = PolynomialRing(K,'Y'); Y=R.gen()
    phi=R.one(); F=R.one()
    for k in range(1, int(N)+1):
        phi *= (Y-K((k-1)*k))/K(k*k); F += phi^2
    psi=R.one()
    for a in range(int(N)): psi *= Y-K(a*(a+1))
    assert psi == (K(4)*Y+1)^N-1
    assert K(factorial(N))^2 == K((-1)^(N+1))
    assert phi^2 == psi^2
    assert F.degree()==p-1 and F.leading_coefficient()==1 and F[0]==1
    return K,R,Y,N,F,psi


def symbolic_WZ():
    R=PolynomialRing(QQ,names=('z','k')); z,k=R.gens(); y=z*(z+1); x=2*z+1
    den=(k+z)^2*(k-z-1)^2; PP=34*z^3+51*z^2+27*z+5
    recurrence=((z+1)^3*(z+k+1)^2*(z+k)^2-PP*den
                +z^3*(z-k)^2*(k-z-1)^2)
    telescope=(4*x*(2*k^2+k-1-4*y)*den
               -4*k^4*x*(2*k^2-3*k-4*y))
    assert recurrence == telescope


def audit_p(p):
    K,R,Y,N,F,psi=build(p); p=ZZ(p); S=PolynomialRing(K,'z'); z=S.gen(); u=z*(z+1)
    B=compose(F,u,S); PP=34*z^3+51*z^2+27*z+5
    lhs=(z+1)^3*B(z+1)-PP*B+z^3*B(z-1)
    rhs=-16*(2*z+1)*(z^p-z)^2
    assert lhs==rhs
    x=2*z+1; prod=S.one()
    for a in range(int(N)+1): prod *= u-K(a*(a+1))
    assert K(4)^(N+1)*prod == x*(x^p-x)
    boundary=(prod/K(factorial(N+1))^2)^2
    assert boundary == x^2*(x^p-x)^2
    assert -4*x*(x^p-x)^2 == rhs
    Fp=F.derivative()
    for n in range(int(N)+1):
        lam=K(n*(n+1))
        assert F(lam)==K(b(n))
        assert K(2*n+1)/K(2)*Fp(lam)==K(qmod(D(n),p)), (p,n)
    G=F.gcd(psi).monic(); q,rem=F.quo_rem(G); assert rem==0
    C=psi*Fp-K(2)*psi.derivative()*F; c,rem=C.quo_rem(G); assert rem==0
    assert q.gcd(c).monic()==F.gcd(Fp).monic()
    return F.gcd(Fp).monic()


def Eval(n,B):
    if n==1: return ZZ(27)
    return ZZ(-3 + sum(Pd(j)*B[j]^2 for j in range(1,n))
              -6*sum(ZZ(j+1)^2*B[j]*B[j+1] for j in range(1,n-1)))


def fixed_audit(limit, prime_limit):
    M=max(limit+1,(prime_limit-1)//2+1)
    B=[b(n) for n in range(M+1)]; DD=[D(n) for n in range(M+1)]; g=[2*x for x in DD]
    assert (B[0],B[1],g[0],g[1])==(1,5,0,12)
    for n in range(1,M):
        source=Pd(n)*B[n]-3*(n+1)^2*B[n+1]-3*n^2*B[n-1]
        assert (n+1)^3*g[n+1]-P(n)*g[n]+n^3*g[n-1]==source
        W=n^3*(B[n-1]*g[n]-g[n-1]*B[n])
        assert W+3*n^2*B[n-1]*B[n]==Eval(n,B)
    # Hostile regression: the old float-based code failed here.
    assert DD[4]==QQ(104825)/2 and qmod(DD[4],11)==3
    K,R,Y,N,F,psi=build(11)
    assert K(9)/K(2)*F.derivative()(K(20))==K(3)
    for p in prime_range(3,prime_limit):
        for n in range(1,min(M,(p-2)//2+1)):
            if p>2*n+1 and B[n]%p==0:
                assert (Eval(n,B)%p==0)==(qmod(DD[n],p)==0)


def main():
    ap=ArgumentParser(); ap.add_argument('--prime-limit',type=int,default=150)
    ap.add_argument('--index-limit',type=int,default=200); a=ap.parse_args()
    symbolic_WZ(); fixed_audit(a.index_limit,a.prime_limit)
    failures=[]
    for p in prime_range(3,a.prime_limit):
        defect=audit_p(p)
        if defect.degree()>0:
            failures.append((ZZ(p),defect)); print('COUNTEREXAMPLE',failures[0]); break
    K,R,Y,N,F,psi=build(3); assert F==Y^2+1 and F.gcd(F.derivative())==1
    if failures: raise AssertionError(failures[0])
    print('PASS: symbolic WZ, exact seam/boundary, exact node normalization,')
    print('      E_n identities, saturated gcd, and finite squarefree regression')
    print('Finite evidence only: no squarefreeness failure for p <',a.prime_limit)


if __name__=='__main__': main()
