"""Measure the tail order of the CRIT-2H certificate objects:
   how fast do the cell-coordinate roots of J_h converge as h grows?
   (decisive for whether an effective h0 exists: need o(h^-2) relative tail)"""
import sympy as sp, mpmath as mp
from sympy import Poly, ZZ
mp.mp.dps=40
X=sp.Symbol('X'); U=sp.Symbol('U')
def A_(t): return 34*t**3-51*t**2+27*t-5
def build(H):
    N={0:Poly(0,X,domain=ZZ),1:Poly(1,X,domain=ZZ)}
    for h in range(2,H+1):
        N[h]=Poly(A_(X+h),X,domain=ZZ)*N[h-1]-Poly((X+h-1)**6,X,domain=ZZ)*N[h-2]
    return {h:N[h].primitive()[1] for h in range(2,H+1)}
NN=build(30)
def J_of(h):
    Nh=NN[h].as_expr()
    q=sp.prod([X+j for j in range(1,h+1)])
    A=sp.expand(q*sp.diff(Nh,X)-3*sp.diff(q,X)*Nh)
    # A is odd/even under X -> -(h+1)-X ; in s=2X+h+1 it should be a poly in u=s^2 times s^eps
    s=2*X+h+1
    # substitute X = (s-h-1)/2 and check parity in s
    S=sp.Symbol('S')
    As=sp.expand(A.subs(X,(S-h-1)/2))
    p=Poly(As,S)
    co=p.all_coeffs()[::-1]  # ascending
    odd=[c for i,c in enumerate(co) if i%2==1 and c!=0]
    even=[c for i,c in enumerate(co) if i%2==0 and c!=0]
    if not odd: pol=[co[i] for i in range(0,len(co),2)]; shift=0
    elif not even: pol=[co[i] for i in range(1,len(co),2)]; shift=1
    else: return None,None
    return Poly(list(reversed(pol)),U), shift
for h in [8,12,16,20,24]:
    J,shift=J_of(h)
    if J is None: print(h,"parity FAIL"); continue
    rts=J.nroots(n=30,maxsteps=400)
    reals=sorted([sp.re(r) for r in rts if abs(sp.im(r))<1e-20])
    print(f"h={h}: deg J={J.degree()} (expect {h-1}) shift={shift} #real u-roots={len(reals)} smallest3={[sp.N(r,10) for r in reals[:3]]}")
