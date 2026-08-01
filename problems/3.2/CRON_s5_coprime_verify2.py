import sympy as sp, math, time
from sympy import Poly, ZZ
X = sp.Symbol('X')
t0=time.time()
def A_(Y): return 34*Y**3 - 51*Y**2 + 27*Y - 5
HMAX=12
Np={0:Poly(0,X,domain=ZZ),1:Poly(1,X,domain=ZZ)}
for h in range(2,HMAX+1):
    Np[h]=Poly(A_(X+h),X,domain=ZZ)*Np[h-1]-Poly(((X+h-1))**6,X,domain=ZZ)*Np[h-2]
Ns={h:Np[h].primitive()[1] for h in range(2,HMAX+1)}
print("degrees:",{h:Ns[h].degree() for h in Ns},"expect 3h-3", flush=True)
N3ref=Poly(1155*X**6+13860*X**5+68535*X**4+178680*X**3+259059*X**2+198156*X+62531,X,domain=ZZ)
print("N3 == banked:",Ns[3]==N3ref or Ns[3]==-N3ref, flush=True)
if Ns[3]!=N3ref and Ns[3]!=-N3ref: print("  our N3:",Ns[3].as_expr())
# reflection
for h in list(Ns)[:6]:
    L=Poly(sp.expand(Ns[h].as_expr().subs(X,-(h+1)-X)),X,domain=ZZ)
    s=(-1)**(h+1)
    ok = (L==Ns[h]*s) or (L==Ns[h]*(-s))
    print(f"reflection h={h}: matches (-1)^(h+1) up to sign: {ok}", flush=True)
# center root even h
print("center roots:",{h:Ns[h].as_expr().subs(X,sp.Rational(-(h+1),2))==0 for h in Ns if h%2==0}, flush=True)
# all-pairs gcd
bad=[]
for h1 in Ns:
    for h2 in Ns:
        if h1<h2:
            g=Ns[h1].gcd(Ns[h2])
            if g.degree()>0: bad.append((h1,h2))
print("all-pairs gcd (h<=12):", bad if bad else "ALL COPRIME", flush=True)
# descent-shift gcds
badsh=[]
for (a,b,c) in [(2,2,2),(3,2,3),(2,3,2),(4,3,4),(3,4,3),(5,3,5),(6,4,6)]:
    sh=Poly(sp.expand(Ns[b].as_expr().subs(X,X+c)),X,domain=ZZ)
    g=Ns[a].gcd(sh)
    if g.degree()>0: badsh.append((a,b,c))
print("descent-shift gcd:", badsh if badsh else "ALL COPRIME", flush=True)
# irreducibility
irr={}
for h in Ns:
    fl=sp.factor_list(Ns[h].as_expr())[1]
    irr[h]=(len(fl)==1 and fl[0][1]==1)
    print(f"irred h={h}: {irr[h]} ({time.time()-t0:.0f}s)", flush=True)
# resultant heights
for (a,b) in [(2,3),(2,4),(3,4),(4,5)]:
    r=int(sp.resultant(Ns[a].as_expr(),Ns[b].as_expr(),X))
    print(f"log|Res(N_{a},N_{b})|={math.log(abs(r)):.1f}", flush=True)
print("total",time.time()-t0)
