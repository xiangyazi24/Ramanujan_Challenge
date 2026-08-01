import sympy as sp
from sympy import Poly, gcd, resultant, factor_list, I
X = sp.Symbol('X')

def A_(Y): return 34*Y**3 - 51*Y**2 + 27*Y - 5

def build_N(hmax):
    # u_{n} = [A(n) u_{n-1} - (n-1)^3 u_{n-2}]/n^3, n = X+2 .. X+h
    # u_{X+h} = P_h(X) u_{X+1} + Q_h(X) u_X ; N_h = P_h * prod_{j=2..h}(X+j)^3 (up to content)
    P = {0: sp.Integer(0), 1: sp.Integer(1)}   # u_{X+0}=u_X -> P=0; u_{X+1} -> P=1
    Q = {0: sp.Integer(1), 1: sp.Integer(0)}
    Ns = {}
    for h in range(2, hmax+1):
        n = X + h
        P[h] = sp.expand((A_(n)*P[h-1] - (n-1)**3*P[h-2])/n**3)
        Q[h] = sp.expand((A_(n)*Q[h-1] - (n-1)**3*Q[h-2])/n**3)
        num = sp.together(P[h])
        num = sp.expand(sp.numer(num))
        # clear content
        pol = Poly(num, X)
        cont = pol.content()
        Ns[h] = Poly(pol / cont, X)
    return Ns

Ns = build_N(12)
print("deg check:", all(Ns[h].degree() == 3*h-3 for h in Ns), [Ns[h].degree() for h in sorted(Ns)])
N3ref = Poly(1155*X**6+13860*X**5+68535*X**4+178680*X**3+259059*X**2+198156*X+62531, X)
print("N3 matches banked:", Ns[3] == N3ref or Ns[3] == -N3ref)
# reflection identity
refl_ok = []
for h in Ns:
    lhs = Poly(sp.expand(Ns[h].as_expr().subs(X, -(h+1)-X)), X)
    sign = (-1)**(h+1)
    refl_ok.append(lhs == Poly(sp.expand(sign*Ns[h].as_expr()), X) or lhs == Poly(sp.expand(-sign*Ns[h].as_expr()), X))
print("reflection identity (up to overall sign):", all(refl_ok))
# center root even h
cent = {h: Ns[h].as_expr().subs(X, sp.Rational(-(h+1),2)) for h in Ns if h%2==0}
print("center root even h (all zero?):", {h: v == 0 for h,v in cent.items()})
# all-pairs gcd
bad = []
for h1 in Ns:
    for h2 in Ns:
        if h1 < h2:
            g = gcd(Ns[h1].as_expr(), Ns[h2].as_expr())
            if sp.degree(g, X) > 0: bad.append((h1,h2,g))
print("all-pairs gcd nontrivial:", bad if bad else "NONE (all coprime)")
# descent-shift gcd examples
badsh = []
for (a,b,c) in [(2,2,2),(3,2,3),(2,3,2),(4,3,4),(3,4,3),(5,3,5)]:
    g = gcd(Ns[a].as_expr(), Ns[b].as_expr().subs(X, X+c))
    if sp.degree(g, X) > 0: badsh.append((a,b,c,g))
print("descent-shift gcd nontrivial:", badsh if badsh else "NONE")
# irreducibility
irr = {}
for h in Ns:
    fl = factor_list(Ns[h].as_expr())[1]
    irr[h] = len(fl) == 1 and fl[0][1] == 1
print("N_h irreducible over Q:", irr)
# resultant sizes (height law)
import math
for (a,b) in [(2,3),(2,4),(3,4),(4,5),(5,6)]:
    r = resultant(Ns[a].as_expr(), Ns[b].as_expr(), X)
    print(f"log|Res(N_{a},N_{b})| = {math.log(abs(int(r))):.1f}")
