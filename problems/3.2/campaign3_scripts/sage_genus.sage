R.<r,s> = QQ[]
def Ppoly(u): return 34*u^3 + 51*u^2 + 27*u + 5
def Npoly(h, v):
    if h == 1: return R(1)
    N1, N2 = R(1), Ppoly(v+1)
    if h == 2: return N2
    a, b = N1, N2
    for m in range(2, h):
        a, b = b, Ppoly(v+m)*b - (v+m)^6*a
    return b
def Dpoly(h, v): return prod([(v+j)^3 for j in range(1, h+1)])

for (h,k) in [(2,3),(2,4),(3,4),(2,5)]:
    F = Npoly(h,r)*Dpoly(k,s) - Npoly(k,s)*Dpoly(h,r)
    A2.<x,y> = AffineSpace(QQ,2)
    Faff = F.subs(r=x, s=y)
    C = Curve(Faff)
    try:
        g = C.geometric_genus()
    except Exception as e:
        g = f"ERR {e}"
    # geometric irreducibility
    Fbar = F.change_ring(QQbar)
    nfac = "?"
    print(f"X_{h},{k}: bidegree=({F.degree(r)},{F.degree(s)}) geometric_genus={g}")
