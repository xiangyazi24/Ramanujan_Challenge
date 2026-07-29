R.<x> = PolynomialRing(ZZ)

def P(t):
    return (2*t + 1)*(17*t^2 + 17*t + 5)

def Ns(H):
    N = [R(0), R(1)]
    for r in range(1, H):
        N.append(P(x+r)*N[r] - (x+r)^6*N[r-1])
    return N

def Q(h):
    if h <= 1:
        return R(1)
    return prod((x+r)^3 for r in range(2, h+1))

def D(h,k,N):
    return (
        N[h]*(Q(k)//Q(h)) - N[k]
        + (x+1)^3*Q(h)*N[k-h](x=x+h)
    )

def primitive_tilde(h,k,N):
    f = D(h,k,N)
    forced = []
    if h % 2 == 0:
        forced.append(2*x+h+1)
    if k % 2 == 0:
        forced.append(2*x+k+1)
    if (k-h) % 2 == 0:
        forced.append(2*x+h+k+1)
    for L in forced:
        q,r = f.quo_rem(L)
        assert r == 0
        f = q
    f = f // gcd(f.list())
    if f.leading_coefficient() < 0:
        f = -f
    return f

def Rcert(h,j,k,N):
    A = primitive_tilde(h,j,N)
    B0 = primitive_tilde(j-h,k-h,N)
    B = B0(x=x+h)
    return A.resultant(B)

N = Ns(7)
for T in [(1,3,4),(1,3,5),(2,3,5),(1,4,6)]:
    r = Rcert(*T,N)
    print(T)
    print(r)
    print(factor(abs(r)))
