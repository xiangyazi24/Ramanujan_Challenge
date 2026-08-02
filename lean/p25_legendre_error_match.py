#!/usr/bin/env python3
"""Compare the target backward solution with Legendre-square Catalan errors."""
from fractions import Fraction as Q


def add(x, y): return (x[0] + y[0], x[1] + y[1])
def mul(c, x): return (c*x[0], c*x[1])


def catalan_moment(s):
    """Integral x^(2s) (-log x)/(1+x^2), as const + coeff*G."""
    partial = sum((Q((-1)**k, (2*k+1)**2) for k in range(s)), Q(0))
    return (Q((-1)**(s+1))*partial, Q((-1)**s))


def poly_mul(a, b):
    c = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): c[i+j] += x*y
    return c


def legendre_x2(m):
    """P_m(sqrt(1-8x^2))^2 as coefficients in X=x^2."""
    # P_0=1, P_1=t; recurrence, store parity-adjusted polynomial in X after t^2=1-8X.
    # First build P_m(t) in t.
    p0 = [Q(1)]
    if m == 0: p = p0
    else:
        p1 = [Q(0), Q(1)]
        if m == 1: p = p1
        else:
            for n in range(1, m):
                tp = [Q(0)] + [Q(2*n+1)*x for x in p1]
                p2 = [(tp[i] if i < len(tp) else 0) -
                      (Q(n)*(p0[i] if i < len(p0) else 0))
                      for i in range(max(len(tp), len(p0)))]
                p2 = [x/Q(n+1) for x in p2]
                while p2 and p2[-1] == 0: p2.pop()
                p0, p1 = p1, p2
            p = p1
    sq = poly_mul(p, p)
    # sq has only even t powers. Substitute t^(2r)=(1-8X)^r.
    out = [Q(0)]*(m+1)
    import math
    for e, c in enumerate(sq):
        if not c: continue
        assert e % 2 == 0
        r=e//2
        for j in range(r+1): out[j] += c*Q(math.comb(r,j))*(-8)**j
    return out


def beukers(m):
    # Integral x^2 weight times square: exponent x^(2(j+1)).
    ans=(Q(0),Q(0))
    for j,c in enumerate(legendre_x2(m)):
        ans=add(ans,mul(c,catalan_moment(j+1)))
    return ans


def target(n):
    n=Q(n)
    return [
        [(-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
         384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
         -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)],
        [(n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
         (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
         (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)],
        [(-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
         (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
         (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)],
    ]


def inv(a):
    x=[row[:] + [Q(i==j) for j in range(3)] for i,row in enumerate(a)]
    for j in range(3):
        p=next(i for i in range(j,3) if x[i][j]);x[j],x[p]=x[p],x[j]
        q=x[j][j];x[j]=[z/q for z in x[j]]
        for i in range(3):
            if i != j:
                q=x[i][j];x[i]=[x[i][k]-q*x[j][k] for k in range(6)]
    return [row[3:] for row in x]


def mv(a,v):
    out=[]
    for i in range(3):
        value=(Q(0),Q(0))
        for j in range(3): value=add(value,mul(a[i][j],v[j]))
        out.append(value)
    return out


w=[(Q(-23809,96),Q(4305,16)),
   (Q(-443527,384),Q(80535,64)),
   (Q(-7086881,1920),Q(257565,64))]
bs=[beukers(m) for m in range(30)]
for N in range(10):
    print('N',N)
    for j,x in enumerate(w):
        matches=[]
        for m,b in enumerate(bs):
            if b[0] and b[1] and x[0]*b[1] == x[1]*b[0]:
                matches.append((m,x[1]/b[1]))
        print(j, 'threshold', -x[0]/x[1], 'matches', matches)
    w=mv(inv(target(N)),w)
