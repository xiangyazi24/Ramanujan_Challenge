from fractions import Fraction as F
import itertools, sys

def ch_coeffs(n):
    return (-8*n**3-51*n**2-105*n-68,
            24*n**5+337*n**4+1833*n**3+4818*n**2+6092*n+2928,
            -(n+2)*(n+3)*(24*n**5+273*n**4+1150*n**3+2154*n**2+1635*n+268),
            (n+1)*(n+2)**4*(n+3)*(8*n**3+75*n**2+231*n+232))

def ch_solve(init, N):
    u={-3:F(init[0]),-2:F(init[1]),-1:F(init[2])}
    for n in range(0,N+1):
        c0,c1,c2,c3=ch_coeffs(n)
        u[n]=-(c1*u[n-1]+c2*u[n-2]+c3*u[n-3])/c0
    return u

N=45
Q = ch_solve([1,12,306],N); P = ch_solve([0,7,179],N)
# index by m = n+3
def q(m): return Q[m-3]
def p(m): return P[m-3]

def find_order2(seq, DEG, MMIN=2, NEQ=None):
    """Look for a(m) s(m) + b(m) s(m-1) + c(m) s(m-2) = 0 with deg <= DEG."""
    nun = 3*(DEG+1)
    if NEQ is None: NEQ = nun + 6
    rows=[]
    for m in range(MMIN, MMIN+NEQ):
        row=[]
        for shift,val in ((0,seq(m)),(1,seq(m-1)),(2,seq(m-2))):
            for d in range(DEG+1):
                row.append(F(m)**d * val)
        rows.append(row)
    # nullspace over Q by gaussian elimination
    import copy
    A=[r[:] for r in rows]
    ncols=nun
    piv=[]
    r=0
    for c in range(ncols):
        pr=None
        for i in range(r,len(A)):
            if A[i][c]!=0: pr=i;break
        if pr is None: continue
        A[r],A[pr]=A[pr],A[r]
        pv=A[r][c]
        A[r]=[x/pv for x in A[r]]
        for i in range(len(A)):
            if i!=r and A[i][c]!=0:
                f=A[i][c]
                A[i]=[a-f*b for a,b in zip(A[i],A[r])]
        piv.append(c); r+=1
        if r==len(A): break
    free=[c for c in range(ncols) if c not in piv]
    return len(free), piv, A, free

for name,seq in (("q",q),("p",p)):
    print("=== sequence %s ===" % name)
    for DEG in range(1,11):
        nfree,_,_,_ = find_order2(seq,DEG)
        print("   order-2 ansatz, coeff deg <= %2d : nullspace dim = %d" % (DEG,nfree))
        if nfree>0: break
