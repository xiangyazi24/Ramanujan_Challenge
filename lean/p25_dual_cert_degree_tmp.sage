from sage.all import *

R.<p,q,v> = PolynomialRing(QQ)
D = p*q*(1+v^2) + 2*v
Dp, Dq, Dv = D.derivative(p), D.derivative(q), D.derivative(v)

def matrix_at(n):
    return matrix(QQ, [
      [(2*n+5)*(n+3)^2*(136*n^4+1424*n^3+5548*n^2+9551*n+6141),
       384*n^6+6384*n^5+44168*n^4+162698*n^3+336377*n^2+369933*n+169011,
       480*n^4+4980*n^3+19210*n^2+32690*n+20730],
      [(n+2)^2*(n+3)^2*(4*n+10)*(48*n^3+386*n^2+1017*n+879),
       (n+2)^2*(272*n^5+3848*n^4+21732*n^3+61184*n^2+85761*n+47808),
       (n+2)^2*(320*n^3+2540*n^2+6610*n+5640)],
      [(4*n+10)*(n+2)^2*(n+3)^2*(32*n^4+302*n^3+1037*n^2+1530*n+813),
       (n+2)^2*(192*n^6+2984*n^5+19116*n^4+64452*n^3+120256*n^2+117279*n+46476),
       (n+2)^2*(16*n^5+408*n^4+2912*n^3+8884*n^2+12254*n+6240)]])

def ops(n):
    def op_p(P):
        return D*(p*(1-p^2)*P.derivative(p) +
                  ((2*n+7)-(4*n+9)*p^2)*P) - \
               (2*n+7)*p*(1-p^2)*Dp*P
    def op_q(P):
        return D*(q*(1-q^2)*P.derivative(q) +
                  ((2*n+6)-(4*n+10)*q^2)*P) - \
               (2*n+7)*q*(1-q^2)*Dq*P
    def op_v(P):
        return D*(v*(1-v^2)*P.derivative(v) +
                  ((2*n+4)-(2*n+6)*v^2)*P) - \
               (2*n+7)*v*(1-v^2)*Dv*P
    return op_p,op_q,op_v

def target(n,i):
    M=matrix_at(n)
    S=p^2*q^2*(1-p^2)*(1-q^2)*v^2
    nxt=[S*D^2,S*2*(n+3)*v*D,
         S*(-(n+3)*v*D+2*(n+3)*(2*n+7)*v^2)]
    cur=[D^4,2*(n+2)*v*D^3,
         (-(n+2)*v*D+2*(n+2)*(2*n+5)*v^2)*D^2]
    lam=(n+1)*(n+2)^2*(n+3)^2*(2*n+7)^2
    return sum(M[i,j]*nxt[j] for j in range(3))-lam*cur[i]

def solvable(n,i,d,which=(0,1,2)):
    op=ops(n)
    basis=[p^a*q^b*v^c for a in range(d+1)
           for b in range(d+1-a) for c in range(d+1-a-b)
           if (a+b)%2==0]
    cols=[]
    for j in which: cols += [op[j](m) for m in basis]
    rhs=target(n,i)
    exps=set(rhs.dict())
    for e in cols: exps.update(e.dict())
    exps=sorted(exps)
    idx={e:k for k,e in enumerate(exps)}
    ent={}
    for j,e in enumerate(cols):
        for exp,c in e.dict().items(): ent[idx[exp],j]=c
    A=matrix(QQ,len(exps),len(cols),ent,sparse=True)
    b=vector(QQ,[rhs.dict().get(e,0) for e in exps])
    try:
        sol=A.solve_right(b)
        nz=sum(x!=0 for x in sol)
        if d == 12:
            print('support', i, [(j//len(basis), next(iter(basis[j % len(basis)].dict())))
                                 for j,x in enumerate(sol) if x != 0], flush=True)
        return (True,A.nrows(),A.ncols(),nz)
    except ValueError:
        return (False,A.nrows(),A.ncols(),None)

for d in range(12,13):
    print('degree',d,flush=True)
    for i in range(3):
        print(i,solvable(2,i,d),flush=True)
