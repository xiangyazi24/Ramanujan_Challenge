from sage.all import *

Kn.<n> = PolynomialRing(QQ)
K = Kn.fraction_field()
R.<p,q,v> = PolynomialRing(K)
D = p*q*(1+v^2) + 2*v
Dp, Dq, Dv = D.derivative(p), D.derivative(q), D.derivative(v)

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

M = matrix(R, [
 [(2*n+5)*(n+3)^2*(136*n^4+1424*n^3+5548*n^2+9551*n+6141),
  384*n^6+6384*n^5+44168*n^4+162698*n^3+336377*n^2+369933*n+169011,
  480*n^4+4980*n^3+19210*n^2+32690*n+20730],
 [(n+2)^2*(n+3)^2*(4*n+10)*(48*n^3+386*n^2+1017*n+879),
  (n+2)^2*(272*n^5+3848*n^4+21732*n^3+61184*n^2+85761*n+47808),
  (n+2)^2*(320*n^3+2540*n^2+6610*n+5640)],
 [(4*n+10)*(n+2)^2*(n+3)^2*(32*n^4+302*n^3+1037*n^2+1530*n+813),
  (n+2)^2*(192*n^6+2984*n^5+19116*n^4+64452*n^3+120256*n^2+117279*n+46476),
  (n+2)^2*(16*n^5+408*n^4+2912*n^3+8884*n^2+12254*n+6240)]])

Snum = p^2*q^2*(1-p^2)*(1-q^2)*v^2
next_num = [Snum*D^2,
            Snum*2*(n+3)*v*D,
            Snum*((-(n+3)*v)*D + 2*(n+3)*(2*n+7)*v^2)]
cur_num = [D^4,
           2*(n+2)*v*D^3,
           (-(n+2)*v*D + 2*(n+2)*(2*n+5)*v^2)*D^2]

degree = 12
basis = [p^a*q^b*v^c for a in range(degree+1)
         for b in range(degree+1-a)
         for c in range(degree+1-a-b) if (a+b) % 2 == 0]
images = [[op_p(m) for m in basis],
          [op_q(m) for m in basis],
          [op_v(m) for m in basis]]
cols = images[0] + images[1] + images[2]
lam = (n+1)*(n+2)^2*(n+3)^2*(2*n+7)^2
rhs = [sum(M[i,j]*next_num[j] for j in range(3)) - lam*cur_num[i]
       for i in range(3)]
exps = set()
for e in cols + rhs:
    exps.update(e.dict().keys())
exps = sorted(exps)
idx = {e:k for k,e in enumerate(exps)}
entries = {}
for j,e in enumerate(cols):
    for exp,c in e.dict().items():
        if c:
            entries[(idx[exp],j)] = c
A = matrix(K, len(exps), len(cols), entries, sparse=True)

def ev2(c):
    return QQ(c.numerator()(2)) / QQ(c.denominator()(2))

A2 = matrix(QQ, A.nrows(), A.ncols(),
            {(i,j):ev2(c) for (i,j),c in A.dict().items()}, sparse=True)

def solve_row(i):
    bvec = vector(K, [rhs[i].dict().get(exp,0) for exp in exps])
    b2 = vector(QQ, [ev2(c) for c in bvec])
    print('row', i, 'shape', A.nrows(), A.ncols(), flush=True)
    sol2 = A2.solve_right(b2)
    support = [j for j,x in enumerate(sol2) if x != 0]
    B2 = A2.matrix_from_columns(support)
    rows = list(B2.transpose().pivots())
    assert len(rows) == len(support)
    B = matrix(K, [[A[r,c] for c in support] for r in rows])
    small_b = vector(K, [bvec[r] for r in rows])
    small_sol = B.solve_right(small_b)
    sol = vector(K, A.ncols())
    for j,x in zip(support,small_sol): sol[j] = x
    assert A*sol == bvec
    certs = [sum(sol[t*len(basis)+j]*basis[j] for j in range(len(basis)))
             for t in range(3)]
    assert op_p(certs[0])+op_q(certs[1])+op_v(certs[2]) == rhs[i]
    print('lambda =', factor(lam), flush=True)
    for name,P in zip(['Pp','Pq','Pv'], certs):
        print(name, 'terms', len(P.dict()), 'degree', P.degree())
        print(P)
    return certs

for i in range(3):
    solve_row(i)
