from sage.all import *


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


def solve(nn, row, lower, degree):
    L = LaurentPolynomialRing(QQ, names=("p", "q", "v"))
    p, q, v = L.gens()
    D = p*q*(1+v^2)+2*v
    exponent = 2*nn+8

    def op(variable, boundary, logpart, poly):
        return D*(boundary*poly.derivative(variable)+logpart*poly) - \
            exponent*boundary*D.derivative(variable)*poly

    ops = [
        lambda P: op(p, p*(1-p^2), (2*nn+7)-(4*nn+9)*p^2, P),
        lambda P: op(q, q*(1-q^2), (2*nn+6)-(4*nn+10)*q^2, P),
        lambda P: op(v, v*(1-v^2), (2*nn+4)-(2*nn+6)*v^2, P),
    ]
    M = matrix_at(nn)
    S = p^2*q^2*(1-p^2)*(1-q^2)*v^2
    nxt = [S*D^2, S*2*(nn+3)*v*D,
           S*(-(nn+3)*v*D+2*(nn+3)*(2*nn+7)*v^2)]
    cur = [D^4, 2*(nn+2)*v*D^3,
           (-(nn+2)*v*D+2*(nn+2)*(2*nn+5)*v^2)*D^2]
    lam = (nn+1)*(nn+2)^2*(nn+3)^2*(2*nn+7)^2
    rhs = sum(M[row,j]*nxt[j] for j in range(3))-lam*cur[row]

    exponents = [(a,b,c)
                 for a in range(lower, degree+1)
                 for b in range(lower, degree+1)
                 for c in range(lower, degree+1)
                 if a+b+c <= degree and (a+b) % 2 == 0]
    basis = [p^a*q^b*v^c for a,b,c in exponents]
    cols = [operator(monomial) for operator in ops for monomial in basis]
    rows = set(rhs.dict())
    for col in cols:
        rows.update(col.dict())
    rows = sorted(rows)
    index = {exp:i for i,exp in enumerate(rows)}
    entries = {}
    for j,col in enumerate(cols):
        for exp,coefficient in col.dict().items():
            if coefficient:
                entries[index[exp],j] = coefficient
    A = matrix(QQ,len(rows),len(cols),entries,sparse=True)
    bvec = vector(QQ,[rhs.dict().get(exp,0) for exp in rows])
    try:
        solution=A.solve_right(bvec)
    except ValueError:
        return None,(A.nrows(),A.ncols())
    return solution,(A.nrows(),A.ncols())


for lower in [-1,-2,-3]:
    for degree in [6,8,10,12]:
        answer,shape=solve(2,0,lower,degree)
        print("lower",lower,"degree",degree,"found",answer is not None,
              "shape",shape,
              "support",None if answer is None else sum(x != 0 for x in answer),
              flush=True)
        if answer is not None:
            raise SystemExit
