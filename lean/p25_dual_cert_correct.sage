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


def operators(n):
    # All three moments use the common denominator D^(2*n+8).
    exponent = 2*n + 8

    def op_p(P):
        return D*(p*(1-p^2)*P.derivative(p) +
                  ((2*n+7)-(4*n+9)*p^2)*P) - \
               exponent*p*(1-p^2)*Dp*P

    def op_q(P):
        return D*(q*(1-q^2)*P.derivative(q) +
                  ((2*n+6)-(4*n+10)*q^2)*P) - \
               exponent*q*(1-q^2)*Dq*P

    def op_v(P):
        return D*(v*(1-v^2)*P.derivative(v) +
                  ((2*n+4)-(2*n+6)*v^2)*P) - \
               exponent*v*(1-v^2)*Dv*P

    return op_p, op_q, op_v


def residual(n, i):
    M = matrix_at(n)
    S = p^2*q^2*(1-p^2)*(1-q^2)*v^2
    nxt = [S*D^2,
           S*2*(n+3)*v*D,
           S*(-(n+3)*v*D+2*(n+3)*(2*n+7)*v^2)]
    cur = [D^4,
           2*(n+2)*v*D^3,
           (-(n+2)*v*D+2*(n+2)*(2*n+5)*v^2)*D^2]
    # The scalar is determined after solving and should be independent of i.
    return sum(M[i,j]*nxt[j] for j in range(3)), cur[i]


def solve_fixed(n, i, degree, use_lambda=True):
    ops = operators(n)
    basis = [p^a*q^b*v^c for a in range(degree+1)
             for b in range(degree+1-a)
             for c in range(degree+1-a-b)
             if (a+b) % 2 == 0]
    rhs, cur = residual(n, i)
    cols = [op(m) for op in ops for m in basis]
    if use_lambda:
        cols.append(cur)
    exps = set(rhs.dict())
    for value in cols:
        exps.update(value.dict())
    exps = sorted(exps)
    index = {exp: row for row, exp in enumerate(exps)}
    entries = {}
    for column, value in enumerate(cols):
        for exp, coefficient in value.dict().items():
            entries[index[exp], column] = coefficient
    A = matrix(QQ, len(exps), len(cols), entries, sparse=True)
    b = vector(QQ, [rhs.dict().get(exp, 0) for exp in exps])
    try:
        solution = A.solve_right(b)
    except ValueError:
        return None, (A.nrows(), A.ncols())
    return solution, (A.nrows(), A.ncols())


if __name__ == "__main__":
    for degree in [15, 16, 18, 20]:
        print("degree", degree, flush=True)
        for i in range(1):
            solution, shape = solve_fixed(2, i, degree)
            if solution is None:
                print(i, False, shape, flush=True)
            else:
                print(i, True, shape, "lambda", factor(solution[-1]),
                      "nonzero", sum(x != 0 for x in solution), flush=True)
