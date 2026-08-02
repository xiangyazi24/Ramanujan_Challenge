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


def solve(nn, row, degree):
    R = PolynomialRing(QQ, names=("p", "q", "v"))
    p, q, v = R.gens()
    K = R.fraction_field()
    D = p*q*(1+v^2)+2*v
    exponent = 2*nn+8

    log_derivatives = [
        (2*nn+6)/p - 2*nn*p/(1-p^2),
        (2*nn+5)/q - 2*(nn+1)*q/(1-q^2),
        (2*nn+3)/v,
    ]
    variables = [p, q, v]
    ops = [
        lambda P, variable=variable, logarithmic=logarithmic:
          K(D*P.derivative(variable) + D*logarithmic*P -
            exponent*D.derivative(variable)*P)
        for variable, logarithmic in zip(variables, log_derivatives)
    ]

    M = matrix_at(nn)
    S = p^2*q^2*(1-p^2)*(1-q^2)*v^2
    nxt = [S*D^2, S*2*(nn+3)*v*D,
           S*(-(nn+3)*v*D+2*(nn+3)*(2*nn+7)*v^2)]
    cur = [D^4, 2*(nn+2)*v*D^3,
           (-(nn+2)*v*D+2*(nn+2)*(2*nn+5)*v^2)*D^2]
    lam = (nn+1)*(nn+2)^2*(nn+3)^2*(2*nn+7)^2
    rhs = sum(M[row,j]*nxt[j] for j in range(3))-lam*cur[row]

    # Clear only the logarithmic-derivative denominators.  The parity of
    # p+q is preserved, so the odd sector cannot contribute.
    common = p*q*v*(1-p^2)*(1-q^2)
    basis = [p^a*q^b*v^c for a in range(degree+1)
             for b in range(degree+1-a)
             for c in range(degree+1-a-b)
             if (a+b) % 2 == 0]
    cols = [R(common*operator(monomial))
            for operator in ops for monomial in basis]
    target = R(common*rhs)
    exponents = set(target.dict())
    for col in cols:
        exponents.update(col.dict())
    exponents = sorted(exponents)
    index = {exp:i for i,exp in enumerate(exponents)}
    entries = {}
    for j,col in enumerate(cols):
        for exp,coefficient in col.dict().items():
            if coefficient:
                entries[index[exp],j] = coefficient
    A = matrix(QQ,len(exponents),len(cols),entries,sparse=True)
    bvec = vector(QQ,[target.dict().get(exp,0) for exp in exponents])
    try:
        solution=A.solve_right(bvec)
    except ValueError:
        return None,(A.nrows(),A.ncols()),basis
    return solution,(A.nrows(),A.ncols()),basis


for nn in [2, 1, 0]:
    for degree in [4, 6, 8, 10, 12, 14]:
        answer,shape,basis=solve(nn,0,degree)
        print("n",nn,"degree",degree,"found",answer is not None,
              "shape",shape,
              "support",None if answer is None else sum(x != 0 for x in answer),
              flush=True)
        if answer is not None:
            chunk=len(basis)
            for operator in range(3):
                support=[(basis[i],answer[operator*chunk+i])
                         for i in range(chunk) if answer[operator*chunk+i]]
                print("operator",operator,"support",support,flush=True)
            raise SystemExit
