from sage.all import *

exec(open("p25_dual_flux_search.sage").read().split("for nn in")[0])


def solve_shift(nn, row, shift, degree):
    R = PolynomialRing(QQ, names=("p", "q", "v"))
    p, q, v = R.gens()
    K = R.fraction_field()
    D = p*q*(1+v^2)+2*v
    exponent = 2*nn+8+shift
    logarithmic = [
        (2*nn+6)/p - 2*nn*p/(1-p^2),
        (2*nn+5)/q - 2*(nn+1)*q/(1-q^2),
        (2*nn+3)/v,
    ]
    variables = [p,q,v]
    ops = [lambda P, x=x, ell=ell:
             K(D*P.derivative(x)+D*ell*P-exponent*D.derivative(x)*P)
           for x,ell in zip(variables,logarithmic)]
    left,current=residual_polynomials(nn,row,R)
    target=(left-current)*D^shift
    common=p*q*v*(1-p^2)*(1-q^2)
    basis=[p^a*q^b*v^c for a in range(degree+1)
           for b in range(degree+1-a)
           for c in range(degree+1-a-b) if (a+b)%2==0]
    cols=[R(common*op(mon)) for op in ops for mon in basis]
    target=R(common*target)
    exponents=set(target.dict())
    for col in cols: exponents.update(col.dict())
    exponents=sorted(exponents); index={e:i for i,e in enumerate(exponents)}
    entries={}
    for j,col in enumerate(cols):
        for e,c in col.dict().items():
            if c: entries[index[e],j]=c
    A=matrix(QQ,len(exponents),len(cols),entries,sparse=True)
    b=vector(QQ,[target.dict().get(e,0) for e in exponents])
    try: answer=A.solve_right(b)
    except ValueError: answer=None
    return answer,(A.nrows(),A.ncols()),basis


def residual_polynomials(nn,row,R):
    p,q,v=R.gens(); D=p*q*(1+v^2)+2*v
    M=matrix_at(nn)
    S=p^2*q^2*(1-p^2)*(1-q^2)*v^2
    nxt=[S*D^2,S*2*(nn+3)*v*D,
         S*(-(nn+3)*v*D+2*(nn+3)*(2*nn+7)*v^2)]
    cur=[D^4,2*(nn+2)*v*D^3,
         (-(nn+2)*v*D+2*(nn+2)*(2*nn+5)*v^2)*D^2]
    lam=(nn+1)*(nn+2)^2*(nn+3)^2*(2*nn+7)^2
    return sum(M[row,j]*nxt[j] for j in range(3)),lam*cur[row]


for shift in range(1,7):
    for degree in [8,10,12,14,16,18]:
        answer,shape,basis=solve_shift(2,0,shift,degree)
        print("shift",shift,"degree",degree,"found",answer is not None,
              "shape",shape,"support",None if answer is None else
              sum(c != 0 for c in answer),flush=True)
        if answer is not None:
            chunk=len(basis)
            for op in range(3):
                support=[(basis[i],answer[op*chunk+i]) for i in range(chunk)
                         if answer[op*chunk+i]]
                print("operator",op,"support",support,flush=True)
            raise SystemExit
