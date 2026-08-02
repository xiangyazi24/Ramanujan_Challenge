from sage.all import *

exec(open("p25_dual_flux_search.sage").read().split("for nn in")[0])


def solve(nn, row, degree):
    R = PolynomialRing(QQ, names=("p", "q", "r", "t"))
    p,q,r,t=R.gens(); K=R.fraction_field()
    h=K(r/(p*q))
    A=K(t+1/t+2/(p*q))
    variables=[p,q,r,t]
    logs=[
      2/p-2*nn*p/(1-p^2)+2*r/(p^2*q),
      1/q-2*(nn+1)*q/(1-q^2)+2*r/(p*q^2),
      (2*nn+3)/r-A,
      -1/t-r*(1-1/t^2),
    ]
    ops=[lambda P,x=x,ell=ell: K(P.derivative(x)+ell*P)
         for x,ell in zip(variables,logs)]
    S=(1-p^2)*(1-q^2)*r^2
    cur=[K(1),h,h^2-h/2]
    M=matrix_at(nn)
    lam=2*(nn+1)*(nn+2)^3*(nn+3)^2*(2*nn+5)*(2*nn+7)^2
    target=sum(M[row,j]*S*cur[j] for j in range(3))-lam*cur[row]
    common=p^2*q^2*r*t^2*(1-p^2)*(1-q^2)
    target=R(K(common)*target)
    basis=[p^a*q^b*r^c*t^d
           for a in range(degree+1)
           for b in range(degree+1-a)
           for c in range(degree+1-a-b)
           for d in range(degree+1-a-b-c)]
    cols=[R(K(common)*op(mon)) for op in ops for mon in basis]
    exponents=set(target.dict())
    for col in cols: exponents.update(col.dict())
    exponents=sorted(exponents); index={e:i for i,e in enumerate(exponents)}
    entries={}
    for j,col in enumerate(cols):
      for e,c in col.dict().items():
        if c: entries[index[e],j]=c
    mat=matrix(QQ,len(exponents),len(cols),entries,sparse=True)
    rhs=vector(QQ,[target.dict().get(e,0) for e in exponents])
    try: answer=mat.solve_right(rhs)
    except ValueError: answer=None
    return answer,(mat.nrows(),mat.ncols()),basis


for degree in range(2,11):
    answer,shape,basis=solve(2,0,degree)
    print("degree",degree,"found",answer is not None,"shape",shape,
          "support",None if answer is None else sum(c != 0 for c in answer),
          flush=True)
    if answer is not None:
      chunk=len(basis)
      for op in range(4):
        support=[(basis[i],answer[op*chunk+i]) for i in range(chunk)
                 if answer[op*chunk+i]]
        print("operator",op,"support",support,flush=True)
      break
