#!/usr/bin/env python3
from fractions import Fraction as F
from math import comb
import sympy as sp

NMAX = 90


def Acoef(n):
    inner = sum(comb(n,k)**2 * comb(n+k,n) * comb(n+2*k,n) for k in range(n+1))
    return F(comb(2*n,n) * inner)


def PA(n):
    n=F(n)
    return F(1024)*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860)

def PB(n):
    n=F(n)
    return F(128)*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)

def PC(n):
    n=F(n)
    return F(16)*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)

def PD(n):
    n=F(n)
    return (n+3)**4*(n+4)**6*(946*n*n+4515*n+5399)


def qterms(N):
    q=[F(-215040420000), F(-167282265043404,905), F(-964185327658080,6071)]
    for n in range(2,N-1):
        q.append(PB(n)/PA(n)*q[n] - PC(n-1)/PA(n-1)*q[n-1] + PD(n-2)/PA(n-2)*q[n-2])
    return q[:N]


def sym(x): return sp.Rational(x.numerator,x.denominator)


def fit_rational_values(vals, maxdeg=12, holdout=12):
    # vals[n] = P(n)/Q(n), normalize Q leading/constant through homogeneous nullspace.
    N=len(vals)
    for dp in range(maxdeg+1):
        for dq in range(maxdeg+1):
            nv=dp+dq+2
            train=min(N-holdout, nv+4)
            if train < nv-1: continue
            rows=[]
            for n in range(train):
                v=sym(vals[n])
                rows.append([sp.Integer(n)**i for i in range(dp+1)] + [-v*sp.Integer(n)**j for j in range(dq+1)])
            ns=sp.Matrix(rows).nullspace()
            for vec in ns:
                P=sum(vec[i]*sp.Symbol('x')**i for i in range(dp+1))
                Q=sum(vec[dp+1+j]*sp.Symbol('x')**j for j in range(dq+1))
                if Q==0: continue
                ok=True
                for n,v0 in enumerate(vals):
                    if sp.expand(P.subs({'x':n})-sym(v0)*Q.subs({'x':n})) != 0:
                        ok=False; break
                if ok:
                    return sp.factor(P),sp.factor(Q)
    return None


def fit_hyper_ratio(vals,maxdeg=12,holdout=12):
    ratios=[vals[n+1]/vals[n] for n in range(len(vals)-1) if vals[n] != 0]
    return fit_rational_values(ratios,maxdeg,holdout)


def guess_prec(vals,max_order=6,max_degree=14,holdout=15):
    N=len(vals)
    x=sp.Symbol('n')
    for r in range(1,max_order+1):
        for d in range(max_degree+1):
            nv=(r+1)*(d+1)
            avail=N-r
            train=min(avail-holdout,nv+5)
            if train < nv-1: continue
            rows=[]
            for n in range(train):
                row=[]
                for j in range(r+1):
                    for k in range(d+1):
                        row.append(sp.Integer(n)**k*sym(vals[n+j]))
                rows.append(row)
            ns=sp.Matrix(rows).nullspace()
            for v in ns:
                ps=[]; idx=0
                for j in range(r+1):
                    ps.append(sp.factor(sum(v[idx+k]*x**k for k in range(d+1))))
                    idx+=d+1
                if all(p==0 for p in ps): continue
                ok=True
                for n in range(avail):
                    if sp.simplify(sum(ps[j].subs(x,n)*sym(vals[n+j]) for j in range(r+1))) != 0:
                        ok=False; break
                if ok: return r,d,ps
    return None


def search_shift_post(f,q,max_order=5,max_degree=18,holdout=15):
    N=min(len(f),len(q)); nvar=sp.Symbol('n')
    for r in range(max_order+1):
        for d in range(max_degree+1):
            nv=(r+1)*(d+1); avail=N-r
            train=min(avail-holdout,nv+5)
            if train<nv: continue
            rows=[]; rhs=[]
            for n in range(train):
                rows.append([sp.Integer(n)**k*sym(f[n+j]) for j in range(r+1) for k in range(d+1)])
                rhs.append(sym(q[n]))
            sol=sp.linsolve((sp.Matrix(rows),sp.Matrix(rhs)))
            if sol is sp.EmptySet or sol==sp.EmptySet: continue
            sols=list(sol)
            if not sols: continue
            vec=sols[0]
            params=set().union(*(e.free_symbols for e in vec))
            vec=tuple(sp.simplify(e.subs({p:0 for p in params})) for e in vec)
            ps=[]; idx=0
            for j in range(r+1):
                ps.append(sp.factor(sum(vec[idx+k]*nvar**k for k in range(d+1))))
                idx+=d+1
            ok=True
            for n in range(avail):
                if sp.simplify(sum(ps[j].subs(nvar,n)*sym(f[n+j]) for j in range(r+1))-sym(q[n]))!=0:
                    ok=False; break
            if ok: return r,d,ps
    return None


def search_gf_post(f,q,max_theta=5,max_zdeg=18,holdout=15):
    # q(z)=sum_i p_i(z) theta^i f(z)
    N=min(len(f),len(q)); z=sp.Symbol('z')
    for dth in range(max_theta+1):
        for dz in range(max_zdeg+1):
            nv=(dth+1)*(dz+1)
            train=min(N-holdout,nv+5)
            if train<nv: continue
            rows=[]; rhs=[]
            for n in range(train):
                row=[]
                for i in range(dth+1):
                    for j in range(dz+1):
                        row.append(sp.Integer(n-j)**i*sym(f[n-j]) if n>=j else 0)
                rows.append(row); rhs.append(sym(q[n]))
            sol=sp.linsolve((sp.Matrix(rows),sp.Matrix(rhs)))
            if sol is sp.EmptySet or sol==sp.EmptySet: continue
            sols=list(sol)
            if not sols: continue
            vec=sols[0]
            params=set().union(*(e.free_symbols for e in vec))
            vec=tuple(sp.simplify(e.subs({p:0 for p in params})) for e in vec)
            ps=[]; idx=0
            for i in range(dth+1):
                ps.append(sp.factor(sum(vec[idx+j]*z**j for j in range(dz+1))))
                idx+=dz+1
            ok=True
            for n in range(N):
                val=0
                for i in range(dth+1):
                    for j in range(dz+1):
                        if n>=j: val += vec[i*(dz+1)+j]*sp.Integer(n-j)**i*sym(f[n-j])
                if sp.simplify(val-sym(q[n])) != 0:
                    ok=False; break
            if ok: return dth,dz,ps
    return None


def main():
    q=qterms(NMAX)
    AA=[Acoef(n) for n in range(NMAX)]
    f4=[AA[n]/F(4)**n for n in range(NMAX)]
    f256=[AA[n]/F(256)**n for n in range(NMAX)]
    h4=[q[n]/f4[n] for n in range(NMAX)]
    h256=[q[n]/f256[n] for n in range(NMAX)]
    print('A[0:10]=',AA[:10])
    print('q[0:8]=',q[:8])
    print('h4=q*4^n/A first 16:')
    for n in range(16): print(n,h4[n],float(h4[n]),'ratio',float(h4[n]/h4[n-1]) if n else None)
    print('h256=q*256^n/A normalized by q0 first 16:')
    for n in range(16): print(n,h256[n]/h256[0],float(h256[n]/h256[0]),'ratio',float(h256[n]/h256[n-1]) if n else None)
    print('fit rational h256/q0:',fit_rational_values([x/h256[0] for x in h256],12,15))
    print('fit hyper ratio h4:',fit_hyper_ratio(h4,14,15))
    print('fit hyper ratio h256:',fit_hyper_ratio(h256,14,15))
    print('guess P-rec h256 order<=6 degree<=14:',guess_prec(h256,6,14,15))
    print('shift post f256 -> q order<=5 degree<=18:',search_shift_post(f256,q,5,18,15))
    print('GF post f256 -> q theta<=5 zdeg<=18:',search_gf_post(f256,q,5,18,15))

if __name__=='__main__': main()
