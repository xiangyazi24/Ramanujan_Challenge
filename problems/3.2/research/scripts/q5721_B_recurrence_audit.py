#!/usr/bin/env python3
"""Exact reconstruction/audit for Q5721.

The reconstruction uses small exact linear systems at fixed M and polynomial
interpolation in M.  It then verifies the result well outside the interpolation
range and audits the structural identities used in the answer.
"""
from math import comb
from sympy import Matrix, Poly, Symbol, expand, factor, interpolate, linsolve

M = Symbol('M'); r = Symbol('r')

def C(n,k): return comb(n,k) if 0 <= k <= n else 0

def B(m,s):
    if s < 0 or s > 2*m: return 0
    return sum(C(m,k)**2*C(m+k,s)**2 for k in range(m+1))

def Q(m,s):
    return (-6*m**3+24*m**2*s-27*m*s**2+9*s**3+29*m**2
            -57*m*s+27*s**2-26*m+25*s+7)

p0=-(2*M-r-1)*(r+1)*(2*M-r)**3*Q(M,r+1)
p3=2*(2*r-2*M+5)*(r-M+2)*(r+2)*(r+3)**2*Q(M,r)

# At each fixed M, solve for the 16 coefficients of two degree-7 r-polynomials.
fixed={}
for mv in range(10,27):
    rows=[]; rhs=[]
    for rv in range(0,2*mv-2):
        vals=[B(mv,rv+i) for i in range(4)]
        rows.append([vals[1]*rv**j for j in range(8)]
                    +[vals[2]*rv**j for j in range(8)])
        rhs.append(-int(p0.subs({M:mv,r:rv}))*vals[0]
                   -int(p3.subs({M:mv,r:rv}))*vals[3])
    sol=list(linsolve((Matrix(rows),Matrix(rhs))))
    assert len(sol)==1 and not any(x.free_symbols for x in sol[0]), mv
    fixed[mv]=sol[0]
print('fixed_M_systems',len(fixed),'PASS')

# Total degree <=7 means [r^j] has M-degree <=7-j.
def recover(offset):
    out=0
    for j in range(8):
        deg=7-j
        pts=[(mv,fixed[mv][offset+j]) for mv in list(fixed)[:deg+1]]
        cj=expand(interpolate(pts,M))
        assert Poly(cj,M).degree() <= deg
        assert all(cj.subs(M,mv)==fixed[mv][offset+j] for mv in fixed)
        out += cj*r**j
    return expand(out)
p1=recover(0); p2=recover(8)

print('p0 =',factor(p0))
print('p1 =',factor(p1))
print('p2 =',factor(p2))
print('p3 =',factor(p3))
print('p1_expanded =',p1)
print('p2_expanded =',p2)
print('degrees =',[Poly(x,M,r).total_degree() for x in (p0,p1,p2,p3)])

checks=0
for mv in range(1,51):
    for rv in range(0,2*mv-2):
        lhs=sum(int(p.subs({M:mv,r:rv}))*B(mv,rv+i)
                for i,p in enumerate((p0,p1,p2,p3)))
        assert lhs==0,(mv,rv,lhs)
        checks+=1
print('recurrence_checks',checks)

# Exact coefficient kernel F_r(z).
def Fcoeff(m,s): return [C(m,k)*C(m+k,s) for k in range(m+1)]
def add(out,v,c,sh=0):
    for i,x in enumerate(v):
        if 0 <= i+sh < len(out): out[i+sh]+=c*x
for mv in range(1,18):
    for rv in range(0,2*mv):
        out=[0]*(mv+2)
        add(out,Fcoeff(mv,rv+1),rv+1); add(out,Fcoeff(mv,rv+1),rv+1,1)
        add(out,Fcoeff(mv,rv),rv-mv); add(out,Fcoeff(mv,rv),2*(rv-mv),1)
        add(out,Fcoeff(mv,rv-1),rv-2*mv-1,1)
        assert all(x==0 for x in out),(mv,rv)
print('kernel_contiguity_checks PASS')

lead=[Poly(p,r).LC() for p in (p0,p1,p2,p3)]
lam=Symbol('lam'); chi=expand(sum(lead[i]*lam**i for i in range(4)))
print('r_degrees',[Poly(p,r).degree() for p in (p0,p1,p2,p3)])
print('r_leading_coefficients',lead)
print('r_infinity_characteristic =',factor(chi))
print('chi_at_1 =',factor(chi.subs(lam,1)))

orbits={'r-2M+Z':4,'r+Z':-2,'r-M+Z':-1,
        'r-M+1/2+Z':-1,'Q(M,r+Z)':0}
print('shift_orbit_totals_p0_over_p3',orbits)
print('mod3_obstruction',{k:v%3 for k,v in orbits.items()})
assert any(v%3 for v in orbits.values())

def sh(f,k): return expand(f.subs(r,r+k))
ext_coeffs=(sh(p0,1)*p0,-sh(p0,1)*p2,sh(p1,1)*p3,-sh(p3,1)*p3)
print('exterior_trailing =',factor(ext_coeffs[0]))
print('exterior_leading =',factor(ext_coeffs[3]))

# First finite-difference determinant factors; larger j are checked numerically.
def fd_matrix(j,MM=M,rr=r):
    return Matrix([[expand(sum(p.subs({M:MM,r:rr+u})
                         *(comb(u+i,h) if h<=u+i else 0)
                         for i,p in enumerate((p0,p1,p2,p3))))
                    for h in range(j)] for u in range(j)])
for j in range(1,4):
    D=expand(fd_matrix(j).det(method='domain-ge'))
    assert D!=0
    print('fd_det',j,'degree',Poly(D,M,r).total_degree(),'factor =',factor(D))
for j in range(1,17):
    D=fd_matrix(j,17,3).det(method='domain-ge')
    assert D!=0,j
print('finite_difference_determinants_numeric_nonzero_through_16 PASS')

print('singular_trailing_factors',[
 '2*M-r-1','r+1','2*M-r (multiplicity 3)','Q(M,r+1)'])
print('singular_leading_factors',[
 '2*r-2*M+5','r-M+2','r+2','r+3 (multiplicity 2)','Q(M,r)'])
