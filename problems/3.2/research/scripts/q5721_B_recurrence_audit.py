#!/usr/bin/env python3
"""Exact reconstruction/audit for Q5721.

Uses exact integer arithmetic plus SymPy rational linear algebra.  It reconstructs
p1,p2 from the supplied p0,p3, verifies the recurrence independently, prints
factorizations and leading characteristic data, and audits the finite-difference
exceptional determinants.
"""
from math import comb
from sympy import Matrix, Poly, QQ, Symbol, expand, factor, linsolve, simplify

M = Symbol('M')
r = Symbol('r')


def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def B(m, s):
    if s < 0 or s > 2*m:
        return 0
    return sum(C(m,k)**2 * C(m+k,s)**2 for k in range(m+1))


def Q(m, s):
    return (-6*m**3 + 24*m**2*s - 27*m*s**2 + 9*s**3
            + 29*m**2 - 57*m*s + 27*s**2 - 26*m + 25*s + 7)

Qsym = Q(M,r)
p0 = -(2*M-r-1)*(r+1)*(2*M-r)**3*Q(M,r+1)
p3 = 2*(2*r-2*M+5)*(r-M+2)*(r+2)*(r+3)**2*Qsym

# The user reports total degree seven for p1,p2.
mons = [M**i*r**j for d in range(8) for i in range(d+1) for j in [d-i]]
assert len(mons) == 36
rows=[]
rhs=[]
# A deterministic overdetermined exact system.
for mv in range(5, 23):
    for rv in range(0, min(mv-2, 10)):
        vals = [B(mv,rv+i) for i in range(4)]
        if vals[1] == 0 and vals[2] == 0:
            continue
        ev = [int(mon.subs({M:mv,r:rv})) for mon in mons]
        rows.append([vals[1]*x for x in ev] + [vals[2]*x for x in ev])
        rhs.append(-int(p0.subs({M:mv,r:rv}))*vals[0]
                   -int(p3.subs({M:mv,r:rv}))*vals[3])

A=Matrix(rows)
b=Matrix(rhs)
print('linear_system_shape', A.shape, 'rank', A.rank())
solset=linsolve((A,b))
sols=list(solset)
assert len(sols)==1
sol=sols[0]
assert not any(x.free_symbols for x in sol)
p1=expand(sum(sol[i]*mons[i] for i in range(36)))
p2=expand(sum(sol[36+i]*mons[i] for i in range(36)))

print('p0 =', factor(p0))
print('p1 =', factor(p1))
print('p2 =', factor(p2))
print('p3 =', factor(p3))
print('p1_expanded =', p1)
print('p2_expanded =', p2)
print('degrees =', [Poly(x,M,r).total_degree() for x in (p0,p1,p2,p3)])

# Independent exact verification outside the reconstruction grid.
checks=0
for mv in range(1, 41):
    for rv in range(0, 2*mv-2):
        lhs=sum(int(p.subs({M:mv,r:rv}))*B(mv,rv+i)
                for i,p in enumerate((p0,p1,p2,p3)))
        assert lhs==0, (mv,rv,lhs)
        checks += 1
print('recurrence_checks', checks)

# Hypergeometric/Jacobi kernel coefficient recurrence check.
# F_r(z)=sum_k binom(M,k)binom(M+k,r)z^k.
def Fcoeff(m,s):
    return [C(m,k)*C(m+k,s) for k in range(m+1)]
def add_scaled(out, vec, scale, shift=0):
    for i,x in enumerate(vec):
        j=i+shift
        if 0 <= j < len(out): out[j]+=scale*x
for mv in range(1,15):
    for rv in range(0,2*mv):
        # (r+1)(1+z)F_{r+1}+(r-M)(1+2z)F_r+z(r-2M-1)F_{r-1}=0
        out=[0]*(mv+2)
        add_scaled(out,Fcoeff(mv,rv+1),rv+1,0)
        add_scaled(out,Fcoeff(mv,rv+1),rv+1,1)
        add_scaled(out,Fcoeff(mv,rv),rv-mv,0)
        add_scaled(out,Fcoeff(mv,rv),2*(rv-mv),1)
        add_scaled(out,Fcoeff(mv,rv-1),rv-2*mv-1,1)
        assert all(x==0 for x in out),(mv,rv,out)
print('kernel_contiguity_checks PASS')

# Leading characteristic at r=infinity with M fixed.
lead=[]
for p in (p0,p1,p2,p3):
    P=Poly(p,r)
    lead.append(P.LC())
print('r_degrees', [Poly(p,r).degree() for p in (p0,p1,p2,p3)])
print('r_leading_coefficients', lead)
lam=Symbol('lam')
chi=expand(sum(lead[i]*lam**i for i in range(4)))
print('r_infinity_characteristic =', factor(chi))
print('chi_at_1 =', factor(chi.subs(lam,1)))

# Shift-orbit valuation totals of p0/p3, represented by orbit labels.
# These are exact factor-divisor totals over Q(M)[r] modulo integer shifts.
orbits={
    'r-2M+Z': 4,
    'r+Z': 1-1-2,
    'r-M+Z': -1,
    'r-M+1/2+Z': -1,
    'Q(M,r+Z)': 0,
}
print('shift_orbit_totals_p0_over_p3', orbits)
print('mod3_obstruction', {k:v%3 for k,v in orbits.items()})
assert any(v%3 for v in orbits.values())

# Exterior-square scalar operator coefficients for X_r=u_r v_{r+1}-u_{r+1}v_r.
def sh(f,k): return expand(f.subs(r,r+k))
ext_coeffs=(
    expand(sh(p0,1)*p0),
    expand(-sh(p0,1)*p2),
    expand(sh(p1,1)*p3),
    expand(-sh(p3,1)*p3),
)
print('exterior_trailing =', factor(ext_coeffs[0]))
print('exterior_leading =', factor(ext_coeffs[3]))

# Finite-difference polynomial-segment determinants.
# If Delta^j B vanishes on at least 3 consecutive positions, B agrees with a
# polynomial of degree <j on j+3 points.  Substitution in j recurrence rows
# gives this determinant.
def fd_det(j):
    rows=[]
    for u in range(j):
        row=[]
        for h in range(j):
            e=0
            for i,p in enumerate((p0,p1,p2,p3)):
                # Newton basis binom(n,h), n=u+i.
                e += p.subs(r,r+u) * (comb(u+i,h) if h <= u+i else 0)
            row.append(expand(e))
        rows.append(row)
    return expand(Matrix(rows).det(method='domain-ge'))

for j in range(1,9):
    D=fd_det(j)
    assert D != 0
    Pd=Poly(D,M,r)
    print('fd_det',j,'total_degree',Pd.total_degree(),'terms',len(Pd.terms()))
    if j <= 3:
        print('fd_det_factor',j,'=',factor(D))
print('finite_difference_determinants_nonzero_through_8 PASS')

# Direct recurrence singular factor list.
print('singular_trailing_factors =', [
    '2*M-r-1','r+1','2*M-r (multiplicity 3)','Q(M,r+1)'])
print('singular_leading_factors =', [
    '2*r-2*M+5','r-M+2','r+2','r+3 (multiplicity 2)','Q(M,r)'])
