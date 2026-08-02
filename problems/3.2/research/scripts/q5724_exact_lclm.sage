#!/usr/bin/env sage
"""Exact Q5724 minimizer for the 14-ray first-cell shell.

Run from a checkout of current main, where q32_first_cell_ray_telescopers.sage
is present.  It derives the fourteen exact telescopers, forms the exact Ore
LCLM, reports all intermediate orders and factored endpoint coefficients, and
then proves generic minimality by specialized finite-field guessing at held-out
moments.  The weighted ray sum is evaluated independently from the exact
coefficient formula.
"""
from itertools import combinations
from math import comb
from ore_algebra import *
from ore_algebra import nullspace
from ore_algebra.ore_algebra import OreAlgebra_generic
from ore_algebra.guessing import guess_raw

# Same QQ-safe patches as the source derivation.
def _associated_commutative_algebra(self):
    try:
        return self._commutative_ring
    except AttributeError:
        self._commutative_ring = PolynomialRing(self.base_ring(), self.variable_names())
        return self._commutative_ring
OreAlgebra_generic.associated_commutative_algebra = _associated_commutative_algebra
_original_kronecker = nullspace.kronecker
def _qq_safe_kronecker(subsolver, presolver=None):
    return nullspace.clear(_original_kronecker(subsolver, presolver))
nullspace.kronecker = _qq_safe_kronecker

R.<M,r,t> = QQ[]
OA.<Sr,St> = OreAlgebra(R)
d=M-r
RAY_CLASSES=(
 ((-1,-1,-1),1),((-1,-1,0),2),((-1,-1,1),2),
 ((-1,0,0),1),((-1,0,1),2),((-1,1,1),1),
 ((0,-1,-1),1),((0,-1,0),2),((0,-1,1),2),
 ((0,0,1),2),((0,1,1),1),((1,0,0),1),
 ((1,0,1),2),((1,1,1),1))

def shifted_binomial_ratio(upper,lower,shift):
    if shift==1:return (upper-lower)/(lower+1)
    if shift==-1:return lower/(upper-lower+1)
    return R.one()

def ray_operator(point):
    u,v,w=point;kx=t-u*d;upper=2*M-t;ky=M-v*d;kz=M-w*d
    tr=(M-t)/(t+1)*(M-kx)/(kx+1)*(upper-ky)/upper*(upper-kz)/upper
    rr=shifted_binomial_ratio(M,kx,u)*shifted_binomial_ratio(upper,ky,v)*shifted_binomial_ratio(upper,kz,w)
    ideal=OA.ideal([tr.denominator()*St-tr.numerator(),rr.denominator()*Sr-rr.numerator()])
    tel,cert=ideal.ct(St-1,certificates=True,early_termination=True,iteration_limit=16)
    assert len(tel)==1
    return tel[0].primitive_part()

ops=[]
for point,mult in RAY_CLASSES:
    L=ray_operator(point)
    ops.append(L)
    print('RAY',point,'MULT',mult,'ORDER',L.order())
    print(' TRAIL',factor(L[0]))
    print(' LEAD',factor(L[L.order()]))

# Greedy pair order: at each stage choose the next operator producing the
# smallest coefficient degree, while recording exact gcrd/lclm orders.
remaining=list(range(1,len(ops)))
L=ops[0]
chosen=[0]
while remaining:
    trials=[]
    for j in remaining:
        g=L.gcrd(ops[j])
        C=L.lclm(ops[j]).primitive_part()
        deg=max(c.degree(r) for c in C)
        trials.append((C.order(),deg,j,g.order(),C))
    _,_,j,go,C=min(trials,key=lambda x:(x[0],x[1],x[2]))
    L=C;chosen.append(j);remaining.remove(j)
    print('LCLM_STEP','chosen',chosen,'new',j,'gcrd_order',go,'order',L.order(),
          'max_r_degree',max(c.degree(r) for c in L))

assert all((L % A)==0 for A in ops)
print('GENERIC_LCLM_ORDER',L.order())
print('GENERIC_TRAILING',factor(L[0]))
print('GENERIC_LEADING',factor(L[L.order()]))
print('GENERIC_CONTENT',factor(L.content()))
print('GENERIC_OPERATOR',L)

# Exact weighted sum values.
def C(n,k):return comb(n,k) if 0<=k<=n else 0
def ray_value(moment,residue,point,mod=None):
    u,v,w=point;node=moment-residue
    z=sum(C(moment,i)*C(moment,i-node*u)*C(2*moment-i,moment-node*v)*C(2*moment-i,moment-node*w) for i in range(moment+1))
    return z if mod is None else z%mod
def F_value(moment,residue,mod=None):
    z=sum(mult*ray_value(moment,residue,point,mod) for point,mult in RAY_CLASSES)
    return z if mod is None else z%mod

# Specialize the exact generic LCLM modulo p, then use its r-degree as a
# rigorous search bound.  A lower-order operator for the distinguished sum
# must appear in guess_raw with these degree bounds after desingularizing the
# specialized LCLM; we also run increasing degree up to the exact LCLM degree.
for p,moment in ((1000003,997),(1000033,1009),(1000037,1013)):
    K=GF(p);Pr.<n>=K[];Ar.<Sn>=OreAlgebra(Pr)
    coeffs=[]
    for c in L:
        q=Pr(c(M=moment,r=n,t=0))
        coeffs.append(q)
    Lm=Ar(coeffs).primitive_part()
    assert Lm.order()==L.order()
    data=[K(F_value(moment,j,p)) for j in range((moment-1)//2+1)]
    maxdeg=max(c.degree() for c in Lm)
    found=None
    for order in range(1,Lm.order()+1):
        # Only degrees with enough held-out equations are attempted.
        degree_cap=min(maxdeg,(len(data)-order-32)//(order+1)-1)
        if degree_cap<0:continue
        basis=guess_raw(data,Ar,order=order,degree=degree_cap,lift=K,cut=None,ensure=16)
        nonzero=[A.primitive_part() for A in basis if A]
        if nonzero:
            found=min(nonzero,key=lambda A:(A.order(),A.degree()))
            break
    assert found is not None
    print('SPECIALIZED_MINIMAL',p,moment,'order',found.order(),'degree',found.degree(),
          'lclm_degree',maxdeg)
    print(' SPECIAL_TRAIL',factor(found[0]))
    print(' SPECIAL_LEAD',factor(found[found.order()]))
    # Held-out exact application to all available positions.
    assert all(x==0 for x in found(data))
    assert found.order()==L.order(), ('generic order drops',p,moment,found.order())

print('PASS: exact 14-ray LCLM and generic minimality')
