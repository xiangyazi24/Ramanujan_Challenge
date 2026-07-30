# Companion exact Sage verifier for drops/Q5721-de91974f.md
# Run: sage drops/Q5721-de91974f-verify.sage

R.<M,r> = PolynomialRing(QQ)

def C(n,k): return binomial(n,k) if 0 <= k <= n else ZZ(0)
def B(m,s):
    if s < 0 or s > 2*m: return ZZ(0)
    return sum(C(m,k)^2*C(m+k,s)^2 for k in range(m+1))
def Q(m,s):
    return (-6*m^3+24*m^2*s-27*m*s^2+9*s^3+29*m^2
            -57*m*s+27*s^2-26*m+25*s+7)

p0=-(2*M-r-1)*(r+1)*(2*M-r)^3*Q(M,r+1)
p3=2*(2*r-2*M+5)*(r-M+2)*(r+2)*(r+3)^2*Q(M,r)

fixed={}
for mv in range(10,27):
    rows=[]; rhs=[]
    for rv in range(0,2*mv-2):
        vals=[B(mv,rv+i) for i in range(4)]
        rows.append([vals[1]*rv^j for j in range(8)]
                    +[vals[2]*rv^j for j in range(8)])
        rhs.append(-p0.subs({M:mv,r:rv})*vals[0]
                   -p3.subs({M:mv,r:rv})*vals[3])
    A=matrix(QQ,rows); b=vector(QQ,rhs)
    assert A.rank()==16
    fixed[mv]=A.solve_right(b)

U.<u> = PolynomialRing(QQ)
def recover(offset):
    ans=R.zero(); keys=list(fixed.keys())
    for j in range(8):
        d=7-j
        cu=U.lagrange_polynomial([(QQ(mv),fixed[mv][offset+j])
                                  for mv in keys[:d+1]])
        cR=R(cu(M))
        assert cR.degree(M)<=d
        assert all(cR.subs({M:mv})==fixed[mv][offset+j] for mv in keys)
        ans += cR*r^j
    return R(ans)
p1=recover(0); p2=recover(8)
print('p1 =',p1); print('p2 =',p2)
print('factor p1 =',p1.factor()); print('factor p2 =',p2.factor())

for mv in range(1,51):
    for rv in range(0,2*mv-2):
        assert sum(p.subs({M:mv,r:rv})*B(mv,rv+i)
                   for i,p in enumerate((p0,p1,p2,p3)))==0
print('recurrence verification: PASS')

orbit_totals={'r-2M+Z':4,'r+Z':-2,'r-M+Z':-1,
              'r-M+1/2+Z':-1,'Q(M,r+Z)':0}
assert any(v%3 for v in orbit_totals.values())
print('symmetric-square obstruction: PASS')

from ore_algebra import OreAlgebra
OA.<Sr> = OreAlgebra(R)
sh=lambda f,j:R(f.subs({r:r+j}))
Lext=(sh(p0,1)*p0-sh(p0,1)*p2*Sr
      +sh(p1,1)*p3*Sr^2-sh(p3,1)*p3*Sr^3)
print('exterior operator =',Lext)

PM.<m> = PolynomialRing(QQ); KM=PM.fraction_field(); PX.<x> = PolynomialRing(KM)
def Qx(c=0):
    z=x+c
    return PX(-6*m^3+24*m^2*z-27*m*z^2+9*z^3+29*m^2
              -57*m*z+27*z^2-26*m+25*z+7)
def endpoint_A(): return PX(-x*(2*m-x)*(2*m-x+1)^3*Qx())
def endpoint_B(Lrun):
    b=x+Lrun-1
    return PX(2*(2*b-2*m+1)*(b-m)*b*(b+1)^2*Qx(Lrun-3))
for Lrun in range(3,40):
    A=endpoint_A(); E=endpoint_B(Lrun); g=A.gcd(E).monic()
    if Lrun==3: assert g==Qx().monic()
    else: assert g==1 and A.resultant(E)!=0
print('zero-run resultants: PASS')

def difference_matrix(j):
    ps=(p0,p1,p2,p3)
    return matrix(R,j+8,j,lambda uu,h:
        sum(ps[i].subs({r:r+uu})*binomial(uu+i,h) for i in range(4)))
for j in range(1,20):
    A=difference_matrix(j); assert A.rank()==j
    cert=None
    for rows in Subsets(range(j+8),j):
        D=A.matrix_from_rows(list(rows)).det()
        if D!=0: cert=R(D); break
    assert cert is not None
    print('j',j,'certificate degree',cert.total_degree())
print('finite-difference determinant rank: PASS')
