from sage.all import *

p = 17
N = 18
Fp = GF(p)
R = PowerSeriesRing(Fp, 'q', default_prec=N+1)
q = R.gen()

# eta products with q-powers already separated.
def prod_eta(d, exponent):
    z = R(1)
    for n in range(1, N//d + 2):
        z *= (1-q^(d*n))^exponent
        z = z.add_bigoh(N+1)
    return z

tq = (q * prod_eta(1,12)*prod_eta(6,12) / (prod_eta(2,12)*prod_eta(3,12))).add_bigoh(N+1)
Eq = (prod_eta(2,7)*prod_eta(3,7) / (prod_eta(1,5)*prod_eta(6,5))).add_bigoh(N+1)
assert tq[1] == 1 and Eq[0] == 1

# Revert t(q) recursively: q(t)=t+...
T = PowerSeriesRing(Fp, 't', default_prec=N+1)
t = T.gen()
tcoeff = [Fp(tq[n]) for n in range(N+1)]
qoft = T(0)
for n in range(1,N+1):
    # coefficient of t^n in t(q(t)) is affine with coefficient 1 in q_n.
    base = qoft
    comp0 = sum(tcoeff[k]*(base^k) for k in range(1,N+1)).add_bigoh(N+1)
    target = Fp(1) if n == 1 else Fp(0)
    need = target - Fp(comp0[n])
    qoft += need*t^n
    comp = sum(tcoeff[k]*(qoft^k) for k in range(1,N+1)).add_bigoh(N+1)
    assert comp[n] == target
assert (sum(tcoeff[k]*(qoft^k) for k in range(1,N+1)).add_bigoh(N+1) - t).add_bigoh(N+1) == 0

# E4(d tau) coefficients.
def sigma3(n):
    return sum(ZZ(d)^3 for d in divisors(n))

def E4d(d):
    a = [Fp(0)]*(N+1)
    a[0] = Fp(1)
    for n in range(1,N+1):
        if n % d == 0:
            a[n] = Fp(240*sigma3(n//d))
    return R(a).add_bigoh(N+1)

# K_d(t(q))=E(q)*U_d(q), delta^3 U_d = E4(d tau)-1.
def Kd_t(d):
    Md = E4d(d)
    U = R(0)
    for n in range(1,N+1):
        U += Md[n]/Fp(n^3)*q^n
    Kq = (Eq*U).add_bigoh(N+1)
    Kt = T(0)
    for n in range(1,N+1):
        Kt += Fp(Kq[n])*(qoft^n)
    return Kt.add_bigoh(N+1)

Ks = [Kd_t(d) for d in [1,2,3,6]]
for d,K in zip([1,2,3,6],Ks):
    print('K',d,'k1,k3,k13',int(K[1]),int(K[3]),int(K[13]))

# Canonical coefficients (-3,4,-9,108)/20.
cc = vector(Fp,[Fp(-3)/Fp(20),Fp(4)/Fp(20),Fp(-9)/Fp(20),Fp(108)/Fp(20)])
for r in [1,3,13]:
    val=sum(cc[i]*Ks[i][r] for i in range(4))
    print('CANON',r,int(val))
assert sum(cc[i]*Ks[i][1] for i in range(4)) == Fp(-36)
assert sum(cc[i]*Ks[i][13] for i in range(4)) == 0
assert sum(cc[i]*Ks[i][3] for i in range(4)) != 0

# Search all forms with same initial kappa_1=-36 and zeros at r=3,13.
A = matrix(Fp,3,4,lambda i,j: Ks[j][[1,3,13][i]])
y = vector(Fp,[Fp(-36),Fp(0),Fp(0)])
print('RANK',A.rank(),'AUGRANK',A.augment(y).rank())
assert A.rank() == A.augment(y).rank()
sol = A.solve_right(y)
print('FLEX_SOL',tuple(int(x) for x in sol),'constant',int(sum(sol)))
print('CANON_COEFF',tuple(int(x) for x in cc),'constant',int(sum(cc)))
assert sol != cc
print('Q7684_EISENSTEIN_FLEX=PASS')
