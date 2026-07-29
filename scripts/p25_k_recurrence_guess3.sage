#!/usr/bin/env sage
"""
P2.5 follow-up: operator relations between Lf (annihilator of f_q) and
Lg (annihilator of g_p), and the geometric rate of eps_k = g_p/f_q - G.
"""
import sys, time
from ore_algebra import OreAlgebra, guess

KMAX = 400

def M_entries(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

t0 = time.time()
E = [vector(QQ,[1,0,0]), vector(QQ,[0,1,0]), vector(QQ,[0,0,1])]
first_col = [[QQ(E[j][0]) for j in range(3)]]
cur = [vector(v) for v in E]
for N in range(KMAX):
    M = matrix(QQ, M_entries(N)) / delta_H(N)
    cur = [v*M for v in cur]
    first_col.append([cur[j][0] for j in range(3)])
q = vector(QQ,[33750,-36000,9000]); p = vector(QQ,[30921,-32972,8240])
Q_vals = [sum(q[j]*first_col[N][j] for j in range(3)) for N in range(KMAX+1)]
P_vals = [sum(p[j]*first_col[N][j] for j in range(3)) for N in range(KMAX+1)]

def Bsummand(N,k):
    return 2**k*binomial(2*k,k)*binomial(N,k)*binomial(N+k,k)
def decompose(vals):
    c = []
    for K in range(len(vals)):
        rhs = vals[K]
        for k in range(K):
            rhs -= c[k]*Bsummand(K,k)
        c.append(rhs/Bsummand(K,K))
    return c
f_q = decompose(Q_vals); g_p = decompose(P_vals)
print("data done (%.1fs)" % (time.time()-t0)); sys.stdout.flush()

R = QQ['k']; A = OreAlgebra(R, 'Sk')
Lf = guess(f_q, A); Lg = guess(g_p, A)
print("Lf: order %d degree %d;  Lg: order %d degree %d"
      % (Lf.order(), Lf.degree(), Lg.order(), Lg.degree()))

def apply_op(L, seq):
    r = L.order(); cs = [L[i] for i in range(r+1)]
    for k in range(len(seq)-r):
        if sum(cs[i](k)*seq[k+i] for i in range(r+1)) != 0:
            return k
    return None

print("g_p under Lf: first failure at k =", apply_op(Lf, g_p))
print("f_q under Lg: first failure at k =", apply_op(Lg, f_q))
sys.stdout.flush()

t1 = time.time()
Lcm = Lf.lclm(Lg)
print("LCLM(Lf,Lg): order %d, degree %d  (%.1fs)"
      % (Lcm.order(), Lcm.degree(), time.time()-t1))
print("  -> common operator annihilating BOTH f and g has order", Lcm.order())
badf = apply_op(Lcm, f_q); badg = apply_op(Lcm, g_p)
print("  LCLM verification: f %s, g %s"
      % ("PASS" if badf is None else "FAIL@%d" % badf,
         "PASS" if badg is None else "FAIL@%d" % badg))
t1 = time.time()
try:
    Gcd = Lf.gcrd(Lg)
    print("GCRD(Lf,Lg): order %d  (%.1fs)" % (Gcd.order(), time.time()-t1))
except Exception as e:
    print("GCRD failed:", e)
sys.stdout.flush()

# eps_k decay rate, high precision
from mpmath import mp, mpf, catalan
mp.dps = 400
G = catalan
print("\neps_k = g_p(k)/f_q(k) - G, ratios eps_{k+1}/eps_k:")
eps = []
for k in range(0, 120):
    ratio = g_p[k]/f_q[k]
    e = mpf(ratio.numerator())/mpf(ratio.denominator()) - G
    eps.append(e)
import math
for k in [10, 20, 40, 60, 80, 100, 110, 118]:
    r = eps[k+1]/eps[k]
    print("  k=%3d: eps=%.3e  ratio=%.12f" % (k, float(eps[k]), float(r)))
print("  -1/8 = -0.125;  17-12*sqrt(2) = %.12f" % (17-12*math.sqrt(2)))
print("total %.1fs" % (time.time()-t0))
