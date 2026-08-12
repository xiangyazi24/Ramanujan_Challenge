from sage.all import *

# Q7674: exact Hankel/Padé audit for the Q7621 transverse sequence.
# Everything below is exact except the explicitly labeled logarithmic summaries.

N = 180
NTEST = 150
MMAX = 10


def P(n):
    n = ZZ(n)
    return 34*n^3 + 51*n^2 + 27*n + 5

# Apéry b_n.
b = [ZZ(1), ZZ(5)]
for n in range(1, N):
    num = P(n)*b[n] - ZZ(n)^3*b[n-1]
    den = ZZ(n+1)^3
    assert num % den == 0
    b.append(num // den)

# Reciprocal-period source g_n.
S = PowerSeriesRing(QQ, 't', default_prec=N+3)
t = S.gen()
F = S(sum(QQ(b[n])*t^n for n in range(N+1))).add_bigoh(N+2)
Delta = (1 - 34*t + t^2).add_bigoh(N+2)
G = (1/(F^2 * Delta.sqrt())).add_bigoh(N+2)
g = [ZZ(G[n]) for n in range(N+1)]
assert g[:6] == [1, 7, 192, 5520, 165168, 5068320]

# Q7621 Eichler kappa: kappa_0=0, kappa_1=-36, A kappa=-5g for n>=2.
kap = [QQ(0), QQ(-36)]
for n in range(2, N+1):
    kap.append((P(n-1)*kap[n-1] - ZZ(n-1)^3*kap[n-2] - 5*g[n]) / ZZ(n)^3)

# Integral adjacent Casoratian Xi_n.
Xi = [ZZ(-1)]
for n in range(1, N+1):
    Xi.append(ZZ(Xi[-1] - 5*g[n]*b[n-1]))
    assert QQ(Xi[n]) == ZZ(n)^3*(b[n-1]*kap[n] - b[n]*kap[n-1])

# Ratio and exact first difference; C_CM cancels from every difference.
r = [kap[n]/QQ(b[n]) for n in range(N+1)]
s = [r[n] - r[n+1] for n in range(N)]
for n in range(1, N):
    assert s[n-1] == -QQ(Xi[n])/(ZZ(n)^3*b[n-1]*b[n])


def ell3(m):
    if m <= 0:
        return ZZ(1)
    L = ZZ(1)
    for k in range(1, m+1):
        L = lcm(L, ZZ(k))
    return L^3


def shared_den(indices):
    q = ZZ(1)
    for j in indices:
        q = lcm(q, ZZ(kap[j].denominator()))
    assert ell3(max(indices)) % q == 0
    return q


def pair_det(i, j, q):
    # q * (b_i kappa_j - b_j kappa_i), integral for q clearing both endpoints.
    return ZZ(b[i] * (q*kap[j]) - b[j] * (q*kap[i]))


def hankel(seq, n, M, ring):
    return matrix(ring, M, M, lambda i,j: seq[n+i+j]).det()


def logabs_int(z, RR):
    z = abs(ZZ(z))
    if z == 0:
        return -Infinity
    return log(RR(z))


print('Q7674 exact audit N=', N, 'NTEST=', NTEST, 'MMAX=', MMAX)
print('b0..5', b[:6])
print('kappa0..5', kap[:6])
print('Xi0..5', Xi[:6])

# ------------------------------------------------------------------
# A. Exact Stieltjes-necessary tests for the ratio defect.
# If delta_n=r_n-C_CM were a positive Stieltjes moment sequence and delta_n->0,
# it would be Hausdorff; all (-1)^k Delta^k r_n (k>=1) must be nonnegative.
# Also s_n=r_n-r_(n+1) would itself be a Hausdorff/Stieltjes moment sequence.
# ------------------------------------------------------------------
cm_fail = []
cm_checked = 0
for k in range(1, MMAX+1):
    for n in range(1, NTEST+1):
        if n+k > N:
            continue
        d = sum((-1)^(k-j)*binomial(k,j)*r[n+j] for j in range(k+1))
        val = (-1)^k*d
        cm_checked += 1
        if val <= 0:
            cm_fail.append((n,k,val))
            break
    if cm_fail:
        break
print('RATIO_COMPLETE_MONOTONE_CHECKED', cm_checked, 'FIRST_FAIL', cm_fail[:1])

s_hankel_fail = []
s_hankel_checked = 0
for M in range(1, MMAX+1):
    for n in range(1, NTEST+1):
        if n + 2*M - 2 >= len(s):
            continue
        H = hankel(s, n, M, QQ)
        s_hankel_checked += 1
        if H <= 0:
            s_hankel_fail.append((n,M,H))
            break
    if s_hankel_fail:
        break
print('S_HANKEL_POS_CHECKED', s_hankel_checked, 'FIRST_FAIL', s_hankel_fail[:1])

# ------------------------------------------------------------------
# B. Apéry Hankel positivity and kappa Hankel nonvanishing/sign scan.
# ------------------------------------------------------------------
b_fail = []
k_zero = []
for M in range(1, MMAX+1):
    for n in range(0, NTEST+1):
        if n + 2*M - 2 > N:
            continue
        Hb = hankel(b, n, M, ZZ)
        if Hb <= 0:
            b_fail.append((n,M,Hb)); break
        Hk = hankel(kap, n, M, QQ)
        if Hk == 0:
            k_zero.append((n,M)); break
    if b_fail or k_zero:
        break
print('B_HANKEL_FIRST_NONPOS', b_fail[:1])
print('KAPPA_HANKEL_FIRST_ZERO', k_zero[:1])

# Exact target-safety counterexample for scalar b-Hankel at the known (17,13) row.
Hb_12_2 = hankel(b, 12, 2, ZZ)
print('B_HANKEL_TARGET_COUNTEREXAMPLE', 'H_b(12,2)=', Hb_12_2, 'mod17=', Hb_12_2 % 17)
assert b[13] % 17 == 0 and Hb_12_2 % 17 != 0

# ------------------------------------------------------------------
# C. Pair-wedge block Hankel.  Rank is at most two over any field.
# For M=2, the determinant is exactly a Pluecker product; M>=3 vanishes.
# ------------------------------------------------------------------
for n in [1, 10, 50, 100, 140]:
    for M in range(1, MMAX+1):
        if n + 2*M - 1 > N:
            continue
        inds = list(range(n, n+2*M))
        q = shared_den(inds)
        left = list(range(n, n+M))
        right = list(range(n+M, n+2*M))
        W = matrix(ZZ, M, M, lambda i,j: pair_det(left[i], right[j], q))
        detW = W.det()
        if M >= 3:
            assert detW == 0
        if M == 2:
            rhs = pair_det(n,n+1,q)*pair_det(n+2,n+3,q)
            assert detW == rhs
    print('PAIR_WEDGE_RANK_AUDIT_N', n, 'OK')

# p-adic endpoint preservation at p=17, r=13 using one shared denominator.
q17 = shared_den(range(12,16))
for j in [12,14,15]:
    assert pair_det(13,j,q17) % 17 == 0
print('PAIR_TARGET_17_SHARED_DEN', q17, 'OK')

# ------------------------------------------------------------------
# D. Integral Hankel pencil from consecutive pairs:
# R_{n,M}(z)=det[ U_(n+i+j) - q*z*b_(n+i+j) ] = q^M det[kappa-z b].
# q is the minimal shared denominator for the whole coefficient window.
# Record polynomial content and primitive height.  A single target coefficient
# does NOT force p into the content; n=12,M=2,p=17 is the exact counterexample.
# ------------------------------------------------------------------
Z = PolynomialRing(ZZ, 'z').gen()
RR = RealField(100)

pencil_rows = []
for n in [1, 5, 10, 12, 25, 50, 75, 100, 125, 150]:
    for M in [2,3,5,10]:
        if M > MMAX or n + 2*M - 2 > N:
            continue
        inds = list(range(n, n+2*M-1))
        q = shared_den(inds)
        U = {j: ZZ(q*kap[j]) for j in inds}
        Mat = matrix(PolynomialRing(ZZ,'z'), M, M,
                     lambda i,j: U[n+i+j] - q*Z*b[n+i+j])
        poly = Mat.det()
        coeffs = [ZZ(c) for c in poly.list()]
        cont = gcd([abs(c) for c in coeffs if c != 0]) if any(coeffs) else ZZ(0)
        prim_h = max(abs(c//cont) for c in coeffs) if cont else ZZ(0)
        Hb = hankel(b,n,M,ZZ)
        # leading coefficient check
        assert poly[M] == (-q)^M * Hb
        pencil_rows.append((n,M,q,cont,prim_h))
        if n == 12 and M == 2:
            print('PENCIL_12_2_CONTENT', cont, 'mod17=', cont % 17,
                  'degree=',poly.degree())
            assert cont % 17 != 0
        print('PENCIL',n,M,'qbits',q.nbits(),'content_bits',cont.nbits() if cont else 0,
              'prim_height_bits',prim_h.nbits() if prim_h else 0,
              'log_prim/(Mn)', (log(RR(prim_h))/(M*n) if prim_h>0 and n>0 else 0),
              'log_content/(Mn)', (log(RR(cont))/(M*n) if cont>0 and n>0 else 0))

# ------------------------------------------------------------------
# E. Height profiles: b-Hankel, minimally reduced kappa-Hankel,
# exact rational s-Hankel, and the target-safe Xi 2x2 Hankel.
# ------------------------------------------------------------------
for n in [20, 50, 100, 150]:
    for M in [1,2,3,5,10]:
        if n + 2*M - 2 > N:
            continue
        Hb = hankel(b,n,M,ZZ)
        Hk = hankel(kap,n,M,QQ)
        Hs = hankel(s,n,M,QQ) if n+2*M-2 < len(s) else None
        print('HEIGHT',n,M,
              'logHb/(Mn)',logabs_int(Hb,RR)/(M*n),
              'Hk_num_bits',abs(ZZ(Hk.numerator())).nbits(),
              'Hk_den_bits',ZZ(Hk.denominator()).nbits(),
              'logHknum/(Mn)',logabs_int(Hk.numerator(),RR)/(M*n),
              'Hs_sign', (sign(Hs) if Hs is not None else 0),
              'Hs_num_bits', (abs(ZZ(Hs.numerator())).nbits() if Hs is not None else 0),
              'Hs_den_bits', (ZZ(Hs.denominator()).nbits() if Hs is not None else 0))

for n in [10,20,50,100,150]:
    if n+2 > N: continue
    HX2 = Xi[n]*Xi[n+2] - Xi[n+1]^2
    print('XI_H2',n,'sign',sign(HX2),'bits',abs(HX2).nbits(),
          'logabs/n',logabs_int(HX2,RR)/n)

# Target preservation for Xi-H2 at r=13: target at 13 makes Xi_13,Xi_14 zero mod17,
# so both H2 starting at 12 and at 13 are divisible by 17.
assert Xi[13] % 17 == 0 and Xi[14] % 17 == 0
for start in [12,13]:
    HX2 = Xi[start]*Xi[start+2]-Xi[start+1]^2
    assert HX2 % 17 == 0
print('XI_H2_TARGET_17_OK')

print('DONE')
