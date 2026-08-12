from sage.all import *

N = 600
RMAX = 300

def P(n):
    n = ZZ(n)
    return 34*n^3 + 51*n^2 + 27*n + 5

# Exact Apéry sequence.
b = [ZZ(1), ZZ(5)]
for n in range(1, N):
    num = P(n)*b[n] - ZZ(n)^3*b[n-1]
    den = ZZ(n+1)^3
    assert num % den == 0
    b.append(num // den)

# Exact reciprocal-period source and Q7621 Eichler sequence.
S = PowerSeriesRing(QQ, 't', default_prec=N+3)
t = S.gen()
F = S(sum(QQ(b[n])*t^n for n in range(N+1))).add_bigoh(N+2)
Delta = (1 - 34*t + t^2).add_bigoh(N+2)
G = (1/(F^2 * Delta.sqrt())).add_bigoh(N+2)
g = [ZZ(G[n]) for n in range(N+1)]
assert g[:6] == [1, 7, 192, 5520, 165168, 5068320]

kap = [QQ(0), QQ(-36)]
for r in range(2, N+1):
    kap.append((P(r-1)*kap[r-1] - ZZ(r-1)^3*kap[r-2] - 5*g[r]) / ZZ(r)^3)

Xi = [ZZ(-1)]
for r in range(1, N+1):
    Xi.append(ZZ(Xi[-1] - 5*g[r]*b[r-1]))
    assert QQ(Xi[r]) == ZZ(r)^3*(b[r-1]*kap[r] - b[r]*kap[r-1])

print('Q7672 BFH SAGE EXACT N', N, 'RMAX', RMAX)
print('b0_5', b[:6])
print('kappa0_5', kap[:6])
print('Xi0_5', Xi[:6])

# Factor each row content once.  For p>r, gcd(b_r,Xi_r) is exactly the
# fresh transverse p-support, with valuation retained.
row_fac = {}
for r in range(1, N+1):
    c = gcd(abs(b[r]), abs(Xi[r]))
    if c > 1:
        row_fac[r] = [(ZZ(p), ZZ(e)) for p,e in factor(c)]

fresh_rows = []
for r, fac in row_fac.items():
    ff = [(p,e) for p,e in fac if p > r]
    if ff:
        fresh_rows.append((r, ff))
print('FRESH_ROWS_P_GT_R', fresh_rows)

# Exact fresh content of D_R(Y)=prod_{R<r<=2R}(b_r+Y Xi_r), restricted p>2R.
# Report all nonzero R and compact interval summary.
RR = RealField(100)
nonzero = []
max_ratio = (RR(0), None, ZZ(1))
for R in range(1, RMAX+1):
    fresh = ZZ(1)
    for r in range(R+1, 2*R+1):
        for p,e in row_fac.get(r, []):
            if p > 2*R:
                fresh *= p^e
    lv = RR(log(RR(fresh))) if fresh > 1 else RR(0)
    ratio = lv / RR(R^2)
    if fresh > 1:
        nonzero.append((R, fresh, ratio))
    if ratio > max_ratio[0]:
        max_ratio = (ratio, R, fresh)
print('BFH_NONZERO_COUNT', len(nonzero))
for R, fresh, ratio in nonzero:
    print('BFH_NONZERO', R, fresh, ratio)
print('BFH_MAX_RATIO', max_ratio)

# Compress consecutive R having the same fresh content.
segments=[]
if nonzero:
    a=nonzero[0][0]; z=nonzero[0][0]; val=nonzero[0][1]
    for R,fresh,ratio in nonzero[1:]:
        if R == z+1 and fresh == val:
            z=R
        else:
            segments.append((a,z,val))
            a=z=R; val=fresh
    segments.append((a,z,val))
print('BFH_NONZERO_SEGMENTS', segments)

# Disprove small fixed-order homogeneous P-recursive relations for Xi:
# sum_{j=0}^ord Q_j(n) Xi_{n+j}=0, deg Q_j <= deg.
# Full column rank mod one prime rules out a rational relation in that box.
def test_polyrec(order, degree, prime=1000003, start=5):
    K = GF(prime)
    unk = (order+1)*(degree+1)
    rows = min(N-order-start+1, unk+30)
    A = matrix(K, rows, unk)
    for ii in range(rows):
        n = start + ii
        col = 0
        for j in range(order+1):
            x = K(Xi[n+j])
            ne = K(1)
            nn = K(n)
            for e in range(degree+1):
                A[ii,col] = ne*x
                ne *= nn
                col += 1
    return (A.rank(), unk, rows)

for od in [(4,4),(6,6),(8,8),(10,10)]:
    rk,unk,rows=test_polyrec(*od)
    print('POLYREC_TEST', od[0], od[1], 'rank', rk, 'unknowns', unk, 'rows', rows, 'FULL', rk==unk)

# Fixed-level modular-form / modular-symbol sizes.  These are candidate ambient
# congruence modules only; the BFH proof does not assume that Taylor coefficients
# are modular-symbol coordinates.
for k in [2,4]:
    try:
        MF = ModularForms(Gamma0(6), k)
        E = MF.eisenstein_subspace()
        Cusp = MF.cuspidal_subspace()
        print('MODFORMS', k, 'dim', MF.dimension(), 'eis', E.dimension(), 'cusp', Cusp.dimension(), 'sturm', MF.sturm_bound())
    except Exception as ex:
        print('MODFORMS_ERROR', k, repr(ex))
    try:
        MS = ModularSymbols(Gamma0(6), k, sign=0)
        CS = MS.cuspidal_subspace()
        print('MODSYMS', k, 'dim', MS.dimension(), 'cusp', CS.dimension())
        try:
            T = MS.hecke_algebra()
            print('HECKE', k, 'rank', T.rank())
        except Exception as ex:
            print('HECKE_ERROR', k, repr(ex))
    except Exception as ex:
        print('MODSYMS_ERROR', k, repr(ex))

# A naive fixed-q-expansion direct-sum transition determinant, only as a
# computational candidate for an Eisenstein/cuspidal congruence index.
# It is basis-normalization dependent and is NOT used as a theorem.
for k in [2,4]:
    try:
        MF = ModularForms(Gamma0(6), k)
        d = MF.dimension()
        if d == 0:
            continue
        prec = max(12, MF.sturm_bound()+d+3)
        BM = list(MF.q_echelon_basis(prec=prec))
        BE = list(MF.eisenstein_subspace().q_echelon_basis(prec=prec))
        BS = list(MF.cuspidal_subspace().q_echelon_basis(prec=prec))
        BES = BE + BS
        if len(BES) != d:
            print('QINDEX_SKIP', k, 'basis lengths', len(BM), len(BES))
            continue
        AM = matrix(QQ, [[f[n] for n in range(prec)] for f in BM])
        AES = matrix(QQ, [[f[n] for n in range(prec)] for f in BES])
        piv = AM.echelon_form().pivots()
        # pivots() here are row pivots for transposed convention; use independent
        # columns selected from AM.transpose().pivot rows.
        cols = AM.transpose().echelon_form().pivots()
        cols = list(cols[:d])
        if len(cols) != d:
            print('QINDEX_SKIP', k, 'independent columns', cols)
            continue
        MM = AM.matrix_from_columns(cols)
        MES = AES.matrix_from_columns(cols)
        qidx = abs(MES.det()/MM.det())
        print('QEXP_DIRECTSUM_DET_RATIO', k, qidx, 'cols', cols)
    except Exception as ex:
        print('QINDEX_ERROR', k, repr(ex))

print('DONE')
