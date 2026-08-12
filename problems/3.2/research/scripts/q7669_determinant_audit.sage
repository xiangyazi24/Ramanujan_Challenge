from sage.all import *

N = 500

def P(n):
    n = ZZ(n)
    return 34*n^3 + 51*n^2 + 27*n + 5

# Apéry b_n exactly.
b = [ZZ(1), ZZ(5)]
for n in range(1, N):
    num = P(n)*b[n] - ZZ(n)^3*b[n-1]
    den = ZZ(n+1)^3
    assert num % den == 0
    b.append(num // den)

# Reciprocal-period source g_n = [t^n] 1/(F^2 sqrt(1-34t+t^2)).
S = PowerSeriesRing(QQ, 't', default_prec=N+3)
t = S.gen()
F = S(sum(QQ(b[n])*t^n for n in range(N+1))).add_bigoh(N+2)
Delta = (1 - 34*t + t^2).add_bigoh(N+2)
G = (1/(F^2 * Delta.sqrt())).add_bigoh(N+2)
g = [ZZ(G[n]) for n in range(N+1)]
assert g[:6] == [1, 7, 192, 5520, 165168, 5037696]

# Q7621 Eichler kappa: kappa_0=0, kappa_1=-36 and A kappa = -5 g for r>=2.
kap = [QQ(0), QQ(-36)]
for r in range(2, N+1):
    kap.append((P(r-1)*kap[r-1] - ZZ(r-1)^3*kap[r-2] - 5*g[r]) / ZZ(r)^3)

# Integral Green/Casoratian coordinate.
Xi = [ZZ(-1)]
for r in range(1, N+1):
    Xi.append(ZZ(Xi[-1] - 5*g[r]*b[r-1]))
    assert QQ(Xi[r]) == ZZ(r)^3*(b[r-1]*kap[r] - b[r]*kap[r-1])

# Minimal denominator-cleared determinant numerator.
def D(m,n):
    return QQ(b[m])*kap[n] - QQ(b[n])*kap[m]

def prim_num(m,n):
    return ZZ(D(m,n).numerator())

def prim_den(m,n):
    return ZZ(D(m,n).denominator())

def fresh_part(x, cutoff):
    x = abs(ZZ(x))
    if x == 0:
        return ZZ(0)
    out = ZZ(1)
    for p,e in factor(x):
        if p > cutoff:
            out *= p^e
    return out

def fresh_factor(x, cutoff):
    x = abs(ZZ(x))
    if x in (0,1):
        return str(x)
    return ' * '.join([('%s^%s' % (p,e) if e != 1 else str(p)) for p,e in factor(x) if p > cutoff]) or '1'

# Exact recurrence / adjacent Casoratian checks.
for m in [0,1,2,7,13,100]:
    for n in range(max(2,m+2), min(N,m+20)+1):
        lhs = ZZ(n)^3*D(m,n) - P(n-1)*D(m,n-1) + ZZ(n-1)^3*D(m,n-2)
        assert lhs == -5*b[m]*g[n]
for n in range(1,N+1):
    assert D(n-1,n) == QQ(Xi[n], ZZ(n)^3)

print('Q7669 exact audit N=',N)
print('b[0:6]=',b[:6])
print('kappa[0:6]=',kap[:6])
print('Xi[0:6]=',Xi[:6])

# Factor the true row primitive contents c_r=gcd(b_r,num(kappa_r)).
nontrivial=[]
fresh=[]
for r in range(1,N+1):
    u=abs(ZZ(kap[r].numerator()))
    c=gcd(abs(b[r]),u)
    if c>1:
        nontrivial.append((r,c,factor(c)))
    ff=[(ZZ(p),ZZ(e)) for p,e in factor(c) if p>r]
    if ff:
        fresh.append((r,ff))
print('NONTRIVIAL_ROW_CONTENT_COUNT',len(nontrivial))
for row in nontrivial:
    print('ROW_CONTENT',row[0],row[1],row[2])
print('FRESH_ROW_CONTENTS',fresh)

# Fixed anchor (0,1): exact primitive-minor gcd.  Since D_01=-36, all p>3
# support must equal the row content support.
anchor_mismatches=[]
anchor_nontrivial=[]
for r in range(2,N+1):
    n0=abs(prim_num(0,r)); n1=abs(prim_num(1,r))
    h=gcd(n0,n1)
    c=gcd(abs(b[r]),abs(ZZ(kap[r].numerator())))
    # compare factors above max(r,3), which is the target range relevant here
    hp=fresh_part(h,r); cp=fresh_part(c,r)
    if hp != cp:
        anchor_mismatches.append((r,hp,cp))
    if h>1:
        anchor_nontrivial.append((r,h,factor(h)))
print('ANCHOR01_MISMATCHES',anchor_mismatches)
print('ANCHOR01_NONTRIVIAL_COUNT',len(anchor_nontrivial))
for row in anchor_nontrivial:
    print('ANCHOR01_CONTENT',row[0],row[1],row[2])

# Dyadic moving-anchor gcds G_{a,b}(r)=gcd(num D_{a,r}, num D_{b,r}).
# Summarize exact fresh support and factor the known target rows.
for R in [8,16,32,64,128,256]:
    lo=R+1; hi=min(2*R,N)
    if lo>hi: continue
    a=lo; aa=lo+1
    if aa>hi: continue
    bad=[]; exact=0; extras=0; max_extra_bits=0
    target_rows=[]
    for r in range(lo,hi+1):
        if r in (a,aa): continue
        h=gcd(abs(prim_num(a,r)),abs(prim_num(aa,r)))
        c=gcd(abs(b[r]),abs(ZZ(kap[r].numerator())))
        hf=fresh_part(h,2*R); cf=fresh_part(c,2*R)
        if hf==cf: exact += 1
        else: bad.append((r,hf,cf))
        if cf>1:
            target_rows.append((r,h,cf,factor(h)))
        extra = h // gcd(h,c)
        if extra>1:
            extras += 1
            max_extra_bits=max(max_extra_bits,extra.nbits())
    print('DYADIC',R,'anchors',a,aa,'rows',hi-lo+1,'fresh_exact',exact,'bad_count',len(bad),'extra_rows',extras,'max_extra_bits',max_extra_bits)
    if bad:
        print('DYADIC_BAD_SAMPLE',R,bad[:10])
    for z in target_rows:
        print('DYADIC_TARGET_CONTENT',R,z[0],'raw_gcd',z[1],'true_fresh',z[2],'raw_factor',z[3])

# Sparse edge-product lower-bound data: adjacent matching in each dyadic block.
# Report actual sum log|primitive numerator| / R^2 and analytic raw determinant version.
RR = RealField(100)
for R in [16,32,64,128,250]:
    lo=R+1; hi=min(2*R,N)
    verts=list(range(lo,hi+1))
    if len(verts)%2==1: verts=verts[:-1]
    edges=[(verts[i],verts[i+1]) for i in range(0,len(verts),2)]
    slog=RR(0); rawlog=RR(0); denlog=RR(0)
    for m,n in edges:
        q=prim_den(m,n); z=abs(prim_num(m,n))
        assert z>0
        slog += log(RR(z))
        rawlog += log(abs(RR(D(m,n))))
        denlog += log(RR(q))
    print('MATCH_HEIGHT',R,'edges',len(edges),'log_num/R2',slog/RR(R^2),'log_raw/R2',rawlog/RR(R^2),'log_den/R2',denlog/RR(R^2))

# Selected exact determinant primitive factorizations where feasible.
for m,n in [(0,1),(0,13),(1,13),(12,13),(13,14),(491,492),(492,493)]:
    z=prim_num(m,n); q=prim_den(m,n)
    print('DET',m,n,'den',q,'num_bits',abs(z).nbits())
    print('DET_FACTOR',m,n,factor(abs(z)))

print('DONE')
