from sage.all import *
import hashlib
import json
import time

T = ZZ(3670)
J = ZZ(46)
A_T = T//7 + 1
B_T = T//3
H_T = T//J


def Q(n, d=1):
    return QQ(n) / QQ(d)


def qfmt(x):
    x = QQ(x)
    return "%s/%s" % (x.numerator(), x.denominator())


def sha(lines):
    h = hashlib.sha256()
    for s in lines:
        h.update((str(s) + "\n").encode("utf-8"))
    return h.hexdigest()


def apery_P(n):
    n = ZZ(n)
    return 34*n**3 + 51*n**2 + 27*n + 5


def apery_zeros(p):
    F = GF(p)
    b = [F(1), F(5)]
    for n in range(1, p-1):
        b.append((F(apery_P(n))*b[n] - F(n)**3*b[n-1]) / F(n+1)**3)
    assert len(b) == p and b[0] == 1 and b[p-1] == 1
    z = tuple(ZZ(i) for i in range(p) if b[i] == 0)
    zs = set(z)
    for i in range(p):
        assert ((i in zs) == ((p-1-i) in zs))
    return z, zs


def norm_interval(lo, hi):
    lo, hi = ZZ(lo), ZZ(hi)
    return None if hi < lo else (lo, hi)


def ilen(I):
    return ZZ(0) if I is None else ZZ(I[1]-I[0]+1)


def ivals(I):
    return tuple() if I is None else tuple(range(I[0], I[1]+1))


def ifmt(I):
    return "EMPTY" if I is None else "%s:%s" % I


def icap(I, K):
    return None if I is None or K is None else norm_interval(max(I[0],K[0]), min(I[1],K[1]))


def row_interval(q):
    return norm_interval(max(ZZ(0), T-6*q+1), min(q-1, 2*T-6*q))


def literal_masks(p, q, t):
    if not (A_T <= p <= B_T and 0 <= t < q and T < 6*q+t <= 2*T):
        return False, False
    d = q-p
    if J*d <= T:
        return False, False
    ell = max(A_T, (6*q+t)//7 + 1)
    cap = min(B_T, q-H_T-1)
    return (ell <= p <= min(cap,q-t-1)), (ell <= p <= min(cap,(12*q+t)//13))


def rank_intervals(p, q):
    I0 = row_interval(q)
    if I0 is None or not (A_T <= p <= B_T) or J*(q-p) <= T:
        return None, None
    lo, hi = I0
    d = q-p
    Im = norm_interval(lo, min(hi, d-1, p-6*d-1))
    Ip = norm_interval(max(lo,p-12*d), min(hi,p-6*d-1))
    return Im, Ip


def pos(x):
    return x if x > 0 else Q(0)


t0 = time.time()
primes = tuple(ZZ(p) for p in prime_range(A_T, B_T+1))
prime_set = set(primes)
Z = {}
ZS = {}
N = {}
for p in primes:
    Z[p], ZS[p] = apery_zeros(p)
    N[p] = ZZ(len(Z[p]))

T3 = {}
for p in primes:
    inv7 = inverse_mod(ZZ(7),p)
    T3[p] = ZZ(sum(1 for a in Z[p] for r in Z[p]
                   if ZZ((inv7*(6*a+13*r+6)) % p) in ZS[p]))
    defect = N[p]**2-T3[p]
    assert 0 <= T3[p] <= N[p]**2
    assert N[p]**2//4 <= defect <= N[p]*(N[p]-1)

# ---------- literal q,t,p path ----------
rowsD = {}
bitsD = {}
pqD = {}
geomD = []
count = dict(prime_count=ZZ(len(primes)),q_rows=ZZ(0),q_geometry=ZZ(0),
             qtp_tests=ZZ(0),minus_masks=ZZ(0),plus_masks=ZZ(0),cap_masks=ZZ(0),
             minus_hits=ZZ(0),plus_hits=ZZ(0),union_hits=ZZ(0))

for q in primes:
    I0 = row_interval(q)
    if I0 is None:
        continue
    count["q_rows"] += 1
    lo, hi = I0
    Nq, Zq = N[q], ZS[q]
    pst = {}
    sumA = ZZ(0)
    hA = ZZ(0)
    dq_literal = Q(0)
    any_geom = False
    for t in range(lo,hi+1):
        zq = ZZ(t in Zq)
        gnum = q*zq-Nq
        At = ZZ(0)
        Rt = Q(0)
        for p in primes:
            count["qtp_tests"] += 1
            cm, cp = literal_masks(p,q,ZZ(t))
            if not (cm or cp):
                continue
            any_geom = True
            d = q-p
            r = 6*d+t
            assert 0 <= r < p
            st = pst.setdefault(p,dict(M=[],P=[],C=[],Y=set(),gY=ZZ(0),gR=ZZ(0)))
            mb = 0
            if cm:
                count["minus_masks"] += 1
                st["M"].append(ZZ(t)); mb |= 1
            if cp:
                count["plus_masks"] += 1
                st["P"].append(ZZ(t)); mb |= 2
            if cm and cp:
                count["cap_masks"] += 1
                st["C"].append(ZZ(t))
            geomD.append("%s|%s|%s|%s" % (q,t,p,mb))
            xm = False
            if cm:
                a = d-t-1
                assert 0 <= a < p
                xm = (a in ZS[p] and r in ZS[p])
            xp = False
            if cp:
                w = 12*d+t-p
                assert 0 <= w < p and w == ((12*d+t) % p)
                xp = (r in ZS[p] and w in ZS[p])
            eb = (1 if xm else 0) | (2 if xp else 0)
            count["minus_hits"] += ZZ(xm)
            count["plus_hits"] += ZZ(xp)
            y = ZZ(bool(eb))
            if y:
                count["union_hits"] += 1
                key = (q,ZZ(t),p)
                assert key not in bitsD
                bitsD[key] = eb
                st["Y"].add(ZZ(t))
            cnum = N[p]**2*(ZZ(cm)+ZZ(cp))-T3[p]*ZZ(cm and cp)
            At += y
            Rt += Q(cnum,p**2)
            st["gY"] += gnum*y
            st["gR"] += gnum*cnum
        sumA += At
        hA += zq*At
        dq_literal += Q(gnum,q)*(Q(At)-Rt)
    if any_geom:
        count["q_geometry"] += 1
    event_cov = Q(0); SR = Q(0); dq = Q(0)
    for p,st in pst.items():
        event_cov += Q(st["gY"],q)
        SR += Q(st["gR"],q*p**2)
        dq += Q(p**2*st["gY"]-st["gR"],q*p**2)
        pqD[(q,p)] = st
    assert dq == dq_literal == event_cov-SR
    g2 = Q(sum((q*ZZ(t in Zq)-Nq)**2 for t in range(q)),q**2)
    G = Q(Nq*(q-Nq),q)
    assert g2 == G
    if G == 0:
        assert dq == 0
        Pq = Q(0)
    else:
        Pq = pos(dq)**2/G
    rowsD[q] = dict(I0=I0,Nq=Nq,G=G,sumA=sumA,hA=hA,event=event_cov,SR=SR,dq=dq,Pq=Pq,k=ZZ(len(pst)))

# ---------- inverse zero-pair event path ----------
bitsS = {}
genM = ZZ(0); genP = ZZ(0); dup = ZZ(0)

def add_event(q,t,p,b):
    global dup
    k = (ZZ(q),ZZ(t),ZZ(p))
    old = bitsS.get(k,0)
    if old & b:
        dup += 1
    bitsS[k] = old | b

for p in primes:
    for a in Z[p]:
        for r in Z[p]:
            s = a+r+1
            if s % 7 == 0:
                d = ZZ(s//7); t = ZZ(d-a-1); q = ZZ(p+d)
                if q in prime_set and literal_masks(p,q,t)[0]:
                    assert r == 6*d+t
                    add_event(q,t,p,1); genM += 1
    for z in Z[p]:
        for w in Z[p]:
            s = z+w+1
            if s % 6 == 0:
                d = ZZ(s//6); t = ZZ(p-1-z-6*d); q = ZZ(p+d)
                if q in prime_set and literal_masks(p,q,t)[1]:
                    assert z == p-1-(6*d+t) and w == 12*d+t-p
                    add_event(q,t,p,2); genP += 1

ediff = sorted(set(bitsD)^set(bitsS))
ebit = sorted(k for k in set(bitsD)&set(bitsS) if bitsD[k] != bitsS[k])
first_event_q = min([k[0] for k in ediff+ebit]) if ediff or ebit else None
byPQ = {}
for (q,t,p),b in bitsS.items():
    byPQ.setdefault((q,p),set()).add(t)

# ---------- interval endpoint + upper-zero prefix-rank path ----------
prefix = {}
for q in rowsD:
    a = [ZZ(0)]*(q+1); s = ZZ(0)
    for t in range(q):
        s += ZZ(t in ZS[q]); a[t+1] = s
    assert s == N[q]
    prefix[q] = a

def irank(q,I):
    return ZZ(0) if I is None else prefix[q][I[1]+1]-prefix[q][I[0]]

rowsS = {}
geomR = []
metaR = []
imismatch = []
staticPQ = ZZ(0); geomPQ = ZZ(0)
for q in sorted(rowsD):
    Nq = N[q]
    sumA = ZZ(0); hA = ZZ(0); event = Q(0); SR = Q(0); dq = Q(0); kq = ZZ(0)
    for p in primes:
        if J*(q-p) <= T:
            continue
        staticPQ += 1
        Im,Ip = rank_intervals(p,q); Ic = icap(Im,Ip)
        dm = tuple(pqD.get((q,p),{}).get("M",[]))
        dp = tuple(pqD.get((q,p),{}).get("P",[]))
        dc = tuple(pqD.get((q,p),{}).get("C",[]))
        rm,rp,rc = ivals(Im),ivals(Ip),ivals(Ic)
        if dm != rm or dp != rp or dc != rc:
            imismatch.append((q,p,dm,rm,dp,rp,dc,rc))
        for t in sorted(set(rm)|set(rp)):
            geomR.append("%s|%s|%s|%s" % (q,t,p,(1 if t in set(rm) else 0)|(2 if t in set(rp) else 0)))
        if Im is None and Ip is None:
            assert not byPQ.get((q,p),set())
            continue
        geomPQ += 1; kq += 1
        nm,np,nc = ilen(Im),ilen(Ip),ilen(Ic)
        mm,mp,mc = irank(q,Im),irank(q,Ip),irank(q,Ic)
        Sm,Sp,Sc = q*mm-Nq*nm, q*mp-Nq*np, q*mc-Nq*nc
        ts = byPQ.get((q,p),set())
        m = ZZ(len(ts)); h = ZZ(sum(t in ZS[q] for t in ts))
        sumA += m; hA += h
        A = q*h-Nq*m
        Rnum = N[p]**2*(Sm+Sp)-T3[p]*Sc
        event += Q(A,q); SR += Q(Rnum,q*p**2); dq += Q(p**2*A-Rnum,q*p**2)
        metaR.append("%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s" %
                     (q,p,ifmt(Im),ifmt(Ip),ifmt(Ic),nm,np,nc,mm,mp,mc,Sm,Sp,Sc,N[p],T3[p],m))
    assert dq == event-SR
    G = Q(Nq*(q-Nq),q)
    if G == 0:
        assert dq == 0
        Pq = Q(0)
    else:
        Pq = pos(dq)**2/G
    rowsS[q] = dict(I0=rowsD[q]["I0"],Nq=Nq,G=G,sumA=sumA,hA=hA,event=event,SR=SR,dq=dq,Pq=Pq,k=kq)

fields = ("Nq","G","sumA","hA","event","SR","dq","Pq","k")
rdiff = [(q,f,rowsD[q][f],rowsS[q][f]) for q in sorted(rowsD) for f in fields if rowsD[q][f] != rowsS[q][f]]
first_row_q = rdiff[0][0] if rdiff else None
first_interval_q = imismatch[0][0] if imismatch else None
all_mq = [q for q in (first_event_q,first_interval_q,first_row_q) if q is not None]
smallest_mismatch = min(all_mq) if all_mq else None

def totals(rows):
    G = sum((r["G"] for r in rows.values()),Q(0))
    Pplus = sum((r["Pq"] for r in rows.values()),Q(0))
    Dplus = sum((pos(r["dq"]) for r in rows.values()),Q(0))
    Daff = sum((r["dq"] for r in rows.values()),Q(0))
    Mmark = ZZ(sum(r["hA"] for r in rows.values()))
    Mq0 = sum((Q(r["Nq"]*r["sumA"],q) for q,r in rows.items()),Q(0))
    SR = sum((r["SR"] for r in rows.values()),Q(0))
    residual = Q(Mmark)-Mq0-SR-Daff
    slack = G*Pplus-Dplus**2
    assert residual == 0 and slack >= 0
    return dict(G=G,Pplus=Pplus,Dplus=Dplus,Daff=Daff,Mmark=Mmark,Mq0=Mq0,SR=SR,residual=residual,slack=slack)

TD,TS = totals(rowsD),totals(rowsS)
assert TD == TS

def event_lines(bits):
    return ["%s|%s|%s|%s" % (q,t,p,bits[(q,t,p)]) for q,t,p in sorted(bits)]

def row_lines(rows):
    return ["%s|%s|%s|%s|%s|%s|%s|%s|%s|%s" %
            (q,r["Nq"],qfmt(r["G"]),r["sumA"],r["hA"],qfmt(r["event"]),qfmt(r["SR"]),qfmt(r["dq"]),qfmt(r["Pq"]),r["k"])
            for q,r in sorted(rows.items())]

hashes = dict(
    zero=sha(["%s|%s" % (p,",".join(map(str,Z[p]))) for p in primes]),
    T3=sha(["%s|%s|%s" % (p,N[p],T3[p]) for p in primes]),
    direct_geometry=sha(sorted(geomD)),
    rank_geometry=sha(sorted(geomR)),
    direct_events=sha(event_lines(bitsD)),
    stream_events=sha(event_lines(bitsS)),
    direct_rows=sha(row_lines(rowsD)),
    stream_rows=sha(row_lines(rowsS)),
    interval_meta=sha(sorted(metaR)),
)
assert hashes["direct_geometry"] == hashes["rank_geometry"]
assert hashes["direct_events"] == hashes["stream_events"]
assert hashes["direct_rows"] == hashes["stream_rows"]

count.update(dict(direct_event_keys=ZZ(len(bitsD)),stream_event_keys=ZZ(len(bitsS)),
                  generated_minus=genM,generated_plus=genP,duplicate_same_sign=dup,
                  event_key_differences=ZZ(len(ediff)),event_bit_differences=ZZ(len(ebit)),
                  interval_mismatches=ZZ(len(imismatch)),row_mismatches=ZZ(len(rdiff)),
                  static_pq=staticPQ,geometry_pq=geomPQ,
                  q_G_positive=ZZ(sum(r["G"]>0 for r in rowsD.values())),
                  q_dq_nonzero=ZZ(sum(r["dq"]!=0 for r in rowsD.values())),
                  q_dq_positive=ZZ(sum(r["dq"]>0 for r in rowsD.values())),
                  q_union_nonzero=ZZ(sum(r["sumA"]>0 for r in rowsD.values()))))

geom_q = [q for q,r in rowsD.items() if r["k"]>0]
ggeom_q = [q for q in geom_q if rowsD[q]["G"]>0]
union_q = [q for q,r in rowsD.items() if r["sumA"]>0]
ndq_q = [q for q,r in rowsD.items() if r["dq"]!=0]
pdq_q = [q for q,r in rowsD.items() if r["dq"]>0]
focus = []
for x in [min(geom_q) if geom_q else None,min(ggeom_q) if ggeom_q else None,
          min(union_q) if union_q else None,min(ndq_q) if ndq_q else None,
          min(pdq_q) if pdq_q else None,ZZ(1129) if ZZ(1129) in rowsD else None]:
    if x is not None and x not in focus:
        focus.append(x)

out = []
out.append("Q7503_PAUC_REPRODUCTION_V2")
out.append("PARAM T=%s J=%s A_T=%s B_T=%s H_T=%s" % (T,J,A_T,B_T,H_T))
out.append("SMALLEST_MISMATCH_Q %s" % ("NONE" if smallest_mismatch is None else smallest_mismatch))
for k in sorted(count): out.append("COUNTER %s=%s" % (k,count[k]))
for k in sorted(hashes): out.append("HASH %s=%s" % (k,hashes[k]))
for k in ("G","Pplus","Dplus","Daff"): out.append("GLOBAL %s=%s" % (k,qfmt(TD[k])))
out.append("IDENTITY Mmark=%s Mq0=%s SR=%s Daff=%s residual=%s" %
           (TD["Mmark"],qfmt(TD["Mq0"]),qfmt(TD["SR"]),qfmt(TD["Daff"]),qfmt(TD["residual"])))
out.append("CAUCHY slack=G*Pplus-Dplus^2=%s" % qfmt(TD["slack"]))
out.append("FOCUS_QS %s" % (",".join(map(str,focus)) if focus else "NONE"))
for q in focus:
    r = rowsD[q]
    out.append("ROW q=%s I0=%s Zq=%s Nq=%s Gq=%s sumA=%s hA=%s event=%s SR=%s dq=%s Pplus_q=%s geometry_p=%s" %
               (q,ifmt(r["I0"]),",".join(map(str,Z[q])),r["Nq"],qfmt(r["G"]),r["sumA"],r["hA"],qfmt(r["event"]),qfmt(r["SR"]),qfmt(r["dq"]),qfmt(r["Pq"]),r["k"]))
    gl = []
    for p in primes:
        Im,Ip = rank_intervals(p,q)
        if Im is not None or Ip is not None:
            gl.append("p=%s,Np=%s,T3=%s,Im=%s,Ip=%s,Ic=%s" % (p,N[p],T3[p],ifmt(Im),ifmt(Ip),ifmt(icap(Im,Ip))))
    out.append("ROW_GEOMETRY q=%s %s" % (q,";".join(gl)))
out.append("SECONDS %.6f" % (time.time()-t0))

with open("q7503_pauc_reproduction.txt","w") as f:
    f.write("\n".join(out)+"\n")
print("\n".join(out))
