from sage.all import *
import hashlib
import json
import time

# Independent reproduction point requested in Q7503.
T = ZZ(3670)
J = ZZ(46)
A_T = T // 7 + 1
B_T = T // 3
H_T = T // J


def apery_coeff_poly(n):
    n = ZZ(n)
    return 34*n**3 + 51*n**2 + 27*n + 5


def apery_zero_set(p):
    """Exact Apéry recurrence in GF(p), including both endpoints."""
    F = GF(p)
    b = [F(1), F(5)]
    for n in range(1, p-1):
        nxt = (F(apery_coeff_poly(n))*b[n] - F(n)**3*b[n-1]) / F(n+1)**3
        b.append(nxt)
    assert len(b) == p
    assert b[0] == 1 and b[p-1] == 1
    Z = tuple(ZZ(r) for r in range(p) if b[r] == 0)
    Zset = set(Z)
    # Apéry reflection is a validation, not an input to the recurrence.
    for r in range(p):
        assert ((r in Zset) == ((p-1-r) in Zset))
    return Z, Zset


def sha256_lines(lines):
    h = hashlib.sha256()
    for line in lines:
        h.update((str(line) + "\n").encode("utf-8"))
    return h.hexdigest()


def qfmt(x):
    x = QQ(x)
    return "%s/%s" % (x.numerator(), x.denominator())


def interval_fmt(I):
    if I is None:
        return "EMPTY"
    return "%s:%s" % (I[0], I[1])


def normalized_interval(lo, hi):
    lo = ZZ(lo)
    hi = ZZ(hi)
    if hi < lo:
        return None
    return (lo, hi)


def interval_length(I):
    if I is None:
        return ZZ(0)
    return ZZ(I[1] - I[0] + 1)


def interval_values(I):
    if I is None:
        return tuple()
    return tuple(range(ZZ(I[0]), ZZ(I[1]) + 1))


def intersect_intervals(I, K):
    if I is None or K is None:
        return None
    return normalized_interval(max(I[0], K[0]), min(I[1], K[1]))


def row_interval(q):
    lo = max(ZZ(0), T - 6*q + 1)       # T < 6q+t
    hi = min(q-1, 2*T - 6*q)           # 6q+t <= 2T
    return normalized_interval(lo, hi)


def masks(p, q, t):
    """Literal floor-aware quotient-six p-window masks."""
    if not (A_T <= p <= B_T):
        return False, False
    if not (0 <= t < q and T < 6*q+t <= 2*T):
        return False, False
    d = q-p
    if not (J*d > T):                  # strict d > T/J
        return False, False
    ell = max(A_T, (6*q+t)//7 + 1)      # strict p > (6q+t)/7
    cap = min(B_T, q-H_T-1)             # exact integer long-gap cap
    minus_hi = min(cap, q-t-1)
    plus_hi = min(cap, (12*q+t)//13)
    return (ell <= p <= minus_hi), (ell <= p <= plus_hi)


def mask_intervals(p, q):
    """Independent endpoint form used by the streamed interval-rank path."""
    Irow = row_interval(q)
    if Irow is None or not (A_T <= p <= B_T) or not (J*(q-p) > T):
        return None, None
    lo, hi = Irow
    common_hi = min(hi, 7*p-6*q-1)
    Iminus = normalized_interval(lo, min(common_hi, q-p-1))
    Iplus = normalized_interval(max(lo, 13*p-12*q), common_hi)
    return Iminus, Iplus


def positive_part(x):
    x = QQ(x)
    return x if x > 0 else QQ(0)


t0 = time.time()
primes = tuple(ZZ(p) for p in prime_range(A_T, B_T+1))
prime_set = set(primes)

zeros = {}
zero_sets = {}
N = {}
for p in primes:
    Z, Zset = apery_zero_set(p)
    zeros[p] = Z
    zero_sets[p] = Zset
    N[p] = ZZ(len(Z))

# Exact complete-plane three-zero overlap coefficient.
T3 = {}
for p in primes:
    inv7 = inverse_mod(ZZ(7), p)
    count = ZZ(0)
    Zp = zero_sets[p]
    for a in zeros[p]:
        for r in zeros[p]:
            b = ZZ((inv7*(6*a + 13*r + 6)) % p)
            if b in Zp:
                count += 1
    T3[p] = count
    assert 0 <= T3[p] <= N[p]**2

zero_lines = ["%s|%s" % (p, ",".join(map(str, zeros[p]))) for p in primes]
t3_lines = ["%s|%s|%s" % (p, N[p], T3[p]) for p in primes]

# ---------------------------------------------------------------------------
# Path A: literal q,t,p loop.  No interval endpoint or inverse-pair formula is
# used to decide a mask or an actual event.
# ---------------------------------------------------------------------------
rows_direct = {}
direct_bits = {}
direct_geom_lines = []
direct_pq = {}

counters = {
    "prime_count": ZZ(len(primes)),
    "q_row_count": ZZ(0),
    "q_geometry_count": ZZ(0),
    "qtp_tests": ZZ(0),
    "minus_mask_leaves": ZZ(0),
    "plus_mask_leaves": ZZ(0),
    "intersection_mask_leaves": ZZ(0),
    "minus_actual_hits": ZZ(0),
    "plus_actual_hits": ZZ(0),
    "union_actual_hits": ZZ(0),
}

for q in primes:
    Irow = row_interval(q)
    if Irow is None:
        continue
    counters["q_row_count"] += 1
    lo, hi = Irow
    Nq = N[q]
    Zq = zero_sets[q]

    pqstats = {}
    sumA = ZZ(0)
    hA = ZZ(0)
    dq_literal = QQ(0)
    row_has_geometry = False

    for t in range(lo, hi+1):
        zqt = ZZ(1 if t in Zq else 0)
        gnum = q*zqt - Nq
        At = ZZ(0)
        Rt = QQ(0)

        for p in primes:
            counters["qtp_tests"] += 1
            cm, cp = masks(p, q, ZZ(t))
            if not (cm or cp):
                continue
            row_has_geometry = True
            d = q-p
            r = 6*d+t
            assert 0 <= r < p

            st = pqstats.setdefault(p, {
                "minus_ts": [], "plus_ts": [], "inter_ts": [],
                "union_ts": set(), "union_upper_ts": set(),
                "gnum_y": ZZ(0), "gnum_rnum": ZZ(0),
            })
            maskbit = 0
            if cm:
                counters["minus_mask_leaves"] += 1
                st["minus_ts"].append(ZZ(t))
                maskbit |= 1
            if cp:
                counters["plus_mask_leaves"] += 1
                st["plus_ts"].append(ZZ(t))
                maskbit |= 2
            if cm and cp:
                counters["intersection_mask_leaves"] += 1
                st["inter_ts"].append(ZZ(t))
            direct_geom_lines.append("%s|%s|%s|%s" % (q, t, p, maskbit))

            a = d-t-1
            xm = False
            if cm:
                assert 0 <= a < p
                xm = (a in zero_sets[p]) and (r in zero_sets[p])

            xp = False
            if cp:
                w = 12*d+t-p
                assert 0 <= w < p
                assert w == ((12*d+t) % p)
                xp = (r in zero_sets[p]) and (w in zero_sets[p])

            bits = (1 if xm else 0) | (2 if xp else 0)
            if xm:
                counters["minus_actual_hits"] += 1
            if xp:
                counters["plus_actual_hits"] += 1
            y = ZZ(1 if bits else 0)
            if y:
                counters["union_actual_hits"] += 1
                key = (q, ZZ(t), p)
                assert key not in direct_bits
                direct_bits[key] = bits
                st["union_ts"].add(ZZ(t))
                if zqt:
                    st["union_upper_ts"].add(ZZ(t))

            coeff_num = N[p]**2 * (ZZ(1 if cm else 0) + ZZ(1 if cp else 0))
            if cm and cp:
                coeff_num -= T3[p]

            At += y
            Rt += QQ(coeff_num, p**2)
            st["gnum_y"] += gnum*y
            st["gnum_rnum"] += gnum*coeff_num

        sumA += At
        if zqt:
            hA += At
        dq_literal += QQ(gnum, q) * (QQ(At) - Rt)

    if row_has_geometry:
        counters["q_geometry_count"] += 1

    event_cov = QQ(0)
    SRq = QQ(0)
    dq_aggregate = QQ(0)
    for p, st in pqstats.items():
        event_cov += QQ(st["gnum_y"], q)
        SRq += QQ(st["gnum_rnum"], q*p**2)
        dq_aggregate += QQ(p**2*st["gnum_y"] - st["gnum_rnum"], q*p**2)
        direct_pq[(q,p)] = st

    assert dq_literal == dq_aggregate
    assert event_cov - SRq == dq_aggregate

    g2num = sum((q*ZZ(1 if t in Zq else 0)-Nq)**2 for t in range(q))
    G_literal = QQ(g2num, q**2)
    G_closed = QQ(Nq*(q-Nq), q)
    assert G_literal == G_closed
    if G_closed == 0:
        assert dq_aggregate == 0
        Pq = QQ(0)
    else:
        Pq = positive_part(dq_aggregate)**2 / G_closed

    rows_direct[q] = {
        "row": Irow,
        "Nq": Nq,
        "G": G_closed,
        "sumA": sumA,
        "hA": hA,
        "event_cov": event_cov,
        "SR": SRq,
        "dq": dq_aggregate,
        "Pq": Pq,
        "geometry_p_count": ZZ(len(pqstats)),
    }

# ---------------------------------------------------------------------------
# Path B1: inverse ordered lower-zero pairs, independently producing actual
# sign events and then Boolean-deduplicating (p,q,t).
# ---------------------------------------------------------------------------
stream_bits = {}
stream_duplicate_same_sign = ZZ(0)
stream_generated_minus = ZZ(0)
stream_generated_plus = ZZ(0)


def add_stream_event(q, t, p, bit):
    global stream_duplicate_same_sign
    key = (ZZ(q), ZZ(t), ZZ(p))
    old = stream_bits.get(key, 0)
    if old & bit:
        stream_duplicate_same_sign += 1
    stream_bits[key] = old | bit


for p in primes:
    Zp = zeros[p]
    # Minus inverse: d=(a+r+1)/7, t=d-a-1, q=p+d.
    for a in Zp:
        for r in Zp:
            s = a+r+1
            if s % 7 != 0:
                continue
            d = ZZ(s//7)
            t = ZZ(d-a-1)
            q = ZZ(p+d)
            if q not in prime_set:
                continue
            cm, cp = masks(p, q, t)
            if cm:
                assert r == 6*d+t
                add_stream_event(q, t, p, 1)
                stream_generated_minus += 1

    # Plus inverse in the unified reflected-upper coordinate:
    # d=(z+w+1)/6, t=p-1-z-6d, q=p+d.
    for z in Zp:
        for w in Zp:
            s = z+w+1
            if s % 6 != 0:
                continue
            d = ZZ(s//6)
            t = ZZ(p-1-z-6*d)
            q = ZZ(p+d)
            if q not in prime_set:
                continue
            cm, cp = masks(p, q, t)
            if cp:
                assert z == p-1-(6*d+t)
                assert w == 12*d+t-p
                add_stream_event(q, t, p, 2)
                stream_generated_plus += 1

# Event comparison and first witness.
event_symmetric_difference = sorted(set(direct_bits) ^ set(stream_bits))
event_bit_differences = sorted(k for k in set(direct_bits) & set(stream_bits)
                               if direct_bits[k] != stream_bits[k])
event_difference_keys = sorted(set(event_symmetric_difference) | set(event_bit_differences))
first_event_mismatch_q = event_difference_keys[0][0] if event_difference_keys else None

# Canonical event hashes.
direct_event_lines = ["%s|%s|%s|%s" % (q,t,p,direct_bits[(q,t,p)])
                      for q,t,p in sorted(direct_bits)]
stream_event_lines = ["%s|%s|%s|%s" % (q,t,p,stream_bits[(q,t,p)])
                      for q,t,p in sorted(stream_bits)]

# Group inverse-stream events by (q,p).
stream_union_ts = {}
for (q,t,p), bits in stream_bits.items():
    stream_union_ts.setdefault((q,p), set()).add(t)

# ---------------------------------------------------------------------------
# Path B2: independent interval endpoints + prefix ranks in the upper zero set.
# This produces the collapsed streamed dq row formula.
# ---------------------------------------------------------------------------
prefix = {}
for q in rows_direct:
    pref = [ZZ(0)]*(q+1)
    running = ZZ(0)
    Zq = zero_sets[q]
    for t in range(q):
        if t in Zq:
            running += 1
        pref[t+1] = running
    assert pref[q] == N[q]
    prefix[q] = pref


def rank_interval(q, I):
    if I is None:
        return ZZ(0)
    return prefix[q][I[1]+1] - prefix[q][I[0]]


rows_stream = {}
interval_mismatches = []
interval_direct_lines = []
interval_rank_lines = []
pq_stream_counter = ZZ(0)
pq_nonempty_counter = ZZ(0)

for q in sorted(rows_direct):
    Nq = N[q]
    Zq = zero_sets[q]
    sumA = ZZ(0)
    hA = ZZ(0)
    dq = QQ(0)
    SRq = QQ(0)
    event_cov = QQ(0)
    geometry_p_count = ZZ(0)

    for p in primes:
        if not (J*(q-p) > T):
            continue
        pq_stream_counter += 1
        Im, Ip = mask_intervals(p, q)
        Ii = intersect_intervals(Im, Ip)

        d_st = direct_pq.get((q,p), None)
        dminus = tuple(d_st["minus_ts"]) if d_st is not None else tuple()
        dplus = tuple(d_st["plus_ts"]) if d_st is not None else tuple()
        dinter = tuple(d_st["inter_ts"]) if d_st is not None else tuple()
        rminus = interval_values(Im)
        rplus = interval_values(Ip)
        rinter = interval_values(Ii)

        # Hash both the directly visited mask sets and the independent endpoint sets.
        for t in sorted(set(dminus) | set(dplus)):
            bit = (1 if t in set(dminus) else 0) | (2 if t in set(dplus) else 0)
            interval_direct_lines.append("%s|%s|%s|%s" % (q,t,p,bit))
        for t in sorted(set(rminus) | set(rplus)):
            bit = (1 if t in set(rminus) else 0) | (2 if t in set(rplus) else 0)
            interval_rank_lines.append("%s|%s|%s|%s" % (q,t,p,bit))

        if dminus != rminus or dplus != rplus or dinter != rinter:
            interval_mismatches.append((q,p,dminus,rminus,dplus,rplus,dinter,rinter))

        if Im is None and Ip is None:
            # No geometric coefficient and, by exact event comparison, no actual event.
            assert not stream_union_ts.get((q,p), set())
            continue

        geometry_p_count += 1
        pq_nonempty_counter += 1
        nm = interval_length(Im)
        np_ = interval_length(Ip)
        ni = interval_length(Ii)
        mm = rank_interval(q, Im)
        mp = rank_interval(q, Ip)
        mi = rank_interval(q, Ii)
        Sm = q*mm - Nq*nm
        Sp = q*mp - Nq*np_
        Si = q*mi - Nq*ni

        uts = stream_union_ts.get((q,p), set())
        mA = ZZ(len(uts))
        hAp = ZZ(sum(1 for t in uts if t in Zq))
        sumA += mA
        hA += hAp

        event_num = q*hAp - Nq*mA
        struct_num = N[p]**2*(Sm+Sp) - T3[p]*Si
        event_cov += QQ(event_num, q)
        SRq += QQ(struct_num, q*p**2)
        dq += QQ(p**2*event_num - struct_num, q*p**2)

        interval_rank_lines.append(
            "META|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s" %
            (q,p,interval_fmt(Im),interval_fmt(Ip),interval_fmt(Ii),
             nm,np_,ni,mm,mp,mi,Sm,Sp,Si,N[p],T3[p],mA)
        )

    assert event_cov - SRq == dq
    Gq = QQ(Nq*(q-Nq), q)
    if Gq == 0:
        assert dq == 0
        Pq = QQ(0)
    else:
        Pq = positive_part(dq)**2 / Gq
    rows_stream[q] = {
        "row": rows_direct[q]["row"],
        "Nq": Nq,
        "G": Gq,
        "sumA": sumA,
        "hA": hA,
        "event_cov": event_cov,
        "SR": SRq,
        "dq": dq,
        "Pq": Pq,
        "geometry_p_count": geometry_p_count,
    }

# Compare every q row and locate the first discrepancy.
row_fields = ("Nq","G","sumA","hA","event_cov","SR","dq","Pq","geometry_p_count")
row_mismatches = []
for q in sorted(rows_direct):
    for field in row_fields:
        if rows_direct[q][field] != rows_stream[q][field]:
            row_mismatches.append((q, field, rows_direct[q][field], rows_stream[q][field]))
first_row_mismatch_q = row_mismatches[0][0] if row_mismatches else None
first_interval_mismatch_q = interval_mismatches[0][0] if interval_mismatches else None
mismatch_qs = [q for q in (first_event_mismatch_q, first_interval_mismatch_q, first_row_mismatch_q)
               if q is not None]
smallest_mismatch_q = min(mismatch_qs) if mismatch_qs else None

# Global exact PAUC and marked-incidence identities.
def global_totals(rows):
    G = sum((r["G"] for r in rows.values()), QQ(0))
    Pplus = sum((r["Pq"] for r in rows.values()), QQ(0))
    Dplus = sum((positive_part(r["dq"]) for r in rows.values()), QQ(0))
    Daff = sum((r["dq"] for r in rows.values()), QQ(0))
    Mmark = sum((r["hA"] for r in rows.values()), ZZ(0))
    Mq0 = sum((QQ(r["Nq"]*r["sumA"], q) for q,r in rows.items()), QQ(0))
    SR = sum((r["SR"] for r in rows.values()), QQ(0))
    residual = QQ(Mmark) - Mq0 - SR - Daff
    return {"G":G,"Pplus":Pplus,"Dplus":Dplus,"Daff":Daff,
            "Mmark":ZZ(Mmark),"Mq0":Mq0,"SR":SR,"residual":residual}

GD = global_totals(rows_direct)
GS = global_totals(rows_stream)
assert GD == GS
assert GD["residual"] == 0

# Canonical row hashes.
def canonical_row_lines(rows):
    out = []
    for q in sorted(rows):
        r = rows[q]
        out.append("%s|%s|%s|%s|%s|%s|%s|%s|%s|%s" %
                   (q,r["Nq"],qfmt(r["G"]),r["sumA"],r["hA"],
                    qfmt(r["event_cov"]),qfmt(r["SR"]),qfmt(r["dq"]),
                    qfmt(r["Pq"]),r["geometry_p_count"]))
    return out

row_lines_direct = canonical_row_lines(rows_direct)
row_lines_stream = canonical_row_lines(rows_stream)

# Select informative rows for human reproduction.
geometry_qs = [q for q,r in rows_direct.items() if r["geometry_p_count"] > 0]
gpositive_geometry_qs = [q for q in geometry_qs if rows_direct[q]["G"] > 0]
union_qs = [q for q,r in rows_direct.items() if r["sumA"] > 0]
nonzero_dq_qs = [q for q,r in rows_direct.items() if r["dq"] != 0]
positive_dq_qs = [q for q,r in rows_direct.items() if r["dq"] > 0]
focus_qs = []
for candidate in [
    min(geometry_qs) if geometry_qs else None,
    min(gpositive_geometry_qs) if gpositive_geometry_qs else None,
    min(union_qs) if union_qs else None,
    min(nonzero_dq_qs) if nonzero_dq_qs else None,
    min(positive_dq_qs) if positive_dq_qs else None,
    ZZ(1129) if ZZ(1129) in rows_direct else None,
]:
    if candidate is not None and candidate not in focus_qs:
        focus_qs.append(candidate)

# Exact counters and hashes.
counters.update({
    "direct_event_keys": ZZ(len(direct_bits)),
    "stream_event_keys": ZZ(len(stream_bits)),
    "stream_generated_minus": stream_generated_minus,
    "stream_generated_plus": stream_generated_plus,
    "stream_duplicate_same_sign": stream_duplicate_same_sign,
    "event_symmetric_difference": ZZ(len(event_symmetric_difference)),
    "event_bit_differences": ZZ(len(event_bit_differences)),
    "interval_mismatches": ZZ(len(interval_mismatches)),
    "row_mismatches": ZZ(len(row_mismatches)),
    "pq_stream_static": pq_stream_counter,
    "pq_nonempty_geometry": pq_nonempty_counter,
    "q_G_positive": ZZ(sum(1 for r in rows_direct.values() if r["G"] > 0)),
    "q_dq_nonzero": ZZ(len(nonzero_dq_qs)),
    "q_dq_positive": ZZ(len(positive_dq_qs)),
    "q_actual_union": ZZ(len(union_qs)),
})

hashes = {
    "zero_sha256": sha256_lines(zero_lines),
    "t3_sha256": sha256_lines(t3_lines),
    "direct_geometry_sha256": sha256_lines(sorted(direct_geom_lines)),
    "interval_geometry_sha256": sha256_lines(sorted(interval_rank_lines)),
    "direct_interval_expansion_sha256": sha256_lines(sorted(interval_direct_lines)),
    "direct_event_sha256": sha256_lines(direct_event_lines),
    "stream_event_sha256": sha256_lines(stream_event_lines),
    "direct_row_sha256": sha256_lines(row_lines_direct),
    "stream_row_sha256": sha256_lines(row_lines_stream),
}

# direct_geom_lines and interval_direct_lines are the same directly visited masks.
assert hashes["direct_geometry_sha256"] == hashes["direct_interval_expansion_sha256"]
assert hashes["direct_event_sha256"] == hashes["stream_event_sha256"]
assert hashes["direct_row_sha256"] == hashes["stream_row_sha256"]

out = []
out.append("Q7503_PAUC_REPRODUCTION_V1")
out.append("PARAM T=%s J=%s A_T=%s B_T=%s H_T=%s" % (T,J,A_T,B_T,H_T))
out.append("SMALLEST_MISMATCH_Q %s" % ("NONE" if smallest_mismatch_q is None else smallest_mismatch_q))
for k in sorted(counters):
    out.append("COUNTER %s=%s" % (k,counters[k]))
for k in sorted(hashes):
    out.append("HASH %s=%s" % (k,hashes[k]))

out.append("GLOBAL G=%s" % qfmt(GD["G"]))
out.append("GLOBAL Pplus=%s" % qfmt(GD["Pplus"]))
out.append("GLOBAL Dplus=%s" % qfmt(GD["Dplus"]))
out.append("GLOBAL Daff=%s" % qfmt(GD["Daff"]))
out.append("IDENTITY Mmark=%s Mq0=%s SR=%s Daff=%s residual=%s" %
           (GD["Mmark"],qfmt(GD["Mq0"]),qfmt(GD["SR"]),
            qfmt(GD["Daff"]),qfmt(GD["residual"])))
out.append("FOCUS_QS %s" % (",".join(map(str,focus_qs)) if focus_qs else "NONE"))

for q in focus_qs:
    r = rows_direct[q]
    geom = []
    for p in primes:
        Im, Ip = mask_intervals(p,q)
        if Im is None and Ip is None:
            continue
        Ii = intersect_intervals(Im,Ip)
        geom.append("p=%s:Np=%s:T3=%s:Im=%s:Ip=%s:Ii=%s" %
                    (p,N[p],T3[p],interval_fmt(Im),interval_fmt(Ip),interval_fmt(Ii)))
    out.append("ROW q=%s row=%s Zq=%s Nq=%s Gq=%s sumA=%s hA=%s event_cov=%s SR=%s dq=%s Pplus_q=%s geometry_p_count=%s" %
               (q,interval_fmt(r["row"]),",".join(map(str,zeros[q])),r["Nq"],
                qfmt(r["G"]),r["sumA"],r["hA"],qfmt(r["event_cov"]),
                qfmt(r["SR"]),qfmt(r["dq"]),qfmt(r["Pq"]),r["geometry_p_count"]))
    out.append("ROW_GEOMETRY q=%s %s" % (q,";".join(geom)))

out.append("SECONDS %.6f" % (time.time()-t0))

with open("q7503_pauc_reproduction.txt","w") as fh:
    fh.write("\n".join(out) + "\n")

print("\n".join(out))
