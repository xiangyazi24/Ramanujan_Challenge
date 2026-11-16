import random
import serial16 as S
import etbstar as E
N = 32; REG = [(0, 0), (0, 1), (1, 0), (1, 1)]; NAME = {(0, 0): 'G00', (0, 1): 'G01', (1, 0): 'G10', (1, 1): 'G11'}
BASE_DEMAND = {7, 9}
def loc(p, i): q = (i - 8 * p) % N; return q // 16, q % 16
def glob(p, r, l): return (16 * r + l + 8 * p) % N
def region_slice(X, p, r): return [X[glob(p, r, l)][p] for l in range(16)]
def lane_sw(p, r, sw):
    out = []
    if 'HEADLESS-REGION' in sw and (p, r) == (0, 1): out.append('NOHEAD')
    if 'NO-SEAM-DRAIN' in sw and p == 1: out.append('NO-DRAIN')
    return tuple(out)
def pair32(X, i, sw=()):
    j = (i + 1) % N; ev = []
    for p in (0, 1):
        r, l = loc(p, i)
        if l == 15: continue
        a, b = X[i][p], X[j][p]
        if 'SHARED-LANE' in sw and p == 1: a = a[:4] + (X[i][0][4],); b = b[:4] + (X[j][0][4],)
        dem = set(BASE_DEMAND)
        if 'NO-MIDPOINT-OWNER' in sw and p == 0: dem |= {1, 15}
        if 'REMOTE-BUSY' in sw and any(s[4][2] for s in region_slice(X, p, 1 - r)): dem = set()
        S.DEMAND = dem
        a2, b2, e = S.pair(l, a, b, lane_sw(p, r, sw)); S.DEMAND = set(BASE_DEMAND)
        if 'SHARED-LANE' in sw and p == 1:
            X[i] = (X[i][0][:4] + (a2[4],), X[i][1]); X[j] = (X[j][0][:4] + (b2[4],), X[j][1])
        X[i] = X[i][:p] + (a2,) + X[i][p + 1:]; X[j] = X[j][:p] + (b2,) + X[j][p + 1:]
        ev += [(p, r) + x for x in e]
    return ev
def totality():
    fb = {0: [], 1: []}
    for i in range(N):
        for p in (0, 1):
            r, l = loc(p, i); r2, l2 = loc(p, (i + 1) % N)
            if l == 15: fb[p].append(i); assert l2 == 0 and r2 != r
            else: assert (r2, l2) == (r, l + 1)
    assert fb == {0: [15, 31], 1: [7, 23]}
    c16 = S.totality(stride=7); rnd = random.Random(3); c = 0; HS = set(S.HDR); MS = set(S.MARK)
    for sw in ((), ('SHARED-LANE',), ('NO-SEAM-DRAIN',), ('NO-MIDPOINT-OWNER',), ('REMOTE-BUSY',), ('HEADLESS-REGION',)):
        for _ in range(4000):
            X = garbage(rnd); pair32(X, rnd.randrange(N), sw); c += 1
            for s in X:
                for q in s: assert q[0] in (0, 1) and q[1] in (0, 1) and q[2] in HS and q[3] in range(7) and q[4] in MS
    return dict(frame_boundary_pairs=fb, serial16_factor_inputs=c16, product_samples=c)
def garbage(rnd): return [tuple((rnd.randrange(2), rnd.randrange(2), rnd.choice(S.HDR), rnd.randrange(7), rnd.choice(S.MARK)) for _ in (0, 1)) for _ in range(N)]
def clean(): z = (0, 0, (S.IDLE, 0, 0, 0, 0, 0), S.NONE, (0, 0, 0, 0, 0)); return [(z, z) for _ in range(N)]
def RETB(X, sw=()):
    out = {}
    for p, r in REG:
        Y = region_slice(X, p, r); bad = []; F = S.fronts(Y, lane_sw(p, r, sw))
        opn = [t for t in range(16) if Y[t][4][2]]; dl = [t for t in range(16) if Y[t][4][3]]
        if [S.role(l) for l in range(16)] != ['H'] + ['I'] * 14 + ['F'] or loc(p, glob(p, r, 0))[1] != 0: bad.append('R0 current bounded frame')
        if len(F) != 1: bad.append('R1 local front')
        if len(opn) > 1: bad.append('R2 one owner per regional namespace')
        if len(F) == 1 and len(opn) <= 1 and not E.ETBstar(Y): bad.append('R3 connected present ownership (ETB*)')
        if any(t not in BASE_DEMAND for t in opn): bad.append('R4 static midpoint assignment')
        if any(glob(p, r, t) != glob(p, r, t) or loc(1 - p, glob(p, r, t))[1] < 0 for t in dl): bad.append('R5 bounded seam ownership')
        drain = len(F) == 1 and ((F[0][0] == 'UP' and F[0][2][0] in (S.DRAINA, S.DRAINF)) or (F[0][0] == 'DOWN' and F[0][2] in (S.ACKDRAIN, S.FAILDRAIN)))
        if dl and not drain and (len(opn) != 1 or dl != [S.parent(opn[0])]): bad.append('R6 destination binding')
        if len(F) == 1:
            kind, pos, hd = F[0]; idle = (kind == 'UP' and hd[0] == S.IDLE) or (kind == 'DOWN' and hd in (S.NONE, S.FREE, S.RESETTURN))
            if idle and any(any(Y[t][4][:4]) for t in range(16)): bad.append('R7 drain before regional reuse')
        out[NAME[(p, r)]] = bad
    return out
def run(steps, seed, sw=(), X=None, sample=97):
    rnd = random.Random(seed); X = X or garbage(rnd); acks = {}; wd = [0, 0]; wa = [0, 0]; dels = [0, 0]; begins = {}
    last_bad = {n: -1 for n in NAME.values()}; simul = simul_same = 0; wit = []
    for t in range(steps):
        i = rnd.randrange(N); pre = {(p, r): [l for l in range(16) if X[glob(p, r, l)][p][4][2]] for p, r in REG}
        ev = pair32(X, i, sw); late = int(t > steps // 2)
        for p, r, kind, *rest in ev:
            nm = NAME[(p, r)]
            if kind == 'BEGIN': begins[(nm, glob(p, r, rest[0]))] = begins.get((nm, glob(p, r, rest[0])), 0) + late
            elif kind == 'DELIVER':
                ok = len(pre[(p, r)]) == 1 and S.parent(pre[(p, r)][0]) == rest[0] and pre[(p, r)][0] in BASE_DEMAND
                dels[late] += 1; wd[late] += not ok
                if not ok and late and len(wit) < 3: wit.append((nm, 'DELIVER at global', glob(p, r, rest[0]), 'open local', pre[(p, r)]))
            elif kind == 'ACK':
                Y = region_slice(X, p, r); dl = [l for l in range(16) if Y[l][4][3]]; ok = dl == [S.parent(rest[0])] and rest[0] in BASE_DEMAND
                wa[late] += not ok
                if ok: acks[(nm, glob(p, r, rest[0]))] = acks.get((nm, glob(p, r, rest[0])), 0) + late
                elif late and len(wit) < 3: wit.append((nm, 'ACK at global', glob(p, r, rest[0]), 'mailboxes local', dl))
        if t % sample == 0:
            V = RETB(X, sw)
            for nm, b in V.items():
                if b: last_bad[nm] = t
            act = [(p, r) for p, r in REG if not V[NAME[(p, r)]] and any(s[4][2] for s in region_slice(X, p, r))]
            simul += len(act) >= 2; simul_same += any(a[0] == b[0] for a in act for b in act if a != b)
    dup = {k: v for k, v in begins.items() if k[0] in ('G00', 'G01') and (k[1] % 16) in (1, 15)}
    return dict(retb_last_violation=last_bad, late_good_acks=dict(sorted(acks.items())), deliveries=dels, wrong_deliveries_early_late=wd,
                wrong_acks_early_late=wa, samples_with_2_active_clean_regions=simul, samples_with_2_active_in_one_cover=simul_same,
                late_begins_of_shifted_jobs_in_base_cover=dup, witnesses=wit)
