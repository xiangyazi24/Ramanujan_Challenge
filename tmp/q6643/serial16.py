"""SERIALIZER primitive, cure D of Q6629 (FABLE-EXEC-15): 16-site framed line, sites 0..15, frame 15 | 0 = the region's own
arc boundary (oracle field: roles H I^14 F).  No ring seam: the pair (15, 0) is clock row D0 (stutter) and is never used.
Static labels (fan.py): lev(t) = v2(t + 1) + 1, ct(t) = 1 iff t is a right child; demand(t) = environment request at host t
(persistent at 7 and 9).  Served move: left child c -> parent(c) (to the right).
FRONT = the privilege of ONE row-90 clock on the frame.  It never waits and never leaves the line: an excursion is a
bounded sequence of whole passes of the front, so no walker exists apart from the front.
One job = passes  IDLE (select: first demanding left-child host with turn = 0 claims: open := origin := 1)
                  INIT (cand := 1 everywhere; verify exactly one origin, equal to the open host, no mailbox)
                  PEEL* (delete every other remaining candidate, starting with the first: pass k deletes the sites of
                         level exactly k; when the origin is deleted, the NEXT surviving candidate is its parent: DELIVER)
                  DRAIN (clear cand / origin / mailbox everywhere; ACK and close the open host; turn := 1)
Level comparison is done by the frame (alternation from the frame start), not by a token: no OPEN stack, no walk back.
Per site: clock (x, u) 4; up header (mode 7, par 2, hit 2, found 2, any 2, oc 3) 336; down header 7;
marks (cand, origin, open, dl, turn) 32.  Alphabet 4 * 336 * 7 * 32 = 301,056.
Switches: 'ALLOW-NESTED-BEGIN', 'NO-DRAIN', 'NOHEAD'."""
import itertools, random, sys
N = 16
v2 = lambda x: (x & -x).bit_length() - 1
lev = lambda t: v2(t + 1) + 1
ct = lambda t: ((t + 1) >> lev(t)) & 1
parent = lambda t: t - (1 << (lev(t) - 1)) if ct(t) else t + (1 << (lev(t) - 1))
DEMAND = {7, 9}
IDLE, CLAIMED, DIRTY, INIT, PEEL, DRAINA, DRAINF = range(7); MODES = range(7)
NONE, START, PEELGO, ACKDRAIN, FAILDRAIN, FREE, RESETTURN = range(7)
HDR = [(m, p, h, f, a, oc) for m in MODES for p in (0, 1) for h in (0, 1) for f in (0, 1) for a in (0, 1) for oc in (0, 1, 2)]
MARK = list(itertools.product((0, 1), repeat=5))
def role(t, sw=()): return 'H' if t == 0 else ('I' if 'NOHEAD' in sw else 'F') if t == N - 1 else 'I'
def up(u, g): return 1 if g == 'H' else 0 if g == 'F' else u
def proc(h, mk, t, sw):
    m, par, hit, found, anyc, oc = h; cand, origin, opn, dl, turn = mk; ev = []
    eligible = t in DEMAND and ct(t) == 0
    dirty = (cand or origin or opn or dl) and 'NO-DRAIN' not in sw
    if m == IDLE:
        if dirty: m = DIRTY
        elif eligible and not turn and not opn: opn = origin = 1; m = CLAIMED; ev.append(('BEGIN', t))
    elif m == CLAIMED:
        if 'ALLOW-NESTED-BEGIN' in sw and eligible and not opn: opn = origin = 1; ev.append(('BEGIN', t))
        elif dirty and 'ALLOW-NESTED-BEGIN' not in sw: m = DIRTY
    elif m == INIT:
        cand = 1
        if origin: oc = min(2, oc + 1)
        if (opn != origin or dl) and not ({'NO-DRAIN', 'ALLOW-NESTED-BEGIN'} & set(sw)): oc = 2
    elif m == PEEL:
        if cand:
            anyc = 1
            if par == 0:
                cand = 0; par = 1
                if origin: hit = 1
            else:
                par = 0
                if hit: dl = 1; found = 1; hit = 0; ev.append(('DELIVER', t))
    elif m in (DRAINA, DRAINF):
        if opn:
            if m == DRAINA and origin: ev.append(('ACK', t))
            opn = 0
            if m == DRAINA: turn = 1
        if 'NO-DRAIN' not in sw:
            if dl: ev.append(('CONSUME', t))
            cand = origin = dl = 0
    return (m, par, hit, found, anyc, oc), (cand, origin, opn, dl, turn), ev
def top(h, sw):
    m, par, hit, found, anyc, oc = h
    if m == IDLE: return RESETTURN
    if m == CLAIMED: return START
    if m == DIRTY: return FAILDRAIN
    if m == INIT: return PEELGO if (oc == 1 or (oc == 2 and {'NO-DRAIN', 'ALLOW-NESTED-BEGIN'} & set(sw))) else FAILDRAIN
    if m == PEEL: return ACKDRAIN if found else PEELGO if anyc else FAILDRAIN
    return FREE
BOTTOM = {NONE: IDLE, START: INIT, PEELGO: PEEL, ACKDRAIN: DRAINA, FAILDRAIN: DRAINF, FREE: IDLE, RESETTURN: IDLE}
def pair(i, a, b, sw=()):
    ga, gb = role(i, sw), role(i + 1, sw); ax, au, ah, ad, am = a; bx, bu, bh, bd, bm = b
    if ax != bx:
        h, mk, ev = proc(ah, bm, i + 1, sw); d = bd
        if gb == 'F':
            d = top(h, sw); ev.append(('TOP', d))
            if d == RESETTURN: mk = mk[:4] + (0,)
        return a, (ax, 1 if gb == 'I' else bu, h, d, mk), ev
    if up(au, ga) == 1 and up(bu, gb) == 0:
        d = bd; mk = am
        if d == RESETTURN: mk = mk[:4] + (0,)
        if ga == 'H':
            h, mk, ev = proc((BOTTOM[d], 0, 0, 0, 0, 0), mk, 0, sw)
            return (1 - ax, au, h, NONE, mk), b, ev + [('BOTTOM', BOTTOM[d])]
        return (ax, 0, ah, d, mk), b, []
    return a, b, []
SITE_COUNT = 4 * len(HDR) * 7 * len(MARK)
def totality(stride=1):
    S = set(HDR); M = set(MARK); c = 0
    for sw in ((), ('ALLOW-NESTED-BEGIN',), ('NO-DRAIN',), ('NOHEAD',)):
        for h in HDR:
            assert top(h, sw) in range(7); c += 1
            for mk in MARK:
                for t in range(N):
                    h2, m2, _ = proc(h, mk, t, sw); c += 1; assert h2 in S and m2 in M
        for i in range(N - 1):
            for ax, au, bx, bu in itertools.product((0, 1), repeat=4):
                for d in range(7):
                    for h in HDR[::stride]:
                        for mk in MARK[::5]:
                            a2, b2, _ = pair(i, (ax, au, h, d, mk), (bx, bu, h, d, mk), sw); c += 1
                            for s in (a2, b2): assert s[0] in (0, 1) and s[1] in (0, 1) and s[2] in S and s[3] in range(7) and s[4] in M
    return c
def fronts(X, sw=()):
    out = []
    for k in range(N - 1):
        a, b = X[k], X[k + 1]
        if a[0] != b[0]: out.append(('UP', k, a[2]))
        elif up(a[1], role(k, sw)) == 1 and up(b[1], role(k + 1, sw)) == 0: out.append(('DOWN', k + 1, b[3]))
    return out
def ETB(X, sw=()):
    F = fronts(X, sw); bad = []
    opn = [t for t in range(N) if X[t][4][2]]; org = [t for t in range(N) if X[t][4][1]]
    cand = [t for t in range(N) if X[t][4][0]]; dl = [t for t in range(N) if X[t][4][3]]
    if len(F) != 1 or len(opn) > 1: return ['E1 single active owner']
    kind, pos, hd = F[0]
    draining = (kind == 'UP' and hd[0] in (DRAINA, DRAINF)) or (kind == 'DOWN' and hd in (ACKDRAIN, FAILDRAIN))
    if draining:
        if kind == 'UP' and any(t <= pos for t in org + cand + dl + opn): bad.append('E5 drain before reuse')
        return bad
    idle = (kind == 'UP' and hd[0] == IDLE) or (kind == 'DOWN' and hd in (NONE, FREE, RESETTURN))
    if (org or cand or dl) and not opn: bad.append('E2 connected ownership witness')
    if opn and org != opn: bad.append('E2 connected ownership witness')
    if any(t not in DEMAND or ct(t) for t in opn): bad.append('E3 stationary waiting jobs')
    if dl and (len(dl) > 1 or not opn or dl[0] != parent(opn[0])): bad.append('E4 destination binding')
    if idle and (opn or org or cand or dl): bad.append('E5 drain before reuse')
    if kind == 'UP' and hd[0] == CLAIMED and (cand or dl or len(opn) != 1 or opn[0] > pos): bad.append('E5 drain before reuse')
    return bad
def step(X, i, sw):
    a, b, ev = pair(i, X[i], X[i + 1], sw); X[i], X[i + 1] = a, b; return ev
def garbage(rnd): return [(rnd.randrange(2), rnd.randrange(2), rnd.choice(HDR), rnd.randrange(7), rnd.choice(MARK)) for _ in range(N)]
def clean():
    z = (0, 0, 0, 0, 0)
    return [(0, 0, (IDLE, 0, 0, 0, 0, 0), NONE, z) for _ in range(N)]
def run(steps, seed, sw=(), X=None):
    rnd = random.Random(seed); X = X or garbage(rnd)
    last_bad = -1; acks = {}; wrong_del = [0, 0]; wrong_ack = [0, 0]; dels = [0, 0]; log = []; first = {}
    for t in range(steps):
        i = rnd.randrange(N - 1); pre_open = [s for s in range(N) if X[s][4][2]]
        ev = step(X, i, sw)
        if ETB(X, sw): last_bad = t
        late = int(t > steps // 2)
        for e in ev:
            if e[0] == 'DELIVER':
                ok = len(pre_open) == 1 and parent(pre_open[0]) == e[1]; dels[late] += 1; wrong_del[late] += not ok
                if not ok and len(log) < 4: log.append(('DELIVER at', e[1], 'open', pre_open, 'origin', [s for s in range(N) if X[s][4][1]]))
            elif e[0] == 'ACK':
                dl = [s for s in range(N) if X[s][4][3]]; ok = dl == [parent(e[1])]
                wrong_ack[late] += not ok
                if ok: acks[e[1]] = acks.get(e[1], 0) + late; first.setdefault(e[1], t)
                elif len(log) < 4: log.append(('ACK at', e[1], 'mailboxes', dl))
    return dict(etb_last_violation=last_bad, fronts_final=len(fronts(X, sw)), late_good_acks=acks, deliveries_early_late=dels,
                wrong_deliveries_early_late=wrong_del, wrong_acks_early_late=wrong_ack, witnesses=log)
def headless(sw=('NOHEAD',)):
    X = clean(); X = [(0, 1, s[2], s[3], s[4]) for s in X]
    fixed = all(pair(i, X[i], X[i + 1], sw)[:2] == (X[i], X[i + 1]) for i in range(N - 1))
    return dict(fronts=len(fronts(X, sw)), fixed_point=fixed, pending_demands=sorted(DEMAND))
