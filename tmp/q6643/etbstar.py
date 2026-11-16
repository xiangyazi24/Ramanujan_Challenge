import random, itertools
import serial16 as S
N = S.N; lev, ct, parent = S.lev, S.ct, S.parent
HOSTS = (7, 9)
LVL = lambda k: frozenset(t for t in range(N) if lev(t) >= k)
def first_eligible(c, t7, t9): return (c == 7 and t7 == 0) or (c == 9 and t7 == 1 and t9 == 0)
def canon_all():
    E = frozenset(); Z = (0, 0, 0, 0, 0)
    for t7, t9 in itertools.product((0, 1), repeat=2):
        for p in range(N - 1):
            if all(t for c, t in ((7, t7), (9, t9)) if c <= p):
                yield ('UP', p, (S.IDLE,) + Z), (E, E, E, E), (t7, t9)
        for p in range(1, N):
            for d in (S.NONE, S.FREE): yield ('DOWN', p, d), (E, E, E, E), (t7, t9)
            if all((t == 0) == (c >= p) for c, t in ((7, t7), (9, t9))): yield ('DOWN', p, S.RESETTURN), (E, E, E, E), (t7, t9)
        for c in HOSTS:
            if not first_eligible(c, t7, t9): continue
            O = frozenset([c]); P = parent(c); L = lev(c)
            for p in range(c, N - 1): yield ('UP', p, (S.CLAIMED,) + Z), (E, O, O, E), (t7, t9)
            for p in range(1, N): yield ('DOWN', p, S.START), (E, O, O, E), (t7, t9)
            for p in range(N - 1): yield ('UP', p, (S.INIT, 0, 0, 0, 0, int(c <= p))), (frozenset(range(p + 1)), O, O, E), (t7, t9)
            for k in range(1, L + 1):
                A, B = LVL(k), LVL(k + 1)
                for p in range(1, N): yield ('DOWN', p, S.PEELGO), (A, O, O, E), (t7, t9)
                for p in range(N - 1):
                    seen = [t for t in A if t <= p]; last = k == L
                    hdr = (S.PEEL, len(seen) % 2, int(last and c <= p < P), int(last and P <= p), int(bool(seen)), 0)
                    cand = frozenset(t for t in B if t <= p) | frozenset(t for t in A if t > p)
                    yield ('UP', p, hdr), (cand, O, O, frozenset([P]) if last and P <= p else E), (t7, t9)
            B = LVL(L + 1); D = frozenset([P])
            for p in range(1, N): yield ('DOWN', p, S.ACKDRAIN), (B, O, O, D), (t7, t9)
            for p in range(N - 1):
                tt = (t7, t9) if c > p else ((1, t9) if c == 7 else (t7, 1))
                yield ('UP', p, (S.DRAINA,) + Z), (frozenset(t for t in B if t > p), O if c > p else E, O if c > p else E, D if P > p else E), tt
CANON = set(canon_all())
def project(X):
    F = S.fronts(X)
    if len(F) != 1: return None
    sets = tuple(frozenset(t for t in range(N) if X[t][4][j]) for j in range(4))
    return F[0], sets, (X[7][4][4], X[9][4][4])
def ETBstar(X): return project(X) in CANON
def build(st, rnd=None):
    (kind, p, hd), (cand, org, opn, dl), (t7, t9) = st; X = []
    for t in range(N):
        if kind == 'UP': x, u = (1, 1) if t <= p else (0, 0)
        else: x, u = 0, int(t < p)
        h = rnd.choice(S.HDR) if rnd else (S.IDLE, 0, 0, 0, 0, 0); d = rnd.randrange(7) if rnd else S.NONE
        turn = t7 if t == 7 else t9 if t == 9 else (rnd.randrange(2) if rnd else 0)
        if kind == 'UP' and t == p: h = hd
        if kind == 'DOWN' and t == p: d = hd
        X.append((x, u, h, d, (int(t in cand), int(t in org), int(t in opn), int(t in dl), turn)))
    return X
def check(variants=3, seed=1):
    rnd = random.Random(seed); steps = bad_self = bad_close = bad_stut = bad_del = bad_ack = dels = acks = 0
    for st in CANON:
        for v in range(variants):
            X = build(st, rnd if v else None)
            if project(X) != st: bad_self += 1; continue
            kind, p, _ = st[0]; e = p if kind == 'UP' else p - 1
            for i in range(N - 1):
                if i != e and S.pair(i, X[i], X[i + 1])[:2] != (X[i], X[i + 1]): bad_stut += 1
            opn = sorted(st[1][2]); a, b, ev = S.pair(e, X[e], X[e + 1]); X[e], X[e + 1] = a, b; steps += 1
            if not ETBstar(X): bad_close += 1
            for x in ev:
                if x[0] == 'DELIVER': dels += 1; bad_del += not (len(opn) == 1 and parent(opn[0]) == x[1])
                if x[0] == 'ACK': acks += 1; bad_ack += [t for t in range(N) if X[t][4][3]] != [parent(x[1])]
    return dict(canonical_states=len(CANON), one_step_checks=steps, self_failures=bad_self, closure_failures=bad_close,
                non_front_pair_moves=bad_stut, deliveries=dels, wrong_deliveries=bad_del, acks=acks, wrong_acks=bad_ack)
def U(p, h):
    X = [((1, 1) if i <= p else (0, 0)) + ((S.IDLE, 0, 0, 0, 0, 0), S.NONE, (0, 0, 0, 0, 0)) for i in range(N)]
    X[p] = X[p][:2] + (h,) + X[p][3:]; return X
def setm(X, t, mk): X[t] = X[t][:4] + (mk,)
def instances():
    r = {}
    X = S.clean(); X[1] = X[1][:3] + (S.START,) + X[1][4:]; r['Q6633 2.3'] = (S.ETB(X) == [], ETBstar(X))
    X = U(8, (S.PEEL, 1, 1, 0, 0, 0)); setm(X, 7, (0, 1, 1, 0, 0)); setm(X, 9, (1, 0, 0, 0, 0)); r['Q6633 3.1'] = (S.ETB(X) == [], ETBstar(X))
    X = U(6, (S.DRAINA, 0, 0, 0, 0, 0)); setm(X, 7, (0, 1, 1, 0, 0)); r['Q6633 3.2'] = (S.ETB(X) == [], ETBstar(X))
    setm(X, 15, (0, 0, 0, 1, 0)); r['Q6633 3.3'] = (S.ETB(X) == [], ETBstar(X))
    setm(X, 15, (1, 0, 0, 1, 0)); r['canonical mid-drain copy'] = (S.ETB(X) == [], ETBstar(X))
    return r
if __name__ == '__main__':
    print('check:', check(), flush=True)
    print('instances:', instances(), flush=True)
    print('canonical states violating old ETB:', sum(1 for st in CANON if S.ETB(build(st))), flush=True)
