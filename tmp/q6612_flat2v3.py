"""FLAT-2SCALE v3 (FABLE-EXEC-12): local pair-interval arcs only, no seam, no births. Ring n = 64.
Arc family (s, p), s in {3, 4}, p in {0, 1}: pair intervals [B, B + 2r), r = 2^s, B = p r mod 2r; each is an H I^(2r-2) F arc
(H at B, F at B + 2r - 1) carrying its own row-90 Dijkstra clock; its privilege is the service head for the block pair
(B, B + r). Evidence = payload of the upward pass (hop-copied): val = E_s at the sel-th address cell of the left block,
compared with A at the sel-th address cell of the right block; sel advances at each bottom move.
Site state: one factor (x, u, sel, val, bad) per family (dense: four factors) + R.
Alphabet: 48^2 * 72^2 * 3 = 35,831,808 (scale 3: sel in {0,1}: 48 states; scale 4: sel in {0,1,2}: 72 states).
Switches: 'NOCHECKPOINT', 'LOCKSTEP', 'TERMREUSE', 'NODRAIN', 'SHAREDCLOCK'."""
import itertools, random, sys
N = 64; FAM = [(3, 0), (3, 1), (4, 0), (4, 1)]; OFF = {3: (3, 6), 4: (6, 12, 13)}
CS = {s: [(x, u, sel, v, bd) for x in (0, 1) for u in (0, 1) for sel in range(len(OFF[s])) for v in ('N', 0, 1) for bd in (0, 1)] for s in (3, 4)}
RS = [None, 'ACC', 'REJ']
def up(c, g): return 1 if g == 'H' else 0 if g == 'F' else c[1]
def clock(s, a, b, ga, gb, ob, coin, sw):
    if not (ga in 'HI' and gb in 'IF'): return a, b, None
    ax, au, asel, av, abd = a; bx, bu, bsel, bv, bbd = b; w = len(OFF[s])
    if ax != bx:
        sel, v, bd = asel, av, abd
        if 'NODRAIN' in sw and coin: sel = (sel + 1) % w
        key = bsel if 'LOCKSTEP' in sw else sel
        if ob['idx'][s] is not None and ob['idx'][s] == key:
            if v != 'N' and v != ob['A']: bd = 1
            v = ob['E'][s]
        na = (ax, au, asel, 'N', abd) if 'NOCHECKPOINT' in sw else a
        return na, (ax, 1 if gb == 'I' else bu, sel, v, bd), (('TOP', sel, bd) if gb == 'F' else ('UP',))
    if up(a, ga) == 1 and up(b, gb) == 0:
        if ga == 'H': return (1 - ax, au, (asel + 1) % w, 'N', 0), b, ('BOTTOM',)
        return (ax, 0, asel, av, abd), b, ('DOWN',)
    return a, b, None
def pair(a, b, oa, ob, coin=0, sw=()):
    ev = []
    if a[4] is not None:
        ev.append(('PUBLISH', a[4])); return a[:4] + ((a[4] if 'TERMREUSE' in sw else None),), b, ev
    na, nb, R = list(a[:4]), list(b[:4]), b[4]
    for k, (s, p) in enumerate(FAM):
        slot = 0 if ('SHAREDCLOCK' in sw and k == 1) else k
        na[slot], nb[slot], e = clock(s, na[slot], nb[slot], oa['g'][k], ob['g'][k], ob, coin, sw)
        if e:
            ev.append((k,) + e)
            if e[0] == 'TOP': R = 'REJ' if (e[2] or R == 'REJ') else 'ACC'
    return tuple(na) + (a[4],), tuple(nb) + (R,), ev
def oracle_ring(digits):
    n = len(digits); O = []
    for i in range(n):
        o = dict(g=[], idx={}, E={}, A=digits[i])
        for s, p in FAM:
            r = 1 << s; q = (i - p * r) % (2 * r); o['g'].append('H' if q == 0 else 'F' if q == 2 * r - 1 else 'I')
        for s in (3, 4):
            r = 1 << s; off = i % r; B = i - off
            o['idx'][s] = OFF[s].index(off) if off in OFF[s] else None
            o['E'][s] = digits[i] ^ int(all(digits[(B + x) % n] == 1 for x in OFF[s] if x < off))
        O.append(o)
    return O
def reference(digits, O, s, top, sel):
    r = 1 << s; n = len(digits); B = (top - (2 * r - 1)) % n; x = OFF[s][sel]
    return int(O[(B + x) % n]['E'][s] != digits[(B + r + x) % n])
def arc_privileges(X, O, k):
    s, p = FAM[k]; r = 1 << s; n = len(X); out = []
    for B in range(p * r, n + p * r, 2 * r):
        c = 0
        for t in range(2 * r - 1):
            i = (B + t) % n; j = (i + 1) % n; a, b = X[i][k], X[j][k]; ga, gb = O[i]['g'][k], O[j]['g'][k]
            c += int(a[0] != b[0] or (up(a, ga) == 1 and up(b, gb) == 0))
        out.append(c)
    return out
def totality():
    cnt = 0
    for s in (3, 4):
        S = set(CS[s])
        for ga, gb in itertools.product('HIF', repeat=2):
            for idx in [None] + list(range(len(OFF[s]))):
                for A, E in itertools.product((0, 1), repeat=2):
                    ob = dict(idx={s: idx}, E={s: E}, A=A)
                    for sw, coin in (((), 0), (('NOCHECKPOINT',), 0), (('LOCKSTEP',), 0), (('NODRAIN',), 0), (('NODRAIN',), 1)):
                        for a in CS[s]:
                            for b in CS[s]:
                                x, y, _ = clock(s, a, b, ga, gb, ob, coin, sw); cnt += 1
                                assert x in S and y in S
    return cnt
def run(digits, steps, seed, sw=(), warm=0.5):
    rnd = random.Random(seed); n = len(digits); O = oracle_ring(digits)
    X = [tuple(rnd.choice(CS[s]) for s, _ in FAM) + (rnd.choice(RS),) for _ in range(n)]
    minpriv = 99; tops = wrong = pubs = late_pubs = 0; seen = set()
    for t in range(steps):
        i = rnd.randrange(n); j = (i + 1) % n
        X[i], X[j], ev = pair(X[i], X[j], O[i], O[j], rnd.randrange(2), sw)
        if t % 211 == 0: minpriv = min([minpriv] + [mm for k in range(4) for mm in arc_privileges(X, O, k)])
        for e in ev:
            if e[0] == 'PUBLISH': pubs += 1; late_pubs += t > warm * steps
            elif e[1] == 'TOP' and t > warm * steps:
                s = FAM[e[0]][0]; tops += 1; wrong += e[3] != reference(digits, O, s, j, e[2]); seen.add((e[0], j, e[2]))
    fin = [mm for k in range(4) for mm in arc_privileges(X, O, k)]
    return dict(min_priv=minpriv, final_priv=(min(fin), max(fin)), late_passes=tops, late_wrong=wrong, jobs_seen=len(seen), pubs=pubs, late_pubs=late_pubs)
if __name__ == '__main__':
    print('alphabet', len(CS[3]), '^2 x', len(CS[4]), '^2 x 3 =', len(CS[3]) ** 2 * len(CS[4]) ** 2 * 3, '| factor totality (exhaustive per factor):', totality(), 'inputs ok', flush=True)
    rnd = random.Random(7); dg = [rnd.randrange(2) for _ in range(N)]
    print('jobs in total (arc, item):', 2 * 4 * 2 + 2 * 2 * 3, flush=True)
    for sw in ((), ('NOCHECKPOINT',), ('LOCKSTEP',), ('NODRAIN',), ('SHAREDCLOCK',), ('TERMREUSE',)):
        print('sw', sw, run(dg, int(sys.argv[1]) if len(sys.argv) > 1 else 250000, 11, sw), flush=True)
