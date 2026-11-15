import itertools, time, sys
sys.path.insert(0, 'tmp')
import q6596_tm5 as tm5

tm5.CODES, tm5.INV = tm5.tail_codes(600)

def fixed_data(name, c, mut=None):
    PH = tm5.phases(c)
    bad, R = tm5.reference_bad(c)
    flush, S = tm5.cartesian(c, mut=mut, T=400)
    unique = all(len(s) == 1 for s in S)
    reject = sorted(i for i,s in enumerate(S) if any(x[5] for x in s))
    fp = tuple(next(iter(s)) for s in S) if unique else None
    print(f'CASE {name} mut={mut} n={len(c)}')
    print('  phases', PH)
    print('  reference_bad', bad)
    print('  residues', R)
    print('  flush', flush, 'unique', unique, 'reject', reject)
    if fp is not None:
        print('  fixed_point', fp)
        for i in reject:
            print('  reject_detail site', i, 'phase', PH[i], 'left_phase', PH[i-1], 'label', c[i], 'state', fp[i])
    return PH, bad, R, flush, unique, reject, fp

A16 = [int(ch) for ch in '4554544154434512']
splice = [tm5.alpha(k) for k in range(100,164)] + [tm5.alpha(k) for k in range(300,364)]
PH16,bad16,R16,f16,u16,rej16,fp16 = fixed_data('alpha16', A16)
PHsp,badsp,Rsp,fsp,usp,rejsp,fpsp = fixed_data('splice128', splice)
fixed_data('splice128', splice, 'M_NOPOS')
fixed_data('alpha16', A16, 'M_NOLAUNCH')

print('OFFSET_ANALYSIS')
for s in bad16:
    print('  bad_start', s, 'R[s]', R16[s], 'R[s+1]', R16[(s+1)%len(A16)],
          'expected', None if R16[s] is None else (R16[s]+1)%32)
for i in rej16:
    print('  reject_site', i, 'candidate_bit3_bad_start', (i-(8+3))%len(A16),
          'candidate_bit4_bad_start', (i-(8+6))%len(A16))

def verdict(c):
    n = len(c)
    PH = tm5.phases(c)
    defect_sites = [i for i in range(n)
                    if PH[i] is None or PH[i-1] is None or PH[i] != (PH[i-1]+1)%8]
    X = tuple((0,)*6 for _ in c)
    for _ in range(3*n+40):
        X = tm5.step(X, c, PH)
    reject_sites = [i for i,x in enumerate(X) if x[5]]
    return defect_sites, reject_sites, X

print('EQUIVALENCE_EXHAUSTIVE')
grand_tot = grand_dis = 0
first = None
t0 = time.time()
for n in range(3,7):
    tot = agree = disagree = 0
    category = {}
    for ctuple in itertools.product(range(6), repeat=n):
        c = list(ctuple)
        bad, R = tm5.reference_bad(c)
        ds, rs, X = verdict(c)
        lhs = bool(bad); rhs = bool(ds or rs)
        key = (lhs, bool(ds), bool(rs))
        category[key] = category.get(key,0)+1
        tot += 1
        if lhs == rhs:
            agree += 1
        else:
            disagree += 1
            if first is None:
                first = (n, c, bad, R, ds, rs, tm5.phases(c), X)
    grand_tot += tot; grand_dis += disagree
    print('  n', n, 'words', tot, 'agreements', agree, 'disagreements', disagree,
          'categories', sorted(category.items()))
print('  total_words', grand_tot, 'total_disagreements', grand_dis,
      'elapsed_seconds', round(time.time()-t0,3))
if first is None:
    print('  FIRST_DISAGREEMENT NONE')
else:
    n,c,bad,R,ds,rs,PH,X = first
    print('  FIRST_DISAGREEMENT n',n,'word',c,'reference_bad',bad,'residues',R,
          'phase_defect_sites',ds,'reject_sites',rs,'phases',PH,'fixed',X)
