#!/usr/bin/env python3
import random, time
import q6612_check as c
import q6612_flat2v3 as m
print('QUICK_CLOCK_THEOREM',211,931,flush=True)
t=time.time(); w16=c.smt_worst_witness(16,211,180000); print('QUICK_SMT16',{k:(len(v) if k=='moves' else v) for k,v in w16.items()},'SEC',round(time.time()-t,3),flush=True)
w32={'status':'fallback'}
ap=c.solve_allpass_digits(); O=m.oracle_ring(ap); refs=[m.reference(ap,O,s,top,sel) for k,s,p,B,top,sites in c.ARCS for sel in range(len(m.OFF[s]))]
print('QUICK_ALLPASS',''.join(map(str,ap)),'ONES',sum(refs),flush=True)
print('QUICK_INTERACTION',c.interaction(ap),flush=True)
print('QUICK_WRONG',c.execute_wrong(ap)[:2],flush=True)
rr=random.Random(7); dg=[rr.randrange(2) for _ in range(c.N)]
for seed in range(4): print('QUICK_RANDOM',seed,c.simulate(dg,c.random_state(seed),20000+seed),flush=True)
print('QUICK_FALLBACK_WORST',c.simulate(dg,c.embed(w16,w32),515151,1200000),flush=True)
X,_=c.wrong_start(ap); print('QUICK_WRONG_RANDOM',c.simulate(ap,X,616161),flush=True)
print('QUICK_ALLTERM',c.simulate(ap,[c.neutral('REJ') for _ in range(c.N)],717171),flush=True)
print('QUICK_NOCHECKPOINT',c.nocheckpoint_diff(dg,500,200),flush=True)
