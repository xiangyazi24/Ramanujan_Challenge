#!/usr/bin/env python3
import random
import q6612_check as c
import q6612_flat2v3 as m
rr=random.Random(7); dg=[rr.randrange(2) for _ in range(c.N)]; O=m.oracle_ring(dg)
refs=[m.reference(dg,O,s,top,sel) for k,s,p,B,top,sites in c.ARCS for sel in range(len(m.OFF[s]))]
print('NOSMT_DIGITS',''.join(map(str,dg)),'REFERENCE_ONES',sum(refs),'JOBS',len(refs),flush=True)
print('NOSMT_INTERACTION',c.interaction(dg),flush=True)
print('NOSMT_WRONG',c.execute_wrong(dg)[:2],flush=True)
print('NOSMT_CLEAN_MERGE',c.clean_merge(dg),flush=True)
for seed in range(8): print('NOSMT_RANDOM',seed,c.simulate(dg,c.random_state(seed),30000+seed),flush=True)
print('NOSMT_ALTERNATING_WORST',c.simulate(dg,c.embed({'status':'fallback'},{'status':'fallback'}),818181,1200000),flush=True)
X,_=c.wrong_start(dg); print('NOSMT_WRONG_RANDOM',c.simulate(dg,X,919191),flush=True)
print('NOSMT_ALLTERM',c.simulate(dg,[c.neutral('REJ') for _ in range(c.N)],101010),flush=True)
print('NOSMT_NOCHECKPOINT',c.nocheckpoint_diff(dg,2000,300),flush=True)
