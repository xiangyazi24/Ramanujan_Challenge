#!/usr/bin/env python3
import random
import q6612_check as c
import q6612_flat2v3 as m
rr=random.Random(7); dg=[rr.randrange(2) for _ in range(c.N)]; O=m.oracle_ring(dg)
refs=[m.reference(dg,O,s,top,sel) for k,s,p,B,top,sites in c.ARCS for sel in range(len(m.OFF[s]))]
print('CORE_DIGITS',''.join(map(str,dg)),'REFERENCE_ONES',sum(refs),'JOBS',len(refs),flush=True)
print('CORE_INTERACTION',c.interaction(dg),flush=True)
tops,pubs,fp=c.execute_wrong(dg); print('CORE_WRONG top_count',len(tops),'wrong',sum(z[3] for z in tops),'pub_count',len(pubs),'tops',tops,'pubs',pubs,'final_priv',sorted(set(fp.values())),flush=True)
print('CORE_CLEAN_MERGE',c.clean_merge(dg),flush=True)
print('CORE_NOCHECKPOINT',c.nocheckpoint_diff(dg,2000,300),flush=True)
# Literal permanent TERMREUSE block at single-top site 7.
X=[c.neutral() for _ in range(c.N)]; X[7]=c.neutral('ACC'); evs=[]
for _ in range(5): X[7],X[8],ev=m.pair(X[7],X[8],O[7],O[8],0,('TERMREUSE',)); evs.append(ev)
print('CORE_TERMREUSE_FOREVER',evs,'R7',X[7][4],flush=True)
