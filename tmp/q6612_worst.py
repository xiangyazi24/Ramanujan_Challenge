#!/usr/bin/env python3
import random
import q6612_check as c
import q6612_fast as f
import q6612_flat2v3 as m

def ewd_worst_line(L):
    N=L-1
    assert N%2==1
    # Q5302/EWD merged-root worst state: root=(1,1), follower 1=(1,0),
    # followers i>=2 alternate (0,1),(1,1).  Split root into H.x and F.x.
    out=[(1,1,0,1,1)]
    for i in range(1,N):
        if i==1: x,u=1,0
        elif i%2==0: x,u=0,1
        else: x,u=1,1
        out.append((x,u,i%3,1,1))
    out.append((1,0,0,1,1))
    assert len(out)==L
    return out

def all_arcs_worst():
    X=[list(c.neutral()) for _ in range(c.N)]
    for k,s,p,B,top,sites in c.ARCS:
        W=ewd_worst_line(len(sites))
        for q,site in enumerate(sites):
            x,u,_,v,bd=W[q]; X[site][k]=(x,u,q%len(m.OFF[s]),v,bd)
    return [tuple(z) for z in X]
rr=random.Random(7); dg=[rr.randrange(2) for _ in range(c.N)]; X=all_arcs_worst(); O=m.oracle_ring(dg)
print('WORST_INITIAL_PRIV',c.privilege_map(X,O),flush=True)
print('WORST_SERVICE',f.fast(dg,X,123456,2500000),flush=True)
