#!/usr/bin/env python3
import collections, random
import q6612_check as c
import q6612_flat2v3 as m

def event_arc(k, kind, edge, j):
    s=m.FAM[k][0]; r=1<<s
    if kind=='TOP': return (k,j)
    if kind=='BOTTOM': return (k,(edge+2*r-1)%c.N)
    return None

def fast(digits,X0,seed,cap=1500000,sw=()):
    rnd=random.Random(seed); O=m.oracle_ring(digits); X=list(X0); p=c.privilege_map(X,O)
    au={a:(0 if p[a]==1 else None) for a in c.ARC_KEYS}; gu=0 if all(v==1 for v in p.values()) else None
    started={a:False for a in c.ARC_KEYS}; jc=collections.Counter(); early=collections.Counter(); earlybad=collections.Counter()
    ghost=[{'ca':False,'cr':False,'qa':False,'qr':False} for _ in range(c.N)]
    tops=correct=postbad=pubs=pwrong=punqual=pafter=pwrongafter=punqualafter=0; merges=collections.Counter(); dmerges=collections.Counter(); tw=None
    ones={a:1 for a in c.ARC_KEYS}
    for t in range(cap):
        if gu is None:
            pb=c.privilege_map(X,O)
            if all(v==1 for v in pb.values()): gu=t; pb=ones
            for a,v in pb.items():
                if au[a] is None and v==1: au[a]=t
        else: pb=ones
        i=rnd.randrange(c.N); j=(i+1)%c.N; oldR=X[j][4]; X[i],X[j],ev=m.pair(X[i],X[j],O[i],O[j],rnd.randrange(2),sw)
        if ev and ev[0][0]=='PUBLISH':
            val=ev[0][1]; pubs+=1; sound=ghost[i]['cr' if val=='REJ' else 'ca']; qual=ghost[i]['qr' if val=='REJ' else 'qa']; pwrong+=not sound; punqual+=not qual
            if gu is not None and t>=gu: pafter+=1; pwrongafter+=not sound; punqualafter+=not qual
            if 'TERMREUSE' not in sw: ghost[i]={'ca':False,'cr':False,'qa':False,'qr':False}
        else:
            rv=oldR
            for z in ev:
                if len(z)<2: continue
                arc=event_arc(z[0],z[1],i,j)
                if z[1]=='BOTTOM':
                    if pb.get(arc)==1: started[arc]=True
                elif z[1]=='TOP':
                    tops+=1; s=m.FAM[z[0]][0]; ref=m.reference(digits,O,s,j,z[2]); ok=z[3]==ref; correct+=ok
                    if au[arc] is None or t<au[arc]: early[arc]+=1; earlybad[arc]+=not ok
                    if gu is not None and t>=gu and not ok: postbad+=1
                    inc='REJ' if z[3] else 'ACC'; nr='REJ' if inc=='REJ' or rv=='REJ' else 'ACC'
                    if rv is not None: merges[(rv,inc,nr)]+=1; dmerges[(rv,inc,nr)]+=int(j in c.DOUBLE_TOPS)
                    if ok: ghost[j]['cr' if inc=='REJ' else 'ca']=True
                    qual=ok and pb.get(arc)==1 and started[arc]
                    if qual: ghost[j]['qr' if inc=='REJ' else 'qa']=True
                    if gu is not None and t>=gu and qual: jc[(z[0],j,z[2])]+=1
                    started[arc]=False; rv=nr
            assert X[j][4]==rv
        if gu is None:
            pa=c.privilege_map(X,O)
            if all(v==1 for v in pa.values()):
                gu=t+1
                for a in c.ARC_KEYS:
                    if au[a] is None: au[a]=t+1
        if gu is not None and tw is None and all(jc[j]>=2 for j in c.JOBS):
            tw=t+1; break
    fp=c.privilege_map(X,O)
    return {'stop':t+1,'global_unique':gu,'arc_unique_max':max(v for v in au.values() if v is not None),'ununique':sum(v is None for v in au.values()),'twice':tw,'jobs_twice':sum(jc[j]>=2 for j in c.JOBS),'tops':tops,'correct':correct,'post_unique_wrong_top':postbad,'early_top':sum(early.values()),'early_wrong':sum(earlybad.values()),'max_early':max(early.values(),default=0),'max_early_wrong':max(earlybad.values(),default=0),'pubs':pubs,'pub_after_unique':pafter,'pub_wrong':pwrong,'pub_unqualified':punqual,'pub_wrong_after_unique':pwrongafter,'pub_unqualified_after_unique':punqualafter,'merges':dict(merges),'double_merges':dict(dmerges),'final_priv':(min(fp.values()),max(fp.values()))}

def main():
    rr=random.Random(7); dg=[rr.randrange(2) for _ in range(c.N)]
    print('FAST_DIGITS',''.join(map(str,dg)),flush=True)
    for seed in range(8): print('FAST_RANDOM',seed,fast(dg,c.random_state(seed),40000+seed),flush=True)
    print('FAST_ALTERNATING',fast(dg,c.embed({'status':'fallback'},{'status':'fallback'}),515151,2000000),flush=True)
    X,_=c.wrong_start(dg); print('FAST_WRONG_START',fast(dg,X,616161),flush=True)
    print('FAST_ALL_TERMINALS',fast(dg,[c.neutral('REJ') for _ in range(c.N)],717171),flush=True)
if __name__=='__main__': main()
