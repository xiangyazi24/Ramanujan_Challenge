#!/usr/bin/env python3
from __future__ import annotations
import collections, hashlib, json, random, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import q6612_flat2v3 as m
N=m.N; DOUBLE_TOPS={15,31,47,63}

def arcs():
    out=[]
    for k,(s,p) in enumerate(m.FAM):
        r=1<<s
        for B0 in range(p*r,N+p*r,2*r):
            out.append((k,s,p,B0%N,(B0+2*r-1)%N,tuple((B0+t)%N for t in range(2*r))))
    return out
ARCS=arcs(); ARC_KEYS=[(k,top) for k,s,p,B,top,sites in ARCS]
JOBS=[(k,top,sel) for k,s,p,B,top,sites in ARCS for sel in range(len(m.OFF[s]))]

def privilege_map(X,O):
    ans={}
    for k,(s,p) in enumerate(m.FAM):
        r=1<<s; vals=m.arc_privileges(X,O,k); Bs=list(range(p*r,len(X)+p*r,2*r))
        for B0,val in zip(Bs,vals): ans[(k,(B0+2*r-1)%len(X))]=val
    return ans

def solve_allpass_digits():
    from z3 import And,Bool,Solver,Xor,is_true,sat
    d=[Bool(f'd{i}') for i in range(N)]; z=Solver()
    for k,sc,p,B,top,sites in ARCS:
        r=1<<sc
        for x in m.OFF[sc]:
            left=(B+x)%N; right=(B+r+x)%N; earlier=[d[(B+y)%N] for y in m.OFF[sc] if y<x]
            carry=And(*earlier) if earlier else True
            z.add(Xor(d[left],carry)==d[right])
    assert z.check()==sat; M=z.model()
    return [1 if is_true(M.eval(x,model_completion=True)) else 0 for x in d]

def clock_privs(x,u):
    L=len(x); out=[]
    def up(i): return 1 if i==0 else 0 if i==L-1 else u[i]
    for i in range(L-1):
        d1=x[i]!=x[i+1]; d2=(not d1) and up(i)==1 and up(i+1)==0
        if d1 or d2: out.append((i,'D1' if d1 else 'D2'))
    return out

def clock_step(x,u,e):
    x=list(x); u=list(u); L=len(x)
    def up(i): return 1 if i==0 else 0 if i==L-1 else u[i]
    if x[e]!=x[e+1]:
        x[e+1]=x[e]
        if e+1<L-1: u[e+1]=1
    else:
        assert up(e)==1 and up(e+1)==0
        if e==0: x[0]^=1
        else: u[e]=0
    return tuple(x),tuple(u)

def smt_worst_witness(L,K,timeout_ms):
    from z3 import And,Bool,If,Implies,Int,Not,Or,Solver,Sum,sat
    x=[[Bool(f'x_{L}_{t}_{i}') for i in range(L)] for t in range(K+1)]
    u=[[None]+[Bool(f'u_{L}_{t}_{i}') for i in range(1,L-1)]+[None] for t in range(K+1)]
    edge=[Int(f'e_{L}_{t}') for t in range(K)]; z=Solver(); z.set(timeout=timeout_ms); z.add(Not(x[0][0]))
    def UP(t,i): return True if i==0 else False if i==L-1 else u[t][i]
    def P(t,i):
        d1=x[t][i]!=x[t][i+1]; d2=And(x[t][i]==x[t][i+1],UP(t,i),Not(UP(t,i+1)))
        return Or(d1,d2),d1,d2
    for t in range(K):
        ps=[P(t,i) for i in range(L-1)]; z.add(Sum([If(q[0],1,0) for q in ps])!=1); z.add(edge[t]>=0,edge[t]<L-1)
        for i,q in enumerate(ps): z.add(Implies(edge[t]==i,q[0]))
        z.add(x[t+1][0]==If(And(edge[t]==0,ps[0][2]),Not(x[t][0]),x[t][0]))
        for j in range(1,L): z.add(x[t+1][j]==If(And(edge[t]==j-1,ps[j-1][1]),x[t][j-1],x[t][j]))
        for j in range(1,L-1): z.add(u[t+1][j]==If(And(edge[t]==j-1,ps[j-1][1]),True,If(And(edge[t]==j,ps[j][2]),False,u[t][j])))
    z.add(Sum([If(P(K,i)[0],1,0) for i in range(L-1)])==1)
    st=z.check()
    if st!=sat: return {'status':str(st)}
    M=z.model(); xx=tuple(1 if bool(M.eval(x[0][i],model_completion=True)) else 0 for i in range(L)); uu=[0]*L
    for i in range(1,L-1): uu[i]=1 if bool(M.eval(u[0][i],model_completion=True)) else 0
    moves=[M.eval(edge[t]).as_long() for t in range(K)]; a=xx; b=tuple(uu); first=None
    for t,e in enumerate(moves):
        assert len(clock_privs(a,b))!=1; a,b=clock_step(a,b,e)
        if len(clock_privs(a,b))==1 and first is None: first=t+1
    assert first==K
    return {'status':'sat','x':xx,'u':tuple(uu),'moves':moves,'first_legit':first,'final_priv':clock_privs(a,b)}

def neutral(R=None): return tuple((0,1,0,'N',0) for _ in m.FAM)+(R,)

def interaction(digits):
    O=m.oracle_ring(digits); X=[neutral() for _ in range(N)]; q=[list(z) for z in X]
    q[6][1]=(0,1,0,'N',0); q[7][1]=(1,0,0,'N',0); q[7][0]=(0,1,0,'N',0); q[8][0]=(1,1,0,'N',0); X=[tuple(z) for z in q]; tr=[]
    for e in (6,7,7):
        j=(e+1)%N; before=(X[e],X[j]); X[e],X[j],ev=m.pair(X[e],X[j],O[e],O[j],0,()); tr.append((e,ev,before,(X[e],X[j])))
    return tr

def wrong_start(digits):
    O=m.oracle_ring(digits); X=[list(neutral()) for _ in range(N)]; exp={}
    for k,s,p,B,top,sites in ARCS:
        ref=m.reference(digits,O,s,top,0); pred=(top-1)%N; X[top][k]=(1,1,0,'N',0); X[pred][k]=(0,1,0,'N',1-ref); exp[(k,top)]=ref
    return [tuple(z) for z in X],exp

def execute_wrong(digits):
    O=m.oracle_ring(digits); X,exp=wrong_start(digits); assert set(privilege_map(X,O).values())=={1}; tops=[]; pubs=[]
    for e in sorted({(top-1)%N for k,s,p,B,top,sites in ARCS}):
        j=(e+1)%N; X[e],X[j],ev=m.pair(X[e],X[j],O[e],O[j],0,())
        for z in ev:
            if len(z)>1 and z[1]=='TOP': tops.append((e,z,exp[(z[0],j)],z[3]!=exp[(z[0],j)]))
    for top in sorted({top for k,s,p,B,top,sites in ARCS}):
        j=(top+1)%N; X[top],X[j],ev=m.pair(X[top],X[j],O[top],O[j],0,()); pubs += [(top,z) for z in ev if z[0]=='PUBLISH']
    return tops,pubs,privilege_map(X,O)

def random_state(seed):
    r=random.Random(seed); return [tuple(r.choice(m.CS[s]) for s,_ in m.FAM)+(r.choice(m.RS),) for _ in range(N)]

def simulate(digits,X0,seed,cap=1200000,sw=()):
    r=random.Random(seed); O=m.oracle_ring(digits); X=list(X0); p=privilege_map(X,O); au={a:(0 if p[a]==1 else None) for a in ARC_KEYS}; gu=0 if all(v==1 for v in p.values()) else None
    started={a:False for a in ARC_KEYS}; jc=collections.Counter(); early=collections.Counter(); earlybad=collections.Counter(); ghost=[{'ca':False,'cr':False,'qa':False,'qr':False} for _ in range(N)]
    tops=correct=postbad=pubs=pwrong=punqual=pafter=pwrongafter=punqualafter=0; merges=collections.Counter(); dmerges=collections.Counter(); tw=None
    for t in range(cap):
        pb=privilege_map(X,O)
        if gu is None and all(v==1 for v in pb.values()): gu=t
        for a,v in pb.items():
            if au[a] is None and v==1: au[a]=t
        i=r.randrange(N); j=(i+1)%N; oldR=X[j][4]; X[i],X[j],ev=m.pair(X[i],X[j],O[i],O[j],r.randrange(2),sw)
        if ev and ev[0][0]=='PUBLISH':
            val=ev[0][1]; pubs+=1; sound=ghost[i]['cr' if val=='REJ' else 'ca']; qual=ghost[i]['qr' if val=='REJ' else 'qa']; pwrong+=not sound; punqual+=not qual
            if gu is not None and t>=gu: pafter+=1; pwrongafter+=not sound; punqualafter+=not qual
            if 'TERMREUSE' not in sw: ghost[i]={'ca':False,'cr':False,'qa':False,'qr':False}
        else:
            rv=oldR
            for z in ev:
                if len(z)<2: continue
                arc=(z[0],j)
                if z[1]=='BOTTOM':
                    if pb.get(arc)==1: started[arc]=True
                elif z[1]=='TOP':
                    tops+=1; s=m.FAM[z[0]][0]; ref=m.reference(digits,O,s,j,z[2]); ok=z[3]==ref; correct+=ok
                    if au[arc] is None or t<au[arc]: early[arc]+=1; earlybad[arc]+=not ok
                    if gu is not None and t>=gu and not ok: postbad+=1
                    inc='REJ' if z[3] else 'ACC'; nr='REJ' if inc=='REJ' or rv=='REJ' else 'ACC'
                    if rv is not None: merges[(rv,inc,nr)]+=1; dmerges[(rv,inc,nr)]+=int(j in DOUBLE_TOPS)
                    if ok: ghost[j]['cr' if inc=='REJ' else 'ca']=True
                    qual=ok and pb.get(arc)==1 and started[arc]
                    if qual: ghost[j]['qr' if inc=='REJ' else 'qa']=True
                    if gu is not None and t>=gu and qual: jc[(z[0],j,z[2])]+=1
                    started[arc]=False; rv=nr
            assert X[j][4]==rv
        if gu is None:
            pa=privilege_map(X,O)
            if all(v==1 for v in pa.values()): gu=t+1; au={a:(t+1 if au[a] is None else au[a]) for a in ARC_KEYS}
        if gu is not None and tw is None and all(jc[j]>=2 for j in JOBS): tw=t+1
        if tw is not None and not sw: break
    fp=privilege_map(X,O)
    return {'stop':t+1,'global_unique':gu,'arc_unique_max':max(v for v in au.values() if v is not None),'ununique':sum(v is None for v in au.values()),'twice':tw,'jobs_twice':sum(jc[j]>=2 for j in JOBS),'tops':tops,'correct':correct,'post_unique_wrong_top':postbad,'early_top':sum(early.values()),'early_wrong':sum(earlybad.values()),'max_early':max(early.values(),default=0),'max_early_wrong':max(earlybad.values(),default=0),'pubs':pubs,'pub_after_unique':pafter,'pub_wrong':pwrong,'pub_unqualified':punqual,'pub_wrong_after_unique':pwrongafter,'pub_unqualified_after_unique':punqualafter,'merges':dict(merges),'double_merges':dict(dmerges),'final_priv':(min(fp.values()),max(fp.values()))}

def embed(w16,w32):
    X=[list(neutral()) for _ in range(N)]
    for k,s,p,B,top,sites in ARCS:
        w=w16 if s==3 else w32
        xs=w['x'] if w.get('status')=='sat' else tuple(t&1 for t in range(len(sites))); us=w['u'] if w.get('status')=='sat' else tuple(0 if t%3==0 else 1 for t in range(len(sites)))
        for t,site in enumerate(sites): X[site][k]=(xs[t],us[t],t%len(m.OFF[s]),1,1)
    return [tuple(z) for z in X]

def clean_merge(digits):
    O=m.oracle_ring(digits); X=[list(neutral()) for _ in range(N)]; targets=[]
    for top in sorted(DOUBLE_TOPS):
        for k,s,p,B,tt,sites in [a for a in ARCS if a[4]==top]:
            ref=m.reference(digits,O,s,top,0); pred=(top-1)%N; X[top][k]=(1,1,0,'N',0); X[pred][k]=(0,1,0,'N',ref); targets.append((k,top,ref))
    X=[tuple(z) for z in X]; te=[]; pu=[]
    for e in sorted({top-1 for top in DOUBLE_TOPS}):
        j=e+1; X[e],X[j],ev=m.pair(X[e],X[j],O[e],O[j],0,()); te += [(j,z) for z in ev if len(z)>1 and z[1]=='TOP']
    for top in sorted(DOUBLE_TOPS):
        j=(top+1)%N; X[top],X[j],ev=m.pair(X[top],X[j],O[top],O[j],0,()); pu += [(top,z) for z in ev if z[0]=='PUBLISH']
    return targets,te,pu

def nocheckpoint_diff(digits,trials=2000,steps=300):
    O=m.oracle_ring(digits); r=random.Random(6612); first=None; topdiff=None
    for tr in range(trials):
        A=[r.choice(m.CS[3]) for _ in range(16)]; B=list(A); word=[]
        for st in range(steps):
            e=r.randrange(15); word.append(e); ga='H' if e==0 else 'I'; gb='F' if e+1==15 else 'I'
            a1,b1,v1=m.clock(3,A[e],A[e+1],ga,gb,O[e+1],0,()); a2,b2,v2=m.clock(3,B[e],B[e+1],ga,gb,O[e+1],0,('NOCHECKPOINT',)); A[e],A[e+1]=a1,b1; B[e],B[e+1]=a2,b2
            if first is None and A!=B: first=(tr,st+1,word.copy(),tuple(A),tuple(B))
            s1=v1 if v1 and v1[0]=='TOP' else None; s2=v2 if v2 and v2[0]=='TOP' else None
            if s1!=s2: topdiff=(tr,st+1,word.copy(),s1,s2); return {'first_state':first,'topdiff':topdiff,'trials':tr+1}
    return {'first_state':first,'topdiff':topdiff,'trials':trials}

def main():
    print('SOURCE_SHA256',hashlib.sha256(Path(__file__).with_name('q6612_flat2v3.py').read_bytes()).hexdigest(),flush=True)
    print('TOTALITY',m.totality(),flush=True); print('ARCS',[(k,s,p,B,t) for k,s,p,B,t,_ in ARCS],flush=True)
    print('CLOCK_EXACT 16',211,'BOUND',256,'CLOCK_EXACT 32',931,'BOUND',1024,flush=True)
    t=time.time(); w16=smt_worst_witness(16,211,240000); print('SMT16',{k:(len(v) if k=='moves' else v) for k,v in w16.items()},'SEC',round(time.time()-t,3),flush=True)
    t=time.time(); w32=smt_worst_witness(32,931,420000); print('SMT32',{k:(len(v) if k=='moves' else v) for k,v in w32.items()},'SEC',round(time.time()-t,3),flush=True)
    ap=solve_allpass_digits(); O=m.oracle_ring(ap); refs=[m.reference(ap,O,s,top,sel) for k,s,p,B,top,sites in ARCS for sel in range(len(m.OFF[s]))]; print('ALLPASS',''.join(map(str,ap)),'ONES',sum(refs),'JOBS',len(refs),flush=True)
    print('INTERACTION',interaction(ap),flush=True)
    tops,pubs,fp=execute_wrong(ap); print('WRONG_FIRST top_count',len(tops),'wrong',sum(z[3] for z in tops),'pub_count',len(pubs),'tops',tops,'pubs',pubs,'final_priv',sorted(set(fp.values())),flush=True)
    print('CLEAN_MERGE',clean_merge(ap),flush=True)
    rr=random.Random(7); dg=[rr.randrange(2) for _ in range(N)]
    for sw in ((),('NOCHECKPOINT',),('LOCKSTEP',),('NODRAIN',),('SHAREDCLOCK',),('TERMREUSE',)): print('AUTHOR',sw,m.run(dg,250000,11,sw),flush=True)
    for seed in range(8): print('RANDOM',seed,simulate(dg,random_state(seed),10000+seed),flush=True)
    print('WORST',simulate(dg,embed(w16,w32),424242,2000000),flush=True)
    X,_=wrong_start(ap); print('WRONG_RANDOM',simulate(ap,X,31337),flush=True)
    print('ALL_TERMINALS',simulate(ap,[neutral('REJ') for _ in range(N)],271828),flush=True)
    O=m.oracle_ring(ap); X=[neutral() for _ in range(N)]; X[7]=neutral('ACC'); evs=[]
    for _ in range(5): X[7],X[8],ev=m.pair(X[7],X[8],O[7],O[8],0,('TERMREUSE',)); evs.append(ev)
    print('TERMREUSE_FOREVER',evs,'R7',X[7][4],flush=True)
    print('NOCHECKPOINT_DIFF',nocheckpoint_diff(dg),flush=True)
if __name__=='__main__': main()
