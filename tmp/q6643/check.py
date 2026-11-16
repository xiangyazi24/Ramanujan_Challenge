import itertools, json, random
import serial16 as S
import etbstar as E
import serial32 as P

def unique_clock_dp():
    d = {}
    for x0 in (0, 1): d[(x0, 1, 0, None)] = 2
    for i in range(1, 16):
        nd = {}; us = (0, 1) if i < 15 else (0,); end_mult = 1 if i < 15 else 2
        for (px, pu, cnt, lab), mult in d.items():
            for x in (0, 1):
                for u in us:
                    ncnt, nlab = cnt, lab
                    if px != x: ncnt += 1; nlab = ('UP', i - 1)
                    elif pu == 1 and u == 0: ncnt += 1; nlab = ('DOWN', i)
                    if ncnt <= 1: nd[(x, u, ncnt, nlab)] = nd.get((x, u, ncnt, nlab), 0) + mult * end_mult
        d = nd
    by = {}
    for (_, _, cnt, lab), mult in d.items():
        if cnt == 1: by[lab] = by.get(lab, 0) + mult
    return by

def clock_variants(kind, p):
    out = []
    for left_x, u0, u15 in itertools.product((0, 1), repeat=3):
        z = []
        for t in range(16):
            if kind == 'UP': x, semu = (left_x if t <= p else 1-left_x), (1 if t <= p else 0)
            else: x, semu = left_x, (1 if t < p else 0)
            z.append((x, u0 if t == 0 else u15 if t == 15 else semu))
        out.append(z)
    return out

def exhaustive_unique_front_scope():
    by = unique_clock_dp(); labels = {('UP',p) for p in range(15)} | {('DOWN',p) for p in range(1,16)}
    assert set(by) == labels and set(by.values()) == {8} and sum(by.values()) == 240
    checks = stutters = closes = self_fail = wrong_del = wrong_ack = dels = acks = 0
    for st in E.CANON:
        kind, p, _ = st[0]
        for cv in clock_variants(kind, p):
            X = E.build(st); X = [(cv[t][0], cv[t][1]) + X[t][2:] for t in range(16)]
            if E.project(X) != st: self_fail += 1; continue
            e = p if kind == 'UP' else p-1
            for i in range(15):
                if i != e:
                    a,b,ev = S.pair(i,X[i],X[i+1]); stutters += 1; assert (a,b,ev)==(X[i],X[i+1],[])
            opn = sorted(st[1][2]); a,b,ev=S.pair(e,X[e],X[e+1]); X[e],X[e+1]=a,b; checks+=1
            closes += E.ETBstar(X)
            for q in ev:
                if q[0]=='DELIVER': dels+=1; wrong_del += not(len(opn)==1 and S.parent(opn[0])==q[1])
                elif q[0]=='ACK': acks+=1; wrong_ack += [t for t in range(16) if X[t][4][3]] != [S.parent(q[1])]
    return dict(unique_front_raw_clocks=sum(by.values()),front_labels=len(by),raw_clocks_per_label_minmax=[min(by.values()),max(by.values())],canonical_tuples=len(E.CANON),tuple_clock_checks=checks,nonfront_stutter_checks=stutters,projection_failures=self_fail,closure_successes=closes,closure_failures=checks-closes,deliveries=dels,wrong_deliveries=wrong_del,acks=acks,wrong_acks=wrong_ack)

def exact_cycles():
    rows=[]
    for u0,u15 in itertools.product((0,1),repeat=2):
        X=S.clean(); X[0]=(0,u0)+X[0][2:]; X[15]=(0,u15)+X[15][2:]
        seen={}; events=[]; projections=[]; bad=wd=wa=0; t=0
        while tuple(X) not in seen:
            seen[tuple(X)]=t; bad += not E.ETBstar(X); projections.append(E.project(X)); F=S.fronts(X); assert len(F)==1
            kind,p,_=F[0]; edge=p if kind=='UP' else p-1; pre=[s for s in range(16) if X[s][4][2]]; ev=S.step(X,edge,())
            for q in ev:
                if q[0]=='DELIVER': wd += not(len(pre)==1 and S.parent(pre[0])==q[1])
                elif q[0]=='ACK': wa += [s for s in range(16) if X[s][4][3]] != [S.parent(q[1])]
            events.append(ev); t+=1; assert t<5000
        mu=seen[tuple(X)]; lam=t-mu; cev=[q for ev in events[mu:] for q in ev]
        rows.append(dict(u0=u0,u15=u15,transient=mu,cycle=lam,projected_cycle_states=len(set(projections[mu:])),etbstar_violations_before_repeat=bad,deliveries_in_cycle=sum(q[0]=='DELIVER' for q in cev),acks_in_cycle=sum(q[0]=='ACK' for q in cev),wrong_deliveries=wd,wrong_acks=wa))
    return rows

def product_witness(limit=100000,seed=0):
    rnd=random.Random(seed); X=P.clean(); first_two=first_same=None; shifted={}
    for t in range(limit):
        i=rnd.randrange(32); ev=P.pair32(X,i,()); V=P.RETB(X)
        act=[(p,r) for p,r in P.REG if not V[P.NAME[(p,r)]] and any(s[4][2] for s in P.region_slice(X,p,r))]
        if first_two is None and len(act)>=2: first_two={'service':t+1,'edge':i,'active':[P.NAME[a] for a in act]}
        if first_same is None:
            for p in (0,1):
                same=[a for a in act if a[0]==p]
                if len(same)>=2: first_same={'service':t+1,'edge':i,'active':[P.NAME[a] for a in same]}; break
        for p,r,kind,*rest in ev:
            if p==1 and kind in ('DELIVER','ACK') and rest:
                shifted.setdefault((P.NAME[(p,r)],kind,rest[0]),{'service':t+1,'global':P.glob(p,r,rest[0])})
        if first_two and first_same and all((name,kind,local) in shifted for name in ('G10','G11') for kind,local in (('DELIVER',15),('ACK',7))): break
    return {'first_two_active':first_two,'first_same_cover_active':first_same,'shifted_events':{str(k):v for k,v in sorted(shifted.items())},'services_examined':t+1}

def structural():
    r5=False
    for p,r in P.REG:
        for t in range(16): r5 |= P.glob(p,r,t)!=P.glob(p,r,t) or P.loc(1-p,P.glob(p,r,t))[1]<0
    return {'frame_boundaries':{p:[i for i in range(32) if P.loc(p,i)[1]==15] for p in (0,1)},'demands':{P.NAME[(p,r)]:{'7':[P.glob(p,r,7),P.glob(p,r,15)],'9':[P.glob(p,r,9),P.glob(p,r,11)]} for p,r in P.REG},'r5_runtime_expression_can_fail':r5,'r8_runtime_field_or_test_present':False,'per_site_alphabet':S.SITE_COUNT**2}

print('ETBSTAR_LITERAL',json.dumps(E.check(),sort_keys=True))
print('UNIQUE_FRONT_SCOPE',json.dumps(exhaustive_unique_front_scope(),sort_keys=True))
print('BOTTOM_CYCLES',json.dumps(exact_cycles(),sort_keys=True))
print('DRAIN_CUTS_IN_ETBSTAR',sum(1 for st in E.CANON if st[0][0]=='DOWN' and st[0][2]==S.FREE and st[1]==(frozenset(),)*4))
print('INSTANCES',json.dumps(E.instances(),sort_keys=True))
print('PRODUCT_WITNESS',json.dumps(product_witness(),sort_keys=True))
print('STRUCTURAL',json.dumps(structural(),sort_keys=True))
print('TOTALITY',json.dumps(P.totality(),sort_keys=True))
print('BASE_CLEAN',json.dumps(P.run(150000,0,(),P.clean()),sort_keys=True))
print('BASE_GARBAGE1',json.dumps(P.run(150000,1),sort_keys=True))
print('BASE_GARBAGE2',json.dumps(P.run(150000,2),sort_keys=True))
for sw in (('SHARED-LANE',),('NO-SEAM-DRAIN',),('NO-MIDPOINT-OWNER',),('REMOTE-BUSY',),('HEADLESS-REGION',),('REMOTE-BUSY','HEADLESS-REGION')):
    X=P.clean()
    if sw==('REMOTE-BUSY','HEADLESS-REGION'):
        g=P.glob(0,1,7);z=X[g][0];X[g]=(z[:4]+((0,1,1,0,0),),X[g][1])
        for l in range(16):
            g=P.glob(0,1,l);z=X[g][0];X[g]=((0,1)+z[2:],X[g][1])
    print('CONTROL',sw,json.dumps(P.run(100000,4,sw,X),sort_keys=True))
