from sage.all import *
import os, time

P_LO = ZZ(os.environ['P_LO'])
P_HI = ZZ(os.environ['P_HI'])
P_GLOBAL = 5000
DS = [1,2,3,6]
OUT = os.environ.get('OUTFILE', 'q7691_samples_%s_%s.tsv' % (P_LO,P_HI))


def P(n):
    n=ZZ(n)
    return 34*n^3+51*n^2+27*n+5

sig3=[ZZ(0)]*(P_GLOBAL+1)
for d in range(1,P_GLOBAL+1):
    dd=ZZ(d)^3
    for m in range(d,P_GLOBAL+1,d):
        sig3[m]+=dd


def apery_mod_targets(p):
    Fp=GF(p)
    b=[Fp(1),Fp(5)]
    for n in range(1,p-1):
        nxt=(Fp(P(n))*b[n]-Fp(n)^3*b[n-1])/Fp(n+1)^3
        b.append(nxt)
    assert b[0]==1 and b[p-1]==1
    return b,[r for r in range(1,p-1) if b[r]==0]


def q_objects_mod(p,prec):
    Fp=GF(p)
    PS=PowerSeriesRing(Fp,'q',default_prec=prec)
    q=PS.gen()
    A=PS(1)
    for a in range(1,prec):
        if a%3 != 0:
            A *= (1+q^a)^12
            A=A.add_bigoh(prec)
    E=PS(1)
    for m in range(1,prec):
        e=-5+7*(m%2==0)+7*(m%3==0)-5*(m%6==0)
        if e>0:
            E *= (1-q^m)^e
        elif e<0:
            E /= (1-q^m)^(-e)
        E=E.add_bigoh(prec)
    H=(1-q*A.derivative()/A).add_bigoh(prec)
    EH=(E*H).add_bigoh(prec)
    return A,EH


def target_state(p,r,bprev,A,EH):
    Fp=GF(p)
    base=(EH*(A^r)).add_bigoh(r+1)
    W={}
    for d in DS:
        kap=Fp(0)
        scale=Fp(240)/Fp(d)^3
        for m in range(1,r//d+1):
            kap += scale*Fp(sig3[m])/Fp(m)^3*base[r-d*m]
        W[d]=Fp(r)^3*bprev*kap
    assert W[1]-28*W[2]+63*W[3]-36*W[6]==Fp(240)
    return [ZZ(W[d]) for d in DS]


t0=time.time()
rows=[]
prime_counts=[]
for p in prime_range(max(ZZ(7),P_LO),P_HI+1):
    p=ZZ(p)
    b,targets=apery_mod_targets(p)
    if not targets:
        continue
    A,EH=q_objects_mod(p,max(targets)+1)
    prime_counts.append((p,len(targets)))
    for r in targets:
        W=target_state(p,ZZ(r),b[r-1],A,EH)
        rows.append([p,ZZ(r),ZZ(b[r-1])]+W)

with open(OUT,'w') as fh:
    fh.write('p\tr\tbprev\tW1\tW2\tW3\tW6\n')
    for row in rows:
        fh.write('\t'.join(map(str,row))+'\n')

print('CHUNK',P_LO,P_HI,'TARGET_PAIRS',len(rows),'TARGET_PRIMES',len(prime_counts),
      'MAX_TARGETS',max([c for _,c in prime_counts] or [0]),
      'SECONDS',time.time()-t0,'OUT',OUT)
