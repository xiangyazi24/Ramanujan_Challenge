from sage.all import *

def P(n): return 34*n^3+51*n^2+27*n+5
p=17; m=16; Fp=GF(p)
b=[Fp(1),Fp(5)]
for n in range(2,m+1):
    b.append((Fp(P(n-1))*b[n-1]-Fp((n-1)^3)*b[n-2])/Fp(n^3))
R=PowerSeriesRing(Fp,'t',default_prec=m+1); t=R.gen()
F=R(b).add_bigoh(m+1)
Delta=(1-34*t+t^2).add_bigoh(m+1)
H=(F*Delta.sqrt()).add_bigoh(m+1)
A=(10*(5-26*t)*F^2).add_bigoh(m+1) # ++
B=(24*F^2).add_bigoh(m+1)            # --
C=(-40*(1+t)*F*H).add_bigoh(m+1)      # +-
D=(30*(t-1)*F*H).add_bigoh(m+1)       # -+
E1=((A+B+C+D)/4).add_bigoh(m+1)
E2=((A-B+C-D)/16).add_bigoh(m+1)
E3=((A-B-C+D)/36).add_bigoh(m+1)
E6=((A+B-C-D)/144).add_bigoh(m+1)
assert [E1[0],E2[0],E3[0],E6[0]] == [1,1,1,1]

def solve(src):
    z=[Fp(0)]*(m+1)
    for n in range(1,m+1):
        rhs=Fp(P(n-1))*z[n-1]+src[n]
        if n>=2: rhs-=Fp((n-1)^3)*z[n-2]
        z[n]=rhs/Fp(n^3)
    return z
ks=[]
for Ed in [E1,E2,E3,E6]:
    src=[Fp(0)]+[Fp(Ed[n]) for n in range(1,m+1)]
    ks.append(solve(src))
v=tuple(int(z[13]) for z in ks)
print('PF_VECTOR_17_13',v)
print('TR',sum(c*x for c,x in zip((-3,4,-9,108),v))%17)
print('SEP',sum(c*x for c,x in zip((1,-28,63,-36),v))%17)
