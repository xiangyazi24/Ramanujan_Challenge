#!/usr/bin/env python3
"""Control: same measurement on RANDOM rational maps of degree 6 and 9.
If control shows std ~ sqrt(19)~4.4 and ours shows 1.46, collapse is real."""
import sympy, random
from collections import Counter
random.seed(7)
# fixed random numerator/denominator coeffs over Z
num2=[random.randint(1,50) for _ in range(4)]   # deg3
den2=[random.randint(1,50) for _ in range(7)]   # deg6
num3=[random.randint(1,50) for _ in range(7)]   # deg6
den3=[random.randint(1,50) for _ in range(10)]  # deg9
def ev(cs,x,p):
    v=0
    for c in reversed(cs): v=(v*x+c)%p
    return v
devs=[]
import math
for p in [q for q in range(300,4000) if sympy.isprime(q)][:150]:
    c2=Counter(); c3=Counter()
    for r in range(p):
        d=ev(den2,r,p)
        if d: c2[ev(num2,r,p)*pow(d,p-2,p)%p]+=1
        d=ev(den3,r,p)
        if d: c3[ev(num3,r,p)*pow(d,p-2,p)%p]+=1
    cnt=sum(c2[a]*c3[a] for a in c2 if a in c3)
    devs.append((cnt-p)/math.isqrt(p))
import statistics
print(f"random deg(6,9) control over {len(devs)} primes: mean={statistics.mean(devs):.2f} std={statistics.stdev(devs):.2f} max={max(abs(d) for d in devs):.2f}")
print("(generic fiber product of deg 6,9 maps: genus ~ (6-1)(9-1)=40-ish arithmetic; expect std ~ 4-6)")
