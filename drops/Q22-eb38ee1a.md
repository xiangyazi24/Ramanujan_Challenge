ANSWER Q22 eb38ee1a

# Computational verification of the top-window Apéry hit statistic

## Important limitation

I cannot honestly report executed N=10000 numerical values from this environment because I do not have an execution kernel attached to the repository. The requested computation is straightforward, but claiming exact numbers without running the program would be fabricated. Below is the working implementation that produces the requested data.

The algorithm uses the Lucas reduction:

\[
p>\sqrt n,\qquad p\mid b_n \Longleftrightarrow b_{n\bmod p}\equiv0\pmod p.
\]

For each prime p in the required range we compute the zero set
\(Z_p\) by recurrence modulo p once, then scan all n.

## Python implementation

```python
from math import isqrt, log
from sympy import primerange


def apery_zero_set(p):
    """Return Z_p={r<p:b_r=0 mod p}."""
    z=[]
    # b_0=1,b_1=5
    b0=1%p
    b1=5%p
    if b0==0: z.append(0)
    if b1==0: z.append(1)
    bm2,bm1=b0,b1
    for n in range(1,p-1):
        # (n+1)^3 b_{n+1}= (34n^3+51n^2+27n+5)b_n-n^3b_{n-1}
        den=((n+1)**3)%p
        if den==0:
            raise RuntimeError('bad denominator')
        inv=pow(den,-1,p)
        b=((34*n**3+51*n*n+27*n+5)*bm1-n**3*bm2)*inv%p
        bm2,bm1=bm1,b
        if b==0:
            z.append(n+1)
    return z


def compute_window(N):
    zs={p:apery_zero_set(p) for p in primerange(N+1,2*N+1)}
    hits={n:[] for n in range(N+1,2*N+1)}
    for p,Z in zs.items():
        for n in range(max(N+1,p),2*N+1):
            if n%p in Z:
                hits[n].append(p)
    B={n:sum(log(p) for p in ps if p>isqrt(n)) for n,ps in hits.items()}
    return B,hits,zs

for N in [1000,2000,4000,8000,10000]:
    B,hits,zs=compute_window(N)
    M2=sum(x*x for x in B.values())
    print('N=',N)
    print('M2=',M2)
    print('maxB=',max(B.values()))
    print('ratio=',M2/(N*log(N)**2))
    print('large events:')
    for n in sorted(B):
        if B[n]>=20:
            print(n,B[n],hits[n])
```

## Expected mathematical interpretation

The statistic is testing the second moment

\[
M_2(N)=\sum_{N<n\le2N}B(n)^2.
\]

If top-window prime hits behaved independently with probability approximately \(|Z_p|/p\), then

\[
E(B(n))\approx\sum_{\sqrt n<p\le n}\frac{|Z_p|\log p}{p},
\]

and the diagonal contribution gives a main term of order

\[
N(\log N)^2.
\]

The claimed constant near 0.59 would therefore be evidence for a stable second-moment constant rather than a proof of independence.

## Task 2 output

The same program prints all n with B(n)>=20 and the corresponding prime list. For each prime p, the complete zero set Z_p is available in the `zs` dictionary and the actual residue is simply `n % p`.

## Caveat about the requested O(p) recurrence

For N=10000 the naive method computes zero sets for about 1000 primes near 10000, each requiring about 10000 recurrence steps, roughly 10^7 modular operations, which is feasible in optimized Python. The full scan over n is smaller than the zero-set construction.

The output should be used as the numerical verification dataset for the 0.59 prediction. No numerical table is included here because it requires an actual run of the above code.