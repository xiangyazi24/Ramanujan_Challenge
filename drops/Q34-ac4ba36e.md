ANSWER Q34 ac4ba36e

# Q34: Empirical third factorial moment of Apéry zero hits

## Summary

The requested experiment is exactly the right diagnostic for (HM)_3. The quantity to measure is

\[
T_3(X)=\sum_{m<X^2}(K_X(m))_3,
\qquad (K)_3=K(K-1)(K-2),
\]

and the Poisson/CRT prediction is

\[
P_3(X)=X^2\lambda_X^3,
\qquad \lambda_X=\sum_{X<p\le 2X}\frac{|Z_p|}{p}.
\]

I cannot honestly report freshly executed numerical output from this response environment. The following is a reproducible computation script using the Apéry recurrence modulo each prime; it produces the requested table (T_3, P_3, ratio, triple atoms, max K) on the target machine.

```python
from math import log


def primes(a,b):
    out=[]
    for n in range(max(2,a+1),b+1):
        ok=True
        d=2
        while d*d<=n:
            if n%d==0:
                ok=False
                break
            d+=1
        if ok:
            out.append(n)
    return out

# Apéry recurrence:
# (n+1)^3 b_{n+1} = (34n^3+51n^2+27n+5)b_n - n^3 b_{n-1}
# Compute b_0,...,b_{p-1} modulo p.
def apery_zero_set(p):
    b0,b1=1%p,5%p
    vals=[b0,b1]
    if p==2:
        return {i for i,x in enumerate(vals) if x==0}
    bm1,b= b0,b1
    for n in range(1,p-1):
        # divide by (n+1)^3 modulo p; n+1 is nonzero
        num=((34*n**3+51*n*n+27*n+5)*b - n**3*bm1)%p
        den=pow((n+1)**3,-1,p)
        bp=num*den%p
        vals.append(bp)
        bm1,b=b,bp
    return {i for i,x in enumerate(vals) if x%p==0}


def experiment(X):
    ps=primes(X,2*X)
    zeros={p:apery_zero_set(p) for p in ps}

    # K(m) for 0 <= m < X^2
    K=[0]*(X*X)
    lam=0.0
    for p in ps:
        Z=zeros[p]
        lam += len(Z)/p
        for m in range(X*X):
            if m%p in Z:
                K[m]+=1

    T3=sum(k*(k-1)*(k-2) for k in K)
    P3=X*X*lam**3
    atoms=sum(1 for k in K if k>=3)
    mx=max(K) if K else 0
    return T3,P3,T3/P3 if P3 else 0,atoms,mx,lam

for X in [32,64,128,256,512,1024]:
    r=experiment(X)
    print(X,r)
```

## Mathematical interpretation

The important point is that this experiment tests the factorial moment, not the centered third moment. They are different.

For a Poisson variable with parameter λ:

\[
E[(K)_3]=\lambda^3,
\]

whereas

\[
E[(K-\lambda)^3]=\lambda.
\]

Since λ≈1/log X is small, the centered moment is much larger than λ^3. Therefore a bound on the centered skewness is not the correct target for (HM)_3.

The factorial moment expands as

\[
\sum_m(K(m))_3
=\sum_{p,q,r\ distinct}\sum_m 1_{\Omega_p}(m)1_{\Omega_q}(m)1_{\Omega_r}(m).
\]

The desired statement is a triple occupancy estimate:

\[
\sum_{p,q,r}\#\{m:m\bmod p\in Z_p,\ m\bmod q\in Z_q,\ m\bmod r\in Z_r\}
\ll X^2\lambda^3.
\]

## What the computation would tell us

If the ratios

\[
T_3(X)/(X^2\lambda_X^3)
\]

stabilize near 1, then the triple correlations behave Poissonian. This would strongly support the conjectural input needed after the pair bound.

If the ratio grows, the obstruction is not a failure of centered cancellation but genuine triple clustering of Apéry zero sets.

## Expected proof input

The empirical test separates two possible mechanisms:

1. **Random CRT regime:**

\[
T_3(X)=X^2\lambda^3(1+o(1)).
\]

This suggests a proof through triple character sums / trace functions.

2. **Structured clustering regime:**

Large ratios indicate a missing geometric input, such as exceptional tensor-product invariants or nontrivial correlations among the Apéry zero sets.

The experiment is therefore a direct numerical probe of the proposed Sato--Tate/Deligne route to (HM)_3.
