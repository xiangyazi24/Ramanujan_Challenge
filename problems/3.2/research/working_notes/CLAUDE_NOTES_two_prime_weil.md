# Two-prime Weil correlation — the potential proof path (2026-07-31)

## Discovery
Numerical computation reveals that both the VERTICAL complete exponential sum
and the TWO-PRIME shifted correlation of the Apéry sequence satisfy Weil-type bounds.

### Data 1: Vertical complete sum
```
C_p(h) = Σ_{a=0}^{p-1} e(h·b_a/p)
|C_p(1)| ≤ 3.27·√p   for ALL 166 primes p ≤ 1000
Average |C_p(1)|/√p = 1.26
```

### Data 2: Two-prime shifted correlation
```
Corr(p,q,d) = Σ_{m=0}^{M-1} e(b_m/p - b_{m+d}/q),   d = |p-q|, M = min(p,q)-|d|-1
|Corr(p,q,d)| ≤ 2.09·√M   for ALL 127 tested prime pairs (p,q ≤ 709)
102 pairs in [200,600] tested, max ratio 2.09, 1 pair over 2.0
```

## The proof chain (if the two-prime bound is proved)

1. **Two-prime Weil bound** (DATA 2): |Σ_m ψ_p(b_m)·ψ̄_q(b_{m+d})| ≤ C√M
2. → **4th-moment bound**: Σ_{N<n≤2N} |S_h(n)|^4 ≤ C'·N^3/log^2 N
   (cross terms in the 4th-moment expansion are O(N^{3/2}) per pair)
3. → **Pointwise bound**: max |S_h(n)| = O(N^{3/4}/log^{1/2}N) = o(N/log N)
4. → **Fejér kernel**: H(n) ≤ P_n/K + (1/K)Σ |S_j| = o(P_n)
5. → **The GCD conjecture**: log G_n = o(n)

## The precise missing theorem (from Q6261 + Q6262)

**Theorem (to prove).** For every pair of distinct large primes p, q and every fixed
shift d, there exists a lisse sheaf G_{p,q,d} on A^1 such that:

(i) Tr(Frob_m | G_{p,q,d}) = ψ_p(b_m)·ψ_q(-b_{m+d})

(ii) G is pointwise pure of weight 0

(iii) G is geometrically irreducible (monodromy acts irreducibly after base change)

(iv) cond(G) = rank + Σ_x(1 + Swan_x) = O(1) independent of p, q

Then Deligne's Riemann Hypothesis gives the needed bound.

## The construction route

The Apéry numbers b_n = Σ C(n,k)^2 C(n+k,k)^2 are the diagonal of a rational function
(Furstenberg-Deligne theory). The associated Picard-Fuchs/Gauss-Manin sheaf provides the
ℓ-adic object. The two-prime version would be a tensor product of two independent
Apéry period sheaves (one mod p, one mod q), with the shift d handled by a Tate twist or
translation functor.

The decisive step: prove the tensor product has no trivial component
(geometric irreducibility = "two-prime shifted independence of the Apéry motive").

Relevant technology: Katz "Exponential Sums and Differential Equations" (1990),
Katz "Gauss Sums, Kloosterman Sums, and Monodromy Groups" (1988),
Adolphson-Sperber toric exponential sums, Fu-Wan polynomial recursion sums.
