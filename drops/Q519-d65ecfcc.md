ANSWER Q519 d65ecfcc

# Kloosterman completion and bilinear inverse sums

## (1) Short interval inverse sums and Iwaniec--Kowalski

The stated bound

\[
\left|\sum_{X<n\le 2X,(n,p)=1} e_p(cn^{-1})\right|\ll \sqrt p\log p
\]

is a standard consequence of the completion method, but it is not usually stated as a standalone theorem in exactly this notation in Iwaniec--Kowalski.

The relevant ingredients in *Analytic Number Theory* (Iwaniec--Kowalski, AMS Colloquium Publications 53, 2004) are:

* the Weil bound for Kloosterman sums
\[
|S(a,b;p)|\le 2\sqrt p,
\]
covered in the discussion of Kloosterman sums and exponential sums in Chapter 11 (especially the sections on Kloosterman sums and the Weil bound);
* the completion method for incomplete exponential sums, discussed in the same chapter in the sections on incomplete sums.

The argument is:

\[
1_{(X,2X]}(n)=\frac1p\sum_{m\bmod p}\widehat w(m)e_p(mn),
\]

so after substitution

\[
\sum_n w(n)e_p(cn^{-1})
=\frac1p\sum_m\widehat w(m)
\sum_{n\bmod p^*}e_p(cn^{-1}+mn).
\]

The inner sum is the Kloosterman sum
\[
S(m,c;p),
\]
so Weil gives \(\ll\sqrt p\). The Fourier coefficients of the interval satisfy
\[
\frac1p\sum_m|\widehat w(m)|\ll \log p,
\]
which yields
\[
\ll \sqrt p\log p.
\]

So the exact statement is a corollary of the IK completion machinery rather than a named theorem with this precise short interval.

---

## (2) Bilinear prime inverse sum

Consider

\[
B=\sum_{q,r\in(X,2X]\atop q,r\ \mathrm{prime}}
 e_p(c(qr)^{-1}).
\]

The naive factorization idea is not quite correct because

\[
e_p(c(qr)^{-1})=e_p(cq^{-1}r^{-1})
\]

couples the two variables multiplicatively. However it is a standard bilinear Kloosterman-type sum.

A direct bound from the one-variable completion estimate is not enough. The trivial estimate gives
\[
|B|\le P^2,
\]
and applying the short interval estimate separately loses the prime structure.

The right framework is:

* bilinear forms with Kloosterman fractions;
* dispersion method;
* large sieve for trace functions.

References:

1. Deshouillers and Iwaniec, *Kloosterman sums and Fourier coefficients of cusp forms*, Invent. Math. 70 (1982). This develops bilinear forms involving Kloosterman sums.

2. Kowalski, *The Large Sieve and its Applications*, especially the chapters on trace functions and bilinear sums.

3. Fouvry, Kowalski, Michel papers on bilinear sums with trace functions provide the modern sheaf-theoretic version.

For arbitrary sets of primes q,r in a dyadic interval, the desired statement

\[
B=o(P^2)
\]

is plausible but is not an automatic consequence of Weil. One needs cancellation after restricting both variables to primes.

---

## Relevance to P3.2

For the Apéry GCD problem, the obstacle is stronger than proving a single bilinear estimate. The second moment of the Weyl error requires uniform control over many coupled prime variables.

The useful target would be a prime bilinear estimate of the form

\[
\sum_{q,r\sim X} \Lambda(q)\Lambda(r)
 e_p(c(qr)^{-1})
\ll X^{2-\delta}
\]

or an analogous estimate after Vaughan/Heath-Brown decomposition.

A proof would likely combine:

1. completion/Kloosterman trace interpretation;
2. prime decomposition via Vaughan identity;
3. bilinear trace-function estimates.

The short-interval IK bound is therefore a necessary local input, but it does not by itself solve the prime bilinear problem.

## Bottom line

* The \(2\sqrt p\log p\) completion bound is standard IK material: Weil Kloosterman bound + Fourier completion.
* It is not usually quoted as a separate theorem for inverse sums over intervals.
* The prime bilinear sum is a genuinely harder object; the correct references are the Deshouillers--Iwaniec dispersion/Kloosterman literature and modern Kowalski--Michel trace-function bilinear estimates.
* For P3.2, the missing ingredient is not the one-variable Kloosterman bound but a uniform bilinear (and eventually multilinear) trace-function estimate over the prime variables.
