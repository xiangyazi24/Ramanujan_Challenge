ANSWER Q47 efc86dd9

# p-adic analysis of the Apéry GCD conjecture

## Executive summary

The p-adic viewpoint is the correct local language, but the naive idea that the base-p digits of n independently control all primes does not directly prove the GCD conjecture. The obstruction is that the profiles

\[
(v_p(b_n))_p
\]

live in different p-adic worlds and there is no known independence theorem between them.

The useful output of p-adic analysis is instead:

1. a precise description of single-prime divisibility;
2. lifting laws for valuations;
3. strong restrictions on repeated zeros of the Apéry orbit.

These give local sparsity, but the missing ingredient is still cross-prime anti-concentration.

## 1. The p-adic valuation mechanism

For p>n, Lucas gives

\[
b_n \equiv b_{n\bmod p}b_{\lfloor n/p\rfloor}\pmod p.
\]

Since \(\lfloor n/p\rfloor=0\),

\[
p\mid b_n \iff n\bmod p\in Z_p.
\]

Thus the top-window primes are exactly the primes whose first base-p digit lands in the Apéry zero set.

For p<n, iterating Lucas-type congruences gives a digit decomposition. Schematically,

\[
v_p(b_n) \approx 3\sum_i 1_{d_i(p)\in Z_p},
\]

where d_i(p) are the base-p digits of n.

The factor 3 comes from the Apéry supercongruence depth:

\[
b_{mp^r}\equiv b_{mp^{r-1}}\pmod {p^{3r}}.
\]

However, this is only a same-prime statement.

## 2. Can different p-adic profiles be coupled?

Suppose

\[
n\equiv r_i\pmod {p_i}.
\]

Then

\[
p_i\mid b_n\iff p_i\mid b_{r_i}.
\]

A k-atom therefore corresponds to simultaneous conditions

\[
r_i\in Z_{p_i}.
\]

The product

\[
\prod_i b_{r_i}
\]

does not help because

\[
\log |b_r|\sim r\log(17+12\sqrt2),
\]

so the available size is exponentially larger than the product of the forcing primes.

The right object would be a gcd or resultant between the different Apéry values, not the product.

## 3. Bounds on gcd(b_r,b_s)

The Casoratian identity is

\[
a_rb_{r+1}-a_{r+1}b_r=-\frac6{(r+1)^3}.
\]

Therefore if a prime p divides both b_r and b_{r+1}, then p divides the denominator obstruction, giving only finitely many possibilities:

\[
p\mid 6(r+1).
\]

So adjacent Apéry values have extremely small gcd.

For general distance h=s-r, the analogue is obtained by propagating the recurrence matrix:

\[
\binom{b_s}{b_{s-1}}=M_{r+h}\cdots M_r\binom{b_r}{b_{r-1}}.
\]

If p divides both b_r and b_s, then the projective orbit starting from the zero point returns to the zero divisor after h steps. Hence

\[
gcd(b_r,b_s)
\]

is controlled by periodicity of the Apéry recurrence orbit modulo p.

There is no known uniform bound like

\[
gcd(b_r,b_{r+h})\ll_h 1
\]

for arbitrary h.

The best available replacement is:

- fixed h gives an algebraic condition of bounded degree;
- large h requires orbit anti-concentration.

This is exactly the same obstruction appearing in the atom-tail formulation.

## 4. p-adic digits and simultaneous bases

A fixed integer n has incompatible digit expansions:

\[
n=\sum_i d_i^{(p)}p^i.
\]

The condition

\[
d_i^{(p)}\in Z_p
\]

is a different dynamical condition for every p. There is currently no theorem saying that these digit conditions behave independently over varying p.

The desired estimate

\[
\sum_{\sqrt n<p\le n}v_p(b_n)\log p=o(n)
\]

would follow from a cross-prime large sieve estimate of the form

\[
\#\{p_1,\dots,p_k:p_i\mid b_n\}
\]

having Poisson-type tails.

This is not supplied by the p-adic supercongruences alone.

## 5. What p-adic analysis does prove

The p-adic approach gives three powerful constraints:

### (a) Valuation depth is rare

For a fixed p, high valuation requires repeated p-adic digit hits. Thus

\[
v_p(b_n)\ge 3j
\]

forces j separate digit conditions.

### (b) Exceptional primes are structured

Primes causing high valuations are linked to the zero set Z_p and recurrence orbit returns. They are not random divisors.

### (c) Repeated divisibility forces orbit collisions

If many primes divide b_n, then many independent projective recurrence orbits hit the zero divisor at prescribed times. This converts the GCD problem into an orbit anti-concentration problem.

## Final conclusion

The p-adic digit perspective is valuable, but the missing step is not a single-prime valuation estimate. The conjecture requires controlling the simultaneous profile

\[
(v_p(b_n))_p.
\]

The natural next theorem would be a cross-prime p-adic large sieve:

\[
\sum_{n<N}\prod_{p\in S}1_{p\mid b_n}
\]

must factor approximately as the product of local densities.

The known supercongruences provide the local factors, but not the independence. A bound on \(gcd(b_r,b_s)\) for separated indices would be an important component, but the full Apéry GCD conjecture still requires controlling simultaneous recurrence-orbit hits across different primes.
