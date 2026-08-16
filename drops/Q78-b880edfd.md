ANSWER Q78 b880edfd

# Goldbach-type reformulation and cross-integer independence

## Main conclusion

This is a genuinely different direction from the original cross-prime problem. The reformulation is correct and important, but the last step is still a correlation theorem about Apéry values. It replaces a hard *same integer / many primes* problem by a potentially easier *many integers / one prime each* problem.

For the top window:

\[
p\in(m/2,m],\qquad p\mid b_m \iff p\mid b_{m-p}.
\]

Putting \(r=m-p\) gives

\[
\#\{p\in(m/2,m]:p\mid b_m\}
=
\#\{0\le r<m/2:m-r\text{ prime and }m-r\mid b_r\}.
\]

This is a Goldbach-type diagonal problem.

## Expected size

For a random integer of size about \(\alpha^r\),

\[
\Pr(m-r\mid b_r)\approx \frac1{m-r}.
\]

Therefore

\[
\sum_{r<m/2}\frac1{m-r}=\log 2+O(1/m),
\]

so the Poisson heuristic predicts O(1) hits in the top window.

The target bound is therefore not unreasonable: one would like

\[
\sum_{r<m/2}1_{m-r\mid b_r}\ll m^{o(1)}.
\]

## Why this is easier than the original problem

The original LGTQ obstruction asks about many primes dividing the same integer b_m. Here each event is attached to a different Apéry number:

\[
X_r=1_{m-r\mid b_r}.
\]

The pair correlation is

\[
X_rX_s=1_{m-r\mid b_r}1_{m-s\mid b_s}.
\]

The moduli are different and the underlying integers are different. This removes the strongest obstruction: simultaneous zeros of one recurrence orbit modulo many characteristics.

## Casoratian input

The identity

\[
a_nb_{n+1}-a_{n+1}b_n=-\frac6{(n+1)^3}
\]

implies

\[
\gcd(b_n,b_{n+1})\mid 6.
\]

Thus adjacent values have complete independence for primes \(>3\). This is a real structural advantage.

For distant pairs, however, a general theorem of the form

\[
\gcd(b_r,b_s)=\exp(o(r+s))
\]

or a stronger uniform bound is not currently available. Without such an input, a large common factor could create correlations.

## Second moment strategy

Define

\[
S(m)=\sum_{r<m/2}X_r.
\]

To prove \(S(m)\ll \log m\) by a second moment method, it is enough to control

\[
\mathbb E S^2
=\sum_rX_r+2\sum_{r<s}X_rX_s.
\]

The desired estimate is

\[
\sum_{r<s}X_rX_s\ll 1+o(\log^2m).
\]

A probabilistic model gives

\[
\Pr(X_rX_s=1)\approx\frac1{(m-r)(m-s)},
\]

and the double sum converges to a constant. The missing theorem is a rigorous replacement for this independence.

## What sieve theory can and cannot do

A Selberg sieve normally needs a family of residue conditions on one parameter set. Here the modulus itself changes with the parameter:

\[
r\mapsto m-r.
\]

So this is not a standard sieve sequence.

The right formulation is closer to a large-sieve estimate for the moving-modulus recurrence:

\[
\sum_{r\le M}1_{p_r\mid b_r},\qquad p_r=m-r.
\]

One would need cancellation in the sequence

\[
b_r\pmod{m-r}.
\]

This is precisely the kind of estimate not supplied by Bombieri--Vinogradov.

## Possible route through recurrence geometry

A promising approach is to study the joint distribution

\[
(b_r\bmod p, b_s\bmod q)
\]

for p,q comparable to m. Since the characteristics differ, one may hope for a product theorem using:

1. the transfer matrix description of the Apéry recurrence;
2. the Casoratian non-degeneracy;
3. monodromy/equidistribution for independent characteristics.

This is closer to a two-family independence statement than the original Mellin sheaf problem.

## Verdict

The Goldbach reformulation is one of the better bypasses discovered:

- It avoids the defining-characteristic condition \(T_p(r)=0\) for a fixed r across many p.
- It converts the problem into divisibility of many different integers.
- The adjacent correlation is killed completely by the Casoratian.

However, the remaining theorem needed is still nontrivial:

\[
\sum_{r<s\le m/2}
1_{m-r\mid b_r}1_{m-s\mid b_s}
\ll m^{o(1)}.
\]

This is not currently a known consequence of standard sieve theory. It would require a new cross-integer independence theorem for Apéry numbers, but it appears substantially more accessible than the original cross-prime Mellin problem.
