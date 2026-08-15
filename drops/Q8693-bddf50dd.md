ANSWER Q8693 bddf50dd

# Audit of the no-isolated-exception propagation route

## Verdict

The two unconditional inputs do not by themselves force a no-isolated-exception theorem. A row with large

\[
W(n)=\sum_{\sqrt n<p\le n,\ p\mid B_{n\bmod p}}\log p
\]

would have to create many nearby rows carrying a positive proportion of the same weight. The available Apéry identities give several exact transport maps, but none preserves enough of the large-prime support. The missing statement is a genuine cross-row correlation estimate.

## 1. What propagation would require

Fix
\[
S(n)=\{p:\sqrt n<p\le n,\ p\mid B_{n\bmod p}\}.
\]
Suppose \(W(n)\ge \epsilon n\). A propagation theorem with window \(H\) would need to prove that for some \(c>0\),

\[
\#\{1\le h\le H:W(n+h)\ge c\epsilon n\}\gg H.
\]

Since each prime contributes at most \(\log n\), this requires that a positive proportion of primes in \(S(n)\) survive under a controlled map
\[
p\mapsto p\mid B_{(n+h)\bmod p}.
\]

The required modulus-side statement is therefore:

For many \(h\),
\[
B_{n\bmod p}\equiv0\pmod p\implies B_{(n+h)\bmod p}\equiv0\pmod p
\]
for a positive proportion of the primes contributing to \(W(n)\).

No such implication follows from the recurrence.

## 2. The maps \(n\mapsto n+h\)

For fixed p, write \(r=n\bmod p\). Then
\[
(n+h)\bmod p=r+h\pmod p.
\]
The Apéry recurrence modulo p is a second-order recurrence, but knowing one zero value
\[
B_r=0\pmod p
\]
does not determine the neighboring values. The only universal restrictions are the known reflection and non-consecutive-root facts.

Reflection gives
\[
B_{p-1-r}=B_r\pmod p,
\]
so one root transports to one reflected root. It does not produce a positive-density cluster of roots.

Thus a large set of primes contributing at n may move to unrelated residues at n+h.

## 3. The digit lift \(n\mapsto n+tp\)

This is the strongest available exact transport. If
\[
n=pq+r,
\]
Lucas gives
\[
B_n\equiv B_qB_r\pmod p.
\]
Changing n by multiples of p fixes r:
\[
n+tp=p(q+t)+r.
\]
Therefore
\[
B_{n+tp}\equiv B_{q+t}B_r\pmod p.
\]
If \(B_r=0\), then every digit lift remains a zero modulo p.

However, these are rows separated by the prime-dependent spacing p. They do not create many ordinary nearby integers m. For a fixed interval length L there are only O(L/p) such lifts per prime, and summing over different p gives no guaranteed multiplicity because the lift lattices have different moduli.

Hence digit lifts preserve isolated support rather than destroy it.

## 4. Reflection transport

The reflection map sends
\[
r\longmapsto p-1-r.
\]
For a row n and prime p this corresponds to another row residue
\[
m\equiv p-1-(n\bmod p)\pmod p.
\]
The resulting rows depend on p. There is no common integer m receiving many reflected primes unless an additional coincidence theorem is proved.

This is exactly the same obstruction as the original problem: one needs simultaneous alignment of many prime-dependent affine conditions.

## 5. Why the codegree theorem does not imply propagation

The exceptional-set proof uses:

\[
\#\{p>N^\delta:p\mid B_{m\bmod p},\ p\mid B_{n\bmod p}\}=O_\delta(|m-n|).
\]

This is a pair-correlation upper bound. It says that two rows do not share too many bad primes.

Propagation requires the opposite type of statement: a large W(n) forces many rows to share many of those primes. The two statements are logically compatible. A sparse collection of rows may each have a large private prime set while still satisfying all pairwise codegree bounds.

## 6. Actual Apéry-compatible obstruction

The obstruction is not a recurrence-free incidence construction. It is the exact prime-dependent digit geometry of the Apéry congruence:

For each prime p, a root r of B modulo p gives a family
\[
\{n:n\equiv r\pmod p\}
\]
and its reflected family
\[
\{n:n\equiv p-1-r\pmod p\}.
\]
Different primes define different lattices. The recurrence supplies these lattices but supplies no alignment among them.

A hypothetical large W(n) can therefore arise from many primes whose root lattices all meet at the single integer n, while their next intersections occur at unrelated locations. The recurrence allows this because Lucas transport preserves each individual lattice independently.

## 7. Narrowest missing correlation theorem

To obtain the desired contradiction one needs a statement of the following form.

There exist constants c,C>0 such that for every n and every set of primes
\[
P\subseteq S(n),\qquad \sum_{p\in P}\log p\ge \epsilon n,
\]
there exists
\[
1\le h\le C
\]
(or a controlled range of h) with
\[
\sum_{p\in P\cap S(n+h)}\log p\ge c\epsilon n.
\]

Equivalently, one needs a positive-density common-root transfer theorem for the actual Apéry zero sets
\(Z_p\).

The existing reflection, Lucas, and recurrence identities do not imply it.

## Conclusion

No-isolated-exception propagation is not forced by the known exceptional-set theorem plus the Apéry recurrence. The obstruction is the absence of cross-prime synchronization: the recurrence controls each modulus p separately, while propagation requires correlation between the families of residue lattices attached to many different primes. The narrow missing input is an actual Apéry zero-set correlation theorem strong enough to transfer a positive fraction of \(W(n)\) to nearby rows. Without that input, isolated large-W rows remain compatible with all currently available identities.
