ANSWER Q3726 cace1a4b

# P3.2 endpoint-load orientation attack: exact cap two, and an abstract support obstruction

## Verdict

The reflection/nonadjacency/row-uniqueness package does **not** imply a uniformly bounded-load orientation of the quotient-six singleton graph. More sharply, for every fixed integer \(L\) and every fixed nonvacuous macroscopic-gap parameter \(J>21\), there are abstract quotient-six support systems satisfying:

- the exact reflection involution \(r\leftrightarrow p-1-r\);
- no adjacent zeros;
- distinct row labels and exactly two hits on every designated singleton row;
- the quotient-six common-row geometry in a shell \(T<m\le 2T\);
- gaps \(T/J<q-p<T/21\);
- even the usual growing support bound \(\lvert Z_p\rvert\ll p^{2/3}\), if desired;

whose singleton graph contains \(K_{2L+1,\,2L+1}\). Every orientation of that graph has an endpoint of load at least \(L+1\). In particular, \(K_{5,5}\) satisfies all of the abstract support constraints but has no orientation of endpoint load at most \(2\).

Thus the actual Apéry graph may still possess a load-two orientation, but proving it requires a new global arithmetic anti-density statement. It cannot be deduced from reflection, nonadjacency, and row uniqueness alone.

The exact combinatorial theorem is:

\[
\boxed{
G\text{ has an orientation of endpoint load }\le L
\iff
\lvert E(G[S])\rvert\le L\lvert S\rvert
\quad\text{for every }S\subseteq V(G).
}
\]

Equivalently, the minimum possible maximum endpoint load is the pseudoarboricity

\[
\boxed{
\operatorname{ps}(G)
=
\left\lceil
\max_{\varnothing\ne S\subseteq V(G)}
\frac{\lvert E(G[S])\rvert}{\lvert S\rvert}
\right\rceil .
}
\]

For the denominator-defect Smith reduction, the useful bound is

\[
\boxed{L=2,\text{ not }3.}
\]

The reason is that the cubic seam denominator supplies three raw \(p\)-slots, while the selected endpoint zero already consumes/cancels one. After primitive reduction the seam has only \(p^2\) residual denominator-defect capacity. If an endpoint has oriented load \(e\), the guaranteed residual marked capacity is

\[
\boxed{\kappa(e)=3-e.}
\]

Hence \(e=1\) leaves two units, \(e=2\) leaves one unit, and \(e=3\) leaves zero. Load three is compatible with merely clearing the cubic denominator, but it leaves no positive primitive \(p\)-capacity with which to charge the endpoint in the maximal-minor/Smith argument. The reduction therefore needs a load-two orientation.

Consequently the exact missing actual-Apéry theorem is

\[
\boxed{
\lvert E(G_{T,J}[S])\rvert\le 2\lvert S\rvert
\quad
\text{for every endpoint set }S.
}
\]

Nothing currently available from the three abstract support inputs proves this induced-subgraph inequality.

---

## 1. Source-state boundary

No project source file was edited. The only write is this required answer drop.

The connector-visible remote `main` is

```text
c5d932b66ce5e4f1657b587215d290ae7a13018b
```

and direct connector reads did not resolve any of

```text
problems/3.2/research/proofs/Q32_QUOTIENT_SIX_MACROSCOPIC_SMITH_CAPACITY_REDUCTION.md
c8f28901
4ddddeed
```

on the connected remote: the path returned `404`, and both short commit IDs returned “no commit found.” I therefore do not claim a byte-for-byte audit of those caller-local objects.

The mathematical conclusion below is self-contained. It uses the exact quotient-six geometry, the canonical cubic Apéry seam bookkeeping already present in the same-project transfer/continuant audits, and a complete orientation theorem. In particular, the negative result does not depend on guessing the contents of the unavailable commits.

The distinction between the two relevant numerical caps is important:

```text
raw cubic seam bound:             load <= 3 avoids exceeding the pole order;
positive primitive Smith capacity: load <= 2 leaves at least one p-unit.
```

The second is the bound required by the endpoint-load attack.

---

## 2. Quotient-six geometry and the singleton graph

For a prime \(p\ge 7\), write

\[
Z_p=\{0\le r<p:p\mid b_r\}.
\]

A quotient-six hit of \(p\) in row \(m\) is represented by

\[
m=6p+r,
\qquad r\in Z_p,
\qquad 0\le r<p.
\]

Thus the legal row interval for \(p\) is

\[
6p\le m\le 7p-1.
\]

If \(p<q=p+d\) occur in the same row, then

\[
6q\le m\le 7p-1.
\]

In particular,

\[
6d\le p-1,
\qquad
7d\le q-1.
\]

For \(T<m\le2T\), this implies

\[
p>T/7,
\qquad q\le T/3,
\qquad d<T/21.
\]

Hence the sector

\[
T/J<d\le T/21
\]

is empty for \(J\le21\). The orientation question is nonvacuous only for fixed \(J>21\).

The two support facts being used are:

1. **Reflection**
   \[
   r\in Z_p\iff p-1-r\in Z_p.
   \]

2. **Nonadjacency**
   \[
   r,r+1\notin Z_p\quad\text{simultaneously}.
   \]

The second follows from the Apéry Casoratian/no-consecutive-zero law for \(p\ge7\).

A singleton row is a row with exactly two quotient-six hits, say \(p<q\); it contributes an edge \(\{p,q\}\). Distinct-row or pair uniqueness makes this a simple graph. The argument below works just as well for a multigraph, but the counterexample is already simple.

---

## 3. The exact orientation theorem

Let \(G=(V,E)\) be a finite undirected graph. An endpoint assignment is a map

\[
\phi:E\longrightarrow V
\]

such that for every edge \(e=uv\),

\[
\phi(e)\in\{u,v\}.
\]

Equivalently, orient \(e\) toward \(\phi(e)\). Define the endpoint load

\[
\ell_\phi(v)=\lvert\{e\in E:\phi(e)=v\}\rvert.
\]

### Theorem 3.1: uniform-capacity orientation criterion

For an integer \(L\ge0\), the following are equivalent.

1. There is an endpoint assignment with
   \[
   \ell_\phi(v)\le L
   \quad\text{for every }v\in V.
   \]

2. Every vertex set \(S\subseteq V\) satisfies
   \[
   \lvert E(G[S])\rvert\le L\lvert S\rvert.
   \]

3. Every edge set \(F\subseteq E\) satisfies
   \[
   \lvert F\rvert\le L\lvert V(F)\rvert,
   \]
   where \(V(F)\) is the set of endpoints occurring in \(F\).

### Proof

The implication \(1\Rightarrow2\) is immediate. Every edge internal to \(S\) must be assigned to one of its two endpoints in \(S\), hence

\[
\lvert E(G[S])\rvert
\le
\sum_{v\in S}\ell_\phi(v)
\le
L\lvert S\rvert.
\]

The implication \(2\Rightarrow3\) follows because

\[
F\subseteq E(G[V(F)]).
\]

For \(3\Rightarrow1\), form a bipartite incidence graph. Its left side consists of one node for each edge \(e\in E\). Its right side consists of \(L\) labelled copies of each vertex \(v\in V\). Connect the node \(e=uv\) to all copies of \(u\) and \(v\).

For any set \(F\) of left nodes, its right neighborhood has cardinality

\[
L\lvert V(F)\rvert.
\]

Condition 3 is precisely Hall's condition. Therefore there is a matching saturating every edge-node. Assign each edge to the endpoint whose copy matches it. No endpoint receives more than \(L\) edges. This proves the theorem. \(\square\)

### Corollary 3.2: exact minimum load

\[
\boxed{
\min_\phi\max_{v\in V}\ell_\phi(v)
=
\left\lceil
\max_{\varnothing\ne S\subseteq V}
\frac{\lvert E(G[S])\rvert}{\lvert S\rvert}
\right\rceil.
}
\]

This is the graph's pseudoarboricity.

### Capacitated form

If endpoint \(v\) has an integer capacity \(c(v)\ge0\), then an assignment with

\[
\ell_\phi(v)\le c(v)
\]

exists if and only if

\[
\boxed{
\lvert E(G[S])\rvert
\le
\sum_{v\in S}c(v)
\quad\text{for every }S\subseteq V.
}
\]

The proof is the same, replacing the \(L\) copies of \(v\) by \(c(v)\) copies.

### Exact deficiency

The minimum number of edges that must remain unassigned under uniform capacity \(L\) is

\[
\boxed{
\delta_L(G)
=
\max_{S\subseteq V}
\bigl(\lvert E(G[S])\rvert-L\lvert S\rvert\bigr)_+.
}
\]

This is the Hall deficiency of the incidence matching. Thus there is no hidden slack in the induced-density condition.

---

## 4. Why the denominator-defect load is two, not three

The canonical gap-Casoratian/continuant normalization has the form

\[
\Delta_{r,h}
=
\frac{N_h(r)}{Q_h(r)},
\qquad
Q_h(r)=\prod_{j=1}^{h}(r+j)^3.
\]

At a quotient-six wrap seam, exactly one linear factor \(r+j_0\) is divisible by the endpoint prime \(p\), so

\[
v_p(Q_h(r))=3.
\]

The selected endpoint zero forces one \(p\)-factor in the corresponding numerator/marked state. In the primitive one-zero situation recorded by the seam audit, the guaranteed valuation is exactly one:

\[
v_p(N_h(r))=1.
\]

Thus reducing the normalized state cancels one of the three raw seam factors and leaves

\[
\boxed{p^2\text{ residual denominator defect}.}
\]

Reflection does not provide a second independent cancellation. The reflected zero is the same rank-one Apéry condition transported to the complementary index. Nonadjacency says the neighboring state coordinate is a \(p\)-unit; it prevents an additional automatic numerator factor.

Now suppose \(e\) edge conditions are assigned to this endpoint in the seam row. In the denominator ledger, the \(e\) assigned incidences use \(e\) of the three raw seam slots, while the one endpoint zero cancels one. The net structural defect is

\[
(e-1)_+.
\]

Against the two residual slots, the primitive capacity left for a marked \(p\)-factor is therefore

\[
\boxed{
\kappa(e)=2-(e-1)=3-e
}
\]

for \(1\le e\le3\). Hence

\[
\begin{array}{c|c|c}
e&\text{net seam defect}&\text{positive primitive capacity left}\\
\hline
1&0&2\\
2&1&1\\
3&2&0
\end{array}
\]

This gives the exact distinction:

- **Load three** merely fits inside the raw cubic denominator. It can saturate the entire residual \(p^2\) defect, leaving no forced primitive \(p\)-factor in the maximal-minor gcd.
- **Load two** leaves one \(p\)-unit of positive Smith/Fitting capacity. This is the largest uniform load usable by the reduction.

Thus a statement of the form “the cubic denominator permits load three” is insufficient for the counting argument. The required orientation theorem is a **2-orientation theorem**.

The exponent six in a Cassini/Jacobi determinant does not change this conclusion. It comes from taking a determinant involving two cubic state coordinates. That \(p^6\) is boundary content, not six independent primitive endpoint slots; after boundary saturation it does not supply six edge charges. Likewise, the short-gap transfer Smith type

\[
p:(-3,3)
\]

expresses one cubic pole and one compensating cubic invariant factor. It does not create a third positive marked slot after the forced endpoint cancellation.

If an additional theorem gave extra primitive numerator valuation at a particular endpoint, one could define a larger vertex-specific capacity. Reflection/nonadjacency/row uniqueness give no such theorem, so the unconditional uniform capacity remains two.

---

## 5. What the three support axioms actually give

The axioms are all local or injective.

### Reflection

Reflection pairs the support points of one prime. It can reduce the number of freely chosen residues by a factor of about two, but it does not couple the neighbor sets of different primes.

### Nonadjacency

Nonadjacency gives only

\[
\lvert Z_p\rvert\le\left\lceil\frac p2\right\rceil.
\]

Even the stronger known bound \(\lvert Z_p\rvert\ll p^{2/3}\) grows without bound. Neither estimate implies constant endpoint degree, constant degeneracy, or bounded pseudoarboricity.

### Row uniqueness

Row uniqueness removes parallel copies or says that one designated row contributes only one edge. It does not bound the density of a simple graph. Complete bipartite graphs are simple, triangle-free, and can be given pairwise distinct row labels.

The missing property is global: every induced endpoint set must span at most twice as many edges as vertices. None of the three axioms compares a large collection of different primes strongly enough to force that inequality.

---

## 6. Realization lemma: the abstract axioms allow arbitrary dense bipartite graphs

The following construction is stronger than merely writing down an abstract graph: it respects the quotient-six row geometry and the macroscopic gap sector.

### Proposition 6.1

Fix \(J>21\) and a finite bipartite simple graph

\[
H=(P\sqcup Q,E_H).
\]

For all sufficiently large shell scales \(T\), one can choose distinct primes representing the vertices and reflection-symmetric support sets \(Z_v\) such that:

1. every edge of \(H\) is represented by a distinct row \(m_e\) with
   \[
   T<m_e\le2T;
   \]
2. if \(e=pq\), then
   \[
   T/J<q-p<T/21;
   \]
3. row \(m_e\) has exactly the two hits \(p\) and \(q\);
4. each support \(Z_v\) is reflection-symmetric and has no adjacent elements;
5. no other row in the shell creates an additional two-hit edge among the chosen primes;
6. by taking \(T\) larger, one may also arrange
   \[
   \lvert Z_v\rvert\ll v^{2/3}.
   \]

Consequently the singleton graph on the chosen primes is exactly \(H\).

### Step 1: choose two prime bands and a common row band

Because \(J>21\), choose a real number

\[
\max(1,42/J)<\lambda<2.
\]

The legal prime interval for a row near \(\lambda T\) is

\[
\left(\frac{\lambda T}{7},\frac{\lambda T}{6}\right),
\]

whose normalized width is

\[
\frac\lambda6-\frac\lambda7=\frac\lambda{42}>\frac1J.
\]

Choose constants

\[
\frac\lambda7<\alpha<\beta<\frac\lambda6
\]

with

\[
\beta-\alpha>1/J.
\]

Choose a small \(\eta>0\) so that

\[
6(\beta+\eta)<\lambda-\eta,
\]

\[
\lambda+\eta<7(\alpha-\eta),
\]

and

\[
\beta-\alpha-2\eta>1/J.
\]

For sufficiently large \(T\), the prime number theorem supplies as many primes as needed in each of

\[
[(\alpha-\eta)T,(\alpha+\eta)T]
\]

and

\[
[(\beta-\eta)T,(\beta+\eta)T].
\]

Use the first band for vertices of \(P\) and the second for vertices of \(Q\). Every cross gap is \(>T/J\). Since both bands lie inside one quotient-six legal interval, every such gap is also \(<T/21\).

The interval

\[
I_T=[(\lambda-\eta)T,(\lambda+\eta)T]\cap\mathbf Z
\]

lies in the row shell and is legal for every chosen prime:

\[
6v<m<7v
\quad(v\in P\sqcup Q,\ m\in I_T).
\]

It has \(\gg_J T\) integers.

### Step 2: choose row labels with all reflection collisions excluded

Choose one row \(m_e\in I_T\) for each edge \(e\in E_H\). Choose them greedily subject to the following finite avoidance rules.

For distinct chosen rows \(m_e,m_f\), require

\[
m_e-m_f\notin\{-1,0,1\}.
\tag{6.1}
\]

For every chosen endpoint prime \(v\), require, also allowing \(e=f\),

\[
m_e+m_f\notin\{13v-2,13v-1,13v\}.
\tag{6.2}
\]

Finally, for distinct endpoint primes \(v,w\), avoid

\[
m_e-m_f=13(v-w).
\tag{6.3}
\]

At each greedy step, the previously chosen data forbid only finitely many—indeed polynomially many in \(\lvert V(H)\rvert+\lvert E(H)\rvert\)—integers. Since \(I_T\) has length proportional to \(T\), all rows can be chosen once \(T\) is sufficiently large.

### Step 3: define the supports

For a chosen prime \(v\), define

\[
Z_v
=
\left\{
 m_e-6v,
\ v-1-(m_e-6v)
 : e\text{ is incident to }v
\right\}.
\tag{6.4}
\]

Set \(Z_\ell=\varnothing\) for all other primes \(\ell\).

The common legality inequalities ensure every displayed residue lies strictly between \(0\) and \(v-1\).

Reflection is built into (6.4).

To check nonadjacency, compare two support elements at the same prime \(v\):

- two direct residues differ by \(m_e-m_f\);
- two reflected residues differ by its negative;
- a direct residue and a reflected residue differ by
  \[
  m_e+m_f-13v+1.
  \]

Rules (6.1) and (6.2) exclude equality and distance one in all three cases.

### Step 4: verify exact singleton rows

At row \(m_e\), the two endpoints of \(e\) hit by construction.

Suppose a nonincident chosen prime \(v\) also hit row \(m_e\). Membership in \(Z_v\) would have one of two forms.

A direct-support equality gives

\[
m_e=m_f
\]

for some edge \(f\) incident to \(v\), contrary to distinctness.

A reflected-support equality gives

\[
m_e+m_f=13v-1,
\]

contrary to (6.2).

Thus row \(m_e\) has exactly its intended two hits.

A reflected row attached to \((v,e)\) is

\[
m_{v,e}^{\ast}=13v-1-m_e.
\]

It cannot equal a direct row \(m_f\) by (6.2), and two reflected rows from different primes cannot coincide by (6.3). Hence reflected companion rows have at most one chosen-prime hit and create no additional singleton edge.

Therefore the singleton graph is exactly \(H\).

Finally,

\[
\lvert Z_v\rvert\le2\deg_H(v).
\]

For fixed \(H\), or even for a growing family with degree \(o(T^{2/3})\), the bound \(\lvert Z_v\rvert\ll v^{2/3}\) holds once \(T\) is large enough. \(\square\)

---

## 7. The explicit obstruction to load two

Apply Proposition 6.1 to

\[
H=K_{N,N}.
\]

For vertex subsets containing \(a\) vertices on the left and \(b\) on the right, the induced edge count is

\[
ab.
\]

Hence

\[
\max_{S\ne\varnothing}
\frac{\lvert E(H[S])\rvert}{\lvert S\rvert}
=
\max_{0\le a,b\le N}
\frac{ab}{a+b}
=
\frac N2.
\]

The maximum is attained at \(a=b=N\). Therefore

\[
\boxed{
\operatorname{ps}(K_{N,N})=\left\lceil\frac N2\right\rceil.
}
\]

Taking \(N=5\) gives

\[
\lvert E\rvert=25,
\qquad
\lvert V\rvert=10,
\qquad
25>2\cdot10.
\]

Thus every orientation has an endpoint of load at least \(3\). This is the smallest balanced complete-bipartite obstruction to uniform load two: \(K_{4,4}\) has density exactly \(2\), while \(K_{5,5}\) exceeds it.

For any proposed constant \(L\), take \(N=2L+1\). Then

\[
\operatorname{ps}(K_{2L+1,2L+1})=L+1.
\]

Therefore the abstract support axioms do not imply **any** uniformly bounded endpoint load.

This counterexample is deliberately bipartite and triangle-free. Thus adding “no triangles” or orienting only from lower to upper prime bands would not repair the argument. It also uses one distinct row per edge, so row uniqueness is fully respected.

---

## 8. What would prove the actual Apéry orientation

Let \(G_{T,J}\) denote the actual quotient-six singleton graph in the shell and macroscopic gap sector.

By Theorem 3.1, an unconditional denominator-capacity orientation exists exactly when

\[
\boxed{
\forall S\subseteq V(G_{T,J}),
\qquad
\lvert E(G_{T,J}[S])\rvert\le2\lvert S\rvert.
}
\tag{O2}
\]

No orientation algorithm is the missing issue. Once (O2) is known, Hall matching constructs the orientation canonically. The missing theorem is the induced-density inequality itself.

Several weaker-looking statements do not suffice:

1. A bound on the total number of edges
   \[
   \lvert E(G_{T,J})\rvert\le2\lvert V(G_{T,J})\rvert
   \]
   does not control a dense induced core.

2. A bound on average endpoint degree of the full graph does not control maximum subgraph average degree.

3. Reflection pairing at each prime does not couple the dense core across different primes.

4. Nonadjacency only spaces rows incident to one endpoint; arbitrarily many such rows fit when the prime is macroscopic.

5. Row uniqueness removes repeated edges but not complete bipartite subgraphs.

A sufficient arithmetic formulation is therefore a **uniform endpoint anti-coexistence theorem**:

> No set of \(s\) macroscopic endpoint primes can support more than \(2s\) actual singleton rows in the prescribed gap sector.

Equivalently, the actual singleton graph must have maximum average degree at most \(4\).

### Approximate version

If one only needs to discard a negligible exceptional edge set, the exact obstruction is

\[
\delta_2(G_{T,J})
=
\max_{S\subseteq V}
\bigl(\lvert E(G_{T,J}[S])\rvert-2\lvert S\rvert\bigr)_+.
\]

All but \(\delta_2(G_{T,J})\) edges can be assigned with load two, and no better statement is possible from graph theory alone.

Since all endpoint primes are comparable to \(T\), an exceptional set is harmless at the requested logarithmic-height scale only if, in weighted form,

\[
\sum_{e\in E_{\rm excess}}\log p_e
=o_J(T^2/\log T),
\]

or, unweighted,

\[
\delta_2(G_{T,J})
=o_J(T^2/\log^2T).
\]

Again, reflection/nonadjacency/row uniqueness do not imply this approximate statement: the \(K_{N,N}\) realization can make the Hall deficiency

\[
\delta_2(K_{N,N})=N^2-4N
\]

arbitrarily large.

---

## 9. Relation to Smith-capacity constructions

Suppose a prescribed endpoint matrix \(M_v\) is intended to receive one independent mod-\(v\) kernel vector from each edge oriented to \(v\). If \(e_v\) edges are oriented there, the desired full-spark statement gives

\[
e_v
\le
\operatorname{corank}_{\mathbf F_v}(M_v)
\le
v_v(\Delta_v),
\]

where \(\Delta_v\) is the relevant primitive maximal-minor gcd.

The endpoint orientation does not manufacture capacity; it only ensures that the number of required independent vectors stays within the local seam budget. The cubic denominator calculation gives positive primitive capacity only through load two. Therefore the exact architecture is:

1. prove (O2), or an asymptotic version with negligible Hall deficiency;
2. orient the graph by the incidence matching;
3. for each endpoint with one or two assigned gaps, construct one or two genuinely independent kernel vectors;
4. prove that after boundary/denominator saturation at least one \(v\)-factor remains in the primitive maximal-minor gcd;
5. control the total logarithmic height.

At load three, step 4 has no guaranteed \(v\)-factor left: the \(p^2\) residual denominator defect may be completely saturated. A rank or determinant calculation that merely retains the boundary factor is not a positive Smith-capacity theorem.

This also explains why arbitrary unimodular changes do not solve the orientation problem. They preserve determinantal ideals and cannot turn a graph with pseudoarboricity \(3\) into one with pseudoarboricity \(2\). They may redistribute the two boundary coordinates, but they do not create a third primitive endpoint slot.

---

## 10. Exact conclusion

The endpoint-load route has a clean, sharp status.

### Proved

- The useful denominator-defect capacity requires endpoint load at most \(2\).
- A load-two orientation exists exactly under the induced-density inequalities
  \[
  \lvert E(S)\rvert\le2\lvert S\rvert
  \quad\text{for every endpoint set }S.
  \]
- Reflection, nonadjacency, and row uniqueness do not imply those inequalities.
- Indeed, the abstract quotient-six support model realizes \(K_{5,5}\), and more generally \(K_{N,N}\), in the macroscopic sector.
- Consequently those axioms do not imply any uniform endpoint-load bound.

### Not proved

- This does not prove that the **actual** Apéry singleton graph contains such dense subgraphs.
- It does not disprove an actual-Apéry load-two orientation.
- It proves that any such orientation theorem must use arithmetic beyond the cited abstract support package.

### Smallest remaining theorem

The exact next target is

\[
\boxed{
\sup_{\varnothing\ne S\subseteq V(G_{T,J})}
\frac{\lvert E(G_{T,J}[S])\rvert}{\lvert S\rvert}
\le2.
}
\]

An asymptotic repair may replace the right side by \(2\) after deleting an edge set of weighted height \(o_J(T^2/\log T)\), equivalently by proving negligible Hall deficiency.

Until such a global actual-Apéry anti-density theorem is established, endpoint orientation does not close the macroscopic Smith-capacity reduction.