# Q8369 BRJ4 hostile audit

## Verdict

**BRJ is not a new obstruction beyond the existing one-`R`/CRT wall.**  On every chart on which the first-order block congruence is valid, the BRJ congruence is exactly the quotient of that congruence by the defining prime, followed by the tautological Vandermonde replacement

\[
\frac{R}{p_i}\equiv \prod_{j\ne i}(p_j-p_i)=V_i\pmod{p_i}.
\]

For the direct `g=1` chart, if the project p-adic lift is

\[
b_{p+h}\equiv 5b_h+10pD_h\pmod{p^2},
\]

then the sign in Q8360 is internally correct:

\[
\Gamma_i:=\frac{b_m-5b_{h_i}}{p_i}\equiv +10D_{h_i}\pmod{p_i},
\qquad
BV_i-5z_i-\Gamma_i\equiv0\pmod{p_i}.
\]

The small height of `V_i` does not help by itself: `B`, `z_i`, and the unreduced quotient represented by `Gamma_i` have no corresponding small-height control.

The proposed pre-CRT face/Boolean elimination also saturates in the natural local model.  Away from primes dividing the chart coefficient (`p_i != 5` for `g=1`), the BRJ equation is linear with a unit coefficient in `z_i`; consequently `z_i` can absorb every value of `B` and every recurrence/face value of `Gamma_i`.  This remains true after adjoining arbitrary canonical `(b,D)` recurrence and face equations **provided those equations do not themselves constrain the first p-adic quotient `z_i=b_{h_i}/p_i` (or the global quotient `B`) modulo `p_i`**.  A new theorem must control exactly that missing quotient digit.

## 0. Provenance / live-tree limitation

The requested files

- `/tmp/gpt/rc/Q8360.md`,
- `problems/3.2/isolated_five_star.tex`,
- `problems/3.2/ORACLE_COMM/codex_cube_boundary_attack.md`,
- `problems/3.2/ORACLE_COMM/codex_canonical_augmented_rectangle_audit.md`, and
- `problems/3.2/ORACLE_COMM/codex_topstrip_full_packet_tesseract.md`

are not exposed by the allowed GitHub connector in the connector-visible repository state.  The visible `main` used for this audit is commit
`734a5a84c1e4fd8703a811aadaa2b4c7f532b20e`; it contains `problems/3.2/ORACLE_COMM/` and the established P3.2 recurrence/Lucas material, but not the named newer BRJ files.  Library/conversation search also did not recover the exact Q8360 text.  I therefore do **not** claim to have source-verified any newer convention that appears only in those local/unpushed files.

What is audited below is (a) exact algebra from the formulas supplied in Q8369, (b) stable repository-visible Apéry recurrence/Lucas facts, and (c) finite exact recurrence evidence.  In particular, `Gamma=+10D` is certified below as the consequence of the displayed `+10pD` direct lift; its exact `D` convention could not be independently checked against the missing live file.

## 1. Direct `g=1` BRJ: every sign and denominator

Let the packet have distinct defining primes `p_i`, common direct-chart value

\[
m=p_i+h_i,
\]

and let

\[
R=\prod_i p_i,\qquad B=b_m/R,\qquad z_i=b_{h_i}/p_i.
\]

### 1.1 Integrality

At a defining node, `p_i | b_{h_i}`, so `z_i` is an integer.  The Apéry-Lucas direct block relation gives

\[
b_{p_i+h_i}\equiv b_1b_{h_i}=5b_{h_i}\equiv0\pmod{p_i}.
\]

Hence every `p_i | b_m`.  Distinct primes are pairwise coprime, so `R | b_m` and `B` is an integer.  Also `p_i | b_m-5b_{h_i}`, so the integer quotient

\[
Q_i:=\frac{b_m-5b_{h_i}}{p_i}
\]

exists and `Gamma_i` is its residue modulo `p_i`.

No inverse of `p_i` is being taken modulo `p_i`; the division is in `Z` **before** reduction.

### 1.2 Exact BRJ derivation

Using `b_m=RB` and `b_{h_i}=p_i z_i`, divide the integer identity

\[
b_m-5b_{h_i}=p_iQ_i
\]

by `p_i`:

\[
\frac{R}{p_i}B-5z_i=Q_i.
\]

Reducing modulo `p_i` gives

\[
\frac{R}{p_i}B-5z_i-\Gamma_i\equiv0\pmod{p_i}.
\]

Now, with Q8369's orientation

\[
V_i=\prod_{j\ne i}(p_j-p_i),
\]

each factor satisfies `p_j-p_i == p_j (mod p_i)`, hence

\[
V_i\equiv\prod_{j\ne i}p_j=R/p_i\pmod{p_i}.
\]

Therefore

\[
\boxed{BV_i-5z_i-\Gamma_i\equiv0\pmod{p_i}}.
\]

There is no additional arithmetic content in this last replacement.

### 1.3 Sign audit

With the stated orientation `p_j-p_i`, **there is no sign** in `R/p_i == V_i (mod p_i)`.

On the direct `g=1` chart, `p_i=m-h_i`; therefore

\[
p_j-p_i=(m-h_j)-(m-h_i)=h_i-h_j,
\]

so

\[
V_i=\prod_{j\ne i}(h_i-h_j).
\]

If a face formula instead uses `prod(h_j-h_i)`, it differs by `(-1)^(N-1)`.  For a 16-node packet, `N-1=15`, so that alternate orientation is the **negative** of Q8369's `V_i`.

For the p-adic correction, the sign is equally rigid.  If

\[
b_{p+h}-5b_h\equiv +10pD_h\pmod{p^2},
\]

then dividing by `p` gives `Gamma == +10D_h (mod p)`.  A source theorem with a minus sign would flip `Gamma`; it cannot be repaired by changing `V_i` under the definition used here.

### 1.4 Unit audit

Distinct packet primes imply `p_i` does not divide `p_j-p_i`, so `V_i` is a unit modulo `p_i`.  The coefficient of `z_i` is `5`; it is a unit unless `p_i=5`.  Thus the local freeness/no-go below applies to all `g=1` top-strip primes except the isolated exceptional prime `5`.

### 1.5 Height audit

If all nodes lie in an interval of diameter `H`, then on the direct chart

\[
|V_i|=\prod_{j\ne i}|h_i-h_j|\le H^{N-1};
\]

for `N=16`, this is `H^15`.

That is a correct small-height statement about `V_i`, but it does **not** make BRJ a small integer identity.  The congruence also contains

- `B=b_m/R`, whose archimedean size is not bounded by the node diameter;
- `z_i=b_{h_i}/p_i`, likewise uncontrolled by the Vandermonde height; and
- `Gamma_i`, which is only a residue of the potentially huge integer quotient `(b_m-5b_h)/p_i`.

Choosing least residues for those terms gives numbers `<p_i`, not a target-selective integer of height `<R`.  That is exactly the old CRT obstruction in different coordinates.

## 2. Precise pre-CRT saturation/no-go

The right statement is local and algebraic, and it includes arbitrary face eliminants and the canonical recurrence so long as they do not contain the missing first quotient digit.

### Theorem (BRJ unit-linear saturation)

Fix a node `i`, put `k=F_{p_i}`, and let `A` be any commutative `k`-algebra generated by the local packet/face/recurrence variables.  Thus `A` may already include all canonical `(b,D)` recurrence relations, Boolean face relations, and the direct relation `Gamma_i=10D_{h_i}`.  Let `Gamma` be the resulting element of `A`.  Let `c` be the chart coefficient (`c=5` for direct `g=1`) and assume `c != 0` in `k`.

Then for any `V in A`,

\[
A[B,z]/(VB-cz-\Gamma)\;\cong\;A[B],
\]

via

\[
z\mapsto c^{-1}(VB-\Gamma).
\]

Consequently the BRJ relation adds **no eliminant in `A[B]` at all**.  Equivalently, if `E` is the ideal of all face/recurrence equations before BRJ, then

\[
(E,VB-cz-\Gamma)\cap k[\text{aux},B]=E\cap k[\text{aux},B].
\]

#### Proof

Because `c` is a unit, the relation is monic linear in `z` after multiplication by `-c^{-1}`.  Substitution `z=c^{-1}(VB-Gamma)` defines inverse homomorphisms between the displayed quotient and `A[B]`.  This is an isomorphism, not a dimension count or heuristic.

### CRT corollary

For distinct packet primes,

\[
\mathbb Z/R\mathbb Z\cong\prod_i\mathbb F_{p_i}.
\]

Thus the residue of the single integer `B` is an arbitrary tuple `(B_i)_i` modulo `R`.  In each factor, if `c` is a unit, `z_i` can be chosen uniquely as

\[
z_i=c^{-1}(V_iB_i-\Gamma_i).
\]

So the phrase “the same `B` occurs on every face” does not create a modular coupling: CRT is precisely the statement that those residues are independent.

### What the canonical `(b,D)` recurrence does and does not change

The recurrence can determine the zeroth p-adic digit `b_h mod p` and auxiliary data such as `D_h mod p`; on the direct chart it may therefore determine `Gamma_i`.  But the quotient definition

\[
b_{h_i}=p_i z_i
\]

reduces modulo `p_i` to `b_{h_i}=0` and loses `z_i`.  To constrain `z_i mod p_i`, one needs one more p-adic digit: a relation modulo `p_i^2` (or an equivalent integral identity) in which `z_i` survives with nonzero coefficient.

Therefore adjoining the canonical recurrence **at residue level** does not defeat the theorem.  Adjoining a genuinely p-adic lifted recurrence that contains `z_i`, or an exact integral elimination performed before reduction that produces such a relation, is outside the no-go—and is exactly the kind of new theorem needed.

### Scope and exceptions

This is deliberately not a global impossibility theorem.

It does **not** cover:

1. primes `p_i | c` (for `g=1`, `p_i=5`), where the `z_i` coefficient vanishes;
2. any new p-adic relation involving `z_i` or `B` itself;
3. an archimedean height theorem coupling the CRT representatives;
4. a cross-prime integral identity formed before reduction whose content is not generated by the local residue ideals.

Within the stated local/face/recurrence scope, however, it is exact: every face eliminant that leaves `z_i` out is saturated away by the unit-linear BRJ equation.

## 3. Full fixed-`g`, fixed-branch packet, including reflection

The robust form uses the block coefficient

\[
c_g:=b_g.
\]

### 3.1 Direct branch

For

\[
m=g p_i+h_i,\qquad 0\le h_i<p_i,
\]

Apéry-Lucas gives

\[
b_m\equiv b_gb_{h_i}=c_gb_{h_i}\pmod{p_i}.
\]

At a defining node `p_i | b_{h_i}`, hence `p_i | b_m`.  Define

\[
\Gamma_i^{(g,+)}:=\frac{b_m-c_gb_{h_i}}{p_i}\pmod{p_i}.
\]

Exactly the same division gives

\[
\boxed{BV_i-c_gz_i-\Gamma_i^{(g,+)}\equiv0\pmod{p_i}}.
\]

The exact node-difference identity is

\[
p_i=\frac{m-h_i}{g},\qquad
p_j-p_i=\frac{h_i-h_j}{g},
\]

hence, for an `N`-node packet,

\[
\boxed{g^{N-1}V_i=\prod_{j\ne i}(h_i-h_j)}.
\]

This cross-multiplied form is safer than writing `g^{-(N-1)}`: it has no hidden modular denominator and remains a literal integer identity.

### 3.2 Reflected branch

Write the reflected base digit as `p_i-1-h_i`.  Then

\[
m=g p_i+(p_i-1-h_i)=(g+1)p_i-1-h_i.
\]

Using the project palindrome `b_{p-1-h} == b_h (mod p)` together with Lucas,

\[
b_m\equiv b_gb_{p_i-1-h_i}\equiv c_gb_{h_i}\pmod{p_i}.
\]

Define

\[
\Gamma_i^{(g,-)}:=\frac{b_m-c_gb_{h_i}}{p_i}\pmod{p_i}.
\]

Then the same quotient algebra survives:

\[
\boxed{BV_i-c_gz_i-\Gamma_i^{(g,-)}\equiv0\pmod{p_i}}.
\]

But the geometry orientation changes:

\[
p_i=\frac{m+h_i+1}{g+1},\qquad
p_j-p_i=\frac{h_j-h_i}{g+1},
\]

so

\[
\boxed{(g+1)^{N-1}V_i=\prod_{j\ne i}(h_j-h_i)}.
\]

For `N=16`, this is the negative of `prod(h_i-h_j)` after the common positive scaling.

### 3.3 What does **not** survive automatically

`Gamma_i=10D_{h_i}` is a **direct `g=1` p-adic correction**, not a chart-free identity.  For general `g`, and especially for the reflected digit, the first p-adic correction can contain additional `g`-dependent and reflection-dependent terms.  The BRJ skeleton survives with `Gamma_i^{(g,branch)}` defined by the quotient; the explicit `10D` replacement must be re-proved from the appropriate `p^2` block theorem.

Similarly, the unit-linear no-go uses `p_i not dividing c_g=b_g`.  Primes dividing `b_g` are exceptional local factors and must be split off before applying the saturation theorem.

### 3.4 Height in the two branches

If the `h_i` have diameter at most `H`, then

\[
|g^{N-1}V_i|\le H^{N-1}
\]

on the direct branch and

\[
|(g+1)^{N-1}V_i|\le H^{N-1}
\]

on the reflected branch.  Again, this is only a Vandermonde-height statement; it does not bound `B` or the first quotient `z_i`.

## 4. The smallest additional theorem that would close the packet

BRJ with unit `V_i` is equivalent to

\[
B\equiv V_i^{-1}(c_gz_i+\Gamma_i)\pmod{p_i}.
\]

So the exact missing information is the common target of the first-quotient combination `c_g z_i+Gamma_i`.

### Conditional target `[WFQA16]` (weighted first-quotient alignment)

For every actual 16-node fixed-`g`, fixed-branch top-strip packet, produce an integer `Q_T`, defined from the packet/face data independently of choosing CRT representatives, such that for every node

\[
\boxed{c_gz_i+\Gamma_i\equiv V_iQ_T\pmod{p_i}}.
\]

Then BRJ gives `B == Q_T (mod p_i)` for every `i`, hence

\[
\boxed{R\mid(B-Q_T)}.
\]

A quantified height statement

\[
\boxed{0<|B-Q_T|<R}
\]

then contradicts that divisibility and kills the packet.

This is the smallest target at the BRJ level because it specifies exactly the CRT residues that are currently free and nothing else.  A face theorem that merely makes `V_i` small does not do this.

### Zero-target special case

Taking `Q_T=0` asks for

\[
c_gz_i+\Gamma_i\equiv0\pmod{p_i}\quad\text{for every }i.
\]

Then `p_i | B` for all nodes, so `R | B` and therefore

\[
\boxed{R^2\mid b_m}.
\]

On direct `g=1`, if `Gamma_i=10D_{h_i}` and `p_i!=5`, this local condition is

\[
5z_i+10D_{h_i}\equiv0\pmod{p_i},
\]

equivalently

\[
z_i+2D_{h_i}\equiv0\pmod{p_i},
\]

or, in integral divisibility form,

\[
\boxed{p_i^2\mid b_{h_i}+2p_iD_{h_i}}.
\]

If one also proves `0<B<R` (equivalently `0<b_m<R^2`), this immediately contradicts `R|B`.

### TP16 gateway formulation

For the weighted full-packet gateway it is enough to prove, for every actual tesseract counted by the combinatorial TP16 theorem, a packet scalar

\[
S_T:=B-Q_T
\]

such that

\[
R_T\mid S_T,\qquad 0<|S_T|<R_T.
\]

`[WFQA16]` is a direct BRJ-level sufficient theorem producing precisely that scalar.  If the desired TP16 argument permits an exceptional subfamily, the same statement can be weakened quantitatively by proving it for all but `o(L^16/X^15)` weighted packets; the combinatorial lower bound then still leaves a contradiction.  The exact exceptional weight depends on the final live TP16 normalization, which was not connector-visible here.

## 5. Finite evidence (actual Apéry values, not a countermodel)

An independent exact-integer computation from the Apéry recurrence, for primes `5 <= p <= 500` and all `0 <= h < p`, found **95** defining pairs `p | b_h`.  Grouping them by the direct common value `m=p+h` gives four groups of size at least two, with maximum group size three:

- `m=200`: `(p,h)=(139,61),(181,19)`;
- `m=321`: `(179,142),(193,128),(211,110)`;
- `m=272`: `(191,81),(233,39)`;
- `m=300`: `(191,109),(227,73)`.

For every node in those groups, exact integer arithmetic verifies both

\[
R/p_i-V_i\equiv0\pmod{p_i}
\]

and

\[
BV_i-5z_i-\Gamma_i\equiv0\pmod{p_i}.
\]

The residue rows `(m,p,h,B,z,Gamma,V,BRJ)` are:

```text
(200,139, 61,  17, 29,  13,  42,0)
(200,181, 19, 170,149,  79, 139,0)
(321,179,142,  33, 97, 158,  90,0)
(321,193,128, 187,178,  43, 134,0)
(321,211,110, 118, 44,  17, 154,0)
(272,191, 81,  47,134, 158,  42,0)
(272,233, 39, 195, 13, 133, 191,0)
(300,191,109, 169, 26,  33,  36,0)
(300,227, 73, 186,  4,  94, 191,0)
```

All entries except `m,p,h` are reduced modulo the displayed `p`.  This is finite evidence for the signs and the `R/p -> V` rewrite only; it is not evidence for asymptotic TP16 closure, and it is not a synthetic Apéry counterexample.

The companion stdlib verifier `chatgpt_q32_brj4_audit_verify.py` reconstructs this census from scratch and also checks the formal unit-linear freeness theorem on finite fields.

## 6. Bottom line

1. **BRJ novelty:** no.  It is the `p^2` quotient relation rewritten with `R/p_i == V_i (mod p_i)`.
2. **Signs:** Q8369's `V_i=prod(p_j-p_i)` has no modular sign; direct `g=1` gives `V_i=prod(h_i-h_j)`.  Reflected geometry reverses that node-difference orientation for even packet size 16.
3. **Denominators:** all displayed divisions are integral only after the defining-prime/Lucas divisibilities are established.  Modular inversion is used only for unit coefficients such as `5`, `c_g`, or `V_i` when explicitly solving a congruence.
4. **Face elimination:** scoped no-go proved.  If the face/recurrence ideal does not contain the first quotient `z_i` or `B`, BRJ eliminates `z_i` and imposes no new condition on `B` at primes not dividing `c_g`.
5. **Full fixed-`g` packet:** the quotient BRJ survives on direct and reflected branches with coefficient `c_g=b_g`; `Gamma=10D` does not automatically survive outside direct `g=1`.
6. **Smallest missing theorem:** a p-adic first-quotient alignment such as `[WFQA16]`, plus `0<|B-Q_T|<R`; zero target gives the concrete square-divisibility route `R^2 | b_m`.

That is the exact place where a new theorem would be genuinely beyond the one-`R`/CRT obstruction.
