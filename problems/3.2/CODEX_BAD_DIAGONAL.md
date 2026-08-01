# The bad-diagonal Mellin inverse conjecture: reconstruction and stress test

## Verdict

The archive contains a precise conjectural framework, not a theorem.  The
full ten-item table was recovered from `chatgpt-answers/Q6413_full.md:75-90`;
the three-line `chatgpt-answers/Q6413.md` is only its summary.  The strongest
archive formulation is the (p^\epsilon)-match conjecture for specializations
of two *fixed* compatible systems, after excluding low-order characters and
Kummer/punctual/negligible constituents (`chatgpt-answers/Q6413_full.md:51-57`).
The fixed-constant threshold is explicitly described as strictly stronger and
unsupported (`chatgpt-answers/Q6413_full.md:56`).

The numerical stress test finds a statement-level defect in the motivating
interpretation.  A Kummer twist shifts the Mellin index; it does not normally
produce equality at the *same* character.  Likewise, the Apéry (q(t))-graph
companion produces no exact same-character matches in the tested centered
integral lift, and a literal rank-three Frobenius graph pair produces none
either.  Thus the ten entries are plausible classifications of non-product
joint monodromy, but they are not automatically ten mechanisms producing many
same-index equalities.  This does not disprove the one-way implication
"many matches implies an exceptional relation"; it disproves the stronger
diagnostic expectation that every listed relation should itself yield many
same-index matches.

## 1. Conventions and pinned quantifiers

Let (c\geq 1).  A **fixed system of complexity at most (c)** means a
geometrically semisimple compatible system (\mathcal P), fixed before (p)
varies, whose good specializations are middle-extension objects on
(\mathbf G_m/\mathbf F_p), with

- generic rank (1\leq \operatorname{rank}(\mathcal P)\leq3);
- Artin conductor at most (c), using rank plus local Artin drops;
- coefficient field degree, weight, and number of excluded primes at most
  (c);
- an integral structure at every good (p).

The rank-three Apéry descents and the archive's conductor convention are
recorded at `chatgpt-answers/Q6413_full.md:102-104`.  Bounding the coefficient
field and weight is **[our pinning]**: the source says "complexity \(\leq c\)"
but does not expand every component of that complexity.

For a good prime (p), put (N=p-1), let
(X_p=\operatorname{Hom}(\mathbf F_p^\times,\mu_N)), and define

\[
 t_{i,p}(x)=\operatorname{Tr}(\operatorname{Frob}_x\mid\mathcal P_{i,p}),
 \qquad
 S_{i,p}(\chi)=\sum_{x\in\mathbf F_p^\times}t_{i,p}(x)\chi(x).
\]

The equality (S_{1,p}(\chi)=S_{2,p}(\chi)) is equality in the compositum of
the coefficient field and (\mathbf Q(\mu_N)), not equality after choosing a
complex approximation and not congruence modulo a prime above (p).

For (0<\eta<1), define the high-order character set

\[
 X_p(\eta)=\{\chi\in X_p:\operatorname{ord}(\chi)\geq p^\eta\}.
\]

The source includes this restriction but leaves (\eta) outside the displayed
name `MI(c,eps)` (`chatgpt-answers/Q6413_full.md:55`).  Making (\eta) a fixed
quantified parameter is **[our pinning]**.

### 1.1 Trace-function/object version: `MI(c,epsilon,eta)`

**[source-claimed conjecture, quantifiers pinned]** For every integer (c\geq1)
and every (0<\epsilon,\eta<1), there is
(p_0=p_0(c,\epsilon,\eta)) such that the following holds.  Let
(\mathcal P_1,\mathcal P_2) be two fixed systems of complexity at most (c).
Assume neither semisimplification contains a punctual, negligible, constant,
Tate, or Kummer constituent.  For every good prime (p>p_0), if

\[
 \#\{\chi\in X_p(\eta):S_{1,p}(\chi)=S_{2,p}(\chi)\}\geq p^\epsilon,
 \tag{MI}
\]

then, after semisimplification and finite coefficient extension, the pair is
explained by at least one of the ten geometric mechanisms in Table 1.

This is the source's **cleaned** formulation: items 2--4 have already been
assumed absent.  The equivalent raw formulation drops those exclusions and
allows items 2--4 as conclusions.  Keeping both clauses explicit avoids the
logical redundancy in saying that an excluded Kummer constituent is later an
exception.

This is the closest faithful formalization of
`chatgpt-answers/Q6413_full.md:53-56`.  It is not proved: the missing step is
precisely the conversion of (p^\epsilon) Frobenius incidences on a
trace-equality hypersurface into a tensor invariant
(`chatgpt-answers/Q6413_full.md:61-64`).

### 1.2 Exact-value version and its one rigorous endpoint

For arbitrary functions (f_1,f_2:\mathbf F_p^\times\to K), where (K) has
characteristic zero and contains (\mu_N), write

\[
 \widehat f_i(\chi)=\sum_x f_i(x)\chi(x),\qquad
 B=\{\chi:\widehat f_1(\chi)\ne\widehat f_2(\chi)\}.
\]

Fourier inversion gives the exact identity

\[
 f_1(x)-f_2(x)=\frac1N\sum_{\chi\in B}
   (\widehat f_1(\chi)-\widehat f_2(\chi))\chi(x)^{-1}. \tag{1}
\]

Consequently, if (|B|\leq C), the pointwise difference is a sum of at most
(C) Kummer characters.  For trace functions, semisimplification and
Chebotarev turn this into the Kummer class.  This near-total-agreement endpoint
is the sole rigorous inverse statement claimed in the source
(`chatgpt-answers/Q6413_full.md:57`).  Without fixed-system geometric
provenance, no conclusion stronger than (1) is possible.

### 1.3 Threshold ledger

| Threshold | Status | Exact source reading |
|---|---|---|
| (\geq\delta(p-1)), fixed (\delta>0) | **[source-claimed candidate]** | Called "positive density" in `CRON_FRESH_EYES_pointwise.md:491` and `chatgpt-answers/questions/CRON_qC_inverse_theorem.txt:3`. |
| (\geq p^\epsilon) among characters of order (\geq p^\eta) | **[source-claimed strongest conjecture]** | The displayed `MI(c,eps)` formulation, `chatgpt-answers/Q6413_full.md:55`. |
| (>C(c)) matches | **[our-guess / explicitly stronger]** | The archive says the fixed-constant version is much stronger and unsupported, `chatgpt-answers/Q6413_full.md:56`. |
| equality for all but (C) characters | **[proved algebraic endpoint]** | Fourier inversion and Kummer classification, `chatgpt-answers/Q6413_full.md:57`. |
| mod-(\mathfrak p\) equality | **not covered by MI** | Exact equality is full-Galois stable; mod-(\mathfrak p\) equality is only decomposition-group stable, `chatgpt-answers/Q6413_full.md:29-33`. |

No archive source proves a (p^{1/2}), (p/\log p), or
(p^{2/3-\delta}) exact match bound (`chatgpt-answers/Q6413_full.md:96-98`).

## 2. The reconstructed ten-item exception table

Table 1 reproduces all ten archive items.  The "pair" column makes explicit a
representative object-level shape; it must not be read as a proof that every
such pair has many same-index matches.

| # | Representative pair ( (\mathcal P_1,\mathcal P_2) ) | Mechanism | Why it evades a product-monodromy inverse statement | Archive source |
|---:|---|---|---|---|
| 1 | (C\oplus U,\ C\oplus V) | Common semisimple constituent | The common summand cancels from the difference and puts the joint group on a diagonal factor. | `chatgpt-answers/Q6413_full.md:77` |
| 2 | A generic object paired with a punctual/negligible summand modification | Zero generic Mellin fiber or finite Fourier support | Generic Mellin data cannot see the added object, so object equality cannot be reconstructed without quotienting it out. | `chatgpt-answers/Q6413_full.md:78` |
| 3 | Kummer constituents or finite Kummer packets | Mellin support at finitely many characters | Sparse Fourier support creates exceptional characters without a large connected monodromy relation. | `chatgpt-answers/Q6413_full.md:79` |
| 4 | Constant/Tate line added to either layer | Endpoint or trivial-character correction | It changes only special characters and can fake or remove coincidences at the endpoints. | `chatgpt-answers/Q6413_full.md:80` |
| 5 | (\mathcal P,\mathcal P\otimes\mathcal L_\eta), with a finite self-twist | Periodicity in character space | Tensoring by (\eta) shifts (S(\chi)) to (S(\chi\eta)); a genuine self-twist makes that shift periodic. | `chatgpt-answers/Q6413_full.md:81` |
| 6 | (\mathcal P,\alpha^*\mathcal P), (\alpha(t)=at,a/t), or a composition | Group automorphism | Characters are relabeled, possibly together with inversion; the joint group is an automorphism graph rather than a product. | `chatgpt-answers/Q6413_full.md:82` |
| 7 | (\mathcal P, D(\mathcal P)\otimes\mathcal L_\eta) (or conjugate dual) | Duality/conjugate duality | It forces paired or real traces and can identify (chi) with (chi^{-1}\eta). | `chatgpt-answers/Q6413_full.md:83` |
| 8 | Two objects induced through (t\mapsto t^d) | Power-map induction | The Mellin spectrum is confined to character congruence classes; primitive projectors may vanish for structural reasons. | `chatgpt-answers/Q6413_full.md:84` |
| 9 | Two representations sharing a simple quotient | Graph subgroup of joint monodromy | Goursat gives a proper subgroup projecting onto both factors; trace equality can become an identity on that graph. | `chatgpt-answers/Q6413_full.md:85` |
| 10 | Two arithmetic representations with a special component | Disconnected arithmetic components | Equality can hold identically on one component while remaining a genuine hypersurface on the identity component. | `chatgpt-answers/Q6413_full.md:86` |

The archive designates the Apéry pair
(\mathcal A_- = \mathcal A_+\otimes\mathcal K_q),
(q(t)=t^2-34t+1), as the first nontrivial graph case
(`chatgpt-answers/Q6413_full.md:88-90`).  Here (\mathcal K_q) is the sign
local system of (y^2=q(t)), not a Kummer sheaf of (t).  The finite-field
audit independently checks that (t\mapsto\chi_2(q(t))) is nonmultiplicative
for every requested prime (`CRON_pushforward_check_report.md:185-218`).

This graph designation must also be separated from the withdrawn claim that
the companion is automatically a self-twist.  The archive explicitly retracts
that conflation (`chatgpt-answers/Q6438.md:18-40`).

## 3. Numerical experiment

### 3.1 Exact arithmetic

The executable is `research/scripts/q32_bad_diagonal_stress.py`.  It uses only
the Python standard library.

For (g) the least primitive root and
(t=g^j), an integer trace vector (f) has

\[
 S_f(k)=\sum_{j=0}^{N-1}f(g^j)\zeta_N^{-kj}.
\]

The script reduces the polynomial
(\sum_jf(g^j)X^{-kj}\) modulo the monic cyclotomic polynomial
(\Phi_N(X)).  Two remainder vectors are equal if and only if the algebraic
integers are equal.  There is no floating point evaluation.  Separately, it
maps (\zeta_N\mapsto g\in\mathbf F_p), reduces each Mellin sum modulo (p),
and compares centered representatives.  The two counts are never merged.

### 3.2 Zoo

For each (p\in\{29,37,41,53,61,73,89,101\}), the full core zoo contains
14--30 objects:

- `A*...`: the centered integral lift of the Apéry Hasse--Witt polynomial and
  every Kummer twist of exact order at most six;
- `Aq*...`: its (\chi_2(q(t))) graph companion and the same Kummer twists;
- `Leg2`, `LegS2`: the Legendre ( {}_2F_1(1/2,1/2) ) family and its
  rank-three symmetric square;
- `Fra2`, `FraS2`: the Franel/Beauville-IV elliptic family and its symmetric
  square;
- `FraS2q`: a literal rank-three Frobenius trace function twisted by
  (\chi_2(q(t))), included to test the graph mechanism without the Apéry
  integral-comparison gap;
- `Hes2`: the Hesse ( {}_2F_1(1/3,2/3) ) elliptic family.

There is one necessary limitation.  `A` is an exact integer test vector, but
the archive has not printed an exact characteristic-zero compatible-system
trace at every Apéry middle-extension stalk; it has only the mod-(p)
comparison (`CODEX_FRANEL_MELLIN.md:229-251`).  Therefore the exact `A/Aq`
column is a centered-lift stress test, not a claim about the missing
characteristic-zero Apéry trace.  `FraS2/FraS2q` is the literal rank-three
Frobenius control experiment.

The counterexample hunt adds eight deterministic rank-two elliptic families

\[
 y^2=x^3+(a_0+a_1t)x+(b_0+b_1t)
\]

with coefficient tuples

```text
(3,5,5,1), (-2,4,-1,5), (2,3,-3,-3), (1,0,4,-4),
(-2,0,1,-1), (1,-1,5,4), (-5,2,4,0), (-3,5,-3,-4).
```

Singular fibers are extended by zero.

### 3.3 Representative match matrices

The executable prints the full 14--30 by 14--30 matrix for every prime.  The
following common 8 by 8 submatrices make the cross-prime comparison readable.
Each cell is `exact/mod-p` and the diagonal is (p-1).

#### p = 29

| |A*1|Aq*1|Leg2|LegS2|Fra2|FraS2|FraS2q|Hes2|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|A*1|28/28|0/0|0/1|0/0|0/1|0/1|0/1|0/2|
|Aq*1|0/0|28/28|0/1|0/1|0/0|0/3|0/0|0/1|
|Leg2|0/1|0/1|28/28|0/1|0/1|0/0|0/1|0/0|
|LegS2|0/0|0/1|0/1|28/28|0/0|0/1|0/0|0/3|
|Fra2|0/1|0/0|0/1|0/0|28/28|1/1|0/1|1/1|
|FraS2|0/1|0/3|0/0|0/1|1/1|28/28|0/0|0/0|
|FraS2q|0/1|0/0|0/1|0/0|0/1|0/0|28/28|1/3|
|Hes2|0/2|0/1|0/0|0/3|1/1|0/0|1/3|28/28|

#### p = 37

| |A*1|Aq*1|Leg2|LegS2|Fra2|FraS2|FraS2q|Hes2|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|A*1|36/36|0/0|0/1|0/0|0/1|0/0|0/0|0/2|
|Aq*1|0/0|36/36|0/0|0/0|0/0|0/2|0/1|0/1|
|Leg2|0/1|0/0|36/36|0/0|0/0|0/2|0/0|0/1|
|LegS2|0/0|0/0|0/0|36/36|0/0|0/0|0/0|0/0|
|Fra2|0/1|0/0|0/0|0/0|36/36|3/3|0/1|0/0|
|FraS2|0/0|0/2|0/2|0/0|3/3|36/36|0/1|0/4|
|FraS2q|0/0|0/1|0/0|0/0|0/1|0/1|36/36|0/0|
|Hes2|0/2|0/1|0/1|0/0|0/0|0/4|0/0|36/36|

#### p = 41

| |A*1|Aq*1|Leg2|LegS2|Fra2|FraS2|FraS2q|Hes2|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|A*1|40/40|0/0|0/2|0/2|0/3|0/0|0/0|0/0|
|Aq*1|0/0|40/40|0/1|0/3|0/0|0/0|0/0|0/2|
|Leg2|0/2|0/1|40/40|0/0|0/2|0/1|0/2|0/0|
|LegS2|0/2|0/3|0/0|40/40|0/1|0/1|0/1|0/0|
|Fra2|0/3|0/0|0/2|0/1|40/40|0/0|0/0|1/1|
|FraS2|0/0|0/0|0/1|0/1|0/0|40/40|0/0|0/0|
|FraS2q|0/0|0/0|0/2|0/1|0/0|0/0|40/40|0/0|
|Hes2|0/0|0/2|0/0|0/0|1/1|0/0|0/0|40/40|

#### p = 53

| |A*1|Aq*1|Leg2|LegS2|Fra2|FraS2|FraS2q|Hes2|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|A*1|52/52|0/2|0/2|0/0|0/3|0/3|0/1|0/0|
|Aq*1|0/2|52/52|0/1|0/0|0/0|0/0|0/0|1/4|
|Leg2|0/2|0/1|52/52|0/1|0/0|0/0|0/1|0/0|
|LegS2|0/0|0/0|0/1|52/52|0/0|0/2|0/0|0/0|
|Fra2|0/3|0/0|0/0|0/0|52/52|1/3|0/2|1/2|
|FraS2|0/3|0/0|0/0|0/2|1/3|52/52|0/1|0/0|
|FraS2q|0/1|0/0|0/1|0/0|0/2|0/1|52/52|0/2|
|Hes2|0/0|1/4|0/0|0/0|1/2|0/0|0/2|52/52|

#### p = 61

| |A*1|Aq*1|Leg2|LegS2|Fra2|FraS2|FraS2q|Hes2|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|A*1|60/60|0/2|0/2|0/2|0/0|0/1|0/1|0/3|
|Aq*1|0/2|60/60|0/0|0/0|0/1|0/1|0/1|0/2|
|Leg2|0/2|0/0|60/60|0/1|0/0|0/0|0/0|0/0|
|LegS2|0/2|0/0|0/1|60/60|0/2|0/1|0/2|1/1|
|Fra2|0/0|0/1|0/0|0/2|60/60|1/1|0/2|2/2|
|FraS2|0/1|0/1|0/0|0/1|1/1|60/60|0/0|0/0|
|FraS2q|0/1|0/1|0/0|0/2|0/2|0/0|60/60|0/0|
|Hes2|0/3|0/2|0/0|1/1|2/2|0/0|0/0|60/60|

#### p = 73

| |A*1|Aq*1|Leg2|LegS2|Fra2|FraS2|FraS2q|Hes2|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|A*1|72/72|0/0|0/2|0/0|0/0|0/2|0/0|0/2|
|Aq*1|0/0|72/72|0/0|0/0|0/2|0/0|0/0|0/2|
|Leg2|0/2|0/0|72/72|0/0|0/0|0/2|0/0|0/2|
|LegS2|0/0|0/0|0/0|72/72|0/0|0/1|0/0|0/0|
|Fra2|0/0|0/2|0/0|0/0|72/72|0/0|0/3|0/1|
|FraS2|0/2|0/0|0/2|0/1|0/0|72/72|0/4|0/1|
|FraS2q|0/0|0/0|0/0|0/0|0/3|0/4|72/72|0/2|
|Hes2|0/2|0/2|0/2|0/0|0/1|0/1|0/2|72/72|

#### p = 89

| |A*1|Aq*1|Leg2|LegS2|Fra2|FraS2|FraS2q|Hes2|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|A*1|88/88|0/2|0/4|0/0|0/3|0/1|0/1|0/2|
|Aq*1|0/2|88/88|0/2|0/0|0/1|0/1|0/0|0/0|
|Leg2|0/4|0/2|88/88|0/1|0/0|0/2|0/1|0/3|
|LegS2|0/0|0/0|0/1|88/88|0/1|0/0|0/1|0/0|
|Fra2|0/3|0/1|0/0|0/1|88/88|0/0|0/0|1/1|
|FraS2|0/1|0/1|0/2|0/0|0/0|88/88|0/1|0/1|
|FraS2q|0/1|0/0|0/1|0/1|0/0|0/1|88/88|0/2|
|Hes2|0/2|0/0|0/3|0/0|1/1|0/1|0/2|88/88|

#### p = 101

| |A*1|Aq*1|Leg2|LegS2|Fra2|FraS2|FraS2q|Hes2|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|A*1|100/100|0/0|0/3|0/2|0/1|0/0|0/0|0/0|
|Aq*1|0/0|100/100|0/2|0/2|0/0|0/0|0/2|0/0|
|Leg2|0/3|0/2|100/100|0/0|0/0|0/1|0/0|0/0|
|LegS2|0/2|0/2|0/0|100/100|0/0|0/3|0/0|0/0|
|Fra2|0/1|0/0|0/0|0/0|100/100|1/1|0/1|1/2|
|FraS2|0/0|0/0|0/1|0/3|1/1|100/100|0/2|0/0|
|FraS2q|0/0|0/2|0/0|0/0|0/1|0/2|100/100|0/1|
|Hes2|0/0|0/0|0/0|0/0|1/2|0/0|0/1|100/100|

### 3.4 Kummer and graph verdicts

For the untwisted Apéry lift and its companion, the counts over the eight
primes are

| comparison | 29 | 37 | 41 | 53 | 61 | 73 | 89 | 101 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| exact | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| mod (p) | 0 | 0 | 0 | 2 | 2 | 0 | 2 | 0 |

For the literal rank-three `FraS2/FraS2q` graph pair:

| comparison | 29 | 37 | 41 | 53 | 61 | 73 | 89 | 101 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| exact | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| mod (p) | 0 | 1 | 0 | 1 | 0 | 4 | 1 | 2 |

Therefore the graph case behaves like an unrelated pair under the literal
same-(\chi) match statistic.  It is genuinely graph-type at the sheaf and
joint-monodromy level, but it does **not** break the naive same-character
inverse theorem by producing a positive proportion of matches.

Distinct Kummer twists within one Apéry layer have exact match ranges

```text
(p,min,max) = (29,0,2), (37,0,6), (41,2,2), (53,0,2),
              (61,0,6), (73,2,6), (89,2,2), (101,0,2).
```

The identity is (S_{f\otimes\eta}(\chi)=S_f(\chi\eta)), not
(S_{f\otimes\eta}(\chi)=S_f(\chi)).  Hence the requested expectation that
all twist-related pairs give (\asymp p) same-index matches is false unless a
separate self-twist/periodicity hypothesis is imposed.

## 4. Counterexample hunt

The stable cross-prime search set consists of 18 objects: the untwisted and
quadratic Kummer twists of each Apéry layer, the eight named hypergeometric or
symmetric-square objects, and eight deterministic random rank-two elliptic
families.  After excluding same-family Kummer pairs, the two graph pairs, and
obvious Legendre/Franel symmetric-power relations, 143 unordered pairs remain.

The pre-registered violation criterion was

\[
 \#\{\chi:S_1(\chi)=S_2(\chi)\}\geq
 \left\lceil\frac{p-1}{8}\right\rceil
 \quad\text{for each }p=73,89,101.
\]

No exact pair passes it.  In fact, among the 143 unrelated pairs the largest
exact count at any tested prime is 2; the leading sequences include

```text
Fra2 vs Hes2: [1,0,1,1,2,0,1,1]
Fra2 vs Rnd3: [0,1,0,2,0,0,0,1]
Fra2 vs Rnd4: [0,1,2,0,0,0,0,0]
Hes2 vs Rnd2: [0,0,0,0,2,0,0,0]
Rnd0 vs Rnd1: [2,0,0,0,0,0,0,0].
```

Thus the exact counterexample hunt is negative on the stated finite search
space.  It is evidence for an (O(1)) unrelated-match regime, not a proof of
uniform boundedness.

The mod-(p) column is qualitatively different.  For example, `Rnd6/Rnd7`
has mod-(p) counts `[2,0,1,0,0,54,1,0]`, and `Rnd0/Rnd7` has a count 30 at
(p=41), while their exact counts remain at most one.  These isolated large
residual coincidences arise because reductions of elliptic trace functions
can have truncated Hasse--Witt Mellin support.  They do not persist across the
last three primes, but they are a concrete warning that exact and
defining-characteristic equality cannot share one inverse theorem.  This is
the numerical form of the archive's Galois-orbit divorce
(`CODEX_LT_MELLIN.md:13-29` and
`chatgpt-answers/Q6413_full.md:29-33`).

## 5. What a proof would need

1. **Correct statistic and symmetry quotient.**  Decide whether the theorem
   counts literal same-character equality or equality after the relabeling
   naturally induced by Kummer twists, automorphisms, inversion, and duality.
   The numerical test shows these are different statements.

2. **Multiplicative-translation independence.**  After quotienting common,
   Kummer, punctual, and power-induced pieces, prove that distinct translates
   have product geometric monodromy.  Goursat--Kolchin--Ribet should reduce a
   proper joint group to the automorphism/duality/power-map/graph cases listed
   in Table 1; the archive assigns this to the joint Tannakian package
   (`chatgpt-answers/Q6413_full.md:59-64` and
   `CRON_FRESH_EYES_pointwise.md:519`).

3. **A primitive-projector theorem.**  Exact equality is stable under the full
   cyclotomic Galois orbit.  One must expand the corresponding primitive
   projector using Ramanujan sums, then prove that it cannot annihilate a
   bounded-conductor noninduced object except through Kummer packets or
   power-map induction (`chatgpt-answers/Q6413_full.md:64`).

4. **An incidence-to-invariant theorem.**  Even with the joint group known,
   there is currently no result turning only (p^\epsilon) Frobenius points on
   the trace-equality hypersurface into a tensor invariant.  This is the exact
   missing theorem, not a consequence of generic Mellin equidistribution
   (`chatgpt-answers/Q6413_full.md:61-63`).

5. **Bilinear ceiling.**  Pairwise Deligne cancellation calibrates only to
   (p^{3/4}); a shift amplifier by itself does not cross the existing
   (2/3) threshold (`CRON_FRESH_EYES_pointwise.md:519` and
   `chatgpt-answers/Q6413_full.md:66-73`).  Higher intersections, a spectral
   large sieve, or another bounded-complexity detector would be necessary.

6. **Separate residual theorem.**  Even a complete exact inverse theorem says
   nothing about nonzero algebraic integers that vanish modulo the selected
   prime above (p).  The archive isolates this as (DRS)/(RLL), with
   (\alpha<1/2) needed at the (k=2) interface
   (`chatgpt-answers/Q6413_full.md:43-49`).  Existing exact orbit norms fail
   because their Archimedean scale is (p^{3/2}), and even orbitwise AM--GM is
   nontrivial only below RMS (p) (`CODEX_LT_MELLIN.md:158-225`).

The correct final status is therefore: a well-posed exact-side research
program after the quantifier and statistic corrections above; no proof at the
(p^\epsilon) threshold; no numerical exact counterexample in this zoo; and a
clear numerical refutation of treating either Kummer or Apéry graph relations
as automatic sources of many literal same-character matches.

---

## Adjudication addendum (life session, 2026-08-01 09:40)

Post-delivery cross-check against cron's Q6445 verdict (appendix AF, commit
1d9ecee), which KILLED the original Q6339 route "positive density zeros =>
O(1/delta) cosets => bounded self-twist" by an explicit dilation-family
counterexample (coset covering lower bound >= q, exhaustively verified):

1. The main formalization above, MI(c,epsilon,eta), is the Q6413 sec.II
   shape and does NOT rely on the killed coset-covering step. It stands as
   stated (still unproved; the incidence-to-invariant step remains the gap).
2. Threshold-ledger row "positive density / [source-claimed candidate]" is
   re-graded: the positive-density hypothesis is now known NOT to imply
   bounded self-twist by the naive projector-covering argument. Its salvage
   is cron's C0 unconditional theorem (zero set = union of complete
   order-packets <=> Phi_d | A <=> projector = 0; all-but-K => difference =
   <= K Ramanujan waves of order <= 2K^2) plus Prop 3 conditional theorem
   (equidistribution + positive-density coincidence => joint-group component
   trace identity => Goursat/Laurent/Mann). See cron appendix AF.
3. The empirical matrices (sec. 3) are data and unaffected.
4. Route consequence (cron concurring): the C-line's next step is B's three
   inputs (joint group / component classification / horizontal
   equidistribution), not another cyclotomic inverse theorem.
