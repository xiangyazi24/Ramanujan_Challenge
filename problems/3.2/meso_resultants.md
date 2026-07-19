# Separated-block resultants at height 24

This is the computational record for CODEX_SPEC_meso.md. Every polynomial
and resultant below is computed over the integers. The calculations are
reconnaissance and exact verification of displayed identities; none of the
finite data is used as proof of a uniform estimate.

Run

    python3 problems/3.2/meso_explore.py
    python3 problems/3.2/meso_verify.py

Both wrappers automatically re-execute under SageMath. The first command
prints the certified partial factorisation of all 121 symmetric
representatives and the prime statistics. The second independently
recomputes the data and checks the symbolic identities used in
meso_result.tex.

## Four different support counts

For $d,r\geq2$, put

\[
X_{d,r}(p)=\{x\in\mathbf F_p:N_d(x)=N_r(x+d)=0\},
\qquad m_{d,r}(p)=|X_{d,r}(p)|.
\]

The computation keeps the following quantities separate:

- $V_p(H)$: ordered pairs with $p\mid\mathcal R_{d,r}$, where
  $\mathcal R_{d,r}=\operatorname{Res}(N_d(X),N_r(X+d))$;
- $A_p(H)$: ordered pairs whose reduced polynomials have a common affine
  root over $\overline{\mathbf F}_p$;
- $W_p(H)$: ordered pairs with an actual common root in $\mathbf F_p$;
- $\mathcal E_p(H)=\sum m_{d,r}(p)$: the root-weighted collision energy.

Always

\[
\mathcal E_p(H)\geq W_p(H),\qquad
W_p(H)\leq A_p(H)\leq V_p(H),
\]

but equality need not hold at either step. In particular, the sentence in
the original specification identifying $p\mid\mathcal R_{d,r}$ with an
$\mathbf F_p$ root is false.

## Exact resultant census

For all $2\leq d,r$ with $d+r\leq24$:

- ordered pairs: 231;
- representatives with $d\leq r$: 121;
- all 121 resultants are positive and nonzero;
- largest: $\mathcal R_{12,12}$, with 2246 decimal digits;
- SHA-256 of the newline-separated records d,r,R in lexicographic order:

    c50c84589bc8426e978b8299e82aa415f4fb047de283083fe2549810f52d672f

The exact integers are deterministically regenerated before the digest is
computed. The default factor table removes every prime at most 5000 and
records a fingerprint of the remaining exact cofactor: its decimal length
and the first 16 hexadecimal digits of its SHA-256. It never labels an
unresolved cofactor as prime.

### Certified complete factorisations for $d+r\leq8$

Every displayed factor was individually primality-proved.

| $(d,r)$ | complete factorisation of $\mathcal R_{d,r}$ |
|---|---|
| $(2,2)$ | $2^5 5^4 11^2 17^3 71$ |
| $(2,3)$ | $5^3\cdot23\cdot103\cdot129763\cdot685784151628061$ |
| $(2,4)$ | $2^5 3^6 5^7 11\cdot13\cdot17^3 41\cdot653\cdot52385933\cdot6024613671641$ |
| $(3,3)$ | $3^{17}5^5 7^8 11^3 89153\cdot110629^2\cdot315735829^2$ |
| $(2,5)$ | $3^6 5^8 43\cdot113\cdot397\cdot11119\cdot30757\cdot34283\cdot233591\cdot21869834584034319215521$ |
| $(3,4)$ | $2^{14}5^5 13\cdot73\cdot131\cdot57131\cdot75577\cdot15260305013831376290046094059172675684932975108889335505111957$ |
| $(2,6)$ | $2^8 3^6 5^7 17^3 41\cdot1091\cdot3919541\cdot13365198871\cdot5998798267201\cdot511088437530213575388589$ |
| $(3,5)$ | $2^{18}5^7 59\cdot129763\cdot5784767\cdot98608507\cdot521339821\cdot170262135457\cdot111600466182360778561621\cdot22456093031506928723038485585657593428847$ |
| $(4,4)$ | $2^{42}5^{11}13\cdot17^3 79\cdot101^2 577^3 4787^2 8431^2 56003\cdot271030933^2\cdot45032941417^2\cdot16802972635975249^2$ |

Full factorisation ceases to be routine already at $d+r=9$. After removing
the factors at most 5000 from $\mathcal R_{3,6}$ and then the certified
prime factors 21187 and 3605629, the following 104-digit composite remains:

    82331328057159957094082641924652164496892111499574525894473111366919585965389051882753610906780816656199

Thus all 121 resultants were computed exactly and bounded-factorised, but a
claim of 121 complete prime factorisations would be false.

### Repetition and known sequences

Among prime factors at most 5000 there are 132 distinct primes. Fourteen
also divide an Apéry value $b_j$, $1\leq j\leq22$, and 37 also divide a
Pell coefficient $\ell_j$, $1\leq j\leq22$; 90 occur in neither list.
The most frequent factors among the 231 ordered pairs are:

| prime | number of pairs | maximum valuation |
|---:|---:|---:|
| 5 | 231 | 146 |
| 2 | 210 | 618 |
| 3 | 203 | 260 |
| 11 | 151 | 38 |
| 7 | 140 | 30 |
| 17 | 120 | 8 |
| 19 | 84 | 8 |
| 13 | 79 | 9 |

Fully factoring the $b_j$ and $\ell_j$ through $j=22$ gives 109 known
prime supports. Removing those primes accounts for only 17.4972% of the
aggregate logarithmic height of the 121 representatives, and 120 of the
121 remaining cofactors are still composite. The overlap is therefore not
explained by a product formula using the tested values $b_j,\ell_j$ with
$j\leq22$. The reproducible structural factors found here are instead the
half-integer center factors and diagonal squares proved in meso_result.tex.

## Prime statistics through 5000

There are 564 primes in the mesoscopic test range
$24^2\leq p\leq5000$. The distributions below list value:number of primes.

| quantity | distribution |
|---|---|
| $V_p(24)$ | $0:506,\ 1:3,\ 2:38,\ 3:2,\ 4:10,\ 6:3,\ 15:1,\ 16:1$ |
| $A_p(24)$ | $0:508,\ 1:3,\ 2:38,\ 3:2,\ 4:10,\ 6:2,\ 16:1$ |
| $W_p(24)$ | $0:509,\ 1:2,\ 2:39,\ 3:1,\ 4:10,\ 6:2,\ 16:1$ |
| $\mathcal E_p(24)$ | $0:509,\ 1:2,\ 2:39,\ 4:11,\ 6:2,\ 16:1$ |

The maximum of both $V_p(24)$ and $\mathcal E_p(24)$ is 16 at $p=653$:

\[
\frac{V_{653}(24)}{24}=\frac23,\qquad
\frac{\mathcal E_{653}(24)}{24^{3/2}}=0.1360827\ldots.
\]

All 16 supported pairs have one simple $\mathbf F_{653}$ common root:

| $(d,r)$ | $x$ | $(d,r)$ | $x$ |
|---|---:|---|---:|
| $(2,4)$ | 313 | $(2,18)$ | 313 |
| $(2,22)$ | 313 | $(4,2)$ | 333 |
| $(4,14)$ | 315 | $(4,18)$ | 315 |
| $(4,20)$ | 315 | $(6,14)$ | 313 |
| $(6,18)$ | 313 | $(14,4)$ | 319 |
| $(14,6)$ | 319 | $(18,2)$ | 319 |
| $(18,4)$ | 315 | $(18,6)$ | 315 |
| $(20,4)$ | 313 | $(22,2)$ | 315 |

The distinctions between the four counts already occur in this range:

- $p=577$: $V=15$, but $A=W=\mathcal E=0$. The indices satisfy
  $d,r\in\{4,8,12,16,20\}$; both leading coefficients vanish and the
  nominal degrees drop by exactly two.
- $p=1153$: the same phenomenon occurs for
  $d,r\in\{6,12,18\}$, giving $V=6$ and no affine edge.
- $p=4787$: $V=A=1$, but $W=\mathcal E=0$; for $(4,4)$ the affine gcd
  is an irreducible quadratic.
- $p=797$: there are three algebraic edges but only two split over
  $\mathbf F_p$.
- $p=3109$: there are three $\mathbf F_p$ edges and energy four, because
  the $(8,8)$ gcd has two distinct roots.

For all 669 primes $p\leq5000$, including the singular
small-characteristic range, the largest raw support is
$V_5(24)=231$, and the largest energy is $\mathcal E_5(24)=747$.
These values are not relevant to $H\leq\sqrt p$.

## Center factors seen in the table

For even $a$, define

\[
T_b^{(a)}=2^{3(b-1)}N_b\left(\frac{a-1}{2}\right)\in\mathbf Z.
\]

These values satisfy an integer recurrence and their odd parts divide both
$\mathcal R_{a,b}$ and $\mathcal R_{b,a+b}$. This explains large repeated
factors invisible in a cutoff-5000 table. For example, $T_{10}^{(4)}$
contains the 35-digit prime

    30681778082168266499711406058345639

which divides both $\mathcal R_{4,10}$ and $\mathcal R_{10,14}$.
The symmetric propagation applies when the second block is even.

## STALL REPORT

The computation supports sparse affine collisions, but it does not prove
G1 or G2. More importantly, the proposed implication chain has two logical
gaps.

First, $V_p(H)$ is only a simple support count, whereas

\[
\mathcal E_p(H)=\sum_{d+r\leq H}m_{d,r}(p)
\]

is weighted. The existing bound
$m_{d,r}\leq3(\min(d,r)-1)$ gives only

\[
V_p(H)\ll H\quad\Longrightarrow\quad\mathcal E_p(H)\ll H^2.
\]

The correct sufficient arithmetic target is the degree-weighted support

\[
\sum_{d+r\leq H}(\min(d,r)-1)
\mathbf1_{p\mid\mathcal R_{d,r}}\ll H^{3/2},
\]

or an equivalent average-gcd-multiplicity or tail estimate.
Kővári--Sós--Turán on the simple support graph cannot control parallel root
witnesses.

Second, even $\mathcal E_p(H)\ll H^{3/2}$ does not control the singleton
columns of the gap-root incidence graph, so it does not by itself imply
$R_p(H)\ll H$. A sufficient missing amplification statement is

\[
R_p(H)\ll H+H^{-1/2}\mathcal E_p(H),
\]

together with a separate bound for cut-edge columns. Neither statement is
provided by the current resultant identities.

Finally, raw $V_p(H)$ is contaminated by a leading-coefficient lattice. If
$\rho_p$ is the first positive index with $p\mid\ell_{\rho_p}$, then all
pairs with $\rho_p\mid d,r$ are automatic projective-infinity edges. At
$p=665857$, $\rho_p=8$ and $H=\lfloor\sqrt p\rfloor=816$, so the raw
support already contains

\[
\frac{102\cdot101}{2}=5151
\]

such pairs, irrespective of affine collisions. Any future counting theorem
should use affine collision support or first saturate the common point at
infinity.

The rigorous output of this attack is therefore G3: root strips, centered
norms, center-factor propagation, and the diagonal square law. The G4
obstruction is also precise: the missing inputs are a weighted affine
resultant estimate and a low-fiber collision-amplification lemma.
