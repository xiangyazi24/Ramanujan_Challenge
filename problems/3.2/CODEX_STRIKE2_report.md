NO UNBOUNDED L PROVED.

# Second deep strike: the (4.8)/(4.9) escape hatches

## 1. Terminal verdict

Let

\[
N=p-2,\qquad D=\sqrt N\,L,
\]

with all collisions nonwrapping.  Neither

\[
P_DQ_D\ll N^2                                                   \tag{1.1}
\]

nor

\[
P_D\ll \frac{N}{L^2\log D}                                    \tag{1.2}
\]

was proved for any unbounded $L=L(p)$.  Consequently this strike does not
improve the $3/2$ energy exponent.

The strike did, however, change the exact diagnosis of this escape hatch.

1. The first-return, cascade, and inclusion--exclusion identities can all be
   made exact, but they do not force a useful union saving.  In the requested
   data, the top dyadic shell saves only $2.3\%--8.3\%$, and most lower
   shells have no repeated base at all.
2. A detector using the first two factorial moments gives no new upper bound
   for $P_D$.  With the restart cap $k_r\le M_D$, the optimal quadratic
   detector is
   
   \[
   P_D\le S_D-\frac{2Q_D}{M_D},                                \tag{1.3}
   \]
   
   which still contains the unknown first moment and does not close.
3. From the stated unconditional inputs alone, the best product bound at the
   wall remains
   
   \[
   P_DQ_D\le 66ND^2(1+\log D)
             =66N^2L^2(1+\log D).                              \tag{1.4}
   \]
4. There is a correction to the claimed calibration of the reflection word
   in `CODEX_LASTWALL_report.md` Section 4.4.  For that literal word,
   
   \[
   P_H=q^3+2q^2-q-2,\qquad Q_H=2\binom{q^2}{3},                 \tag{1.5}
   \]
   
   so $P_HQ_H\asymp q^9$, not $q^8$.  It still satisfies
   $P_HQ_H/N_q^2\asymp q^{-1}\to0$, so the original word does not refute
   (1.1).
5. More decisively, the same triangle-rich word can be padded by
   primitive-only reflected pairs.  The padded words retain reflection, no
   adjacent equality, $C_d\le3(d-1)$, the full window bound
   $m_J(v)-1\le4\operatorname{span}(J)^{2/3}$, the exact $h=2$ root-count
   law, and every word-level renewal identity, but satisfy
   
   \[
   \frac{D_q}{\sqrt{N_q}}\longrightarrow\infty,qquad
   P_{D_q}\asymp N_q,qquad Q_{D_q}\asymp q^6,qquad
   \frac{P_{D_q}Q_{D_q}}{N_q^2}\asymp q\longrightarrow\infty. \tag{1.6}
   \]
   
   This is not an Apéry counterexample.  It proves that (1.1) does not follow
   from the currently listed word/clique inputs.  A successful proof must use
   arithmetic that rules out coexistence of a triangle-rich core with a large
   $k_r=1$ primitive stratum.

The persistent verifier is `CODEX_STRIKE2_verify.py`.

## 2. V4 numerical calibration -- done first

### 2.1 Conventions and independent construction gate

The two integer cutoffs are

\[
D=\lceil\sqrt p\log p\rceil,qquad D=\lceil p^{3/5}\rceil.
\]

The natural logarithm is used.  The $p^{3/5}$ ceiling is selected by exact
integer fifth-power comparisons; the first ceiling uses 80-digit Decimal
arithmetic and is checked to lie more than $10^{-50}$ from an integer.

For each prime, the verifier constructs the collision set in two independent
ways:

1. equality of the normalized projective Apéry orbit values;
2. simultaneous propagation of every $N_d(r)$ by
   
   \[
   N_{d+1}(r)=P(r+d)N_d(r)-(r+d)^6N_{d-1}(r).
   \]

All $8,782,308$ admissible orbit/continuant comparisons agree.

### 2.2 Main exact table

Here **SL** denotes $\lceil\sqrt p\log p\rceil$, and **3/5** denotes
$\lceil p^{3/5}\rceil$.  The fraction in the penultimate column is exact.
The last column is the relative global union saving $(S_D-P_D)/S_D$.

| $p$ | scale | $D$ | $S_D$ | $P_D$ | $Q_D$ | $P_DQ_D/N^2$ | union saving |
|---:|:---:|---:|---:|---:|---:|:---|---:|
| 997  | SL  | 219 | 313  | 266  | 58  | $15428/990025=0.015583444863$ | 15.015974% |
| 997  | 3/5 | 63  | 87   | 87   | 0   | $0$ | 0% |
| 1999 | SL  | 340 | 498  | 431  | 84  | $36204/3988009=0.009078214217$ | 13.453815% |
| 1999 | 3/5 | 96  | 144  | 137  | 8   | $1096/3988009=0.000274823853$ | 4.861111% |
| 4001 | SL  | 525 | 754  | 695  | 78  | $18070/5330667=0.003389819698$ | 7.824934% |
| 4001 | 3/5 | 145 | 198  | 196  | 2   | $392/15992001=0.000024512255$ | 1.010101% |
| 7919 | SL  | 799 | 1229 | 1114 | 142 | $158188/62678889=0.002523784364$ | 9.357201% |
| 7919 | 3/5 | 219 | 349  | 329  | 32  | $1504/8954127=0.000167967240$ | 5.730659% |

The ratios $P_D/D$ in the same order are

\[
1.2146, 1.3810, 1.2676, 1.4271, 1.3238, 1.3517, 1.3942, 1.5023.
\]

Thus the primitive support is empirically linear in $D$, but nearly every
observed collision is already primitive.  The finite-size diagnostic

\[
\Xi_D:=\frac{P_DD^2\log D}{N^2}
       =\frac{P_D}{N/(L^2\log D)}
\]

is $69.44,72.82,75.03,75.83$ at the **SL** cutoffs and
$1.445,1.445,1.282,1.357$ at the **3/5** cutoffs.  This does not contradict
the asymptotic heuristic $P_D\asymp D$, but it shows that these primes are
not numerically close to (1.2) with constant one.

### 2.3 Exact first-return spacing distributions

For a base with a return, let $\tau(r)$ be its first-return lag.  The next
table gives exact counts in the ten equal-width bins

\[
\left(\left\lfloor\frac{jD}{10}\right\rfloor,
      \left\lfloor\frac{(j+1)D}{10}\right\rfloor\right],
\qquad 0\le j<10.
\]

The entries in each row sum exactly to $P_D$.

| $p$ | $D$ | exact bin counts for $\tau(r)$ |
|---:|---:|:---|
| 997  | 219 | 26, 29, 33, 29, 25, 23, 21, 33, 22, 25 |
| 997  | 63  | 3, 9, 11, 9, 9, 7, 10, 5, 9, 15 |
| 1999 | 340 | 63, 42, 39, 43, 41, 32, 61, 40, 36, 34 |
| 1999 | 96  | 10, 29, 13, 13, 10, 10, 19, 11, 7, 15 |
| 4001 | 525 | 74, 68, 76, 83, 71, 76, 70, 51, 57, 69 |
| 4001 | 145 | 21, 23, 15, 24, 13, 15, 23, 18, 25, 19 |
| 7919 | 799 | 128, 111, 112, 111, 112, 118, 99, 101, 115, 107 |
| 7919 | 219 | 46, 32, 33, 35, 19, 39, 26, 32, 39, 28 |

The exact nearest-rank deciles $q_0,q_{0.1},\ldots,q_1$ are:

| $p$ | $D$ | first-return deciles |
|---:|---:|:---|
| 997  | 219 | 2, 22, 43, 62, 78, 102, 125, 154, 170, 194, 218 |
| 997  | 63  | 2, 10, 18, 22, 26, 33, 42, 46, 56, 62, 63 |
| 1999 | 340 | 2, 24, 58, 90, 117, 163, 204, 225, 254, 301, 340 |
| 1999 | 96  | 2, 10, 16, 21, 30, 44, 56, 65, 74, 90, 96 |
| 4001 | 525 | 2, 50, 104, 153, 197, 244, 293, 344, 407, 472, 524 |
| 4001 | 145 | 2, 12, 25, 42, 57, 75, 91, 104, 120, 130, 144 |
| 7919 | 799 | 2, 66, 148, 226, 308, 389, 460, 551, 639, 717, 799 |
| 7919 | 219 | 2, 16, 34, 63, 82, 108, 128, 154, 176, 193, 218 |

The lag distributions are broad and approximately flat.  There is no
visible concentration at a small collection of restart spacings which could
be removed separately.

### 2.4 Full dyadic shell union ledger

For a shell $J$, write

\[
S_J=\sum_r k_J(r),\qquad
U_J=\#\{r:k_J(r)>0\},\qquad
F_J=\#\{r:\tau(r)\in J\}.
\]

Each item below is `inclusive shell: S_J/U_J/F_J`.  Every number is exact.

```text
p=997, D=219:
110-219:161/149/124, 55-109:85/80/75, 28-54:30/30/30,
14-27:25/25/25, 7-13:9/9/9, 4-6:2/2/2, 2-3:1/1/1, 1:0/0/0.
sum U_J=296, P_D=266, cross-shell repetition=30.

p=997, D=63:
32-63:46/46/46, 16-31:28/28/28, 8-15:8/8/8,
4-7:4/4/4, 2-3:1/1/1, 1:0/0/0.
sum U_J=87, P_D=87, cross-shell repetition=0.

p=1999, D=340:
171-340:253/232/203, 86-170:121/115/107, 43-85:53/53/53,
22-42:29/28/26, 11-21:27/27/27, 6-10:9/9/9,
3-5:5/5/5, 2:1/1/1, 1:0/0/0.
sum U_J=470, P_D=431, cross-shell repetition=39.

p=1999, D=96:
49-96:66/64/62, 25-48:32/32/31, 13-24:30/29/28,
7-12:9/9/9, 4-6:6/6/6, 2-3:1/1/1, 1:0/0/0.
sum U_J=141, P_D=137, cross-shell repetition=4.

p=4001, D=525:
263-525:371/350/323, 132-262:204/202/195, 66-131:89/89/89,
33-65:42/41/40, 17-32:22/22/22, 9-16:10/10/10,
5-8:10/10/10, 3-4:3/3/3, 2:3/3/3, 1:0/0/0.
sum U_J=730, P_D=695, cross-shell repetition=35.

p=4001, D=145:
73-145:100/100/100, 37-72:44/43/42, 19-36:27/27/27,
10-18:9/9/9, 5-9:12/12/12, 3-4:3/3/3,
2:3/3/3, 1:0/0/0.
sum U_J=197, P_D=196, cross-shell repetition=1.

p=7919, D=799:
400-799:616/590/540, 200-399:290/283/271, 100-199:158/154/147,
50-99:77/74/71, 25-49:36/36/35, 13-24:28/28/26,
7-12:9/9/9, 4-6:8/8/8, 2-3:7/7/7, 1:0/0/0.
sum U_J=1189, P_D=1114, cross-shell repetition=75.

p=7919, D=219:
110-219:175/171/164, 55-109:77/74/71, 28-54:44/44/43,
14-27:27/27/25, 7-13:11/11/11, 4-6:8/8/8,
2-3:7/7/7, 1:0/0/0.
sum U_J=342, P_D=329, cross-shell repetition=13.
```

In the four **SL** top shells, $1-U_J/S_J$ is respectively
$7.45\%,8.30\%,5.66\%,4.22\%$.  At the **3/5** cutoffs it is
$0\%,3.03\%,0\%,2.29\%$.  The numerical instruction to aim V1 therefore
has a sharp answer: the hoped-for factor is not coming from within-shell
cascade overlap.  The observed room is in the unexpectedly small first
moment $S_J\asymp |J|$, not in removing multiplicity from a degree-sized
sum.

## 3. Exact audit of the Section 4.4 reflection word

Let $p_q$ be the first prime at least $q^5$, let $N_q=p_q-2$, and put

\[
a_n=2qn+((n\bmod q)^2\bmod q),\qquad 0\le n<q^2,
\]

with $H=2q^3-2q+1$.  The two special clusters each contain $q^2$
positions and have span exactly $H$.  Cross-cluster gaps exceed $H$.

Each cluster contributes $q^2-1$ active bases.  The remaining private
reflection pairs contribute one active base for every even gap
$2\le d\le H$, hence

\[
\frac{H-1}{2}=q^3-q
\]

further bases.  Therefore

\[
P_H=2(q^2-1)+(q^3-q)=q^3+2q^2-q-2.                \tag{3.1}
\]

Only the two special clusters make triangles, so

\[
Q_H=2\binom{q^2}{3}.                               \tag{3.2}
\]

The exact edge count is

\[
S_H=2\binom{q^2}{2}+q^3-q.                         \tag{3.3}
\]

The phrase “$Q$ concentrates on $\asymp q^2$ bases” is correct, but the
$P$ in (1.1) is total primitive support and includes the
$q^3-q$ isolated mirror bases.  Thus

\[
P_HQ_H\sim\frac13q^9,qquad N_q^2\asymp q^{10}.    \tag{3.4}
\]

Literal word construction, not just formula evaluation, gives:

| $q$ | $p_q$ | $N_q$ | $H$ | $S_H$ | $P_H$ | $Q_H$ | $P_HQ_H/N_q^2$ |
|---:|---:|---:|---:|---:|---:|---:|:---|
| 5  | 3137   | 3135   | 241  | 720   | 168  | 4600   | $10304/131043=0.078630678480$ |
| 7  | 16811  | 16809  | 673  | 2688  | 432  | 36848  | $1768704/31393609=0.056339619953$ |
| 11 | 161053 | 161051 | 2641 | 15840 | 1560 | 575960 | $7425600/214358881=0.034640972025$ |

All differences, reflection pairs, primitive counts, and row inequalities
$C_d\le3(d-1)$ pass exactly.  The corrected word remains compatible with
(1.1), but with only a $q^{-1}$, rather than $q^{-2}$, margin.

## 4. V1 -- exact first-return renewal and its terminal obstruction

### 4.1 Successor description -- PROVED

For each projective value $v$, list its positions in increasing order:

\[
x_{v,1}<x_{v,2}<\cdots<x_{v,m_v}.
\]

A primitive gap-$d$ collision is exactly a consecutive pair in this list
whose spacing is $d\le D$.  Hence

\[
P_d^{\rm prim}
=\#\{(v,i):x_{v,i+1}-x_{v,i}=d\},                  \tag{4.1}
\]

and

\[
P_D=\sum_v\sum_{i<m_v}
       \mathbf 1_{x_{v,i+1}-x_{v,i}\le D}.          \tag{4.2}
\]

This proves at once that both primitive bases and primitive endpoints are
injective.  Primitive arrows form increasing paths inside the projective
fibres, broken wherever a consecutive spacing exceeds $D$.

### 4.2 Exact cascade intersections -- PROVED

Put

\[
Z_d=\{r:1\le r\le N-d,\ N_d(r)=0\}.
\]

For $a,g\ge1$, the addition law and $N_{a+1}(r)\ne0$ at a physical root
of $N_a$ give

\[
r\in Z_a\cap Z_{a+g}
\quad\Longleftrightarrow\quad
r\in Z_a\ \hbox{ and }\ r+a\in Z_g.                \tag{4.3}
\]

Iterating, for $d_1<\cdots<d_t$,

\[
r\in\bigcap_{i=1}^t Z_{d_i}                         \tag{4.4}
\]

is equivalent to the restarted chain

\[
N_{d_1}(r)=0,\quad
N_{d_2-d_1}(r+d_1)=0,\quad\ldots,\quad
N_{d_t-d_{t-1}}(r+d_{t-1})=0.                       \tag{4.5}
\]

Thus the cascade is an exact description of every higher intersection; it
does not create an inequality in a favorable direction.

Let

\[
k_r=\#\{d\le D:r\in Z_d\},
\]

and define the total $j$-fold intersection mass

\[
I_j(D)=\sum_{d_1<\cdots<d_j}
       \left|Z_{d_1}\cap\cdots\cap Z_{d_j}\right|.
\]

Then exactly

\[
I_j(D)=\sum_r\binom{k_r}{j},                         \tag{4.6}
\]

so $I_1=S_D$, $I_2=Q_D$, and full inclusion--exclusion is

\[
P_D=I_1-I_2+I_3-I_4+\cdots.                          \tag{4.7}
\]

The verifier checks every one of the 404 two-return cascade instances in the
eight requested cells and checks the complete finite sum (4.7), through the
observed maximum multiplicity, in each cell.

### 4.3 Exact shell identities -- PROVED

For a lag shell $J$, put

\[
k_J(r)=\#\{d\in J:r\in Z_d\},\qquad
I_{j,J}=\sum_r\binom{k_J(r)}j.
\]

Let $\mathcal U_J=\bigcup_{d\in J}Z_d$ and $U_J=|\mathcal U_J|$.  Then

\[
U_J=\left|\mathcal U_J\right|
    =\sum_{j\ge1}(-1)^{j+1}I_{j,J}.                  \tag{4.8}
\]

In particular, Bonferroni gives

\[
S_J-I_{2,J}\le U_J\le S_J.                          \tag{4.9}
\]

For $J=(A,B]$, the number whose first return lies in the shell is

\[
F_J=P_B-P_A
   =\left|\mathcal U_J\setminus\mathcal U_{[1,A]}\right|.      \tag{4.10}
\]

Consequently $P_D=\sum_JF_J\le\sum_JU_J$ for any disjoint shell
partition.  This is the exact dyadic mechanism requested in V1.

Its sign is the obstruction.  An upper estimate for cascade intersections
does not upper-bound a union below its row sum.  A useful union saving would
need a lower bound for overlap, whereas the current arithmetic supplies upper
bounds for $Q_D=I_2(D)$.  Adjacent coprimality can make intersections empty,
which makes $U_J=S_J$ and is maximally adverse for this purpose.

### 4.4 V1 verdict and the thinner conditional residual

The numerical shell ledger shows that (4.9) is usually close to its upper,
not lower, endpoint.  Therefore:

**V1 via cascade overlap is DEAD.**  Renewal describes multiplicity exactly
but does not force multiplicity.  What remains is a quenched first-return
census, not a further inclusion--exclusion identity.

A useful residual weaker than (1.2) is the following.  For any fixed
$\eta>0$, suppose one could prove, uniformly in the relevant range,

\[
[\mathrm{FR}_\eta]\qquad P_D\ll D^{2-\eta}.          \tag{4.11}
\]

Then the banked pair bound gives

\[
P_DQ_D\ll D^{4-\eta}(1+\log D)
 =N^{2-\eta/2}L^{4-\eta}(1+\log D).                 \tag{4.12}
\]

Taking, for example, $L=\log\log N$ makes (4.12) $o(N^2)$.  Thus any
power saving over the degree-sum scale for primitive support would already
break the $3/2$ record.  It would suffice to prove the dyadic version

\[
\left|\bigcup_{D/2<d\le D}Z_d\right|\ll D^{2-\eta}  \tag{4.13}
\]

uniformly, since summing its rescalings is geometric.  Neither (4.11) nor
(4.13) is currently banked.

## 5. V2 -- moment and detector route

### 5.1 The first two factorial moments have the wrong one-sided content

The exact zero detector is

\[
\mathbf 1_{k_r>0}
=\sum_{j\ge1}(-1)^{j+1}\binom{k_r}{j}.               \tag{5.1}
\]

Summing (5.1) is (4.7).  If only $S_D$ and $Q_D$ are known, Bonferroni
gives

\[
S_D-Q_D\le P_D\le S_D.                               \tag{5.2}
\]

The second moment improves only the lower bound.  The next upper bound is

\[
P_D\le S_D-Q_D+I_3(D),                               \tag{5.3}
\]

and $I_3(D)$ is a new, uncontrolled cross-gap moment.

There is one exact refinement using the restart cap.  If
$M_D=\max_r k_r$, then for every $1\le k\le M_D$,

\[
1\le k-\frac{2}{M_D}\binom{k}{2}
  =\frac{k(M_D-k+1)}{M_D}.                            \tag{5.4}
\]

Summing proves (1.3).  It is optimal among quadratic majorants
$ak+b\binom{k}{2}$ with $a=1$: the cases $k=1$ and $k=M_D$ force
$a\ge1$ and $b\ge-2/M_D$.

With $M_D\le4D^{2/3}$, (1.3) becomes

\[
P_D\le S_D-\frac{Q_D}{2D^{2/3}}.                    \tag{5.5}
\]

This cannot be used with an upper bound for $Q_D$ to guarantee a
subtraction, and $S_D$ is precisely the unknown first moment.  Bases with
$k_r=1$ attain equality in the raw detector $1\le k$ and are invisible to
every higher factorial moment.

### 5.2 Why the fixed-gap moment theorems do not supply the missing input

Additive orthogonality for one rational map $\delta_d$ writes its zero count
as

\[
R_d=\frac1p\left(|U_d|+
       \sum_{t\ne0}\sum_{x\in U_d}e_p(t\delta_d(x))\right).    \tag{5.6}
\]

Using the banked $(4d-1)\sqrt p$ bound term by term gives only

\[
R_d\le1+O(d\sqrt p),                                 \tag{5.7}
\]

which is weaker than the algebraic degree bound $R_d\le3(d-1)$.  The
Sp-full and Morse certificates through $d=32$ give strong fixed-$d$
information, but only for finitely many rows.  They do not control the OR
detector across the growing family $d\le D$, nor the cross-gap terms
$I_j(D)$.  Fixed-$d$ Chebotarev averages over $p$; V2 needs one fixed
$p$ and growing $d$.

**V2 verdict: DEAD with the available moments.**  Any mollifier that removes
multiplicity either reproduces $P_D\le S_D$, introduces an uncontrolled
higher cross-gap moment, or uses the negative $Q_D$ term in (5.5) without a
lower bound for it.

## 6. V3 -- clique geometry and the primitive-padding obstruction

### 6.1 Best bound from the stated scalar inputs -- PROVED

In the bounded-gap collision graph, $S_D$ is the edge count and $Q_D$ is
the triangle count.  Put

\[
B_D=\frac32D(D-1),\qquad
W_D=66D^2(1+\log D),\qquad c=\frac{\sqrt2}{3}.
\]

The stated inputs are

\[
P_D\le\min(N,S_D),\quad S_D\le B_D,
\quad Q_D\le\min(W_D,cS_D^{3/2}).                    \tag{6.1}
\]

Therefore

\[
P_DQ_D\le
\min\{NW_D,\;NcB_D^{3/2},\;B_DW_D,\;cB_D^{5/2}\}.  \tag{6.2}
\]

The constants in the two spectral branches simplify asymptotically to

\[
NcB_D^{3/2}\sim\frac{\sqrt3}{2}ND^3,qquad
cB_D^{5/2}\sim\frac{3\sqrt3}{4}D^5.                 \tag{6.3}
\]

At $D=\sqrt N L$ with $L\ge1$, the first branch of (6.2) is the
asymptotically smallest one and gives exactly (1.4).  Thus the listed clique
and scalar inputs miss (1.1) by $L^2\log D$.

The structural reason is simple.  A disjoint two-vertex component adds one
to $P_D$ and $S_D$, but adds nothing to $Q_D$ or to any intersection
mass $I_j$, $j\ge2$.  Kruskal--Katona constrains the triangle-rich core;
it cannot prevent arbitrarily many triangle-free primitive components from
being placed beside that core.  In the $k_r$ language, it has no control of

\[
P_{D,1}=\#\{r:k_r=1\}.                                \tag{6.4}
\]

### 6.2 Primitive-padding theorem -- PROVED abstractly

There is an explicit family showing that this is not merely a defect of the
inequality (6.2).

**Theorem.**  There are reflection-symmetric finite words $w^{(q)}$, with
lengths $N_q$ and cutoffs $D_q$, such that

1. $D_q/\sqrt{N_q}\to\infty$;
2. $w_n=w_{N_q+1-n}$ and $w_n\ne w_{n+1}$;
3. $C_d\le3(d-1)$ for every $d\le D_q$;
4. every fibre and every interval $J$ satisfy
   $m_J(v)-1\le4\operatorname{span}(J)^{2/3}$;
5. the exact Apéry $h=2$ count
   $C_2=1+2\mathbf 1_{\left(\frac{-51}{p_q}\right)=1}$ can be imposed;
6. all primitive, cascade, split, and renewal identities hold;
7. nevertheless $P_{D_q}Q_{D_q}/N_q^2\to\infty$.

**Proof.**  Start with the Section 4.4 word at an odd prime $q$, with
$p_q$ the first prime at least $q^5$, $N_q=p_q-2$, and

\[
D_q=2q^3-2q+1.
\]

Bertrand's postulate gives $q^5\le p_q<2q^5$, so

\[
\frac{D_q}{\sqrt{N_q}}\asymp q^{1/2}\longrightarrow\infty.   \tag{6.5}
\]

The left special cluster occupies $[1,D_q+1]$.  Stop all new left-hand
positions at

\[
L_q=\frac{N_q-D_q}{2};                                \tag{6.6}
\]

then every cross-reflection gap from the new region exceeds $D_q$.  The
available length after the special cluster is

\[
A_q=L_q-D_q-1=\frac{N_q-3D_q-2}{2}\asymp q^5.         \tag{6.7}
\]

Let $t_q=\lfloor A_q/(2D_q)\rfloor\asymp q^2$, and choose distinct gaps

\[
d_j=D_q-j,\qquad 0\le j<t_q.                          \tag{6.8}
\]

The $j$-th block has length $2d_j$.  Pair its first $d_j$ positions
with its last $d_j$ positions in order, give each pair a new private
colour, and give the reflected pair the same colour.  The total block length
is at most $2t_qD_q\le A_q$, so all blocks fit before (6.6).

Each block contributes $2d_j$ active bases and edges, but no triangle:
each new colour has one gap-$d_j$ pair on each side, and the two sides are
more than $D_q$ apart.  Thus

\[
P_D^{\rm pad}=2\sum_{j<t_q}d_j\asymp q^5\asymp N_q,  \tag{6.9}
\]

while $Q_D$ is unchanged from (3.2).

For the row bound, the unpadded word has the uniform estimate

\[
C_d^{(0)}\le2q^2+1.                                   \tag{6.10}
\]

On a selected row the padding adds exactly $2d$.  Since
$d\ge D_q-t_q=2q^3-O(q^2)$, for all sufficiently large $q$,

\[
C_d^{(0)}+2d\le3(d-1).                                \tag{6.11}
\]

All unselected rows retain their old bound.  If
$\left(\frac{-51}{p_q}\right)=1$, reserve one further reflected pair of gap
two; it adds two
to the forced $C_2=1$, giving the exact value three.  This $O(1)$ patch
does not affect (6.7)--(6.9).

Reflection and no-adjacent-equality are built into the blocks.  For the
window bound, consecutive special-cluster positions are at least $q+1$
apart.  A subwindow of span $H\le D_q$ therefore contains at most
$1+H/(q+1)\le1+4H^{2/3}$ special positions.  A window meeting both special
clusters has span $\asymp q^5$, making the bound still weaker.  Every new
or private fibre has at most two positions on one side and at most four in
total, so it also satisfies the bound.  Renewal identities hold for every
word.

Finally, (3.2), (6.9), and $N_q\asymp q^5$ give

\[
\frac{P_{D_q}Q_{D_q}}{N_q^2}\asymp
\frac{q^5q^6}{q^{10}}\asymp q\longrightarrow\infty. \tag{6.12}
\]

This proves the theorem. \(\square\)

### 6.3 Finite hostile gates

The verifier constructs the padded words literally.  It does not infer their
statistics from the proof.  Here $D^2/N=L^2$.

| $q$ | $N$ | $D$ | $D^2/N$ | padded primitive bases | total $P_D$ | $Q_D$ | $P_DQ_D/N^2$ |
|---:|---:|---:|---:|---:|---:|---:|:---|
| 5  | 3135   | 241  | 18.526634768740 | 1204  | 1372  | 4600   | $252448/393129=0.642150540917$ |
| 7  | 16809  | 673  | 26.945624367898 | 7394  | 7826  | 36848  | $22182496/21734037=1.020633948493$ |
| 11 | 161051 | 2641 | 43.308523387002 | 76562 | 78122 | 575960 | $33805520/19487171=1.734757702901$ |

The literal runs check every row, all reflections, all adjacent positions,
the $h=2$ law, and 295,760 tight fibre windows.  The maximum multiplicity
remains the original special-cluster value; padding only creates $k_r=1$
bases.  The finite verifier also uses any leftover left-hand capacity after
the top-gap blocks (3, 6, and 15 blocks respectively); this only increases
$P_D$ and is not needed for the asymptotic proof.

**V3 verdict: DEAD from word/clique geometry.**  The padded theorem is a new
abstract obstruction class.  It does not disprove (1.1) for the Apéry orbit;
it identifies the arithmetic coexistence statement that an Apéry proof must
establish.

## 7. Audit of the banked arithmetic

The requested arithmetic inputs were checked for whether they can rule out
primitive padding in the actual orbit.

| Input | What it controls | Why it does not close (1.1) or (1.2) |
|:---|:---|:---|
| Reflection and parity | One forced root in every even row; pairs all other row roots | Gives a lower baseline and isolated primitive edges, not an upper bound for noncentral first returns |
| Exact $h=2$, $-51$ split | $C_2=1$ or $3$ | One fixed row; the padded obstruction can satisfy it exactly |
| Leading-coefficient apparition | Identifies rows where the nominal degree drops | Does not bound roots of the remaining polynomial or any nonapparition row |
| Adjacent and shifted resultants | Detect same-base common roots and hence cascade intersections | These are $Q$-type events; absence of a common root increases distinctness and leaves $P$ uncontrolled |
| $N_h$ squarefree for all $h$ | Removes repeated characteristic-zero roots | Squarefree polynomials may have many distinct roots modulo one fixed prime |
| Morse/irreducibility/Sp certificates through $h=32$ | Fixed-height geometric monodromy and moments | A finite set of rows is negligible when $D\to\infty$; no growing-height fixed-prime sieve follows |
| Per-height $(4h-1)\sqrt p$ sums | Additive cancellation for one $\delta_h$ | The zero detector gives (5.7), worse than the degree bound |
| Saturated mesopair diagnostics | Empirically isolate the true balanced $Q_D$ support | Finite evidence only; the missing uniform balanced pair theorem would itself solve the wall |
| Transfer codegree theorem | Annealed two-clock mixing | First returns use the ordered, locked clock; the known quenched-order gap remains |

The exact missing input can be attacked in either of two useful strategic
forms.

1. **First-return census:** prove (4.11), or ideally
   $P_D\ll D\operatorname{polylog}D$, for one supercritical range.
2. **Simple/core coexistence:** with $P_{D,1}$ from (6.4), prove
   
   \[
   P_{D,1}Q_D\ll N^2                                      \tag{7.1}
   \]
   
   together with the analogous estimate on the $k_r\ge2$ core.

Neither is a consequence of the current resultants or fixed-height moment
certificates.  They are quenched growing-gap arithmetic statements.

## 8. Exact machine audit

Run:

```sh
PYTHONPYCACHEPREFIX=/tmp/codex-strike2-pycache \
  python3 -m py_compile CODEX_STRIKE2_verify.py
python3 CODEX_STRIKE2_verify.py
```

The final verifier SHA-256 is

```text
dc43e18a1900469e330d3ce7ad54f612d4dff3cd0ecbcd221c3cae7b141d54a7
```

The run ends with

```text
ALL STRIKE2 GATES PASS
```

The persistent gates include:

- exact orbit/continuant agreement on 8,782,308 physical cells;
- reflection at every orbit position;
- injectivity and consecutiveness of every primitive endpoint;
- every two-return cascade, with total count exactly $Q_D$;
- complete finite inclusion--exclusion in every requested cell;
- the optimal quadratic detector and squared spectral triangle bound;
- every dyadic shell sum, union, pair intersection, and first-return count;
- literal construction of the three original reflection words;
- literal construction of the three padded hostile words;
- all row bounds, $h=2$ laws, and 295,760 tight-window multiplicity checks.

## 9. Status ledger

| Statement or route | Status | Terminal reason |
|:---|:---|:---|
| V4 requested calibration | `VERIFIED-N` | Two independent collision constructions; all exact tables above |
| Original Section 4.4 $P,Q$ count | `PROVED`, `VERIFIED-N` | $P\asymp q^3$, $PQ\asymp q^9$, still $o(N^2)$ |
| Primitive = fibre successor | `PROVED`, `VERIFIED-N` | Exact ordered-fibre description (4.1)--(4.2) |
| Higher cascade/intersection identities | `PROVED`, `VERIFIED-N` | Iterated restart gives (4.5)--(4.7) |
| V1 union saving from cascade | `DEAD` | Overlap is not forced and has the wrong one-sided sign |
| $[\mathrm{FR}_\eta]$ | `CONDITIONAL` | Any $\eta>0$ breaks $3/2$ for a sufficiently slow unbounded $L$ |
| V2 fixed moments/mollified detector | `DEAD` | Best upper detector is (1.3), still containing unknown $S_D$ |
| Scalar/clique product envelope | `PROVED` | Gives only (1.4), off by $L^2\log D$ |
| Primitive-padding word theorem | `PROVED`, `VERIFIED-N` | Refutes word-level implication to (1.1) |
| Current arithmetic rules out padding in Apéry | `OPEN` | Requires a quenched first-return or simple/core coexistence theorem |
| (4.8) $P_DQ_D\ll N^2$ for Apéry | **NOT PROVED** | No growing-gap arithmetic coupling is available |
| (4.9) primitive-support bound | **NOT PROVED** | Empirical linearity in $D$ has no uniform fixed-prime proof |
| Any unbounded-$L$ improvement of $3/2$ | **NOT PROVED** | The last arithmetic wall remains |
