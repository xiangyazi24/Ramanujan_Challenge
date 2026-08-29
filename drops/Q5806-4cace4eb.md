ANSWER Q5806 4cace4eb

# Problem 3.2 — cold audit of the claimed unconditional \(H^3T^2=N^{12/5}\) affine-mass bound

## Verdict

\[
oxed{	extbf{FAIL}.}
\]

The local affine Fourier inequality used by Q5767 is correct, and the elementary bound

\[
\sqrt{E_
u}\le M_
u
\]

is correct. For the **two-gap** family that actually underlies Q5605's localized cubic-moment argument, the pointwise estimate

\[
M_
u\le H^{2+o(1)}
\]

is also legitimate **provided that** \(
u\) fixes the characteristic, the rational residual root/local factor, the primary level, and every unbounded outer coordinate.

The claimed global consequence does not follow. Its first false numerical line is

\[
oxed{\sum_
u w_
u\ll H^{1+o(1)}.}
	ag{F}
\]

In Q5767, \(
u=(p,\mathfrak q,r,	au_
u)\), with \(\mathfrak q=X-x\) for a rational root and \(w_
u=\log p\) at degree one. Thus \(
u\) ranges over the prime shell, rational residual roots, primary levels, and bounded chart data. None of the canonical sources proves that the total weight of those cells is \(H^{1+o(1)}\). The only source-valid unconditional aggregate supplied by the Q5605 Sylvester/Smith ledger is

\[
\sum_
u w_
u M_
u\ll H^{4+o(1)},
	ag{C1}
\]

not \(H^{3+o(1)}\). Consequently the trivial affine-energy majorant gives only

\[
\sum_
u w_
u\sqrt{E_
u}\le\sum_
u w_
u M_
u\ll H^{4+o(1)},
	ag{C2}
\]

which recovers Q5605's existing dyadic line

\[
H^4T^2=N^4T^{-2}\le N^{14/5+o(1)}\qquad\left(N^{3/5}\le T\le N^{2/3}ight).
	ag{C3}
\]

The algebraic substitution

\[
H^3T^2=N^3/T\le N^{12/5}
\]

is arithmetically correct **conditional on a new aggregate theorem** of strength \(H^{3+o(1)}\). Q5767 did not prove that aggregate theorem; it replaced it by the incompatible pair of assertions “\(M_
u\le H^2\)” and “\(\sum_
u w_
u\le H\).” Fixing the residual root makes the first assertion available but forces every prime/root/level cell into the second sum. Omitting the root from \(
u\) would avoid that exact sum only by putting the root multiplicity back into \(M_
u\).

There is a second, logically earlier source gap. Q5674's literal distinguished-state carrier tests \(B_x=0\) at a rational common root. Q5739/Q5767 instead introduce an abstract two-coordinate terminal projection \(P_\omega T_\omega(1,5)^t\) and test that vector for zero. No canonical Problem 3.2 source currently states the one-to-one, multiplicity-preserving theorem identifying these vector-zero occurrences with Q5605's selected actual Apéry triple starts. The affine origin inequality is valid as an abstract finite-group identity, but its exact downstream consumer binding is not yet a canonical theorem.

Accordingly:

```text
new unconditional N^(12/5) cubic theorem:           FAIL
local Fourier/origin inequality:                     PASS
sqrt(E_nu) <= M_nu:                                  PASS
M_nu <= H^(2+o(1)):                                  PASS, two-gap/fixed-root scope only
sum_nu w_nu <= H^(1+o(1)):                           FAIL / not proved
sum_nu w_nu M_nu <= H^(3+o(1)):                      FAIL / not proved
source-valid unconditional aggregate:                H^(4+o(1))
source-valid cubic consequence:                       N^(14/5+o(1))
raw three-gap H^3 family covered by Q5767 ledger:     NO
canonical insertion of an N^(12/5) theorem:           NOT JUSTIFIED
```

No authoritative TeX was edited.

---

## 0. Source boundary and pins

This audit uses only Problem 3.2 in the canonical repository and the exact connected same-project records named below. No workspace/Zinan repository was used.

### Canonical repository

```text
repository: xiangyazi24/Ramanujan_Challenge
branch: main
HEAD: c5d932b66ce5e4f1657b587215d290ae7a13018b
tree: 1abbf9eccd470122f782d2ab82d7788595558a7a
```

Pinned canonical files:

```text
problems/3.2/proof.tex
  blob efbede7ea8ac6e040a5d380860ef5009a564fb01

problems/3.2/atom_tail_section.tex
  blob 96051be65d3be3d01d6cf42cb6552f95662f963e

problems/3.2/projective_variance_reduction.tex
  blob 4ac4e34709be27db1ea7d8378442d45779fedca1

problems/3.2/energy_result.tex
  canonical at the same HEAD

problems/3.2/aligned_corank_verify.py
  blob 340ebbb82a11d464e49b69951d6e020a7276aaab

problems/3.2/fully_deflated_corank_verify.py
  blob d325a3b96757ea85cd8575b98e1d1cddd23baf77

problems/3.2/primitive_fd_candidate_verify.py
  blob 1f200891ae67ce901460071bf759e8004893a5c1

problems/3.2/centered_residual_verify.py
  blob d7cd445cde12dc491ff107c0d9f1ae832b53d3bd

problems/3.2/primitive_projective_prime_scan.cpp
  blob 0b1c106ad5b449affc44c7a851aa141e73c6d11d
```

`proof.tex` includes `oracleA_result`, `hm3_result`, and `atom_tail_section`; it does not include Q5605, Q5739, or Q5767 as authoritative source modules.

### Connected research records audited

```text
Q5605 548ceea6
  Notion page 3cb7a6c4-fa84-812d-b609-e83934dbcbd5
  localized triple-resultant / cubic-moment ledger

Q5674 82e67fb3
  Notion page 3cb7a6c4-fa84-8131-a48bf0b8792189f5
  actual distinguished-state carrier and primitive-return masks

Q5739 02a18711
  Google Doc 1MFvjs1pjm3WNloXASVnrqA9x6H078WJ8AJom9iZbxfI
  module A_Y/(B_Y,C_Y), local length, and abstract affine variance

Q5767 41375002
  Google Doc 1T8h-IdP234meXEztUXzklHYKW4pIPgRVYrUCzpekpg0
  claimed H^3 affine fallback and N^(12/5) consequence
```

Q5739 contains an important scope warning: it explicitly says that the literal Q5605 normalization was not exposed to that audit and that it was using the interface stated in its prompt. Its phrase “in a shell with \(H^{1+o(1)}\) local levels” is therefore a conditional template, not a theorem about Q5605's actual prime/root/level index set.

---

## 1. The literal Q5605 quantity

### 1.1 Dyadic prime and zero-count shells

For a fixed dyadic characteristic block, put

\[
\mathcal P(N)=\{p	ext{ prime}:N<p\le2N\},\qquad \mathcal Z_p=\{0\le x<p:b_x\equiv0\pmod p\},\qquad Z(p)=|\mathcal Z_p|.
\]

For a dyadic zero-count height \(T\), Q5605 uses

\[
\mathcal P_T=\{p\in\mathcal P(N):T<Z(p)\le2T\},
	ag{1.1}
\]

in the licensed high range

\[
N^{3/5}\lesssim T\lesssim N^{2/3}.
	ag{1.2}
\]

It chooses

\[
H=\left\lceilrac{16N}{T}ightceil,\qquad H\asymprac NT.
	ag{1.3}
\]

The factor \(16\) is irrelevant to exponents but ensures that the block loss below is a fixed fraction of \(T\).

### 1.2 Selected actual zero triples

Partition \([0,p-1]\) into consecutive intervals of length at most \(H\). If one interval contains actual Apéry zeros

\[
z_1<z_2<\cdots<z_m,
\]

select its consecutive triples \((z_i,z_{i+1},z_{i+2})\). Their total number \(Q_p(H)\) satisfies

\[
Q_p(H)\ge Z(p)-2\left\lceilrac pHightceil.
	ag{1.4}
\]

For \(p\in\mathcal P_T\) and the choice (1.3),

\[
Q_p(H)\gg T.
	ag{1.5}
\]

Every selected triple has the unique representation

\[
r<r+a<r+a+b<p,\qquad a,b\ge2,\qquad a+b<H.
	ag{1.6}
\]

Define

\[
t_p^{\mathrm{Ap}}(a,b):=\#\{r:0\le r<r+a<r+a+b<p,\ r,r+a,r+a+b\in\mathcal Z_p,\ 	ext{and the Q5605 block-selection conditions hold}\}.
	ag{1.7}
\]

Then

\[
Q_p(H)=\sum_{a,b\ge2,\ a+b<H}t_p^{\mathrm{Ap}}(a,b).
	ag{1.8}
\]

This is a **two-gap** family. The start \(r\) is an occurrence/root variable, not a third gap parameter.

### 1.3 The positive primary mass

The logarithmically weighted mass is

\[
\mathfrak M_{\mathrm{Ap}}(N,H):=\sum_{a,b\ge2,\ a+b<H}\ \sum_{N<p\le2N}t_p^{\mathrm{Ap}}(a,b)\log p.
	ag{1.9}
\]

By (1.5)–(1.8),

\[
T|\mathcal P_T|\log N\ll\mathfrak M_{\mathrm{Ap}}(N,H).
	ag{1.10}
\]

This is the first place where the prime-shell weight is paid. The later unweighted cubic moment loses the single \(\log N\) again; it does not create an additional prime sum.

### 1.4 The source-valid Sylvester/Smith envelope

A selected actual triple forces

\[
N_a(r)\equiv0\pmod p,\qquad N_b(r+a)\equiv0\pmod p.
	ag{1.11}
\]

Let

\[
S_{a,b}=\operatorname{Res}_Xigl(N_a(X),N_b(X+a)igr).
	ag{1.12}
\]

The canonical root-strip theorem makes \(S_{a,b}
e0\). Define

\[
t_p^{\mathrm{form}}(a,b):=\#\{x\in\mathbf F_p:N_a(x)=N_b(x+a)=0\},
	ag{1.13}
\]

with the relevant safe central/structural restrictions inserted when the deflated form is used. Then

\[
t_p^{\mathrm{Ap}}(a,b)\le t_p^{\mathrm{form}}(a,b)\le v_p(S_{a,b}).
	ag{1.14}
\]

The second inequality is the fixed-size Sylvester/Smith corank argument. It pays distinct common roots and repeated local multiplicity; no root-simplicity assumption is needed.

For a fixed gap pair,

\[
\sum_{N<p\le2N}v_p(S_{a,b})\log p\le\log|S_{a,b}|.
	ag{1.15}
\]

The canonical height estimate is

\[
\log|S_{a,b}|\ll(a+b)^2\log(a+b).
	ag{1.16}
\]

Consequently

\[
\mathfrak M_{\mathrm{Ap}}(N,H)\le\sum_{a,b\ge2,\ a+b<H}\log|S_{a,b}|\ll\sum_{s\le H}s^3\log s\ll H^{4+o(1)}.
	ag{1.17}
\]

This \(H^4\) mass is what Q5605 actually banks unconditionally.

### 1.5 Why the cubic contribution has a factor \(T^2\)

On \(\mathcal P_T\), \(Z(p)^3\ll T^3\). From (1.10) and (1.17),

\[
|\mathcal P_T|\llrac{H^{4+o(1)}}{T\log N}.
	ag{1.18}
\]

Therefore

\[
\sum_{p\in\mathcal P_T}Z(p)^3\ll T^3|\mathcal P_T|\ll H^{4+o(1)}T^2.
	ag{1.19}
\]

The \(T^2\) is exactly

```text
T^3 from the cubic zero count
minus
one factor T supplied by the selected-triple lower bound Q_p(H) >> T.
```

Using \(H\asymp N/T\),

\[
H^4T^2=N^4T^{-2}.
	ag{1.20}
\]

This decreases throughout (1.2), so the worst endpoint is \(T=N^{3/5}\):

\[
N^4N^{-6/5}=N^{14/5}.
	ag{1.21}
\]

Thus the prior \(N^{14/5}\) term counted the heavy-prime contribution to the **third zero moment**, via actual consecutive zero triples and their two-gap Smith-primary envelope. It did not count a free raw three-gap family.

---

## 2. Reindexing the formal occurrences correctly

Let

\[
\mathcal Y_H=\{(a,b):a,b\ge2,\ a+b<H\},\qquad |\mathcal Y_H|=H^{2+o(1)}.
	ag{2.1}
\]

For \(Y=(a,b)\), a characteristic \(p\), and a rational common root \(x\), let \(e_{p,Y,x}\) be the local multiplicity retained after the literal safe masks. A fully expanded formal primary occurrence has the shape

\[
\iota=(Y,p,x,r,	au),\qquad 1\le r\le e_{p,Y,x},
	ag{2.2}
\]

where \(	au\) records bounded orientation, source, chart, and physical labels.

The formal occurrence set must retain:

```text
p in the fixed dyadic characteristic shell;
x in F_p, excluding extension roots;
Y=(a,b), with a,b>=2 and a+b<H;
nonwrapping;
the two shell/gap equations at the same x;
primitive/first-return masks when the consumer uses them;
center and structural exceptions split exactly as in the source;
primary level r;
all bounded occurrence labels exactly once.
```

Let

\[
s(\iota)=P_\iota T_\iota(1,5)^t
	ag{2.3}
\]

be the proposed affine state.

To perform Fourier analysis, occurrences may be grouped only when their values lie in the same finite additive group. For degree-one factors a favorable grouping is by

\[

u=(p,x,r,	au_
u),\qquad R_
u\simeq\mathbf F_p[\epsilon]/(\epsilon^r),
	ag{2.4}
\]

with every local polynomial expressed in its own parameter \(X-x\) and then transported to the canonical model. Put

\[
\Omega_{
u,H}:=\{(Y,	au):(Y,p,x,r,	au)	ext{ is a formal occurrence in cell }
u\},\qquad M_
u:=|\Omega_{
u,H}|.
	ag{2.5}
\]

At degree one the layer weight is, up to bounded normalization,

\[
w_
u=\log p.
	ag{2.6}
\]

The exact incidence identity is

\[
oxed{\sum_
u w_
u M_
u=\sum_{Y\in\mathcal Y_H}\ \sum_{N<p\le2N}\ \sum_{x\in\mathbf F_p}\ \sum_{1\le r\le e_{p,Y,x}}\ \sum_	au(\log p)\,\mathbf1_{\{	ext{all formal masks}\}}.}
	ag{2.7}
\]

The right side is the total formal prime/root/level occurrence mass. Regrouping cannot change it.

If the shell equations in (2.7) are the same two common-root equations paid by the Q5605 Sylvester matrix, then Smith normal form gives the unconditional envelope

\[
\sum_
u w_
u M_
u\ll H^{4+o(1)}.
	ag{2.8}
\]

If the proposed terminal projection introduces a different local scheme not yet identified with that common-root scheme, then even (2.8) for the new object requires a bridge theorem. In neither interpretation does the source give \(H^3\).

---

## 3. Audit of the pointwise bounds

### 3.1 \(\sqrt{E_
u}\le M_
u\): PASS

Let

\[
n_
u(z):=\#\{\omega\in\Omega_{
u,H}:s(\omega)=z\},\qquad z\in R_
u^2,
	ag{3.1}
\]

and

\[
C_
u:=\sum_z n_
u(z)^2,\qquad E_
u:=C_
u-rac{M_
u^2}{|R_
u|^2}.
	ag{3.2}
\]

Then

\[
0\le E_
u\le C_
u\le\left(\sum_z n_
u(z)ight)^2=M_
u^2,
	ag{3.3}
\]

so

\[
\sqrt{E_
u}\le M_
u.
	ag{3.4}
\]

This step is unconditional once the finite occurrence set and common local ring are correctly defined.

### 3.2 \(M_
u\le H^{2+o(1)}\): PASS only in fixed-root two-gap scope

If \(
u\) fixes \(p,x,r\) and every unbounded non-gap variable, then the only free large parameters in Q5605 are

\[
(a,b)\in\mathcal Y_H.
\]

Hence

\[
M_
u\le|\mathcal Y_H|H^{o(1)}\le H^{2+o(1)}.
	ag{3.5}
\]

The \(H^{o(1)}\) may absorb only labels already proved to have bounded/subpolynomial multiplicity. It may not absorb a new root, characteristic, host, or gap parameter.

### 3.3 Why the root cannot be omitted for free

If \(x\) is omitted from \(
u\), then \(\Omega_
u\) must include every rational root belonging to each gap pair. Since a gap polynomial has degree \(O(H)\), the elementary bound becomes

\[
M_
u\le H^{3+o(1)},
	ag{3.6}
\]

not \(H^{2+o(1)}\). More importantly, states at different roots live in differently based local rings until a canonical transport is specified.

Thus there is an unavoidable bookkeeping dichotomy:

```text
include x in nu:
    M_nu <= H^(2+o(1)),
    but every active prime/root/level cell appears in sum_nu w_nu;

omit x from nu:
    root multiplicity reappears inside M_nu,
    and a common local additive group is no longer automatic.
```

Q5767 uses the first choice for the pointwise bound and then reasons as though it had used the second choice in the level count.

### 3.4 Raw three-gap family: not covered

The canonical fully deflated four-return family has

\[
a,b,c\ge2,\qquad a+b+c\le H,
	ag{3.7}
\]

and hence \(H^{3+o(1)}\) raw gap triples. The canonical aligned-pencil definition also retains three gap symbols; one gets only two effective parameters after a **specific proved alignment relation** fixes one combination.

Q5605's localized cubic consumer is the two-gap triple-start family, so (3.5) is correct for that consumer. It cannot be exported to the raw three-gap/four-return source. For the raw family the corresponding fixed-root trivial bound is

\[
M_
u\le H^{3+o(1)}.
	ag{3.8}
\]

Q5739's phrase “aligned three-gap system with two free parameters” is an abstract scope statement. It is not a canonical theorem supplying a deterministic \(c=c(a,b,\ldots)\) relation for every raw three-gap consumer.

---

## 4. The fatal total-weight error

Q5767 writes

\[
\sum_
u w_
u\ll H^{1+o(1)},\qquad M_
u\le H^{2+o(1)},
	ag{4.1}
\]

and concludes

\[
\sum_
u w_
u M_
u\ll H^{3+o(1)}.
	ag{4.2}
\]

The second implication would be valid if both premises were proved. The first premise is not a Q5605 fact.

### 4.1 Literal expansion of the alleged level weight

For the degree-one actual-root decomposition,

\[
\sum_
u w_
u=\sum_{N<p\le2N}\ \sum_{x\in\mathbf F_p}\ \sum_{r\ge1}\ \sum_{	au_
u:\exists Y	ext{ formal in }(p,x,r,	au_
u)}\log p.
	ag{4.3}
\]

This sum includes:

```text
all active characteristics p in (N,2N];
all active rational residual roots x;
all primary layers r;
all separately retained charts/branches.
```

There is no canonical estimate reducing (4.3) to the degree \(O(H)\) of one gap polynomial. Degree \(O(H)\) is a **per-system** statement; Q5605 has \(H^2\) gap systems and an unrestricted dyadic prime shell.

Since every active cell has \(M_
u\ge1\), one has merely

\[
\sum_
u w_
u\le\sum_
u w_
u M_
u.
	ag{4.4}
\]

The source bounds the right side by \(H^{4+o(1)}\), not by \(H\).

### 4.2 The Smith ledger shows exactly where the missing factors live

For one pair \(Y=(a,b)\), local prime/root/level multiplicity is paid by

\[
\sum_{N<p\le2N}\ \sum_x\ \sum_{r\le e_{p,Y,x}}\log p\le\sum_{N<p\le2N}v_p(S_{a,b})\log p\le\log|S_{a,b}|.
	ag{4.5}
\]

The right side is \(H^{2+o(1)}\) in the worst pair range. Summing over the \(H^{2+o(1)}\) pairs gives \(H^{4+o(1)}\).

Q5767 effectively kept the \(H^2\) pair count inside \(M_
u\), kept one \(H\) as an alleged number of levels, and dropped the remaining prime/root/valuation height. That dropped factor is not a logarithm; in the source ledger it is a full power of \(H\) at aggregate scale.

### 4.3 Prime-shell and residue-degree audit

- **Prime shell.** The audit is already at one fixed dyadic \(N\)-shell, so summing over dyadic shells costs only \(N^{o(1)}\). But the primes *inside* \((N,2N]\) are not one level. They are paid through (4.5).
- **Rational roots.** Each \(x\in\mathbf F_p\) is a separate degree-one local factor. Fixing \(x\) is what makes (3.5) true; it cannot then be omitted from (4.3).
- **Primary levels.** Summing \(r=1,\ldots,e_{p,Y,x}\) gives the local length. Smith valuation pays this sum. Replacing all levels by one cell would lose the exact primary multiplicity.
- **Residue degree.** Q5605's actual starts use degree-one factors. If higher-degree factors are included, \(w_
u=(\deg\mathfrak q)\log p\) and extension-root pollution is introduced; this cannot improve (4.3).
- **Charts and physical labels.** Source-proved bounded multiplicities may be absorbed in \(H^{o(1)}\). They do not repair the missing prime/root mass.

### 4.4 A legal interface obstruction

The claimed implication is refuted even at the abstract deterministic interface. Consider a permitted incidence profile with \(H^{4-o(1)}\) active cells \(
u\), each having

\[
M_
u=1,\qquad |R_
u|	o\infty.
\]

Then

\[
M_
u\le H^2,\qquad E_
u=1-|R_
u|^{-2},\qquad\sqrt{E_
u}=1+o(1),
\]

but

\[
\sum_
u w_
u\sqrt{E_
u}=H^{4-o(1)}
\]

when logarithmic weights are normalized by \(\log N\). This profile respects the Q5605 aggregate height ceiling and every pointwise inequality used by Q5767. It is not asserted to be an Apéry realization; it proves that the displayed scalar inputs do not imply the claimed \(H^3\) conclusion.

The canonical sources contain no covariance/overlap theorem forcing those active cells to collapse into only \(H^{1+o(1)}\) common local levels.

---

## 5. The affine origin identity and the downstream consumer

### 5.1 The finite-group identity: PASS

For a correctly defined finite occurrence set in one ring \(R_
u\), two-dimensional character orthogonality gives

\[
E_
u=rac1{|R_
u|^2}\sum_{\lambda
e0}|S_
u(\lambda)|^2,
	ag{5.1}
\]

and

\[
n_
u(0)\lerac{M_
u}{|R_
u|^2}+\sqrt{E_
u}.
	ag{5.2}
\]

Defining \(\Omega_{
u,H}\) **before** imposing \(s(\omega)=0\) is the correct noncircular order. If \(\Omega\) were pre-restricted to the distinguished origin, every surviving occurrence would lie in one state fibre and the energy statement would assume the conclusion.

Restricting \(x\) to \(\mathbf F_p\) also correctly removes extension-root pollution.

### 5.2 The mean term must remain in the consumer

The aggregate local-length mass has the form

\[
\mathcal M_{\mathrm{aff}}(N,H)=\sum_
u w_
u n_
u(0).
	ag{5.3}
\]

Equation (5.2) gives

\[
\mathcal M_{\mathrm{aff}}(N,H)\le\sum_
u w_
urac{M_
u}{|R_
u|^2}+\sum_
u w_
u\sqrt{E_
u}.
	ag{5.4}
\]

Q5767's \(H^3\) headline discusses only the second term. Under the source-valid formal envelope (2.8), the mean is harmless but must be shown:

\[
|R_
u|\ge p\asymp N\quad\Longrightarrow\quad\sum_
u w_
urac{M_
u}{|R_
u|^2}\le N^{-2}\sum_
u w_
u M_
u\ll N^{-2}H^{4+o(1)}.
	ag{5.5}
\]

After the outer \(T^2\), this is

\[
N^{-2}H^4T^2=N^2T^{-2}\le N^{4/5+o(1)},
	ag{5.6}
\]

so it is not the obstruction. The obstruction is the energy term's aggregate formal mass.

### 5.3 Exact consumer identification: not established in canonical source

Q5605's literal selected occurrence is an actual zero triple:

\[
b_x\equiv b_{x+a}\equiv b_{x+a+b}\equiv0\pmod p.
	ag{5.7}
\]

Q5674 refines a formal common root by retaining the distinguished initial slope. In its literal normalization, with a companion solution \(C_x\), the state resolvent is

\[
\Theta_{p,X}(Y,Z)=\prod_{x\in X(\mathbf F_p)}(C_xY-B_xZ),
	ag{5.8}
\]

and the actual Apéry start is selected by \(B_x=0\). The companion coordinate is generally a unit rather than zero.

Q5739 uses different symbols \(B_Y,C_Y\): they are asserted to be two **terminal residual coordinates after a projection** \(P_\omega\), and the desired event is

\[
P_\omega T_\omega(1,5)^t=(B_Y,C_Y)^t=0.
	ag{5.9}
\]

These are not the same two coordinates as the fundamental Apéry/companion row in Q5674. That distinction is necessary—an invertible raw transfer cannot send the primitive initial vector to zero—but it creates a proof obligation.

The canonical repository currently contains no multiplicity-preserving theorem asserting

```text
Q5605 selected actual start (p,x,a,b)
if and only if
one fully masked Q5739 occurrence omega has s(omega)=0.
```

Q5739 itself records that it was not reading the literal Q5605 normalization.

Therefore (5.2) is a valid abstract reduction, but Q5767 does not yet establish that its left side is exactly the Q5605 mass (1.9). This is an additional reason that no canonical \(N^{12/5}\) theorem may be inserted.

---

## 6. Canonical variance inputs do not supply the missing level count

### 6.1 Projective variance

`projective_variance_reduction.tex` proves, for four-consecutive windows in one projective fibre, statements including \(M_p\le p\) when the first index determines the window, and the short-chain energy bound

\[
\sum_{q\in\mathbf P^1(\mathbf F_p)}C_p(q)^2\le30pH.
	ag{6.1}
\]

It also proves a mass-sensitive variant

\[
\sum_q C_p(q)^2\le20p\sqrt{M_p}.
	ag{6.2}
\]

These estimates do not imply (F):

1. they are per-characteristic and retain a factor \(p\);
2. they concern nonzero projective directions;
3. the distinguished zero vector has no projective image;
4. they forget the radial primary level \(r\);
5. they concern four-return projective windows, not the complete two-gap affine-origin mass.

At the endpoint \(H\asymp p^{1/3}\), the canonical bound \(pH\) is \(H^4\), exactly the trivial square scale for an \(H^2\)-sized family. It cannot be converted into an \(H\)-sized set of affine local levels.

### 6.2 Canonical continuant collision energy

`energy_result.tex` proves an unconditional formal-column bound

\[
\mathcal E_p(H)\ll H^{8/3}
\]

for a different energy built from pairs of zero levels \(N_h(x)=0\) in a continuant column. It explicitly distinguishes actual affine witnesses from raw resultant support and warns that boundary/resultant-supported rows need not have a common rational witness.

That theorem is not an affine-state collision estimate for \(P_\omega T_\omega(1,5)^t\). It supplies neither (F) nor a primary-level/root count for the Q5605 mass.

### 6.3 Aligned three-gap content

`atom_tail_section.tex` defines the raw aligned pencil with three gaps and proves a valuation statement for common roots of three polynomials. Its desired aggregate \(H^{3+o(1)}\) aligned-content mass is stated as a hypothesis, not an unconditional theorem. Exact computations support the hypothesis at small heights but are not asymptotic proof.

Thus no canonical projective, continuant-energy, or aligned-content statement fills the missing line in Q5767.

---

## 7. Correct replacement ledger

### 7.1 Local affine statement

The following remains correct:

\[
n_
u(0)\lerac{M_
u}{|R_
u|^2}+\sqrt{E_
u},\qquad\sqrt{E_
u}\le M_
u.
	ag{7.1}
\]

### 7.2 Aggregate statement available from current source

Under the most favorable source-aligned interpretation of the shell equations,

\[
\sum_
u w_
u M_
u\ll H^{4+o(1)}.
	ag{7.2}
\]

Therefore

\[
\sum_
u w_
u\sqrt{E_
u}\ll H^{4+o(1)},
	ag{7.3}
\]

and

\[
\mathcal M_{\mathrm{aff}}(N,H)\ll H^{4+o(1)}.
	ag{7.4}
\]

This is not a new estimate; it is the existing formal Smith envelope rewritten in affine-energy notation.

### 7.3 Downstream cubic moment

The correct unconditional Q5605 chain remains

\[
T|\mathcal P_T|\log N\ll\mathfrak M_{\mathrm{Ap}}(N,H)\ll H^{4+o(1)},
	ag{7.5}
\]

and

\[
\sum_{p\in\mathcal P_T}Z(p)^3\ll H^{4+o(1)}T^2=N^4T^{-2+o(1)}\le N^{14/5+o(1)}.
	ag{7.6}
\]

No fixed power below \(14/5\) follows from the trivial affine-energy inequality alone.

### 7.4 What would be needed for \(N^{12/5}\)

A genuinely new theorem of either of the following forms would suffice:

\[
\sum_
u w_
u M_
u\ll H^{3+o(1)},
	ag{7.7}
\]

or, more directly,

\[
\mathfrak M_{\mathrm{Ap}}(N,H)\ll H^{3+o(1)}.
	ag{7.8}
\]

Together with the exact consumer bridge, either gives

\[
\sum_{p\in\mathcal P_T}Z(p)^3\ll H^3T^2=N^3/T\le N^{12/5+o(1)}.
	ag{7.9}
\]

But (7.7) is not a consequence of \(M_
u\le H^2\); it is itself a global prime/root/level sparsity theorem. Equation (7.8) is an actual-primary theorem of the same strength. Neither is currently canonical.

---

## 8. Full exponent ledger and sigma thresholds

Let a valid aggregate primary-mass estimate be

\[
\mathfrak M(N,H)\ll H^{\mu+o(1)}.
	ag{8.1}
\]

Then the Q5605 heavy-tail contribution is

\[
H^\mu T^2=N^\mu T^{2-\mu+o(1)}.
	ag{8.2}
\]

For \(T=N^t\), \(3/5\le t\le2/3\), its exponent is

\[
f_\mu(t)=\mu+(2-\mu)t.
	ag{8.3}
\]

| Mass exponent | Worst endpoint | Cubic exponent | Status |
|---|---:|---:|---|
| \(\mu=4\) | \(t=3/5\) | \(14/5\) | proved by current Smith ledger |
| \(\mu=3\) | \(t=3/5\) | \(12/5\) | conditional; Q5767 did not prove it |
| \(\mu=2\) | all \(t\) | \(2\) | desired \(H^2\) primary-mass scale |
| \(\mu=2-\delta\) | \(t=2/3\) | \(2-\delta/3\) | requires a subquadratic mass theorem |

The staged cubic targets recorded in Q5605 are

\[
N^{29/15-6\sigma+o(1)}
	ag{8.4}
\]

and the stronger

\[
N^{2-6\sigma+o(1)}.
	ag{8.5}
\]

For a subquadratic mass \(H^{2-\delta+o(1)}\), (8.2) reaches (8.4) exactly when

\[
2-rac{\delta}{3}\lerac{29}{15}-6\sigma\quad\Longleftrightarrow\quadoxed{\delta\gerac15+18\sigma},
	ag{8.6}
\]

and reaches (8.5) exactly when

\[
2-rac{\delta}{3}\le2-6\sigma\quad\Longleftrightarrow\quadoxed{\delta\ge18\sigma}.
	ag{8.7}
\]

For comparison, the incorrectly claimed \(N^{12/5}\) line would still exceed \(N^{29/15-6\sigma}\) by exponent \(7/15+6\sigma\), exceed \(N^{2-6\sigma}\) by \(2/5+6\sigma\), and exceed the RR-style target \(N^{11/5-2\sigma}\) by \(1/5+2\sigma\).

The source-valid \(N^{14/5}\) line exceeds those three targets by, respectively,

\[
rac{13}{15}+6\sigma,\qquadrac45+6\sigma,\qquadrac35+2\sigma.
	ag{8.8}
\]

Thus even a future proof of the conditional \(N^{12/5}\) fallback would be progress, not full closure.

---

## 9. PASS/FAIL table

| Audit item | Verdict | Reason |
|---|---|---|
| Canonical repository and include order pinned | PASS | `proof.tex` and included Problem 3.2 files were read at `c5d932b...` |
| Q5605 prior quantity reconstructed | PASS | heavy-prime third moment via selected actual zero triples |
| Original family has two free gaps | PASS | selected triple \(r,r+a,r+a+b\) |
| \(T^2\) outer factor | PASS | \(T^3\) cubic weight divided by \(T\) selected triples per heavy prime |
| \(H^4T^2=N^4/T^2\) ledger | PASS | direct from Smith mass \(H^4\) |
| Worst endpoint \(T=N^{3/5}\) for \(\mu=4\) | PASS | exponent decreases with \(T\) |
| Local affine Fourier identity | PASS | finite two-coordinate orthogonality |
| \(\sqrt{E_
u}\le M_
u\) | PASS | \(E_
u\le M_
u^2\) |
| \(M_
u\le H^{2+o(1)}\) | PASS, scoped | only when \(p,x,r\) are fixed and family is two-gap |
| Same \(M_
u\) bound for raw three-gap family | FAIL | raw family has \(H^3\) gap triples |
| \(\sum_
u w_
u\ll H^{1+o(1)}\) | FAIL | missing prime/root/level aggregate; no source theorem |
| \(\sum_
u w_
u M_
u\ll H^{3+o(1)}\) | FAIL | depends on preceding false premise |
| Source-valid formal aggregate | PASS at \(H^{4+o(1)}\) | Sylvester/Smith plus resultant height |
| Mean term retained | PARTIAL in Q5767 | omitted in headline; harmless after the \(H^4\) envelope |
| Degree-one restriction avoids extension roots | PASS | if \(x\in\mathbf F_p\) is literal in \(\Omega\) |
| Same root used in both affine coordinates | PASS in abstract setup | must retain ordered root/local-ring data |
| Vector origin equals Q5605 actual start | NOT ESTABLISHED | no canonical multiplicity-preserving consumer theorem |
| Projective variance supplies affine radial mass | FAIL | loses zero state and primary radius |
| New unconditional \(N^{12/5}\) theorem | FAIL | \(H^3\) aggregate is unproved |
| Canonical theorem insertion justified | NO | replacement ledger remains \(N^{14/5+o(1)}\) |

---

## 10. First false line and required withdrawals

The first exact false/unproved inequality in Q5767's exponent proof is

\[
oxed{\sum_{
u=(p,\mathfrak q,r,	au_
u)}w_
u\ll H^{1+o(1)}.}
\]

It silently treats the \(O(H)\) local degree/height of one system as the total number of weighted local cells across all \(H^2\) systems and the entire prime shell.

The following Q5767 claims must therefore be withdrawn or relabeled as conditional:

```text
“Q5739's Q5605 normalization has sum_nu w_nu << H^(1+o(1)).”

“Combining sqrt(E_nu)<=M_nu and M_nu<=H^2 proves
 sum_nu w_nu sqrt(E_nu)<<H^3.”

“The current canonical inputs unconditionally improve the Q5605
 N^(14/5) cubic line to N^(12/5).”
```

The local sector decompositions in Q5767—diagonal, reflection, overlapping windows, and bounded labels—remain useful cardinality observations. Their advertised global \(H^2\) weighted consequences also used the same unproved level-weight normalization and therefore should not be treated as established Q5605 exponents without a corrected aggregate weight theorem.

Q5739's conditional statement

```text
“in a shell with H^(1+o(1)) local levels ...”
```

need not be withdrawn; it must simply remain a hypothesis template rather than being cited as a proved fact about the literal Q5605 shell.

---

## 11. No canonical insertion theorem

Because the verdict is FAIL, there is no insertion-ready unconditional \(N^{12/5}\) theorem for `proof.tex`, `oracleA_result.tex`, or `atom_tail_section.tex`.

The strongest safe research-note replacement is:

```latex
\begin{remark}[Trivial affine-energy majorant does not improve the localized cubic ledger]
Let $\Omega_{\nu,H}$ be a source-aligned fully masked two-gap occurrence family, grouped by characteristic, rational residual root, and primary level, and let $M_\nu=|\Omega_{\nu,H}|$. Then $\sqrt{E_\nu}\le M_\nu$, whence
\[
 \sum_\nu w_\nu\sqrt{E_\nu}\le \sum_\nu w_\nu M_\nu.
\]
The current Sylvester--Smith/resultant-height argument gives only
\[
 \sum_\nu w_\nu M_\nu\ll H^{4+o(1)},
\]
not $H^{3+o(1)}$. Consequently this trivial affine majorant recovers $H^4T^2\ll N^{14/5+o(1)}$ on $N^{3/5}\le T\le N^{2/3}$. A bound $H^{3+o(1)}$ would be a new global prime/root/level sparsity theorem and would imply the conditional fallback $N^{12/5+o(1)}$.
\end{remark}
```

This text is a correction note, not a new theorem, and was not inserted into authoritative source.

---

## 12. Final conclusion

The affine-state framework itself survives the audit:

\[
n_
u(0)\lerac{M_
u}{|R_
u|^2}+\sqrt{E_
u}
\]

is the right noncircular local inequality, and it preserves the same residual root in both coordinates. What fails is Q5767's globalization of the trivial bound. The prime/root/primary-level cells are precisely the arithmetic mass that Q5605's Smith/resultant argument was designed to pay. They do not form an \(H^{1+o(1)}\)-weight index set merely because one recurrence transfer has length \(O(H)\).

The corrected unconditional endpoint is therefore

\[
oxed{\sum_{N<p\le2N}Z(p)^3\ll N^{14/5+o(1)}.}
\]

There is no new \(N^{-2/5}\) saving from the trivial affine collision bound. Obtaining \(N^{12/5+o(1)}\) requires a genuinely new \(H^3\) aggregate theorem plus an exact source-level consumer identification; obtaining the staged sigma-dependent targets requires still stronger, subquadratic primary-mass savings.