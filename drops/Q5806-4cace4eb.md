ANSWER Q5806 4cace4eb

# Problem 3.2 — cold audit of the claimed unconditional H³T² = N^(12/5) affine-mass bound

## Verdict

**FAIL.**

The local affine Fourier inequality used by Q5767 is correct, and the elementary bound

$$
√E_ν ≤ M_ν
$$

is correct. For the **two-gap** family that actually underlies Q5605's localized cubic-moment argument, the pointwise estimate

$$
M_ν ≤ H^{2+o(1)}
$$

is also legitimate, provided that ν fixes the characteristic, the rational residual root or local factor, the primary level, and every unbounded outer coordinate.

The claimed global consequence does not follow. Its first false numerical line is

$$
∑_ν w_ν ≪ H^{1+o(1)}.
\tag{F}
$$

In Q5767, ν = (p, 𝔮, r, τ_ν), with 𝔮 = X − x for a rational root and w_ν = log p at degree one. Thus ν ranges over the prime shell, rational residual roots, primary levels, and bounded chart data. None of the canonical sources proves that the total weight of those cells is H^(1+o(1)). The only source-valid unconditional aggregate supplied by the Q5605 Sylvester/Smith ledger is

$$
∑_ν w_ν M_ν ≪ H^{4+o(1)},
\tag{C1}
$$

not H^(3+o(1)). Consequently the trivial affine-energy majorant gives only

$$
∑_ν w_ν √E_ν ≤ ∑_ν w_ν M_ν ≪ H^{4+o(1)},
\tag{C2}
$$

which recovers Q5605's existing dyadic line

$$
H^4 T^2 = N^4 T^{−2} ≤ N^{14/5+o(1)}
$$

throughout N^(3/5) ≤ T ≤ N^(2/3).

The algebraic substitution

$$
H^3 T^2 = N^3/T ≤ N^{12/5}
$$

is arithmetically correct **conditional on a new aggregate theorem** of strength H^(3+o(1)). Q5767 did not prove that aggregate theorem; it replaced it by the incompatible pair of assertions “M_ν ≤ H²” and “∑_ν w_ν ≤ H.” Fixing the residual root makes the first assertion available but forces every prime/root/level cell into the second sum. Omitting the root from ν would avoid that exact sum only by putting the root multiplicity back into M_ν.

There is a second, logically earlier source gap. Q5674's literal distinguished-state carrier tests B_x = 0 at a rational common root. Q5739/Q5767 instead introduce an abstract two-coordinate terminal projection P_ω T_ω (1,5)ᵗ and test that vector for zero. No canonical Problem 3.2 source currently states the one-to-one, multiplicity-preserving theorem identifying these vector-zero occurrences with Q5605's selected actual Apéry triple starts. The affine origin inequality is valid as an abstract finite-group identity, but its exact downstream consumer binding is not yet a canonical theorem.

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

This audit uses only Problem 3.2 in the canonical repository and the exact connected same-project records named below. No workspace or Zinan repository was used.

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
  blob 05ad067fd11d75a51e456c85bad5d97f0c9514ee

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
  Notion page 3cb7a6c4-fa84-8131-a48b-f0b8792189f5
  actual distinguished-state carrier and primitive-return masks

Q5739 02a18711
  Google Doc 1MFvjs1pjm3WNloXASVnrqA9x6H078WJ8AJom9iZbxfI
  module A_Y/(B_Y,C_Y), local length, and abstract affine variance

Q5767 41375002
  Google Doc 1T8h-IdP234meXEztUXzklHYKW4pIPgRVYrUCzpekpg0
  claimed H^3 affine fallback and N^(12/5) consequence
```

Q5739 contains an important scope warning: it explicitly says that the literal Q5605 normalization was not exposed to that audit and that it was using the interface stated in its prompt. Its phrase “in a shell with H^(1+o(1)) local levels” is therefore a conditional template, not a theorem about Q5605's actual prime/root/level index set.

---

## 1. The literal Q5605 quantity

### 1.1 Dyadic prime and zero-count shells

For a fixed dyadic characteristic block, put

$$
P(N) = {p prime : N < p ≤ 2N},
$$

$$
Z_p = {0 ≤ x < p : b_x ≡ 0 mod p},
$$

and Z(p) = |Z_p|.

For a dyadic zero-count height T, Q5605 uses

$$
P_T = {p in P(N) : T < Z(p) ≤ 2T}
\tag{1.1}
$$

in the licensed high range

$$
N^{3/5} ≲ T ≲ N^{2/3}.
\tag{1.2}
$$

It chooses

$$
H = ceil(16N/T),
$$

so H is asymptotic to N/T. The factor 16 is irrelevant to exponents but ensures that the block loss below is a fixed fraction of T.

### 1.2 Selected actual zero triples

Partition [0,p−1] into consecutive intervals of length at most H. If one interval contains actual Apéry zeros

$$
z_1 < z_2 < ... < z_m,
$$

select its consecutive triples (z_i,z_(i+1),z_(i+2)). Their total number Q_p(H) satisfies

$$
Q_p(H) ≥ Z(p) − 2 ceil(p/H).
\tag{1.3}
$$

For p in P_T and the chosen H,

$$
Q_p(H) ≫ T.
\tag{1.4}
$$

Every selected triple has the unique representation

$$
r < r+a < r+a+b < p,
$$

where a,b ≥ 2 and a+b < H.

Define t_p^Ap(a,b) to be the number of selected starts r with these exact gaps. Then

$$
Q_p(H) = ∑_{a,b≥2, a+b<H} t_p^Ap(a,b).
\tag{1.5}
$$

This is a **two-gap** family. The start r is an occurrence/root variable, not a third gap parameter.

### 1.3 The positive primary mass

The logarithmically weighted mass is

$$
M_Ap(N,H) = ∑_{a,b≥2, a+b<H} ∑_{N<p≤2N} t_p^Ap(a,b) log p.
\tag{1.6}
$$

Equations (1.4) and (1.5) give

$$
T |P_T| log N ≪ M_Ap(N,H).
\tag{1.7}
$$

This is where the prime-shell weight is paid. The later unweighted cubic moment loses the single log N again; it does not create an additional prime sum.

### 1.4 The source-valid Sylvester/Smith envelope

A selected actual triple forces

$$
N_a(r) ≡ 0 mod p,
$$

and

$$
N_b(r+a) ≡ 0 mod p.
$$

Let

$$
S_(a,b) = Res_X(N_a(X), N_b(X+a)).
\tag{1.8}
$$

The canonical root-strip theorem makes S_(a,b) nonzero. Define t_p^form(a,b) as the number of rational common roots x in F_p of those two polynomials, with the relevant safe central and structural restrictions inserted when the deflated form is used. Then

$$
t_p^Ap(a,b) ≤ t_p^form(a,b) ≤ v_p(S_(a,b)).
\tag{1.9}
$$

The second inequality is the fixed-size Sylvester/Smith corank argument. It pays distinct common roots and repeated local multiplicity; no root-simplicity assumption is needed.

For a fixed gap pair,

$$
∑_{N<p≤2N} v_p(S_(a,b)) log p ≤ log |S_(a,b)|.
\tag{1.10}
$$

The canonical height estimate is

$$
log |S_(a,b)| ≪ (a+b)^2 log(a+b).
\tag{1.11}
$$

Consequently

$$
M_Ap(N,H) ≤ ∑_{a,b≥2, a+b<H} log |S_(a,b)| ≪ ∑_{s≤H} s^3 log s ≪ H^{4+o(1)}.
\tag{1.12}
$$

This H⁴ mass is what Q5605 actually banks unconditionally.

### 1.5 Why the cubic contribution has a factor T²

On P_T, Z(p)³ ≪ T³. From (1.7) and (1.12),

$$
|P_T| ≪ H^{4+o(1)} / (T log N).
\tag{1.13}
$$

Therefore

$$
∑_{p in P_T} Z(p)^3 ≪ T^3 |P_T| ≪ H^{4+o(1)} T^2.
\tag{1.14}
$$

The T² is exactly

```text
T^3 from the cubic zero count
minus
one factor T supplied by the selected-triple lower bound Q_p(H) >> T.
```

Using H asymptotic to N/T,

$$
H^4 T^2 = N^4 T^{−2}.
\tag{1.15}
$$

This decreases throughout the licensed T range, so the worst endpoint is T = N^(3/5):

$$
N^4 N^{−6/5} = N^{14/5}.
\tag{1.16}
$$

Thus the prior N^(14/5) term counted the heavy-prime contribution to the **third zero moment**, via actual consecutive zero triples and their two-gap Smith-primary envelope. It did not count a free raw three-gap family.

---

## 2. Reindexing the formal occurrences correctly

Let

$$
Y_H = {(a,b) : a,b ≥ 2 and a+b < H}.
$$

Then |Y_H| = H^(2+o(1)).

For Y = (a,b), a characteristic p, and a rational common root x, let e_(p,Y,x) be the local multiplicity retained after the literal safe masks. A fully expanded formal primary occurrence has the shape

$$
iota = (Y,p,x,r,tau),
$$

where 1 ≤ r ≤ e_(p,Y,x), and tau records bounded orientation, source, chart, and physical labels.

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

$$
s(iota) = P_iota T_iota (1,5)^t
\tag{2.1}
$$

be the proposed affine state.

To perform Fourier analysis, occurrences may be grouped only when their values lie in the same finite additive group. For degree-one factors a favorable grouping is by

$$
nu = (p,x,r,tau_nu),
$$

with local ring R_nu isomorphic to F_p[epsilon]/(epsilon^r). Put

$$
Omega_(nu,H) = {(Y,tau) : (Y,p,x,r,tau) is a formal occurrence in cell nu},
$$

and M_nu = |Omega_(nu,H)|.

At degree one the layer weight is, up to bounded normalization,

$$
w_nu = log p.
\tag{2.2}
$$

The exact incidence identity is

$$
∑_nu w_nu M_nu
= ∑_{Y in Y_H} ∑_{N<p≤2N} ∑_{x in F_p} ∑_{1≤r≤e_(p,Y,x)} ∑_tau (log p) 1_(all formal masks).
\tag{2.3}
$$

The right side is the total formal prime/root/level occurrence mass. Regrouping cannot change it.

If the shell equations in (2.3) are the same two common-root equations paid by the Q5605 Sylvester matrix, then Smith normal form gives the unconditional envelope

$$
∑_nu w_nu M_nu ≪ H^{4+o(1)}.
\tag{2.4}
$$

If the proposed terminal projection introduces a different local scheme not yet identified with that common-root scheme, then even (2.4) for the new object requires a bridge theorem. In neither interpretation does the source give H³.

---

## 3. Audit of the pointwise bounds

### 3.1 √E_nu ≤ M_nu: PASS

Let n_nu(z) count occurrences in Omega_(nu,H) whose state is z in R_nu². Define

$$
C_nu = ∑_z n_nu(z)^2,
$$

and

$$
E_nu = C_nu − M_nu^2 / |R_nu|^2.
$$

Then

$$
0 ≤ E_nu ≤ C_nu ≤ (∑_z n_nu(z))^2 = M_nu^2.
$$

Therefore

$$
√E_nu ≤ M_nu.
\tag{3.1}
$$

This step is unconditional once the finite occurrence set and common local ring are correctly defined.

### 3.2 M_nu ≤ H^(2+o(1)): PASS only in fixed-root two-gap scope

If nu fixes p,x,r and every unbounded non-gap variable, then the only free large parameters in Q5605 are (a,b) in Y_H. Hence

$$
M_nu ≤ |Y_H| H^{o(1)} ≤ H^{2+o(1)}.
\tag{3.2}
$$

The H^(o(1)) may absorb only labels already proved to have bounded or subpolynomial multiplicity. It may not absorb a new root, characteristic, host, or gap parameter.

### 3.3 Why the root cannot be omitted for free

If x is omitted from nu, then Omega_nu must include every rational root belonging to each gap pair. Since a gap polynomial has degree O(H), the elementary bound becomes

$$
M_nu ≤ H^{3+o(1)},
\tag{3.3}
$$

not H^(2+o(1)). More importantly, states at different roots live in differently based local rings until a canonical transport is specified.

There is an unavoidable bookkeeping dichotomy:

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

The canonical fully deflated four-return family has a,b,c ≥ 2 and a+b+c ≤ H, hence H^(3+o(1)) raw gap triples. The canonical aligned-pencil definition also retains three gap symbols; one gets only two effective parameters after a **specific proved alignment relation** fixes one combination.

Q5605's localized cubic consumer is the two-gap triple-start family, so (3.2) is correct for that consumer. It cannot be exported to the raw three-gap/four-return source. For the raw family the corresponding fixed-root trivial bound is

$$
M_nu ≤ H^{3+o(1)}.
\tag{3.4}
$$

Q5739's phrase “aligned three-gap system with two free parameters” is an abstract scope statement. It is not a canonical theorem supplying a deterministic c = c(a,b,...) relation for every raw three-gap consumer.

---

## 4. The fatal total-weight error

Q5767 writes

$$
∑_nu w_nu ≪ H^{1+o(1)},
$$

and

$$
M_nu ≤ H^{2+o(1)},
$$

then concludes

$$
∑_nu w_nu M_nu ≪ H^{3+o(1)}.
\tag{4.1}
$$

The last implication would be valid if both premises were proved. The first premise is not a Q5605 fact.

### 4.1 Literal expansion of the alleged level weight

For the degree-one actual-root decomposition,

$$
∑_nu w_nu
= ∑_{N<p≤2N} ∑_{x in F_p} ∑_{r≥1} ∑_{tau_nu: some Y is formal in (p,x,r,tau_nu)} log p.
\tag{4.2}
$$

This sum includes:

```text
all active characteristics p in (N,2N];
all active rational residual roots x;
all primary layers r;
all separately retained charts/branches.
```

There is no canonical estimate reducing (4.2) to the degree O(H) of one gap polynomial. Degree O(H) is a **per-system** statement; Q5605 has H² gap systems and an unrestricted dyadic prime shell.

Since every active cell has M_nu ≥ 1, one has merely

$$
∑_nu w_nu ≤ ∑_nu w_nu M_nu.
\tag{4.3}
$$

The source bounds the right side by H^(4+o(1)), not by H.

### 4.2 The Smith ledger shows exactly where the missing factors live

For one pair Y = (a,b), local prime/root/level multiplicity is paid by

$$
∑_{N<p≤2N} ∑_x ∑_{r≤e_(p,Y,x)} log p
≤ ∑_{N<p≤2N} v_p(S_(a,b)) log p
≤ log |S_(a,b)|.
\tag{4.4}
$$

The right side is H^(2+o(1)) in the worst pair range. Summing over the H^(2+o(1)) pairs gives H^(4+o(1)).

Q5767 effectively kept the H² pair count inside M_nu, kept one H as an alleged number of levels, and dropped the remaining prime/root/valuation height. That dropped factor is not a logarithm; in the source ledger it is a full power of H at aggregate scale.

### 4.3 Prime-shell and residue-degree audit

- **Prime shell.** The audit is already at one fixed dyadic N-shell, so summing over dyadic shells costs only N^(o(1)). But the primes *inside* (N,2N] are not one level. They are paid through (4.4).
- **Rational roots.** Each x in F_p is a separate degree-one local factor. Fixing x is what makes (3.2) true; it cannot then be omitted from (4.2).
- **Primary levels.** Summing r = 1,...,e_(p,Y,x) gives the local length. Smith valuation pays this sum. Replacing all levels by one cell would lose the exact primary multiplicity.
- **Residue degree.** Q5605's actual starts use degree-one factors. If higher-degree factors are included, w_nu = (deg 𝔮) log p and extension-root pollution is introduced; this cannot improve (4.2).
- **Charts and physical labels.** Source-proved bounded multiplicities may be absorbed in H^(o(1)). They do not repair the missing prime/root mass.

### 4.4 A legal interface obstruction

The claimed implication is refuted even at the abstract deterministic interface. Consider a permitted incidence profile with H^(4−o(1)) active cells nu, each having M_nu = 1 and |R_nu| tending to infinity. Then

$$
M_nu ≤ H^2,
$$

$$
E_nu = 1 − |R_nu|^{−2},
$$

and √E_nu = 1+o(1), but

$$
∑_nu w_nu √E_nu = H^{4−o(1)}
$$

when logarithmic weights are normalized by log N. This profile respects the Q5605 aggregate height ceiling and every pointwise inequality used by Q5767. It is not asserted to be an Apéry realization; it proves that the displayed scalar inputs do not imply the claimed H³ conclusion.

The canonical sources contain no covariance or overlap theorem forcing those active cells to collapse into only H^(1+o(1)) common local levels.

---

## 5. The affine origin identity and the downstream consumer

### 5.1 The finite-group identity: PASS

For a correctly defined finite occurrence set in one ring R_nu, two-dimensional character orthogonality gives

$$
E_nu = |R_nu|^{−2} ∑_{lambda≠0} |S_nu(lambda)|^2,
\tag{5.1}
$$

and

$$
n_nu(0) ≤ M_nu / |R_nu|^2 + √E_nu.
\tag{5.2}
$$

Defining Omega_(nu,H) **before** imposing s(omega) = 0 is the correct noncircular order. If Omega were pre-restricted to the distinguished origin, every surviving occurrence would lie in one state fibre and the energy statement would assume the conclusion.

Restricting x to F_p also correctly removes extension-root pollution.

### 5.2 The mean term must remain in the consumer

The aggregate local-length mass has the form

$$
M_aff(N,H) = ∑_nu w_nu n_nu(0).
\tag{5.3}
$$

Equation (5.2) gives

$$
M_aff(N,H)
≤ ∑_nu w_nu M_nu / |R_nu|^2
+ ∑_nu w_nu √E_nu.
\tag{5.4}
$$

Q5767's H³ headline discusses only the second term. Under the source-valid formal envelope (2.4), the mean is harmless but must be shown. Since |R_nu| ≥ p is asymptotic to N,

$$
∑_nu w_nu M_nu / |R_nu|^2
≤ N^{−2} ∑_nu w_nu M_nu
≪ N^{−2} H^{4+o(1)}.
\tag{5.5}
$$

After the outer T², this is

$$
N^{−2} H^4 T^2 = N^2 T^{−2} ≤ N^{4/5+o(1)},
\tag{5.6}
$$

so it is not the obstruction. The obstruction is the energy term's aggregate formal mass.

### 5.3 Exact consumer identification: not established in canonical source

Q5605's literal selected occurrence is an actual zero triple:

$$
b_x ≡ b_(x+a) ≡ b_(x+a+b) ≡ 0 mod p.
\tag{5.7}
$$

Q5674 refines a formal common root by retaining the distinguished initial slope. In its literal normalization, with a companion solution C_x, the state resolvent is

$$
Theta_(p,X)(Y,Z) = ∏_{x in X(F_p)} (C_x Y − B_x Z),
\tag{5.8}
$$

and the actual Apéry start is selected by B_x = 0. The companion coordinate is generally a unit rather than zero.

Q5739 uses different symbols B_Y,C_Y: they are asserted to be two **terminal residual coordinates after a projection** P_omega, and the desired event is

$$
P_omega T_omega (1,5)^t = (B_Y,C_Y)^t = 0.
\tag{5.9}
$$

These are not the same two coordinates as the fundamental Apéry/companion row in Q5674. That distinction is necessary—an invertible raw transfer cannot send the primitive initial vector to zero—but it creates a proof obligation.

The canonical repository currently contains no multiplicity-preserving theorem asserting

```text
Q5605 selected actual start (p,x,a,b)
if and only if
one fully masked Q5739 occurrence omega has s(omega)=0.
```

Q5739 itself records that it was not reading the literal Q5605 normalization.

Therefore (5.2) is a valid abstract reduction, but Q5767 does not yet establish that its left side is exactly the Q5605 mass. This is an additional reason that no canonical N^(12/5) theorem may be inserted.

---

## 6. Canonical variance inputs do not supply the missing level count

### 6.1 Projective variance

`projective_variance_reduction.tex` proves, for four-consecutive windows in one projective fibre, statements including M_p ≤ p when the first index determines the window, and the short-chain energy bound

$$
∑_{q in P^1(F_p)} C_p(q)^2 ≤ 30 p H.
\tag{6.1}
$$

It also proves a mass-sensitive variant

$$
∑_q C_p(q)^2 ≤ 20 p √M_p.
\tag{6.2}
$$

These estimates do not imply (F):

1. they are per-characteristic and retain a factor p;
2. they concern nonzero projective directions;
3. the distinguished zero vector has no projective image;
4. they forget the radial primary level r;
5. they concern four-return projective windows, not the complete two-gap affine-origin mass.

At the endpoint H asymptotic to p^(1/3), the canonical bound pH is H⁴, exactly the trivial square scale for an H²-sized family. It cannot be converted into an H-sized set of affine local levels.

### 6.2 Canonical continuant collision energy

`energy_result.tex` proves an unconditional formal-column bound

$$
E_p(H) ≪ H^{8/3}
$$

for a different energy built from pairs of zero levels N_h(x) = 0 in a continuant column. It explicitly distinguishes actual affine witnesses from raw resultant support and warns that boundary or resultant-supported rows need not have a common rational witness.

That theorem is not an affine-state collision estimate for P_omega T_omega (1,5)^t. It supplies neither (F) nor a primary-level/root count for the Q5605 mass.

### 6.3 Aligned three-gap content

`atom_tail_section.tex` defines the raw aligned pencil with three gaps and proves a valuation statement for common roots of three polynomials. Its desired aggregate H^(3+o(1)) aligned-content mass is stated as a hypothesis, not an unconditional theorem. Exact computations support the hypothesis at small heights but are not asymptotic proof.

Thus no canonical projective, continuant-energy, or aligned-content statement fills the missing line in Q5767.

---

## 7. Correct replacement ledger

### 7.1 Local affine statement

The following remains correct:

$$
n_nu(0) ≤ M_nu / |R_nu|^2 + √E_nu,
$$

and

$$
√E_nu ≤ M_nu.
\tag{7.1}
$$

### 7.2 Aggregate statement available from current source

Under the most favorable source-aligned interpretation of the shell equations,

$$
∑_nu w_nu M_nu ≪ H^{4+o(1)}.
\tag{7.2}
$$

Therefore

$$
∑_nu w_nu √E_nu ≪ H^{4+o(1)},
\tag{7.3}
$$

and

$$
M_aff(N,H) ≪ H^{4+o(1)}.
\tag{7.4}
$$

This is not a new estimate; it is the existing formal Smith envelope rewritten in affine-energy notation.

### 7.3 Downstream cubic moment

The correct unconditional Q5605 chain remains

$$
T |P_T| log N ≪ M_Ap(N,H) ≪ H^{4+o(1)},
\tag{7.5}
$$

and

$$
∑_{p in P_T} Z(p)^3 ≪ H^{4+o(1)} T^2 = N^4 T^{−2+o(1)} ≤ N^{14/5+o(1)}.
\tag{7.6}
$$

No fixed power below 14/5 follows from the trivial affine-energy inequality alone.

### 7.4 What would be needed for N^(12/5)

A genuinely new theorem of either of the following forms would suffice:

$$
∑_nu w_nu M_nu ≪ H^{3+o(1)},
\tag{7.7}
$$

or, more directly,

$$
M_Ap(N,H) ≪ H^{3+o(1)}.
\tag{7.8}
$$

Together with the exact consumer bridge, either gives

$$
∑_{p in P_T} Z(p)^3 ≪ H^3 T^2 = N^3/T ≤ N^{12/5+o(1)}.
\tag{7.9}
$$

But (7.7) is not a consequence of M_nu ≤ H²; it is itself a global prime/root/level sparsity theorem. Equation (7.8) is an actual-primary theorem of the same strength. Neither is currently canonical.

---

## 8. Full exponent ledger and sigma thresholds

Let a valid aggregate primary-mass estimate be

$$
M(N,H) ≪ H^{mu+o(1)}.
\tag{8.1}
$$

Then the Q5605 heavy-tail contribution is

$$
H^mu T^2 = N^mu T^{2−mu+o(1)}.
\tag{8.2}
$$

For T = N^t, with 3/5 ≤ t ≤ 2/3, its exponent is

$$
f_mu(t) = mu + (2−mu)t.
\tag{8.3}
$$

| Mass exponent | Worst endpoint | Cubic exponent | Status |
|---|---:|---:|---|
| mu = 4 | t = 3/5 | 14/5 | proved by current Smith ledger |
| mu = 3 | t = 3/5 | 12/5 | conditional; Q5767 did not prove it |
| mu = 2 | all t | 2 | desired H² primary-mass scale |
| mu = 2−delta | t = 2/3 | 2−delta/3 | requires a subquadratic mass theorem |

The staged cubic targets recorded in Q5605 are

$$
N^{29/15−6 sigma+o(1)}
\tag{8.4}
$$

and the stronger

$$
N^{2−6 sigma+o(1)}.
\tag{8.5}
$$

For a subquadratic mass H^(2−delta+o(1)), equation (8.2) reaches the first target exactly when

$$
delta ≥ 1/5 + 18 sigma,
\tag{8.6}
$$

and reaches the second exactly when

$$
delta ≥ 18 sigma.
\tag{8.7}
$$

For comparison, the incorrectly claimed N^(12/5) line would still exceed N^(29/15−6 sigma) by exponent 7/15+6 sigma, exceed N^(2−6 sigma) by 2/5+6 sigma, and exceed the RR-style target N^(11/5−2 sigma) by 1/5+2 sigma.

The source-valid N^(14/5) line exceeds those three targets by, respectively,

$$
13/15+6 sigma,
$$

$$
4/5+6 sigma,
$$

and

$$
3/5+2 sigma.
\tag{8.8}
$$

Thus even a future proof of the conditional N^(12/5) fallback would be progress, not full closure.

---

## 9. PASS/FAIL table

| Audit item | Verdict | Reason |
|---|---|---|
| Canonical repository and include order pinned | PASS | `proof.tex` and included Problem 3.2 files were read at `c5d932b...` |
| Q5605 prior quantity reconstructed | PASS | heavy-prime third moment via selected actual zero triples |
| Original family has two free gaps | PASS | selected triple r,r+a,r+a+b |
| T² outer factor | PASS | T³ cubic weight divided by T selected triples per heavy prime |
| H⁴T² = N⁴/T² ledger | PASS | direct from Smith mass H⁴ |
| Worst endpoint T = N^(3/5) for mu = 4 | PASS | exponent decreases with T |
| Local affine Fourier identity | PASS | finite two-coordinate orthogonality |
| √E_nu ≤ M_nu | PASS | E_nu ≤ M_nu² |
| M_nu ≤ H^(2+o(1)) | PASS, scoped | only when p,x,r are fixed and family is two-gap |
| Same M_nu bound for raw three-gap family | FAIL | raw family has H³ gap triples |
| ∑_nu w_nu ≪ H^(1+o(1)) | FAIL | missing prime/root/level aggregate; no source theorem |
| ∑_nu w_nu M_nu ≪ H^(3+o(1)) | FAIL | depends on preceding false premise |
| Source-valid formal aggregate | PASS at H^(4+o(1)) | Sylvester/Smith plus resultant height |
| Mean term retained | PARTIAL in Q5767 | omitted in headline; harmless after the H⁴ envelope |
| Degree-one restriction avoids extension roots | PASS | if x in F_p is literal in Omega |
| Same root used in both affine coordinates | PASS in abstract setup | must retain ordered root/local-ring data |
| Vector origin equals Q5605 actual start | NOT ESTABLISHED | no canonical multiplicity-preserving consumer theorem |
| Projective variance supplies affine radial mass | FAIL | loses zero state and primary radius |
| New unconditional N^(12/5) theorem | FAIL | H³ aggregate is unproved |
| Canonical theorem insertion justified | NO | replacement ledger remains N^(14/5+o(1)) |

---

## 10. First false line and required withdrawals

The first exact false or unproved inequality in Q5767's exponent proof is

$$
∑_{nu=(p,q,r,tau_nu)} w_nu ≪ H^{1+o(1)}.
$$

It silently treats the O(H) local degree or height of one system as the total number of weighted local cells across all H² systems and the entire prime shell.

The following Q5767 claims must therefore be withdrawn or relabeled as conditional:

```text
“Q5739's Q5605 normalization has sum_nu w_nu << H^(1+o(1)).”

“Combining sqrt(E_nu)<=M_nu and M_nu<=H^2 proves
 sum_nu w_nu sqrt(E_nu)<<H^3.”

“The current canonical inputs unconditionally improve the Q5605
 N^(14/5) cubic line to N^(12/5).”
```

The local sector decompositions in Q5767—diagonal, reflection, overlapping windows, and bounded labels—remain useful cardinality observations. Their advertised global H² weighted consequences also used the same unproved level-weight normalization and therefore should not be treated as established Q5605 exponents without a corrected aggregate weight theorem.

Q5739's conditional statement

```text
“in a shell with H^(1+o(1)) local levels ...”
```

need not be withdrawn; it must simply remain a hypothesis template rather than being cited as a proved fact about the literal Q5605 shell.

---

## 11. No canonical insertion theorem

Because the verdict is FAIL, there is no insertion-ready unconditional N^(12/5) theorem for `proof.tex`, `oracleA_result.tex`, or `atom_tail_section.tex`.

The strongest safe correction is the following research-note statement:

```text
Trivial affine-energy majorant.

For a source-aligned fully masked two-gap occurrence family, grouped by
characteristic, rational residual root, and primary level, let M_nu be
the cell cardinality. Then sqrt(E_nu) <= M_nu, so

    sum_nu w_nu sqrt(E_nu) <= sum_nu w_nu M_nu.

The current Sylvester–Smith/resultant-height argument gives only

    sum_nu w_nu M_nu << H^(4+o(1)),

not H^(3+o(1)). Hence the trivial affine majorant recovers
H^4 T^2 << N^(14/5+o(1)) on N^(3/5) <= T <= N^(2/3).
A bound H^(3+o(1)) would be a new global prime/root/level sparsity
theorem and would imply the conditional fallback N^(12/5+o(1)).
```

This is a correction note, not a new theorem, and was not inserted into authoritative source.

---

## 12. Final conclusion

The affine-state framework itself survives the audit:

$$
n_nu(0) ≤ M_nu / |R_nu|^2 + √E_nu
$$

is the right noncircular local inequality, and it preserves the same residual root in both coordinates. What fails is Q5767's globalization of the trivial bound. The prime/root/primary-level cells are precisely the arithmetic mass that Q5605's Smith/resultant argument was designed to pay. They do not form an H^(1+o(1))-weight index set merely because one recurrence transfer has length O(H).

The corrected unconditional endpoint is therefore

$$
∑_{N<p≤2N} Z(p)^3 ≪ N^{14/5+o(1)}.
$$

There is no new N^(−2/5) saving from the trivial affine collision bound. Obtaining N^(12/5+o(1)) requires a genuinely new H³ aggregate theorem plus an exact source-level consumer identification; obtaining the staged sigma-dependent targets requires still stronger, subquadratic primary-mass savings.