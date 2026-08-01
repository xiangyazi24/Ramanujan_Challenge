ANSWER Q6405 0ff4350f

# Lang–Trotter for defining-characteristic Mellin values: exact moments, norm obstruction, and the real gap

## Executive verdict

The sheaf-theoretic normalization is now strong enough to make the characteristic-zero moment problem completely legitimate, but it does not by itself cross the defining-characteristic local-limit barrier.

The strongest conclusions I can justify are:

1. [PROVED-PARSEVAL] There is an exact Parseval identity. Because the degree-two map

is the Fricke w_6 quotient of the Beauville-IV modular pencil, its two generic points are 6-isogenous and have the same Frobenius trace. This supplies an extra factor of two that is absent from the naive “(p-1)\sum_x a_x^4” estimate.

1. [PROVED-MOMENT] For all characters, including the trivial character,

For nontrivial characters only,

An exact finite-p cohomological formula is given below. The Sato–Tate fourth moment 2 is the input behind \sum_xa_{p,x}^4\sim2p^3; the Fricke double cover changes the Mellin second-moment constant from 2 to 4.

1. [NEGATIVE-NORM] The proposed cyclotomic norm argument is valid after one important correction: Galois conjugation sends the chosen prime \mathfrak p above p to another prime above p. Same-\mathfrak p zeros in one Galois orbit do imply divisibility of a norm by p^k, but the archimedean bound |S_r|\ll p^{3/2} yields only

which is weaker than the trivial orbit-size bound. Parseval does not improve this exponent. The norm method has a structural threshold mismatch: divisibility costs p, while a typical conjugate has size p^{3/2}.

1. [NEGATIVE-KATZ-LOCAL] Katz’s Mellin equidistribution is fully applicable to the characteristic-zero normalized Weil traces, once the certified monodromy hypotheses are inserted. It gives continuous Sato–Tate distribution and every fixed moment. It does not give equidistribution of the reductions in \mathbf F_p, nor exact zero counts. Even hypothetical individual square-root bounds for the relevant additive Fourier sums would give only O(p^{1/2}) zeros.

1. [GAP-LT-MELLIN] The clean deterministic new target is

A uniform Z_p=O(1) statement should be separated as the stronger [GAP-LT-MELLIN+]. Contrary to the wording “the truth is O(1),” an independent random/Poisson model predicts bounded moments and mean 1, but it predicts that the maximum over infinitely many primes is slowly unbounded. A uniform absolute bound would therefore reflect extra arithmetic rigidity, not ordinary Lang–Trotter randomness.

1. The best unconditional pointwise bound currently banked in the workspace remains

from the recurrence/continuant argument in working_notes/Q3.2_density_theorem.md. The new sheaf interpretation explains the expected scale and gives exact characteristic-zero moments, but it does not yet improve that exponent.

---

## 1. Setup and scope of the audit

I used the explicit state in the question as certified input:

- the Beauville-IV elliptic pencil E_x and its compatible rank-two system;

- H_p(x) as its Hasse scalar and H_p(x)^2\equiv a_{p,x}^2\pmod p;

- the integral rank-three companion \mathcal A_- of absolute conductor at most 11;

- the equality

I also checked the current transfer repository at commit 831e147, especially:

- working_notes/Q3.2_density_theorem.md, which contains the unconditional 3p^{2/3}+O(1) zero bound;

- working_notes/Q32_KatzMellin_vertical.md, which already records the essential warning that Sato–Tate equidistribution does not imply exact characteristic-p vanishing.

Write

and, for x\in U(\mathbf F_p),

For t\in\mathbf F_p^\times, define the pushforward coefficient

With \chi_r=\omega^{-r},

The index range 1\le r\le p-1 is one complete set of all N multiplicative characters: r=p-1 is the trivial character. When I say “nontrivial characters,” I mean 1\le r\le p-2.

---

## 2. The hidden factor two: \phi is the Fricke quotient

This is the cleanest new structural computation in the audit.

Zagier identifies the Franel sequence

with Beauville family IV

of modular group \Gamma_1(6); see Don Zagier, “Integral solutions of Apéry-like recurrence equations”, DOI 10.1090/crmp/047/22, PDF pp. 5 and 12–13, especially the Beauville table and the correspondence table.

For the project model

one has

Maier’s Hauptmodul x_6 on X_0(6) satisfies the Fricke formula

see Robert Maier, “Algebraic Hypergeometric Transformations of Modular Origin”, DOI 10.1090/S0002-9947-07-04128-1, Appendix A, PDF p. 24, formula (A.7), together with the modular interpretation on PDF p. 2. Comparing Maier’s and Zagier’s displayed j-formulas gives

Therefore w_6:x_6\mapsto72/x_6 becomes, in the project coordinate,

Moreover,

and the last expression is invariant under x_6\mapsto72/x_6. Hence

This is precisely the nontrivial deck involution of the quadratic map \phi.

For p>3 and away from the finitely many bad or branch fibers, w_6 replaces an elliptic curve with its degree-six quotient. Thus E_x and E_{\iota(x)} are isogenous over \mathbf F_p, and

This conclusion is not a heuristic use of Sato–Tate: it is exact isogeny invariance of the Frobenius polynomial.

The fixed points of \iota satisfy

so there are at most two branch points. Their contribution to any fourth-moment sum is O(p^2) because |a_{p,x}|\le2\sqrt p.

---

## 3. Exact Parseval and the fourth moment

### 3.1 Orthogonality

Multiplicative-character orthogonality gives the exact identity

Because a generic fiber is \{x,\iota(x)\} and the two traces agree, each generic fiber contributes

The sum of a^4 over the two points of that fiber is only 2a_{p,x}^4. Consequently

If a slightly different convention removes additional exceptional points, (3.2) acquires only an explicitly finite O(p^2) correction. In particular,

So the leading constant is twice the naive one-fiber fourth moment.

### 3.2 Exact cohomological fourth-moment identity

Let

be the rank-two elliptic local system. At x, let \alpha_x,\beta_x be its Frobenius eigenvalues, so

Set

Then, exactly,

and

For the non-isotrivial Beauville family, the geometric monodromy is the standard \mathrm{SL}_2 type; hence \operatorname{Sym}^2\mathcal V and \operatorname{Sym}^4\mathcal V have no geometrically constant subobject. The Grothendieck–Lefschetz trace formula therefore gives

Since |U(\mathbf F_p)|=p-3, (3.4) yields the exact finite-p formula

Deligne’s weight theorem gives O(p^{3/2}) for the first cohomological trace and O(p^{5/2}) for the second. Thus

The constant 2 is the fourth moment of the standard \mathrm{SU}(2) trace law. A primary source for the weight estimate is Pierre Deligne, “La conjecture de Weil II”, Publ. Math. IHÉS 52 (1980), 137–252, DOI 10.1007/BF02684780.

### 3.3 Final all-character and nontrivial-character moments

Combining (3.2) and (3.5) gives an exact expression:

Therefore

The trivial-character value is also explicit:

so

Subtracting it from (3.8) gives

These formulas validate the weight-three scale |M_p(r)|\asymp p^{3/2} for a generic nontrivial character.

### 3.4 What Parseval gives for S_p=M_p-T_p

Let B_p(t) be the coefficient trace function whose Mellin transform is T_p(r), and put

Then one still has the exact identity

For a pure semisimple compatible system of weight two, the leading constant in the right side is the dimension of the geometrically constant part of

equivalently, it is the sum of squares of the multiplicities of the geometrically irreducible constituents of \mathcal D. Hence

for a computable nonnegative integer c_{\mathcal D}, after the usual finite exceptional corrections.

The conductor bound \operatorname{cond}(\mathcal A_-)\le11 is enough for a uniform O(p^4) estimate, but not enough to identify c_{\mathcal D}. Computing that constant requires the actual semisimplified constituents and the mutual-twist calculation between the Fricke pushforward and \mathcal A_-. This is a finite monodromy calculation, not a local-limit theorem.

Most importantly, an upper moment does not upper-bound the number of divisible values. Zeros contribute nothing to (3.10). At best, a nonzero second moment plus the Weil bound proves that a positive proportion of the characteristic-zero values are nonzero. It does not distinguish units from multiples of \mathfrak p.

---

## 4. Cyclotomic fields, Galois orbits, and why the norm is too weak

Let r\bmod N have character order

Ignoring a fixed coefficient field for a moment,

If the companion has a fixed coefficient field E, replace K_m by the compositum EK_m and use a relative norm over E; none of the conclusions below changes except by constants depending on E.

### 4.1 Exact Galois action

For a\in(\mathbf Z/m\mathbf Z)^\times,

The primitive orbit has size

The genuine field norm is

If one instead multiplies over all a\in(\mathbf Z/N\mathbf Z)^\times, each primitive conjugate is repeated \varphi(N)/\varphi(m) times.

### 4.2 The prime-ideal correction

Because m\mid p-1, the rational prime p splits completely in K_m. If

then

This is generally a different prime above p. Thus one zero at the chosen reduction map does not force the whole Galois orbit to vanish at that same reduction map.

Conversely, if k distinct members S_p(ar) of the orbit vanish at the same chosen \mathfrak p, then S_p(r) is divisible by k distinct conjugate primes above p. For nonzero S_p(r) this gives the valid implication

This is the strongest correct form of the proposed orbit-product observation.

### 4.3 Archimedean size gives no nontrivial count

Suppose every conjugate satisfies

Then

Together with (4.3), this yields only

Since trivially k\le d, (4.4) says nothing.

Using AM–GM and a hypothetical optimal orbit-level Parseval estimate does not help:

The natural orbit average is of order p^3, so the right side is again p^{(3/2+o(1))d}.

This is not a technical loss. It is the core obstruction:

- one same-prime divisibility condition contributes a factor p to the norm;

- one typical archimedean conjugate already has size p^{3/2}.

A norm argument could become useful only after a genuinely new estimate such as

or an archimedean geometric mean <p^{1-o(1)}. Neither follows from Weil bounds or Sato–Tate.

Also, if S_p(r)=0 as an algebraic integer, its norm is zero and the norm method gives no information at all.

### 4.4 There is no verified “small value field” shortcut

For order m, the natural field is already the smaller cyclotomic field K_m, not the full K_N. A further strict subfield would mean that S_p(r) is fixed by a nontrivial subgroup of (\mathbf Z/m\mathbf Z)^\times, hence that several values S_p(ar) coincide exactly.

Such coincidences may occur at exceptional characters, but no theorem in the current state makes the generic degree bounded independently of p. On the contrary, the no-unbounded-Kummer-self-twist picture predicts that a generic primitive value has orbit degree close to \varphi(m).

For bounded m the degree is bounded, but the number of such character indices is itself small:

That separates low-order exceptional characters, but it gives no control over the bulk.

### 4.5 Resultant reformulation

Choose a primitive root g\bmod p and write t=g^j. Reduction at the chosen prime sends the Teichmüller value to g^j\in\mathbf F_p. Define the finite Fourier polynomial

Then

Consequently

Since p\nmid N, Y^N-1 is squarefree. On the order-m orbit, the corresponding condition is a common factor with the reduction of \Phi_m(Y); before reduction, the norm in (4.2) is a cyclotomic resultant.

This produces a precise intermediate problem:

[GAP-ORBIT-NORM] Prove a uniform bound, or even p^{o(1)} bound, for

v_poperatorname{Res}(P_p,Phi_m)

uniformly in m\mid p-1, after removing the certified bounded list of self-twist characters.

A bounded valuation would imply O(1) same-prime hits in each orbit, hence at best p^{o(1)} total after summing over the p^{o(1)} divisors of p-1. It would not by itself imply an absolute O(1) total.

---

## 5. The correct distribution heuristic

### 5.1 Uniform residue model

A random element of the residue field at a degree-one prime above p is zero with probability 1/p, regardless of the degree of the ambient number field. There are N=p-1 character indices, so the first-order prediction is

Without forced symmetries, the natural model is

If the exact reflection symmetry r\leftrightarrow N-r forces generic zeros to occur in pairs, remove the finitely many fixed characters and let Y_p count zero pairs. Then the sharper model is

It predicts

This is entirely consistent with mean 1.01 and maximum 8 below 25{,}000.

### 5.2 Important correction: Poisson does not predict a uniform bound

A Poisson or pair-Poisson sequence is tight for each prime and has bounded moments, but over infinitely many independent trials its maximum grows slowly. Heuristically, among primes up to X the maximum should be on the order of

up to the factor-two pairing convention.

Therefore:

- “mean approximately 1” and “usually at most a small constant” are the Lang–Trotter/random predictions;

- a theorem Z_p=O(1) for every prime would be stronger and would require special rigidity not present in the random model.

The data currently distinguish neither possibility.

### 5.3 Why archimedean Sato–Tate does not determine the factor 1/p

Katz’s Sato–Tate law controls

in complex embeddings. It describes the mass of values in fixed real or complex regions. The event

is a residue-class condition whose scale changes with p. It is not the event that the normalized complex trace lies in a small interval around zero.

This distinction is especially stark here because |S_p(r)| is naturally of size p^{3/2}. If S_p(r) were an ordinary integer, p\mid S_p(r) would allow roughly O(\sqrt p) possible multiples of p inside the Weil range. In the actual cyclotomic field there is not even a lower bound |S_p(r)|\ge p in one chosen complex embedding: divisibility is an ideal-theoretic condition, and units can redistribute archimedean sizes among conjugates.

For an elliptic curve, by contrast, |a_p|\le2\sqrt p<p for large p, so

That is why many “a_p\bmod p” questions for weight-two objects collapse to exact-trace Lang–Trotter. The present weight-three Mellin value does not enjoy that collapse. A closer classical analogue is the nonordinary-prime problem for higher-weight modular forms; see Fernando Gouvêa, “Non-Ordinary Primes: A Story”, Experimental Mathematics 6 (1997), 195–205.

---

## 6. What Katz proves—and exactly where it stops

Nicholas Katz’s primary reference is Convolution and Equidistribution: Sato–Tate Theorems for Finite-Field Mellin Transforms, Annals of Mathematics Studies 180, Princeton University Press, 2012, DOI 10.23943/princeton/9780691153308.001.0001.

The relevant statements are:

- Theorem 7.3 and Corollary 7.4, PDF p. 47: equidistribution of the normalized Mellin Frobenius conjugacy classes as the character varies;

- Remark 7.5, PDF pp. 47–48: effective O(q^{-1/2}) error for each fixed irreducible representation/test character;

- Theorem 28.1 and the estimates immediately following it, PDF pp. 178–179: uniform versions along sequences of finite fields, with constants growing with the tensor/representation degree.

Inserted into the now-certified fixed sheaves, these theorems legitimately give:

- the characteristic-zero Sato–Tate distribution of normalized M_p(r) or S_p(r);

- every fixed tensor moment;

- O(p^{-1/2}) discrepancy against a fixed continuous class function.

They do not give:

- distribution of the reductions \overline{S_p(r)}\in\mathbf F_p;

- the probability of one residue class of width 1/p;

- a test function whose complexity grows fast enough with p to isolate divisibility by \mathfrak p.

The constants in Katz’s fixed-moment estimates grow exponentially with tensor degree. Reaching a local scale of 1/p would require test complexity growing with p, outside the fixed-representation theorem.

Katz’s “Lang–Trotter Revisited” makes the same separation in another form. The mod-N_0 representations used there require N_0 prime to the characteristic; see PDF p. 15. The defining-characteristic ordinary unit root is treated separately. The function-field Lang–Trotter “Hope” statements on PDF pp. 17–18 concern exact characteristic-zero traces, not reduction of a weight-three trace modulo the same prime.

Katz’s “Wieferich Past and Future”, PDF p. 2, explicitly presents the 1/p random-residue heuristic and compares it with Lang–Trotter. That is heuristic evidence for (5.1), not a theorem deriving local residue statistics from Sato–Tate.

Finally, Katz’s “On a Question of Rudnick: Do We Have Square Root Cancellation for Error Terms in Moment Calculations?”, especially PDF pp. 1–3, emphasizes that even strengthening the standard O(q^{-1/2}) moment error is a separate and delicate question. This is much weaker than the local limit needed here.

---

## 7. The exact finite-level requirement and the square-root barrier

Define the residue-fiber counts

Then Z_p=N_p(0). For a nontrivial additive character e_p, Fourier inversion gives

This formula quantifies the missing theorem.

Even if one could construct a bounded-conductor sheaf in the r-variable and prove the individual Weil bound

uniformly in u\ne0, (7.1) would give only

Thus “one more application of Deligne” is not enough. To get Z_p=O(1) from (7.1), one needs the averaged cancellation

but (7.3) is essentially equivalent to the target itself.

There is also a normalization trap in the phrase “equidistributed with error o(1).” Since

is the probability of zero, a bound N_p(0)=O(1) requires probability O(1/p). Ordinary weak equidistribution with error o(1) is far too weak. Even an error O(p^{-1/2}) in probability permits O(p^{1/2}) points in one fiber. An asymptotic N_p(0)=1+o(1) would require an error o(1/p) in probability and, because N_p(0) is integral and often paired, would be unrealistically rigid.

A more realistic local-limit statement is distributional across primes, or Poisson factorial moments after averaging over primes.

---

## 8. Comparison with Fermat-quotient and additive-combinatorial results

The closest cited literature does not provide a general fixed-conductor theorem of the needed type.

### 8.1 Fermat quotients

For

one has the exceptional homomorphism

Equivalently, additive phases of q_p become multiplicative characters modulo p^2. The main estimates exploit precisely this p^2-level group structure.

Primary comparisons:

- Jean Bourgain, Kevin Ford, Sergei Konyagin and Igor Shparlinski, “On the Divisibility of Fermat Quotients”, Michigan Math. J. 59 (2010), 313–328. Its theorems use smooth numbers, multiplicative subgroups in residue rings, and Heilbronn sums to find small nonvanishing Fermat quotients; they do not prove singleton-fiber equidistribution for a generic trace-family Fourier transform.

- Igor Shparlinski, “Fermat Quotients: Exponential Sums, Value Set and Primitive Roots”, Bull. London Math. Soc. 43 (2011), 1228–1238, DOI 10.1112/blms/bdr058. The introduction states that Heath-Brown’s fixed-p estimate becomes nontrivial for intervals of length at least p^{1/2+\varepsilon}; Shparlinski reaches p^\varepsilon only on average over p. Lemma 5 is where the special character-mod-p^2 representation enters, and Theorem 8 is an average-over-primes estimate.

Our Mellin trace has no known analogue of q_p(uv)=q_p(u)+q_p(v).

### 8.2 Bourgain–Garaev–Konyagin–Shparlinski methods

Sum-product and subgroup methods can prove strong cancellation when the input has explicit additive/multiplicative growth, or when variables range over intervals/subgroups whose energy can be controlled. They do not imply that the defining-characteristic reductions of an arbitrary bounded-conductor sheaf Mellin transform have singleton fibers of bounded size.

The relevant incidence estimates, such as those for

are powerful because the phase has an explicit low-complexity product structure. For the present r\mapsto\overline{S_p(r)}, the exponent r ranges through all p-1 tame characters and the output is itself a complete global sum. No reduction to the standard BGKS setup is currently certified.

### 8.3 Targeted literature conclusion

After targeted searches through August 1, 2026, I did not find a theorem of the following general form:

for every fixed bounded-conductor compatible system \mathcal F, the defining-characteristic reductions of its complete Mellin values are equidistributed in \mathbf F_p as the multiplicative character varies, with O(1) or p^{o(1)} singleton fibers.

Such a statement would in fact be false without strong nondegeneracy assumptions: constant or Kummer constituents can make almost all Mellin values vanish. Katz’s Mellin-nondegeneracy and monodromy hypotheses repair the characteristic-zero Sato–Tate problem, but no published implication from them to a defining-characteristic local limit was located.

---

## 9. What a proof would actually require

### 9.1 Route A: a defining-characteristic Mellin local-limit theorem

One would need more than an \ell-adic compatible system. A plausible package would be:

1. an integral overconvergent F-isocrystal or Dwork module realizing all reductions \overline{S_p(r)};

1. a genuine parameter object for the character index r whose rank and conductor remain bounded as p grows;

1. exclusion of unit-root/Kummer degeneracies and all unbounded self-twists;

1. a local-limit or large-sieve theorem at residue scale 1/p, not merely a purity theorem.

Kiran Kedlaya’s “Fourier Transforms and p-adic ‘Weil II’”, Compositio Math. 142 (2006), proves purity in rigid cohomology for overconvergent F-isocrystals. That supplies a p-adic analogue of weight control; it does not supply the fourth item.

There is also a geometry problem: the set of all characters of \mathbf F_p^\times is the degree-(p-1) group scheme \mu_{p-1}. Its size grows with p. Katz’s Tannakian Mellin formalism handles this elegantly for characteristic-zero Frobenius classes, but it does not automatically produce a fixed finite-type, bounded-complexity r-parameter space in defining characteristic.

### 9.2 Route B: growing moments and incidence varieties

For fixed k, character orthogonality expands

into a weighted count on the incidence variety

For k=2, this is exactly the four-variable variety mentioned in the question. Its top-dimensional components include:

- permutation diagonals;

- Fricke pairings x\leftrightarrow\iota(x);

- components forced by any bounded Kummer self-twist;

- mixed components involving the companion system.

Classifying these components and applying Deligne off the diagonals would rigorously compute the fourth moment and the characteristic-zero collision structure. That is worthwhile and finite.

However:

- the four-variable variety controls a fourth archimedean moment, not divisibility by \mathfrak p;

- a singleton local limit normally requires moments of order growing with p;

- to approximate residue-scale delta functions one expects at least k\asymp\log p/\log\log p;

- one would need Betti-number/conductor bounds uniform in growing k, while Katz’s constants already grow exponentially with tensor degree;

- an additional p-adic or additive-character input is still needed to detect \mathfrak p\mid S_p(r).

Thus a bounded cohomology dimension for the k=2 incidence variety is not enough for [GAP-LT-MELLIN]. The genuinely useful growing-moment theorem would have to classify all top-dimensional components uniformly in k and keep the off-diagonal cohomology under control.

### 9.3 Route C: cyclotomic resultants and Cartier slopes

The exact reformulation (4.5) suggests a more arithmetic route:

1. construct an integral Fourier polynomial P_p before reduction;

1. identify the p-adic slopes of the resultants \operatorname{Res}(P_p,\Phi_m) through Cartier/Dwork theory;

1. prove that each primitive character orbit contributes only O(1) chosen-prime zeros;

1. combine with a bound on the number of exceptional orders or self-twist classes.

This is the most concrete way the norm idea could become nontrivial. It asks for a valuation theorem, not a size theorem.

### 9.4 Route D: exploit the Apéry recurrence directly

The current 3p^{2/3}+O(1) theorem already exploits special structure absent from generic trace functions: the order-two Apéry recurrence, nonconsecutive zeros, and polynomial continuants. Any pointwise improvement obtained from additional multi-gap constraints would be unconditional and could precede the full local-limit theory.

The sheaf picture can help identify which gap configurations are genuinely exceptional, but a recurrence/incidence hybrid may be more realistic than a universal defining-characteristic equidistribution theorem.

---

## 10. Clean formal gap statements

### [GAP-LT-MELLIN] — robust deterministic target

For every \varepsilon>0, uniformly over good primes,

Equivalent strengthened form:

This is already enough to beat the current p^{2/3} bound decisively and is compatible with random maxima.

### [GAP-LT-MELLIN-AVG] — realistic Lang–Trotter/local-limit theorem

After removing and recording the finitely many reflection-fixed characters, let Y_p be the number of generic zero pairs. Prove, for every fixed k,

This is the factorial-moment formulation of

and gives mean Z_p\to1 after averaging over primes. It is closer to the actual random/Lang–Trotter prediction than a uniform bound.

### [GAP-LT-MELLIN+] — extra-rigid pointwise conjecture

There is an absolute constant C such that

This is consistent with current computations but is not predicted by an independent Poisson model. It should be advertised as a stronger Apéry-specific rigidity conjecture.

### [GAP-ORBIT-NORM] — concrete intermediate arithmetic target

For every primitive order m\mid p-1, outside a bounded exceptional self-twist list,

This would give bounded same-prime hits per Galois orbit and a total p^{o(1)} bound after summing over orders.

---

## 11. Bankable conclusions for the project ledger

- [VERIFIED-THEOREM] \phi is the Fricke w_6 quotient; its deck partners are 6-isogenous and have equal a_{p,x}.

- [VERIFIED-THEOREM] Exact Parseval is (3.1), and the Fricke symmetry gives (3.2).

- [VERIFIED-THEOREM] The exact fourth-moment identity is (3.5).

- [VERIFIED-THEOREM] All-character second moment: 4p^4+O(p^{7/2}).

- [VERIFIED-THEOREM] Nontrivial-character second moment: 3p^4+O(p^{7/2}).

- [NEGATIVE-NORM] Norm plus Weil size plus Parseval gives no nontrivial zero-count bound.

- [NEGATIVE-KATZ-LOCAL] Katz gives characteristic-zero Sato–Tate and fixed moments, not defining-characteristic singleton fibers.

- [NEGATIVE-SQRT] Even an individual O(\sqrt p) additive-Fourier estimate yields only Z_p=1+O(\sqrt p).

- [OPEN] Compute c_{\mathcal D} for S=M-T from the semisimplified fixed pair and mutual twists.

- [GAP-LT-MELLIN] Prove Z_p\ll_\varepsilon p^\varepsilon.

- [GAP-LT-MELLIN-AVG] Prove the pair-Poisson factorial moments across primes.

- [GAP-LT-MELLIN+] Uniform O(1) requires extra Apéry-specific rigidity.

## Least-confident step

The least-confident part is not any of the moment or norm calculations; it is the proposed Route A. I do not know that the full character set r\bmod p-1 can be organized into a fixed-rank, fixed-conductor overconvergent F-isocrystal family whose trace is \overline{S_p(r)}. The natural character scheme \mu_{p-1} has degree growing with p, and even a successful construction would still need a new local-limit theorem beyond p-adic Weil II. The recurrence/resultant routes therefore look more concrete than invoking a hypothetical “vertical Weil II” alone.