# Q6413 FULL TEXT (retrieved from Notion page 3af7a6c4fa84817e893aeb50d4ae2b0f, 2026-08-01 08:10)

ANSWER Q6413 5f9bc927
R4-b — two-sheaf Mellin inverse program and defining-characteristic transfer audit

## Executive verdict

The audit objection is correct, and the previous FFK claim must be withdrawn more strongly than stated in the question.
Forey–Fresán–Kowalski (FFK) prove that, on a one-dimensional torus, only O(1) characters are cohomologically nongeneric for a fixed perverse sheaf. They do not prove that only O(1) generic Mellin fibers have Frobenius trace exactly zero. They certainly do not control reduction of a nonzero Mellin value modulo a prime above the defining characteristic.
I found no published theorem that, for a fixed bounded-complexity sheaf and one fixed prime p, bounds the number of tame characters with trace in a chosen prime P | p.
The Apéry/Franel geometry does have both étale and crystalline/rigid realizations at good primes. However, p | b_r is naturally a trace-zero congruence for the Frobenius on the r-twisted Mellin cohomology, not the assertion that a fixed Cartier operator is noninvertible. Packaging all r requires p-1 tame blocks, or a Kummer cover whose degree grows with p.
An exact-vanishing inverse theorem would still leave a separate defining-characteristic residue theorem. The exact-zero set and the nonzero-but-P-divisible set are disjoint and require different inputs.
Current general technology gives no bankable p^{1/2}, p/log p, or even p^{2/3-delta} count for exact equality of two Mellin transforms in this horizontal character aspect. Pairwise Deligne correlation, even if an exact-zero incidence sheaf were available, naturally stops at p^{3/4}. A square-root count would require a genuinely bounded-complexity zero detector; no such detector is presently known.
The Franel local system is rank 2, not rank 3. The primitive rank-3 object is Sym^2 of the Franel elliptic system, equivalently the Apéry rank-3 system after descent. The notation 3F2(1/3,2/3,1;1,1) has a cancelling numerator/denominator parameter and represents the same scalar series as 2F1(1/3,2/3;1); it is not a canonical irreducible rank-3 Franel sheaf.
There is also a necessary convention correction. Let pi be the quadratic Apéry–Franel cover, and let K_q be its sign local system. If t_! means the full degree-two direct image pi_!, then the projection formula gives pi_!(S) ⊗ K_q ≅ pi_!(S) because pi^* K_q is trivial. Thus the literal difference of the two Mellin transforms would be identically zero. Consequently, the machine-verified two-sheaf identity must use one eigendescent, a projector, or a virtual constituent — not the unprojected full pushforward. Below the primitive constituent is A_+ and its companion A_-.
This notation should be fixed before any monodromy or conductor constant is banked.

## I. Audit objection

### I.1 What FFK actually controls

FFK's generic-vanishing theorem controls the exceptional set on which the cohomology of M ⊗ L_chi is not concentrated in its generic degree; on G_m this is O_M(1) characters. That statement concerns cohomological jumps, not the value of the Frobenius trace on the generic cohomology group. For a generic character one can have a fixed-dimensional nonzero vector space and nevertheless a zero Frobenius trace. The three conflated notions:
1. Cohomological exception (FFK controls; O(1) characters on G_m).
2. Exact scalar vanishing (FFK does not count).
3. Defining-characteristic vanishing: a generally nonzero algebraic integer trace lies in P | p (FFK does not address).
FFK's equidistribution tests normalized Mellin traces against bounded continuous functions; an exact point mass at zero is not such a test. Their general horizontal quantitative theory requires an unestablished quantitative stratified-vanishing theorem.
The earlier FFK-based claim is RETRACTED; the mod-p inference was a second, independent error.

### I.2 Published defining-characteristic zero counts: verdict

No published general theorem controls #{tame chi : trace in P_p} for one fixed bounded-conductor compatible sheaf as the defining prime p varies. Open in this generality. Neighboring but insufficient: Perret-Gentil's large sieve (auxiliary coefficient primes lambda != p; algebraic parameter, not the character set); T-adic / slope theories (Liu–Wan, Davis–Wan–Xiao, Liu–Wan–Xiao, Ren–Wan–Xiao–Yu — vary wild/additive character or weight inside ONE tame component; Teichmüller exponent held fixed); Gross–Koblitz/Stickelberger (special hypergeometric cases only); Hasse–Witt ordinarity (fixed algebraic Hasse polynomial in a connected parameter, not per-tame-component traces).
Archimedean size O(p^{3/2}) does not help: divisibility by one prime above p is one local condition; in Q(mu_{p-1}) the prime p splits completely, so a P-condition imposes ONE coordinate, not all conjugates.
Exact equality is stable under the full Galois orbit of a character; vanishing modulo one P is stable only under the decomposition group of P (which can be trivial in the cyclotomic part). This is the exact-to-mod-p divorce.

### I.3 Crystalline companion, Cartier, and Z_p

Rational crystalline/rigid companions exist (Kedlaya). Not automatic: a canonical Frobenius-stable integral lattice whose reduction mod p matches the truncation coefficients. Each tame twist chi = omega^{-r} has a rank-one unit-root overconvergent F-isocrystal; a rigid Lefschetz trace formula expresses the Mellin sum as an alternating Frobenius trace, endpoint terms to be fixed. The missing item: integral comparison b_r ≡ (trace on integral rigid Mellin cohomology) mod P for every generic r, with torsion control — a concrete comparison theorem still to be printed.
Z_p is NOT literally a Cartier kernel: the condition is trace ≡ 0, not Frobenius noninvertible; generic Mellin cohomology of A_+ has dimension 2, so the one-dimensional simplification is unavailable generically. Packaging all p-1 characters needs rank ~p objects or a growing Kummer cover — the moving-complexity barrier.
Quarter-point formula (p ≡ 1 mod 24, m=(p-1)/4): naturally the explicit description of the order-four tame component (2A = Gauss/CM factor, U_m = residual universal recurrence). Not a template for all r.
Templates and stopping points: Adolphson–Sperber (p-adic Hasse–Witt refinement), Beukers–Vlasenko Dwork Crystals I–III (closest computational template; no zero-block counting theorem), Katz/FFK (no defining-char lattice), Kedlaya (rational polynomials only).
Verdict: constructing the integral trace comparison is feasible; counting its zero residues is the genuinely new theorem.

### I.4 The exact-to-mod-p interface lemma (DRS)/(RLL)

Fix compatible systems P_1, P_2 with integral structures; for good p choose P_p | p; D_p(chi) = difference of the two Mellin traces. Assume: outside E_p of size O_c(1), integral comparison holds; torsion-free in range; endpoint corrections in E_p; complexity bounded by c. Decompose the zero set into the exact part and the residual part. Required additional theorem, either:
- (DRS) direct anti-concentration: #{chi : D_p(chi) != 0 exactly but D_p(chi) in P_p} << p^alpha, or
- (RLL) the stronger residue local limit.
Thresholds: alpha < 2/3 improves the current Apéry fiber exponent; alpha < 1/2 is the k=2 interface threshold; alpha = 0 gives bounded row count.
(DRS) is logically independent of the exact inverse theorem; an exact proof likely uses the full Galois orbit, which the local event D_p(chi) in P_p destroys (Galois moves P_p).

## II. Exact two-sheaf inverse problem: correct formulation

Bounded conductor at one prime alone is FALSE without a long exception list (Kummer trace functions have Mellin transform supported at one character). Realistic strong target:

**Conjectural Mellin inverse theorem MI(c, eps):** for every complexity bound c and eps > 0 there is p_0(c,eps) such that: P_1, P_2 specializations at p > p_0 of fixed compatible systems of complexity <= c, no negligible/punctual/Kummer constituents, restricted to characters of order >= p^eta; if #{chi : S_1(chi) = S_2(chi)} >= p^eps then the pair lies in the finite geometric exception list: common constituent, Kummer/duality/automorphism relation, power-map induction, or graph relation between the two Tannakian representations.
The fixed-constant version (>= C matches forces relation) is much stronger; no current theorem implies it.
Rigorous endpoint: if D_p(chi) = 0 for all but at most C characters, Fourier inversion makes the difference trace function a sum of <= C Kummer characters; semisimplification + Chebotarev put the pair in the Kummer class. The hard range is p^eps zeros.

## III. Tannakian formulation

Work in Katz/Gabber–Loeser convolution category (FFK quantitative form). Generic chi gives fiber functor via H_c^1; arithmetic Frobenius gives conjugacy class Fr_{p,chi} in the arithmetic Tannaka group; traces recover S_i(chi). For the pair: joint group H <= G_1 x G_2 with the two defining representations; exact equality = trace-equality hypersurface F = 0. If H = G_1 x G_2 and F != 0, zero density follows conceptually — but NO lattice-scale count.
Steps where Goursat–Kolchin–Ribet, Larsen–Pink, Sawin enter: compute connected derived groups; Goursat for proper subgroups projecting onto both factors (must come from isomorphism of common simple quotients); translate graphs into sheaf relations (isomorphism, duality, finite-order twist, automorphism pullback, power-map induction); Larsen–Pink invariant tensors for lambda-independence; quantitative sheaf theory (Sawin/FFKS) for conductor control; FKMS robust Goursat + Xu stratification for multi-correlation diagonal classification.
MISSING: no theorem converts a sparse exact incidence set (p^eps Frobenius classes on the hypersurface) into a tensor invariant. That incidence-to-invariant step IS the new inverse theorem.
Exact arithmetic amplifier: one high-order exact zero propagates over the full Galois orbit (index <= [E:Q]) — a large cyclotomic packet. Full-orbit vanishing = vanishing of the primitive Fourier projector; via Ramanujan sums this is an alternating combination of subgroup averages = direct image/pullback under power maps. Program: prove a bounded-conductor noninduced sheaf cannot be annihilated by a high-order primitive projector; classify failures as power-map induction or finite Kummer packet. Limitation: power-map degree/subgroup index may grow with m.

## IV. Bilinear forms and amplifiers

Moment unfolding on the constraint torus; Deligne square-root cancellation when no geometrically trivial constituent; multiplicity m of the trivial rep controls the main term. Four amplifiers:
A. Shift amplifier: pairwise bilinear alone gives |S| << p^{3/4} — does NOT beat 2/3. Higher intersections or spectral large sieve essential.
B. Galois-orbit/primitive-projector amplifier: strongest exact-regime tool (see III).
C. Auxiliary-prime sieve (Perret-Gentil model): needs integral Mellin monodromy + strong approximation + finite-level equidistribution; plausible p^{1-delta} with small delta; research program, not off-the-shelf.
D. Polynomial/small-ball amplifier: degree must grow with p; conductor constants blow up; at best logarithmic density under unproved exponential moment bounds. Do not bank.
Deligne supplies p^{1/2} only AFTER a bounded-complexity zero detector exists; none is known. alpha < 1/2 is stronger than the one-shot Deligne barrier.

## V. Diagonal and exceptional pairs (THE 10-ITEM EXCEPTION LIST)

1. Common semisimple constituents (cancel from the difference).
2. Punctual and negligible objects (zero generic Mellin fiber / finite Fourier support).
3. Kummer constituents and finite Kummer packets (Mellin support at finitely many characters).
4. Constant and Tate lines (endpoint and trivial-character corrections).
5. Finite self-twists: P ≅ P ⊗ L_eta forces periodicity in character space.
6. Group automorphisms: pullback by t -> at, inversion t -> a/t, compositions (relabel characters, may combine with duality).
7. Duality/conjugate duality: P_2 ≅ D(P_1) ⊗ L_eta forces real or paired traces.
8. Power-map induction: induced from t -> t^d has Mellin spectrum supported on character congruence classes.
9. Graph subgroups of the joint monodromy (common simple quotient of G_1, G_2 gives nonproduct joint group).
10. Disconnected arithmetic components (equality identically true on one component while a hypersurface on the identity component).

### V.1 The Apéry quadratic companion is an exceptional graph case

The pair (A_+, A_- = A_+ ⊗ K_q as trace functions via chi_2(q(t)) twist) must be treated separately: K_q is the sign local system of y^2 = q(t), NOT a multiplicative Kummer sheaf (q(t) not monomial-times-square), so the twist cannot be absorbed into chi -> chi eta. Both systems come from the same descended rank-3 source, so product monodromy is false or unavailable without further calculation. Required: the joint Mellin Tannaka group and the deck-involution action on the two generic Mellin cohomology spaces. First nontrivial test case for the two-sheaf inverse theorem.

### V.2 Full-pushforward convention check

pi_! Sym^2 decomposes as A_+ ⊕ A_-; twisting the full direct image by K_q swaps the summands; a nonzero difference must use a projector/eigendescent/virtual combination. Record the projector explicitly or all rank/monodromy/constant claims are ambiguous by a factor of two.

## VI. Explicit current bounds (calibration)

Neither p^{1/2} nor p/log p is proved for exact two-transform equality at one prime. Only automatic: trivial p-1 with O_c(1) removable (cohomological). Under horizontal equidistribution assumptions: qualitative o(p), no rate competitive with the project's p^{2/3} fiber theorem. An exact inverse theorem by itself would NOT improve |Z_p| (addresses Z_p^exact only; actual zeros dominated by Z_p^res interface).

## VII. Exact Franel object, rank, conductor

Franel system rank 2 (elliptic Picard–Fuchs; CFVZ gauge of 2F1(1/3,2/3;1)); Sym^2 rank 3 primitive; A_± rank-3 descents; full pi_!Sym^2 rank 6; pi_!(E⊗E) rank 8 before removing determinant/Tate. Mod p the determinant line contributes a Tate trace divisible by p (why Hasse polynomial sees rank 3).
Local monodromy (p > 3): E lisse away from {0, 1/8?, -1, infinity} (four points on x-line); Sym^2 tame J_3 there. On t-line with alpha, beta roots of q(t) = t^2 - 34t + 1: A_+ has J_3 at 0, infinity and diag(1,1,-1) at alpha, beta; A_- (K_q twist) has diag(-1,-1,1) at alpha, beta. Geometric identity component SO_3; A_+ = disconnected O_3 with det = K_q; A_- connected.
Conductor conventions (two mixed in prior rounds): FKM analytic conductor vs Artin/drop complexity. Artin: cond(E) = 6, cond(Sym^2 E) = 11, cond(A_+) = 9, cond(A_-) = 11. The quoted 11 is a safe coarse/Artin bound, NOT the FKM analytic conductor — every quantitative argument must state its convention. Generic Mellin cohomology dimension of A_+ is 2. Weight-2 input => unnormalized Mellin traces weight 3, size O(p^{3/2}).

## VIII. Dispatchable research program (Packages A–F)

A — notation + integral comparison: replace ambiguous pushforward by A_±; print sign/Tate twist/endpoint/unit conventions; construct Frobenius-stable integral lattice in the rigid realization; prove b_r ≡ trace mod P; torsion-freeness; list nongeneric characters. (Dwork-crystal machinery, not equidistribution.)
B — joint Mellin monodromy: convolution Tannaka groups of A_± separately; joint group; deck involution action on generic Mellin cohomology; product vs graph vs disconnected extension; classify all Kummer/inversion/duality/power-map self-relations. (= [GAP-2] formal construction sheet.)
C — exact cyclotomic inverse theorem: Galois-orbit stability; high-order zero orbit => primitive-projector vanishing; Ramanujan-sum expansion; prove bounded-conductor noninduced sheaf cannot be annihilated (else power-map induced / Kummer constituent). Cleanest exact route.
D — finite-level auxiliary-prime sieve: integral representations mod lambda != p; strong approximation; finite-level equidistribution; large sieve on the trace-equality hypersurface. Maybe first explicit p^{1-delta}, likely not competitive with 2/3.
E — defining-characteristic residue local limit: prove (DRS) or (RLL). The ONLY package touching |Z_p|. Candidate inputs: integral Cartier matrices for tame-twisted rigid complexes; cross-component complexity theorem; explicit Gross–Koblitz in hypergeometric specializations; a new p-adic large sieve / local limit over the finite étale tame-character scheme. NO existing reference completes this.
F — special order-four component: quarter-point identity isolates the order-4 block (CM/Gauss factor + universal recurrence U_m as Frobenius entry of a rank-2 crystal). Tractable test of A; do not extrapolate to all r. [Our Q6360 already delivered this — see appendix X / FABLE_NOTES.]

## IX. References checked

FFK Arithmetic Fourier transforms (Thm 1.1 generic/stratified vanishing, Thm 1.3 + §5 equidistribution, Rmk 1.4(3) horizontal caveat, §§9–11); Katz, Convolution and Equidistribution (PUP 2012); Perret-Gentil IMRN 2020 (DOI 10.1093/imrn/rny202, Prop 4.1, §4) and Trans. AMS 371 (2019) (DOI 10.1090/tran/7333); Fouvry–Kowalski–Michel–Sawin bilinear forms with trace functions (2025/26 preprint; Xu stratification; robust GKR); Sawin–Forey–Fresán–Kowalski, Quantitative sheaf theory (JAMS 36, 2023); Kedlaya, Étale and crystalline companions II; Adolphson–Sperber (DOI 10.1007/s12188-021-00243-1); Beukers–Vlasenko Dwork Crystals I–III (DOI 10.1093/imrn/rnad101); Caruso–Fürnsinn–Vargas-Montoya–Zudilin, Galois groups of Apéry-like series modulo primes (DOI 10.1017/S0004972725100932, publ. 2026-02-13; Franel equation = their (2.1)); Gabber–Loeser, Duke 83 (1996).

## Final calibration

The two-sheaf inverse theorem is a legitimate exact-char-0 research program; strongest exploitable feature = Galois-orbit amplification of high-order exact zeros. Not presently a theorem; FFK does not supply its base case; and even if proved it would not count the actual Apéry zeros without the independent (DRS)/(RLL) residue-separation lemma. The crystalline companion supplies the right operator on each tame block; the missing theorem is a uniform local limit across all p-1 blocks.
