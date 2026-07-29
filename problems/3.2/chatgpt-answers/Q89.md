ANSWER Q89 a7fe1cf7

# Literature check: Apéry truncations, Hasse invariants, and the mod-p² complete sum

Search completed 2026-07-19. I use “not found” literally: I did not locate a source stating the claimed result after checking the named papers/authors and targeted searches. I distinguish statements actually written in a source from deductions that appear mathematically natural but still need a proof or a citable reference.

## Executive findings

1. Roots of B_p = parameters above supersingular j: not found as an explicit literature statement. Stienstra–Beukers do give the formal-Brauer/Hasse-Witt interpretation of the truncated Apéry period: on the smooth Apéry K3 pencil, nonvanishing of the first truncation is the height-one/ordinary criterion. But vanishing says nonordinary K3 (formal Brauer height at least 2), not by itself “Artin-supersingular K3.” The stronger identification with the pullback of the elliptic supersingular divisor should follow from a Frobenius-compatible symmetric-square/Shioda–Inose identification, with care at branch and singular points, but I did not find that statement written for the CFVZ polynomial B_p.

1. Truncated period = Hasse-Witt invariant is definitely in the literature. A very clean general reference is Huang–Lian–Yau–Yu, Theorem 1.2, equation (1.11); their equation (1.6) gives the Dwork-family prototype. However, the precise level-6 modular-form formula

1. Sun’s mod-p² conjecture was proved. Chen Wang and Zhi-Wei Sun, Theorem 1.1, prove exactly the conjectural congruence for Σ_{k=0}^{p−1} A_k. The preprint is arXiv:1910.06856; the published version is Nanjing University Journal of Mathematical Biquarterly 41 (2024), no. 1, 34–56.

---

## 1. Does the literature identify the roots of B_p with supersingular parameters?

### 1.1 What CFVZ actually prove

Xavier Caruso, Florian Fürnsinn, Daniel Vargas-Montoya, Wadim Zudilin, “Galois Groups of Apéry-like Series Modulo Primes,” Bulletin of the Australian Mathematical Society 114 (2026), no. 1, 65–78; published online 13 February 2026. DOI: 10.1017/S0004972725100932. Preprint: arXiv:2510.23298.

In their notation

A_p(t) = Σ_{n=0}^{p−1} b_n t^n.

Theorem 1.2 states that there is B_p ∈ 𝔽_p[t] such that

- A_p = B_p² for p ≡ 1,5,7,11 (mod 24);

- A_p = (t²−34t+1)B_p² for p ≡ 13,17,19,23 (mod 24).

The paper derives this from the rational substitution relating the Apéry and Franel generating functions and from the quadratic extension with discriminant t²−34t+1. The published full text contains no occurrence of “Hasse” or “supersingular.” Therefore:

CFVZ themselves do not state that the roots of B_p are supersingular parameters.

The factor t²−34t+1 is both the discriminant of their quadratic substitution and the finite singular/discriminant locus of the Apéry Picard–Fuchs equation. It must be separated from any assertion about the smooth supersingular locus.

### 1.2 What Stienstra–Beukers safely gives

Jan Stienstra and Frits Beukers, “On the Picard-Fuchs Equation and the Formal Brauer Group of Certain Elliptic K3-Surfaces,” Mathematische Annalen 271 (1985), no. 2, 269–304. DOI: 10.1007/BF01455990.

A companion exposition is:

Jan Stienstra, “Les groupes formels d’Artin–Mazur et les congruences d’Atkin–Swinnerton-Dyer,” Groupe d’étude d’analyse ultramétrique 12 (1984–1985), no. 2, exposé 18, 1–13; especially §§5–6. Numdam record.

The safe formal-group content is as follows. For a one-dimensional formal Brauer group whose logarithm is obtained from the normalized holomorphic period, the coefficient controlling the first Frobenius/Hasse-Witt map is the period truncated through degree p−1. For the Apéry K3 pencil, after the usual normalization at the maximally unipotent point, this scalar is

F_p(t) = Σ_{n=0}^{p−1} b_n t^n  (mod p).

Thus, at a smooth fiber in characteristic p, one has the first-height criterion

- F_p(t₀) ≠ 0  ⇔ the formal Brauer group has height 1 ⇔ the K3 fiber is ordinary;

- F_p(t₀) = 0  ⇒ the formal Brauer group has height at least 2, possibly infinite ⇔ the K3 fiber is nonordinary.

This is the precise point at which terminology matters:

F_p(t₀)=0 is a nonordinary-locus criterion. It does not, from the first Hasse invariant alone, imply that the K3 surface is supersingular in Artin’s sense (formal Brauer height ∞, equivalently Picard rank 22 under the standard conjectural/theorem context).

One needs higher Frobenius/Hasse-Witt data to distinguish finite heights 2,…,10 from height ∞.

I could not reliably recover a theorem number in the 1985 Math. Ann. scan for the exact sentence “ordinary iff this truncated Apéry series is nonzero.” I therefore would cite the paper and its formal-Brauer calculation, but not invent a theorem number. The modern theorem below in §2 supplies an unambiguous theorem-numbered truncation/Hasse-Witt statement.

### 1.3 The desired supersingular-j conclusion is plausible, but it is an additional theorem

The expected argument is:

1. Beukers’ modular parametrization realizes the third-order Apéry Picard–Fuchs equation as a symmetric square of a second-order elliptic Picard–Fuchs equation.

1. If this lifts integrally mod p to a Frobenius-compatible identification of Hodge lines/crystals,

1. Igusa’s theorem says Ha_ell has simple zeros exactly at supersingular elliptic points. Therefore Ha_K3 has double zeros over that divisor.

1. In a chosen t-coordinate and period trivialization, Ha_K3 is represented by F_p(t). Hence, on the smooth unramified locus, the square root B_p should cut out the reduced pullback of the elliptic supersingular divisor.

This explains exactly why the CFVZ double roots invite the supersingular interpretation. But several nontrivial details must be checked before citing it as a theorem:

- the symmetric-square relation must be integral and Frobenius-compatible, not merely an equality of complex differential equations;

- one must specify the modular curve before or after the Fricke quotient;

- ramification can change multiplicities in the Hauptmodul t;

- cusps, elliptic fixed points, and the singular fibers t²−34t+1=0 must be excluded or treated separately;

- p=2,3 and primes of bad reduction must be excluded.

Therefore my answer to the literal question is:

Not found: I did not find a paper stating that the roots of the CFVZ polynomial B_p are exactly the t-parameters lying over supersingular j-invariants for the Apéry K3 family.

A defensible formulation for a new proposition would be set-theoretic and restricted to the good smooth locus:

V(B_p) = image/pullback in the t-line of the supersingular locus of X₀(6)

for p>3, with multiplicities handled separately and the two discriminant points removed.

### 1.4 Candidate literature checked

- Frits Beukers and Masha Vlasenko, “Dwork Crystals I,” International Mathematics Research Notices 2021, no. 12, 8807–8844, DOI 10.1093/imrn/rnaa119; “Dwork Crystals II,” IMRN 2021, no. 6, 4427–4444, DOI 10.1093/imrn/rnaa120; “Dwork Crystals III: From Excellent Frobenius Lifts Towards Supercongruences,” IMRN 2023, no. 23, 20433–20483, DOI 10.1093/imrn/rnad101. These provide general Cartier/Dwork-crystal and truncated-period machinery. The exact B_p/supersingular-j statement was not found.

- Jeng-Daw Yu, “Local Structure of the Moduli Space of K3 Surfaces in Positive Characteristic,” International Mathematics Research Notices 2009, no. 23, 4480–4495, DOI 10.1093/imrn/rnp096, and Jeng-Daw Yu and Noriko Yui, “K3 Surfaces of Finite Height over Finite Fields,” Journal of Mathematics of Kyoto University 48 (2008), no. 3, 499–519, DOI 10.1215/kjm/1250271381. These concern K3 height/ordinary strata and examples, but I found no theorem identifying the Apéry B_p locus with supersingular elliptic j-values.

- Marie-José Bertin and Odile Lecacheux, “Apéry–Fermi Pencil of K3-Surfaces and 2-Isogenies,” Journal of the Mathematical Society of Japan 72 (2020), no. 2, 599–637, DOI 10.2969/jmsj/80638063. This is useful for the geometry/isogenies of the pencil, but the desired mod-p B_p statement was not found.

- Targeted searches of work by Matthias Schütt and Duco van Straten likewise produced no exact statement. Not found.

---

## 2. Truncated period times a period form as a Hasse invariant

### 2.1 A theorem-numbered general reference

An Huang, Bong Lian, Shing-Tung Yau, Chenglong Yu, “Hasse–Witt Matrices, Unit Roots and Period Integrals,” Mathematische Annalen 387 (2023), 145–173. DOI: 10.1007/s00208-022-02464-y. Preprint: arXiv:1801.01189.

Two exact places answer the general question:

- Equation (1.6): for the Dwork family, the Hasse-Witt polynomial is the degree-p−1 truncation of the normalized hypergeometric period.

- Theorem 1.2, equation (1.11): for Calabi–Yau hypersurfaces in the stated toric/flag setting, if the period near the distinguished degeneration is

So the answer to “is this written for any family?” is an unequivocal yes.

For the Legendre family, a convenient modern reference is:

Alan Adolphson and Steven Sperber, “Hasse Invariants and Mod p Solutions of A-Hypergeometric Systems,” Journal of Number Theory 142 (2014), 183–210. DOI: 10.1016/j.jnt.2014.02.010. Preprint: arXiv:1209.2448.

They explicitly recall Igusa’s observation that the Legendre Hasse invariant is a mod-p hypergeometric solution. The traditional Legendre Hasse polynomial is equivalently the coefficient of x^{p−1} in [x(x−1)(x−λ)]^{(p−1)/2}.

### 2.2 The exact Γ₀(6)+6 formula

For the Apéry family, let E(τ) be the chosen holomorphic-period modular form and t=t(τ) the Hauptmodul. In a local trivialization of the K3 Hodge line by E, the general geometric transformation rule gives an expression of the form

Ha_K3(τ) = u · F_p(t(τ)) · E(τ)^(p−1)

where u ∈ 𝔽_p^× is fixed by a normalization, normally by the q-expansion at the cusp. This is the natural global-section version of “truncation equals Hasse-Witt scalar.”

However:

Not found: I did not locate this exact formula, with Γ₀(6)+6, Beukers’ eta-quotient t, and the Apéry truncation F_p, stated as a published theorem.

This should be regarded as a short proposition to prove from the integral modular parametrization plus the general Hasse-Witt/truncation theorem, rather than as something presently backed by a precise citation.

### 2.3 Weight correction

If E(τ) has modular weight 2, then

F_p(t(τ)) E(τ)^(p−1)

has weight 2(p−1). That is appropriate for the Hasse invariant of a K3 Hodge line whose automorphic realization has weight 2.

The classical elliptic Hasse invariant is a section of ω^(p−1) and has weight p−1. Thus:

- if A_{p−1} denotes the K3 Hasse invariant, the displayed weight is consistent;

- if A_{p−1} denotes the elliptic Hasse invariant, the expected symmetric-square identity is instead

A putative square-root expression involving B_p would have weight p−1, but its global formulation can acquire a character or a branch/discriminant factor. That is precisely where the residue classes mod 24 and the factor t²−34t+1 may enter.

### 2.4 Simple zeros exactly at supersingular points

The standard reference is:

Nicholas M. Katz and Barry Mazur, Arithmetic Moduli of Elliptic Curves, Annals of Mathematics Studies 108, Princeton University Press, 1985; see §12.4.3, the Igusa theorem on the Hasse invariant.

The result, on a representable prime-to-p level cover, is:

- the Hasse invariant vanishes exactly on the supersingular locus;

- every zero is simple, so its divisor is the reduced supersingular divisor.

On a stack this is the clean statement. When pushed to a coarse modular curve/Hauptmodul, elliptic automorphisms and ramification can modify the apparent multiplicity in the coarse coordinate. This caveat matters when comparing the simple roots of an elliptic Hasse polynomial with the double roots of the K3 polynomial F_p.

### 2.5 Supersingular points of X₀(N) are 𝔽_{p²}-rational

For p ∤ N, every supersingular point of X₀(N) is defined over 𝔽_{p²}. An explicit source that states this is:

János A. Csirik, Joseph L. Wetherell, Michael E. Zieve, “On the Genera of X₀(N),” preprint arXiv:math/0006096, introduction/motivation.

A direct proof is also standard. Every supersingular j lies in 𝔽_{p²}; one can choose an 𝔽_{p²}-model whose p²-Frobenius acts as a scalar ±p. For p∤N, that scalar action preserves every cyclic subgroup of order N, so the pair (E,C) defining an X₀(N) point descends to 𝔽_{p²}. Consequently, the image on the Fricke quotient X₀(N)^+ is also 𝔽_{p²}-rational.

For the Apéry level N=6, this applies for p>3.

---

## 3. Sun’s complete Apéry sum modulo p²

### 3.1 Original mod-p reference

Zhi-Wei Sun, “On Sums of Apéry Polynomials and Related Congruences,” Journal of Number Theory 132 (2012), no. 11, 2673–2699. DOI: 10.1016/j.jnt.2012.05.014. Preprint: arXiv:1101.1946.

As in the question, Corollary 1.2 gives the mod-p evaluation using (-2/p) and the representation p=x²+2y²; Remark 1.2, following equation (1.9), conjectures the strengthening modulo p².

### 3.2 The conjectured mod-p² formula was proved

Chen Wang and Zhi-Wei Sun, “p-adic Analogues of Hypergeometric Identities and Their Applications,” Nanjing University Journal of Mathematical Biquarterly 41 (2024), no. 1, 34–56. Preprint: arXiv:1910.06856.

Theorem 1.1 proves, for every odd prime p,

Σ_{k=0}^{p−1} A_k ≡ 4x²−2p (mod p²) if p=x²+2y²,

and

Σ_{k=0}^{p−1} A_k ≡ 0 (mod p²) if p≡5,7 (mod 8).

In Wang–Sun this target is displayed as equation (1.1). Their Remark 1.1 identifies Sun’s earlier result as the mod-p predecessor. This is the exact affirmative answer to the question.

### 3.3 Relation to Mortenson and whether there is a cleaner modular-form proof

The hypergeometric/newform input cited in this circle is:

Eric Mortenson, “Supercongruences for Truncated ${}_{n+1}F_n Hypergeometric Series with Applications to Certain Weight Three Newforms,” Proceedings of the American Mathematical Society 133 (2005), no. 2, 321–330. DOI: 10.1090/S0002-9939-04-07697-X.

Wang–Sun is the cleanest single later reference I found because it proves the full p² assertion, not merely the mod-p evaluation. Its method is framed through p-adic analogues of hypergeometric identities and uses the relevant truncated-hypergeometric supercongruences.

I also checked:

Robert Osburn, Brundaban Sahu, Armin Straub, “Supercongruences for Sporadic Sequences,” Proceedings of the Edinburgh Mathematical Society 59 (2016), no. 2, 503–518. DOI: 10.1017/S0013091515000255. Preprint: arXiv:1312.2195.

That paper proves two-term supercongruences for sporadic sequences and surveys related directions; I did not find the exact complete-sum congruence Σ_{k<p} A_k presented there as a simpler modular-form theorem.

Not found: I did not locate a later Kilbourn/Osburn–Sahu–Straub-style modular-form proof that supersedes Wang–Sun as a cleaner proof of this exact complete-sum evaluation. The modular-form interpretation of Mortenson’s ${}_4F_3 value explains the CM quantity, but the direct reference proving Sun’s full equation (1.9) modulo p² is Wang–Sun, Theorem 1.1.

---

## Recommended citation-level conclusions

For a paper or formalization note, I would state the conclusions this way:

1. Proved/citable: The Apéry period truncation is the first Hasse-Witt/formal-Brauer invariant; nonzero is equivalent to ordinary height 1. Cite Stienstra–Beukers (1985) together with Huang–Lian–Yau–Yu, Theorem 1.2, equation (1.11), for a theorem-numbered modern formulation.

1. Not yet located as a citation: V(B_p) is exactly the reduced supersingular divisor pulled back through Beukers’ modular parametrization. Present it as a proposition to prove, not as a known quoted theorem, unless a more specific source is found.

1. Citable standard facts: The elliptic Hasse invariant has reduced zero divisor equal to the supersingular locus—Katz–Mazur §12.4.3—and supersingular X₀(N) points are 𝔽_{p²}-rational for p∤N.

1. Solved conjecture: Sun’s Remark 1.2 mod-p² conjecture is Wang–Sun, Theorem 1.1, arXiv:1910.06856 / Nanjing Univ. J. Math. Biquarterly 41 (2024), 34–56.