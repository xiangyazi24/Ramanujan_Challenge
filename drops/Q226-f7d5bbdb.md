ANSWER Q226 f7d5bbdb

# Verdict

**No presently available theorem on Picard-rank jumps, ordinary/supersingular reduction, Artin--Tate, or K3 moduli proves**

\[
\sum_{p\le X} Z(p)\ll \pi(X),
\qquad
Z(p)=\#\{0\le j<p:b_j\equiv0\pmod p\},
\]

or even the weaker target

\[
\sum_{p\le X}Z(p)=o(X).
\]

The obstruction is more fundamental than a missing quantitative estimate. The proposed argument identifies two different zero loci:

- Picard rank, height, ordinarity, and supersingularity concern **values of the K3 family at base points** \(a\in\mathbb F_p\);
- \(Z(p)\) concerns **zero coordinates / multiplicative Mellin coefficients** indexed by \(j\in\mathbb Z/(p-1)\mathbb Z\).

These are Fourier-dual pieces of information. Even complete knowledge of the ordinary/supersingular locus, or even of the Picard rank of every fiber \(X_a\), does not control how many Mellin coefficients vanish.

There is one important positive statement. Because the Apéry K3 family has a generic rank-\(19\) algebraic lattice, its smooth reductions have a particularly simple height dichotomy: a finite-height fiber must be ordinary, while a nonordinary smooth fiber is supersingular. Thus the original family geometry gives an unusually clean description of the **evaluation Hasse locus**. But the repository's exact Mellin identities show that this still does not bound \(Z(p)\).

The correct geometric route would have to construct and analyze a **horizontal, defining-characteristic Mellin transform** of the primitive K3 sheaf. That would be a genuinely new theorem, not an application of current Picard-rank distribution results.

This audit is pinned to repository `main` commit

```text
f914983116049a10ced0844172d7e7f99a0eb238
```

and in particular to:

- `problems/3.2/oracleB_result.tex`;
- `problems/3.2/toric_mellin_square.tex`;
- `problems/3.2/toric_fiber_k3.tex`;
- `problems/3.2/toric_fiber_sym2.tex`;
- `problems/3.2/hasse_franel_descent.tex`;
- `problems/3.2/campaign3_questions/q81_avgzero.txt`;
- `problems/3.2/ERRATA.md`.

# 1. Three distinct arithmetic problems are being conflated

It is essential to distinguish three axes.

| Object | Variables | What existing geometry controls |
|---|---|---|
| A fixed K3 surface \(X_{a_0}/K\) | vary the prime \(p\) | ordinary/supersingular primes, Frobenius, Picard-rank jumps, Chebotarev |
| The K3 family in characteristic \(p\) | vary \(a\in\mathbb F_p\) | Hasse divisor, height strata, Picard ranks of the evaluated fibers \(X_a\) |
| The Apéry zero set | vary the character/index \(j\in\mathbb Z/(p-1)\mathbb Z\) | zeros of a defining-characteristic Mellin transform |

The fixed-surface theorems of Bogomolov--Zarhin, Charles, Hui, and Shankar--Shankar--Tang--Tayou concern the first row. The height-stratification theorems of van der Geer--Katsura concern the second row. The desired average \(Z\)-bound concerns the third row.

There is no formal implication from either of the first two rows to the third.

# 2. The repository already contains the exact separation theorem

Put

\[
A_p(T)=\sum_{j=0}^{p-1}b_jT^j\in\mathbb F_p[T].
\]

Then

\[
\mathcal Z_p=\{j:[T^j]A_p(T)=0\}.
\]

For \(1\le j\le p-2\), multiplicative orthogonality gives the exact identity recorded in `oracleB_result.tex`:

\[
\boxed{
 b_j=-\sum_{t\in\mathbb F_p^\times}A_p(t)t^{-j}.
}
\]

Thus \(b_j\) is a finite Mellin coefficient of the **evaluation function** \(t\mapsto A_p(t)\). Consequently,

\[
\boxed{
\mathcal Z_p
\text{ is a coefficient/Mellin zero set, not the root locus of }A_p.
}
\]

The toric K3 realization makes the same distinction geometrically. If \(\mu_p(a)\) is the toric fiber-count function, then the repository proves

\[
\boxed{
 b_j=-\sum_{a\in\mathbb F_p^\times}\mu_p(a)a^j
}
\]

for the nontrivial character range. After removing the explicit algebraic baseline, \(b_j\) is the defining-characteristic Mellin coefficient of the primitive K3 trace function.

Every \(b_j\) therefore depends on **all fibers simultaneously**, ordinary and supersingular alike. It is not attached to one fiber whose Picard rank could be inspected.

## A concrete repository counterexample

At \(p=7\), `oracleB_result.tex` computes

\[
A_7(T)
 =1+5T+3T^2+3T^3+3T^4+5T^5+T^6
 =(T-1)^2(T^2+1)^2\pmod7.
\]

All seven coefficients are nonzero, so

\[
Z(7)=0,
\]

while \(A_7\) has a nonempty evaluation root locus. Thus evaluation degeneration and coefficient vanishing already disagree at the smallest relevant scale.

# 3. The Franel/Hasse description does not repair the mismatch

The repository proves the exact cover-adapted identity

\[
K_p(x)^2=(1+x)^{p-1}A_p\!\left(\frac{x(1-8x)}{1+x}\right)
       =\sum_{m=0}^{p-1}b_m\Psi_{p,m}(x),
\]

where the \(\Psi_{p,m}\) form a basis of the relevant polynomial space. Hence

\[
Z(p)=\#\{m:[K_p^2]_{\Psi_{p,m}}=0\}.
\]

Even if \(K_p\), or equivalently \(A_p\) after the displayed pullback and gauge, is identified with a Hasse section of the K3/elliptic square-root family, \(Z(p)\) is the number of vanishing **coordinates of that section in a special basis**. It is not the number of zeros of the section.

The CFVMZ factorization in the repository is

\[
A_p(a)=\Delta(a)^{\varepsilon_p}B_p(a)^2,
\qquad
\Delta(a)=a^2-34a+1,
\]

with

\[
\deg B_p=\frac{p-1-2\varepsilon_p}{2}=\Theta(p).
\]

This is itself a warning against the proposed bounded-complexity argument. A fixed geometric family can produce a new Hasse section in every characteristic whose degree grows linearly with \(p\). Geometry of the original family does not turn that moving section into a fixed bounded-degree divisor.

The height-stratum formulas of van der Geer--Katsura are compatible with this behavior: the Hasse/height locus is a geometric divisor whose class depends on \(p\). Pulling it back to a one-dimensional family controls an **evaluation divisor of degree growing with \(p\)**. It does not control which coordinates of the defining equation vanish.

# 4. What the rank-\(19\) geometry really gives

The toric calculation in `toric_fiber_k3.tex` exhibits a generic Frobenius-stable algebraic subspace of dimension \(19\). The primitive/transcendental part is generically rank \(3\), and `toric_fiber_sym2.tex` identifies it with a symmetric-square/Asai construction from elliptic data.

Let \(Y/\overline{\mathbb F}_p\) be a smooth fiber with formal Brauer height \(h\). If \(h<\infty\), the Newton slopes on \(H^2_{\mathrm{cris}}\) are

\[
1-\frac1h\quad(h\text{ times}),
\qquad
1\quad(22-2h\text{ times}),
\qquad
1+\frac1h\quad(h\text{ times}).
\]

Algebraic divisor classes have slope \(1\), so, using the Tate conjecture in odd characteristic,

\[
\rho(Y)\le 22-2h.
\]

For an Apéry fiber carrying the rank-\(19\) lattice,

\[
19\le \rho(Y)\le22-2h.
\]

If \(h<\infty\), this forces \(h=1\). Therefore, away from the singular parameters and small characteristics,

\[
\boxed{
\text{a smooth Apéry K3 fiber is either ordinary }(h=1)
\text{ or supersingular }(h=\infty).
}
\]

This is a strong and useful geometric simplification. It means that, after the Hasse-section identification is made carefully, the evaluation equation “Hasse \(=0\)” detects supersingular fibers rather than intermediate finite-height fibers.

It still does not imply anything like

\[
Z(p)\le 1\quad\text{at ordinary primes}.
\]

There is no single “ordinary prime” for the family: at a fixed prime \(p\), different parameters \(a\in\mathbb F_p\) can produce different fibers. More importantly, \(Z(p)\) is a statistic of the Mellin transform of the whole fiber function.

# 5. Ordinary reduction does not mean Picard rank is preserved

One premise in the question should be corrected.

For a fixed K3 surface \(X/K\), specialization injects the characteristic-zero Néron--Severi group into that of a good reduction, so

\[
\rho(X_{\overline{\mathbb F}_p})\ge \rho(X_{\overline K}).
\]

Ordinarity does **not** force equality. For an ordinary K3, the slope-\(1\) part has dimension \(20\), so it only gives

\[
\rho(X_{\overline{\mathbb F}_p})\le20.
\]

Thus a rank-\(19\) K3 can have ordinary reductions of Picard rank \(19\) or \(20\), and supersingular reductions of Picard rank \(22\). Picard-rank preservation is governed by the transcendental endomorphism field and Frobenius, not by ordinarity alone.

Charles's density-one theorem describes the generic specialization rank in terms of the endomorphism field of the transcendental Hodge structure. Shankar--Shankar--Tang--Tayou prove that exceptional Picard-rank jumps nevertheless occur at infinitely many primes under potential good reduction. Neither result supplies a bound for Mellin coefficient zeros.

# 6. Picard-rank data is information-theoretically insufficient

There is a simple finite Fourier obstruction which does not depend on any unproved geometry.

Let \(G=\mathbb F_p^\times\), and let

\[
\widehat f(\chi)=\sum_{a\in G}f(a)\chi(a).
\]

Consider, for \(p>3\),

\[
f_1(a)=1,
\qquad
f_2(a)=1+2\mathbf1_{a=1}.
\]

Both functions are nonzero at every point of \(G\), so their evaluation zero loci are identical: both are empty. But

- \(\widehat f_1(\chi)=0\) for every nontrivial character;
- every Mellin coefficient of \(f_2\) is nonzero.

Indeed, for nontrivial \(\chi\), \(\widehat f_2(\chi)=2\), while for the trivial character it equals \((p-1)+2=1\pmod p\).

Therefore even the **complete evaluation zero pattern** does not bound the number of zero Mellin coefficients. The multiset of Picard ranks contains still less information than the complete evaluation function.

This gives a rigorous no-go statement:

> Any argument using only the distribution of ordinary, supersingular, or Picard-rank strata of the original K3 fibers cannot prove an upper bound for \(Z(p)\). One needs arithmetic information about the actual Frobenius trace values and their character-weighted cancellations.

A marked coordinate and the full trace function can of course contain more information than the Picard ranks. The point is not that geometry is irrelevant, but that **Picard-rank distribution by itself is the wrong invariant**.

# 7. Audit of the relevant state of the art

## 7.1 Bogomolov--Zarhin: density-one ordinary reduction

Bogomolov and Zarhin prove that, after a finite extension of the ground number field, a fixed K3 surface has ordinary reduction at a density-one set of places.

This is qualitative, fixed-surface, prime-aspect information. It neither counts fibers \(a\in\mathbb F_p\) of a family nor controls the character-indexed Mellin coefficients \(b_j\).

It also gives no quantitative exceptional-set estimate strong enough to absorb the current \(p^{2/3}\) bound, even if one postulated a connection between supersingularity and large \(Z(p)\).

## 7.2 Hui 2025: the strongest directly relevant fixed-surface bound

For a **fixed non-CM** abelian variety or K3 surface over a number field, Chun-Yin Hui proves that the supersingular primes have density zero. Writing \(G_\ell\) for the algebraic monodromy group, the paper gives, for every \(\epsilon>0\), bounds of the form

\[
\pi_{SS_X}(x)
 \ll_\epsilon
 \frac{x}{(\log x)^{1+1/\dim G_\ell-\epsilon}},
\]

and, under GRH,

\[
\pi_{SS_X}(x)
 \ll_\epsilon
 x^{1-1/(2\dim G_\ell)+\epsilon}.
\]

When \(G_\ell\) is connected, the codimension can be sharpened using
\((\operatorname{rk}G_\ell-1)/\dim G_\ell\). Hui also proves convergence of

\[
\sum_{v\in SS_X}\frac1{\#\mathbb F_v}.
\]

This is a genuine advance in the state of the art, but it still concerns one fixed \(X\) and one Frobenius conjugacy class per prime. It does not concern the \(p-1\) Kummer characters that occur in \(Z(p)\).

Moreover, the non-CM assumption matters. For CM/singular K3 surfaces, supersingular primes can have positive density, analogous to inert primes for a CM elliptic curve. Thus no uniform “supersingular primes are very rare for every K3” principle is available.

## 7.3 Maulik--Poonen: specialization loci, not prime counts

Maulik--Poonen prove two kinds of statements:

1. in a characteristic-zero smooth proper family there exists a closed fiber with the same Picard rank as the geometric generic fiber;
2. for a family with good reduction at a fixed \(p\), the Picard-rank jumping locus is nowhere dense in the relevant \(p\)-adic analytic topology.

“Nowhere \(p\)-adically dense” is not a quantitative estimate for

\[
\#\{a\in\mathbb F_p:\rho(X_a)>\rho_{\rm gen}\},
\]

and certainly is not a count over varying rational primes. It allows no estimate for the coordinate zeros of a moving Hasse polynomial. There is no uniformity in the moving characteristic and no Mellin transform in the theorem.

## 7.4 Shankar--Shankar--Tang--Tayou: infinitely many jumps

The March 22, 2026 version of their theorem proves, under potential good reduction, that a fixed K3 has infinitely many primes at which the geometric Picard rank jumps.

This is a lower-bound/existence theorem for exceptional primes, not an upper bound. The arithmetic-special-divisor technology in that work is extremely powerful, but its counted intersections encode extra algebraic classes. They do not encode vanishing of the defining-characteristic Mellin coefficient \(b_j\).

## 7.5 van der Geer--Katsura and Ekedahl--van der Geer

Their height and Artin-invariant stratifications compute cycle classes of loci in the moduli of polarized K3 surfaces. These results are the right tools for questions such as:

- what is the class of the nonordinary locus?
- what is the expected intersection number of a family curve with a height stratum?
- how large is the supersingular locus in K3 moduli?

After pulling back along the Apéry family, they can constrain the evaluation Hasse divisor. They do not constrain the zero coordinates of the pulled-back section. In fact, the linear-in-\(p\) degree of \(B_p\) is exactly the sort of moving degree that these moduli formulas allow.

## 7.6 Tate and Artin--Tate

For a K3 surface \(Y/\mathbb F_q\), the Tate conjecture identifies \(ho(Y)\) with the multiplicity of Frobenius eigenvalues of the form \(q\zeta\), with \(\zeta\) a root of unity. The Artin--Tate formula relates the leading term of the zeta function at the Tate pole to

- the discriminant of \(\operatorname{NS}(Y)\);
- the order of the Brauer group;
- elementary powers of \(q\).

These are pointwise identities for one finite-field surface. They do not give an inequality coupling different parameters \(a\), different characters \(j\), or different residue characteristics \(p\). In particular:

- they do not control the phases of ordinary Frobenius traces;
- they do not control cancellation in
  
  \[
  \sum_a\theta_p(a)a^j;
  \]
- they do not turn a Picard-rank jump count into a Mellin anti-concentration estimate.

Artin--Tate is therefore not the missing averaging theorem.

# 8. The elliptic-curve analogy is not the right analogue

For a fixed non-CM elliptic curve \(E/\mathbb Q\), there is one Frobenius trace \(a_p(E)\) per prime. For \(p>3\),

\[
E\bmod p\text{ is supersingular}
\iff a_p(E)=0.
\]

That is one binary event per prime. By contrast, \(Z(p)\) counts up to \(p\) character-indexed events at the same prime.

Two corrections are useful:

1. Elkies's celebrated theorem proves **infinitude** of supersingular primes for elliptic curves over \(\mathbb Q\); the classical \(x^{3/4}\)-type upper bound is a GRH/Chebotarev bound, not Elkies's infinitude theorem.
2. If one varies the elliptic curve inside a universal family at a fixed \(p\), there are on the order of \(p/12\) supersingular \(j\)-invariants over \(\overline{\mathbb F}_p\). Thus even in dimension one, “fixed curve across primes” and “family fibers at one prime” are entirely different counting problems.

The closer elliptic analogue of the Apéry question would be:

> Take the Deuring/Hasse polynomial of an elliptic family, multiplicatively Fourier transform its value function, and ask how many Fourier coefficients vanish modulo the same characteristic.

Standard theorems on supersingular primes of a fixed elliptic curve do not answer that question either.

# 9. Even a hypothetical ordinary/supersingular decomposition is quantitatively too weak

Suppose, much more strongly than is known, that one could find an exceptional set \(E\) of primes such that

\[
Z(p)\le C\quad(p\notin E),
\qquad
Z(p)\ll p^{2/3}\quad(p\in E).
\]

Then

\[
\sum_{p\le X}Z(p)
 \ll \pi(X)+X^{2/3}E(X).
\]

To obtain the expected-strength estimate \(O(\pi(X))\), one would need

\[
E(X)\ll \frac{X^{1/3}}{\log X}.
\]

To obtain merely \(o(X)\), one still needs

\[
E(X)=o(X^{1/3}).
\]

Density zero,

\[
E(X)=o(\pi(X)),
\]

is vastly too weak. Hui's unconditional logarithmic saving, inserted into this fantasy decomposition, still leaves a contribution near \(X^{5/3}\) up to logarithms. Even the GRH power-saving exponent furnished by the dimension of a fixed monodromy group is nowhere close to the \(>2/3\) saving needed after multiplying by \(X^{2/3}\).

This calculation is secondary, because the decomposition itself is unjustified: ordinary primes have not been shown to contribute \(O(1)\) to \(Z(p)\). But it demonstrates that qualitative Picard-rank rarity would not close the average even if the basic geometric bridge existed.

# 10. Why global \(H^2\) data does not automatically control the mod-\(p\) zeros

The global compatible system on \(H^2\) does impose strong restrictions:

- purity and determinant;
- local monodromy;
- Hodge and Newton polygons;
- Sato--Tate/Chebotarev distribution of fixed Frobenius conjugacy classes;
- Tate classes and Picard-rank jumps.

The event \(b_j\equiv0\pmod p\) has two extra moving features:

1. the Kummer character \(\chi_j\) varies through the full character group of \(\mathbb F_p^\times\);
2. the prime used to reduce the Mellin value is the same prime \(p\) that defines the finite field and the character group.

This is a **same-characteristic, moving-character divisibility event**. Ordinary \(\ell\)-adic equidistribution fixes \(\ell\) and controls complex or \(\ell\)-adic Frobenius classes. It does not control divisibility by the varying residual characteristic of the resulting algebraic Mellin sums.

The central index already displays the difficulty. The repository records

\[
b_{(p-1)/2}\equiv a_p(f)\pmod p,
\qquad
f=\eta(2\tau)^4\eta(4\tau)^4.
\]

Thus even one moving index contains the nonordinary-prime problem for a non-CM weight-\(4\) modular form. Ordinary Sato--Tate controls fixed real intervals for normalized \(a_p(f)\); it does not settle the shrinking divisibility condition \(p\mid a_p(f)\).

# 11. The exact new theorem that would make geometry work

Let \(\mathcal T\) denote the primitive rank-\(3\) K3 trace object on the parameter torus, after removing the explicit algebraic classes. For a multiplicative character \(\chi\), the natural cohomological Mellin fiber is schematically

\[
\mathcal M_{p,\chi}
 =R\Gamma_c\!\left(\mathbb G_{m/\mathbb F_p},
   \mathcal T_p\otimes\mathcal L_\chi\right).
\]

What is needed is not another theorem about \(ho(X_a)\). It is an integral/crystalline theorem identifying a marked scalar in \(\mathcal M_{p,\chi_j}\) whose reduction is \(b_j\), together with prime-aspect anti-concentration.

A clean sufficient statement would be:

```text
HORIZONTAL CRYSTALLINE MELLIN ANTI-CONCENTRATION.
There is an integral compatible Mellin object M for the primitive Apéry
K3 sheaf such that, for every good prime p and every Teichmüller character
χ_j of F_p^×, a marked Frobenius/unit-root coordinate u_{p,j} satisfies

    u_{p,j} mod p = 0  <=>  b_j mod p = 0,

and uniformly in X,

    Σ_{p≤X} #{j mod p-1 : u_{p,j} mod p = 0}  <<  π(X).
```

A proof would require at least the following inputs.

1. **A genuine horizontal parameter space for characters.** The finite sets
   \(\widehat{\mathbb F_p^\times}\) must be packaged compatibly as \(p\) varies. An ordinary algebraic \(j\)-line does not do this.
2. **Integral/crystalline control.** A bounded-rank \(\ell\)-adic Mellin fiber is not enough; one needs a lattice and a marked coordinate whose reduction modulo the varying prime is the Apéry coefficient.
3. **Uniform conductor and monodromy.** One must exclude Kummer self-twists and prove that the Mellin family has no exceptional components causing systematic zeros.
4. **A prime-character large sieve or effective Chebotarev theorem.** It must be uniform simultaneously in \(p\) and in characters of order up to \(p-1\).
5. **Defining-characteristic anti-concentration.** Complex Sato--Tate and Deligne bounds control sizes of trace sums, not their divisibility by \(p\).

No theorem currently combines these five properties for the Apéry K3 family.

A weaker but still sufficient target is the first-moment statement itself,

\[
\sum_{p\le X}\sum_{j=0}^{p-1}
 \mathbf1_{b_j\equiv0\ (p)}\ll\pi(X),
\]

proved by a mixed prime-character large sieve. Another plausible route is a uniform second-factorial-moment estimate for the Mellin zero sets followed by a separate first-moment calculation. Both are new horizontal inputs, not consequences of Picard-rank distribution.

# 12. What the existing K3 geometry can still contribute

The geometry is not useless. It has already supplied several indispensable pieces:

- a fixed toric K3 family realizing the fiber-count function;
- an explicit rank-\(19\) algebraic lattice and rank-\(3\) primitive part;
- an exact symmetric-square/Asai description;
- explicit discriminant and exceptional parameters;
- a defining-characteristic Mellin identity for \(b_j\);
- strong local monodromy and candidate large geometric monodromy;
- a clean ordinary/supersingular dichotomy for smooth rank-\(19\) fibers.

These are precisely the ingredients from which a horizontal Mellin theorem might eventually be built. But they stop one logical step before the target. The missing step is not “count supersingular primes more accurately”; it is “control reduction modulo \(p\) of all character Mellin fibers uniformly in the moving characteristic.”

# 13. Final answer to the proposed routes

## Can Picard-rank distribution prove average \(Z\)?

**Not by itself.** Picard ranks concern evaluated fibers; \(Z(p)\) concerns Mellin coefficients. The finite Fourier counterexample above shows that no argument based only on rank/ordinary strata can control the number of zero coefficients.

## Does Artin--Tate constrain \(\sum_p Z(p)\)?

**No known implication.** Artin--Tate is pointwise in a finite-field K3 and relates a zeta-function leading term to Néron--Severi and Brauer invariants. It contains no cross-character or cross-prime anti-concentration statement.

## Do Maulik--Poonen-type theorems give the required average?

**No.** Their smallness notion is geometric or \(p\)-adic nowhere density for Picard-jump loci. It is not a uniform finite-field point count, and it does not see Mellin coordinates.

## Do recent supersingular-prime bounds help?

**They are the wrong quantifier and, even under a hypothetical bridge, quantitatively insufficient.** Hui's 2025 theorem is for one fixed non-CM surface; \(Z(p)\) uses the whole family and all \(p-1\) moving characters.

## Is there a viable non-circular geometric program?

**Yes, but it must change objects.** The promising object is the integral/crystalline multiplicative Mellin transform of the primitive Apéry K3 sheaf, together with a prime-character large sieve. The Picard lattice and K3 trace formulas would be inputs to that program, not the final averaging mechanism.

# Bottom line

The proposed “ordinary primes contribute \(O(1)\), supersingular primes are rare” strategy does not work. Its first implication is false, its exceptional-prime estimate is too weak, and it targets evaluation geometry rather than coefficient/Mellin zeros.

The sharp current conclusion remains:

\[
\boxed{
\text{No unconditional improvement of }
\sum_{p\le X}Z(p)
\ll \frac{X^{5/3}}{\log X}
\text{ follows from known Picard-rank or supersingularity results.}
}
\]

The genuinely new theorem to pursue is a **horizontal defining-characteristic Mellin anti-concentration theorem**. That is the exact point where the global K3 geometry must be converted into an average-\(Z\) bound.

# References

## Repository sources

- `problems/3.2/oracleB_result.tex`, especially “Coefficient zeros are not Hasse-polynomial roots” and “Fixed family versus moving Hasse section.”
- `problems/3.2/toric_mellin_square.tex`, exact toric Mellin identities and the warning that the collapse is not an anti-concentration estimate.
- `problems/3.2/toric_fiber_k3.tex`, explicit rank-\(19\) lattice and rank-\(3\) primitive trace.
- `problems/3.2/toric_fiber_sym2.tex`, symmetric-square/Asai realization.
- `problems/3.2/hasse_franel_descent.tex`, cover-adapted coordinate identity and CFVMZ factorization.
- `problems/3.2/campaign3_questions/q81_avgzero.txt`, exact `[AVG-ZERO]` target and current unconditional \(p^{2/3}\) baseline.
- `problems/3.2/ERRATA.md`, defining-characteristic and rank-distinction cautions.

## External sources

1. F. A. Bogomolov and Yu. G. Zarhin, [*Ordinary reduction of K3 surfaces*](https://arxiv.org/abs/0902.1548).
2. Chun-Yin Hui, [*On distribution of supersingular primes of abelian varieties and K3 surfaces*](https://arxiv.org/abs/2504.08088), 2025.
3. D. Maulik and B. Poonen, [*Néron--Severi groups under specialization*](https://arxiv.org/abs/0907.4781), Duke Math. J. 161 (2012).
4. F. Charles, [*On the Picard number of K3 surfaces over number fields*](https://arxiv.org/abs/1111.4117), Algebra Number Theory 8 (2014).
5. A. N. Shankar, A. Shankar, Y. Tang, and S. Tayou, [*Exceptional jumps of Picard ranks of reductions of K3 surfaces over number fields*](https://arxiv.org/abs/1909.07473), version dated March 22, 2026.
6. G. van der Geer and T. Katsura, [*On a stratification of the moduli of K3 surfaces*](https://arxiv.org/abs/math/9910061).
7. G. van der Geer, [*A stratification on the moduli of K3 surfaces in positive characteristic*](https://arxiv.org/abs/1511.00345).
8. E. Costa and Yu. Tschinkel, [*Variation of Néron--Severi ranks of reductions of K3 surfaces*](https://arxiv.org/abs/1405.2265).
9. F. Charles, [*The Tate conjecture for K3 surfaces over finite fields*](https://arxiv.org/abs/1206.4002).
10. K. Madapusi Pera, [*The Tate conjecture for K3 surfaces in odd characteristic*](https://arxiv.org/abs/1301.6326).
