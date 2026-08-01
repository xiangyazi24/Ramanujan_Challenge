ANSWER Q6393 eb41c725

# Executive verdict

The hypergeometric-motive reframing is mathematically correct but does not produce a new counting bound.

For every nontrivial character

the Greene value in the question is the trace at the balanced point \lambda=1 of the Tate-normalized Katz hypergeometric system

Generically in \lambda this is an irreducible rank-four, weight-three, symplectically self-dual hypergeometric local system of Calabi–Yau type, with Hodge numbers

However, \lambda=1 is not a lisse fiber. It is the balanced conifold point. The local monodromy is a rank-one unipotent pseudoreflection, and the middle-extension stalk has rank three, not four. In the normalization matching the Apéry formula, its Frobenius polynomial is

where

and, at the chosen Teichmüller prime \mathfrak p\mid p,

Consequently

This equivalence is exact. It is also a Hasse-invariant tautology: the Hasse scalar is precisely b_r modulo p. It does not replace b_r by a lower-complexity function.

There is a universal Hasse section over the finite character scheme

but \mathscr X_p is a zero-dimensional scheme of length p-1, not a fixed algebraic curve. Its coordinate ring contains every function on the character set. The Hasse section has a unique representative of degree at most p-2, and the exact Gross–Koblitz/Stickelberger calculation shows that its unit part contains 1+\min(r,p-1-r) summands—generically \asymp p—whose reduction is exactly the original Apéry binomial sum. No bounded-degree polynomial emerges.

Thus the proposed degree argument does not give O(1), O(\log p), or even p^{1-\delta} zeros. The current unconditional record \#\{r:b_r=0\}\ll p^{3/4} is not improved by this reframing.

---

# 1. From Greene’s 4F_3 to Katz’s balanced hypergeometric sheaf

Let N=p-1, let \psi be a character of order N, and put

The three displayed lower characters in Greene’s {}_4F_3 conceal the usual fourth trivial character in Katz’s balanced notation. With the project’s Gauss-sum convention, define over \mathbf F_q, q=p^f,

This is the trace convention for the Tate-normalized balanced object

The opposite Katz ordering exchanges the two character lists and sends \lambda to \lambda^{-1}; at \lambda=1 it gives the same stalk.

The factor q^{-2} in (1.1) is important. The raw balanced (4,4) convolution has Weil weight 7; the Tate twist by 2 lowers the weight by 4, giving the motivic normalization of weight 3.

For A\ne\varepsilon, the two character multisets are disjoint, hence the generic hypergeometric system is irreducible of rank four. The inverse-stable parameter set makes it self-dual with alternating polarization. In the Beukers–Heckman/Katz classification it is the primitive symplectic rank-four case; its connected geometric monodromy is \mathrm{Sp}_4. The local data below are enough for all later arguments even if one chooses not to invoke the full monodromy-group classification.

## 1.1 Classical differential operator

Write

The corresponding complex hypergeometric equation is

Its local exponent data are

At 0 the monodromy is maximally unipotent. For A^2\ne1, the monodromy at infinity has two size-two Jordan blocks, one with eigenvalue A and one with eigenvalue A^{-1}. When A is quadratic, the two eigenvalues coincide and the infinity block becomes a single size-four block. At 1 the monodromy is a nontrivial unipotent pseudoreflection.

## 1.2 Hodge numbers

In Fedorov’s convention take

Every \beta_k is strictly positive, including a=1/2. Hence

Therefore the generic Hodge polynomial is

After the standard shift, this is a pure weight-three variation with

This computation is independent of the order of A, as long as A is nontrivial.

## 1.3 Field of definition and low-order characters

Let

The single rank-four motive is naturally defined over the real cyclotomic coefficient field

It is not, for high M, one fixed rank-four motive over \mathbf Q. To obtain a \mathbf Q-motive one must take the Galois orbit or restriction of scalars, whose rank grows with [E_M:\mathbf Q]. This is another reason the collection indexed by all r is not a fixed bounded-rank \mathbf Q-family.

The low-order cases are as follows.

- M=1: A=\varepsilon. The upper and lower parameters cancel. The rank-four description degenerates completely, and the endpoint values r=0,N require the separate corrections already present in the verified formula.

- M=2: A=\phi. The system remains irreducible of rank four; its generic parameters are (1/2,1/2,1/2,1/2;1,1,1). At \lambda=1 the weight-three rank-two quotient is the classical weight-four, level-eight modular motive, and

- M=3,4,6: the coefficient field is still \mathbf Q; there is no rank drop. The infinity monodromy has the two size-two blocks described above.

- Higher M: the same rank, weight, Hodge numbers, and conifold behavior persist, but the coefficient field varies with M.

No low-order exception changes the slope-zero calculation below. The only true parameter cancellation is A=\varepsilon.

---

# 2. The point \lambda=1: conifold nearby cycles, not a smooth rank-four fiber

The exponent list 0,1,1,2 at 1 means that

The generic rank-four variation has a limit mixed Hodge structure whose monodromy filtration, centered at weight 3, has

with

The middle-extension stalk is

so its rank is three. Its boundary Hodge polynomial is

Equivalently, the boundary object consists of

- one Tate line of type (1,1) and weight 2;

- one rank-two pure weight-three quotient of types (3,0)+(0,3).

Thus it is inaccurate to call the value at 1 a pure rank-four motive. The finite-field hypergeometric value is the trace of this rank-three middle-extension/nearby-cycle object.

## 2.1 Exact Frobenius factor

The project’s arrangement model gives the factor without guessing from numerics. The exceptional divisor at the unique non-normal-crossing point contributes the Tate eigenvalue p. An involution sending the Kummer character A to A^{-1}, combined with Poincaré duality, gives an alternating perfect pairing on the remaining rank-two quotient with multiplier p^3. Hence that quotient has determinant p^3.

If S_{p,r} is the total effective trace, the quotient trace is S_{p,r}-p. Therefore

Set

The verified Jacobi/Greene congruence gives

and therefore

This is the precise bridge from the hypergeometric trace to the Apéry residue.

---

# 3. Exact ordinarity statement

Let

normalized by v_{\mathfrak p}(p)=1. The Newton polygon of

gives the two quotient slopes

Adding the Tate slope 1 gives

The Hodge slopes of the boundary object are (0,1,3). Hence it is ordinary exactly when a_{p,r} is a unit. By (2.4),

The first Hasse invariant is the reduction of Frobenius on the unique Hodge-slope-zero line. In the effective normalization it is

So the answer to Question (2) is yes, with two qualifications:

1. It is the Hasse invariant of the rank-three middle-extension object at the conifold, or equivalently of its rank-two weight-three quotient—not of a smooth rank-four fiber at 1.

1. The equivalence gives no simplification: the Hasse scalar is exactly the original Apéry residue.

There can never be two unit roots: the boundary Hodge polynomial has only one slope-zero slot. Thus the hoped-for mechanism “one explicit unit root plus one moving unit root, with zeros arising from collisions” is structurally impossible here.

---

# 4. The character parameter is a finite étale scheme, not a fixed curve

For fixed p, the tame character group is represented over \mathbf F_p by

Because N is prime to p, this is split finite étale of length N. After choosing a generator g of \mathbf F_p^\times, the point

corresponds to A_r.

One may form the universal Kummer twist over the finite product U\times\mathscr X_p, push forward the relevant cohomology, and take Frobenius on the unique Hodge-slope-zero line. This produces a genuine Hasse section

with

Thus

This is a correct algebraic zero-locus statement. It does not imply a useful degree bound.

## 4.1 Explicit polynomial representative

There is a unique representative of degree below N. Define the inverse Fourier coefficients

and put

Then character orthogonality gives

If one uses the fixed trace-function/Mellin representation

then, writing x=g^m, the same section is simply

Formula (4.6) is the concrete universal Hasse section requested in the question.

## 4.2 Why algebraicity gives no bound

The scheme \mathscr X_p itself has degree p-1. Every function on its p-1 points is regular. A nonzero regular function can vanish on p-2 points. For example,

vanishes at every nontrivial $N$th root of unity.

Moreover \operatorname{Pic}(\mathscr X_p)=0, so presenting the Hasse invariant as a section of a line bundle does not add a geometric degree. The length of its zero scheme is the quantity one is trying to bound.

The known reflection

implies that H_p^{\mathrm{char}} is reciprocal and descends to the inversion quotient through Y=Z+Z^{-1}. This can reduce a natural representative degree from N-1 to about N/2, but that remains linear in p.

Therefore a theorem of the following form would indeed finish the vertical problem:

There exist a unit u_p on \mathscr X_p and a Laurent polynomial Q_p of degree D_p=o(p)—ideally O(1) or O(\log p)—such that

operatorname{Ha}_p(A_r)=u_p(A_r)Q_p(A_r).

It would give

But this is a new Apéry-specific theorem, not a consequence of hypergeometric rank, conductor, purity, or ordinarity. The project already found raw interpolation degree p-2 at p=11,13, so bounded degree is not automatic.

Also, CFVZ’s factorization

lives in the generating/fiber variable T. It does not factor the character-space section H_p^{\mathrm{char}}(Z), which is a discrete Fourier transform of the coefficient vector.

---

# 5. Gross–Koblitz and Stickelberger: the exact calculation is circular, not compressive

The verified Jacobi formula can be reindexed into the equivalent form

Let T_{r,k} denote the product of the two squared Jacobi sums. Gross–Koblitz and the Jacobi factorization

produce the exact Stickelberger valuation

Thus the valuation-zero indices are

and their number is

For a generic character this is \asymp p, not O(1).

After reduction modulo \mathfrak p, the unit terms are exactly the Apéry summands. If

then

and hence

Equation (5.4) is strong confirmation that the Jacobi formula is normalized correctly. It is also the obstruction: cancellation among the valuation-zero Jacobi terms is literally the original Apéry cancellation.

The p-adic gamma expression therefore has O(p) moving fractional-part breakpoints and O(p) unit summands. It does not yield a polynomial of bounded degree in a Teichmüller scalar. Any claimed Gross–Koblitz compression must explain an additional cancellation or factorization beyond (5.2)–(5.4).

A general counterexample shows that bounded sheaf complexity cannot supply it. Fix a character \xi of order d. The rank-one Kummer sheaf with Mellin values

has fixed rank and conductor, but Stickelberger gives p\mid M_j for a positive proportion—about 1/d—of all tame characters. Thus fixed cohomological complexity is compatible with \gg p nonordinary twists.

---

# 6. Why the usual supersingular-divisor analogy does not give bounded degree

The statement “the nonordinary locus in a family over a curve is a finite divisor” is true for each fixed characteristic. The phrase “bounded degree” must be qualified.

For a family of elliptic curves over a curve C in characteristic p, the Hasse invariant is a section of

Its divisor has degree

which generally grows linearly with p.

The Legendre family is the elementary test. Let

The Hasse invariant is the coefficient of x^{p-1} in

Direct expansion gives

which has degree

So even a fixed non-isotrivial rank-two family with full \mathrm{SL}_2 monodromy has a Hasse polynomial of degree \Theta(p).

Function-field Lang–Trotter or supersingular-count theorems often fix the characteristic p and let the constant-field extension grow. In that regime the Hasse divisor is fixed. Here p itself varies, and the character parameter scheme has length p-1. That is a different quantifier order.

Our situation is less favorable than a fixed curve family:

1. the geometric point \lambda is frozen at the singular value 1;

1. the local monodromy character A_r changes with r;

1. the coefficient field E_{\operatorname{ord}(A_r)} changes;

1. the parameter space \mu_{p-1} grows with p.

Therefore no fixed-curve divisor theorem supplies a uniform O(1) degree.

---

# 7. Precise theorem that can be banked

For every prime p\ge5 and every nontrivial character A_r, the following statements hold in the project’s normalization.

## Theorem

Let

Then:

1. the generic \lambda-system is the rank-four, weight-three, symplectically self-dual hypergeometric variation with Hodge polynomial 1+T+T^2+T^3;

1. the middle-extension stalk at \lambda=1 has Hodge polynomial 1+T+T^3 and Frobenius polynomial

1. at the Teichmüller prime,

1. consequently,

1. the resulting Hasse invariant over the character scheme is the regular function

The theorem identifies the exact crystalline meaning of Z_p. It does not estimate |Z_p|.

---

# 8. Gaps and the actual missing lemma

[GAP-1: integral nearby-cycle comparison in publication-ready form.] The project’s arrangement computation proves the semisimplified factor (2.3), and the congruence (2.4) is independently verified by the Jacobi and Greene formulas. For a fully external arithmetic-geometry paper, the comparison between the \ell-adic middle-extension stalk at the conifold and a canonical integral overconvergent F-crystal should be written with all lattice and normalization choices explicit. This is the technically delicate bridge underlying the phrase “the Hasse invariant is b_r.” It does not affect the finite-field congruence or the Newton-polygon conclusion once (2.3) is accepted.

[GAP-2: low-order component-group audit.] The connected symplectic classification is clear from the primitive reciprocal hypergeometric data. I have not reproduced Katz’s complete finite component-group/exceptional-list check for every small order M in this report. This does not affect rank, Hodge numbers, the conifold factor, or the Hasse equivalence.

[GAP-3: Apéry-specific tame-Hasse compression.] The decisive missing statement would be a nontrivial factorization

on \mu_{p-1} with u_p nowhere zero and

uniformly. Degree O(\log p) would give the empirically expected scale. Nothing in Katz’s generic classification, BCM point-count realization, purity, self-duality, or Gross–Koblitz provides this. Equations (5.2)–(5.4) show that the naive Jacobi/Gamma route simply reconstructs the original Apéry sum.

Accordingly, the honest status is:

The current 8p^{3/4} bound remains the unconditional record in the project.

---

# Least-confident step

The least-confident step is not the Hodge-number calculation or the Stickelberger obstruction; those are explicit. It is the integral nearby-cycle normalization at the singular point \lambda=1: promoting the semisimplified \ell-adic factor

to a canonical integral crystalline object whose first Hasse map is literally the chosen reduction of a_{p,r}. The repo’s arrangement calculation and three independent finite-field formulas make the conclusion highly constrained, but this is the place where a referee will demand the most careful comparison theorem and lattice bookkeeping.

---

# References and project anchors

- N. M. Katz, Exponential Sums and Differential Equations, Annals of Mathematics Studies 124, Princeton University Press, 1990.

- F. Beukers and G. Heckman, “Monodromy for the hypergeometric function {}_nF_{n-1},” Inventiones Mathematicae 95 (1989), 325–354.

- R. Fedorov, “Variations of Hodge Structures for Hypergeometric Differential Operators and Parabolic Higgs Bundles,” IMRN 2018, 5583–5608. arXiv

- F. Beukers, H. Cohen, A. Mellit, “Finite hypergeometric functions,” Pure and Applied Mathematics Quarterly 11 (2015), 559–589. arXiv

- Project commit f438bd9bd535368104e389c5039f95f55c238af6: independent confirmation of the Jacobi/Greene formula.

- Project file problems/3.2/lemmaBprime_result.tex at commit 42f4ec1b324029f55593f3d8b444c764a1aaa6c4: Hodge calculation, conifold factor, Newton slopes, and the G1/Hasse-tautology verdict.

- Project commit 64227dc48543afe17660c309dcb1fdabd787bc2e: exact Stickelberger valuation and proof that the unit Jacobi terms reconstruct the Apéry sum.

- Project commits 230d3af3a48bc8ec007e49b423d39190a5f7b3af and 1bb78487556800ed81c4c126dbbe5e96e891e704: bounded-Hasse criterion, Jacobi counterexample, and the distinction between fiberwise and character-twisted Hasse invariants.