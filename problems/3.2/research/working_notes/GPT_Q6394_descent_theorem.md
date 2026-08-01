ANSWER Q6394 34452441

# Descent of the Franel symmetric square and the quadratic companion

## Executive verdict

The geometric descent does work, but one must formulate it in the correct category and at the correct rank.

- The Franel period is carried by a rank-two elliptic Gauss–Manin system on the x-line.

- Its literal tensor square has rank four, but it splits as

- On the quadratic cover

- The + descent is the transcendental Apéry rank-three system. The - descent is exactly its quadratic companion.

- Consequently, for every good prime p,

This gives the requested theorem:

Theorem (quadratic companion as a bounded-conductor Mellin object). There is an integral compatible rank-three system \mathcal G_2=\mathcal A_- on \mathbf G_m with absolute conductor bounded by 11 such that, for 1\le r\le p-2,

T(r)=sum_{tinmathbf F_p^times}chi_2(q(t))A_p(t)omega^{-r}(t)

equiv

sum_{tinmathbf F_p^times}

operatorname{tr}(operatorname{Frob}_tmidmathcal G_2)omega^{-r}(t)

pmod{mathfrak p}.

For a generic character \chi, the right side is minus the Frobenius trace on a four-dimensional space

H_c^1(mathbf G_m,mathcal G_2otimesmathcal L_chi).

Even better, the two-term decomposition collapses:

The first Franel sum is the Mellin transform of this rank-six pushforward, while T(r) is the Mellin transform of \mathcal A_-. Therefore

For generic r, this is the reduction of a rank-two Mellin cohomology trace. Thus the apparent “rank 4+4” object reduces, modulo p, to rank 3+3, and the final difference reduces to one rank-three input sheaf whose generic Mellin transform has dimension two.

The geometric monodromy of \mathcal A_+ is O_3 with identity component SO_3; its determinant is \mathcal K_q. The twist \mathcal A_- has connected geometric monodromy SO_3. The direct sum does not have product monodromy: its image is the diagonal copy

What remains open is not descent. It is the defining-characteristic problem already isolated in Q6371: Katz’s complex Mellin equidistribution controls normalized complex Frobenius traces, but does not bound the number of characters for which their integral values reduce to zero modulo the same prime p.

---

## 1. The quadratic cover

The equation \phi(x)=t is

Thus the two roots satisfy

and the deck involution is

The discriminant is

Let \alpha,\beta be the two roots of q. Put

After deleting the two ramification points upstairs, \pi:V\to U is a finite étale Galois cover with group C_2.

The associated quadratic rank-one sheaf is

characterized by

Its geometric local monodromy is

- trivial at t=0;

- -1 at t=\alpha;

- -1 at t=\beta;

- trivial at t=\infty.

The degree of q is even, so there is no geometric ramification at infinity.

---

## 2. An explicit Franel elliptic family

Let

Consider the toric family

in (\mathbf G_m)^2, with its smooth projective genus-one compactification. The period at x=0 is

because

Let

be the rank-two Gauss–Manin sheaf of this elliptic family. Its Picard–Fuchs operator is

Equivalently,

The singular points are

The Riemann scheme is

and the local monodromy at each point is a nontrivial unipotent J_2. The family is non-isotrivial, hence

### Hasse invariant

The invariant differential on C_x is represented by the toric residue. Its Cartier coefficient is

Using

we get

If

then, with this differential normalization,

Therefore

This is why the correct effective object has rank three rather than four.

---

## 3. The actual descent datum

Put

The characteristic-zero transformation

is the scalar-period shadow of the symmetric-square pullback identity

At the differential-module level this can be verified without modular forms. Let

The involution satisfies

A direct operator calculation gives

Thus multiplication by g gives an involutive descent datum on the rank-three connection.

The arithmetic realization comes from the modular interpretation of the Franel family: the involution \iota is the index-two modular/Atkin–Lehner involution underlying the Apéry–Franel transformation. Its isogeny correspondence acts on R^1g_*\mathbf Q_\ell; after taking \operatorname{Sym}^2 and dividing by the isogeny degree, the correspondence squares to the identity. Hence it supplies a Frobenius-compatible étale and crystalline descent datum.

[GAP-1: explicit isogeny formula]. The existence of the modular isogeny correspondence is standard in the modular parametrization used for the Apéry–Franel identity, but this report does not write a Weierstrass model and the rational maps of the isogeny. For a completely stand-alone algebraic proof, that map should be printed and checked. The connection-level cocycle and all monodromy consequences below are independent of this missing display.

Given the equivariant structure, define

Then

The projection formula and

give

This is the desired invariant/anti-invariant decomposition.

---

## 4. Local monodromy on the t-line

The Apéry differential operator is

Its Riemann scheme is

The Fuchs sum is 6, as required for an order-three equation with four regular singularities.

The ramified pullback doubles the exponents at \alpha,\beta, so the pullback is regular there. This forces the local monodromy at each root of q to have order two and hence to be semisimple.

### The invariant descent \mathcal A_+

Its geometric local monodromy is

Therefore

The identity component of its geometric monodromy is

and the two reflections enlarge it to

### The anti-invariant descent \mathcal A_-

Twisting by \mathcal K_q changes only the two branch-point monodromies:

Its determinant is trivial and

All monodromy is tame. The sums of the local Artin drops are

and

Hence, under the convention \operatorname{cond}=\operatorname{rank}+\sum a_x, one has

Under Katz’s rank-plus-number-of-singularities convention, both are at most 7.

---

## 5. The trace congruences

The exact polynomial identity

shows that the symmetric-square Hasse invariant is deck invariant and descends. Since pullback of functions along the dominant map \pi is injective, the descended Hasse invariant is exactly A_p(t).

For the compatible crystalline/étale system this gives

Twisting gives

At q(t)=0, the middle-extension stalk of \mathcal K_q has trace zero, agreeing with the convention \chi_2(0)=0.

Equation (5.2) proves the requested bounded-conductor interpretation of the companion.

Category warning. Equations (5.1)–(5.2) are congruences after reducing an integral compatible-system trace modulo a prime above p. They are not equalities between A_p(t) and a complex trace. This distinction is exactly why complex equidistribution does not solve the divisibility problem.

---

## 6. Mellin transforms and the exact cohomological ranks

Let j:U\hookrightarrow\mathbf G_m. For a multiplicative character \chi, write \mathcal L_\chi for its Kummer sheaf.

The Grothendieck trace formula gives, for generic \chi,

There is no H_c^0 or H_c^2 for generic \chi.

On

one has

Middle extension across \alpha,\beta adds the dimensions of inertia invariants.

For \mathcal A_+:

so

For \mathcal A_-:

so

This explains the observed architecture:

- Franel pushforward: Mellin rank 2+4=6;

- quadratic companion T(r): Mellin rank 4;

- their difference, b_r: Mellin rank 2.

Indeed, if

then

is the reduction of the Mellin transform of

Since T(r) is the Mellin transform of \mathcal A_-,

is minus the Mellin transform of \mathcal A_+. Therefore

for all generic r, with finitely many low-order-character corrections handled by the usual H_c^0/H_c^2 terms.

This is the strongest non-circular descent statement: b_r is the defining-characteristic reduction of a two-dimensional Mellin cohomology trace of a fixed rank-three, bounded-conductor sheaf.

---

## 7. What happens to the rank-four tensor square

If one insists on the literal square of a rank-two trace, use

Since

its trace is p at an \mathbf F_p-point. The full pushforward decomposes as

These are the anticipated rank-4+4 invariant and anti-invariant objects. After reduction modulo p, both Tate summands disappear. Thus the effective mod-p decomposition is rank 3+3.

---

## 8. Monodromy of the combined object

Because

the direct-sum representation is

Therefore its image is a single copy of O_3, not O_3\times SO_3.

The virtual difference relevant to b_r is simply -\mathcal A_+. Hence its geometric monodromy is

This is an important correction to a possible “rank-eight product monodromy” interpretation: the two pieces are maximally related by a quadratic twist.

---

## 9. What Katz Mellin equidistribution gives

Let

The sheaf j_*\mathcal A_+[1] is a perverse Mellin object of generic fibre dimension two. Katz’s convolution formalism, and the more general arithmetic Fourier-transform formalism of Forey–Fresán–Kowalski, attach a reductive Tannakian group G_{\mathrm{Mell}} to it. After normalization by p^{3/2}, the conjugacy classes governing M_p(\chi) equidistribute in a maximal compact subgroup of G_{\mathrm{Mell}} as the finite field or extension degree varies, once the finite exceptional set of characters is removed.

The input local system has O_3 monodromy, but G_{\mathrm{Mell}} is the monodromy of the two-dimensional cohomology family, not the same O_3. Its standard representation has dimension two.

[GAP-2: exact Mellin Tannakian group]. The local monodromy and the absence of an obvious nontrivial Kummer self-twist strongly suggest that the connected Mellin group is SL_2, with a possible finite central determinant factor; equivalently, the compact distribution should be a central phase times an SU_2 trace distribution. A proof requires the Katz convolution irreducibility, determinant, and autoduality-sign calculation. This report does not assert the exact full group.

What is unconditional from the general theorem is complex equidistribution after the Tannakian group is computed. It predicts square-root-scale cancellation and bounded moments for normalized complex Mellin values.

What it does not give is

Reduction to zero modulo the defining prime is invisible to the archimedean conjugacy-class distribution. Thus the descent solves the object problem but leaves [GAP-DCM], defining-characteristic Mellin zero-density, unchanged.

---

## 10. Does this explain the census fact F_4=0?

Not by itself.

The descent explains three structural facts:

1. the quadratic companion is not an unrelated error term; it is the determinant twist of the same O_3 system;

1. after subtraction, the Apéry coefficient is governed by a two-dimensional Mellin cohomology space;

1. low-order characters are exactly where extra cohomology or self-twist phenomena may occur.

This makes the observed random-scale distribution plausible. It also explains why special low-order quarter-point laws can occur.

But F_4=0 is a statement that four different primes do not align on one moving integer column. The present theorem is vertical, one prime at a time. It supplies no cross-prime independence theorem. In particular, neither O_3 monodromy nor complex Mellin equidistribution forbids a defining-characteristic zero at a high-order character.

---

## 11. Machine-checkable predictions

### 11.1 Operator descent

Symbolically verify

and

Then compute the symmetric-square scalar operator of L_F and check

This is an exact rational-function identity over \mathbf Q.

### 11.2 Riemann schemes

Indicial calculations should print:

```plain text
Franel rank 2:
  x=0       [0,0]
  x=1/8     [0,0]
  x=-1      [0,0]
  x=infinity[1,1]

Apery rank 3:
  t=0       [0,0,0]
  t=alpha   [0,1/2,1]
  t=beta    [0,1/2,1]
  t=infinity[1,1,1]
```

### 11.3 Toric Hasse invariant

For several primes, compute

and verify coefficientwise that it equals H_p(x).

### 11.4 Point-count congruence

For nonsingular x\in\mathbf F_p, count the compactified curve C_x and verify

### 11.5 Companion trace

For t away from 0,\alpha,\beta,\infty, the predicted reductions are

An explicit K3/elliptic-correspondence implementation should verify this at split and nonsplit fibres.

### 11.6 Mellin degrees

Compute the local L-function of the character-twisted sums over \mathbf F_{p^n} for n=1,2,3,4. Newton reconstruction predicts:

- degree 2 for the \mathcal A_+ Mellin transform;

- degree 4 for the \mathcal A_- Mellin transform;

- degree 6 for the Franel pushforward.

This is the most discriminating same-day test of the descent architecture.

---

## 12. References actually used

- X. Caruso, F. Fürnsinn, D. Vargas-Montoya, W. Zudilin, Galois Groups of Apéry-like Series Modulo Primes, arXiv:2510.23298; especially the quadratic cover, its involution, the Franel differential equation, and the Apéry–Franel transformation.

- F. Beukers and C. A. M. Peters, A family of K3 surfaces and \zeta(3), J. reine angew. Math. 351 (1984), 42–54, for the Apéry K3 realization.

- N. Katz, Convolution and Equidistribution: Sato–Tate Theorems for Finite-Field Mellin Transforms, Annals of Mathematics Studies 180, 2012.

- A. Forey, J. Fresán, E. Kowalski, Arithmetic Fourier transforms over finite fields: generic vanishing, convolution, and equidistribution, arXiv:2109.11961.

---

## Least-confident step

The least-confident step is the arithmetic upgrade in Section 3 from the explicitly verified connection cocycle to a written Frobenius-compatible isogeny correspondence for the toric elliptic family. The modular interpretation makes the correspondence standard, and all downstream local-monodromy and Hasse-invariant formulas are consistent with it, but a referee-proof self-contained version should print the Weierstrass model of C_x, the isogeny C_x\to C_{\iota(x)}, and verify that its normalized symmetric square realizes the cocycle g(x)=\frac89(1+x)^2. This is precisely [GAP-1]; it is narrower than the original [GAP-DESCENT] and is directly machine-checkable.