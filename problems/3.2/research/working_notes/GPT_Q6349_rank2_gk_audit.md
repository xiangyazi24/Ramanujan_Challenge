ANSWER Q6349 3f94d009

# Executive verdict

The proposed implication

is false as a general principle.

Rank two controls the dimension of the space of horizontal sections and the number of Frobenius eigenvalues at a fixed geometric fiber. It does not control the complexity of the (j)-th Taylor jet of a chosen local section after nonlinear algebraic reversion. Gross--Koblitz converts the multiplicative Mellin coefficients of a finite-field hypergeometric trace into bounded products of (p)-adic Gamma values; it does not collapse the inverse Mellin transform, and it does not commute with Lagrange coefficient extraction.

For the present Apéry pullback, the strongest conclusions I can derive are:

1. There is a clean exact Lagrange formula, and it simplifies substantially:

1. Substituting the Franel expansion gives exact finite hypergeometric sums. Modulo (p), at the quarter indices, every summand is a bounded product of Morita Gamma values, equivalently a bounded product of Jacobi sums. The number of summands still grows polynomially with (n): in the most direct Gamma/Jacobi lattice formula it is (Theta(n^3)).

1. By expanding the Jacobi sums and summing the triangular index set geometrically, one can compress the formula to a bounded-dimensional complete character sum with finitely many lower-dimensional boundary corrections. This is a real compression, but it is not an (O(1))-term Jacobi-monomial evaluation. Moreover, the equality is a (p)-adic/reduction statement; a complex Weil--Deligne estimate for the complete sum does not by itself decide whether the reduction is zero modulo (p).

1. At the special quarter indices, some powers become low-order Kummer characters. This makes a CM/Cartier evaluation plausible and is consistent with the observed quarter-zero law. But one still needs an additional identity converting the local coefficient functional into a finite-field fiber trace or a controlled Cartier matrix entry. Rank two alone does not supply that identity. This is the precise remaining opening, marked [GAP-CARTIER] below.

1. Replacing (b_r=(tau*tau)_r) by the two rank-2 formulas does not create a new cross-prime Deligne object. It expands and then recombines to the already-known terminating

There are bounded-dimensional character-sum presentations, but the condition (b_requiv0pmod p) is a (p)-adic divisibility condition on a trace value, not a complex cancellation estimate. Applying an additive character to the value of a trace is not a tensor operation on the original rank-2 sheaf.

So the skeptical verdict is:

A special bounded evaluation at the quarter point remains possible, but it would be a new arithmetic identity, not a formal consequence of rank two or of the quadratic pullback.

# 1. Exact reversion algebra

Put

so that (t=phi(x)=x/Phi(x)) and the inverse branch satisfies

The inverse is explicitly quadratic:

hence

This degree-two algebraicity explains why the two square-root branches and the polynomial (q(t)=1-34t+t^2) occur. It does not imply that coefficients after composition are hypergeometric terms.

For any formal series (G(x)), Lagrange--Bürmann gives

An integration-by-parts form is more useful. Since

we obtain

A direct calculation gives

Now use

Substitution into (1.3) proves (0.1). This cancellation is the main exact simplification contributed by the degree-two pullback.

# 2. Exact finite sums

The hypergeometric representation of (h) gives

The second equality follows from

Alternatively,

## 2.1 A one-layer coefficient kernel

Define

Then

It is also one terminating Gauss function:

Here the symbol ({}_2F_1) compresses a length-(L+1) terminating sum; it is not a Gamma-product evaluation.

Using (2.2) in (0.1) gives

Let

Then

Equations (2.5)--(2.6) are already enough to refute the naive implication “rank two means one Gamma monomial.” They are exact finite convolutions whose length grows with (n).

## 2.2 Appell form

For completeness, define

Then

and equivalently

Thus nonlinear reversion naturally produces a terminating Appell-type object. The differential equation may have rank two in its geometric parameter, while the coefficient extraction has acquired a second summation direction.

# 3. Exact quarter-index formulas modulo (p)

Let (p>3) and suppose (nle(p-3)/4) or (n=(p-1)/4), as in the quarter-point applications. Set

These are the representatives in ([0,p-1]) of (n-1/2) and (n+1/2) modulo (p). Since generalized binomial coefficients are polynomial in their top parameter for lower index (<p),

In the quarter range every factorial occurring below is (<p). Expanding (f_m) and (K_L) therefore gives the completely elementary formula

Similarly,

The number of lattice terms in (3.3) is

The number in (3.4) is

Thus at a quarter index (nasymp p), the direct bounded-product Gamma/Jacobi expansion has (Theta(p^3)) terms.

# 4. Morita Gamma reformulation

Use Morita's convention

For (0le Ble A<p),

This is an exact identity in (mathbf Z_p), not merely a congruence, because no factorial contains (p).

Define

Then (3.3) becomes

with the same triangular ranges as (3.3). Formula (3.4) has the identical replacement.

This is a valid rank-independent (p)-adic Gamma formula. Each summand contains only a bounded number of Gamma factors, and every argument is affine in the indices. The unresolved operation is the growing sum over (m,k,a).

This distinction matters:

- “a bounded product of Gamma values per summand” is immediate factorial bookkeeping;

- “a bounded number of Gamma monomials for the whole coefficient” is a strong evaluation theorem and is not implied by the first statement.

# 5. Jacobi-sum reformulation

Let (omega:mathbf F_p^timestomu_{p-1}) be the Teichmüller character, reduced at a fixed prime above (p), and define

For (0<B<A<p-1), direct reduction of the Jacobi sum gives

Indeed, after reduction, only the exponent (p-1) survives in the sum over (xinmathbf F_p), and its coefficient is (-binom AB). Boundary cases (B=0,A) are treated separately.

For the last binomial in (3.3), the second character is especially simple:

Hence every interior summand of (3.3) is a product of five Jacobi sums:

- three copies arising from (binom mk^3);

- one from (binom{A_sigma}{a});

- one from (5.2).

Gross--Koblitz then rewrites each Gauss sum occurring in these Jacobi sums as

Together with (J(A,B)=g(A)g(B)/g(AB)) in the nondegenerate cases, this is the exact Gross--Koblitz reformulation.

Again, (5.3) converts each Mellin coefficient into a bounded Gamma product. It does not evaluate the sum over the growing triangular region.

# 6. Fiber traces versus coefficients

This is the conceptual center of the audit.

## 6.1 Finite-field hypergeometric traces

Greene's normalization is

and

When cubic characters exist, the finite-field hypergeometric object associated with the parameters (1/3,2/3;1) is obtained by taking (A=chi_3), (B=chi_3^2), (C=varepsilon), up to the normalization and twist dictated by the chosen geometric model.

Equation (6.1) already displays the correct lesson:

- the Mellin transform has a bounded product of two Jacobi sums;

- the inverse Mellin transform has (p-1) character terms.

The sheaf has rank two because it has two numerator characters, not because its trace is a sum of two Jacobi monomials.

For (pequiv2pmod3), there is no cubic character on (mathbf F_p). One must either pass to (mathbf F_{p^2}), use a descent, or use the (p)-adic hypergeometric Gamma-sum formalism. None of these operations removes the inverse Mellin summation automatically.

## 6.2 Taylor coefficients after algebraic reversion

The coefficient (tau_j) is the (j)-th local jet at (t=0) of

Lagrange inversion inserts (Phi(x)^j), so the (j)-th coefficient is a growing moment of the original local solution. It is not the Frobenius trace of the rank-2 family at one point (xinmathbf F_p).

Rank is preserved under algebraic pullback of the differential equation. Jet complexity is not.

A useful exact constant-term form is

To replace (6.2) by a complete finite-field trace sum, one must replace the local truncated expansion by a globally defined finite-field function without changing this constant term. Agreement of rational functions on (mathbf F_p)-points is not enough: formal Taylor coefficients are not determined by values on (mathbf F_p).

This missing comparison is exactly:

[GAP-CARTIER] Construct a bounded-conductor Frobenius/Cartier object whose selected matrix entry is the coefficient in (6.2), including the nonlinear reversion and the correct local branch.

Without [GAP-CARTIER], the rank-2 trace and the quarter coefficient are different arithmetic functionals.

## 6.3 Morita Gamma reductions

Morita Gamma formulas such as (4.1) are exact factorial identities in the safe range. Gross--Koblitz formulas such as (5.3) concern Gamma values at fractions (c/(p-1)), after a choice of Gauss sum and uniformizer.

Passing from (4.1) to (5.3) is legitimate term by term. It does not eliminate the summation indices. Confusing these two stages is the central logical error in the proposed “rank-2 horizontalization.”

# 7. What the quadratic pullback can and cannot compress

The degree-two inverse (1.1) gives two sheets. The second sheet tends to (x=1/8) when (tto0); it is not another copy of the same local Taylor branch at (x=0). Algebraic conjugation therefore does not express the desired coefficient as the sum of two local hypergeometric coefficients.

The pullback does provide:

- the order-two recurrence;

- the discriminant (1-34t+t^2);

- the exact Lagrange cancellation (0.1);

- the finite (K_L), Appell, Gamma, and Jacobi formulas above.

It does not provide an (O(1))-term evaluation of (K_L), (C_L), or the outer Franel convolution.

At a special index such as (n=(p-1)/4), some endpoint powers become quartic Kummer characters on (mathbf F_p^times). This is extra arithmetic structure, not a rank argument. For (pequiv3pmod4), there is no quartic character over (mathbf F_p), so any quarter evaluation in the (pequiv23pmod{24}) class must involve an extension-field character or a descent identity. A direct (mathbf F_p) quartic-Jacobi formula is therefore invalid in that class.

The observed quarter-zero law may well come from an exceptional CM/Cartier evaluation. The present derivation reduces the problem to proving such an exceptional evaluation; it does not derive it.

# 8. A bounded-dimensional complete sum does exist

Although the Jacobi-monomial sum is long, its index polytope has fixed dimension. This permits a different compression.

Consider the interior part of (3.3). Expand the five Jacobi sums using variables

The dependence on (m,k,a) is through three multiplicative monomials

The triangular index sum becomes

It has the closed form

with the evident limiting interpretation on the denominator strata.

Consequently the principal interior block of (sigma_n) is one complete five-variable character sum with kernel (8.3), plus a bounded number of lower-dimensional boundary sums. The (tau_n) formula is a sum of the same construction for (n,n-1,n-2).

This is the strongest genuine compression I see:

But this does not prove the mission claim:

- it is not an (O(1))-term product/short sum of Jacobi sums;

- the rational kernel has singular strata requiring separate treatment;

- the equality relevant to (tau_nbmod p) is (p)-adic/reductive;

- a complex Weil bound for the complete sum does not determine its reduction modulo the chosen prime above (p).

A successful rank-2 Gross--Koblitz proof would have to do more: it would need a (p)-adic cohomological evaluation or stationary-phase argument that collapses (8.3) at the quarter index.

# 9. What to check at (p=13) and (p=29)

The following checks are exact and discriminating.

## 9.1 Direct coefficients

From the recurrences,

Therefore

For (p=29), the tau recurrence gives

Thus

## 9.2 Formula verification

At (p=13,n=3):

- formula (3.3) has (T(3)=20) lattice terms;

- formula (3.4) has (T(3)+T(2)+T(1)=34) terms.

At (p=29,n=7):

- the sigma formula has (T(7)=120) terms;

- the tau formula has

terms.

A machine audit should verify, in this order:

```plain text
1. Construct b_m and f_m exactly.
2. Revert t=x(1-8x)/(1+x) as a formal series.
3. Verify both pullback identities through degree n.
4. Compute tau_n and sigma_n by square-root convolution.
5. Verify (0.1), then (2.5)-(2.6), then (3.3)-(3.4).
6. Replace every binomial by the Morita Gamma ratio (4.1).
7. Verify the Jacobi congruence (5.1) for every interior binomial.
8. Expand the Jacobi sums and verify the five-variable formula plus boundary strata.
9. At p=29, print all nonzero term contributions to confirm that tau_7=0 is a cancellation, not a termwise Gamma zero.
```

## 9.3 What a finite search cannot prove

A search at (p=13,29) can refute a specific proposed formula. It cannot prove that no bounded short formula exists. With (p)-dependent constants, any finite data set can be fitted.

# 10. What would constitute a proof of impossibility?

One must first define the allowed formula class.

## 10.1 Characteristic-zero hypergeometric terms

For a claim that (tau_n) or (sigma_n) is a fixed sum of classical Gamma/Pochhammer monomials, a proof of impossibility would be an Ore-factorization certificate:

- derive the second-order recurrence operator;

- prove it has no first-order right factor over (overline{mathbf Q}(n));

- equivalently, prove its rational Riccati equation has no rational solution.

That would rule out a basis of hypergeometric terms and fixed two-term decompositions in the classical sense.

I do not use this as a load-bearing step here; the present audit is about the stronger, (p)-dependent Jacobi/Gross--Koblitz claim.

## 10.2 Uniform finite-field Jacobi formulas

To rule out uniformly an expression made from at most (K) Jacobi monomials with character exponents affine in (j), one needs a structural invariant. A plausible route is:

1. formalize the ansatz as a bounded-conductor Mellin sheaf on the character torus;

1. prove that the actual coefficient-index function (jmapstotau_jbmod p) has Mellin conductor, singular support, or geometric monodromy growing with (p);

1. derive a contradiction for fixed (K).

A large interpolation degree is not enough: bounded-conductor trace functions can have full interpolation degree as ordinary polynomials on (mathbf F_p).

No such conductor lower bound is proved here. [GAP-IMPOSSIBILITY]

## 10.3 The special quarter coefficient

Even a uniform impossibility theorem for generic (j) would not exclude an exceptional formula at (j=(p-1)/4) or (j=(p-3)/4). Special CM points often admit Gauss/Jacobi evaluations that generic points do not.

Therefore the quarter-point problem remains open to a special evaluation theorem. Failed bounded-ansatz searches are only evidence.

# 11. Returning to (b_r=(tau*tau)_r)

The direct coefficient formula is already

For (r<p), every summand is a bounded Morita-Gamma product and, away from boundary cases, a bounded Jacobi product. The sum still has (r+1) terms.

Expanding the Jacobi factors and summing (k) geometrically yields a bounded-dimensional complete character sum. Thus bounded-dimensional presentations are not new consequences of splitting (F=tau^2); they are already implicit in the terminating ({}_4F_3).

Substituting the two Lagrange formulas for (tau_jtau_{r-j}) makes the intermediate expression much larger and then collapses algebraically back to (11.1). It does not lower the arithmetic complexity.

# 12. Why this is not yet a Weil--Deligne cross-prime tool

There are four separate obstructions.

1. Complex size versus (p)-adic zero. Deligne bounds the complex absolute values of trace sums. The event (b_requiv0pmod p) is divisibility at a chosen (p)-adic place. A square-root complex bound does not decide this divisibility unless one also has a sufficiently small global norm/height statement.

1. The character parameter moves. In a finite-field hypergeometric representation of (11.1), the characters depend on (r). In the horizontal problem (r=nbmod p), so both the field and the Mellin parameter move with (p).

1. Additive character of a trace value is nonlinear. Even if (b_r) is itself a trace of a bounded-rank sheaf, the phase

is an additive character applied to the numerical trace value. This is not obtained by tensoring the original sheaf with a fixed Artin--Schreier sheaf unless (b_r) is a bounded-degree regular function of the geometric parameter.

1. Convolution does not preserve the needed horizontal form. The identity (b=tau*tau) is coefficient convolution. Sheaf convolution applies to trace functions in a field variable; it does not automatically identify coefficient convolution after reversion with a bounded-conductor sheaf in the index (r).

A genuine positive result would require a uniform Mellin-sheaf realization of the coefficient-index function (rmapstotau_r), with controlled conductor as (p) varies, followed by a theorem transporting convolution and the moving specialization (r=nbmod p). [GAP-HORIZONTAL-SHEAF]

# 13. Final boundary map

The exact status is:

- Valid: rank-2 finite-field hypergeometric fiber traces have bounded Gauss-product Mellin transforms.

- Valid: Gross--Koblitz rewrites each Mellin coefficient as a bounded product of (p)-adic Gamma values.

- Valid: the reverted coefficients admit exact length-growing Gamma/Jacobi sums, equations (3.3)--(5.3).

- Valid: those long sums can be reorganized as bounded-dimensional complete sums, equation (8.3), with boundary corrections.

- Invalid implication: rank two alone does not turn a Taylor coefficient after nonlinear reversion into (O(1)) Gamma/Jacobi monomials.

- Open special possibility: the quarter coefficient may satisfy an exceptional CM/Cartier evaluation. [GAP-CARTIER]

- No present cross-prime consequence: the convolution formula does not yield a Deligne estimate for the moving divisibility condition (pmid b_{nbmod p}).

The strongest honest research target is therefore not “apply Gross--Koblitz because the system has rank two.” It is:

Construct and evaluate the specific Cartier/Mellin functional defined by (6.2), at the quarter character, including the algebraic reversion and the local branch. Show that this functional reduces to a bounded Jacobi expression, or prove that its Mellin conductor grows.

That is a precise theorem-shaped fork. Everything before it is exact; everything after it is currently [GAP].

# Bibliographic anchors checked

- Benedict H. Gross and Neal Koblitz, Gauss sums and the p-adic Gamma-function, Annals of Mathematics 109 (1979), 569--581, DOI 10.2307/1971226.

- John Greene, Hypergeometric functions over finite fields, Transactions of the American Mathematical Society 301 (1987), 77--101, DOI 10.1090/S0002-9947-1987-0879564-8.