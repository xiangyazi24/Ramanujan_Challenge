ANSWER Q8797 0dfff2df

# Q8797: Short-window reflection principle for the Apéry zero set

## Executive conclusion

The short-window reflection principle is a plausible Apéry-specific route to improving

\(|Z_p|=O(p^{2/3})\)

but it is not currently a consequence of the known recurrence, reflection symmetry, or square factorization. The principle is stronger than the available orbit information: it requires a new rigidity statement about the marked Apéry solution, not just the underlying projective dynamical system.

The most promising interpretation is not Lyapunov theory but arithmetic rigidity coming from the integral Apéry motive: the Beukers integral structure, Dwork congruences, and the K3 Picard-Fuchs realization. However, no currently known theorem converts those structures into the required four-zero local exclusion.

---

# 1. What distinguishes the true Apéry solution from phantom returns?

The transfer matrix recurrence alone is insufficient.

The recurrence

\[
v_{d+1}=M(d)v_d
\]

defines a projective orbit in \(\mathbf P^1(\mathbf F_p)\). The phantom construction shows that an arbitrary initial vector can produce many artificial returns to the zero point. Therefore the property cannot be:

- recurrence order,
- determinant one structure,
- projective orbit size,
- no-adjacent-zero condition.

The extra arithmetic data are likely one of:

1. The Apéry integrality condition.

The pair \((a_n,b_n)\) satisfies

\[
\frac{a_n}{b_n}\to \zeta(3)
\]

and the denominator cancellation is exceptional. A phantom solution does not correspond to a globally integral approximant.

2. Dwork congruence structure.

The congruences

\[
b_{qp+r}\equiv b_qb_r\pmod p
\]

and higher-depth congruences constrain the marked orbit across p-adic scales.

3. The K3 motive.

The Apéry recurrence is not an arbitrary recurrence; it is a Picard-Fuchs equation. The true solution is a period vector selected by geometry.

A possible theorem would therefore have the form:

"A short cluster of Frobenius zeros of the distinguished period vector forces a Hodge-theoretic symmetry."

Such a theorem is not presently available.

---

# 2. Why four short-window zeros should force a reflected pair

Assume zeros

\[
d_1<d_2<d_3<d_4
\]

lie in an interval of length \(\eta\sqrt p\).

Reflection sends

\[
d\mapsto p-1-d.
\]

A centered adjacent pair means two zeros satisfy

\[
d_i+d_j=p-1.
\]

The heuristic mechanism is:

- zeros are rare because the orbit intersects a codimension-one condition;
- reflection doubles every non-central zero;
- four returns in a very short interval are too many degrees of freedom unless one return is forced by the involution.

This resembles a geometric incidence statement, not a simple recurrence statement.

The K3 interpretation suggests that the relevant object is the intersection of a Frobenius orbit with a special divisor. The reflection involution comes from the functional equation. A proof would probably require a local monodromy or Picard lattice argument.

However, the current K3 information gives global symmetry, not local spacing.

---

# 3. Can the square factorization \(A_p=B_p^2\) prove it?

This is the most interesting algebraic route, but there is a major obstacle.

Write

\[
B_p(t)=\sum c_i t^i.
\]

Then

\[
b_d=[t^d]B_p(t)^2=\sum_i c_i c_{d-i}.
\]

A zero gives one quadratic relation.

Four nearby zeros give four quadratic equations. The hope is that the coefficient variables are too constrained.

The problem is that the number of unknown coefficients is roughly half the degree of \(B_p\), namely \(\asymp p\). Four equations do not create over-determination.

The square structure becomes useful only if combined with additional restrictions:

- low complexity of \(B_p\),
- bounded rank of the coefficient Hankel matrix,
- strong sparsity of the support of \(B_p\),
- extra Hasse-Witt constraints.

At present the square factorization alone does not beat the gap-polynomial bound.

---

# 4. Lyapunov exponent viewpoint

The archimedean Lyapunov exponent explains the exponential growth of solutions:

\[
\lambda=\frac12\log(17+12\sqrt2).
\]

It does not directly control modular returns.

A finite-field orbit is not governed by the archimedean norm. A large real Lyapunov exponent can coexist with many finite-field coincidences.

The relevant quantity would be a p-adic/Frobenius Lyapunov exponent for the associated overconvergent F-isocrystal.

Even then, a Lyapunov exponent controls average growth of Frobenius eigenvalues, not spacing of zeros of one coordinate of one vector.

Thus this route is suggestive but not presently a proof mechanism.

---

# 5. Quantitative value of eta

The data suggest that \(\eta\) should be a small absolute constant.

A random model gives a heuristic. If

\[
|Z_p|\approx O(1),
\]

then four zeros should have typical span comparable with p. The smallest four-zero span among many primes should decrease slowly by extreme-value statistics.

The observation that the first examples with four zeros have minimum span about \(p/5\) is consistent with a constant such as

\[
\eta=10^{-2}\text{ to }10^{-1}
\]

but this is only heuristic.

No theoretical value is known.

A realistic first target would be:

\[
\eta=10^{-3}.
\]

Any fixed positive eta would already imply the desired \(O(\sqrt p)\) bound.

---

# 6. What theorem would actually finish the route?

The missing theorem can be formulated cleanly:

For some absolute \(\eta>0\), if

\[
d_1<d_2<d_3<d_4,
\]

are Apéry zeros modulo p and

\[
d_4-d_1<\eta\sqrt p,
\]

then

\[
d_i+d_j=p-1
\]

for some pair.

This is a four-point incidence theorem for the marked Apéry Frobenius orbit.

It is stronger than:

- gap polynomial bounds,
- reflection symmetry,
- no consecutive zeros,
- Lucas congruences.

The likely ingredients are:

1. arithmetic of the integral period vector;
2. Frobenius/K3 monodromy;
3. a quantitative incidence theorem for the corresponding F-crystal.

---

# Final assessment

Among the proposed routes:

| Method | Expected power |
|---|---|
| Reflection symmetry alone | insufficient |
| Recurrence/no adjacent zeros | insufficient |
| Lyapunov exponent | heuristic only |
| Square factorization | promising but incomplete |
| Hasse-Witt/F-crystal | potentially relevant |
| Short-window reflection principle | sufficient if proved |

The short-window reflection principle is probably the cleanest formulation of the missing Apéry-specific theorem. It is not currently proved, but it identifies the exact local statement needed to replace the generic gap-polynomial argument and reach \(|Z_p|=O(\sqrt p)\).