ANSWER Q6424 c742fc04

# Adversarial referee report: the Franel–Mellin “bypass”

## Executive verdict

The finite-field algebra in links **(1)–(3) is correct** for primes \(p\ge 7\) and interior characters \(1\le r\le p-2\). The boundary bookkeeping can be closed completely, and no hidden correction occurs at branch points or at nonsplit fibers.

The claimed consequence **(4) does not follow from (1)–(3)**. There are three distinct assertions being conflated:

1. the two kernels are fixed **as functions of the Mellin character \(r\)** for a fixed prime \(p\);
2. those kernels are reductions of a compatible family of fixed-rank, uniformly bounded-conductor sheaves/crystals;
3. bounded conductor controls exact vanishing modulo the defining prime.

Only (1) follows from the displayed identities. Assertion (2) requires an independent geometric descent theorem, including an integral Frobenius-compatible realization. Assertion (3) is false in general, even for rank-one Kummer sheaves of fixed conductor.

More sharply, if one writes

\[
M_r=\sum_{t\in\mathbf F_p^\times}A_p(t)t^{-r},\qquad
C_r=\sum_{t\in\mathbf F_p^\times}\chi_2(q(t))A_p(t)t^{-r},
\]

and

\[
P_r=\sum_{\substack{x\in\mathbf F_p\\x\ne-1,\ \phi(x)\ne0}}
H_p(x)^2\phi(x)^{-r},
\]

then the entire derivation is

\[
P_r=M_r+C_r,\qquad M_r=-b_r,
\]

hence

\[
b_r=-P_r+C_r.
\]

Thus, **as a displayed two-sum identity, the “bypass” is the tautology \(P_r=M_r+C_r\)** coming from the fiber multiplicity \(1+\chi_2(q)\). Equivalently,

\[
C_r=b_r+P_r.
\]

The correction term contains exactly the information lost on the nonsplit fibers. Calling it an independent bounded-complexity companion is unjustified until a separate descent theorem constructs it from Franel geometry.

There is nevertheless real, non-tautological content: the rational pullback reconstructs \(A_p(t)\) from the Franel polynomial on the quadratic cover, including inert fibers over \(\mathbf F_{p^2}\). This removes coefficient-level circularity. It does **not** by itself produce a bounded-conductor trace function, and it gives no zero-counting theorem.

The GK/Jacobi convolution in link **(5)** is also an exact interior-character theorem once conventions are fixed. It does not rescue (4): it has \(p-1\) terms, and its Stickelberger unit block reduces term-by-term to the original Apéry binomial sum. It is a second exact coordinate system, not a complexity reduction for defining-characteristic zeros.

**Overall referee verdict:**

- **THEOREM:** links (1), (2), and (3), with the exact domains stated below.
- **THEOREM, but separate:** the interior GK convolution, after writing all character conventions and exceptional Jacobi cases.
- **BOOKKEEPING OWED:** endpoint corrections, branch-point tests in the scripts, and middle-extension conventions if the formula is promoted to sheaves.
- **UNPROVED IN THIS CHAIN:** “fixed bounded-conductor pair.”
- **WISHFUL:** obtaining zero density from complex/\(\ell\)-adic equidistribution.
- **REFORMULATION, not a quantitative reduction:** “the sole remaining target is defining-characteristic Mellin zero density.”

---

# 1. Notation and exact range

Let

\[
N=p-1,
\qquad
A_p(t)=\sum_{n=0}^{N}b_nt^n,
\qquad
H_p(x)=\sum_{n=0}^{N}f_nx^n,
\]

where \(f_n=\sum_k\binom nk^3\) is the Franel sequence, and put

\[
\phi(x)=\frac{x(1-8x)}{1+x},
\qquad
q(t)=t^2-34t+1.
\]

Throughout the finite-field derivation, \(p\ge7\), so \(2,3,8,9\) are units. Extend the quadratic character by

\[
\chi_2(0)=0.
\]

The Mellin identity is naturally indexed by characters of \(\mathbf F_p^\times\), hence by \(r\bmod N\). The clean coefficient-extraction range is

\[
1\le r\le N-1=p-2.
\]

At \(r\equiv0\pmod N\), coefficient \(0\) and coefficient \(N\) alias. Any statement claiming the individual value \(b_0\) or \(b_N\) from that Mellin character requires an endpoint correction.

---

# 2. Link (1): derivation of the rational and pointwise pullback

This link can be proved formally; it need not remain a numerical observation.

Let

\[
\mathcal F(t)=\sum_{n\ge0}b_nt^n,
\qquad
h(x)=\sum_{n\ge0}f_nx^n.
\]

Use the characteristic-zero Apéry–Franel transformation

\[
\mathcal F(\phi(x))=(1+x)h(x)^2.
\tag{2.1}
\]

The Dwork/Lucas congruences are

\[
\mathcal F(t)\equiv A_p(t)\mathcal F(t^p)\pmod p,
\tag{2.2}
\]

and

\[
h(x)\equiv H_p(x)h(x^p)\pmod p.
\tag{2.3}
\]

Substitute \(t=\phi(x)\) in (2.2). Since \(\phi\) has coefficients in \(\mathbf F_p\),

\[
\phi(x)^p=\phi(x^p).
\]

Using (2.1) at \(x\) and \(x^p\), and (2.3), gives

\[
(1+x)H_p(x)^2h(x^p)^2
=A_p(\phi(x))(1+x^p)h(x^p)^2.
\]

The series \(h(x^p)\) has constant term \(1\), so it is a unit and cancels. Since

\[
1+x^p=(1+x)^p
\]

in characteristic \(p\), one obtains

\[
\boxed{
A_p(\phi(x))(1+x)^{p-1}=H_p(x)^2.
}
\tag{2.4}
\]

Although derived in \(\mathbf F_p[[x]]\), (2.4) is a polynomial identity after clearing the denominator:

\[
(1+x)^N A_p\!\left(\frac{x(1-8x)}{1+x}\right)
=\sum_{n=0}^{N}b_nx^n(1-8x)^n(1+x)^{N-n}.
\]

Therefore (2.4) holds in \(\mathbf F_p[x]\), not merely to a finite jet.

For every \(x\in\mathbf F_p\setminus\{-1\}\), Fermat gives

\[
(1+x)^{p-1}=1,
\]

hence

\[
\boxed{A_p(\phi(x))=H_p(x)^2.}
\tag{2.5}
\]

### Referee verdict on link (1)

**Holds**, assuming the three named inputs (2.1)–(2.3) are already proved with exactly these normalizations. The pointwise script is only a regression test; the proof is the four-line cancellation above.

The point \(x=-1\) must remain excluded. It is a pole of \(\phi\), and Fermat cannot be invoked there.

---

# 3. Link (2): the fiber discriminant and every exceptional fiber

The equation \(\phi(x)=t\) is

\[
t(1+x)=x-8x^2,
\]

or

\[
8x^2+(t-1)x+t=0.
\tag{3.1}
\]

Its discriminant is exactly

\[
(t-1)^2-32t=t^2-34t+1=q(t).
\tag{3.2}
\]

Because the leading coefficient \(8\) is a unit, the number of affine roots is

\[
\boxed{
\nu(t)=1+\chi_2(q(t)).
}
\tag{3.3}
\]

This formula includes all three cases:

- \(q(t)\) a nonzero square: two roots;
- \(q(t)=0\): one double root;
- \(q(t)\) a nonsquare: no roots.

## 3.1 The excluded pole \(x=-1\)

Substituting \(x=-1\) in the cleared fiber polynomial (3.1) gives

\[
8-(t-1)+t=9.
\]

For \(p\ge7\), this is nonzero. Therefore \(x=-1\) is not secretly a root of any finite fiber. Removing it from the domain loses no finite \(t\)-point.

The rational map on \(\mathbf P^1\) has two points above \(t=\infty\), namely \(x=-1\) and \(x=\infty\); neither belongs to the affine Mellin sum over \(t\in\mathbf F_p^\times\).

## 3.2 The zero fiber

At \(t=0\), equation (3.1) becomes

\[
x(8x-1)=0.
\]

Thus the two roots are

\[
x=0,\qquad x=\frac18.
\]

Also \(q(0)=1\), so (3.3) correctly gives two roots. Since the Mellin character \(t^{-r}\) is defined only for \(t\ne0\), the \(t=0\) term is omitted, and both \(x\)-points with \(\phi(x)=0\) must be omitted. This is exactly the exclusion in the claimed formula. No correction remains.

## 3.3 Branch points

If \(q(t)=0\), then \(\chi_2(q(t))=0\) and \(\nu(t)=1\). The pullback sum contains one copy of \(A_p(t)t^{-r}\), while the companion sum contains zero copies. The final coefficient is therefore \(-1\), exactly as required.

The current script `franel_mellin_mult_test.py` does **not** test this case: its fiber assertion is guarded by

```python
if (t*t - 34*t + 1) % p
```

and therefore skips every branch point. The code should be repaired, even though the omitted case is settled by (3.1)–(3.3).

## 3.4 Nonsplit fibers

If \(q(t)\) is a nonsquare, then \(\nu(t)=0\). There is no \(x\in\mathbf F_p\) above \(t\), and the first sum contributes nothing. The correction contributes

\[
- A_p(t)t^{-r}.
\]

This is exactly the entire desired Mellin contribution at that \(t\). Consequently the nonsplit locus is not a harmless error set: it is where the correction term carries all of the missing data.

### Referee verdict on link (2)

**Holds.** The numerical verification is incomplete at branch points, but the quadratic calculation is decisive.

---

# 4. Link (3): exact Mellin bookkeeping and the endpoint alias

For every integer \(m\), multiplicative orthogonality in \(\mathbf F_p\) gives

\[
\sum_{t\in\mathbf F_p^\times}t^m
=
\begin{cases}
-1,&N\mid m,\\
0,&N\nmid m.
\end{cases}
\tag{4.1}
\]

Define

\[
M_r:=\sum_{t\in\mathbf F_p^\times}A_p(t)t^{-r}.
\]

Expanding \(A_p\),

\[
M_r
=\sum_{n=0}^{N}b_n
\sum_{t\in\mathbf F_p^\times}t^{n-r}.
\tag{4.2}
\]

If \(1\le r\le N-1\), then

\[
-(N-1)\le n-r\le N-1.
\]

The only multiple of \(N\) in this interval is \(0\), so the only contributing index is \(n=r\). Therefore

\[
\boxed{M_r=-b_r\qquad(1\le r\le p-2).}
\tag{4.3}
\]

There is no mysterious cancellation between \(b_0\) and \(b_{p-1}\). In the interior range, each endpoint term vanishes separately:

\[
N\nmid -r,
\qquad
N\nmid N-r.
\]

At \(r\equiv0\pmod N\), both endpoints contribute:

\[
M_0=-(b_0+b_N)=-2,
\tag{4.4}
\]

because \(b_0=b_N=1\). Thus a character-space Mellin transform cannot distinguish \(b_0\) from \(b_N\). The formula below returns \(2\), not either endpoint value \(1\). This is the exact endpoint correction.

Now put

\[
P_r:=\sum_{\substack{x\in\mathbf F_p\\x\ne-1,\ \phi(x)\ne0}}
H_p(x)^2\phi(x)^{-r},
\]

and

\[
C_r:=\sum_{t\in\mathbf F_p^\times}
\chi_2(q(t))A_p(t)t^{-r}.
\]

Using the pointwise identity and grouping by fibers,

\[
\begin{aligned}
P_r
&=\sum_{t\in\mathbf F_p^\times}
\nu(t)A_p(t)t^{-r}\\
&=\sum_{t\in\mathbf F_p^\times}
(1+\chi_2(q(t)))A_p(t)t^{-r}\\
&=M_r+C_r.
\end{aligned}
\tag{4.5}
\]

Combining (4.3) and (4.5),

\[
\boxed{
 b_r
 =-\sum_{\substack{x\in\mathbf F_p\\x\ne-1,\ \phi(x)\ne0}}
 H_p(x)^2\phi(x)^{-r}
 +\sum_{t\in\mathbf F_p^\times}
 \chi_2(q(t))A_p(t)t^{-r}
}
\tag{4.6}
\]

for every \(1\le r\le p-2\).

A per-fiber audit makes the absence of corrections completely transparent:

| \(\chi_2(q(t))\) | \(\nu(t)\) | coefficient in \(-P_r+C_r\) |
|---:|---:|---:|
| \(1\) | \(2\) | \(-2+1=-1\) |
| \(0\) | \(1\) | \(-1+0=-1\) |
| \(-1\) | \(0\) | \(0-1=-1\) |

Thus every nonzero \(t\) contributes exactly

\[
-A_p(t)t^{-r}.
\]

At \(r\equiv0\pmod N\), the same calculation gives

\[
-P_0+C_0=-M_0=2,
\]

so one must subtract \(1\) to recover either conventionally selected endpoint coefficient.

### Referee verdict on link (3)

**Holds exactly in the stated interior range.** The phrase “the wraparound terms cancel” should be deleted. They do not occur in the interior and they add at the endpoint.

---

# 5. The precise circularity in link (4)

Equation (4.6) looks like a decomposition of \(b_r\) into two new Mellin values. Algebraically it contains only

\[
P_r=M_r+C_r,
\qquad
M_r=-b_r.
\]

In particular,

\[
\boxed{C_r=b_r+P_r.}
\tag{5.1}
\]

This is the exact circularity. The correction term is not an independently controlled error. It is the original coefficient plus the Franel pullback term.

The reason is visible before summing. Since

\[
\nu(t)=1+\chi_2(q(t)),
\]

the coefficient identity behind (4.6) is simply

\[
-\nu(t)+\chi_2(q(t))=-1.
\tag{5.2}
\]

Therefore the two-sum formula, by itself, is a change-of-variables tautology. It cannot yield a new estimate unless one proves that \(C_r\) has an independent realization with useful structure.

The split locus does not provide that realization. On split fibers, \(A_p(t)\) is indeed visible as \(H_p(x)^2\) at either rational preimage. On nonsplit fibers there is no \(\mathbf F_p\)-point upstairs, and the entire value \(-A_p(t)\) is restored by the companion. The “error” is exactly the inaccessible inert information.

### Referee verdict

The sentence

> “The term explosion obstruction is bypassed because \(b_r\) is the Mellin value of a fixed bounded-conductor pair”

is **not proved by links (1)–(3)**. At this stage one has a fixed-in-\(r\) pair of finite-field functions, one member of which is explicitly built from the original Apéry Hasse polynomial.

---

# 6. What genuinely non-circular content survives

The situation is not completely tautological. The rational identity reconstructs the inert values from Franel data over the quadratic cover.

## 6.1 The deck involution

If \(x\) is one root of the fiber equation, the other root is

\[
\iota(x)=\frac{1-8x}{8(1+x)}.
\tag{6.1}
\]

Indeed, the product of the two roots of (3.1) is \(t/8\). One checks

\[
\phi\circ\iota=\phi,
\qquad
\iota^2=1,
\tag{6.2}
\]

and

\[
q(\phi(x))
=\left(\frac{1-16x-8x^2}{1+x}\right)^2.
\tag{6.3}
\]

Thus the \(x\)-line is the fixed quadratic cover of the \(t\)-line with branch divisor \(q(t)=0\).

## 6.2 Scalar descent of the Franel Hasse section

Let \(e=N/2\) and define the rational function

\[
G_p(x)=\frac{H_p(x)}{(1+x)^e}.
\tag{6.4}
\]

Equation (2.4) says

\[
G_p(x)^2=A_p(\phi(x)).
\tag{6.5}
\]

Since \(\phi\circ\iota=\phi\),

\[
G_p(\iota(x))^2=G_p(x)^2.
\]

In the field \(\mathbf F_p(x)\), characteristic is not \(2\), so

\[
\boxed{G_p\circ\iota=\epsilon_pG_p}
\qquad
(\epsilon_p\in\{1,-1\}).
\tag{6.6}
\]

For every finite \(t\), choose either geometric preimage \(x\in\mathbf F_{p^2}\). Then

\[
\boxed{A_p(t)=G_p(x)^2.}
\tag{6.7}
\]

This is independent of the chosen root by (6.6). If the fiber is inert, then \(x^p=\iota(x)\), so

\[
G_p(x)^p=\epsilon_pG_p(x),
\]

and \(G_p(x)^2\in\mathbf F_p\), as required.

Therefore the companion can be rewritten without using the coefficients \(b_n\):

\[
C_r
=\sum_{t\in\mathbf F_p^\times}
\chi_2(q(t))G_p(x_t)^2t^{-r},
\tag{6.8}
\]

where \(x_t\) is any geometric point above \(t\).

This is genuine non-circular finite-field content. It says that both terms arise from the same Franel Hasse section on a fixed quadratic cover, with invariant and anti-invariant descent.

The CFVZ factorization

\[
A_p=S_p^2
\quad\text{or}\quad
A_p=qS_p^2
\]

and the identification of \(S_p\) with a truncation of a fixed \(\tau\)- or \(\sigma\)-series give another coefficient-level non-circular description. They do not change the next obstruction: the degrees of these Hasse sections grow linearly with \(p\).

## 6.3 What (6.7) does not prove

Equation (6.7) is a statement about scalar Hasse functions. It does not automatically construct:

1. an \(\ell\)-adic or crystalline compatible system whose integral Frobenius trace reduces to \(G_p(x)^2\);
2. a Frobenius-compatible equivariant structure under \(\iota\);
3. the invariant and anti-invariant descents on the \(t\)-line;
4. bounded conductor at \(0\), the roots of \(q\), and infinity;
5. the exact middle-extension trace values at the branch points.

Those are separate geometric theorems. The later note `GPT_Q6394_descent_theorem.md` attempts precisely this upgrade and explicitly retains `[GAP-1: explicit isogeny formula]`. Until that correspondence or another integral descent is written, “bounded-conductor pair” remains conditional.

### Referee verdict on non-circularity

- **The two-sum equation alone:** tautological.
- **The Franel reconstruction on the quadratic cover:** real and non-circular.
- **Bounded-conductor compatible-system descent:** extra theorem, not supplied here.
- **Any zero-density consequence:** still absent.

---

# 7. “Fixed function” is not “fixed bounded-conductor sheaf”

For a fixed prime \(p\), both

\[
x\longmapsto H_p(x)^2
\]

and

\[
t\longmapsto \chi_2(q(t))A_p(t)
\]

are independent of \(r\). In that weak sense they are fixed Mellin kernels.

That is not the sense required by Deligne, Katz, Fouvry–Kowalski–Michel, or Forey–Fresán–Kowalski. The required statement is a compatible family of middle-extension sheaves or overconvergent crystals with ranks, Swan conductors, and singular loci uniformly bounded as \(p\) varies, together with an integral trace congruence.

The polynomial data do not imply this:

- \(A_p\) has degree \(p-1\) and contains every Apéry residue as a coefficient;
- \(H_p\) has degree \(p-1\);
- \(S_p\) has degree about \(p/2\);
- an arbitrary function on \(\mathbf F_p\) has an interpolation polynomial of degree at most \(p-1\).

Thus polynomial presentation alone gives conductor on the scale of \(p\), not \(O(1)\).

It may be true that \(A_p\) and \(H_p\) are Hasse–Witt reductions of fixed geometric families. But a Hasse invariant is a Cartier/Frobenius scalar in characteristic \(p\), not automatically a complex trace function. A proof must supply the full integral cohomological object and the reduction of its Frobenius trace. This is the category mismatch already recorded in the project as `[GAP-CARTIER]` or `[GAP-1]`.

At the two branch points, the finite-field formula uses \(\chi_2(0)=0\). A sheaf formulation must use the correct middle extension of the quadratic Kummer sheaf, whose branch-point trace is zero. This is straightforward once the sheaf is constructed, but it is not optional bookkeeping.

---

# 8. Link (5): the GK/Jacobi convolution

The formula

\[
b_r=-\frac1N\sum_{t=0}^{N-1}
J(\varphi_2\psi^{t+r},\psi^{N/2-t})^2
J(\varphi_2\psi^{t-r},\psi^{N/2-t})^2
\pmod{\mathfrak p_\psi}
\tag{8.1}
\]

is consistent with the Greene hypergeometric formula and with the direct mod-\(p\) checks. Here \(\varphi_2\) denotes the quadratic character, to avoid confusing it with the rational map \(\phi(x)\).

The repository script `q6356_gk_convolution_check.py` checks

\[
1\le r\le p-2
\]

for \(p=29,37,41\). It does not check \(r=0\) or \(r=N\). Its Jacobi convention omits \(0\) and \(1\), which is equivalent to extending every multiplicative character—including the trivial character—by zero at \(0\). Any written theorem must state this convention, because the exceptional identities

\[
J(\varepsilon,\chi),\qquad
J(\chi,\chi^{-1}),\qquad
J(\varepsilon,\varepsilon)
\]

depend on it.

For interior \(r\), (8.1) is a legitimate theorem after one writes the Greene-to-Jacobi normalization and all inverse-character degeneracies. At the endpoints \(A_r=\varepsilon\), the generic Greene conversion degenerates and a separate correction is required.

## 8.1 Why the convolution does not lower complexity

The sum in (8.1) has \(N=p-1\) terms. The exact Stickelberger ledger, after reindexing, gives the valuation of the \(k\)-th product as

\[
2\mathbf 1_{k>r}+2\mathbf 1_{k>N-r}.
\tag{8.2}
\]

Hence the slope-zero block contains

\[
1+\min(r,N-r)
\]

unit terms. Their reductions are precisely the Apéry summands

\[
\binom rk^2\binom{r+k}{k}^2,
\]

and their sum is \(b_r\), using reflection when \(r>N/2\).

Therefore the GK formula exposes, rather than removes, the defining-characteristic cancellation:

\[
p\mid b_r
\]

means cancellation among a block of equally minimal \(p\)-adic valuation. There is no unique minimal term and no bounded-size Stickelberger container.

### Referee verdict on link (5)

**Exact and useful as a formula; irrelevant as evidence that the Franel–Mellin correction is bounded conductor.** It is an independent coordinate description of the same residue and retains the \(O(p)\) term explosion.

---

# 9. The defining-characteristic caveat is fatal to the proposed counting consequence

Assume the best possible geometric upgrade: suppose \(b_r\) is the reduction modulo \(\mathfrak p\mid p\) of the Mellin Frobenius trace of one fixed bounded-conductor sheaf.

Complex or \(\ell\)-adic equidistribution controls quantities such as

\[
\frac{M_p(\chi)}{p^{w/2}}
\]

under complex embeddings, and fixed continuous class functions of the associated compact monodromy conjugacy class. It does not control

\[
M_p(\chi)\equiv0\pmod{\mathfrak p}.
\]

The latter is a \(p\)-adic valuation event depending on the integral lattice and on the same prime that defines the field.

## 9.1 A decisive bounded-conductor counterexample

Fix \(d\ge2\). For every prime \(p\equiv1\pmod d\), let

\[
\psi=\omega_p^{-(p-1)/d}
\]

be a character of order \(d\). The rank-one Kummer sheaf attached to \(\psi\) has fixed rank and conductor. Its Mellin transforms are

\[
M_j=J(\omega_p^{-j},\psi).
\]

Gross–Koblitz/Stickelberger gives

\[
v_p(M_j)=1
\]

for every

\[
j\ge \frac{(p-1)(d-1)}d,
\]

up to the standard endpoint exclusions. Thus

\[
\#\{j:p\mid M_j\}\ge \frac{p-1}{d}-O(1).
\]

So a fixed rank-one, fixed-conductor sheaf with square-root Weil bounds can have a positive proportion of its tame Mellin transforms vanish modulo the defining prime.

This refutes any implication of the form

\[
\text{bounded conductor + big monodromy/equidistribution}
\Longrightarrow o(p)\text{ defining-characteristic zeros}.
\]

Apéry-specific arithmetic may still force zero density, but it must be proved separately.

## 9.2 The character parameter is not a fixed curve

The character set is the finite étale scheme

\[
\widehat{\mathbf F_p^\times}\simeq\mu_{p-1}.
\]

Its coordinate ring is

\[
\mathbf F_p[X]/(X^{p-1}-1),
\]

which is the full algebra of functions on \(p-1\) points. Every function on the character set has a polynomial representative of degree at most \(p-2\). Therefore saying that the zero set is the divisor of a Hasse section on character space gives no bound unless one independently proves that the section has degree \(O(1)\) or \(o(p)\).

The raw Hasse section here is exactly \(r\mapsto b_r\). Gross–Koblitz does not lower its degree; it expands it into a growing block of unit terms.

---

# 10. Reduction versus reformulation

The statement

\[
p\mid b_r
\iff
\text{a Mellin value of the Franel/Apéry companion vanishes in }\mathbf F_p
\]

is exact. Quantitatively, it is the same zero set.

It becomes a genuine reduction only if the new side satisfies an independently proved theorem not already equivalent to \(b_r=0\), for example:

- a bounded-degree Hasse polynomial on a genuinely fixed parameter curve;
- a theorem that all but \(O(p^{1-\delta})\) twists have a unique unit root;
- a defining-characteristic local-limit theorem for the specific Mellin compatible system;
- an explicit arithmetic formula reducing the zero event to a bounded-order character condition.

None follows from (1)–(5).

The project’s character-order counting lemma does provide a genuine separate reduction: low-order characters are sparse in the moving-prime window. But replacing the remaining high-order event by the phrase “Mellin non-ordinarity” does not simplify it. It names the original condition in crystalline language.

There is also a benchmark mismatch worth correcting. For the particular zero fiber, the project already records

\[
\#\{r:b_r=0\}\le 3p^{2/3}+O(1).
\]

The bound \(8p^{3/4}\) concerns the uniform all-value multiplicity. Any claimed improvement for the zero set should be compared with the \(2/3\) exponent, not \(3/4\).

---

# 11. Required repairs before this chain is used downstream

1. **Promote link (1) to a written proof.** Use equations (2.1)–(2.4), not only coefficient tests.

2. **Fix the fiber-count regression test.** Test the branch points \(q(t)=0\); the current script skips them.

3. **State the exact Mellin range.** Formula (4.6) is for \(1\le r\le p-2\). At the trivial character it returns \(b_0+b_{p-1}=2\).

4. **Delete the phrase “wraparound cancellation.”** The endpoint coefficients are absent separately in the interior and add at the endpoint.

5. **Separate three meanings of “fixed.”** Write “fixed in \(r\) for each \(p\)” unless a compatible-system theorem has actually been proved.

6. **Label the descent theorem as an extra dependency.** The scalar cover identity is proved. Frobenius-compatible invariant/anti-invariant descent is not a consequence of it. The present project note still marks the explicit isogeny as `[GAP-1]`.

7. **Write the integral trace congruence.** An \(\ell\)-adic complex trace and an \(\mathbf F_p\)-valued Hasse invariant are different objects. Specify the lattice, prime \(\mathfrak p\mid p\), Tate factors, and branch-point middle extensions.

8. **Do not cite FFK/Katz as a zero-density theorem.** Their equidistribution conclusions are archimedean/\(\ell\)-adic. A new defining-characteristic theorem is required.

9. **Keep the GK endpoint conventions separate.** The interior convolution and the trivial-character endpoint are different statements.

10. **Downgrade the strategic conclusion.** The correct current statement is:

> The Franel pullback packages the Apéry coefficient as the defining-characteristic reduction of Mellin data on a fixed quadratic cover. Conditional on a Frobenius-compatible bounded-conductor descent, it removes the length-\(p\) expansion at the sheaf-object level. It does not by itself imply any saving for the number of mod-\(p\) zeros.

---

# 12. Final classification of the five links

## Link (1): rational pullback and pointwise equality

**THEOREM.** Complete proof in Section 2. Exclude \(x=-1\).

## Link (2): discriminant and fiber multiplicity

**THEOREM.** Complete proof in Section 3. Branch fibers are handled by \(\chi_2(0)=0\). The numerical script omitted them.

## Link (3): Mellin-with-multiplicity identity

**THEOREM for \(1\le r\le p-2\).** No corrections. At the trivial character the expression equals \(2\), not \(1\).

## Link (4): fixed bounded-conductor pair and bypass

- “\(r\) appears only through a Mellin character”: **TRUE**.
- “the kernels are fixed finite-field functions for each \(p\)”: **TRUE**.
- “the displayed second kernel is independent of the original Apéry data”: **FALSE as written**; it contains \(A_p\), and equation (5.1) displays the circularity.
- “all values can be reconstructed from Franel data on the quadratic cover”: **TRUE at the scalar Hasse-function level**.
- “there is a fixed bounded-conductor compatible-system pair”: **[GAP-DESCENT] in this chain**.
- “the term explosion is bypassed for purposes of mod-\(p\) zero counting”: **FALSE / WISHFUL**.

## Link (5): GK convolution

**THEOREM for nontrivial characters after conventions are fixed.** Endpoint corrections remain separate. It has \(O(p)\) terms and its unit block is the original Apéry sum, so it does not prove the bypass.

---

# Bottom line

The finite-field identities are solid and worth keeping. The strongest honest interpretation is that they reveal a quadratic-cover descent architecture for the Franel and Apéry Hasse functions.

The sentence “the BGK obstruction is bypassed and only standard Mellin zero density remains” is not acceptable. The displayed two-term identity is tautological until the inert companion is constructed independently; even after a bounded-conductor descent is proved, defining-characteristic zero density is a new theorem of essentially the original difficulty.

So the bypass is **not wholly empty**, but it has not yet paid any quantitative dividend. It changes the geometric packaging, not the arithmetic hardness.

## Least-confident step

The least-confident judgment is the current status of the Frobenius-compatible descent, not any of the finite-field algebra. The later project note `GPT_Q6394_descent_theorem.md` gives a convincing connection-level and monodromy-level construction but explicitly leaves the normalized isogeny correspondence as `[GAP-1]`. It is possible that the required correspondence is already implicit in a standard modular model and can be written down quickly. If that gap is closed, the bounded-conductor **object** claim should be upgraded from conditional to theorem; none of the negative conclusions about defining-characteristic zero counting would change.