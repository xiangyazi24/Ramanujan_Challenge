ANSWER Q3579 dde8bef0

# P3.2 — can geometric monodromy improve the Apéry zero-set bound?

## Verdict

For the project sequence

\[
b_r=A_r=\sum_{k=0}^r \binom rk^2\binom{r+k}{k}^2
      =1,5,73,1445,33001,\ldots,
\]

I do **not** find a published theorem that upgrades

\[
Z_p=\{0\le r<p:A_r\equiv0\pmod p\}
\]

to

\[
|Z_p|\ll_\varepsilon p^{1/2+\varepsilon}
\]

from the geometric origin of the Apéry numbers. In fact, the proposed `Sp4 monodromy -> square-root zero count` route has two separate category errors.

1. For the classical \(\zeta(3)\) Apéry numbers used in this repository, the **minimal Picard–Fuchs system is rank 3**, coming from a one-parameter family of K3 surfaces, and the order-3 Apéry differential equation is a **symmetric square of an order-2 equation**. The relevant monodromy is therefore orthogonal / symmetric-square type, not a rank-4 \(Sp_4\) local system.

2. Much more importantly, \(Z_p\) is a set of **zero coefficients** of a polynomial truncation. It is not the set of \(\mathbf F_p\)-points where a bounded-degree geometric trace function vanishes. Weil/Deligne point-counting and large-monodromy square-root cancellation do not directly control coefficient sparsity.

As a literature audit, I also do **not** find the exact bound \(|Z_p|\ll p^{2/3}\) stated as a theorem in the standard Apéry divisibility literature. Malik–Straub prove Lucas structure and a reflection symmetry and then discuss zero/divisibility statistics heuristically; they do not prove a sublinear zero-count estimate of this form. Thus the repository's \(p^{2/3}\) estimate should be regarded as an **internal lemma whose proof must be audited separately**. Conditional on that lemma being valid, I find no published geometric improvement below exponent \(2/3\), in particular no \(1/2+\varepsilon\) result.

The strongest recent geometric/algebraic information I found actually explains why the naive bounded-genus argument cannot work: Caruso–Fürnsinn–Vargas-Montoya–Zudilin (2026) prove that the mod-\(p\) Apéry generating function generates a cyclic extension of \(\mathbf F_p(t)\) of degree \((p-1)/2\) or \(p-1\). Its algebraic complexity therefore **grows with \(p\)** rather than remaining bounded.

---

## 1. The classical Apéry object here is K3 / rank 3, not rank 4 \(Sp_4\)

Let

\[
F(t)=\sum_{n\ge0}A_n t^n.
\]

The recurrence in the repository is equivalent to the order-3 differential equation

\[
t^2(1-34t+t^2)F'''
+t(3-153t+6t^2)F''
+(1-112t+7t^2)F'
+(t-5)F=0.
\]

Beukers and Peters identified this equation with the Picard–Fuchs equation of a one-parameter family of K3 surfaces; the generic transcendental part has rank 3. See:

- F. Beukers and C. A. M. Peters, **A family of K3 surfaces and \(\zeta(3)\)**, *J. reine angew. Math.* **351** (1984), 42–54, DOI `10.1515/crll.1984.351.42`.
- J. Stienstra and F. Beukers, **On the Picard-Fuchs equation and the formal Brauer group of certain elliptic K3-surfaces**, *Math. Ann.* **271** (1985), 269–304, DOI `10.1007/BF01455990`.

Moreover, the Apéry order-3 operator is the symmetric square of the order-2 operator

\[
t(t^2-34t+1)y''+(2t^2-51t+1)y'+\frac14(t-10)y=0.
\]

This symmetric-square fact goes back to Apéry/Dwork and is explicitly discussed in Beukers' expositions and in later work on the Apéry differential equation. It is also exactly the structure exploited in the 2026 Caruso–Fürnsinn–Vargas-Montoya–Zudilin paper.

So the statement “the Apéry sequence has \(Sp_4\) geometric monodromy because it is a diagonal in four variables” is not correct for the classical sequence \(1,5,73,\ldots\). A four-*variable* rational diagonal is not the same thing as a rank-4 Picard–Fuchs local system.

Armin Straub proves the four-variable diagonal representation

\[
A_n=[x_1^n x_2^n x_3^n x_4^n]
\frac{1}{(1-x_1-x_2)(1-x_3-x_4)-x_1x_2x_3x_4}.
\]

See A. Straub, **Multivariate Apéry numbers and supercongruences of rational functions**, *Algebra & Number Theory* **8** (2014), 1985–2008, DOI `10.2140/ant.2014.8.1985`, arXiv:`1401.0854`.

That is a coefficient-extraction realization. It does not change the minimal order-3 Picard–Fuchs equation governing \(F(t)\).

Even if one embeds the motive into some larger ambient cohomology carrying a symplectic representation, large monodromy of that ambient representation would still not solve the actual zero-coefficient problem below.

---

## 2. The main obstruction: \(Z_p\) is a coefficient-zero set, not a geometric root locus

Define the usual truncation

\[
A_p(t)=\sum_{r=0}^{p-1} A_r t^r\in\mathbf F_p[t].
\]

Then

\[
Z_p=\{r: [t^r]A_p(t)=0\}.
\]

Compare this with the genuinely geometric set

\[
R_p=\{x\in\mathbf F_p:A_p(x)=0\}.
\]

These are completely different statistics.

A Weil/Hasse–Weil/Deligne theorem can control \(R_p\), or more generally the number of field parameters \(x\) satisfying a geometric condition, provided one has a bounded-complexity variety/sheaf over the **parameter \(x\)**. It does not count how many Taylor coefficients of its period polynomial happen to vanish.

In particular, the sentence

> “if the zero locus comes from a curve or surface of bounded genus/degree, apply Weil”

would help only after proving a representation of the **index map**

\[
r\longmapsto A_r\pmod p
\]

as the value/trace of a bounded-conductor family with \(r\in\mathbf F_p\) as geometric parameter. The K3 family uses the generating-function variable \(t\), not the coefficient index \(r\). No such bounded-conductor `index sheaf` for \(r\mapsto A_r\) is supplied by the diagonal/K3 construction.

This distinction is not cosmetic. The literature on truncated hypergeometric polynomials shows that geometry can indeed give strong root bounds in the **variable** of the truncation, but those results address \(\#\{x:Q_p(x)=s\}\), not the number of zero coefficients of \(Q_p\). For example:

- A. Ghosh and K. Ward, **Counting roots of truncated hypergeometric series over finite fields**, arXiv:`1601.06765`, prove a general \(O(p^{11/12})\) bound for roots of certain truncated \({}_2F_1\) polynomials and obtain sharper bounds in special elliptic/K3 geometric cases.

That is exactly the right comparison: strong geometry can improve a **root-counting** problem, but \(Z_p\) is not that problem.

---

## 3. What the best current mod-\(p\) structure actually says

The key modern theorem is:

X. Caruso, F. Fürnsinn, D. Vargas-Montoya, W. Zudilin,
**Galois Groups of Apéry-like Series Modulo Primes**, *Bull. Aust. Math. Soc.* **114** (2026), 65–78, DOI `10.1017/S0004972725100932`, arXiv:`2510.23298`.

Let

\[
f_\alpha(t)=\sum_{n\ge0}A_n t^n.
\]

Gessel's Lucas congruence implies

\[
f_\alpha(t)\equiv A_p(t)f_\alpha(t)^p\pmod p.
\]

Caruso et al. prove the following two very concrete facts.

### Theorem 1.1 of Caruso et al.

Writing \(S=(\mathbf F_p^\times)^2\),

\[
\operatorname{Gal}(\mathbf F_p(t,f_\alpha)/\mathbf F_p(t))
=\begin{cases}
S,&p\equiv1,5,7,11\pmod{24},\\
\mathbf F_p^\times,&p\equiv13,17,19,23\pmod{24}.
\end{cases}
\]

Thus the algebraic degree is respectively

\[
\frac{p-1}{2}\quad\text{or}\quad p-1.
\]

### Theorem 1.2 of Caruso et al.

There is \(B_p(t)\in\mathbf F_p[t]\) such that

\[
A_p(t)=\begin{cases}
B_p(t)^2,&p\equiv1,5,7,11\pmod{24},\\
(t^2-34t+1)B_p(t)^2,&p\equiv13,17,19,23\pmod{24}.
\end{cases}
\]

This is a very strong theorem, but it acts in the **polynomial variable \(t\)**. It determines root multiplicities/factorization of \(A_p(t)\); it does not say that the coefficient vector of \(A_p\) is sparse or dense.

Most importantly for the proposed Weil argument, Theorem 1.1 says that the natural algebraic cover attached to the mod-\(p\) generating series has degree \(\asymp p\), not bounded degree. So there is no fixed-genus/fixed-conductor family here to which one can blindly apply a uniform \(O(p^{1/2})\) error term and read off \(|Z_p|\).

Indeed, the factorization \(A_p=B_p^2\) is not a coefficient-support theorem: coefficients of a square are convolutions

\[
[t^r]B_p(t)^2=\sum_{i+j=r}b_i b_j,
\]

and can vanish through cancellation for reasons unrelated to the roots of \(B_p\).

---

## 4. What is actually known directly about divisibility among the first \(p\) Apéry numbers?

The standard direct reference is:

A. Malik and A. Straub, **Divisibility properties of sporadic Apéry-like numbers**, *Research in Number Theory* **2** (2016), article 5, DOI `10.1007/s40993-016-0036-8`, arXiv:`1508.00297`.

For the classical Apéry numbers they use Gessel's Lucas congruences and prove, in particular, the reflection congruence (their Lemma 6.2)

\[
A_r\equiv A_{p-1-r}\pmod p,\qquad 0\le r<p.
\]

Thus \(Z_p\) is symmetric under \(r\leftrightarrow p-1-r\). Together with Lucas congruences, this lets them reduce the question whether a prime ever divides an Apéry number to the first base-\(p\) digit block.

But their discussion of how often zeros occur is explicitly probabilistic/heuristic. They model the values in roughly half the interval as independent uniform residues and from this predict a positive density of primes with no Apéry zero at all. That is evidence that the true \(|Z_p|\) should usually be tiny (random-residue heuristics suggest order 1), but it is not a deterministic upper bound.

For background on the Lucas congruence itself, see I. Gessel, **Some congruences for Apéry numbers**, *J. Number Theory* **14** (1982), 362–368, DOI `10.1016/0022-314X(82)90071-3`.

The older K3/modular work gives remarkable congruences at special indices — for example Beukers' congruences and the modular-form relation at \((p-1)/2\) — but again does not bound the number of all \(r<p\) for which \(A_r\equiv0\pmod p\).

---

## 5. Why large monodromy alone would not give a zero-count even if \(Sp_4\) were present

Suppose, hypothetically, one had a bounded-conductor lisse sheaf \(\mathcal F\) on \(\mathbf A^1/\mathbf F_p\) with a trace function

\[
T_p(x)=\operatorname{tr}(\operatorname{Frob}_x\mid\mathcal F_x).
\]

Large geometric monodromy plus Deligne's Riemann hypothesis can give square-root cancellation in **linear sums** such as

\[
\sum_{x\in\mathbf F_p}T_p(x)\psi(g(x)).
\]

It does not automatically imply

\[
\#\{x:T_p(x)=0\}\ll p^{1/2+\varepsilon}.
\]

The indicator of an exact trace value is nonlinear. To use additive-character orthogonality one would need uniform square-root control of

\[
\sum_x \psi(hT_p(x))
\]

for all \(h\ne0\), with bounded complexity of the composed objects. Such a conclusion is not a formal consequence of \(Sp_4\) monodromy of \(\mathcal F\).

For the actual Apéry problem we are one step further away: \(A_r\) is not currently represented as such a bounded-conductor trace function in the variable \(r\).

So even a genuine \(Sp_4\) theorem would need an additional **index-sheaf / value-distribution theorem** before it could imply the desired zero-set estimate.

---

## 6. Consequence for P3.2

I would record the route as follows.

### Closed as stated

> **`Sp4 monodromy of the four-variable Apéry diagonal => |Z_p| <= p^{1/2+eps}`**

is not a valid deduction.

Reasons:

1. the relevant classical Apéry Picard–Fuchs system is rank 3 / symmetric-square K3, not rank 4 \(Sp_4\);
2. \(Z_p\) counts coefficient zeros, not geometric roots in the parameter \(t\);
3. the mod-\(p\) algebraic cover of the generating series has degree \(\asymp p\), by Caruso et al. 2026, so the hoped-for bounded-degree/bounded-genus hypothesis fails;
4. large monodromy controls trace sums/equidistribution, not exact zero counts without an additional theorem.

### Best defensible bound/status

- I found **no published Apéry-specific theorem** giving \(|Z_p|\ll p^{1/2+\varepsilon}\).
- I found no published theorem specifically for this coefficient-zero set giving an exponent better than \(2/3\).
- In fact, I did not locate the project's claimed \(O(p^{2/3})\) coefficient-zero estimate in the standard literature at all. Therefore it should be retained only to the extent that the repository already has a correct internal proof of it.
- **Conditional on that internal \(O(p^{2/3})\) theorem being correct, it remains the best rigorous bound available for the P3.2 argument from the sources I could verify.** The K3/diagonal/monodromy literature does not currently upgrade it to square root.

### What would genuinely reopen a square-root route

One would need at least one of the following genuinely new inputs:

1. **Index trace function:** construct, uniformly in \(p\), a bounded-rank/bounded-conductor sheaf \(\mathcal G_p\) on the *index line* with
   \[
   A_r\bmod p = T_{\mathcal G_p}(r),\qquad r\in\mathbf F_p,
   \]
   and then prove a strong value-distribution theorem for the trace value 0.

2. **Stepanov/differential coefficient theorem:** exploit the recurrence/differential equation to prove directly that a length-\(p\) solution cannot have more than \(p^{1/2+\varepsilon}\) zero coefficients. Existing Stepanov results for roots of truncated hypergeometric polynomials do not directly apply because they treat zeros in the polynomial variable.

3. **Exploit the 2026 square factorization at coefficient level:** derive a new additive-combinatorial theorem from
   \[
   A_p=B_p^2\quad\text{or}\quad (t^2-34t+1)B_p^2
   \]
   that forces most convolution coefficients to be nonzero. No such theorem is contained in Caruso et al.

4. **A bounded-complexity incidence variety whose projection is exactly \(Z_p\):** not merely a K3 family producing the period sequence, but a variety/sheaf with the discrete coefficient index \(r\) as its \(\mathbf F_p\)-parameter and complexity independent of \(p\).

Until one of these is produced, the geometric origin is extremely useful for congruences and factorization, but it does **not** supply the requested \(p^{1/2+\varepsilon}\) bound.

---

## References

1. F. Beukers and C. A. M. Peters, *A family of K3 surfaces and \(\zeta(3)\)*, J. reine angew. Math. **351** (1984), 42–54. DOI: `10.1515/crll.1984.351.42`.
2. J. Stienstra and F. Beukers, *On the Picard-Fuchs equation and the formal Brauer group of certain elliptic K3-surfaces*, Math. Ann. **271** (1985), 269–304. DOI: `10.1007/BF01455990`.
3. I. Gessel, *Some congruences for Apéry numbers*, J. Number Theory **14** (1982), 362–368. DOI: `10.1016/0022-314X(82)90071-3`.
4. F. Beukers, *Some congruences for the Apéry numbers*, J. Number Theory **21** (1985), 141–155. DOI: `10.1016/0022-314X(85)90047-2`.
5. A. Straub, *Multivariate Apéry numbers and supercongruences of rational functions*, Algebra & Number Theory **8** (2014), 1985–2008. DOI: `10.2140/ant.2014.8.1985`; arXiv:`1401.0854`.
6. A. Malik and A. Straub, *Divisibility properties of sporadic Apéry-like numbers*, Research in Number Theory **2** (2016), article 5. DOI: `10.1007/s40993-016-0036-8`; arXiv:`1508.00297`.
7. A. Ghosh and K. Ward, *Counting roots of truncated hypergeometric series over finite fields*, arXiv:`1601.06765`.
8. X. Caruso, F. Fürnsinn, D. Vargas-Montoya and W. Zudilin, *Galois Groups of Apéry-like Series Modulo Primes*, Bull. Aust. Math. Soc. **114** (2026), 65–78. DOI: `10.1017/S0004972725100932`; arXiv:`2510.23298`.

## Bottom line

**No: the known geometric monodromy of the classical Apéry motive does not currently improve the P3.2 zero-set bound to \(p^{1/2+\varepsilon}\).** The premise is doubly mismatched (rank-3 K3/symmetric-square rather than \(Sp_4\), and coefficient zeros rather than geometric roots). The most recent exact mod-\(p\) geometry gives a degree-\(\asymp p\) algebraic cover and a striking square factorization of the truncation, not a bounded-genus realization of \(Z_p\). Conditional on the repository's \(O(p^{2/3})\) estimate being correct, that remains the best rigorous bound I can verify for this exact P3.2 quantity.