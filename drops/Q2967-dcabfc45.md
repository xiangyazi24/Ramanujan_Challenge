ANSWER Q2967 dcabfc45

# P3.2 exact Gessel-energy height attack

## Verdict

Let

\[
\alpha=17+12\sqrt2,
\qquad
\alpha^{-1}=17-12\sqrt2,
\]

and let

\[
G_n^{(2)}:=\gcd(E_n,b_n^2).
\]

The exact outcome is as follows.

1. The projective-height question is **exactly** a weighted prime-power gcd question:

   \[
   \boxed{
   h_{\mathbf P^1}\!\left(\frac{E_n}{b_n^2}\right)
   =2n\log\alpha-\log G_n^{(2)}+O(\log n).
   }
   \]

   Consequently a positive linear lower bound exists if and only if, for some \(\delta>0\),

   \[
   \boxed{
   \log G_n^{(2)}
   \le (2\log\alpha-\delta)n+O(\log n).
   }
   \tag{HC_\delta}
   \]

   This is the smallest precise missing lemma. It is not a finite-scan statement.

2. There is a useful one-copy reduction:

   \[
   \boxed{
   \gcd(E_n,b_n)\mid G_n^{(2)}\mid
   b_n\gcd(E_n,b_n).
   }
   \]

   Hence the stronger estimate

   \[
   \log\gcd(E_n,b_n)\le(\log\alpha-\delta)n+O(\log n)
   \]

   already implies a positive linear height bound. In particular,
   \(\gcd(E_n,b_n)=\exp(o(n))\) would give slope at least \(\log\alpha\).

3. For every prime \(p>n\) dividing \(b_n\), the first layer of the gcd is exactly a simultaneous Apéry/jet zero:

   \[
   \boxed{
   p\mid E_n\iff p\mid g_n.
   }
   \]

   Away from the central prime \(p=2n+1\), reflection folds this condition to the de-ramified terminal simple-root problem for the finite Racah polynomial. That all-prime simple-root statement remains open in the current project state.

4. The central prime is a genuine exception, not a bookkeeping artifact. If \(p=2n+1\), then \(g_n\equiv0\pmod p\) automatically. The exact first example is \((n,p)=(5,11)\). The banked midpoint supercongruence does, however, show that the central contribution is only \(O(\log n)\) whenever the corresponding weight-4 Fourier coefficient is nonzero. Potentially large central powers are confined to the further exceptional condition \(a_p(f_8)=0\), plus a jet-lift condition.

5. Mod-\(p\) Racah squarefreeness is not the whole height theorem. If a common zero occurs, the exponent in \(\gcd(E_n,b_n^2)\) is governed by an additional \(p\)-adic contact congruence. Primes \(5\le p\le n\) form a separate singular/carry ledger; ordinary Lucas congruences do not control that ledger. The primes \(2\) and \(3\) are completely absent because they never divide an Apéry number.

6. Dividing the target congruence by \(b_{r-1}\) is perfectly legitimate modulo the target prime. It does **not**, however, remove characteristic-zero height. In lowest terms it produces an integer carrier whose size is governed by
   \(h_{\mathbf P^1}(E_{r-1}/b_{r-1}^2)\). Thus one need not multiply back by the raw \(b_{r-1}^2\), but the possible saving is exactly the unresolved gcd in \((HC_\delta)\).

So I do not obtain an unconditional positive linear height theorem from the currently available identities. What is proved below is the full reduction, including prime powers, \(p=2,3\), the central exception, and the exact rational-congruence carrier. It also shows why saying merely “this is equivalent to Racah squarefreeness” would be too strong: squarefreeness resolves only the noncentral clean-prime sector.

## 0. Source boundary

I worked against the connector-visible repository

```text
xiangyazi24/Ramanujan_Challenge
main = 47fa0e653f52c4a71e9a8c26b31ca9f66f6bbe86
```

and the immediately preceding normalized-height audit on `chatgpt-drop`,

```text
drops/Q2961-04938ab9.md
```

at commit `63f4a05e3c387adb24403f2bbe9ee0590c7604f5`.
The caller-local new proof notes are not present on pushed `main`, so the exact recurrence, integrality of \(E_n\), and the displayed target congruence in the question are treated as authenticated current inputs. The pushed `problems/3.2/proof.tex` supplies the Apéry recurrence, Lucas law, reflection, and no-consecutive-zero framework. No project source was edited; only this delivery file was created.

## 1. Exact discrete Lagrange identity

Write

\[
(n+1)^3b_{n+1}=P(n)b_n-n^3b_{n-1},
\qquad
P(n)=34n^3+51n^2+27n+5,
\]

and

\[
(n+1)^3g_{n+1}-P(n)g_n+n^3g_{n-1}=S_n,
\]

where

\[
S_n=P'(n)b_n-3(n+1)^2b_{n+1}-3n^2b_{n-1}.
\]

Define

\[
E_n=n^3(b_{n-1}g_n-b_ng_{n-1}),
\qquad E_1=12.
\]

### Lemma 1 — energy increment

For every \(n\ge1\),

\[
\boxed{E_{n+1}-E_n=b_nS_n.}
\tag{1.1}
\]

### Proof

Subtract \(g_n\) times the homogeneous recurrence from \(b_n\) times the inhomogeneous recurrence:

\[
\begin{aligned}
E_{n+1}
&=(n+1)^3(b_ng_{n+1}-b_{n+1}g_n)\\
&=n^3(b_{n-1}g_n-b_ng_{n-1})+b_nS_n\\
&=E_n+b_nS_n.
\end{aligned}
\]

This proves (1.1). \(\square\)

Eliminating \(b_{n+1}\) gives the smaller exact source

\[
\boxed{
S_n=\frac{3}{n+1}
\left((17n^2+16n+4)b_n-n^2b_{n-1}\right).
}
\tag{1.2}
\]

The current Gessel construction proves that this value is integral. Hence

\[
E_n=12+\sum_{j=1}^{n-1}b_jS_j\in\mathbf Z.
\tag{1.3}
\]

There is one exact adjacent consequence:

\[
\boxed{
\gcd(E_{n+1},b_n)=\gcd(E_n,b_n).
}
\tag{1.4}
\]

This is the strongest immediate adjacent Bézout identity. It does not compare the gcd with \(b_{n+1}\), so it does not descend or bound the height by itself.

It is also useful to put

\[
h_n:=\frac{g_n}{b_n}.
\]

Then

\[
\boxed{
K_n:=\frac{E_n}{b_n^2}
=n^3\frac{b_{n-1}}{b_n}(h_n-h_{n-1}).
}
\tag{1.5}
\]

This separates the ordinary real size of \(K_n\) from its arithmetic cancellation.

## 2. The analytic size is polynomial; the arithmetic height is a gcd

The standard uniform saddle-point expansion of the entire Apéry interpolation gives

\[
b_n=C\alpha^n n^{-3/2}(1+O(n^{-1})),
\tag{2.1}
\]

and, after differentiating the same expansion in the index variable,

\[
\frac{g_n}{b_n}
=\log\alpha-\frac{3}{2n}+O(n^{-2}).
\tag{2.2}
\]

Equations (1.5), (2.1), and (2.2) yield the analytic asymptotic stated in the question:

\[
\boxed{
K_n=\frac32(17-12\sqrt2)n+O(1).
}
\tag{2.3}
\]

Indeed,

\[
\frac{b_{n-1}}{b_n}=\alpha^{-1}(1+O(n^{-1})),
\qquad
h_n-h_{n-1}=\frac{3}{2n(n-1)}+O(n^{-3}).
\]

No arithmetic conclusion follows from (2.3) alone. It only shows that
\(\log^+|K_n|=O(\log n)\).

Now put

\[
G_n^{(2)}=\gcd(E_n,b_n^2),
\]

and write

\[
K_n=\frac{A_n}{B_n},
\qquad
\gcd(A_n,B_n)=1,\quad B_n>0.
\]

Then exactly

\[
A_n=\frac{E_n}{G_n^{(2)}},
\qquad
B_n=\frac{b_n^2}{G_n^{(2)}}.
\tag{2.4}
\]

For the logarithmic projective height

\[
h_{\mathbf P^1}(K_n)=\log\max(|A_n|,B_n),
\]

we therefore have

\[
\boxed{
\begin{aligned}
h_{\mathbf P^1}(K_n)
&=2\log b_n-\log G_n^{(2)}+\log\max(1,|K_n|)\\
&=2n\log\alpha-\log G_n^{(2)}+O(\log n).
\end{aligned}
}
\tag{2.5}
\]

In particular,

\[
\boxed{
\liminf_{n\to\infty}\frac{h_{\mathbf P^1}(K_n)}n
=2\log\alpha-
\limsup_{n\to\infty}\frac{\log G_n^{(2)}}n.
}
\tag{2.6}
\]

Thus a positive linear lower bound is equivalent to \((HC_\delta)\). The empirical numerator and denominator sizes are evidence for such a bound, but they are not part of this proof.

A warning is important here. The approximation

\[
K_n=cn+O(1)
\]

cannot itself force exponential height. Rational numbers of small height can approximate a linear real quantity with bounded absolute error. One would need an exponentially accurate approximation to a fixed irrational or transcendental number before a Liouville/Baker argument could enter; (2.3) supplies nothing of that kind.

## 3. One-copy gcd reduction

Put

\[
G_n^{(1)}:=\gcd(E_n,b_n).
\]

### Lemma 2 — square-gcd sandwich

For every \(n\),

\[
\boxed{
G_n^{(1)}\mid G_n^{(2)}\mid b_nG_n^{(1)}.
}
\tag{3.1}
\]

### Proof

For a prime \(p\), set

\[
a=v_p(b_n),\qquad e=v_p(E_n).
\]

The three exponents in (3.1) are

\[
\min(e,a),\qquad
\min(e,2a),\qquad
a+\min(e,a).
\]

The inequalities

\[
\min(e,a)\le\min(e,2a)\le a+\min(e,a)
\]

prove the divisibilities prime by prime. \(\square\)

Since the reduced denominator is \(b_n^2/G_n^{(2)}\), (3.1) gives

\[
\boxed{
\operatorname{den}(K_n)
\ge\frac{b_n}{G_n^{(1)}}.
}
\tag{3.2}
\]

Consequently

\[
h_{\mathbf P^1}(K_n)
\ge \log b_n-\log G_n^{(1)}.
\tag{3.3}
\]

This proves the useful sufficient criterion

\[
\limsup\frac{\log G_n^{(1)}}n<\log\alpha
\quad\Longrightarrow\quad
h_{\mathbf P^1}(K_n)=\Omega(n).
\tag{3.4}
\]

It does not prove (3.4)'s hypothesis. The problem has merely been reduced from a square gcd to a one-copy common-content gcd.

## 4. Clean primes: exact common-zero criterion

Let \(p>n\) be prime and suppose \(p\mid b_n\).

The Racah realization is

\[
\Phi_k(Y)=\frac1{(k!)^2}
\prod_{s=0}^{k-1}(Y-s(s+1)),
\qquad
R_n(Y)=\sum_{k=0}^n\Phi_k(Y)^2.
\tag{4.1}
\]

At the terminal triangular node \(\lambda_n=n(n+1)\),

\[
\boxed{
R_n(\lambda_n)=b_n,
\qquad
(2n+1)R_n'(\lambda_n)=g_n.
}
\tag{4.2}
\]

Moreover

\[
(n!)^4R_n(Y)\in\mathbf Z[Y].
\tag{4.3}
\]

Thus \(g_n\) and \(g_{n-1}\) are \(p\)-integral for every \(p>n\). The Apéry no-consecutive-zero argument applies because all recurrence coefficients up to row \(n\) are \(p\)-units:

\[
p\mid b_n\quad\Longrightarrow\quad p\nmid b_{n-1}.
\tag{4.4}
\]

Since also \(p\nmid n\), reducing

\[
E_n=n^3(b_{n-1}g_n-b_ng_{n-1})
\]

modulo \(p\) proves

\[
\boxed{
 p\mid b_n
 \quad\Longrightarrow\quad
 \bigl(p\mid E_n\iff p\mid g_n\bigr),
 \qquad p>n.
}
\tag{4.5}
\]

This is the exact clean-prime support theorem.

### Reflection and the de-ramified terminal problem

For \(0\le j<p\), the current reflection laws are

\[
b_{p-1-j}\equiv b_j\pmod p,
\qquad
g_{p-1-j}\equiv-g_j\pmod p.
\tag{4.6}
\]

If \(p>n\), \(p\ne2n+1\), and \(p\mid b_n,g_n\), put

\[
j=\min(n,p-1-n).
\]

Then

\[
0\le j<\frac{p-1}{2},
\qquad p>2j+1,
\qquad p\mid b_j,g_j.
\tag{4.7}
\]

Conversely any such strict common zero produces a clean common zero at its reflected index.

Because \(p>2j+1\), the factor \(2j+1\) and every denominator in (4.2) are \(p\)-units. Hence

\[
\boxed{
 p\mid b_j,\quad p\mid g_j
 \iff
 R_j(\lambda_j)\equiv R_j'(\lambda_j)\equiv0\pmod p.
}
\tag{4.8}
\]

The relevant pointwise theorem is therefore:

> **FJT — de-ramified first-jet transversality.** For every \(j\ge0\) and every prime \(p>2j+1\),
> \[
> p\mid b_j\Longrightarrow p\nmid g_j.
> \]
> Equivalently, the terminal triangular node \(\lambda_j\) is never a multiple root modulo \(p\) of the cleared finite Racah polynomial when it is a root.

FJT is not proved by the current recurrence, reflection, Lucas, or Hasse-polynomial inputs. It is precisely the open noncentral clean-radical seam.

A logical distinction matters:

- FJT would force the entire noncentral clean contribution to \(G_n^{(2)}\) to vanish, including every prime power, because (4.5) would already give \(v_p(E_n)=0\).
- The projective-height theorem does **not** imply FJT; finitely or sparsely many clean common zeros could coexist with positive linear height.
- Conversely FJT alone does **not** imply the height theorem, because it says nothing about primes \(p\le n\) or about the exceptional central lift.

Thus the height question is equivalent to a weighted gcd theorem, while FJT is a strong pointwise solution to one sector of that theorem.

## 5. The central prime and its prime powers

If

\[
p=2n+1
\]

is prime, then (4.2) gives

\[
\boxed{g_n\equiv0\pmod p.}
\tag{5.1}
\]

This is the ramification factor of the triangular coordinate. Hence

\[
p\mid b_n\quad\Longrightarrow\quad p\mid E_n.
\tag{5.2}
\]

The first exact instance is

\[
n=5,\qquad p=11,
\]

with

\[
b_4=33001,\quad b_5=819005,\quad
 g_4=104825,\quad g_5=\frac{13276637}{5}.
\]

Here

\[
E_5
=25\left(33001\cdot13276637
      -5\cdot819005\cdot104825\right)
\equiv66\pmod{121}.
\]

Therefore

\[
\boxed{
v_{11}(b_5)=v_{11}(E_5)=1.
}
\tag{5.3}
\]

This is an exact counterfixture to any claim that all primes above the row are excluded from the gcd.

There is, however, a useful midpoint refinement. Let

\[
f_8(\tau)=\eta(2\tau)^4\eta(4\tau)^4
=\sum_{m\ge1}a_mq^m.
\]

The banked midpoint supercongruence is

\[
 b_{(p-1)/2}\equiv a_p(f_8)\pmod{p^2}.
\tag{5.4}
\]

For the normalized weight-4 newform,

\[
|a_p(f_8)|\le2p^{3/2}<p^2
\qquad(p\ge5).
\tag{5.5}
\]

Suppose \(p=2n+1\), \(p\mid b_n\), and \(a_p(f_8)\ne0\). Then \(p\mid a_p(f_8)\), but \(p^2\nmid a_p(f_8)\) by (5.5). Equation (5.4) therefore gives

\[
\boxed{v_p(b_n)=1.}
\tag{5.6}
\]

Consequently the central contribution to \(\log G_n^{(2)}\) is at most

\[
2\log(2n+1)=O(\log n)
\tag{5.7}
\]

whenever \(a_p(f_8)\ne0\).

The only central case in which the denominator can carry more than \(p^2\) is thus confined to

\[
\boxed{a_p(f_8)=0,}
\tag{5.8}
\]

followed by a higher jet/contact condition. Uniform nonvanishing (5.8) is not a currently available all-prime theorem for this non-CM form, and (5.4) gives no bound on \(v_p(b_n)\) once \(a_p=0\). The earlier statement “there is only one central prime, so it is harmless” is therefore valid at radical level, but not by itself at prime-power height level. Equations (5.4)–(5.7) identify the precise harmless generic case.

## 6. Exact clean-prime valuation formula

Let \(p>n\), \(p\mid b_n\), and put

\[
a=v_p(b_n)>0,
\qquad c=v_p(g_n),
\qquad d=v_p(g_{n-1}).
\]

All these valuations are nonnegative because the jet values are \(p\)-integral. Write

\[
b_{n-1}=u,\quad b_n=p^aB,\quad
 g_n=p^cG,\quad g_{n-1}=p^dH,
\]

where \(u,B,G,H\) are \(p\)-units whenever the corresponding value is nonzero. Since \(n\) is a unit,

\[
v_p(E_n)=v_p\left(up^cG-p^{a+d}BH\right).
\]

Therefore

\[
\boxed{
 v_p(E_n)=
 \begin{cases}
 c,&c<a+d,\\
 a+d,&c>a+d,\\
 c+v_p(uG-BH),&c=a+d.
 \end{cases}
}
\tag{6.1}
\]

The definitive gcd exponent is

\[
\boxed{
 v_p(G_n^{(2)})
 =\min\bigl(v_p(E_n),2a\bigr).
}
\tag{6.2}
\]

Equations (6.1)–(6.2) explicitly separate the two levels:

- if \(c=0\), there is no cancellation at any power of \(p\);
- if \(c>0\), the mod-\(p\) common-zero event has occurred;
- if \(c=a+d\), an additional normalized congruence
  \[
  uG\equiv BH\pmod{p^t}
  \]
  can raise the energy valuation by \(t\).

Thus a count of common zeros only controls the support. A square-gcd height theorem must also control these Bockstein/contact exponents. Formula (6.1) handles arbitrary prime powers; no squarefree assumption has been made.

## 7. Small primes and the singular ledger

The primes \(2\) and \(3\) can be removed exactly.

Gessel's Lucas law gives

\[
b_n\equiv\prod_i b_{n_i}\pmod p
\]

for the base-\(p\) digits \(n_i\). For \(p=2\),

\[
b_0\equiv b_1\equiv1\pmod2,
\]

so every \(b_n\) is odd. For \(p=3\),

\[
b_0\equiv1,\qquad b_1\equiv2,\qquad b_2\equiv1\pmod3,
\]

so every digit factor is nonzero. Hence

\[
\boxed{2\nmid b_n,\qquad3\nmid b_n\quad\text{for every }n.}
\tag{7.1}
\]

Therefore neither prime contributes to \(G_n^{(2)}\), irrespective of jet denominators.

For \(5\le p\le n\), the clean argument fails for structural reasons:

1. recurrence steps at multiples of \(p\) are singular;
2. the rational jet may have a \(p\)-denominator even though \(E_n\) is integral;
3. ordinary Lucas determines whether \(p\mid b_n\), but not \(v_p(b_n)\), \(v_p(E_n)\), or the cancellation between the two terms defining \(E_n\);
4. the mod-\(p^2\) Gessel formula introduces the first jet as a correction term—it does not prove that the correction is nonzero;
5. higher valuation-keyed Lucas lifts are known to encounter genuine memory terms, so one cannot promote a mod-\(p\) digit statement to arbitrary prime powers.

Define the exact singular content

\[
\boxed{
\Sigma_n
:=\sum_{5\le p\le n}
\min\bigl(v_p(E_n),2v_p(b_n)\bigr)\log p.
}
\tag{7.2}
\]

No current Lucas or lcm theorem gives the required uniform linear saving for \(\Sigma_n\). These are not “finitely many harmless primes”: the set grows with \(n\), and its weighted total can in principle be linear.

## 8. Exact decomposition and the smallest missing lemma

For every prime define

\[
\gamma_p(n)=\min\bigl(v_p(E_n),2v_p(b_n)\bigr).
\]

Besides \(\Sigma_n\), put

\[
\Xi_n
:=\sum_{\substack{p>n\\p\ne2n+1}}
\gamma_p(n)\log p,
\tag{8.1}
\]

and

\[
M_n
:=
\begin{cases}
\gamma_{2n+1}(n)\log(2n+1),&2n+1\text{ prime},\\
0,&\text{otherwise}.
\end{cases}
\tag{8.2}
\]

By (7.1),

\[
\boxed{
\log G_n^{(2)}=\Sigma_n+\Xi_n+M_n.
}
\tag{8.3}
\]

Combining (2.5) and (8.3) gives the exact theorem

\[
\boxed{
 h_{\mathbf P^1}(K_n)
 =2n\log\alpha-\Sigma_n-\Xi_n-M_n+O(\log n).
}
\tag{8.4}
\]

This is the complete prime-power reduction.

The smallest missing statement is exactly

> **JET-CONTENT\((\delta)\).** There exists \(\delta>0\) such that
> \[
> \Sigma_n+\Xi_n+M_n
> \le(2\log\alpha-\delta)n+O(\log n)
> \]
> for every \(n\).

By (8.4), JET-CONTENT\((\delta)\) is equivalent to the desired positive linear projective-height bound.

A more structured sufficient package is:

1. de-ramified FJT, which gives \(\Xi_n=0\);
2. the midpoint nonzero-coefficient/contact estimate \(M_n=o(n)\), with (5.7) proving it whenever \(a_{2n+1}(f_8)\ne0\);
3. a singular-prime estimate
   \[
   \Sigma_n\le(2\log\alpha-\delta)n+O(\log n).
   \]

The third item remains even after perfect clean-prime transversality. Conversely, one could replace pointwise FJT by the weaker weighted estimate \(\Xi_n=o(n)\). This is why the exact height theorem is a weighted gcd theorem rather than literally the pointwise common-zero conjecture.

## 9. Audit of the proposed mechanisms

### 9.1 Adjacent energy and continuants

Equation (1.4) is exact, but it fixes the middle Apéry value:

\[
\gcd(E_{n+1},b_n)=\gcd(E_n,b_n).
\]

It does not compare \(\gcd(E_n,b_n)\) with
\(\gcd(E_{n-1},b_{n-1})\) or
\(\gcd(E_{n+1},b_{n+1})\). Hence it supplies no induction on the common content.

The homogeneous continuant controls \(b_n\) and the determinant of a transfer block. The jet is an inhomogeneous extension, and the source term in (1.1) prevents its Wronskian from being constant. A constant-Casoratian Bézout argument therefore does not apply.

### 9.2 Racah Bézout and discriminants

Let

\[
F_n(Y)=(n!)^4R_n(Y)\in\mathbf Z[Y].
\]

A polynomial Bézout identity gives

\[
A_n(Y)F_n(Y)+B_n(Y)F_n'(Y)
=\operatorname{Res}(F_n,F_n').
\]

Evaluating at \(\lambda_n\) shows that any noncentral clean common divisor divides the resultant/discriminant. This is true but insufficient:

- the resultant counts every repeated root, not the distinguished terminal node;
- its generic coefficient-height bound is far larger than the needed linear budget;
- squarefreeness over \(\mathbf Q\) does not exclude new repeated roots after reduction modulo primes depending on \(n\);
- no explicit Bézout evaluation of logarithmic size below \(2n\log\alpha\) is currently available.

The adjacent Racah identity

\[
R_n(Y)=R_{n-1}(Y)+\Phi_n(Y)^2
\tag{9.1}
\]

also does not contradict a terminal multiple root. At \(Y=\lambda_n\),

\[
\Phi_n(\lambda_n)=\binom{2n}{n}.
\]

For \(p>2n\) this is a unit, and a common zero merely forces

\[
R_{n-1}(\lambda_n)=-\Phi_n(\lambda_n)^2,
\qquad
R_{n-1}'(\lambda_n)=-2\Phi_n(\lambda_n)\Phi_n'(\lambda_n),
\]

which is consistent. Thus (9.1) transports the condition but does not prove transversality.

### 9.3 Lucas digit laws

Ordinary Lucas is decisive for (7.1) and for locating the support of small-prime divisibility in \(b_n\). It is silent about the jet and about prime powers. Gessel's mod-\(p^2\) refinement shows that the derivative is exactly the first obstruction to lifting Lucas; it does not make that derivative nonzero. At deeper powers, predecessor/memory terms enter. Hence Lucas identifies the singular ledger rather than bounding it.

### 9.4 Lcm denominators

The banked denominator theorem places the denominator of \(g_n\) in a power of the small lcm (in the usual normalization, it divides \(L_{2n}\)). Clearing this denominator does not control common numerator factors with \(b_n\).

The exact central example already disproves any claim that the gcd is supported only on that lcm:

\[
11\mid\gcd(E_5,b_5^2),
\qquad
11\nmid L_{10}.
\]

Thus denominator clearing is useful for defining reductions modulo clean primes, but it gives no upper bound for \(G_n^{(2)}\).

## 10. The divided target congruence

Let the polynomial denoted by \(D_r\) in the question be the current low-degree source polynomial. In the Q2961 normalization it is

\[
D_r=17r^2-16r+4.
\]

The underlying characteristic-zero identity is

\[
\boxed{
\begin{aligned}
&(r-1)r^3b_{r-1}g_r
-(r-1)E_{r-1}
+3D_rb_{r-1}^2\\
&\qquad=
 b_r\left((r-1)r^3g_{r-1}+3r^2b_{r-1}\right).
\end{aligned}
}
\tag{10.1}
\]

Let \(p>r\) and \(p\mid b_r\). No consecutive zeros give

\[
p\nmid b_{r-1},
\]

and \(p\nmid r(r-1)\). Therefore division in \(\mathbf F_p\) is legitimate and (10.1) gives exactly

\[
\boxed{
\frac{E_{r-1}}{b_{r-1}^2}
\equiv
r^3\frac{g_r}{b_{r-1}}
+\frac{3D_r}{r-1}\pmod p.
}
\tag{10.2}
\]

So the answer to the local question is **yes**: one may use (10.2) without multiplying by the raw exponential integer \(b_{r-1}\).

The global height question is different. Write

\[
\frac{E_{r-1}}{b_{r-1}^2}=\frac{A}{B}
\]

in lowest terms. Since \(p\nmid b_{r-1}\), one also has \(p\nmid B\). On the common-zero subevent \(p\mid g_r\), (10.2) becomes

\[
\boxed{
 p\mid (r-1)A-3D_rB.
}
\tag{10.3}
\]

This is a valid characteristic-zero integer carrier, and it clears only the **reduced** denominator \(B\), not \(b_{r-1}^2\).

The following elementary height lemma shows exactly what has and has not been gained.

### Lemma 3 — rational congruence carrier

Let \(x=A/B\) and \(c=u/v\) be reduced rationals with positive denominators. Let \(\mathcal S\) be a finite set of primes satisfying

\[
p\nmid Bv,\qquad x\equiv c\pmod p.
\]

If \(x\ne c\), then

\[
\boxed{
\prod_{p\in\mathcal S}p\mid vA-uB,
}
\tag{10.4}
\]

and

\[
\boxed{
\sum_{p\in\mathcal S}\log p
\le h_{\mathbf P^1}(x)+h_{\mathbf P^1}(c)+\log2.
}
\tag{10.5}
\]

### Proof

Each congruence gives \(p\mid vA-uB\); distinct primes multiply. Moreover

\[
|vA-uB|
\le2\max(|A|,B)\max(|u|,v),
\]

which proves (10.5). \(\square\)

Applying the lemma to (10.3) shows that the logarithmic cost of the divided carrier is

\[
h_{\mathbf P^1}\!\left(\frac{E_{r-1}}{b_{r-1}^2}\right)+O(\log r).
\tag{10.6}
\]

Thus local inversion does not magically remove global height. It moves the issue from the raw denominator \(b_{r-1}^2\) to the reduced denominator

\[
B=rac{b_{r-1}^2}{\gcd(E_{r-1},b_{r-1}^2)}.
\]

Whether this is exponentially large is exactly \((HC_\delta)\).

If one instead multiplies (10.2) by \(b_{r-1}\), then

\[
b_{r-1}K_{r-1}=rac{E_{r-1}}{b_{r-1}},
\]

so the old raw-height carrier is literally restored. The reduced formulation (10.3) avoids that unnecessary multiplication, but it does not avoid the unresolved projective height.

## 11. Final theorem-quality reduction

Collecting the preceding sections gives the following precise statement.

### Theorem 4 — height/common-content equivalence and scoped obstruction

Assume the analytic expansion (2.1)–(2.3). For every \(n\ge1\), define \(\Sigma_n,\Xi_n,M_n\) by (7.2), (8.1), and (8.2). Then

\[
h_{\mathbf P^1}\!\left(\frac{E_n}{b_n^2}\right)
=2n\log(17+12\sqrt2)-\Sigma_n-\Xi_n-M_n+O(\log n).
\]

Moreover:

1. \(p=2,3\) never contributes.
2. De-ramified FJT implies \(\Xi_n=0\).
3. If \(2n+1=p\) is prime, \(p\mid b_n\), and \(a_p(f_8)\ne0\), then \(M_n\le2\log p\).
4. At every clean prime, the exact prime-power exponent is given by (6.1)–(6.2).
5. The desired positive linear height lower bound is equivalent to JET-CONTENT\((\delta)\) for some \(\delta>0\).

Hence the smallest live proof target is not another finite computation and not merely a mod-\(p\) resultant. It is one of the following equivalent-strength arithmetic advances:

- prove JET-CONTENT\((\delta)\) directly;
- prove a one-copy common-content saving for \(\gcd(E_n,b_n)\) strong enough for (3.4);
- prove de-ramified FJT, control the exceptional central Hecke-zero lift, and prove the singular ledger bound for \(5\le p\le n\).

The current identities prove all reductions to these statements but do not prove any one of the required uniform weighted bounds. In particular, the exact computation through \(n=800\) is fully consistent with a positive linear height theorem, but it cannot replace JET-CONTENT\((\delta)\).
