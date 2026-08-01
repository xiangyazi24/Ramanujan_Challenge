PARTIAL: [NO-RIGHT-2-3] PROVED FOR ALL H; WALL NOT BROKEN, NEW RESIDUAL [ZERO-TAIL-2]

# Breakwall report

## 1. Verdict

No condition W1--W7 is proved here. In particular, no unbounded choice of
\(L(p)\) is proved to improve the \(3/2\) energy record.

There are nevertheless two rigorous new conclusions.

1. **Characteristic-zero progress toward W7.** For every \(h\geq2\), the
   Apéry rational function
   \[
   f_h(X)=\frac{N_h(X)}{q_h(X)^3},\qquad
   q_h(X)=\prod_{j=1}^h(X+j),
   \]
   has no nontrivial decomposition \(f_h=g\circ u\) over
   \(\overline{\mathbb Q}\) with right degree
   \(\deg u\in\{2,3\}\). Equivalently, a block system for the geometric
   monodromy of \(f_h\) cannot have block size \(2\) or \(3\). This is
   [NO-RIGHT-2-3]. It is an all-\(h\) theorem, but it does not exclude block
   sizes at least \(4\), so it does not prove W7.

2. **The proposed L2 census bridge is impossible from row marginals alone.**
   Even if every gap \(h\leq D\) is granted a value histogram of total mass
   \(N-h\), value-collision moment at most \(2(N-h)\), and zero-fibre size at
   most \(3(h-1)\), the zero-fibre total can still be
   \(\Theta(NL)\) at \(D=\sqrt N L\). This remains true if such an L2
   estimate is granted for every \(h\leq D\); hence no growth rate for the
   certificate cutoff \(h_0(p)\leq D\) can by itself tip direction A.

The precise replacement isolated by that obstruction is [ZERO-TAIL-2].
It asks only for a weak-L2 tail bound on the exceptional lags with unusually
large zero fibres. It is far weaker in scale than [MESO-S1-2/3] and the
banked [SAME-LAG-BDH], yet it implies W1.

## 2. Banked facts used

Only the following already proved facts about \(f_h\) enter the
characteristic-zero argument.

- \(\deg f_h=3h\), and its finite poles are exactly
  \(-1,\ldots,-h\), each of order \(3\).
- \(N_h\) is squarefree and coprime to \(q_h\). Thus the zero fibre of
  \(f_h\) consists of \(3h-3\) simple finite zeros and the single point
  \(X=\infty\) with local degree \(3\).
- For \(h=3\), the recurrence gives
  \[
  \begin{split}
  N_3(X)={}&1155X^6+13860X^5+68535X^4+178680X^3\\
           &+259059X^2+198156X+62531.
  \end{split}
  \tag{2.1}
  \]

No critical-value hypothesis, finite certificate, or unproved
irreducibility statement is used below.

## 3. An all-height right-degree gap

### Theorem 3.1 ([NO-RIGHT-2-3])

For every \(h\geq2\), there do not exist rational functions
\(g,u\in\overline{\mathbb Q}(X)\), both of degree greater than one, such that
\[
f_h=g\circ u
\quad\hbox{and}\quad
\deg u\in\{2,3\}.
\tag{3.1}
\]

#### Proof

Assume \(f_h=g\circ u\), and write \(a=\deg u\). Local degrees multiply
under composition.

**Step 1: the zero fibre forces \(e_u(\infty)=3\).**
Every zero of \(g\) is simple. Indeed, if a zero \(y\) of \(g\) had local
degree at least \(2\), every point of \(u^{-1}(y)\) would be a multiple zero
of \(f_h\). The only such zero is \(\infty\), so that fibre would have to be
supported at \(\infty\). It would then give
\[
e_{f_h}(\infty)=e_g(y)e_u(\infty)=e_g(y)a\geq4,
\]
contrary to \(e_{f_h}(\infty)=3\).

Put \(y_0=u(\infty)\). Since \(g(y_0)=0\), the preceding paragraph gives
\(e_g(y_0)=1\), and hence
\[
e_u(\infty)=3.
\tag{3.2}
\]
This already excludes \(a=2\).

It remains to treat \(a=3\). Then \(u\) is totally ramified at infinity.
Since \(\deg f_h=3h\), the outer degree is \(\deg g=h\).

**Step 2: classify the pole fibres of \(g\).**
Let \(\beta\) be a pole of \(g\), of order \(m\). Every point
\(x\in u^{-1}(\beta)\) is one of the poles \(-1,\ldots,-h\) of \(f_h\), and
\[
m e_u(x)=3.
\tag{3.3}
\]
Consequently \(m\in\{1,3\}\). Above a simple pole of \(g\), \(u\) has one
totally ramified point; above an order-three pole of \(g\), \(u\) has three
distinct unramified points.

Let \(k\) and \(\ell\) be the numbers of order-one and order-three poles of
\(g\). The pole divisor of a degree-\(h\) rational function has degree \(h\),
so
\[
k+3\ell=h.
\tag{3.4}
\]
Riemann--Hurwitz gives total ramification defect \(2\cdot3-2=4\) for \(u\).
Infinity already contributes \(2\), and every simple pole of \(g\) forces
another totally ramified point and contributes \(2\). Hence \(k\leq1\).

**Step 3: the case \(h\equiv2\pmod3\).**
Equation (3.4) would require \(k\equiv2\pmod3\), which is incompatible with
\(0\leq k\leq1\). Thus this case is impossible.

**Step 4: the case \(h\equiv1\pmod3\).**
Here \(k=1\). Let \(x_0\in\{-1,\ldots,-h\}\) be the totally ramified point
over the unique simple pole \(\beta_0\) of \(g\). Choose a target Möbius
map \(M\) with \(M(y_0)=0\) and \(M(\beta_0)=\infty\). The degree-three map
\(v=M\circ u\) has divisor
\[
\operatorname{div}(v)=3[\infty]-3[x_0],
\]
and therefore
\[
v(X)=\frac{c}{(X-x_0)^3}\qquad(c\ne0).
\tag{3.5}
\]

Since \(h\geq2\) and \(h\equiv1\pmod3\), in fact \(h\geq4\) and
\(\ell=(h-1)/3\geq1\). Choose an order-three pole \(\beta\) of \(g\).
Its three preimages are distinct members of the real set
\(\{-1,\ldots,-h\}\). On the other hand, (3.5) says that they are the
three roots of
\[
(X-x_0)^3=c/M(\beta).
\]
Three distinct roots of this equation cannot all be real: after subtracting
their centroid \(x_0\), their nonzero differences are proportional to
\(1,\omega,\omega^2\), where \(\omega\) is a nonreal cube root of unity.
This is a contradiction.

**Step 5: the case \(h\equiv0\pmod3\).**
Now \(k=0\) and \(\ell=h/3\). Choose one pole \(\beta_\infty\) of \(g\), and
choose \(M\) with \(M(y_0)=0\), \(M(\beta_\infty)=\infty\). The map
\(v=M\circ u\) has a zero of order \(3\) at infinity. Writing it in lowest
terms and using \(\deg v=3\) gives
\[
v(X)=\frac{c}{B(X)},\qquad \deg B=3,\quad c\ne0.
\tag{3.6}
\]

For each pole \(\beta_i\) of \(g\), its three preimages form one block of
the set \(\{-1,\ldots,-h\}\). For \(\beta_\infty\) the block polynomial is
\(B(X)\); for every other pole it is \(B(X)-c/M(\beta_i)\). After making
these cubics monic, all have the same \(X^2\)- and \(X\)-coefficients.
Changing signs from roots \(-j\) to labels \(j\), it follows that the
\(h\) labels \(1,\ldots,h\) admit a partition into triples, every triple
having the same sum \(\sigma\) and the same sum of squares \(\tau\).
Summing over the \(h/3\) blocks forces
\[
\sigma=\frac{3(h+1)}2,
\qquad
\tau=\frac{(h+1)(2h+1)}2.
\tag{3.7}
\]

If \(h\) is even, \(\sigma\) is not an integer, a contradiction. Suppose
\(h\) is odd. If the block containing \(1\) is \(\{1,r,s\}\), then
\[
r+s=\sigma-1,
\qquad r^2+s^2=\tau-1,
\]
and hence
\[
(r-s)^2
=2(\tau-1)-(\sigma-1)^2
=-\frac{(h-1)(h-5)}4.
\tag{3.8}
\]
This is negative for \(h\geq9\). The only remaining multiple of three is
\(h=3\).

For \(h=3\), (3.6) and the sole pole fibre give, after rescaling,
\[
v(X)=\frac{c}{q_3(X)},\qquad q_3(X)=(X+1)(X+2)(X+3).
\]
The transformed outer map \(g\circ M^{-1}\) has its only pole at infinity,
of order \(3\), and has three simple zeros, one at \(0\). It is therefore a
cubic polynomial \(az(z-\alpha)(z-\beta)\). Clearing \(q_3^3\) in the
composition would force
\[
N_3(X)\in\operatorname{span}_{\overline{\mathbb Q}}
\{1,q_3(X),q_3(X)^2\}.
\tag{3.9}
\]
But the \(X^6\)-coefficient in (2.1) would make the coefficient of
\(q_3^2\) equal to \(1155\), while
\[
[X^4]q_3^2=58,
\qquad 1155\cdot58=66990\ne68535=[X^4]N_3.
\]
This contradicts (3.9) and finishes every case. \(\square\)

### Corollary 3.2

Every nontrivial block system for the geometric monodromy of \(f_h\) has
block size at least \(4\). This follows from the intermediate-field theorem
and Lüroth's theorem: a block system of size \(a\) gives a right component of
degree \(a\).

This is genuine but insufficient progress toward W7. Full primitivity still
requires excluding all block sizes \(a\geq4\).

## 4. Why Theorem 3.1 stops exactly at degree four

The pole-divisor argument above cannot simply be extended to \(a=4\).
There is an exact decomposable comparison cover with all the following Apéry
signatures:

- eight consecutive finite poles, each of order \(3\);
- a squarefree numerator of degree \(21=3\cdot8-3\);
- a zero of order \(3\) at infinity;
- the same odd reflection law around \(-9/2\).

Put
\[
\begin{aligned}
P(X)&=(X+1)(X+4)(X+6)(X+7),\\
Q(X)&=(X+2)(X+3)(X+5)(X+8),\\
A(X)&=P(X)-Q(X)=-16X-72,\\
u(X)&=A(X)/P(X),\\
\phi(z)&=-z/(1-z),\\
g(z)&=\frac{z(z-4)(z-3)(z-2)(2z-3)(3z-4)}
             {16(z-1)^3}.
\end{aligned}
\tag{4.1}
\]
The two four-element blocks have equal first two moments:
\[
1+4+6+7=2+3+5+8=18,
\]
\[
1^2+4^2+6^2+7^2=2^2+3^2+5^2+8^2=102.
\]
This is why \(P-Q\) has degree only one and \(u\) has a triple zero at
infinity.

For \(\iota(X)=-9-X\), direct algebra gives
\[
P\circ\iota=Q,\qquad u\circ\iota=\phi\circ u,
\qquad g\circ\phi=-g.
\tag{4.2}
\]
Consequently \(F=g\circ u\) satisfies \(F\circ\iota=-F\). Exact
cancellation gives
\[
F(X)=\frac{M(X)}{\prod_{j=1}^8(X+j)^3},
\qquad \deg M=21,\qquad \gcd(M,M')=1.
\tag{4.3}
\]
Yet \(F\) has the displayed right component \(u\) of degree \(4\).

This comparison is not the Apéry function \(f_8\), and hence is not a
counterexample to W7. It is a sharp counterexample to any attempted proof
of full indecomposability that uses only the divisor multiplicities,
consecutive pole locations, numerator squarefreeness, and reflection.
Additional Apéry-specific information is unavoidable from degree four on.

## 5. Optimal obstruction to the L2 census route

For this section let \(p=N+1\) be odd, let \(1\leq D\leq N/2\), and put
\(M_h=N-h\). Abstract a gap-\(h\) value distribution as a histogram
\(\nu_h:\mathbb F_p\to\mathbb Z_{\geq0}\) with
\[
\sum_a\nu_h(a)=M_h,
\qquad R_h=\nu_h(0),
\qquad \operatorname{Col}(h)=\sum_a\nu_h(a)^2.
\]

### Theorem 5.1 (sharp marginal countermodel)

There exist such histograms simultaneously for every \(1\leq h\leq D\)
with
\[
R_h\leq3(h-1),
\qquad
\operatorname{Col}(h)\leq2M_h,
\tag{5.1}
\]
but, if \(D=\lfloor\sqrt N L\rfloor\), \(L\to\infty\), and
\(D\leq N/2\), then
\[
\sum_{h\leq D}R_h=\Theta(NL).
\tag{5.2}
\]

#### Proof

Set
\[
k_h=\min\{3(h-1),\lfloor\sqrt{M_h}\rfloor\}.
\tag{5.3}
\]
Put \(k_h\) objects in the zero cell and put each of the remaining
\(M_h-k_h\) objects in a distinct nonzero cell. This is possible because
\(M_h-k_h\leq N=p-1\). Then \(R_h=k_h\leq3(h-1)\), while
\[
\operatorname{Col}(h)=k_h^2+M_h-k_h\leq2M_h.
\tag{5.4}
\]

For
\(h\geq\lceil\sqrt N/3\rceil+1\), the degree cap in (5.3) is at least
\(\sqrt N\). Since \(h\leq D\leq N/2\), also
\(\lfloor\sqrt{M_h}\rfloor\geq\lfloor\sqrt{N/2}\rfloor\). Therefore
\[
\sum_{h\leq D}R_h
\geq
\left(D-\left\lceil\frac{\sqrt N}{3}\right\rceil\right)
\left\lfloor\sqrt{\frac N2}\right\rfloor
=\left(\frac1{\sqrt2}+o(1)\right)NL.
\tag{5.5}
\]
The matching upper bound \(O(NL)\) follows from \(R_h\leq\sqrt{M_h}\leq
\sqrt N\). \(\square\)

The bound in (5.4) is already at the conjectural random L2 scale. Thus even
an all-height version of [L2-FREQ] with a uniform ideal constant cannot
imply W1 or W3 through the row masses, row L2 moments, and degree bounds
alone. A certificate theorem valid only for \(h\leq h_0(p)\) contains still
less information. There is consequently no threshold growth rate for
\(h_0(p)\) that makes the proposed interpolation close: the case \(h_0=D\)
already fails.

The construction is deliberately a family of independent histograms, not an
Apéry orbit. Its conclusion is exactly a non-implication statement about
the proposed marginal statistics; it is not a counterexample to W1 for the
actual orbit.

## 6. The new residual [ZERO-TAIL-2]

The countermodel identifies the missing datum as concentration of the
distinguished value \(0\) across the lag variable, not further flatness
inside each individual value histogram.

For any nonnegative integer sequence \(R_1,\ldots,R_D\), define
\[
T_0=\left\lceil\sqrt{\frac ND}\right\rceil,
\qquad
A_D(t)=\#\{h\leq D:R_h\geq t\},
\tag{6.1}
\]
and
\[
\mathcal Z_D=\max_{t>T_0}t^2A_D(t),
\tag{6.2}
\]
with an empty maximum equal to zero.

### Proposition 6.1 (exact tail bridge)

For every \(N,D\),
\[
S_D=\sum_{h\leq D}R_h
\leq DT_0+\frac{\mathcal Z_D}{T_0}.
\tag{6.3}
\]
Consequently, at \(D=\sqrt N L\), with \(L\to\infty\) and
\(L=N^{o(1)}\), the residual
\[
\boxed{\ [\mathrm{ZERO\text{-}TAIL\text{-}2}]\qquad
\mathcal Z_D\ll N\ }
\tag{6.4}
\]
implies
\[
S_D\ll\sqrt{ND}+D
=N^{3/4}L^{1/2}+N^{1/2}L=o(N).
\tag{6.5}
\]
In particular, (6.4) implies W1 and hence the banked energy improvement.

#### Proof

The layer-cake identity gives
\[
S_D=\sum_{t\geq1}A_D(t).
\]
For \(t\leq T_0\), use \(A_D(t)\leq D\); for \(t>T_0\), definition (6.2)
gives \(A_D(t)\leq\mathcal Z_D/t^2\). Therefore
\[
S_D\leq DT_0+
\mathcal Z_D\sum_{t>T_0}\frac1{t^2}
\leq DT_0+\frac{\mathcal Z_D}{T_0},
\]
where the last step is the integral bound
\(\sum_{t>T_0}t^{-2}\leq\int_{T_0}^{\infty}x^{-2}\,dx\). Since
\(T_0\geq\sqrt{N/D}\), (6.5) follows. \(\square\)

This is not a relabeling of W1 or W3. It controls only the high tail above
the automatic threshold \(T_0\). For example, after discarding the first
\(O(T_0)\) lags to respect \(R_h\leq3(h-1)\), the model
\(R_h=T_0\) has \(\mathcal Z_D=0\), but
\[
S_D\asymp\sqrt{ND},\qquad \sum_{h\leq D}R_h^2\asymp N.
\]
At the mesoscopic scale its first quantity is much larger than
\(N^{2/3}\), and its second is much larger than the conjectural
[SAME-LAG-BDH] scale \(D N^\varepsilon\) for small fixed
\(\varepsilon\). Conversely, [SAME-LAG-BDH] implies (6.4), since
\[
t^2A_D(t)\leq\sum_{h\leq D}R_h^2.
\]
Thus [ZERO-TAIL-2] is a strictly relaxed, zero-fibre-only target. The
per-row L2 information gives only
\(\mathcal Z_D\leq\sum_h\operatorname{Col}(h)=O(ND)\), missing (6.4) by
the full factor \(D\); Theorem 5.1 shows that loss is real.

The subsequently banked residual [AVG-MOMENT-8] is also strictly stronger
than (6.4). Indeed, if
\(\sum_{h\leq D}\sum_a\nu_h(a)^8\ll ND\), then, for every \(t>T_0\),
\[
t^8A_D(t)\leq\sum_{h\leq D}R_h^8\ll ND,
\]
and hence
\[
\mathcal Z_D\ll\frac{ND}{T_0^6}
\leq\frac{D^4}{N^2}=L^4\ll N
\]
at \(D=\sqrt N L\), \(L=N^{o(1)}\). Thus [ZERO-TAIL-2] is the weaker
zero-value projection of that eighth-moment face; it does not ask for any
control at nonzero values.

## 7. Exact machine verification

The verifier is CRON_breakwall_verify.py, SHA-256
b7612ab83a9da7e652cc5eefdca72ac2b46bfdf60e2a59078791666ee293a621.
It performs the following exact gates.

1. It rebuilds \(N_3\) from the recurrence and checks the coefficient
   obstruction in (3.9).
2. It verifies every identity in (4.1)--(4.3), including exact cancellation,
   numerator squarefreeness, denominator \(q_8^3\), and reflection.
3. It verifies the symbolic identity (3.8).
4. It constructs the marginal extremizer row by row and checks mass, degree,
   and \(\operatorname{Col}(h)\leq2M_h\) exactly.
5. It extracts the function orbit(p) verbatim from
   CRON_b1_crosscorr.py by AST, computes the actual nonwrapping root counts,
   and checks (6.3) at all three required primes. The finite test takes
   \(L=\lceil\log\log p\rceil\) and the least integer \(D\) with
   \(D^2\geq NL^2\). All counts and inequalities are exact; displayed ratios
   alone use floating point.

Reproduction:

~~~sh
PYTHONPYCACHEPREFIX=/tmp/codex-breakwall-pycache \
  python3 -m py_compile CRON_breakwall_verify.py
PYTHONPYCACHEPREFIX=/tmp/codex-breakwall-pycache \
  python3 CRON_breakwall_verify.py
~~~

Output:

~~~text
SYMBOLIC_GATES PASS
N3_X4_MISMATCH 68535 != 1155*58
DEGREE4_PROUHET P1-P2=-16*X - 72
DEGREE4_COMPARISON deg_num=21 deg_den=24 squarefree_num=yes reflection_odd=yes
FINITE_FIELD_GATES
p N L D T0 S_D maxR Z_D Z_D/N coarse_bound model_S model_maxR model_Z_D model_maxCollision/mass
1009 1008 2 64 4 90 5 50 0.049603175 268.500000 1791 31 48600 1.967741935
  R_HIST {0: 20, 1: 19, 2: 10, 3: 11, 4: 2, 5: 2}
3001 3000 3 165 5 250 8 98 0.032666667 844.600000 8316 54 412923 1.981481481
  R_HIST {0: 50, 1: 48, 2: 27, 3: 26, 4: 5, 5: 7, 7: 1, 8: 1}
10007 10006 3 301 6 475 8 147 0.014691185 1830.500000 28020 99 2573872 1.989898990
  R_HIST {0: 94, 1: 79, 2: 45, 3: 55, 4: 8, 5: 15, 6: 2, 7: 2, 8: 1}
FINITE_FIELD_GATES PASS
ALL_BREAKWALL_GATES_PASS
~~~

The three actual orbits lie far inside [ZERO-TAIL-2] in these finite
windows: \(\mathcal Z_D/N\leq0.050\). This is evidence only. The abstract
extremizer simultaneously keeps every row below twice its mass and produces
large model_Z_D, exactly as the proof predicts.

## 8. Terminal status

- [NO-RIGHT-2-3]: **PROVED for every \(h\geq2\)**.
- Full indecomposability / W7: **OPEN for right degrees at least \(4\)**.
- Per-row L2 census alone: **DEAD as a route to W1/W3**, even if available
  uniformly for every \(h\leq D\).
- [ZERO-TAIL-2]: **new precise sufficient residual**, empirically supported
  at the three mandated primes, but **not proved asymptotically**.
- W1--W7 and a strict improvement over the \(3/2\) energy record:
  **NOT PROVED**.
