ANSWER Q6268 ddd07fdb

# Problem 3.2 — cold audit of Q6231 and the simple distinguished-state complement

## Bottom line

I get a **strict partial theorem, not the requested power saving**.

1. The max-agent separated-resultant frontier is source-compatible and, in one place, can be strengthened: for the actual affine root multiplicity
   \(m_{a,b}(p)=\#\{x\in\mathbf F_p:N_a(x)=N_b(x+a)=0\}\), Smith normal form of the integer Sylvester matrix gives
   \[
   m_{a,b}(p)\le v_p\!\left(\mathcal R_{a,b}\right),
   \qquad
   \mathcal R_{a,b}=\operatorname{Res}(N_a(X),N_b(X+a))\ne0.
   \]
   Therefore
   \[
   \sum_{X<p\le2X}\sum_{a+b\le H}m_{a,b}(p)
   \ll \frac{H^4\log(2H)}{\log X}.
   \]
   This repairs the full factor-\(H\) support-versus-root-multiplicity loss recorded in `meso_result.tex`; projective roots at infinity and nonsplit factors only consume extra resultant valuation and do not invalidate the one-way bound.

2. This does **not** thin to the distinguished Apéry state. On the safe range there is an exact projective resolvent \(q_x\):
   \[
   x\in R_p(a,b)
   \Longrightarrow
   \pi_p(x)=\pi_p(x+a)=\pi_p(x+a+b)=q_x.
   \]
   The actual Apéry triples are exactly the slice
   \[
   x\in Z_p\cap R_p(a,b)
   \quad\Longleftrightarrow\quad
   x\in R_p(a,b)\ \text{and}\ q_x=q_0,
   \]
   where \(q_0\) is the distinguished initial projective line. Simplicity/transversality of the common root says nothing presently proved about the equality \(q_x=q_0\).

3. The first exact failed implication in the simple, off-centre, nonreflection route is therefore
   \[
   N_a(x)=N_b(x+a)=0,
   \quad x\ \text{simple/transverse}
   \quad\Longrightarrow\quad q_x=q_0.
   \]
   The cited local identities prove only the existence of the unique repeated state \(q_x\). They do not identify it with the distinguished initial line. Q6190's exact local-return countermodels already prove that recurrence/continuant identities alone cannot make that identification; I do **not** use their finite frequency as asymptotic evidence.

4. The normalization in the proposed critical consumer has to be corrected. If \(K=Z(p)=N^{3/5+o(1)}\), \(H=N^{2/5+o(1)}\), and the sum counts **consecutive actual-zero starts**, then one start determines one consecutive gap pair. Its natural mass per rich critical prime is \(\asymp K\), not \(KH\). Consequently
   \[
   I(N,H)\ll H^{4-\varepsilon+o(1)}
   \]
   gives only
   \[
   \#\mathcal P_{3/5}\ll N^{1-2\varepsilon/5+o(1)}.
   \]
   To reach
   \[
   \#\mathcal P_{3/5}\ll N^{3/5-2\sigma},
   \]
   the unweighted consecutive-start incidence must satisfy
   \[
   I(N,H)\ll H^{3-5\sigma+o(1)},
   \]
   equivalently an \(H^{4-\varepsilon}\) theorem with
   \[
   \varepsilon\ge 1+5\sigma.
   \]
   A small unspecified \(\varepsilon>0\) is not enough.

5. Therefore the Q6231 singular-rich population exponent \(\le3/5\) is **not fully source-justified from the connector-visible source**. Its exponent transfer is correct if Q6231 really proved a **root-witness-weighted** singular estimate
   \[
   I_{\rm sing}(N,H)\ll H^{3+o(1)}
   \]
   and a positive proportion of the \(K\) consecutive starts are legal. But the exact Q6231 theorem lives in the caller-authoritative unpushed SHA, and that SHA is not resolvable by GitHub; the visible mesoscopic source explicitly warns that support counts cannot be substituted for root-weighted energy. I cannot certify an \(H^3\) singular witness theorem from stale remote text.

6. For the **simple distinguished-state complement**, the exact state-dependent inputs presently give no positive power \(H^{-\varepsilon}\). The strongest unconditional statement obtained here is the \(H^{4+o(1)}\) root-weighted bound above together with the exact selector \(q_x=q_0\). Thus the state-zero saving is
   \[
   \varepsilon_{\rm state}=0.
   \]
   The desired \(N^{3/5-2\sigma}\) population bound is **not reached**.

The rest of this report gives the theorem-level audit and proofs.

---

# 0. Source boundary

The caller-authoritative state is

```text
/Users/huangx/repos/Ramanujan_Challenge
HEAD 670c5b947d66fa0dacf51df1e112fe5f08433d37
```

The GitHub connector cannot resolve that commit. I therefore cannot honestly claim to have re-read a theorem that exists only in that local, unpushed object, including the literal Q6231 endpoint predicates if they were added there.

I did recheck the same-project source that is connector-visible, principally:

- `problems/3.2/meso_result.tex` on `q6190-census`, blob `661385d846a0326dc2274bbb13ef74686cedc2ae`;
- `problems/3.2/meso_resultants.md` on `q6190-census`, blob `a1b42eeec6a19ad0d1ce8b4beeaabe6c71303ebc`;
- `problems/3.2/projective_variance_reduction.tex` on `q6190-census`, blob `4ac4e34709be27db1ea7d8378442d45779fedca1`;
- `problems/3.2/full_cycle_bridge_reciprocity.tex` on `q6190-census`;
- `problems/3.2/gap_kernel_cartier_defect.tex` on `q6190-census`;
- `drops/Q6111-9c3e62dd.md`, blob `fe2a9669dadd409a1d2d38e908e2e9a47ab30e22`;
- `drops/Q6190-0ed08843.md`, blob `86c54b34c7dc78a11770037c08c6be1c1a8bda96`.

No authoritative theorem/proof file was edited. The only write is this requested delivery drop, after the audit.

---

# 1. Reconstruct the actual critical consumer first

Let the actual Apéry zeros in one characteristic be
\[
0\le z_1<z_2<\cdots<z_K<p,
\qquad K=Z(p).
\]
Write consecutive gaps
\[
g_i=z_{i+1}-z_i.
\]
A consecutive triple start \(z_i\) has the unique pair
\[
(a_i,b_i)=(g_i,g_{i+1}),
\qquad z_{i+2}-z_i=a_i+b_i.
\]
Thus the map
\[
i\longmapsto (z_i,a_i,b_i)
\]
is injective and carries **one witness per start**. Summing over legal \((a,b)\) does not create an additional factor \(H\).

There is also an exact participation lemma, so no probabilistic gap model is needed.

## Lemma 1.1 — deterministic legal-triple participation

For any \(H>0\),
\[
\#\{1\le i\le K-2:z_{i+2}-z_i\le H\}
\ge K-2-\frac{2p}{H}.
\]

### Proof

The two-gap spans satisfy
\[
\sum_{i=1}^{K-2}(z_{i+2}-z_i)
=\sum_{i=1}^{K-2}(g_i+g_{i+1})
\le 2(z_K-z_1)<2p.
\]
Hence fewer than \(2p/H\) spans can exceed \(H\). ∎

Take
\[
H=C\frac pK,
\qquad C>2\ \text{fixed}.
\]
Then
\[
\#\{\text{span-}H\text{ consecutive triple starts}\}
\ge \left(1-\frac2C\right)K-O(1).
\]
At the critical scale \(K=N^{3/5+o(1)}\), this is \(H=N^{2/5+o(1)}\). Removing the single reflection-centred adjacency and finitely many endpoint/centering exceptions costs only \(O(1)\) starts; a singular/simple dichotomy can then split the remaining \(\gg K\) starts.

This is the precise unnormalized mass that a prime-population argument may divide by: **\(K\)**.

If Q6231 has an additional RR occurrence weight or a genuinely independent legal-anchor variable whose total mass is \(\gg KH\), that factor must appear literally in the summed object and be proved there. It cannot be restored by normalizing selected mass or by multiplying the unweighted consecutive-start count by \(H\).

## Corollary 1.2 — exact exponent transfer

Suppose a subpopulation \(\mathcal P\) of critical primes has \(\gg K\) legal starts per prime and
\[
I_{\mathcal P}(N,H)\ll H^{\theta+o(1)}.
\]
Then
\[
\#\mathcal P
\ll \frac{H^{\theta+o(1)}}K
=N^{(2\theta-3)/5+o(1)}.
\]
In particular:

\[
\begin{array}{c|c}
I\text{ bound}&\#\mathcal P\text{ exponent}\\ \hline
H^{4+o(1)}&1\\
H^{4-\varepsilon+o(1)}&1-2\varepsilon/5\\
H^{3+o(1)}&3/5\\
H^{3-5\sigma+o(1)}&3/5-2\sigma.
\end{array}
\]

This is the normalization I use below.

---

# 2. Cold audit of the max-agent frontier

## 2.1 `DA4R_eta`

Status: **still unproved**.

Q6190 states it only as a conditional premise: four consecutive actual Apéry zeros in a fixed \(\eta\sqrt p\) span should force an adjacent reflection pair. It also records exact local-recurrence return carriers showing that such a conclusion is false for arbitrary solutions of the same second-order recurrence. Those finite carriers are useful only to refute a recurrence-only universal implication; I do not extrapolate their frequency or distribution.

If `DA4R_eta` were proved for the distinguished Apéry orbit, Q6190's deterministic covering argument gives \(Z(p)\ll_\eta p^{1/2}\), which would eliminate the \(K=p^{3/5+o(1)}\) critical regime outright. There is no such proof in the source I can verify.

## 2.2 Nonzero separated resultants and triple divisibility

This part rechecks.

The source defines
\[
N_{h+1}(x)=P(x+h)N_h(x)-(x+h)^6N_{h-1}(x),
\qquad N_0=0,\ N_1=1,
\]
with \(P(n)=(2n+1)(17n^2+17n+5)\), and
\[
\mathcal R_{d,r}=\operatorname{Res}_x(N_d(x),N_r(x+d)).
\]

`meso_result.tex` proves that every complex root of \(N_d\) lies in
\[
-d<\Re x<-1,
\]
whereas the roots of \(N_r(x+d)\) lie in
\[
-d-r<\Re x<-d-1.
\]
The strips are disjoint, hence
\[
\mathcal R_{d,r}>0
\]
for all \(d,r\ge2\). This is a uniform proof, not finite-resultant extrapolation.

If
\[
r<r+h<r+h+k<p
\]
are three actual zeros, the exact gap/bridge identity gives
\[
N_h(r)=0,
\qquad N_k(r+h)=0\pmod p,
\]
so
\[
p\mid \mathcal R_{h,k}.
\]
The source's Hadamard/resultant-height estimates give the stated scale
\[
\log \mathcal R_{h,k}\ll hk\log(2(h+k)).
\]
Therefore
\[
\log p\ll hk\log(2(h+k))
\]
for such a triple.

## 2.3 Three actual zeros in a \(p^\alpha\) window

This also rechecks, and follows from the stronger root-weighted estimate proved in §3 below.

Set \(H=(2X)^\alpha\). A prime \(X<p\le2X\) with three actual zeros in an interval of length \(p^\alpha\) contributes at least one actual affine root witness to some \((h,k)\) with \(h+k\le H\). Hence
\[
\#\{X<p\le2X:\text{such a triple exists}\}
\ll \frac{H^4\log(2H)}{\log X}
=X^{4\alpha+o(1)}.
\]
For \(\alpha<1/4\) this is a genuine sparse-prime bound. At the critical \(\alpha=2/5\) it is far above the number of primes and gives no useful population thinning.

## 2.4 Exact obstruction

The obstruction is exactly as stated: ordinary common local-return roots are state-blind. The rigorous form is §4 below.

---

# 3. New strict partial theorem: root multiplicity can be paid by resultant valuation

The visible source correctly warns that
\[
W_p(H)=\#\{(a,b):m_{a,b}(p)>0\}
\]
and
\[
\mathcal E_p(H)=\sum_{a+b\le H}m_{a,b}(p)
\]
are different. `meso_result.tex` therefore only took the crude route
\[
m_{a,b}(p)\ll\min(a,b)
\]
and lost a factor \(H\), producing an \(H^5\)-scale average from an \(H^4\)-scale support estimate.

For the actual affine root count, that loss can be removed by retaining \(p\)-adic valuation rather than only the indicator \(p\mid\mathcal R_{a,b}\).

## Theorem 3.1 — Sylvester valuation pays affine root multiplicity

Let
\[
A(X)=N_a(X),
\qquad B(X)=N_b(X+a),
\qquad R=\operatorname{Res}(A,B)\ne0.
\]
For a prime \(p\ge7\) with \(a,b<p\), let
\[
m_{a,b}(p)
=\#\{x\in\mathbf F_p:\bar A(x)=\bar B(x)=0\}.
\]
Then
\[
\boxed{m_{a,b}(p)\le v_p(R).}
\]

### Proof

The source gives \(\deg N_h=3(h-1)\). If the leading coefficient \(\ell_h\) vanishes modulo \(p\), its centered-coefficient/Cassini calculation shows that, for \(p\ge7\) and \(h<p\), the next relevant coefficient
\[
u_h=\frac{5h\ell_{h-1}}{256}
\]
is nonzero; the degree drops by exactly two. Thus \(\bar N_h\) is not the zero polynomial. Translation preserves nonzeroness, so \(\bar A,\bar B\ne0\).

Let \(S(A,B)\) be the nominal integer Sylvester matrix. Over \(\mathbf F_p\), every distinct common \(\mathbf F_p\)-root contributes a linear factor to
\[
\gcd(\bar A,\bar B),
\]
so
\[
m_{a,b}(p)
\le\deg\gcd(\bar A,\bar B)
\le\operatorname{corank}_{\mathbf F_p}S(A,B).
\]
The last inequality remains safe when nominal degrees drop; degree defects only add possible Sylvester kernel.

Since \(R=\pm\det S(A,B)\ne0\), put the integer Sylvester matrix in Smith normal form. If its reduction modulo \(p\) has corank \(c\), then at least \(c\) Smith invariant factors are divisible by \(p\). Hence
\[
c\le\sum_i v_p(s_i)=v_p(\det S)=v_p(R).
\]
Combining the inequalities proves the claim. ∎

### Why infinity/nonsplit factors do not break this theorem

A common projective point at infinity can increase the Sylvester corank and \(v_p(R)\), but it does not create an extra affine \(\mathbf F_p\)-root on the left. Likewise a nonsplit common factor consumes resultant valuation without being counted by \(m_{a,b}(p)\). Both phenomena make the inequality looser, never false.

## Corollary 3.2 — global root-weighted \(H^4\) bound

For \(2\le H<X\),
\[
\boxed{
\sum_{X<p\le2X}\mathcal E_p(H)
\ll
\frac{H^4\log(2H)}{\log X}.}
\]

Indeed,
\[
\begin{aligned}
(\log X)\sum_{X<p\le2X}\mathcal E_p(H)
&\le
\sum_{a+b\le H}\sum_{X<p\le2X}
  v_p(\mathcal R_{a,b})\log p\\
&\le
\sum_{a+b\le H}\log\mathcal R_{a,b}\\
&\ll
\sum_{a+b\le H}ab\log(2H)\\
&\ll H^4\log(2H).
\end{aligned}
\]

This is an actual improvement over the support-plus-crude-degree route recorded in the visible source. It is also exactly the right way to cold-audit root multiplicity: use valuations, not selected-mass normalization and not the indicator \(p\mid R\) alone.

For the distinguished state one immediately gets the strict partial theorem
\[
\sum_{X<p\le2X}\sum_{a+b\le H}
\#(Z_p\cap R_p(a,b))
\ll
\frac{H^4\log(2H)}{\log X},
\]
but this is only the trivial subset inequality at the state level. It contains **no** \(H^{-\varepsilon}\) distinguished-state saving.

---

# 4. Exact state resolvent: what a simple local return root really gives

Let \(\pi_p(n)\in\mathbf P^1(\mathbf F_p)\) denote the projective state whose recurrence solution vanishes at index \(n\). Equivalently, choose any two fundamental solutions \(B,C\); the evaluation row at \(n\) is nonzero on the safe range and its kernel is a unique projective line \(q_n\). I write
\[
q_n=\pi_p(n).
\]
Let \(q_0\) denote the line of the distinguished Apéry initial solution.

The exact gap criterion used in `projective_variance_reduction.tex` is the projective form of the continuant bridge:
\[
N_h(x)=0
\quad\Longleftrightarrow\quad
\pi_p(x+h)=\pi_p(x)
\]
on a nonwrapping safe segment.

Therefore:

## Theorem 4.1 — state-resolvent decomposition of a common return root

If
\[
0\le x<x+a<x+a+b<p,
\qquad
N_a(x)=N_b(x+a)=0\pmod p,
\]
then there is a unique projective state
\[
q_x=\pi_p(x)
\]
for which the corresponding recurrence solution vanishes at all three indices
\[
x,\quad x+a,\quad x+a+b.
\]
Moreover
\[
\boxed{
x\in Z_p\cap R_p(a,b)
\iff
x\in R_p(a,b)\ \text{and}\ q_x=q_0.\ }
\]

### Proof

The first return gives
\[
\pi_p(x+a)=\pi_p(x),
\]
and the second gives
\[
\pi_p(x+a+b)=\pi_p(x+a).
\]
Thus all three evaluation kernels are the same unique projective line \(q_x\). By definition the distinguished Apéry sequence vanishes at \(x\) exactly when its initial line is that kernel, namely \(q_x=q_0\). Once \(q_x=q_0\), the two return identities propagate the same distinguished zero to the two later indices. ∎

This theorem is the exact state-zero thinning formulation:
\[
R_p(a,b)=\bigsqcup_q R_p^{(q)}(a,b),
\qquad
Z_p\cap R_p(a,b)=R_p^{(q_0)}(a,b).
\]
Ordinary resultants control the union over \(q\). The missing theorem must control one fixed fiber.

## First exact failed implication

The simple branch would need something of the form
\[
\boxed{
\begin{array}{c}
N_a(x)=N_b(x+a)=0,\\
\text{simple/transverse},\ \text{off-centre},\ \text{nonreflection}
\end{array}
\Longrightarrow
q_x=q_0.}
\]
No such implication is present in the audited source. Simplicity says that the local root/intersection is transverse in the \(x\)-geometry. It does not select which one-dimensional projective kernel realizes the return word.

The exact conclusion available before inserting the distinguished initial condition is only
\[
\boxed{\text{there exists a unique repeated state }q_x.}
\]
This is the first place where an ordinary-resultant proof ceases to address the actual Apéry orbit.

---

# 5. Audit of the genuinely distinguished inputs

I checked each of the allowed mechanisms specifically for an additional equation or a power-saving count on the slice \(q_x=q_0\).

## 5.1 Initial projective line

This is necessary and exact: it identifies the target fiber as \(q_0\). But by itself it turns the problem into
\[
\#\{x\in R_p(a,b):\pi_p(x)=q_0\},
\]
which is exactly the unknown state concentration. There is no low-degree polynomial in \((x,a,b)\) supplied by this identification.

## 5.2 Full-period selector and Lucas

The vector Lucas law retains rank two; the relevant block acts by the scalar Apéry factor on a two-dimensional state. At an actual zero the scalar vanishes, but the law does not collapse the ordinary local-return locus onto the distinguished line. In particular, it does not convert
\[
q_x=q_0
\]
into a new \(H\)-degree equation independent of the Apéry zero predicate.

## 5.3 Exact first-Witt/Bockstein transport

Q6111 proves, on a fully actual safe segment, the exact endpoint formula
\[
b_{x+1}A_u(x)^{-1}\frac{N_u(x)}p
=\eta_{x+u}+\lambda_u(x)\eta_x\pmod p,
\qquad \lambda_u(x)\ne0,
\]
and in the literal RR chart
\[
\Theta_c
=\eta_{x+u}+\lambda_u(x)\eta_x+\kappa_cJ_u(x).
\]

This is real extra distinguished information, but on an off-centre open chain it is a **transport relation**, not a closed constraint. Given \(\eta_x\), it propagates the first-Witt coordinate to the next selected zero; there is no source theorem forcing cancellation or restricting the number of admissible \(q_0\)-roots by a power of \(H\). Q6190 already records the corresponding `DA4R_eta`/first-Witt rigidity as unproved.

## 5.4 Reflection and full-cycle reciprocity

`full_cycle_bridge_reciprocity.tex` gives, on a return,
\[
v_{x+h}=\mu_h(x)v_x,
\qquad \mu_h(x)\ne0,
\]
and the reflection law gives the reciprocal scalar relation at the reflected segment. These identities become a closed product only when the geometry actually meets the reflected/full-period segment. The present seam explicitly removes reflection-centred cases. In the nonreflection simple complement they provide no independent state selector.

## 5.5 Cartier/Mellin defect

`gap_kernel_cartier_defect.tex` proves an exact boundary defect whose coefficient is the Apéry value/Hasse section itself; schematically the Cartier extraction recovers
\[
-a^3b_{a-1}z.
\]
The bivariate version stores the full Hasse section
\[
H_p(T)=\sum_{r=0}^{p-1}b_rT^r.
\]
Thus Cartier/Mellin does remember the distinguished orbit, but using that defect as an \(H\)-degree selector for \(q_x=0\) is circular: the coefficient being tested is precisely the Apéry-zero data one is trying to thin. No independent bounded-degree eliminant emerges.

## 5.6 Exact state resolvent

The resolvent \(q_x=\pi_p(x)\) gives the cleanest exact characterization, but its distinguished fiber is exactly \(Z_p\). Turning it into an integer resultant against \(q_0\) requires the global initial-to-\(x\) transfer (or the degree-\(p\) Hasse section), losing the short \(H\)-scale height on which the separated-resultant method depends.

## 5.7 Projective variance

This is the one source theorem that genuinely tries to pass from all states to one fixed state.

For four-consecutive-occurrence windows of span \(\le H\), `projective_variance_reduction.tex` proves
\[
\sum_q C_p(q)^2\le30pH.
\]
It also proves that a separated-energy bound
\[
\sum_{X<p\le2X}E_p^{\rm sep}
\ll \frac{X^{2+o(1)}}{\log X}
\]
would imply for a fixed distinguished state
\[
\sum_{X<p\le2X}C_p(q_0)
\ll \frac{X^{3/2+o(1)}}{\log X}.
\]

Two limitations prevent this from closing the present seam.

First, the needed separated-energy hypothesis is not proved. The unconditional short-chain estimate at \(H=X^{2/5+o(1)}\) gives an aggregate variance scale \(X^{12/5+o(1)}/\log X\), which after Cauchy is no better than about \(X^{17/10+o(1)}/\log X\) for a fixed fiber.

Second, even the **conditional** \(X^{3/2}\) fixed-fiber estimate is a four-window statement. Relative to \(H=X^{2/5}\), it is \(H^{15/4+o(1)}\), only an \(H^{-1/4}\) improvement from \(H^4\), and it does not directly bound every consecutive triple start. If one could transfer it without loss to the triple consumer, division by \(K=X^{3/5}\) would still give exponent \(9/10\), not \(3/5\).

So projective variance is a genuine distinguished route, but the current theorem is quantitatively short of the required \(H^{-1-5\sigma}\) saving.

---

# 6. Q6231 singular-rich branch: exact audit conditions

Because the actual Q6231 theorem text at SHA `670c5b...` is not remotely resolvable, I cannot stamp the claimed \(\le3/5\) exponent as fully source-verified. The exact conditions under which that exponent is valid are nevertheless clear.

Let \(I_{\rm sing}(N,H)\) count **singular root witnesses** \((p,a,b,x)\), not just supported pairs \((p,a,b)\). Suppose:

1. every prime in the singular-rich critical subpopulation has
   \[
   I_{{\rm sing},p}\gg K;
   \]
2. one has the global root-witness estimate
   \[
   I_{\rm sing}(N,H)\ll H^{3+o(1)}.
   \]

Then, and only using the correct denominator \(K\),
\[
\#\mathcal P_{\rm sing}
\ll \frac{H^{3+o(1)}}K
=N^{3/5+o(1)}.
\]
So the **exponent transfer itself is correct**.

What must not be substituted is a theorem only about pair support,
\[
\#\{(p,a,b):\text{a singular eliminant vanishes}\}\ll H^{3+o(1)},
\]
unless the proof also pays the number of distinct starts \(x\) (or the relevant intersection multiplicity). The visible source contains an explicit audit of exactly this mistake for ordinary resultants: \(W_p\), \(V_p\), and \(\mathcal E_p\) are different objects, and a pair can have multiple root witnesses.

The Smith-normal-form theorem in §3 repairs this multiplicity problem for the ordinary resultant at the \(H^4\) scale. It does **not** by itself manufacture an \(H^3\) singular estimate. To certify Q6231's \(3/5\) singular-rich claim one must re-read its actual singular eliminant and verify an analogue such as
\[
m^{\rm sing}_{a,b}(p)\le v_p(D_{a,b}),
\qquad
\sum_{a+b\le H}\log|D_{a,b}|\ll H^{3+o(1)}\log X,
\]
or another genuinely root-weighted \(H^3\) theorem. That object is not present in the source I can resolve.

**Audit verdict:** Q6231's \(3/5\) singular-rich exponent is **conditionally arithmetically consistent but not fully source-certified here**. The two things that must be checked in the inaccessible local theorem are (a) root-witness multiplicity, and (b) the actual \(H^3\), rather than support-only, height/count.

---

# 7. What the simple distinguished-state complement would actually need

Define the literal unweighted critical simple incidence
\[
I_{0,\rm simp}(N,H)
=
\sum_{\substack{N<p\le2N\\p\ \rm critical}}
\sum_{\substack{a,b\ge2\\a+b\le H\\\rm legal}}
\#\left\{
\begin{array}{l}
x:\ x,x+a,x+a+b\text{ are consecutive actual zeros},\\
x\in R_p(a,b),\ \text{common root simple/transverse},\\
\text{off-centre and nonreflection}
\end{array}
\right\}.
\]

For the simple-rich complement, Lemma 1.1 supplies \(I_{0,\rm simp,p}\gg K\) after the singular branch and the finitely many excluded geometric cases have been removed by the Q6231 dichotomy.

The quantitatively correct target is therefore
\[
\boxed{
I_{0,\rm simp}(N,H)
\ll H^{3-5\sigma+o(1)}.}
\]
This implies
\[
\#\mathcal P_{3/5,\rm simp}
\ll
\frac{H^{3-5\sigma+o(1)}}K
=N^{3/5-2\sigma+o(1)}.
\]

Writing the goal instead as
\[
I_{0,\rm simp}\ll H^{4-\varepsilon+o(1)}
\]
is equivalent only if one remembers the threshold
\[
\boxed{\varepsilon\ge1+5\sigma.}
\]
A generic small \(\varepsilon\) does not close the RR consumer.

The source-backed theorem from §3 is only
\[
I_{0,\rm simp}\le I_{\rm all}
\ll H^{4+o(1)}.
\]
Thus the exact simple-state gain currently proved is
\[
\boxed{\varepsilon_{\rm state}=0.}
\]

---

# 8. Required-output checklist

## (i) Theorem-level audit of Q6231

- Critical legal-start participation can be made deterministic: with \(H=Cp/K\), \(C>2\), there are \(\gg K\) span-\(H\) consecutive triple starts.
- Each consecutive start determines one \((a,b)\); there is no automatic factor \(H\).
- Pair support and root-witness energy are different. The visible source explicitly records this gap.
- Ordinary root multiplicity can be repaired by the new exact inequality
  \[
  m_{a,b}(p)\le v_p(\mathcal R_{a,b}),
  \]
  yielding global \(H^{4+o(1)}\) root-weighted energy.
- A singular-rich \(H^3\) witness theorem would indeed imply population exponent \(3/5\), but I cannot verify that theorem in the inaccessible local Q6231 source. Therefore the claim is not fully source-certified.

## (ii) Proof / counterfamily / strict partial theorem for state-zero thinning

Strict partial theorem proved:
\[
Z_p\cap R_p(a,b)
=\{x\in R_p(a,b):q_x=q_0\},
\]
with unique projective resolvent \(q_x\), and
\[
\sum_{X<p\le2X}\sum_{a+b\le H}
\#(Z_p\cap R_p(a,b))
\ll \frac{H^4\log(2H)}{\log X}.
\]
This is exact but has no state power saving.

The source also contains exact recurrence-only local-return countermodels, so local continuant identities alone cannot imply distinguishedness. I use those only as logical counterexamples, not as finite-data extrapolation.

## (iii) First exact failed implication

\[
N_a(x)=N_b(x+a)=0
+\text{ simple/off-centre/nonreflection}
\not\Rightarrow_{\rm proved}
q_x=q_0.
\]
The exact proved conclusion is only a unique repeated state \(q_x\).

## (iv) Exact population exponent gained; does it reach the target?

- New multiplicity repair: improves the previously audited ordinary root-weighted global estimate from the crude \(H^5\) route to \(H^4\). This is a full factor \(H\) improvement in the **ordinary-root accounting**.
- Distinguished-state saving beyond that: \(\varepsilon_{\rm state}=0\).
- For the literal unweighted consecutive-start consumer, \(H^4/K=N^{1+o(1)}\), so there is no nontrivial critical-prime population exponent from this alone.
- A verified singular witness bound \(H^3\) would give \(N^{3/5+o(1)}\).
- The requested \(N^{3/5-2\sigma}\) requires \(H^{3-5\sigma+o(1)}\), or \(\varepsilon\ge1+5\sigma\) in an \(H^{4-\varepsilon}\) formulation.
- **Target not reached.**

## (v) Discipline

- No finite census is extrapolated into an asymptotic theorem.
- No selected/union mass is normalized into a missing factor \(H\).
- No authoritative source was edited before or during the audit.
- The connector-visible repository cannot resolve the caller-authoritative SHA `670c5b947d66fa0dacf51df1e112fe5f08433d37`; statements unique to that local commit remain explicitly uncertified.

---

# 9. Recommended next theorem on this seam

The next useful theorem is not another ordinary resultant. It is a genuinely fixed-fiber estimate at the correct normalization:

\[
\boxed{
\sum_{N<p\le2N}
\sum_{a+b\le H}
\#\{x\in R_p(a,b):q_x=q_0,\ x\text{ simple/off-centre/nonreflection and consecutive}\}
\ll H^{3-\delta}}
\]
for some \(\delta>0\), with \(H\asymp N^{2/5}\). To feed the stated \(2\sigma\) RR target one needs \(\delta\ge5\sigma\).

Equivalently, one needs a fixed-state anti-concentration theorem gaining a full factor \(H\) over the ordinary root-weighted \(H^4\) bound, plus the desired \(H^{-5\sigma}\) margin. Among the audited mechanisms, the only one already formulated as a fixed-fiber theorem is projective variance; its missing separated-energy input and its four-window-to-triple loss are therefore the sharpest source-faithful subproblems. First-Witt, reflection, Lucas, and Cartier currently transport or re-encode the distinguished state but do not provide that anti-concentration.
