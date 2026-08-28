ANSWER Q5224 ec39aebe

# Problem 3.2 — far–far active-pair attack after Q5217

## Verdict

I do **not** close the complete far–far active-degree second moment

\[
\sum_x w_xD_{x,T}^2
\]

at any exponent \(2\rho<145/93\) from the currently proved Problem 3.2 interfaces.  The two activity witnesses do, however, produce one new exact determinant and one genuinely useful unconditional deletion.

The main conclusions are:

1. **Exact two-witness algebra.**  If the two active lower copies are \(\ell_1,\ell_2\) and chosen live targets are \(\lambda_1,\lambda_2\) at oriented \(p\)-gaps \(h_1,h_2\in[T,2T)\), then each activity edge satisfies

   \[
   \varepsilon_i h_i\ell_i\lambda_i
   \equiv A_x(\ell_i-\lambda_i)\pmod p,
   \qquad \varepsilon_i\in\{\pm1\},
   \]

   and hence

   \[
   \boxed{
   p\mid
   \Phi_{12}:=
   \varepsilon_1h_1\ell_1\lambda_1(\ell_2-\lambda_2)
   -\varepsilon_2h_2\ell_2\lambda_2(\ell_1-\lambda_1).}
   \tag{A}
   \]

   In the mesoscopic range left after Q5217, if the four physical moving primes are pairwise distinct, then \(\Phi_{12}\ne0\).  For \(h_1=h_2\) and the same orientation, \(\Phi_{12}/h_1\) is, after relabelling, exactly the reciprocal four-prime cubic already isolated in Q5155.

2. **The hoped-for second CRT stage fails exactly.**  Writing the canonical full-\(P=pq\) completions as \(m_i,m_i^+\), activity controls only their \(p\)-difference:

   \[
   \boxed{m_i^+-m_i=\varepsilon_i h_i+pJ_i,\qquad J_i\in\mathbb Z.}
   \tag{B}
   \]

   If \(e_i=k_i^+-k_i\) is the reciprocal quotient displacement of the activity edge, then

   \[
   \boxed{
   \frac{(\lambda_i-\ell_i)m_i+arepsilon_i\lambda_i h_i}{p}
   +\lambda_iJ_i
   =q e_i.}
   \tag{C}
   \]

   The unselected copy has **no \(q\)-zero test**, so the integer \(J_i\) is not constrained by the activity condition.  Equation (C) is the first exact failed implication: the small actual \(p\)-return gap does not lift to a small/full-\(P\) row displacement, and (A) cannot be promoted to a second independent \(q\)-divisibility statement.

3. **New unconditional far–far deletion.**  Choose one canonical live activity witness for each active copy.  For an ordered pair of active lower copies, use the witness of the left lower label and let \(b\) be the nonwrapping distance from its target \(p\)-label to the other lower \(p\)-label.  For every dyadic \(B\ge2\), the literal pair mass with \(B\le b<2B\) satisfies

   \[
   \boxed{
   \sum_x w_x\,#\{\text{such ordered pairs}\}
   \ll N^{o(1)}T B\min(B,T)M_0.}
   \tag{D}
   \]

   This uses only one activity witness and therefore pays **no** \(T^{2/3}\) witness multiplicity.  The three actual \(p\)-zeros form a genuine \(m_{a,b}(p)\) witness, and the canonical nonvanishing/degree theorem gives

   \[
   m_{a,b}(p)\le3(\min(a,b)-1).
   \]

   No \(H\le\sqrt p\) collision-energy hypothesis is needed.

   Consequently, for \(2\le K\le T\), the cross-near sector \(b\le K\) is

   \[
   \boxed{
   \ll N^{o(1)}T K^2M_0.}
   \tag{E}
   \]

   The endpoint case \(b=0\) contributes only \(N^{2/3+o(1)}M_0\); on the Q4955 far fixed-\(p\) sector it forces \(T\gg N^{3/5}\), hence it is

   \[
   \ll T^{10/9+o(1)}M_0.
   \tag{F}
   \]

4. **Gate exponent.**  Taking

   \[
   K=T^{26/93-\zeta}
   \]

   gives

   \[
   TK^2=T^{145/93-2\zeta}.
   \]

   Thus the entire far–far sub-sector in which one canonical activity target lies within \(T^{26/93-\zeta}\) of the other lower \(p\)-label is already below the canonical AD2 gate.  For the ideal line, the corresponding proved cutoff is

   \[
   K=T^{1/10-\zeta},
   \qquad
   TK^2=T^{6/5-2\zeta}.
   \]

5. **What remains.**  After applying Q4955 and (E), the genuinely new obstruction may be restricted to far–far active lower pairs satisfying, in both directions,

   \[
   \boxed{
   |d|>N^{8/15},\qquad
   g_p^\diamond>N^{3/5},\qquad
   b_{1\to2}>T^{26/93-o(1)},\qquad
   b_{2\to1}>T^{26/93-o(1)}.}
   \tag{G}
   \]

   In this residual, either the two canonical activity targets coalesce, or the four physical moving primes are distinct and carry the nonzero determinant (A).  No current weighted-Palm theorem exploits either alternative with a power saving.

6. **Phase attack also stops at an exact place.**  Expanding all activity witnesses rather than choosing one incurs up to \(T^{2/3+o(1)}\) witnesses per active copy and hence potentially \(T^{4/3+o(1)}\) pair multiplicity.  Choosing one witness avoids that loss but is nonlinear and provides no Fourier averaging variable.  Moreover Q4703's exact Fourier expansion shows that at every successful moving-zero test one Artin–Schreier summand is the trivial sheaf.  Two successful activity targets therefore produce two positive trivial summands, not a new orthogonality mechanism.

So the strongest bankable result in this round is the one-witness cross-gap theorem (D)--(E).  The complete far–far AD2 estimate remains open.

---

# 0. Source boundary

Only `xiangyazi24/Ramanujan_Challenge`, Problem 3.2, is used.

The connected canonical mathematical head remains

```text
c5d932b66ce5e4f1657b587215d290ae7a13018b
```

and I re-read the relevant canonical files through the GitHub connector:

- [`problems/3.2/proof.tex`](https://github.com/xiangyazi24/Ramanujan_Challenge/blob/c5d932b66ce5e4f1657b587215d290ae7a13018b/problems/3.2/proof.tex): `lem:no-consec`, `lem:gap-poly`, `lem:nonvanish`, the restart/subinterval/column facts;
- [`problems/3.2/energy_result.tex`](https://github.com/xiangyazi24/Ramanujan_Challenge/blob/c5d932b66ce5e4f1657b587215d290ae7a13018b/problems/3.2/energy_result.tex):
  \[
  m_{a,b}(p)=\#\{u:N_a(u)=N_b(u+a)=0\},
  \]
  exact column/energy identities and the structured shallow-count argument;
- [`problems/3.2/pairpalm_result.tex`](https://github.com/xiangyazi24/Ramanujan_Challenge/blob/c5d932b66ce5e4f1657b587215d290ae7a13018b/problems/3.2/pairpalm_result.tex): exact pair--Palm factorial hierarchy and the warning that higher Palm extension moments are genuinely new arithmetic;
- [`drops/Q4760-47227b62.md`](https://github.com/xiangyazi24/Ramanujan_Challenge/blob/chatgpt-drop/drops/Q4760-47227b62.md): determinant-counting reciprocal-fibre second moment and its explicit counting-measure/Palm-measure mismatch.

The current local AT continuation is newer than connected `main`.  For that layer I use the exact same-project continuation interfaces already cold-audited in:

- Q4955 `bf086dde`: literal two-moving-copy Palm variables, short reciprocal-displacement and short fixed-\(p\) diamond-gap closures, far--far residual;
- Q5199 `7c002b92`: source-faithful active-degree definition and pointwise compression;
- Q5217 `e627891c`: reduction to AD2 and the exact required exponents;
- Q5155 `f9bf535f`: reciprocal four-prime cubic and the “no second reciprocal divisor” audit;
- Q4703 `e53d1b19`: exact one-moving-selector Fourier/trace normal form and trivial-summand obstruction.

The Q5199 and Q5217 artifacts are current continuation records; no claim is made that their unpushed labels occur in stale GitHub `main`.

No Python, sandbox, finite computation, or unrelated project is used.  No canonical source file is edited.

---

# 1. Literal two-copy Palm variables

Fix one selected Palm root

\[
x=(\mathbf b,\delta,\ell_0),
\qquad
\mathbf b=(\sigma,p,q,p'),
\qquad
P=pq,
\qquad
\Delta_x=\sigma\delta.
\tag{1.1}
\]

Retain the complete nonnegative selected-root atom

\[
\boxed{
w_x=
\Lambda_{\mathbf b}|u_{\ell_0}|^2
 a_{\mathbf b}(\delta,\ell_0)
 J_{\mathbf b}(\delta,\ell_0),
\qquad
M_0=\sum_xw_x.}
\tag{1.2}
\]

Every four-test selected-root condition, source mask, phase-support restriction, branch/orientation mask, shell restriction, high-centred condition, near/reflected-near deletion and four-distinct root condition stays inside \(w_x\).  Nothing below adds a second outer weight.

Put

\[
C_x=\left\langle\Delta_x\overline{p'}\right\rangle_P,
\qquad
 a_x=\frac{p'C_x-\Delta_x}{P}.
\tag{1.3}
\]

For a literal unselected moving shell prime \(\ell\), define the Q4955/Q4904 completed coordinates

\[
 m_x(\ell)
 =\left\langle\Delta_x\overline{p'\ell}\right\rangle_P,
\tag{1.4}
\]

\[
 k_x(\ell)
 =\frac{\ell m_x(\ell)-C_x}{P},
\qquad
 n_x(\ell)=a_x+p'k_x(\ell).
\tag{1.5}
\]

Thus

\[
\boxed{
\Delta_x=p'\ell m_x(\ell)-P n_x(\ell).}
\tag{1.6}
\]

The two literal unselected actual-zero tests are

\[
 r_x(\ell):=m_x(\ell)\bmod p\in\mathcal Z_p,
\tag{1.7}
\]

\[
 s_x(\ell):=n_x(\ell)\bmod\ell
 =\left\langle-\Delta_x\overline P\right\rangle_\ell
 \in\mathcal Z_\ell.
\tag{1.8}
\]

The unselected copy has **no new \(q\)-zero selector**.  Its \(p'\)-reduction is inherited algebraically:

\[
 n_x(\ell)\bmod p'
 =-\Delta_x\overline P\pmod{p'},
\tag{1.9}
\]

which is the one already-selected \(p'\)-zero residue.  This is not a second \(p'\) selector.

Let \(\chi_x(\ell)\in\{0,1\}\) contain the literal moving shell, exclusions, one-copy source/bulk/orientation masks and the tests (1.7)--(1.8).  Every active copy below has \(\chi_x(\ell)=1\).

Now take two lower moving copies \(\ell_1,\ell_2\).  Write

\[
m_i=m_x(\ell_i),\quad
k_i=k_x(\ell_i),\quad
n_i=n_x(\ell_i),\quad
r_i=r_x(\ell_i),\quad
s_i=s_x(\ell_i).
\tag{1.10}
\]

The exact reciprocal quotient displacement is

\[
\boxed{
d:=k_2-k_1,}
\tag{1.11}
\]

and subtraction of (1.5) gives

\[
\boxed{
\ell_2m_2-\ell_1m_1=Pd,}
\tag{1.12}
\]

\[
\boxed{n_2-n_1=p'd.}
\tag{1.13}
\]

Q4955's remaining lower-pair sector imposes

\[
\boxed{|d|>N^{8/15}}
\tag{1.14}
\]

and

\[
\boxed{
g_p^\diamond(r_1,r_2)>N^{3/5},}
\tag{1.15}
\]

where

\[
g_p^\diamond(r_1,r_2)
=\min\bigl(\|r_2-r_1\|_p,
           \|r_2+r_1+1\|_p\bigr).
\tag{1.16}
\]

All pair masks in (1.14)--(1.16) are retained in the exact far--far AD2 residual; when proving an upper bound below I sometimes drop them, which only enlarges a nonnegative count.

---

# 2. Activity witnesses and the multiplicity issue

Fix one of the finite source orientations.  It is convenient to encode ordinary versus reflected motion by

\[
\varepsilon\in\{+1,-1\}.
\]

For an active lower copy \((\varepsilon,r,\ell)\), activity means that there exists a live target prime \(\lambda\), with the same one-copy source package, whose underlying actual \(p\)-label \(r^+\) satisfies

\[
\boxed{
r^+\equiv r+\varepsilon h\pmod p,
\qquad T\le h<2T,}
\tag{2.1}
\]

with the corresponding oriented representatives nonwrapping.  Both

\[
r\in\mathcal Z_p,
\qquad
r^+\in\mathcal Z_p
\tag{2.2}
\]

are actual distinguished-characteristic zeros, and the target's moving-characteristic test

\[
s_x(\lambda)\in\mathcal Z_\lambda
\tag{2.3}
\]

is also literal.

For the two lower copies choose witnesses

\[
(\lambda_i,h_i),
\qquad T\le h_i<2T,
\qquad i=1,2.
\tag{2.4}
\]

Write target completed coordinates

\[
m_i^+=m_x(\lambda_i),\quad
k_i^+=k_x(\lambda_i),\quad
n_i^+=n_x(\lambda_i),\quad
r_i^+=r_x(\lambda_i),\quad
s_i^+=s_x(\lambda_i).
\tag{2.5}
\]

Then

\[
r_i^+\equiv r_i+\varepsilon_i h_i\pmod p,
\qquad
r_i^+\in\mathcal Z_p,
\qquad
s_i^+\in\mathcal Z_{\lambda_i}.
\tag{2.6}
\]

Define the activity-edge reciprocal quotient displacement

\[
\boxed{e_i:=k_i^+-k_i.}
\tag{2.7}
\]

Exactly as in (1.12)--(1.13),

\[
\boxed{
\lambda_i m_i^+-\ell_i m_i=P e_i,}
\tag{2.8}
\]

\[
\boxed{n_i^+-n_i=p'e_i.}
\tag{2.9}
\]

The target-pair quotient displacement is therefore

\[
\boxed{
d^+:=k_2^+-k_1^+=d+e_2-e_1.}
\tag{2.10}
\]

No shortness of \(e_i\) or \(d^+\) follows from activity.

## 2.1 A witness is not unique

For one active lower label \(r\), the possible target \(p\)-labels lie in a consecutive interval of length \(O(T)\).  `cor:subinterval` gives

\[
\#\{r^+\in\mathcal Z_p:T\le |r^+-r|<2T\}
\ll T^{2/3+o(1)}.
\tag{2.11}
\]

Each target \(p\)-label has \(O(1)\) moving-prime lifts by the exact reciprocal fibre.  Hence

\[
\boxed{
\#\{\text{literal activity witnesses of one active copy}\}
\ll T^{2/3+o(1)}.}
\tag{2.12}
\]

Thus blindly attaching witnesses to both copies may cost

\[
T^{4/3+o(1)}.
\tag{2.13}
\]

This is not a hypothetical bookkeeping issue.  A fixed target label also has up to \(T^{2/3+o(1)}\) actual predecessor labels in the backward interval, so reversing a target-pair count back to lower pairs can lose the same factor per coordinate.

## 2.2 Canonical witness selection

For every active oriented copy, order its finite witness set first by \(h\), then by the target prime, and choose the least witness.  Denote it

\[
W^\star(\varepsilon,r,\ell)=(h^\star,\lambda^\star).
\tag{2.14}
\]

This is only a deterministic proof device.  It does not alter \(D_{x,T}\), does not remove an actual-zero predicate, and introduces no new weight.  It lets us use one actual witness **without** paying (2.12).

The price is analytic: a nonlinear least-witness selector is not a summation variable on which the existing Fourier identities provide orthogonality.  This tradeoff is important in Attack III below.

---

# 3. Attack I — eliminate the two activity targets

## 3.1 One activity edge

Reduce the exact reciprocal identity

\[
\ell_i m_i\equiv C_x\pmod P,
\qquad
\lambda_i m_i^+\equiv C_x\pmod P
\tag{3.1}
\]

modulo \(p\).  Put

\[
A_x:=C_x\bmod p.
\tag{3.2}
\]

Then

\[
\ell_i r_i\equiv A_x\pmod p,
\qquad
\lambda_i r_i^+\equiv A_x\pmod p.
\tag{3.3}
\]

Using (2.6),

\[
\lambda_i(r_i+\varepsilon_i h_i)
\equiv\ell_i r_i\pmod p.
\]

Since \(r_i\equiv A_x\ell_i^{-1}\pmod p\), this is exactly

\[
\boxed{
\varepsilon_i h_i\ell_i\lambda_i
\equiv A_x(\ell_i-\lambda_i)\pmod p.}
\tag{3.4}
\]

Equivalently, for fixed \((\ell_i,h_i)\), the target prime lies in the familiar Möbius residue class

\[
\boxed{
\lambda_i
\equiv
\frac{A_x\ell_i}{A_x+\varepsilon_i h_i\ell_i}
\pmod p,}
\tag{3.5}
\]

whenever the denominator is nonzero.  This is exactly the one-stage \(p\)-restriction already implicit in the reciprocal lift.

## 3.2 Four-prime determinant

Eliminate \(A_x\) between (3.4) for \(i=1,2\).  Every literal two-witness occurrence satisfies

\[
\boxed{
 p\mid\Phi_{12},}
\tag{3.6}
\]

where

\[
\boxed{
\Phi_{12}
=
\varepsilon_1h_1\ell_1\lambda_1(\ell_2-\lambda_2)
-
\varepsilon_2h_2\ell_2\lambda_2(\ell_1-\lambda_1).}
\tag{3.7}
\]

This is a genuine characteristic-zero integer determinant.  It uses both activity witnesses and is unavailable for an arbitrary pair of moving hits with no return witnesses.

### Nonvanishing on the four-distinct mesoscopic sector

After Q5217's large-\(T\) deletion, the unresolved gate range has \(T=N^{\tau+o(1)}\) with \(\tau<1\).  Hence, in every fixed shell box and for large \(N\),

\[
2T<\min(\ell_1,\lambda_1,\ell_2,\lambda_2).
\tag{3.8}
\]

Assume the four physical moving primes are pairwise distinct.  If \(\Phi_{12}=0\), then

\[
\varepsilon_1h_1\ell_1\lambda_1(\ell_2-\lambda_2)
=
\varepsilon_2h_2\ell_2\lambda_2(\ell_1-\lambda_1).
\tag{3.9}
\]

The prime \(\ell_1\) cannot divide \(h_2\), \(\ell_2\), or \(\lambda_2\).  Therefore it must divide \(\ell_1-\lambda_1\), hence \(\ell_1\mid\lambda_1\), forcing \(\ell_1=\lambda_1\), contradiction.  Thus

\[
\boxed{
\Phi_{12}\ne0
\quad\text{on the four-distinct mesoscopic sector}.}
\tag{3.10}
\]

Moreover

\[
|\Phi_{12}|\ll T N^3.
\tag{3.11}
\]

Hence for fixed \((h_1,h_2,\ell_1,\lambda_1,\ell_2,\lambda_2)\), only \(O(1)\) shell primes can divide \(\Phi_{12}\).

That observation is **absolute**, not Palm-relative.  The selected mass can concentrate on one of those divisors, so (3.11) does not imply a factor saving relative to \(M_0\).

## 3.3 Equal witness gaps collapse to the known cubic

If the two orientations agree and

\[
h_1=h_2=h,
\]

then

\[
\Phi_{12}
=-h\,\mathcal C_4,
\]

with

\[
\mathcal C_4
=(\ell_1-\lambda_1)\ell_2\lambda_2
 -(\ell_2-\lambda_2)\ell_1\lambda_1.
\tag{3.12}
\]

After the relabelling

\[
(\ell_1,\lambda_1,\ell_2,\lambda_2)
\leftrightarrow
(\ell_1,\ell_2,\ell_3,\ell_4),
\]

this is exactly the reciprocal four-prime cubic audited in Q5155.  Q5155 already shows that the other equal-gap/equal-offset Plücker expression does not supply a second independent reciprocal divisor.

Thus the two-witness determinant is genuinely new only for unequal \(h_1,h_2\); at equal gap it reduces to the known rank-one Plücker carrier.

## 3.4 A second \(p\)-minor is still rank one

Let \(g\) be a chosen signed direct \(p\)-label gap from \(r_1\) to \(r_2\), and let

\[
g^+=g+\varepsilon_2h_2-\varepsilon_1h_1
\tag{3.13}
\]

be the corresponding target-label gap.  The same reciprocal calculation gives

\[
g\ell_1\ell_2
\equiv A_x(\ell_1-\ell_2)\pmod p,
\tag{3.14}
\]

\[
g^+\lambda_1\lambda_2
\equiv A_x(\lambda_1-\lambda_2)\pmod p.
\tag{3.15}
\]

Hence

\[
p\mid
\Psi_{12}:=
 g\ell_1\ell_2(\lambda_1-\lambda_2)
 -g^+\lambda_1\lambda_2(\ell_1-\ell_2).
\tag{3.16}
\]

This looks like a second determinant, but over \(\mathbb F_p\) both (3.7) and (3.16) are minors of the same rank-one reciprocal parametrization

\[
r=A_x/\ell.
\]

Their Plücker combinations recover the three-prime gap equations; they do not increase the CRT rank.  In the equal-gap case they collapse to the same cubic family described above.

So two \(p\)-minors are not two independent moduli.

---

# 4. The exact failed second-stage CRT implication

The full completed labels are not separated by the small physical gap \(h_i\).  Since (2.6) is only a statement modulo \(p\), there is an integer \(J_i\) such that

\[
\boxed{
m_i^+-m_i=\varepsilon_i h_i+pJ_i.}
\tag{4.1}
\]

Insert this into the exact full-\(P\) edge identity (2.8):

\[
(\lambda_i-\ell_i)m_i
+\varepsilon_i\lambda_i h_i
+p\lambda_iJ_i
=pq e_i.
\tag{4.2}
\]

The first two terms are divisible by \(p\) by (3.4).  Therefore

\[
B_i
:=
\frac{(\lambda_i-\ell_i)m_i
      +\varepsilon_i\lambda_i h_i}{p}
\in\mathbb Z
\tag{4.3}
\]

and

\[
\boxed{B_i+\lambda_iJ_i=q e_i.}
\tag{4.4}
\]

This is the exact \(q\)-stage ledger.

The critical point is that activity imposes no condition on

\[
m_i\bmod q
= C_x\ell_i^{-1}\pmod q,
\qquad
m_i^+\bmod q
= C_x\lambda_i^{-1}\pmod q.
\tag{4.5}
\]

Those are algebraic residues, not actual selected \(q\)-zeros.  Consequently no current theorem forces \(J_i=0\), makes \(J_i\) short, or correlates \(J_1,J_2\) with a power saving.

The exact false implication is therefore

\[
\boxed{
\begin{gathered}
r_i,r_i+\varepsilon_i h_i\in\mathcal Z_p,\\
s_i\in\mathcal Z_{\ell_i},\quad
s_i^+\in\mathcal Z_{\lambda_i},\\
T\le h_i<2T
\end{gathered}
\not\Longrightarrow
m_i^+-m_i=\varepsilon_i h_i
\text{ in }\mathbb Z/P\mathbb Z.}
\tag{4.6}
\]

Only the reduction modulo \(p\) is available.

This is why the two activity witnesses do **not** produce a two-stage CRT saving.

---

# 5. Why there is no legal six-zero resultant

It is useful to state the type ledger explicitly.

For each of the four moving copies

\[
\ell_1,\lambda_1,\ell_2,\lambda_2,
\]

activity supplies:

- an actual zero in the common characteristic \(p\):
  \[
  r_1,r_1^+,r_2,r_2^+\in\mathcal Z_p;
  \]
- one actual zero in that copy's own moving characteristic:
  \[
  s_1\in\mathcal Z_{\ell_1},\quad
  s_1^+\in\mathcal Z_{\lambda_1},\quad
  s_2\in\mathcal Z_{\ell_2},\quad
  s_2^+\in\mathcal Z_{\lambda_2};
  \]
- the **same** inherited actual \(p'\)-residue for every completed second row, by (1.9);
- no actual \(q\)-zero on the unselected copies.

Thus \(p\) is the only characteristic receiving several independent new physical zero positions.  The four moving-characteristic zeros live in four different fields, and the \(p'\)-condition is one repeated residue rather than four distinct returns.

Therefore a resultant formed by pretending that, for example, the \(\ell_1\)-zero and the \(\lambda_1\)-zero are roots in one common characteristic is illegal.  The activity witnesses create a four-zero \(p\)-configuration, not a six-zero same-characteristic configuration.

---

# 6. Attack II — one canonical witness gives a three-zero codegree theorem

The determinant route does not close the moment.  A different attack does give an unconditional deletion.

## 6.1 Define the canonical cross gap

Fix one selected root \(x\) and one of the finitely many source orientation charts.  For every active oriented lower copy choose the canonical witness (2.14).

Take two distinct active lower copies.  After exchanging them if necessary, put their oriented actual \(p\)-labels in the order

\[
0\le r<u\le p-1.
\tag{6.1}
\]

Let

\[
W^\star(r,\ell_1)=(h,\lambda),
\qquad
T\le h<2T,
\qquad
r+h\in\mathcal Z_p,
\tag{6.2}
\]

where in the chosen chart the activity edge is nonwrapping.

Define the **one-witness cross gap**

\[
\boxed{b:=|u-(r+h)|.}
\tag{6.3}
\]

The pair is interlaced if \(r<u<r+h\) and noninterlaced if \(r+h<u\).

The target witness is literal.  No target actual-zero condition has been removed, and no target is summed with multiplicity.

## 6.2 Root count for one dyadic cross gap

Assume first

\[
B\le b<2B,
\qquad B\ge2.
\tag{6.4}
\]

### Interlaced case

If

\[
r<u<r+h,
\]

put

\[
a=u-r=h-b.
\tag{6.5}
\]

Then

\[
r,\quad r+a=u,\quad r+a+b=r+h
\]

are three actual nonwrapping \(p\)-zeros.  Hence

\[
N_a(r)=0,
\qquad
N_b(r+a)=0
\pmod p.
\tag{6.6}
\]

For fixed \((a,b)\), the number of such \(r\) is

\[
m_{a,b}(p).
\]

By `lem:nonvanish`, both gap polynomials are nonzero modulo \(p\), with degrees \(3(a-1)\) and \(3(b-1)\).  Therefore

\[
\boxed{m_{a,b}(p)\le3(\min(a,b)-1)\ll\min(T,B).}
\tag{6.7}
\]

There are \(O(T)\) values of \(h\) and \(O(B)\) values of \(b\), so the number of possible label triples is

\[
\ll T B\min(T,B).
\tag{6.8}
\]

### Noninterlaced case

If

\[
r<r+h<u,
\]

then

\[
r,\quad r+h,\quad r+h+b=u
\]

are three actual zeros.  Now

\[
N_h(r)=0,
\qquad
N_b(r+h)=0,
\]

so the same degree argument gives

\[
m_{h,b}(p)\ll\min(T,B).
\tag{6.9}
\]

Again (6.8) follows.

## 6.3 Recovering literal moving copies

For each fixed oriented \(p\)-label, the reciprocal congruence places the moving prime in one residue class modulo \(p\), and the shell has length \(O(N)\) with \(p\asymp N\).  Hence each of

\[
r,\quad u,\quad r+h
\]

has only \(O(1)\) literal moving-prime lifts in a fixed orientation/source box.

Thus passing from the label triple back to

\[
(\ell_1,\ell_2,\lambda)
\]

costs only \(N^{o(1)}\).  Every moving-characteristic zero test, source mask, reciprocal condition, lower-pair far condition and additional cross-copy deletion is a \(0\)--\(1\) restriction and can only decrease the count.

This proves, root by root,

\[
\boxed{
\#\{(\ell_1,\ell_2):
B\le |u-(r+h^\star)|<2B\}
\ll N^{o(1)}T B\min(B,T).}
\tag{6.10}
\]

Multiplying by the unchanged \(w_x\) and summing gives

\[
\boxed{
\mathfrak D_{2,\times B}(T)
\ll N^{o(1)}T B\min(B,T)M_0.}
\tag{6.11}
\]

This is the promised dyadic one-witness cross-gap theorem.

### Why no \(T^{2/3}\) factor appears

The ordered lower pair and the **chosen** target determine the three-zero configuration.  We never sum over all targets of the first lower copy.  Conversely, after fixing the three \(p\)-labels, all three physical prime fibres are \(O(1)\).  The possible \(T^{2/3}\) alternate witnesses are irrelevant.

This is exactly the advantage of canonical witness selection.

---

# 7. Cumulative cross-near deletion and exponent ledger

Let

\[
2\le K\le T.
\]

Summing (6.7) directly over \(2\le b\le K\) gives

\[
\boxed{
\mathfrak D_{2,\times[2,K]}(T)
\ll N^{o(1)}T K^2M_0.}
\tag{7.1}
\]

There are two exceptional tiny cross gaps.

## 7.1 \(b=1\)

This would make the other lower actual \(p\)-zero consecutive to the chosen target actual \(p\)-zero.  `lem:no-consec` therefore makes this sector empty.

## 7.2 \(b=0\)

The other lower copy has the same oriented \(p\)-label as the chosen target.  A fixed label has \(O(1)\) live prime lifts, so root by root this contributes \(O(D_{x,T})\).  Hence

\[
\mathfrak D_{2,\times0}(T)
\ll N^{2/3+o(1)}M_0.
\tag{7.2}
\]

On the Q4955 far fixed-\(p\) sector, \(b=0\) implies that the lower-label gap is the activity gap, so

\[
N^{3/5}<g_p^\diamond(r_1,r_2)\le h<2T.
\]

Thus

\[
N^{2/3}\ll T^{10/9}.
\tag{7.3}
\]

Therefore

\[
\boxed{
\mathfrak D_{2,\times0}^{\rm ff}(T)
\ll T^{10/9+o(1)}M_0.}
\tag{7.4}
\]

This is below both consumers:

\[
\frac65-\frac{10}{9}=\frac4{45}>0,
\]

\[
\frac{145}{93}-\frac{10}{9}=\frac{125}{279}>0.
\]

## 7.3 Canonical \(38/31\) gate

Take

\[
\boxed{K=T^{26/93-\zeta}}
\qquad(\zeta>0).
\tag{7.5}
\]

Then

\[
TK^2
=T^{1+52/93-2\zeta}
=T^{145/93-2\zeta}.
\tag{7.6}
\]

Combining (7.1) and (7.4), the far--far ordered pair sector in which the canonical target of one lower copy lies within \(K\) of the other lower label satisfies

\[
\boxed{
\mathfrak D_{2,\rm ff}^{\times\text{-near}}
\ll
T^{145/93-2\zeta+o(1)}M_0
+T^{10/9+o(1)}M_0.}
\tag{7.7}
\]

The same bound holds for the union of the two directions \(1\to2\) and \(2\to1\), at a constant-factor cost.

Thus this whole sub-sector satisfies AD2 with

\[
\boxed{\rho=\frac{145}{186}-\zeta.}
\tag{7.8}
\]

## 7.4 Ideal \(6/5\) line

Taking instead

\[
\boxed{K=T^{1/10-\zeta}}
\tag{7.9}
\]

gives

\[
TK^2=T^{6/5-2\zeta},
\]

and therefore the same cross-near sector is below the ideal line with a fixed margin.

This is the strongest unconditional new deletion obtained in this round.

---

# 8. What two activity witnesses look like after the deletion

Fix the canonical witness of each lower copy.  After removing the two directional cross-near sectors at the gate scale, we may impose

\[
\boxed{
|r_2-r_1^+|>T^{26/93-o(1)},
\qquad
|r_1-r_2^+|>T^{26/93-o(1)}
}
\tag{8.1}
\]

in every nonwrapping orientation chart, with the corresponding reflected versions treated in the finite reflection alphabet.

Together with Q4955,

\[
|d|>N^{8/15},
\qquad
g_p^\diamond(r_1,r_2)>N^{3/5}.
\tag{8.2}
\]

This is a strictly smaller residual than Q5217's raw rich-degree event.

There are then two structural subcases.

## 8.1 Coalescent target

It is possible that

\[
\lambda_1=\lambda_2.
\tag{8.3}
\]

The two lower primes are still distinct.  In a common orientation, target equality forces

\[
r_1+\varepsilon_1h_1
\equiv
r_2+\varepsilon_2h_2
\pmod p.
\tag{8.4}
\]

This is a genuine V-shaped three-zero configuration.  It is **not** excluded by (8.1), because each lower-to-common-target distance is of order \(T\), much larger than the gate cross-near cutoff.

A common target can have as many as \(T^{2/3+o(1)}\) actual predecessors in the relevant interval, so the present subinterval theorem allows a nontrivial star.  The available large-window collision energy is not applicable in exactly the far range where \(|h_1-h_2|>N^{3/5}\) forces \(T>N^{3/5}\).

Thus target coalescence is one honest residual configuration.

## 8.2 Four distinct moving primes

If the targets do not coalesce, then after the \(b=0\) deletion the four physical moving primes are distinct in the generic orientation sector, and the nonzero determinant (3.7) is available:

\[
p\mid\Phi_{12}\ne0.
\tag{8.5}
\]

The exact q-lift defect (4.4) remains.  No current Palm-relative divisor theorem converts (8.5) into a power saving.

This is the other honest residual configuration.

---

# 9. Attack III — phase/large-sieve use of the activity targets

The activity witnesses also seem to offer new oscillatory variables.  There are two ways to expose them, and both fail for exact reasons.

## 9.1 Sum all witnesses

Replace the active indicator by the positive witness count:

\[
\mathbf1_{\{\alpha\text{ active}\}}
\le
\sum_{(h,\lambda)\in\mathscr W_T(\alpha)}1.
\tag{9.1}
\]

For two lower copies this exposes two target primes and two return gaps.  However (2.12) gives

\[
|\mathscr W_T(\alpha)|\ll T^{2/3+o(1)},
\]

so a two-target expansion can be larger by

\[
T^{4/3+o(1)}.
\tag{9.2}
\]

Reindexing by target pairs does not remove the loss: a fixed target has up to \(T^{2/3+o(1)}\) actual predecessors in the backward interval.

For example, feeding Q4955's \(N^{6/5+o(1)}M_0\) short-pair theorem through this inverse multiplicity would give at best

\[
N^{6/5+o(1)}T^{4/3}M_0,
\]

far above the AD2 gate.

## 9.2 Fourier-expand the target moving-zero tests

For a target characteristic \(\lambda\), the exact canonical identity is

\[
\mathbf1_{\mathcal Z_\lambda}(s_\lambda)
-\frac{Z(\lambda)}\lambda
=
\frac1\lambda
\sum_{a\ne0}F_\lambda(a)e_\lambda(-as_\lambda).
\tag{9.3}
\]

Q4703 audits the corresponding trace function.  If the target actually hits, the term indexed by the actual zero

\[
u=s_\lambda
\]

produces the constant Artin--Schreier trace.  Its complete frequency sum is \(\lambda-1\), not square-root sized.

Thus one successful target contributes a positive trivial summand.  Two successful activity targets contribute the product of two such positive trivial summands.  Local conductor bounds and independent finite-field Cauchy do not produce cancellation across changing characteristics.

Choosing one canonical target avoids (9.2), but then there is no linear target sum left on which to apply a large sieve.

Hence activity does not currently produce a phase-orthogonality theorem unavailable for arbitrary moving-hit pairs.

A genuinely new **horizontal** operator theorem, uniform under the selected-Palm conditioning and able to suppress the trivial-hit contribution after coupling different characteristics, would be required.

---

# 10. Sharp current-interface saturation after the new deletion

The following is an **interface saturation packet**, not an actual Apéry construction.  It is included to show that the presently proved degree, spacing, zero-count, reciprocal-fibre and Palm-pair inequalities still allow the new cross-deep far--far scale.

Take

\[
N\asymp T^2.
\tag{10.1}
\]

Choose

\[
\frac{145}{186}<\kappa<\frac45.
\tag{10.2}
\]

The interval is nonempty because

\[
\frac45-\frac{145}{186}
=rac{19}{930}>0.
\tag{10.3}
\]

Concentrate the selected-Palm weight on one formal legal root and put

\[
D\asymp T^\kappa.
\tag{10.4}
\]

## 10.1 Distinguished-\(p\) geometry

Let

\[
R=N^{3/5}=T^{6/5}.
\]

Place \(D\) lower actual \(p\)-labels in a central nonreflected interval, separated by \(\asymp R\).  This fits because

\[
D R
\asymp T^{\kappa+6/5}
\le T^2\asymp p
\]

precisely when \(\kappa\le4/5\).

Attach one activity target to each lower label at the same fixed gap

\[
h_0\asymp T.
\tag{10.5}
\]

Since \(T\ll R\), all different lower clusters, their targets, and every lower-to-foreign-target cross gap are \(\gg R\).  In particular the Q4955 fixed-\(p\) far condition and the new cross-deep condition

\[
b\gg T^{26/93}
\]

are simultaneously compatible.

The fixed-gap root cap permits this: at \(h_0\) there are only \(D<T\) pairs, while

\[
A_p(h_0)\le3(h_0-1)\asymp T.
\]

A length-\(T\) interval sees only \(O(1)\) packet zeros, well below \(T^{2/3}\).  The total distinguished-\(p\) zero count is

\[
O(D)\ll T^{4/3}=N^{2/3}.
\]

Reflection can be added in a disjoint interval at constant cost.

## 10.2 Reciprocal far coordinate

The Q4955 reciprocal threshold is

\[
N^{8/15}=T^{16/15}.
\]

At the level of the current reciprocal-fibre interfaces there is room for \(D\) lower quotient coordinates separated by this amount, because

\[
D\,T^{16/15}
=T^{\kappa+16/15}
=o(T^2)
\]

for every \(\kappa<14/15\), in particular for (10.2).

Thus the numerical interfaces permit all lower pairs to lie in the reciprocal far sector as well.

## 10.3 AD2 scale

The far--far active pair moment of this packet is

\[
D^2\asymp T^{2\kappa}.
\tag{10.6}
\]

For every fixed small

\[
0<\epsilon<\frac{19}{930},
\qquad
\kappa=\frac{145}{186}+\epsilon,
\]

this is

\[
T^{145/93+2\epsilon},
\]

strictly above the canonical AD2 gate.

The packet also survives the new cross-near deletion because the cluster separation \(T^{6/5}\) is much larger than \(T^{26/93}\).

## 10.4 What this packet does **not** claim

It does not construct actual Apéry zero sets or prove that the required shell primes realize the prescribed two-dimensional \((r,k)\) placement.  In particular, exact simultaneous distribution of the full reciprocal quotient coordinate is precisely one of the horizontal arithmetic issues not controlled by the current theorems.

It also makes no claim about alignment of the original four-copy phases.

The packet proves only the following scoped statement:

> The currently banked **upper-bound interfaces**—global zero count, fixed-gap root degree, subinterval zero count, reflection, finite reciprocal fibres, Q4955 short-pair deletions and Q5199 activity definitions—do not by themselves force the remaining cross-deep far--far AD2 moment below \(T^{145/93}\).

An actual arithmetic proof can still win by using information not represented in those interfaces; the exact q-lift and horizontal moving-zero distributions identified above are the natural places.

---

# 11. First exact failed implication

The earliest load-bearing failure in the hoped-for two-witness closure is not an exponent estimate.  It is the type error

\[
\boxed{
 r_i^+-r_i=\varepsilon_i h_i\pmod p
 \quad\not\Longrightarrow\quad
 m_i^+-m_i=\varepsilon_i h_i\pmod{pq}.}
\tag{11.1}
\]

The exact correction is

\[
 m_i^+-m_i
 =\varepsilon_i h_i+pJ_i.
\]

The free \(J_i\) is the missing q-coordinate.  Therefore:

- two activity witnesses give one common-\(p\) four-zero geometry;
- they do not give a second \(q\)-return geometry;
- they do not give six distinct same-characteristic zeros;
- they do not turn the Q4955 far reciprocal displacement into a short quantity;
- and they do not supply a second outer weight.

This is the first exact implication that fails.

---

# 12. Smallest new arithmetic lemma after the proved deletion

Fix \(\zeta>0\) and put

\[
K_\zeta=T^{26/93-\zeta}.
\tag{12.1}
\]

For each active oriented copy choose the canonical witness (2.14).  Let

\[
\mathfrak D_{2,\rm ff}^{\rm crossdeep}(T;K_\zeta)
\]

denote the literal Palm-weighted ordered pair count with

\[
|d|>N^{8/15},
\qquad
g_p^\diamond>N^{3/5},
\tag{12.2}
\]

and, in both directional orientation charts,

\[
\operatorname{dist}
(\text{target of copy 1},\text{lower label of copy 2})>K_\zeta,
\]

\[
\operatorname{dist}
(\text{target of copy 2},\text{lower label of copy 1})>K_\zeta.
\tag{12.3}
\]

All actual moving-zero tests, reciprocal labels, source masks, lower-copy distinctness and the single selected-root weight are retained.

Sections 6--7 prove that the complement of this cross-deep event is already subcritical for the \(38/31\) gate.

Therefore the smallest clean remaining arithmetic input is:

\[
\boxed{
\mathfrak D_{2,\rm ff}^{\rm crossdeep}(T;K_\zeta)
\ll T^{145/93-\delta+o(1)}M_0
}
\tag{FFAD2-crossdeep}
\]

for some fixed \(\delta>0\).

This is strictly narrower than Q5217's AD2 target: it has already removed Q4955's two short pair sectors and the new activity cross-near sector.  On its four-distinct part it also carries the nonzero determinant \(p\mid\Phi_{12}\ne0\); its only separate collision geometry is canonical-target coalescence.

A proof of `(FFAD2-crossdeep)` would close the canonical active-degree gate.  No current canonical theorem proves it.

For the ideal line, replace \(K_\zeta\) by

\[
T^{1/10-\zeta}
\]

and \(145/93\) by \(6/5\).

---

# 13. Insertion-ready LaTeX — proved statements only

The following blocks use only already-proved canonical facts plus deterministic selection from the literal active witness set.  No open Palm-dispersion premise is included.

```latex
\begin{lemma}[Two active witnesses and the reciprocal $p$-determinant]
\label{lem:at-active-two-witness-determinant}
Fix a selected Palm root
\[
 x=(\mathbf b,\delta,\ell_0),\qquad
 \mathbf b=(\sigma,p,q,p'),\qquad P=pq,
\]
and put
\[
 C_x=\langle\Delta_x\overline{p'}\rangle_P,
 \qquad A_x=C_x\pmod p.
\]
For $i=1,2$, let $\ell_i$ be a literal live moving-prime copy and
let $(h_i,\lambda_i)$ be a literal activity witness.  In a fixed
ordinary/reflected orientation write
\[
 r_i^+\equiv r_i+\varepsilon_i h_i\pmod p,
 \qquad \varepsilon_i\in\{\pm1\},
 \qquad T\le h_i<2T.
\]
Then
\[
 \varepsilon_i h_i\ell_i\lambda_i
 \equiv A_x(\ell_i-\lambda_i)\pmod p,
\]
and consequently
\[
 p\mid
 \Phi_{12}:=
 \varepsilon_1h_1\ell_1\lambda_1(\ell_2-\lambda_2)
 -\varepsilon_2h_2\ell_2\lambda_2(\ell_1-\lambda_1).
\]
If moreover the four physical moving primes are pairwise distinct and
$2T$ is smaller than every moving-prime shell element, then
$\Phi_{12}\ne0$.

Writing
\[
 m_i^+-m_i=\varepsilon_i h_i+pJ_i,
 \qquad e_i=k_i^+-k_i,
\]
one has the exact integer identity
\[
 \frac{(\lambda_i-\ell_i)m_i
       +\varepsilon_i\lambda_i h_i}{p}
 +\lambda_iJ_i
 =q e_i.
\]
\end{lemma}

\begin{proof}
The reciprocal completion gives
$\ell_i r_i\equiv A_x\pmod p$ and
$\lambda_i r_i^+\equiv A_x\pmod p$.  Substituting
$r_i^+=r_i+\varepsilon_i h_i$ and eliminating $r_i$ proves the first
congruence.  Eliminating $A_x$ between the two copies gives the
$\Phi_{12}$ divisibility.

If $\Phi_{12}=0$, the prime $\ell_1$ divides
$h_2\ell_2\lambda_2(\ell_1-\lambda_1)$.  Under the stated size and
four-distinct hypotheses it divides none of the first three factors,
so $\ell_1\mid(\ell_1-\lambda_1)$, forcing
$\ell_1=\lambda_1$, a contradiction.

Finally
$\lambda_i m_i^+-\ell_i m_i=pq e_i$.  Substituting
$m_i^+-m_i=\varepsilon_i h_i+pJ_i$ gives the displayed integer
identity after division by $p$.
\end{proof}

\begin{lemma}[One-witness cross-gap bound]
\label{lem:at-active-one-witness-crossgap}
Fix one selected Palm root, one literal source box and one member of the
finite orientation alphabet.  For every active lower moving-prime copy
choose deterministically one literal activity witness, for example the
lexicographically least pair $(h,\lambda)$.

Let $\mathscr P_{x,T}(B)$ be any subfamily of ordered pairs of distinct
active lower copies with the following property.  After ordering their
oriented actual $p$-labels as
\[
 0\le r<u\le p-1,
\]
use the chosen witness of the left copy,
\[
 T\le h<2T,\qquad r+h\in\mathcal Z_p,
\]
and assume the three labels are nonwrapping and distinct and
\[
 B\le |u-(r+h)|<2B,
 \qquad B\ge2.
\]
All additional literal source, moving-zero, reciprocal, prime-exclusion
and pair masks may be retained in $\mathscr P_{x,T}(B)$.
Then
\[
 \#\mathscr P_{x,T}(B)
 \ll N^{o(1)}T B\min(B,T).
\]
Consequently
\[
 \sum_xw_x\#\mathscr P_{x,T}(B)
 \ll N^{o(1)}T B\min(B,T)M_0.
\]
\end{lemma}

\begin{proof}
If $r<u<r+h$, put $b=r+h-u$ and $a=u-r=h-b$.
Then $r,r+a,r+a+b$ are three actual nonwrapping $p$-zeros, so
\[
 N_a(r)=N_b(r+a)=0\pmod p.
\]
If $r+h<u$, put $b=u-r-h$; then
$r,r+h,r+h+b$ are three actual zeros and
\[
 N_h(r)=N_b(r+h)=0\pmod p.
\]
In either case Lemma~\ref{lem:nonvanish} and
$\deg N_j=3(j-1)$ give at most
$3(\min(B,T)+O(1))$ possible starting residues for each fixed pair of
gap parameters.  There are $O(T)$ choices of $h$ and $O(B)$ choices of
$b$.

A fixed oriented $p$-label has only $O(1)$ moving-prime lifts in the
fixed shell because the reciprocal congruence fixes one residue class
modulo $p\asymp N$.  Thus restoring the two lower primes and the chosen
target costs only $N^{o(1)}$.  Every remaining literal mask is a
deletion.  This proves the result.
\end{proof}

\begin{corollary}[Cross-near deletion at the active-degree gate]
\label{cor:at-active-crossnear-gate}
In the Q4955 far--far lower-pair sector, let
$\mathfrak D_{2,\times\le K}^{\rm ff}(T)$ denote the Palm-weighted
ordered pair mass for which, in at least one of the two directional
orientation charts, the chosen activity target lies at nonwrapping
integer distance at most $K$ from the other lower $p$-label.  Then, for
$2\le K\le T$,
\[
 \mathfrak D_{2,\times\le K}^{\rm ff}(T)
 \ll N^{o(1)}TK^2M_0+N^{2/3+o(1)}M_0.
\]
The distance-one part is empty.  On the distance-zero part the far
fixed-$p$ condition forces $N^{3/5}<2T$, and hence
\[
 N^{2/3}\ll T^{10/9}.
\]
In particular, for every fixed $\zeta>0$, taking
\[
 K=T^{26/93-\zeta}
\]
gives
\[
 \mathfrak D_{2,\times\le K}^{\rm ff}(T)
 \ll
 T^{145/93-2\zeta+o(1)}M_0
 +T^{10/9+o(1)}M_0.
\]
Taking instead $K=T^{1/10-\zeta}$ gives
\[
 \mathfrak D_{2,\times\le K}^{\rm ff}(T)
 \ll
 T^{6/5-2\zeta+o(1)}M_0
 +T^{10/9+o(1)}M_0.
\]
\end{corollary}
```

---

# 14. Final ledger

```text
literal selected root weight w_x:                 RETAINED ONCE
lower moving copies ell_1, ell_2:                 ACTUAL two-test hits
far reciprocal displacement:                      |d| > N^(8/15)
far fixed-p diamond gap:                           > N^(3/5)
activity witness per lower:                        EXISTS, actual target
possible witnesses per lower:                      <= T^(2/3+o(1))
canonical witness selection:                       USED to avoid multiplicity
all-witness pair loss if expanded:                 up to T^(4/3+o(1))

one activity edge:
  eps h ell lambda = A(ell-lambda) mod p

two activity edges:
  p | Phi_12
  Phi_12 nonzero on four-distinct mesoscopic sector
  equal-h specialization:                          Q5155 cubic

full-P lift:
  m^+ - m = eps h + p J
  q-stage identity:                                B + lambda J = q e
  activity controls J:                             NO
  second independent q divisor:                    NOT PROVED

same-characteristic zeros created by activity:
  p:                                                four endpoints
  p':                                               one repeated inherited residue
  q:                                                no unselected zero test
  moving characteristics:                          one zero in each different field
  legal six-zero same-characteristic resultant:    NO

new proved one-witness cross-gap bound:
  dyadic b ~ B:                                    T B min(B,T) M0
  cumulative 2 <= b <= K <= T:                     T K^2 M0
  b = 1:                                            EMPTY
  b = 0 in far-p sector:                            T^(10/9+o(1)) M0

canonical gate cutoff:
  K = T^(26/93-zeta)
  exponent:                                         145/93 - 2 zeta

ideal cutoff:
  K = T^(1/10-zeta)
  exponent:                                         6/5 - 2 zeta

remaining cross-deep residual:
  |d| > N^(8/15)
  g_p^diamond > N^(3/5)
  both directional canonical cross gaps > T^(26/93-o(1))
  target coalescence OR four-distinct nonzero Phi_12

phase route:
  successful target creates trivial trace summand
  all-witness Fourier expansion costs T^(4/3)
  canonical witness has no linear averaging variable
  new horizontal operator theorem required

full far-far AD2 with rho < 145/186:               NOT PROVED
new restricted far-far deletion:                    PROVED
smallest remaining lemma:                            FFAD2-crossdeep
```

## Final conclusion

The additional activity information is not inert.  It has two precise consequences that were absent in Q5217:

- algebraically, two activity edges create the nonzero four-prime determinant \(p\mid\Phi_{12}\), but the missing \(q\)-row lift \(J_i\) blocks a second CRT modulus;
- combinatorially, **one** canonical activity edge turns an arbitrary second active lower label into a three-zero gap configuration.  Whenever the second lower lies within \(K\) of the target, the canonical gap-polynomial degree gives the sharp positive bound \(TK^2M_0\), enough to delete \(K=T^{26/93-\zeta}\) at the canonical gate and \(K=T^{1/10-\zeta}\) at the ideal line.

After those deletions, the obstruction is materially smaller: a Palm-weighted pair of far reciprocal/far spatial active lowers whose two canonical activity edges remain cross-deep.  On the generic four-distinct part it carries a nonzero characteristic-zero determinant; on the exceptional part the two activity edges coalesce at one target.  Neither configuration is controlled with a power saving by the present canonical interfaces.