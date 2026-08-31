ANSWER Q7397 97ecebed

# Problem 3.2 — three lower-density projections of the exact quotient-six C4 pre-upper load

## Verdict

Let

\[
L:=\log T,
\qquad
\mathcal B(T):=\frac{T^2}{L^2},
\]

so that \(\mathcal B(T)\) is the unweighted Boolean-incidence target corresponding to the weighted C4 target \(o_J(T^2/L)\).

For the exact raw Boolean union pre-upper load \(A_q(t)\), there are three natural decompositions:

1. the **global lower-plane projection**, with pair density
   \[
   \pi_p=\frac{(Z(p))_2}{(p)_2}
   \]
   and the indispensable double-sign correction
   \[
   \tau_p=\frac{(Z(p))_3}{(p)_3};
   \]
2. the **exact fixed-\((p,d)\) line projection**, with pair-line means
   \[
   \frac{C_p(7d)}p,\qquad \frac{C_p(6d)}p
   \]
   and exact triple-line correction \(H_p(d)/p\);
3. the **upper \(q\)-fibre projection**, with mean \(q^{-1}\sum_t A_q(t)\) tested against the upper zero density \(Z(q)/q\).

All three give exact centered covariance identities while retaining every prime, shell, row, long-gap, sign-window, and raw-OR filter. Their deterministic mean terms are all already subcritical from the current banks:

\[
M_{\mathrm{plane},0},\ M_{\mathrm{line},0},\ M_{q,0}
\ll_J
T^{-1/3}S_2(T)
\ll_J
T^{28/15}L^{-3/5}
=o_J\!\left(\frac{T^2}{L^2}\right),
\]

where

\[
S_2(T):=
\sum_{p\asymp_J T}Z(p)(Z(p)-1)
\ll_JT^{11/5}L^{-3/5}.
\]

Consequently the three centered covariances differ from each other only by \(o_J(T^2/L^2)\). At the C4 target scale they are equivalent gates.

The exact \(d\)-line mean is the formulation **closest** to the same-characteristic second factorial theorem, because

\[
\sum_{h\ne0}C_p(h)=Z(p)(Z(p)-1).
\]

This helps only with the deterministic mean. It does not control where the starts on a line occur, nor whether those starts align with the moving upper zero set \(\mathcal Z_{p+d}\) inside the exact window. An abstract line-incidence countermodel below keeps every line count and the complete second-factorial budget fixed while placing all occupied starts on upper-marked residues. Its centered line covariance is then of target size. Thus no linewise centering is strictly more accessible to the existing same-characteristic theorem.

In compact form:

```text
plane mean with tau correction:       CLOSED
exact d-line mean:                    CLOSED, most directly by SF2
q-fibre mean:                         CLOSED, cheapest analytically
plane centered covariance:            OPEN
line centered covariance:             OPEN
q-fibre centered covariance:          OPEN
any one centered covariance o(B):     equivalent to the C4 gate
linewise centered gain from SF2 alone: NO
```

---

## 0. Source and scope boundary

The connector-visible public `main` is the older source state

```text
c5d932b66ce5e4f1657b587215d290ae7a13018b.
```

It contains the canonical recurrence, reflection, zero-count, subinterval-zero, and factorial-moment theorems in `problems/3.2/proof.tex`. The late quotient-six C4 notes are not all present on that visible branch. I therefore use the exact formulas in the question together with the synchronized same-project interfaces:

- `Q4156 6dd57a00`: exact raw plus/minus supports, Boolean union, and raw minus-first reconciliation;
- `Q4159 ed2c3efe`: exact pre-upper load, upper Fourier split, and marked covariance obstruction;
- `Q4160 100453d7`: marked cross-lower-prime factorial sharpening and its critical exponent;
- `Q7287 7e36d696`: quotient stripping, full-window fibre one, and Boolean occurrence boundary;
- `Q7374 384705af`: exact plus fixed-\(d\) and projected upper-digit windows;
- the canonical same-characteristic bank
  \[
  S_1(T)\ll_JT^{8/5}L^{-4/5},
  \qquad
  S_2(T)\ll_JT^{11/5}L^{-3/5},
  \qquad
  \#(\mathcal Z_p\cap I)\ll |I|^{2/3}+1.
  \]

No finite computation and no independence assumption is used. This note counts distinct numerical lower primes in the exact Boolean raw union. It does not silently count physical rows, source presentations, or Palm copies.

---

## 1. Exact Boolean pre-upper object

For a prime \(m\), write

\[
z_m(u):=
\mathbf 1_{\mathcal Z_m}(u\bmod m),
\qquad
\mathcal Z_m=\{0\le v<m:m\mid b_v\},
\qquad
Z(m)=|\mathcal Z_m|.
\]

The proved Apéry reflection is

\[
z_p(u)=z_p(-1-u).
\tag{1.1}
\]

Fix a lower prime \(p\), a positive integer \(d\), and put

\[
q:=p+d.
\]

For the upper digit \(t\), put

\[
r=t+6d,
\qquad
a_-=d-t-1,
\qquad
a_+=t+12d-p.
\tag{1.2}
\]

Let

\[
\omega^-_{p,d}(t),\qquad \omega^+_{p,d}(t)\in\{0,1\}
\]

be the exact deterministic raw masks. Each mask retains literally:

- \(p,q=p+d\) prime;
- \(T/7<p<T/3\), \(p<q\), and \(d>T/J\);
- \(0\le t<q\) and \(T<6q+t\le2T\);
- the appropriate sign-specific digit inequalities;
- every fixed exceptional-prime, bounded-lift, shell, and raw-window predicate;
- the numerical-prime OR convention inside that sign.

The masks omit only the two lower Apéry-zero tests and the upper test \(q\mid b_t\). Extra masks can only thin the sets in every estimate below.

Using reflection on the primary lower digit and, in the minus branch, on the secondary lower digit, the two exact lower-pair indicators are

\[
B^-_{p,d}(t)
:=z_p(t+6d)z_p(t-d),
\tag{1.3-}
\]

\[
B^+_{p,d}(t)
:=z_p(t+6d)z_p(t+12d).
\tag{1.3+}
\]

Indeed,

\[
z_p(p-1-r)=z_p(r)=z_p(t+6d),
\]

\[
z_p(a_-)=z_p(-1-a_-)=z_p(t-d),
\]

and \(a_+\equiv t+12d\pmod p\).

The two signs share the first lower zero. Their exact simultaneous indicator is

\[
B^-_{p,d}(t)B^+_{p,d}(t)
=z_p(t+6d)z_p(t-d)z_p(t+12d).
\tag{1.4}
\]

The exact Boolean lower-prime incidence is therefore

\[
U_{p,d}(t)
:=
\omega^-_{p,d}(t)B^-_{p,d}(t)
+\omega^+_{p,d}(t)B^+_{p,d}(t)
-\omega^-_{p,d}(t)\omega^+_{p,d}(t)
 B^-_{p,d}(t)B^+_{p,d}(t).
\tag{1.5}
\]

This is inclusion-exclusion for the numerical-prime union, not a signed multiset count. Equivalently, it is what the raw minus-first tie-break counts once.

The exact pre-upper load is

\[
\boxed{
A_q(t)
:=
\sum_{\substack{p<q\\d=q-p}}U_{p,d}(t).
}
\tag{1.6}
\]

The remaining unweighted marked incidence is

\[
\boxed{
M^\vee(T,J)
:=
\sum_{q\ \mathrm{prime}}
\sum_{0\le t<q}z_q(t)A_q(t).
}
\tag{1.7}
\]

Every identity below is an exact decomposition of (1.7).

### 1.1 Exact fixed-\(d\) window lengths

The deterministic masks have the following basic support bounds before intersecting with the q6 row interval.

For minus,

\[
0\le t<\min\{d,p-6d\},
\tag{1.8-}
\]

so its real length is at most \(d\).

For plus,

\[
\max\{0,p-12d\}\le t<p-6d,
\tag{1.8+}
\]

so its real length is

\[
\min\{6d,p-6d\}\le6d.
\]

Both windows require \(6d<p\). Intersecting with the upper row interval and any additional masks only shortens them.

For later use define the actual upper-marked window occupancies

\[
N^\sigma_{p,d}
:=
\sum_t z_{p+d}(t)\omega^\sigma_{p,d}(t).
\tag{1.9}
\]

The proved subinterval theorem gives uniformly

\[
N^-_{p,d}\ll(d+1)^{2/3},
\qquad
N^+_{p,d}\ll(6d+1)^{2/3}.
\tag{1.10}
\]

If a bounded number of lift charts is retained separately, partition into \(O_J(1)\) intervals; only the implicit fixed-\(J\) constant changes.

---

## 2. Projection I: global lower-plane density and the \(\tau\) correction

For \(k\ge1\), write

\[
(x)_k=x(x-1)\cdots(x-k+1).
\]

On every actual row, \(6d<p\). For all sufficiently large lower primes, the three residues

\[
u,\qquad u-7d,\qquad u+6d
\]

are distinct modulo \(p\): a collision would force \(p\mid6d\), \(p\mid7d\), or \(p\mid13d\), hence only a fixed small exceptional characteristic.

The exact off-diagonal global densities are therefore

\[
\boxed{
\pi_p:=\frac{(Z(p))_2}{(p)_2},
\qquad
\tau_p:=\frac{(Z(p))_3}{(p)_3}.
}
\tag{2.1}
\]

Here \(\pi_p\) is the density of an ordered distinct zero pair in the complete lower plane, while \(\tau_p\) is the density of an ordered distinct zero triple. The latter is the required double-sign correction: the two sign events share one lower-zero slot.

Define

\[
\mu^{\mathrm{pl}}_{p,d}(t)
:=
\omega^-_{p,d}(t)\pi_p
+\omega^+_{p,d}(t)\pi_p
-\omega^-_{p,d}(t)\omega^+_{p,d}(t)\tau_p.
\tag{2.2}
\]

Then, identically,

\[
\begin{aligned}
U_{p,d}(t)
={}&\mu^{\mathrm{pl}}_{p,d}(t)\\
&+\omega^-_{p,d}(t)
  \bigl(B^-_{p,d}(t)-\pi_p\bigr)\\
&+\omega^+_{p,d}(t)
  \bigl(B^+_{p,d}(t)-\pi_p\bigr)\\
&-\omega^-_{p,d}(t)\omega^+_{p,d}(t)
  \bigl(B^-_{p,d}(t)B^+_{p,d}(t)-\tau_p\bigr).
\end{aligned}
\tag{2.3}
\]

Thus

\[
\boxed{
M^\vee=M_{\mathrm{pl},0}+\mathcal C_{\mathrm{pl}},
}
\tag{2.4}
\]

where

\[
M_{\mathrm{pl},0}
:=
\sum_{p,d,t}z_{p+d}(t)\mu^{\mathrm{pl}}_{p,d}(t),
\tag{2.5}
\]

and \(\mathcal C_{\mathrm{pl}}\) is the same sum with the three centered lines on the right of (2.3). Every prime and window filter remains inside \(\omega^\sigma\).

The \(-\tau_p\) term is not optional. Replacing the Boolean union by two independent signed planes changes the counted object on every double-sign row. For an upper bound on the deterministic mean, however, it is favorable:

\[
0\le\mu^{\mathrm{pl}}_{p,d}(t)
\le
\bigl(\omega^-_{p,d}(t)+\omega^+_{p,d}(t)\bigr)\pi_p,
\tag{2.6}
\]

because \(0\le\tau_p\le\pi_p\).

### 2.1 The plane mean is subcritical

Using (1.9), (1.10), and (2.6),

\[
\begin{aligned}
M_{\mathrm{pl},0}
&\le
\sum_p\pi_p\sum_{0<d<p/6}
  \bigl(N^-_{p,d}+N^+_{p,d}\bigr)\\
&\ll
\sum_p\pi_p\sum_{d<p/6}(d+1)^{2/3}\\
&\ll
\sum_p\frac{Z(p)(Z(p)-1)}{p^2}\,p^{5/3}\\
&\ll_J
T^{-1/3}S_2(T).
\end{aligned}
\tag{2.7}
\]

Therefore

\[
\boxed{
M_{\mathrm{pl},0}
\ll_JT^{28/15}L^{-3/5}
=o_J(\mathcal B(T)).
}
\tag{2.8}
\]

This uses the same-characteristic second factorial theorem and the proved upper zero-set subinterval theorem. It does not use any lower/upper independence.

### 2.2 Convention with replacement

If one calls

\[
\rho_p^2=\left(\frac{Z(p)}p\right)^2,
\qquad
\rho_p^3=\left(\frac{Z(p)}p\right)^3
\]

the plane and triple densities, the decomposition remains exact after replacing \(\pi_p,\tau_p\) everywhere by those constants. The finite-population difference is smaller:

\[
|\pi_p-\rho_p^2|\le\frac{Z(p)}{p^2},
\]

and its total contribution is

\[
\ll_JT^{-1/3}S_1(T)=o_J(\mathcal B(T)).
\]

The factorial convention (2.1) is preferable because the actual C4 slots are off-diagonal.

---

## 3. Projection II: exact mean on each \((p,d)\)-line

Define the cyclic lower pair correlations

\[
C_p(h)
:=
\sum_{u\bmod p}z_p(u)z_p(u+h).
\tag{3.1}
\]

Then

\[
C_p(-h)=C_p(h),
\qquad
\sum_{h\ne0}C_p(h)=Z(p)(Z(p)-1).
\tag{3.2}
\]

For the simultaneous two-sign event define

\[
H_p(d)
:=
\sum_{u\bmod p}
 z_p(u)z_p(u-7d)z_p(u+6d).
\tag{3.3}
\]

After the substitution \(u=t+6d\), the exact full-line means are

\[
\lambda^-_{p,d}:=\frac{C_p(7d)}p,
\qquad
\lambda^+_{p,d}:=\frac{C_p(6d)}p,
\qquad
\theta_{p,d}:=\frac{H_p(d)}p.
\tag{3.4}
\]

The exact line-projected Boolean mean is

\[
\mu^{\mathrm{line}}_{p,d}(t)
:=
\omega^-_{p,d}(t)\lambda^-_{p,d}
+\omega^+_{p,d}(t)\lambda^+_{p,d}
-\omega^-_{p,d}(t)\omega^+_{p,d}(t)\theta_{p,d}.
\tag{3.5}
\]

Consequently

\[
\begin{aligned}
U_{p,d}(t)
={}&\mu^{\mathrm{line}}_{p,d}(t)\\
&+\omega^-_{p,d}(t)
  \bigl(B^-_{p,d}(t)-\lambda^-_{p,d}\bigr)\\
&+\omega^+_{p,d}(t)
  \bigl(B^+_{p,d}(t)-\lambda^+_{p,d}\bigr)\\
&-\omega^-_{p,d}(t)\omega^+_{p,d}(t)
  \bigl(B^-_{p,d}(t)B^+_{p,d}(t)-\theta_{p,d}\bigr).
\end{aligned}
\tag{3.6}
\]

Thus

\[
\boxed{
M^\vee=M_{\mathrm{line},0}+\mathcal C_{\mathrm{line}},
}
\tag{3.7}
\]

with

\[
M_{\mathrm{line},0}
:=
\sum_{p,d,t}z_{p+d}(t)\mu^{\mathrm{line}}_{p,d}(t),
\tag{3.8}
\]

and the centered covariance obtained by summing the last three lines of (3.6) against \(z_{p+d}(t)\).

This is an exact mean on the **complete** \(p\)-line. The window and prime filters are not averaged away; they remain as \(\omega^\sigma_{p,d}(t)\) after centering. A “mean” conditioned on the occupied window or on the upper mark would already contain the unknown concentration and would not be supplied by the second factorial theorem.

### 3.1 The line mean is subcritical

Since the triple event is a subset of either pair event,

\[
0\le H_p(d)\le\min\{C_p(6d),C_p(7d)\}.
\]

Therefore

\[
0\le\mu^{\mathrm{line}}_{p,d}(t)
\le
\omega^-_{p,d}(t)\frac{C_p(7d)}p
+\omega^+_{p,d}(t)\frac{C_p(6d)}p.
\tag{3.9}
\]

Using (1.10),

\[
M_{\mathrm{line},0}
\ll
\sum_p\frac1p
\sum_{d<p/6}
(d+1)^{2/3}
\bigl(C_p(6d)+C_p(7d)\bigr).
\tag{3.10}
\]

Since \((d+1)^{2/3}/p\ll p^{-1/3}\), and multiplication by \(6\) or \(7\) is injective on the legal \(d\)-range modulo every sufficiently large \(p\),

\[
\sum_{d<p/6}
\bigl(C_p(6d)+C_p(7d)\bigr)
\le
2\sum_{h\ne0}C_p(h)
=2Z(p)(Z(p)-1).
\tag{3.11}
\]

Hence

\[
\boxed{
M_{\mathrm{line},0}
\ll_JT^{-1/3}S_2(T)
\ll_JT^{28/15}L^{-3/5}
=o_J(\mathcal B(T)).
}
\tag{3.12}
\]

This is the most literal use of the same-characteristic second factorial theorem among the three projections.

---

## 4. Projection III: the upper \(q\)-fibre mean

Keep the exact Boolean pre-upper load \(A_q(t)\) from (1.6). Put

\[
A_q^{\mathrm{tot}}:=\sum_{0\le u<q}A_q(u),
\qquad
\overline A_q:=\frac{A_q^{\mathrm{tot}}}{q},
\tag{4.1}
\]

\[
\rho_q:=\frac{Z(q)}q,
\qquad
g_q(t):=z_q(t)-\rho_q.
\tag{4.2}
\]

Since \(\sum_{t\bmod q}g_q(t)=0\),

\[
\boxed{
M^\vee=M_{q,0}+\mathcal C_q,
}
\tag{4.3}
\]

where

\[
M_{q,0}
:=
\sum_q\rho_qA_q^{\mathrm{tot}}
=
\sum_q Z(q)\overline A_q,
\tag{4.4}
\]

and

\[
\boxed{
\mathcal C_q
:=
\sum_q\sum_{t\bmod q}
 g_q(t)\bigl(A_q(t)-\overline A_q\bigr).
}
\tag{4.5}
\]

Equivalently,

\[
M^\vee
=
\sum_q Z(q)\overline A_q
+
\sum_q\sum_{t\in\mathcal Z_q}
 \bigl(A_q(t)-\overline A_q\bigr).
\tag{4.6}
\]

All lower prime and window filters are already inside \(A_q(t)\). This is the genuine pre-upper split: imposing \(z_q(t)\) inside the definition of \(A_q(t)\) would make the centered term circular.

### 4.1 The q-fibre mean is subcritical

By the fixed-sign ordered-lower-pair inverse, forgetting the upper mark maps every sign-specific candidate injectively to an ordered distinct pair of zeros in one lower characteristic. The Boolean union only decreases the sign sum. Therefore

\[
\sum_qA_q^{\mathrm{tot}}
\ll_J S_2(T).
\tag{4.7}
\]

The pointwise zero theorem gives, uniformly for \(q\asymp_JT\),

\[
\rho_q=\frac{Z(q)}q\ll_JT^{-1/3}.
\tag{4.8}
\]

Thus

\[
\boxed{
M_{q,0}
\ll_JT^{-1/3}S_2(T)
\ll_JT^{28/15}L^{-3/5}
=o_J(\mathcal B(T)).
}
\tag{4.9}
\]

This is analytically the cheapest deterministic mean: it uses the pre-upper lower-pair bound and the pointwise upper zero density, but not the sign-specific upper subinterval theorem.

---

## 5. Comparison of the three covariances

The deterministic terms satisfy

\[
M_{\mathrm{pl},0},\quad
M_{\mathrm{line},0},\quad
M_{q,0}
=o_J(\mathcal B(T)).
\tag{5.1}
\]

Since all three decompositions equal the same exact \(M^\vee\),

\[
\mathcal C_{\mathrm{pl}}-\mathcal C_{\mathrm{line}}
=M_{\mathrm{line},0}-M_{\mathrm{pl},0}
=o_J(\mathcal B(T)),
\tag{5.2}
\]

and

\[
\mathcal C_{\mathrm{line}}-\mathcal C_q
=M_{q,0}-M_{\mathrm{line},0}
=o_J(\mathcal B(T)).
\tag{5.3}
\]

Therefore

\[
\boxed{
M^\vee=o_J(\mathcal B(T))
\iff
\mathcal C_{\mathrm{pl}}=o_J(\mathcal B(T))
\iff
\mathcal C_{\mathrm{line}}=o_J(\mathcal B(T))
\iff
\mathcal C_q=o_J(\mathcal B(T)).
}
\tag{5.4}
\]

The equivalence is at the target scale, not an assertion that the three centered summands agree term by term.

| Projection | Exact deterministic mean | Current input paying it | Bound | Centered status |
|---|---|---|---:|---|
| global plane | \(\pi_p+\pi_p-\tau_p\), with exact masks | lower SF2 + upper subinterval zeros | \(\ll_JT^{-1/3}S_2\) | open |
| exact \(d\)-line | \(C_p(7d)/p+C_p(6d)/p-H_p(d)/p\) | lower SF2 + upper subinterval zeros | \(\ll_JT^{-1/3}S_2\) | open |
| upper q-fibre | \((Z(q)/q)\sum_tA_q(t)\) | lower SF2/fibre one + pointwise \(Z(q)\) | \(\ll_JT^{-1/3}S_2\) | open |

The negative double-sign corrections need no new triple-factorial theorem for the mean: \(\tau_p\le\pi_p\) and \(H_p(d)\le C_p(6d),C_p(7d)\). They must nevertheless remain in the exact identities.

---

## 6. Why the exact line mean does not make the centered gate easier

The same-characteristic second factorial theorem controls the line masses

\[
C_p(h)=\sum_{u\bmod p}B_{p,h}(u)
\]

only after summing over starts, gaps, and lower primes. It does not control the location of the occupied starts on a fixed line.

For one abstract line, let

\[
B(u)\in\{0,1\},
\qquad
m:=\sum_{u\bmod p}B(u),
\]

and let \(E\subseteq\mathbb F_p\) be the exact set selected by the upper mark and all retained windows. The line-centered test is

\[
\boxed{
\sum_{u\in E}\left(B(u)-\frac mp\right)
=|E\cap\operatorname{supp}B|-\frac{m|E|}{p}.
}
\tag{6.1}
\]

The value \(m\), and hence every contribution to the second factorial theorem, is unchanged if the support of \(B\) is moved around the line. It is compatible with

\[
\operatorname{supp}B\subseteq E,
\]

in which case (6.1) equals

\[
m\left(1-\frac{|E|}{p}\right),
\]

which is essentially \(m\) whenever \(|E|=o(p)\). It is equally compatible with disjoint support, producing a negative covariance. Thus a line count supplies no discrepancy estimate against an adaptive marked window.

### 6.1 Aggregate abstract countermodel at the C4 scales

This can be scaled to all current numerical banks.

Put

\[
U_0:=T^{8/5}L^{-4/5},
\qquad
K_0:=\frac{\mathcal B(T)}{U_0}
=T^{2/5}L^{-6/5}.
\]

Construct an abstract plus-only incidence system as follows.

1. Choose \(\asymp U_0\) legal upper marked states in a fixed interior plus region where all prime, row, shell, long-gap, and digit windows are nonempty.
2. Give every such state \(\asymp K_0\) distinct lower-prime tokens. The exact Boolean union then has
   \[
   M^\vee\asymp U_0K_0=\mathcal B(T).
   \]
3. Assign every token to a distinct legal \((p,d)\)-line and put exactly one occupied lower-pair start on that line, at the selected upper residue. Thus \(C_p(6d)=1\) on each used line and zero on unused lines.
4. Use only the plus sign, so both the actual overlap and every \(\tau\)-correction vanish.
5. Distribute the marked upper residues through their macroscopic legal windows so that every pointwise and subinterval upper-zero bound holds.

The total lower pair mass is only

\[
\asymp\mathcal B(T)
=o\!\left(S_2(T)\right),
\]

because

\[
\frac{S_2(T)}{\mathcal B(T)}
\asymp T^{1/5}L^{7/5}\to\infty.
\]

The upper marked-state mass is \(\asymp U_0\), exactly within the banked first-moment scale. There are \(\gg T^2/L\) possible lower prime-gap lines, so using only \(\asymp T^2/L^2\) singleton lines respects the line-capacity interface. The tokens can be distributed over lower characteristics with about \(T/L\) pairs per characteristic, compatible with zero-set size about \((T/L)^{1/2}\ll T^{2/3}\) and with the exported subinterval envelope. This is an interface-level model, not a claim that the actual Apéry zeros realize it.

On every used line, the exact line mean is \(1/p\), while the selected centered contribution is \(1-1/p\). Hence

\[
M_{\mathrm{line},0}=o(\mathcal B(T)),
\qquad
\mathcal C_{\mathrm{line}}\asymp\mathcal B(T).
\]

The q-fibre load is likewise concentrated on the marked residues, so

\[
M_{q,0}=o(\mathcal B(T)),
\qquad
\mathcal C_q\asymp\mathcal B(T).
\]

The plane mean is subcritical by (2.8), and therefore

\[
\mathcal C_{\mathrm{pl}}\asymp\mathcal B(T)
\]

as well. Every deterministic mean theorem and the entire same-characteristic SF2 budget survive, but none of the centered gates closes.

This proves the required logical obstruction:

> The same-characteristic second factorial theorem controls how many lower zero pairs exist. It does not control whether their affine starts are placed on the upper-marked residues selected by \(q=p+d\) and the exact windows.

No actual Apéry counterexample is asserted. The point is that a deduction using only the exported second-factorial and subinterval interfaces cannot rule the model out.

### 6.2 Why window-conditioned line means are circular

One might try to replace \(C_p(cd)/p\) by the average of \(B^\sigma_{p,d}(t)\) over the exact legal window. That average is

\[
\frac{
\sum_t\omega^\sigma_{p,d}(t)B^\sigma_{p,d}(t)
}{
\sum_t\omega^\sigma_{p,d}(t)
}.
\]

The existing second factorial theorem does not estimate it: all line pairs may lie in the window. If the upper mark \(z_{p+d}(t)\) is also included in the conditioning, the resulting “mean term” is simply the marked incidence being estimated. The covariance becomes small only because the target has been moved into the mean. This is not a weaker route.

---

## 7. Final claim boundary

The following statements are proved by the current banks and the exact geometry:

- all three projection identities are exact and occurrence-preserving for the Boolean raw union;
- prime, row, long-gap, and sign-window filters remain literal in every covariance;
- the global-plane \(\tau_p\) term and the exact line triple term are the correct double-sign corrections;
- every deterministic mean is \(o_J(T^2/L^2)\);
- the exact d-line mean is paid most directly by the same-characteristic SF2 theorem;
- the three centered covariances are target-equivalent modulo subcritical deterministic terms;
- SF2 alone does not control any centered covariance, by (6.1) and the aggregate countermodel.

What is not proved is a centered cross-characteristic discrepancy estimate. A genuinely new input must constrain the alignment

\[
\text{lower pair start on its exact }(p,d)\text{-line}
\quad\text{versus}\quad
\mathcal Z_{p+d}	ext{ in the retained upper window}.
\]

Changing the lower-density projection does not remove that alignment problem.