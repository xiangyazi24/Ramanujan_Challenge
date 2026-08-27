ANSWER Q4148 2eab0efa

# Verdict

**The prescribed excess sum is not closed unconditionally by the current P3.2 inputs.**  After retaining the exact q6 geometry, the strongest current shallow/short-window deletions, the genuine selected-upper condition, and every indexed occurrence, the best unconditional bound remains

\[
\boxed{
\mathfrak X_\vee(R)
:=\sum_{(q,t)}\bigl(\kappa_\vee(q,t)-2R\bigr)_+
\ll_J
\frac{T^{11/5}}{(\log T)^{8/5}}
+
\frac{T^2}{\log T}.
}
\tag{V1}
\]

Here and below

\[
L:=\log T,
\qquad
R\asymp_J\frac{T^{2/5}}{L^{6/5}\Omega(T)},
\qquad
\Omega(T)\to\infty,
\qquad
\Omega(T)=T^{o(1)}.
\]

The leading term in `(V1)` exceeds the required scale \(T^2/L^2\) by

\[
\boxed{T^{1/5}L^{2/5},}
\tag{V2}
\]

and the second term exceeds it by \(L\).  The corrected lower-zero second factorial moment is weaker here by one logarithm.

There are genuine positive reductions.

1. The fixed-state carrier simplifies exactly to
   \[
   \boxed{
   G_\vee(q,t)
   =\gcd\!\left(
      \frac{b_n}{q},
      \operatorname{lcm}\!\left(\frac{b_s}{q},\frac{b_N}{q}\right)
   \right),
   }
   \tag{V3}
   \]
   where
   \[
   n=6q+t,\qquad s=q-1-t,\qquad N=12q+t.
   \]
   Thus every counted foreign prime divides the one fixed value \(b_n/q\).

2. Outside the fixed finite set of primes dividing \(5b_6b_{13}\), the sign-tagged support is exact.  If \(d=q-p\), then
   \[
   r=t+6d=n-6p,
   \]
   and the two secondary indices are
   \[
   a_-=d-1-t=s-p,
   \qquad
   a_+=t+12d-p=N-13p.
   \]
   On the exact minus window,
   \[
   p\mid G_-(q,t)\iff p\mid b_r,b_{a_-};
   \]
   on the exact plus window,
   \[
   p\mid G_+(q,t)\iff p\mid b_r,b_{a_+}.
   \]

3. The free quota really is negligible:
   if \(B_1(T)\) is the number of actual selected upper occurrences, then
   \[
   B_1(T)\ll_J T^{8/5}L^{-4/5}
   \]
   and therefore
   \[
   \boxed{
   2R B_1(T)
   \ll_J
   \frac{T^2}{L^2\Omega(T)}
   =o_J(T^2/L^2).
   }
   \tag{V4}
   \]
   Equivalently, one may reserve the first and last \(R\) carrier primes in every fixed state at negligible aggregate cost.  The live excess consists only of interior rank-\(R\) divisors.

4. The proved shallow plus collar and the elementary short-geometric-window collar are subcritical with full state multiplicity.  After they are removed, however, legal boxes with all digits and both window depths comparable with \(T\) remain.  No current local deletion bounds their fixed-state large-prime divisor excess.

The first precise missing arithmetic estimate is the following **deep legal-window radical tail**.  Let \(P_-(T)\asymp_JT\) be the lower edge of the retained lower-prime shell, remove the upper prime and the fixed exceptional primes from `(V3)`, and let

\[
\mathcal R_\vee(q,t)
:=
\prod_{\substack{p\text{ in the exact retained plus-or-minus window}\\
                  p\ne q,\ p\mid G_\vee(q,t)}}p.
\tag{V5}
\]

It would suffice to prove

\[
\boxed{
\sum_{(q,t)}
\left[
  \log \mathcal R_\vee(q,t)
  -2R\log P_-(T)
\right]_+
=o_J\!\left(\frac{T^2}{L}\right).
}
\tag{VEE-RADICAL-EXCESS}
\]

Indeed,

\[
\bigl(\kappa_\vee(q,t)-2R\bigr)_+\log P_-(T)
\le
\left[
  \log \mathcal R_\vee(q,t)
  -2R\log P_-(T)
\right]_+,
\]

and \(\log P_-(T)\asymp L\).

A stronger reusable sufficient theorem is the fixed-state average codegree bound

\[
\boxed{
\sum_{\ell\in\mathscr A_\vee(T)}D_\vee(\ell)
=o_J\!\left(
   \frac{T^{1/5}}{\Omega(T)L^{8/5}}
   \left|\mathscr A_\vee(T)\right|
\right),
}
\tag{VEE-AVG-CODEGREE}
\]

where \(\mathscr A_\vee(T)\) is the complete multiset of exact legal incidences \(\ell=((q,t),p)\), and

\[
D_\vee((q,t),p):=\kappa_\vee(q,t)-1.
\]

This implies the required excess through the second factorial moment.  A fixed-power version

\[
\sum_\ell D_\vee(\ell)
\ll_J T^{1/5-\varepsilon}|\mathscr A_\vee(T)|
\]

for any fixed \(\varepsilon>0\) is more than enough.

The first obstruction is already present in the deep **plus-only** subfamily: the minus carrier may be empty, while one fixed plus gcd can have many legal shell-prime divisors.  Different lower primes live in different lower characteristics, so neither the edge-disjoint C4 packing nor the complete-fibre theorem forces a second fixed integer or a repeated prime power.  The newest fixed-state C4 carrier theorem makes the selected prime product squarefree, but squarefreeness is exactly what `(V5)` already records; it supplies no aggregate tail estimate.

No finite scan, independence assumption, GRH input, or repository-source edit is used.

---

# 0. Source-state boundary

The connector-visible repository heads used for this audit are

```text
main:
  c5d932b66ce5e4f1657b587215d290ae7a13018b

chatgpt-drop before this delivery:
  fa30892fd23515020692e3889fe6f3f337c93d11
```

The newest visible fixed-upper-state carrier audit is

```text
drops/Q4145-e4816a6b.md
```

and the immediately relevant same-project interfaces include

```text
Q3869 23068ed1
  shallow plus collar, fixed-upper short-window count,
  corrected B1/B2 scales, and common-value first moment

Q3970 4f483040
  same-sign fixed-gap and terminal-overlap collars

Q4041 2ec4319b / Q4042 56308b04
  exact fixed-state q-free two-value gcd support

Q4065 e940b0a9 / Q4066 4c8c7181
  exact P/M inverse forms and common-n additive-large-sieve bound

Q4097 10140349
  fixed-fibre pair-energy obstruction and critical average codegree

Q4144 3eebaa08 / Q4145 e4816a6b
  simple upper-state incidence conversion and
  branchwise squarefree packed-cycle radical carrier
```

The connector does not expose the caller's local `/Users/huangx/...` working tree byte-for-byte.  In accordance with the request not to use Python, the sandbox, or local-file materialization, I do not claim unseen local line numbers.  The exact reduction and carrier formulas in Q4148 are audited directly against the current connected project state above.

---

# 1. Exact q-free fixed-state carrier

Fix an actual selected upper occurrence

\[
\sigma=(q,t),
\qquad
q\text{ prime},
\qquad
0\le t<q,
\qquad
q\mid b_t.
\tag{1.1}
\]

Put

\[
n:=6q+t,
\qquad
s:=q-1-t,
\qquad
N:=12q+t.
\tag{1.2}
\]

Since \(b_0=1\), an actual selected occurrence has \(t\ge1\).

## 1.1 All three quotients are integral

For all sufficiently large upper-shell primes, the quotient digits \(6\) and \(12\) are below \(q\).  Gessel--Lucas gives

\[
b_n\equiv b_6b_t\equiv0\pmod q,
\tag{1.3}
\]

\[
b_N\equiv b_{12}b_t\equiv0\pmod q.
\tag{1.4}
\]

Apéry reflection gives

\[
b_s=b_{q-1-t}\equiv b_t\equiv0\pmod q.
\tag{1.5}
\]

Thus

\[
A:=\frac{b_n}{q},
\qquad
B:=\frac{b_s}{q},
\qquad
C:=\frac{b_N}{q}
\tag{1.6}
\]

are integers.  Dividing by one copy of \(q\), or by the full \(q\)-primary part, has exactly the same support at every foreign lower prime \(p\ne q\).

Define

\[
G_-:=\gcd(A,B),
\qquad
G_+:=\gcd(A,C),
\tag{1.7}
\]

and

\[
G_\vee:=\operatorname{lcm}(G_-,G_+).
\tag{1.8}
\]

## 1.2 Exact lcm--gcd simplification

For every prime \(\ell\), if

\[
a=v_\ell(A),\quad b=v_\ell(B),\quad c=v_\ell(C),
\]

then

\[
v_\ell(G_\vee)
=
\max\{\min(a,b),\min(a,c)\}
=
\min\{a,\max(b,c)\}.
\]

Therefore

\[
\boxed{
G_\vee
=
\gcd\bigl(A,\operatorname{lcm}(B,C)\bigr).
}
\tag{1.9}
\]

In particular,

\[
\boxed{G_\vee\mid A=b_n/q.}
\tag{1.10}
\]

This is the strongest immediate characteristic-zero compression: every foreign lower prime in the combined carrier divides one fixed Apéry coefficient.  It does not improve the height, because

\[
\log A\ll_JT.
\tag{1.11}
\]

---

# 2. Exact plus and minus windows

Let \(p<q\) be a foreign lower-shell prime and put

\[
d:=q-p.
\tag{2.1}
\]

The common lower primary digit is

\[
\boxed{r:=t+6d=n-6p.}
\tag{2.2}
\]

Because \(p<q\) and \(t\ge0\), the condition \(r\ge0\) is automatic.  The strict upper digit condition \(r<p\) is

\[
\boxed{p>n/7.}
\tag{2.3}
\]

## 2.1 Minus

The minus secondary digit is

\[
\boxed{a_-:=d-1-t=s-p.}
\tag{2.4}
\]

The conditions \(0\le a_-<p\) are

\[
s/2<p\le s.
\tag{2.5}
\]

But

\[
\frac n7-rac s2
=
\frac{5q+9t+7}{14}>0.
\]

Hence the common primary lower bound dominates, and the exact raw minus prime window is

\[
\boxed{
I_-(q,t)
=
\left\{p\text{ prime}:\frac n7<p\le s\right\}.
}
\tag{2.6}
\]

It is nonempty only when

\[
q-7-8t>0,
\tag{2.7}
\]

and its continuous length is

\[
\boxed{|I_-|_{\rm cont}=\frac{q-7-8t}{7}.}
\tag{2.8}
\]

All lower-shell, parity, primitive, close-shift, shallow, and later residual filters are to be intersected with `(2.6)`.

## 2.2 Plus

The plus secondary digit is

\[
\boxed{a_+:=t+12d-p=N-13p.}
\tag{2.9}
\]

The conditions \(0\le a_+<p\) are

\[
N/14<p\le N/13.
\tag{2.10}
\]

Since \(t\ge1\),

\[
\frac n7-rac N{14}=rac t{14}>0.
\]

Thus the exact raw plus window is

\[
\boxed{
I_+(q,t)
=
\left\{p\text{ prime}:\frac n7<p\le\frac N{13}\right\}.
}
\tag{2.11}
\]

Its continuous length is

\[
\boxed{|I_+|_{\rm cont}=\frac{6(q-t)}{91}.}
\tag{2.12}
\]

Again all inherited residual filters are retained.

Both windows have the same lower endpoint \(n/7\).  Their union is one interval ending at

\[
\max\{s,N/13\},
\tag{2.13}
\]

before the sign-specific filters are imposed.

---

# 3. Exact sign support outside a fixed finite set

Let

\[
E_0:=\operatorname{rad}(5b_6b_{13}).
\tag{3.1}
\]

Every lower-shell prime avoids \(E_0\) for all sufficiently large \(T\).

For \(p\in I_-(q,t)\), one has

\[
n=6p+r,
\qquad
s=p+a_-.
\]

Therefore Lucas gives

\[
b_n\equiv b_6b_r\pmod p,
\qquad
b_s\equiv b_1b_{a_-}=5b_{a_-}\pmod p.
\]

Since \(q\) is a \(p\)-unit,

\[
\boxed{
 p\mid G_-
 \iff
 p\mid b_r\text{ and }p\mid b_{a_-}.
}
\tag{3.2}
\]

For \(p\in I_+(q,t)\),

\[
N=13p+a_+,
\]

and hence

\[
b_N\equiv b_{13}b_{a_+}\pmod p.
\]

Thus

\[
\boxed{
 p\mid G_+
 \iff
 p\mid b_r\text{ and }p\mid b_{a_+}.
}
\tag{3.3}
\]

These are exact occurrence-preserving support statements.  Fixed \((q,t,p)\) determines \(d,r,a_-,a_+\), so there is no hidden fibre.

## 3.1 Sign-tagged versus literal lcm overcarrier

The most favorable exact interpretation of the prescribed carrier is

\[
\mathcal P_-(q,t)
:=\{p\in I_-^{\rm res}(q,t):p\mid G_-\},
\]

\[
\mathcal P_+(q,t)
:=\{p\in I_+^{\rm res}(q,t):p\mid G_+\},
\]

\[
\mathcal P_\vee(q,t)
:=\mathcal P_-(q,t)\cup\mathcal P_+(q,t).
\tag{3.4}
\]

Then

\[
\kappa_\vee(q,t)=|\mathcal P_\vee(q,t)|.
\]

If instead the literal proposed object first takes the lcm and then intersects its support with the untagged union window, a cross term can occur: a prime in the minus-only part of the window may divide \(G_+\), or conversely, without satisfying the corresponding sign's digit legality.  Such an object is a valid overcarrier, but the actual lower-zero-pair moment does not control those cross terms.  The common-\(n\) bound below still controls them because `(1.10)` remains true.

The main negative verdict is proved under the more favorable sign-tagged interpretation `(3.4)`.  It therefore remains valid for the larger literal overcarrier.

---

# 4. Multiplicity audit

Every relevant multiplicity is as follows.

1. **Within one fixed state and one sign.**  A lower prime \(p\) determines one signed candidate.  It contributes at most one support prime.
2. **Both signs at one state.**  The lcm/radical counts the same numerical \(p\) once.  There is no factor two from a double-sign row.
3. **Prime powers.**  \(\kappa_\vee\) counts distinct foreign primes.  Extra \(p\)-adic valuation in a gcd gives no additional occurrence.
4. **The upper prime.**  The numerical prime \(q\) is excluded.  Any residual \(q\)-primary valuation after dividing by one copy of \(q\) is irrelevant.
5. **Repeated lower prime across states.**  The same numerical \(p\) used at two different selected occurrences is counted twice, as it represents two indexed state incidences.
6. **Repeated numerical upper prime.**  Distinct actual zeros \(t\) of one \(q\) are distinct selected occurrences.  They must not be collapsed.
7. **Reflection.**  If both \(t\) and \(q-1-t\) occur in the selected-state sum, the proposed reduction counts both unless it explicitly works with typed reflection orbits.  This report follows the displayed sum over `(q,t)` and does not divide by two.
8. **Fixed exceptional primes.**  They disappear from the macroscopic shell; retaining them would cost only an absolute constant per state, hence \(O(B_1)=o(T^2/L^2)\).

The external factor \(1/2\) in the assumed deficiency reduction changes no exponent.

---

# 5. Exact sign split

Write

\[
\kappa_-(q,t):=|\mathcal P_-(q,t)|,
\qquad
\kappa_+(q,t):=|\mathcal P_+(q,t)|.
\]

Then

\[
\kappa_\vee\le\kappa_-+\kappa_+.
\tag{5.1}
\]

For nonnegative \(x,y\),

\[
(x+y-2R)_+
\le
(x-R)_++(y-R)_+.
\tag{5.2}
\]

Therefore

\[
\boxed{
\mathfrak X_\vee(R)
\le
\sum_{(q,t)}(\kappa_-(q,t)-R)_+
+
\sum_{(q,t)}(\kappa_+(q,t)-R)_+.
}
\tag{5.3}
\]

This is multiplicity-safe even when the same \(p\) supports both signs, because the right side is an upper bound.

An intersection theorem cannot close the union.  A prime supporting both signs divides

\[
\gcd(A,B,C),
\]

and supplies three lower zeros, but

\[
|\mathcal P_-\cup\mathcal P_+|
=|\mathcal P_-|+|\mathcal P_+|-|\mathcal P_-\cap\mathcal P_+|.
\]

A bound on the intersection only decreases the union count.  No current theorem forces a positive proportion of one sign's support to lie in the other.  In particular, the minus set may be empty while the plus set is large.  The plus-only family is already the first live obstruction.

The same warning applies to the reflected state \(q-1-t\): continuant or gcd identities for primes occurring in both reflected branches control only their intersection.

---

# 6. What the current local deletions really remove

Let

\[
H:=\lfloor R\rfloor.
\]

Replacing \(R\) by \(H\) changes the positive part by at most an absolute constant per state, hence by

\[
O(B_1)=o(T^2/L^2).
\]

## 6.1 The free fixed-state quota

Order the combined support in one state:

\[
p_1<p_2<\cdots<p_k,
\qquad k=\kappa_\vee(q,t).
\]

Reserve

\[
p_1,\dots,p_H,
\qquad
p_{k-H+1},\dots,p_k,
\]

when they exist.  The total number reserved over all states is at most

\[
2HB_1(T)
\ll_J
\frac{T^2}{L^2\Omega(T)}.
\tag{6.1}
\]

The remaining indices

\[
H<j\le k-H
\tag{6.2}
\]

are exactly \((k-2H)_+\) interior carrier incidences.  Thus the assumed \(2R\) subtraction is an exact endpoint trim at negligible aggregate cost.

It does not bound the interior.  Every live incidence has at least \(H\) carrier primes on either side in the same fixed state.

## 6.2 Rank-\(H\) prime spacing

If an interval contains \(H\) shell primes, Brun--Titchmarsh gives

\[
H\ll\frac{Y}{\log Y}.
\]

Since \(H=T^{2/5-o(1)}\), this forces

\[
\boxed{
Y\gg H\log H
\asymp
\frac{T^{2/5}}{L^{1/5}\Omega(T)}.
}
\tag{6.3}
\]

Consequently every interior excess prime has rank-\(H\) companions separated on both sides at the scale `(6.3)`.  Through the affine formulas, its lower digits are correspondingly far from the relevant geometric endpoints.

This places the live family in the already-identified deep bulk.  It does not create two zeros in one fixed lower characteristic: the neighboring primes are different characteristics.

## 6.3 Proved shallow plus collar

The current arbitrary-subinterval theorem says

\[
Z_p(I)\le C_A|I|^{2/3}.
\]

With threshold \(H\), put

\[
U_0:=\left(\frac{H}{C_A}\right)^{3/2}
\asymp_J
\frac{T^{3/5}}
     {L^{9/5}\Omega(T)^{3/2}}.
\tag{6.4}
\]

Choose

\[
\eta(T)\to\infty,
\qquad
\eta=o(\Omega^{3/2}),
\qquad
U_*:=U_0\eta.
\]

The proved shallow plus theorem gives, with every selected-state occurrence retained,

\[
\#\{\text{plus incidences with }t\le U_*\}
\ll_J U_*^{2/3}B_1(T).
\]

Since \(U_*^{2/3}\asymp H\eta^{2/3}\),

\[
\boxed{
U_*^{2/3}B_1(T)
\ll_J
\frac{\eta^{2/3}}{\Omega(T)}
\frac{T^2}{L^2}
=o_J(T^2/L^2).
}
\tag{6.5}
\]

Thus the shallow plus sector is genuinely deleted at the present, unamplified target scale.

## 6.4 Short geometric windows

Let \(Y_+(q,t)\) and \(Y_-(q,t)\) denote the continuous lengths `(2.12)` and `(2.8)`.  In a dyadic block \(Y\le Y_\pm<2Y\), ordinary prime interval counting gives at most

\[
O\!\left(\frac{Y}{\log(2+Y)}\right)
\]

candidate lower primes per selected state.  Therefore, summing geometrically over all blocks below

\[
Y_*:=\frac{T^{2/5}}{L^{1/5}\omega_W(T)},
\qquad
\omega_W(T)\to\infty,
\]

gives

\[
\boxed{
\#\{\text{sign incidences with }Y_\pm\le Y_*\}
\ll_J
\frac{Y_*}{L}B_1(T)
\ll_J
\frac1{\omega_W(T)}\frac{T^2}{L^2}.
}
\tag{6.6}
\]

Thus every sufficiently short fixed-state prime window is also subcritical.

Conversely, a state with more than \(R\) sign primes must have

\[
Y_\pm\gg R\log R
\asymp
\frac{T^{2/5}}{L^{1/5}\Omega(T)},
\]

so an excess state automatically lies beyond the natural short-window scale.

## 6.5 Other inherited local deletions

The shallow primary/secondary, low-product, close-shift, fixed prescribed-gap, and terminal-overlap theorems may all be imposed as additional filters in \(I_\pm^{\rm res}\).  They remove genuine subfamilies, but none proves a tail for the remaining macroscopic boxes.

For two same-sign lower primes \(p_1<p_2\), putting \(\Delta=p_2-p_1\), the exact translated lower-zero pairs are

\[
(r_2-r_1,a_{-,2}-a_{-,1})=(-6\Delta,-\Delta)
\]

in the minus cell, and

\[
(r_2-r_1,a_{+,2}-a_{+,1})=(-6\Delta,-13\Delta)
\]

in the plus cell.  A theorem for one prescribed \(\Delta\), or for a sparse set of gaps, does not sum over the unrestricted rank-\(H\) gaps generated by `(6.2)`.

After every current local deletion, boxes with

\[
t\asymp q-t\asymp p-r\asymp p-a_+\asymp T
\]

in the plus sector remain.  The plus-only obstruction can be placed entirely there.

---

# 7. Strongest current first-occurrence bound

Let

\[
K_\vee(T):=\sum_{(q,t)}\kappa_\vee(q,t).
\tag{7.1}
\]

Every counted incidence satisfies

\[
p\mid G_\vee(q,t)\mid b_n/q,
\]

hence

\[
p\mid b_n.
\]

Also \(q\mid b_n\) by `(1.3)`.  Therefore every incidence gives a two-prime common-value tuple

\[
(n;p,q).
\tag{7.2}
\]

The map

\[
((q,t),p)\longmapsto(n;p,q)
\tag{7.3}
\]

has fibre one: from \((n,q)\) recover

\[
t=n-6q.
\]

There is no extra sign multiplicity, because \(\kappa_\vee\) counts the lower prime once.

Let

\[
S_1(T):=\sum_{\ell\asymp_JT}Z(\ell)
\ll_J T^{8/5}L^{-4/5}.
\]

The banked common-integer additive-large-sieve theorem gives

\[
\boxed{
K_\vee(T)
\ll_J
\frac{S_1(T)^2}{T}
+
\frac{T^2}{L}.
}
\tag{7.4}
\]

Thus

\[
\boxed{
K_\vee(T)
\ll_J
T^{11/5}L^{-8/5}+T^2L^{-1}.
}
\tag{7.5}
\]

This remains valid for the literal untagged lcm overcarrier, because every cross term still divides \(b_n\).

Since

\[
\mathfrak X_\vee(R)\le K_\vee(T),
\]

we obtain `(V1)`.  The exact comparisons are

\[
\frac{T^{11/5}L^{-8/5}}{T^2L^{-2}}
=T^{1/5}L^{2/5},
\tag{7.6}
\]

\[
\frac{T^2L^{-1}}{T^2L^{-2}}=L.
\tag{7.7}
\]

The corrected lower ordered-zero-pair moment gives only

\[
K_\vee(T)\ll_J T^{11/5}L^{-3/5}
\]

for the exact signed occurrence subfamily, which is weaker than `(7.5)`.

The pointwise carrier height and prime-interval bounds give

\[
K_\vee(T)
\ll_J
B_1(T)\frac{T}{L}
\ll_J
T^{13/5}L^{-9/5},
\]

which is weaker still.

Hence `(7.5)` is the strongest current unconditional full-bulk ledger.

---

# 8. Why the \(2R\) subtraction does not improve the first moment

Let

\[
\mathscr U_{\rm high}:=\{(q,t):\kappa_\vee(q,t)>2R\}.
\]

Then exactly

\[
\mathfrak X_\vee(R)
=
\sum_{\sigma\in\mathscr U_{\rm high}}\kappa_\vee(\sigma)
-2R|\mathscr U_{\rm high}|.
\tag{8.1}
\]

The current inputs give no useful lower bound for the second term.  Its largest possible total over all selected states is

\[
2RB_1(T)
=o(T^2/L^2)
\]

by `(V4)`.  Therefore subtracting the full quota from every possible state is itself only a subcritical correction.  It cannot turn a potentially supercritical first occurrence mass `(7.5)` into the desired bound.

This is the exact point where a tail estimate, second factorial correlation, or radical-height distribution theorem becomes necessary.

---

# 9. Raw radical height and its exact endpoint

Remove the upper prime and the fixed exceptional set from \(G_\vee\), and restrict its radical to the exact retained lower-prime window.  This is \(\mathcal R_\vee(q,t)\) in `(V5)`.

Since every retained lower prime satisfies

\[
\log p\asymp_JL,
\]

one has statewise

\[
\log\mathcal R_\vee(q,t)
asymp_J L\kappa_\vee(q,t).
\tag{9.1}
\]

Also

\[
\log\mathcal R_\vee(q,t)
\le\log G_\vee(q,t)
\le\log(b_n/q)
\ll_JT.
\tag{9.2}
\]

Thus

\[
\kappa_\vee(q,t)\ll_J T/L.
\tag{9.3}
\]

Summing `(9.1)` and using `(7.5)` gives the strongest aggregate radical first moment

\[
\boxed{
\sum_{(q,t)}\log\mathcal R_\vee(q,t)
\ll_J
T^{11/5}L^{-3/5}+T^2.
}
\tag{9.4}
\]

The radical target corresponding to the desired excess is \(T^2/L\).  The two terms in `(9.4)` miss it by exactly the same factors `(7.6)--(7.7)`.

The threshold subtraction contributes at most

\[
2RLB_1(T)
\ll_J
\frac{T^2}{L\Omega(T)}
=o(T^2/L).
\tag{9.5}
\]

Hence total radical height plus the known quota still cannot prove a positive-part tail.

This proves that `(VEE-RADICAL-EXCESS)` is a genuinely new distribution statement, not a consequence of total height.

---

# 10. Second factorial route and its exact exponent

Put

\[
F_\vee(T)
:=
\sum_{(q,t)}
\kappa_\vee(q,t)
\bigl(\kappa_\vee(q,t)-1\bigr).
\tag{10.1}
\]

For every integer \(k\ge0\) and \(H\ge1\),

\[
(k-2H)_+
\le
\frac{k(k-1)}{2H}.
\tag{10.2}
\]

Therefore

\[
\mathfrak X_\vee(R)
\le
\frac{F_\vee(T)}{2H}
+o(T^2/L^2).
\tag{10.3}
\]

A sufficient second factorial theorem is thus

\[
\boxed{
F_\vee(T)
=o_J\!\left(
R\frac{T^2}{L^2}
\right)
=
o_J\!\left(
\frac{T^{12/5}}
     {\Omega(T)L^{16/5}}
\right).
}
\tag{10.4}
\]

No such theorem is currently banked.

Using only the pointwise bound `(9.3)` and the one-leaf estimate `(7.5)` gives

\[
F_\vee(T)
\ll_J
\frac TL K_\vee(T),
\]

hence

\[
\boxed{
F_\vee(T)
\ll_J
T^{16/5}L^{-13/5}
+
T^3L^{-2}.
}
\tag{10.5}
\]

The ratios of the two terms in `(10.5)` to the required scale `(10.4)` are

\[
\boxed{
\Omega(T)T^{4/5}L^{3/5}
}
\tag{10.6}
\]

and

\[
\boxed{
\Omega(T)T^{3/5}L^{6/5},
}
\tag{10.7}
\]

respectively.  The factorial route therefore does not close without a new fixed-state codegree theorem.

---

# 11. A concrete narrower sufficient theorem

Let

\[
\mathscr A_\vee(T)
:=
\{\ell=((q,t),p):p\in\mathcal P_\vee(q,t)\}
\tag{11.1}
\]

with every exact state, lower prime, sign/window tag, and inherited residual flag retained.  Define

\[
D_\vee(\ell)
:=
\#\{p'\ne p:((q,t),p')\in\mathscr A_\vee(T)\}.
\tag{11.2}
\]

Finite interchange gives

\[
\boxed{
F_\vee(T)
=
\sum_{\ell\in\mathscr A_\vee(T)}D_\vee(\ell).
}
\tag{11.3}
\]

The common-\(n\) theorem gives

\[
|\mathscr A_\vee(T)|
\ll_J
T^{11/5}L^{-8/5}+T^2L^{-1}.
\tag{11.4}
\]

The critical average partner scale is

\[
D_*(T)
:=
\frac{R T^2/L^2}{T^{11/5}L^{-8/5}}
\asymp
\boxed{
\frac{T^{1/5}}
     {\Omega(T)L^{8/5}}.
}
\tag{11.5}
\]

Therefore `(VEE-AVG-CODEGREE)` implies

\[
F_\vee(T)
=o_J(RT^2/L^2),
\]

and then `(10.3)` proves the desired excess estimate.

This theorem is narrower than a full CSRE statement:

```text
- the upper occurrence (q,t) is fixed inside each fibre;
- only foreign primes in the exact residual plus-or-minus window count;
- the common fixed value n=6q+t is retained;
- the two sign-specific secondary gcd conditions are retained;
- all shallow, close-shift, overlap, and product filters remain active;
- no fourth C4 edge or arbitrary upper-state pair is counted.
```

It is exactly the missing cross-characteristic large-prime-divisor codegree.

---

# 12. Why current C4 packing information does not prove the tail

The newest fixed-state carrier theorem gives a useful squarefree product for selected cycles.  In the present reduced problem, \(\kappa_\vee\) already counts a squarefree set of distinct lower primes.  Thus that theorem introduces no additional exponent.

For different lower primes

\[
p_1\ne p_2,
\]

the corresponding support edges lie in different lower-characteristic graphs.  Edge-disjointness is automatic across those graphs.  A maximum packing can therefore use one private C4 for each legal divisor of one fixed-state carrier.

The other three edges of the cycles may have upper states depending on the cycle.  Complete fibres do not force them into a second fixed integer attached only to `(q,t)`.  Consequently there is no implication such as

\[
p^2\mid G_\vee(q,t)
\]

or

\[
p\mid H(q,t)
\]

for a second fixed carrier of smaller aggregate height.

This is an arithmetic obstruction, not merely a graph-theoretic one: the only fixed characteristic-zero object currently common to all such lower primes is `(V3)` itself.

---

# 13. Support-level saturation of all current ledgers

The following is an abstract support/incidence model.  It is not asserted to be realized by the actual Apéry sequence.  Its role is precise: it shows that the current moment, height, local-collar, and common-value inequalities do not logically imply `(VEE-RADICAL-EXCESS)`.

Choose

\[
P_0\asymp\frac TL,
\qquad
\rho\asymp T^{3/5}L^{1/5},
\qquad
m_*\asymp T^{3/5}L^{-4/5},
\qquad
S_*\asymp T^{8/5}L^{-4/5}.
\tag{13.1}
\]

## 13.1 Lower zero masks

Give each of \(P_0\) lower-prime labels a separated \(2/3\)-Frostman zero mask of size \(\rho\), so that

\[
|Z_p\cap I|\ll |I|^{2/3}+1
\]

for every interval.  Then

\[
P_0\rho
\asymp
T^{8/5}L^{-4/5},
\tag{13.2}
\]

\[
P_0\rho^2
\asymp
T^{11/5}L^{-3/5}.
\tag{13.3}
\]

Thus the corrected first and second zero moments and the arbitrary-subinterval bound are respected.

## 13.2 Selected states and plus-only fibres

Choose \(S_*\) selected state labels in interior macroscopic plus boxes and give each state \(m_*\) distinct lower-prime neighbors.  Take the minus support empty.

The one-leaf mass is

\[
S_*m_*
\asymp
T^{11/5}L^{-8/5},
\tag{13.4}
\]

which saturates the leading common-\(n\) scale.

The fixed-state radical cost is

\[
m_*L
\asymp
T^{3/5}L^{1/5}
=o(T),
\tag{13.5}
\]

so it is compatible with the pointwise \(O(T)\) Apéry coefficient height.  Since \(m_*\le\rho\), the lower zero masks have ample ordered-pair capacity.

Place every state with

\[
t\asymp q-t\asymp T
\]

and every lower prime a fixed positive proportion away from both window endpoints.  All shallow and short-window collars are avoided.  The prime neighbors may be spread so that no prescribed sparse set of gaps carries a significant fraction.

Private C4 completions may be attached in the separate lower-characteristic graphs, exactly as in the fixed-state packing obstruction.  This adds no support reuse.

## 13.3 Excess

One has

\[
\frac{m_*}{R}
\asymp
\Omega(T)T^{1/5}L^{2/5}
\longrightarrow\infty.
\tag{13.6}
\]

Therefore

\[
\begin{aligned}
\sum_\sigma(m_*-2R)_+
&\sim S_*m_*\\
&\asymp T^{11/5}L^{-8/5}.
\end{aligned}
\tag{13.7}
\]

Meanwhile the total reserved quota is

\[
2RS_*
\asymp
\frac{T^2}{L^2\Omega(T)}.
\tag{13.8}
\]

The ratio of `(13.7)` to the desired target is

\[
\boxed{T^{1/5}L^{2/5}.}
\tag{13.9}
\]

The model simultaneously respects

```text
selected-state count:              S_* ~ T^(8/5)L^(-4/5)
lower first moment:                P_0 rho ~ T^(8/5)L^(-4/5)
lower second factorial moment:     P_0 rho^2 ~ T^(11/5)L^(-3/5)
arbitrary-subinterval cap:         |Z_p cap I| << |I|^(2/3)
common-n one-leaf bound:           S_* m_* ~ T^(11/5)L^(-8/5)
fixed-state carrier height:        m_* L = o(T)
pointwise prime capacity:           m_* << T/L
free quota:                         2R S_* = o(T^2/L^2)
shallow/short collars:              avoided by interior placement
sign overlap:                       empty; plus-only
occurrence fibre:                   one per ((q,t),p)
```

Again, this is not an actual-Apéry counterexample.  It proves that a new arithmetic theorem about the distribution of legal large prime divisors in the fixed gcds is indispensable.

---

# 14. Final ledger

| Input or route | Exact consequence | Status versus \(T^2/L^2\) |
|---|---:|---:|
| Free quota | \(2RB_1\ll T^2/(L^2\Omega)\) | subcritical |
| Shallow plus collar | \(\ll (\eta^{2/3}/\Omega)T^2/L^2\) | subcritical |
| Short sign window | \(\ll T^2/(L^2\omega_W)\) | subcritical |
| Corrected lower pair moment | \(T^{11/5}L^{-3/5}\) | misses by \(T^{1/5}L^{7/5}\) |
| Fixed carrier height + \(B_1\) | \(T^{13/5}L^{-9/5}\) | misses by \(T^{3/5}L^{1/5}\) |
| Common-\(n\) first moment | \(T^{11/5}L^{-8/5}+T^2/L\) | misses by \(T^{1/5}L^{2/5}\) and \(L\) |
| Current factorial bound | \(T^{16/5}L^{-13/5}+T^3/L^2\) | power-supercritical |
| Plus/minus intersection | controls overlap only | no union bound |
| Reflection intersection | controls overlap only | no union bound |
| Fixed prescribed gaps | one/sparse gap set | cannot sum live gaps |
| Edge-disjoint maximum C4 packing | squarefree selected product | already encoded by radical |
| Required new coefficient theorem | `VEE-RADICAL-EXCESS` | sufficient |
| Required reusable correlation | `VEE-AVG-CODEGREE` | sufficient |

# Bottom line

Even after the exact reduction

\[
\delta_{C4,LU}
\le
\frac12\mathfrak X_\vee(R),
\]

the current arithmetic does not prove

\[
\mathfrak X_\vee(R)=o_J(T^2/L^2).
\]

The local deletion package successfully removes the free quota, shallow plus endpoints, short fixed-state windows, prescribed close shifts, and other boundary sectors.  What remains is an interior fixed-upper-state large-prime divisor tail.  All counted primes divide one fixed coefficient \(b_{6q+t}/q\) and at least one of the two companion coefficients, but the existing theory controls only the total one-leaf incidence, not its concentration in a single state.

The first exact missing estimate is `VEE-RADICAL-EXCESS`; the most reusable narrower correlation theorem is `VEE-AVG-CODEGREE`.  Either would close the proposed Hall-deficiency reduction with every state, prime, sign overlap, and physical occurrence counted correctly.