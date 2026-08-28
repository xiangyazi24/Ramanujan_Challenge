ANSWER Q4980 6ec4acfc

# Q4980 — cold source-consumer audit of the PS14 literal support-size deletion

## Verdict

**The three advertised threshold exponents are mathematically correct, but Q4958's proposed Mellin-`L^2` → endpoint-operator bridge is not a valid source inference.**  The first false/unproved line is the promotion of the banked `N^{161/100}` **actual-row endpoint** to a full high-character coefficient-space operator inequality of the form

```latex
\|T_{N,\omega}x\|_2^2
\ll N^{241/100+o(1)}\|x\|_2^2
\qquad\text{for arbitrary support-projected }x,
\tag{FALSE-OP}
```

or equivalently

```latex
|\langle x,\Gamma_{N,\omega}y\rangle|
\ll N^{241/100+o(1)}\|x\|_2\|y\|_2.
```

The canonical Oracle A source explicitly distinguishes the needed moving-reciprocal estimate as a **scalar Rayleigh bound for the actual Apéry vector** and says that it “neither asks for a full operator norm” (`problems/3.2/oracleA_result.tex`, immediately after `eq:oracleA-sdc`).  Its completed reciprocal endpoint is obtained by the source's fixed-`d` Cauchy/finite-Parseval calculation, not by a uniform operator theorem on arbitrary Mellin coefficient vectors.

There is also a normalization trap.  The canonical additive coefficients are

```latex
c_{p,a}=\frac{(\log p)F_p(a)}p,
```

and therefore the raw additive coefficient square is

```latex
Q_P=\sum_p(\log p)^2
\left(\frac{Z(p)}p-\frac{Z(p)^2}{p^2}\right).
```

Thus the restricted raw additive `L^2` mass is of order **`N E_\alpha`**, not `E_\alpha`.  The estimate

```latex
E_\alpha\ll \min\{N^{-2/5},N^{\alpha-1}\}N^{o(1)}
```

is correct only when `E_\alpha` means the **Mellin/high-character coefficient mass**

```latex
E_\alpha
:=\sum_{p\in\mathcal S_\alpha}|u_p|^2C_p^H,
\qquad
u_p=\frac{\log p}{p}v(p/N),
```

where `C_p^H` is the squared coefficient mass of the high-character projection.  It is not the restricted `Q_P` from the first displayed normalization in `oracleA_result.tex`.

**However, the actual deletion itself is salvageable, with exactly the exponents claimed in Q4958, without any `M_4` bridge and without `(FALSE-OP)`.**  One reruns the source's fixed-`d` endpoint after projecting the *prime supports*.  If

```latex
L_\alpha:=\sum_{p\in\mathcal S_\alpha}|u_p|C_p^H,
```

then

```latex
L_\alpha\ll N^{\alpha+o(1)},
\qquad
L:=\sum_{p\asymp N}|u_p|C_p^H\ll N^{3/5+o(1)}.
```

For a row whose two defining characteristics are restricted to prime sets `X,Y`, the same Cauchy + finite Parseval argument gives

```latex
\sum_{0<|d|\asymp D}|B_{X,Y}(d)|^2
\ll D\,L_XL_YN^{o(1)},
\qquad D=N^{41/100}.
\tag{R}
```

Consequently a Gram block in which `k` of the four characteristic **occurrences** lie in `\mathcal S_\alpha` satisfies

```latex
\boxed{
|\mathfrak G_k(\alpha)|
\ll
N^{161/100+(k/2)(\alpha-3/5)+o(1)}.}
\tag{Gk}
```

This proves exactly:

```text
k >= 1, alpha = 14/25-delta:  <= N^(159/100-delta/2+o(1));
k = 2,  alpha = 29/50-delta:  <= N^(159/100-delta+o(1));
k = 4,  alpha = 59/100-delta:  <= N^(159/100-2delta+o(1)).
```

So the correct classification is

```text
Q4958 E_alpha Mellin coefficient bound:                  VALID
Q4958 pair-tensor L2 bookkeeping:                        VALID
Q4958 abstract endpoint operator promotion:              INVALID / NOT A SOURCE THEOREM
Q4958 numerical Gram support exponent (Gk):              VALID BY A DIFFERENT PROOF
any-small threshold alpha=14/25-delta:                   VALID
one-small-source-each-side threshold 29/50-delta:        VALID
all-four-small threshold 59/100-delta:                   VALID
extra D^2 at this d-shell stage:                          FALSE; only D occurs
one-prime M4/Q4928 bridge used:                           NO
```

The corrected insertion-ready TeX is in §10 below.  It should be inserted **after the current PS14 mask-completion / reciprocal-row endpoint**, not upstream of the frequency-mixing masks.

---

## 0. Source pin and visibility boundary

I used only `xiangyazi24/Ramanujan_Challenge`, `problems/3.2`.

The connector-visible canonical proof head is

```text
main = c5d932b66ce5e4f1657b587215d290ae7a13018b
```

and the current connector-visible canonical Oracle A file is

```text
problems/3.2/oracleA_result.tex
blob f0bfccc441bda22c17658c5586fbb2e6c6431238
```

The delivery branch moved during this audit; immediately before writing this answer it was

```text
chatgpt-drop = dc98d16d8c85e7b43cb5e629402fe5ccf1c9a5f1
```

with message `Add Q4975 actual-orbit core audit`.

The connector-visible `main` file predates the literal local `PS14` labels, so I do not invent line-local PS14 quotations.  The source facts used below are rederived from the canonical equations that are visible:

1. `oracleA_result.tex`, `ssec:oracleA-result`:

   ```latex
   c_{p,a}=\frac{(\log p)F_p(a)}p.
   ```

2. `eq:oracleA-Q-parseval`:

   ```latex
   Q_P
   =\sum_p(\log p)^2
     \left(\frac{Z(p)}p-\frac{Z(p)^2}{p^2}\right).
   ```

3. `eq:oracleA-cross-bilinear` and `eq:oracleA-moving-reciprocal` retain the literal pair coefficient

   ```latex
   \frac{(\log p)(\log q)}{pq}F_p(\cdots)F_q(\cdots),
   ```

   so in the smooth fixed-ratio block one has exactly the later normalization

   ```latex
   u_p=\frac{\log p}{p}v(p/N),
   \qquad |u_p|\ll N^{-1+o(1)}.
   ```

4. Immediately after `eq:oracleA-sdc`, the source explicitly says the desired moving reciprocal estimate is a scalar Rayleigh bound for the actual Apéry vector and **not** a full operator-norm assertion.

5. `atom_tail_section.tex` records the unconditional average

   ```latex
   \sum_{p\le x} Z(p)\ll x^{8/5}/\log x.
   ```

   Hence on a shell `p\asymp N`,

   ```latex
   \sum_{p\asymp N}|u_p|Z(p)\ll N^{3/5+o(1)},
   \qquad
   \sum_{p\asymp N}|u_p|^2Z(p)\ll N^{-2/5+o(1)}.
   ```

No Q4928 one-prime fourth-moment bridge, no numerical scan, and no external project is used.

---

# 1. `u_p` normalization: the factor of `p` that must not be lost

The source additive coefficient is

```latex
c_{p,a}=u_pF_p(a)
```

(up to the bounded smooth cutoff `v(p/N)`).  Finite additive Parseval gives

```latex
\sum_{a\bmod p}|F_p(a)|^2=pZ(p),
```

and on nonzero frequencies

```latex
\sum_{a\ne0}|F_p(a)|^2=pZ(p)-Z(p)^2.
```

Therefore

```latex
\sum_{a\ne0}|c_{p,a}|^2
=|u_p|^2\bigl(pZ(p)-Z(p)^2\bigr).
\tag{1.1}
```

This is exactly `eq:oracleA-Q-parseval`.  In particular the raw additive one-prime `L^2` mass contains **one factor `p`**.

After passing to multiplicative characters on `\mathbb F_p^\times`, write the high-character projection as

```latex
H_p^H(a)=\sum_{\chi\in\mathcal E_p}\gamma_p(\chi)\chi(a),
\qquad a\ne0,
```

and put

```latex
C_p^H:=\sum_{\chi\in\mathcal E_p}|\gamma_p(\chi)|^2.
```

Because this is an orthogonal projection of the nonprincipal/centered character expansion,

```latex
C_p^H\le C_p\ll Z(p).
\tag{1.2}
```

Equivalently, multiplicative Parseval gives

```latex
\sum_{a\ne0}|H_p^H(a)|^2=(p-1)C_p^H.
\tag{1.3}
```

The coefficient mass in the **unnormalized character basis** is therefore

```latex
E:=\sum_p|u_p|^2C_p^H,
\tag{1.4}
```

whereas the corresponding additive evaluation-vector square has size `\asymp N E`.  That factor is carried by the character-evaluation map.  It cannot be discarded and then recovered for free by calling the endpoint an operator norm.

This is the first normalization check on Q4958.

---

# 2. The small-characteristic Mellin coefficient mass is correct

For `0\le\alpha\le2/3`, put

```latex
\mathcal S_\alpha
:=\{p\asymp N:Z(p)\le N^\alpha\},
```

and define precisely

```latex
E_\alpha
:=\sum_{p\in\mathcal S_\alpha}|u_p|^2C_p^H.
\tag{2.1}
```

By (1.2) and the global shell energy,

```latex
E_\alpha\le E
\ll N^{-2/5+o(1)}.
\tag{2.2}
```

On the other hand, `|u_p|\ll N^{-1+o(1)}`, `C_p^H\ll Z(p)\le N^\alpha`, and there are `N^{1+o(1)}` shell primes, so

```latex
E_\alpha
\ll N^{1+o(1)}N^{-2}N^\alpha
=N^{\alpha-1+o(1)}.
\tag{2.3}
```

Hence

```latex
\boxed{
E_\alpha
\ll
\min\{N^{-2/5},N^{\alpha-1}\}N^{o(1)}.}
\tag{2.4}
```

So the **first displayed numerical estimate in the dispatch is valid**, provided `E_\alpha` is defined by (2.1).

If Q4958 instead identifies `E_\alpha` with the restriction of the canonical raw additive `Q_P`, the line is false by a factor `N`: from (1.1),

```latex
Q_{P,\alpha}
\asymp N E_\alpha
```

at power precision.

No fourth moment appears anywhere in (2.1)--(2.4).

---

# 3. High-character projection and pair tensor norm

The high-character projection is harmless for the support-size bookkeeping because it is orthogonal characteristic by characteristic:

```latex
0\le C_p^H\le C_p\ll Z(p).
```

Prime-support projection and character projection commute: the first acts on the `p` label, the second on the `\chi` label within that characteristic.

For prime sets `X,Y`, let the pair coefficient tensor have entries, up to unit phases,

```latex
A_{p,q,\chi,\psi}
=1_X(p)1_Y(q)1_{p\ne q}
 u_pu_q\gamma_p(\chi)\gamma_q(\psi).
```

Then exactly by positivity,

```latex
\begin{aligned}
\|a_{X,Y}^{HH}\|_2^2
&=\sum_{\substack{p\in X,q\in Y\\p\ne q}}
 |u_p|^2|u_q|^2 C_p^HC_q^H\\
&\le
\left(\sum_{p\in X}|u_p|^2C_p^H\right)
\left(\sum_{q\in Y}|u_q|^2C_q^H\right).
\end{aligned}
\tag{3.1}
```

Thus, with `X=\mathcal S_\alpha` or the full shell,

```latex
\|a_{X,Y}^{HH}\|_2^2\le E_XE_Y.
\tag{3.2}
```

This part of Q4958 is sound.  The restriction `p\ne q` only removes nonnegative terms from the coefficient norm.

What does **not** follow is that the output Gram is obtained by multiplying (3.2) by a source-proved uniform `2\to2` norm of size `N^{241/100}`.  The source endpoint is an actual-row finite-Parseval estimate.  That is the first invalid promotion.

---

# 4. First false line in the proposed proof

The first load-bearing false/unproved statement is any line equivalent to

```latex
\|T_{N,\omega}P a\|_2^2
\ll N^{241/100+o(1)}\|Pa\|_2^2
\tag{4.1-FALSE}
```

for an arbitrary prime-support projection `P` of the high-character coefficient space.

Why this is not licensed by the source:

1. The visible canonical Oracle A source itself calls its missing reciprocal estimate a **scalar Rayleigh bound for the actual Apéry vector** and explicitly says it is not a full operator norm.

2. The banked `161/100` endpoint has the form

   ```latex
   \sum_{0<|d|\asymp D}|B_d|^2
   \ll D B_N^2N^{o(1)},
   \qquad
   B_N:=\sum_p|u_p|Z(p),
   \tag{4.2}
   ```

   with `D=N^{41/100}` and `B_N\ll N^{3/5+o(1)}`.  It is obtained by Cauchy and two finite Parsevals on the **actual zero transforms**.

3. Rewriting the numerical exponent as

   ```latex
   N^{161/100}
   =N^{241/100}\,N^{-4/5}
   ```

   does not turn (4.2) into

   ```latex
   \|T\|_{2\to2}^2\ll N^{241/100}.
   ```

   The `N^{-4/5}` is the actual pair-tensor coefficient norm; the source proof of (4.2) uses the stronger row structure and does not factor through a full-space operator norm.

4. The additive and Mellin coordinate normalizations differ by the `p-1` character-basis norm in (1.3).  Applying an additive endpoint operator directly to the Mellin coefficient norm (2.1) silently drops these factors.

Therefore the implication

```text
small E_alpha  +  banked endpoint
        => projected Gram via a universal operator norm
```

is not a valid source-consumer step.

This does **not** refute the desired projected Gram bound itself.  It only refutes Q4958's proposed bridge to it.

---

# 5. Correct source-native replacement: support-projected row endpoint

Define the first-mass quantity

```latex
L_X:=\sum_{p\in X}|u_p|C_p^H.
\tag{5.1}
```

For the full shell, (1.2) and the canonical average zero bound give

```latex
L:=L_{\mathcal P_N}
\ll\sum_{p\asymp N}|u_p|Z(p)
\ll N^{3/5+o(1)}.
\tag{5.2}
```

For `X=\mathcal S_\alpha`,

```latex
L_\alpha
\ll N^{1+o(1)}N^{-1}N^\alpha
=N^{\alpha+o(1)}.
\tag{5.3}
```

Now take one completed PS14 high-high row block with prime supports `X,Y`:

```latex
B_{X,Y}(d)
=\sum_{\substack{p\in X,q\in Y\\p\ne q}}
 u_pu_q\,\omega_{p,q}(d)
 H_p^H(-d\bar q_p)H_q^H(-d\bar p_q),
\tag{5.4}
```

where all remaining smooth/orientation factors are included in `\omega` and satisfy

```latex
|\omega_{p,q}(d)|\ll N^{o(1)}.
```

For fixed nonzero `d` with `|d|\asymp D<N`, weighted Cauchy gives

```latex
|B_{X,Y}(d)|^2
\le N^{o(1)}\mathcal P_X(d)\mathcal P_Y(d),
\tag{5.5}
```

with

```latex
\mathcal P_X(d)
:=\sum_{p\in X}|u_p|
  \sum_{q\in Y}|u_q|
  |H_p^H(-d\bar q_p)|^2,
```

and the symmetric expression for `\mathcal P_Y`.

Fix `p`.  Since all shell primes lie in fixed-ratio intervals of length `O(N)`, reduction modulo `p\asymp N` has bounded multiplicity.  Because `0<|d|<p`, the map

```latex
q\mapsto -d\bar q_p\pmod p
```

hits only nonzero residues and has `O(1)` multiplicity.  Therefore, using

```latex
u_*:=\max_{q\asymp N}|u_q|\ll N^{-1+o(1)}
```

and multiplicative Parseval (1.3),

```latex
\begin{aligned}
\sum_{q\in Y}|u_q|
 |H_p^H(-d\bar q_p)|^2
&\ll N^{-1+o(1)}
 \sum_{a\ne0}|H_p^H(a)|^2\\
&=N^{-1+o(1)}(p-1)C_p^H\\
&\ll C_p^HN^{o(1)}.
\end{aligned}
\tag{5.6}
```

Hence

```latex
\mathcal P_X(d)\ll L_XN^{o(1)},
\qquad
\mathcal P_Y(d)\ll L_YN^{o(1)},
```

and

```latex
\boxed{
|B_{X,Y}(d)|^2
\ll L_XL_YN^{o(1)}.}
\tag{5.7}
```

There are `N^{41/100+o(1)}` integers in the `d` shell.  Thus

```latex
\boxed{
\sum_{0<|d|\asymp D}|B_{X,Y}(d)|^2
\ll D L_XL_YN^{o(1)}.}
\tag{5.8}
```

This is exactly the banked endpoint argument, rerun after the prime-support projection.  It is source-native and does not require a full coefficient-space operator theorem.

---

# 6. Gram block with `k` small characteristic occurrences

Split each of the two characteristic slots in each row according to

```latex
1=1_{\mathcal S_\alpha}+1_{\mathcal S_\alpha^c}.
```

Let `\sigma=(\sigma_1,\sigma_2)` and `\tau=(\tau_1,\tau_2)` be two row signatures, each slot being small or unrestricted/large, and let `k=k(\sigma,\tau)` be the number of the four **occurrences** constrained to `\mathcal S_\alpha`.

By Cauchy in the `d` shell and (5.8),

```latex
\begin{aligned}
|\langle B_\sigma,B_\tau\rangle|
&\le \|B_\sigma\|_2\|B_\tau\|_2\\
&\ll
D\,L_\alpha^{k/2}L^{(4-k)/2}N^{o(1)}.
\end{aligned}
\tag{6.1}
```

Using (5.2)--(5.3),

```latex
\begin{aligned}
|\langle B_\sigma,B_\tau\rangle|
&\ll
N^{41/100}
N^{k\alpha/2}
N^{(4-k)3/10}
N^{o(1)}\\
&=
\boxed{
N^{161/100+(k/2)(\alpha-3/5)+o(1)}}.
\end{aligned}
\tag{6.2}
```

This is exactly the exponent asserted in the dispatch.

Notice that the proof did not use (4.1-FALSE).  The Mellin coefficient bound `E_\alpha` and pair-tensor norm remain correct bookkeeping facts, but the output estimate is furnished by the actual-row endpoint (5.8).

---

# 7. Exact threshold ledger

The critical shell exponent is `159/100`.

## 7.1 Any block containing at least one small characteristic

Take

```latex
\alpha=\frac{14}{25}-\delta
=\frac{56}{100}-\delta.
```

For `k\ge1`, the worst case is `k=1`.  Then

```latex
\begin{aligned}
\frac{161}{100}
+\frac12\left(\alpha-\frac35\right)
&=\frac{161}{100}
+\frac12\left(-\frac4{100}-\delta\right)\\
&=\boxed{\frac{159}{100}-\frac\delta2}.
\end{aligned}
```

Hence every support block with at least one occurrence in

```latex
Z(p)\le N^{14/25-\delta}
```

is strictly subcritical.

## 7.2 Projected self-shell with one small source prime on each side

A self-shell of a row with exactly one small characteristic has two small **occurrences** after taking the square, so `k=2`.  With

```latex
\alpha=\frac{29}{50}-\delta
=\frac{58}{100}-\delta,
```

one gets

```latex
\frac{161}{100}
+\left(\alpha-\frac35\right)
=
\boxed{\frac{159}{100}-\delta}.
```

Equivalently, directly from (5.8),

```latex
D L_\alpha L
\ll
N^{41/100+29/50+3/5-\delta+o(1)}
=N^{159/100-\delta+o(1)}.
```

If the same prime label occurs on both sides of the positive self-shell diagonal, this does not create a multiplicity loss: the diagonal is already included in the positive norm used in (5.8).  In the literal four-distinct cross term such equal-characteristic diagonals are removed by the pre-existing prime-intersection decomposition; see §8.

## 7.3 All four occurrences small

For a row with both characteristics in `\mathcal S_\alpha`, its self-shell has `k=4`.  Take

```latex
\alpha=\frac{59}{100}-\delta.
```

Then

```latex
\frac{161}{100}
+2\left(\alpha-\frac35\right)
=
\boxed{\frac{159}{100}-2\delta}.
```

Equivalently

```latex
D L_\alpha^2
\ll
N^{41/100+118/100-2\delta+o(1)}.
```

Thus all three proposed cutoffs are arithmetically correct.

---

# 8. Mask audit: what is safe and what is not

The support estimate must be inserted at the **completed PS14 reciprocal-row stage**.  The following distinctions are load-bearing.

### 8.1 Prime-support and high-character projections

Safe.  Prime support acts on `p`; high-character projection acts on `\chi`.  They commute, and

```latex
C_p^H\le C_p.
```

No principal/small-character mass is reintroduced.

### 8.2 Smooth factors and pair twists

Safe.  The proof of (5.7) uses only

```latex
|\omega_{p,q}(d)|\ll N^{o(1)}.
```

No factorization of `\omega` in `p` and `q` is required.

### 8.3 High-centered, near/reflected-near, alias, orientation, and circular masks

Do **not** delete these termwise from the original signed form by monotonicity.  Some of them depend on the row zero/shadow coordinate and therefore mix the finite Fourier variable before completion.

The correct consumer order is the one already used by the current PS14 endpoint:

```text
exact partition/completion of closed masks
    -> completed reciprocal row
    -> high-character projection
    -> prime-support projection
    -> endpoint (5.8).
```

Support projection is characteristic-diagonal, so it commutes with the exact finite partitions.  The previously banked complementary sectors remain banked; restricting their prime support cannot enlarge the positive/absolute estimates used for those sectors.

### 8.4 Literal all-four-distinct mask

This is the most important sign warning.  One must **not** claim

```latex
|\text{four-distinct sub-sum}|
\le |\text{unmasked Gram}|
```

by deleting terms.  That is false for a signed Gram.

Use the existing exact prime-intersection inclusion--exclusion:

```text
literal four-distinct
 = completed projected Gram
   - same-row / same-P pieces
   - one-shared-prime pieces
   + finite intersection corrections.
```

The first term is bounded by (6.2).  The non-four-distinct pieces are precisely the strata closed before the residual PS14 high-high form.  The support split does not create a new physical stratum.  Thus the support deletion is valid for the **literal** four-distinct consumer only after this exact decomposition.

If Q4958 instead removed the four-distinct mask by absolute monotonicity, that would be a separate invalid line.  It occurs after the earlier operator-promotion defect identified in §4.

### 8.5 A genuinely uncompleted affine-level selector

If one retains an extra selector depending on the Fourier/shadow variable that has **not** already been completed or banked, then (5.4) is no longer the correct factorized row and the insertion cannot be applied verbatim.  The requested insertion is safe specifically after the current global PS14 completion.

---

# 9. Midpoint, diamond, `D` versus `D^2`, and cross terms

## 9.1 Midpoint corrections

No argument here assumes that high-centeredness eliminates the reflection-fixed midpoint zero.

With

```latex
U_N:=\sum_{p\asymp N}|u_p|\ll N^{o(1)},
\qquad
B_N:=\sum_p|u_p|Z(p)\ll N^{3/5+o(1)},
```

the standard midpoint expansion has, at the `d`-shell level,

```latex
\|B_{00}\|_2^2 \ll D B_N^2          \ll N^{161/100+o(1)},
\|B_{10}\|_2^2+\|B_{01}\|_2^2
                  \ll D U_N^2 B_N    \ll N^{101/100+o(1)},
\|B_{11}\|_2^2 \ll D U_N^4          \ll N^{41/100+o(1)}.
```

Therefore the largest midpoint cross term is

```latex
|\langle B_{00},B_{10}+B_{01}\rangle|
\ll N^{131/100+o(1)},
```

and all further midpoint crosses are smaller.  These are far below `159/100`.  Prime-support restriction only decreases the positive masses entering these estimates.

Hence midpoint terms do not change any of the support thresholds.

## 9.2 Diamond correction

The reflection-compatible endpoint correction has one-prime weighted size

```latex
\mathcal R_N
\ll N^{-2/5+o(1)}.
```

Replacing a canonical centered factor by the diamond factor changes one fixed-`d` two-prime row by at most

```latex
N^{1/5+o(1)},
```

so the correction shell square is

```latex
\ll D N^{2/5+o(1)}
=N^{81/100+o(1)}.
```

Its cross term with the `161/100` main shell is at most

```latex
N^{(161/100+81/100)/2+o(1)}
=N^{121/100+o(1)}.
```

Again this is safely below the target, and support restriction can only reduce the absolute one-prime correction mass.  No central-zero exclusion is needed.

## 9.3 Exactly one `D`, not `D^2`

At the stage of (5.8) the quantity is already the discrete shell norm

```latex
\sum_{0<|d|\asymp D}|B_d|^2.
```

There are `D N^{o(1)}` frequencies, so the fixed-`d` bound is multiplied by **one `D`**:

```latex
D=N^{41/100}.
```

The `D^2` factor that appears when recovering a shell from a Fejér/short-arc integral belongs to the **earlier transform/restoration step**.  It is not multiplied again after the `d`-shell endpoint has been formed.

Paying `D^2` here would double-count the restoration and destroy the `161/100` normalization.  Therefore Q4958 is correct only with one `D` at this stage.

## 9.4 Cross terms among support pieces

The prime split produces only finitely many row pieces.  Every cross term is handled by

```latex
|\langle B_1,B_2\rangle|
\le \|B_1\|_2\|B_2\|_2,
```

which is exactly what yields the factor `k/2` in (6.2).  No sign assumption and no orthogonality between different prime-support blocks is needed.

The already banked mask/midpoint/diamond errors have shell exponents at most the values above; crossing them with a support-projected main term remains subcritical by the same Cauchy estimate.

---

# 10. Corrected insertion-ready TeX for placement after PS14

The following is the safe replacement for the Q4958 insertion.  It deliberately records `E_\alpha` as a coefficient bookkeeping lemma but proves the Gram deletion through the support-projected row endpoint, not through a fictitious global operator norm.

```latex
\subsubsection{Support-size deletion inside the completed PS14 high--high row}

Keep the notation and mask completion of PS14, and put
\[
 D=N^{41/100},\qquad
 u_p=\frac{\log p}{p}v(p/N).
\]
For the high-character projection of the reflection-compatible
nonzero-frequency factor write
\[
 H_p^H(a)=\sum_{\chi\in\mathcal E_p}\gamma_p(\chi)\chi(a),
 \qquad
 C_p^H:=\sum_{\chi\in\mathcal E_p}|\gamma_p(\chi)|^2.
\]
By multiplicative Parseval and orthogonal projection,
\[
 \sum_{a\ne0}|H_p^H(a)|^2=(p-1)C_p^H,
 \qquad
 C_p^H\ll Z(p).
\tag{PS14-S1}
\]
For $0\le\alpha\le3/5$ define
\[
 \mathcal S_\alpha
 :=\{p\asymp N:Z(p)\le N^\alpha\},
\]
\[
 E_\alpha
 :=\sum_{p\in\mathcal S_\alpha}|u_p|^2C_p^H,
 \qquad
 L_\alpha
 :=\sum_{p\in\mathcal S_\alpha}|u_p|C_p^H,
\]
and let $E,L$ denote the corresponding sums over all shell primes.
Then
\[
 E_\alpha
 \ll \min\{N^{-2/5},N^{\alpha-1}\}N^{o(1)},
\tag{PS14-S2}
\]
while
\[
 L_\alpha\ll N^{\alpha+o(1)},
 \qquad
 L\ll N^{3/5+o(1)}.
\tag{PS14-S3}
\]
Indeed $|u_p|\ll N^{-1+o(1)}$, $C_p^H\ll Z(p)$,
and the unconditional average-zero bound gives
$\sum_{p\asymp N}|u_p|Z(p)\ll N^{3/5+o(1)}$ and
$\sum_{p\asymp N}|u_p|^2Z(p)\ll N^{-2/5+o(1)}$.

For prime sets $X,Y$ let
\[
 B_{X,Y}(d)
 :=\sum_{\substack{p\in X,\ q\in Y\\p\ne q}}
 u_pu_q\,\omega_{p,q}(d)
 H_p^H(-d\bar q_p)H_q^H(-d\bar p_q),
\]
where $\omega_{p,q}(d)$ denotes the bounded smooth/orientation
factor in the completed PS14 reciprocal row.  Uniformly for
$0<|d|\asymp D$,
\[
 |B_{X,Y}(d)|^2\ll L_XL_YN^{o(1)}.
\tag{PS14-S4}
\]
To see this, apply weighted Cauchy.  For fixed $p$, the map
$q\mapsto-d\bar q_p$ has bounded multiplicity on the fixed-ratio
prime shell, and hence
\[
 \sum_{q\in Y}|u_q|
 |H_p^H(-d\bar q_p)|^2
 \ll N^{-1+o(1)}
      \sum_{a\ne0}|H_p^H(a)|^2
 \ll C_p^HN^{o(1)}.
\]
The symmetric estimate in $q$ proves~\textup{(PS14-S4)}.  Summing the
single available $d$-shell gives
\[
 \sum_{0<|d|\asymp D}|B_{X,Y}(d)|^2
 \ll D L_XL_YN^{o(1)}.
\tag{PS14-S5}
\]
There is no additional $D^2$ at this stage.

Now split each of the two characteristic slots in each row by
$1_{\mathcal S_\alpha}+1_{\mathcal S_\alpha^c}$.  Consider one Gram
block and let $k$ be the number of its four characteristic occurrences
which are constrained to $\mathcal S_\alpha$.  By Cauchy between the
two row shells and~\textup{(PS14-S3)--(PS14-S5)},
\[
 \boxed{
 |\mathfrak G_k(\alpha)|
 \ll
 N^{161/100+(k/2)(\alpha-3/5)+o(1)}.}
\tag{PS14-S6}
\]
Consequently:
\begin{enumerate}
\item If $\alpha=14/25-\delta$, every block with $k\ge1$ is
\[
 \ll N^{159/100-\delta/2+o(1)}.
\]
\item A projected self-shell with one small source characteristic on
both sides has $k=2$; with $\alpha=29/50-\delta$ it is
\[
 \ll N^{159/100-\delta+o(1)}.
\]
\item If all four characteristic occurrences are small, then $k=4$;
with $\alpha=59/100-\delta$ the block is
\[
 \ll N^{159/100-2\delta+o(1)}.
\]
\end{enumerate}

These estimates are applied after the exact PS14 completion of the
high-centered/near/alias/orientation masks.  The literal four-distinct
condition is retained by the preceding prime-intersection
inclusion--exclusion; it is not removed by monotonicity.  All
same-$P$/one-shared correction strata are already banked.  The midpoint
and reflection-centering correction terms remain below the critical
shell by the preceding PS14 estimates, and their crosses with the
support-projected main term are bounded by Cauchy.

\begin{remark}
The coefficient estimate~\textup{(PS14-S2)} is in the Mellin
coefficient normalization.  The canonical additive coefficient square
$Q_P$ is larger by one factor $p\asymp N$.  Moreover,
\textup{(PS14-S6)} is proved from the actual-row endpoint
\textup{(PS14-S4)--(PS14-S5)}, not from a full high-character
operator inequality.  No one-prime fourth-moment input is used.
\end{remark}
```

---

# 11. Final premise audit

| item requested in Q4980 | audit |
|---|---|
| `u_p` normalization | `u_p=(\log p)/p\,v(p/N)`; **valid**. Raw additive `Q_P` carries an extra factor `p` relative to Mellin `E`. |
| `E_\alpha\ll\min(N^{-2/5},N^{\alpha-1})` | **valid only for Mellin/high-character coefficient mass** `\sum|u_p|^2C_p^H`. |
| high-character projection | **valid/contractive**: `C_p^H\le C_p\ll Z(p)`. |
| pair tensor norm | **valid**: `\|a_{X,Y}\|_2^2\le E_XE_Y`. |
| full operator promotion from pair tensor | **invalid / not a source theorem**; this is the first false line in the proposed proof. |
| Gram exponent `161/100+(k/2)(alpha-3/5)` | **valid after repair**, via support-projected fixed-`d` Cauchy + finite Parseval. |
| `D` versus `D^2` | **one `D`** in the shell endpoint. Another `D^2` would double-count Fejér restoration. |
| all-four-distinct mask | retain by **exact inclusion--exclusion**, not monotone deletion. |
| high-centered / near / alias masks | apply support lemma after exact completion; no new loss. |
| orientations / circular branches | finite exact partition; `N^{o(1)}` only. |
| midpoint terms | largest cross `N^{131/100+o(1)}`; safely subcritical. |
| diamond correction | correction square `N^{81/100+o(1)}`, energy cross `N^{121/100+o(1)}`; safely subcritical. |
| support cross terms | Cauchy gives exactly the `k/2` exponent. |
| `alpha=14/25-delta` any-small deletion | **valid**, margin `delta/2`. |
| `alpha=29/50-delta` one-small-each-side self-shell | **valid**, margin `delta`. |
| `alpha=59/100-delta` all-four-small | **valid**, margin `2delta`. |
| Q4928 one-prime `M_4` bridge | **not used**. |

## Bottom line

Q4958 found the right numerical cutoffs, but its cleanest apparent proof path mixes two different Hilbert-space normalizations.  The small-support Mellin coefficient energy really is `N^{\alpha-1}`, and the pair tensor really has the corresponding `L^2` saving, but the source's `161/100` endpoint is not a uniform full-space operator norm to which one may simply feed that smaller tensor.

The deletion nevertheless survives.  The actual PS14 reciprocal row can be support-projected **before the same fixed-`d` Cauchy/Parseval endpoint is rerun**.  That gives the row bound `D L_XL_Y`, and Cauchy between two projected rows yields exactly

```latex
N^{161/100+(k/2)(\alpha-3/5)+o(1)}.
```

Hence the three support thresholds `14/25`, `29/50`, and `59/100` are source-correct after the proof repair above.