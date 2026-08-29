# Q32 RR68 average repair

## Verdict

The proposed average repair is **proved at the staged RR68 level**, subject only to the literal RR68c hypotheses stated below and the already-recorded physical range \(1\le B\ll_v N\). The four requested estimates are

\[
 w_k\ll N^{-3+o(1)}z_km_k,\qquad
 m_k\ll_v {N^{1+o(1)}\over B},\qquad
 \sum_k z_k\ll N^{8/5+o(1)},\qquad
 W(	heta)\ll {N^{-2/5+o(1)}\over B}.
\]

The decisive point is that, for fixed \(	heta\) and fixed source prime \(p\), the congruence fixes the **integer \(k\)** modulo \(b\), and the physical interval for \(k\) has length \(O_v(B)\) while \(b\asymp_v B\). Hence \(p\) belongs to only \(O_v(1)\) exact-\(k\) packets. The overlap sum is therefore controlled by the canonical first zero moment

\[
 \sum_{p\asymp N}Z(p)\ll N^{8/5+o(1)}.
\]

There is no additional \(Z(p)\) multiplicity in this overlap: the zero-label multiplicity for \(p\) is already the factor \(Z(p)\) inside \(z_k\).

Combining this pointwise-in-\(	heta\) bound with the staged positive \(L^1\) ledger gives

\[
 \sum_	heta W(	heta)^2
 \ll {N^{2/5+o(1)}M\over B},\qquad
 M=\min\!\left(B,{D\over N^{2/5}}ight).
\]

The literal host summation then restores

\[
 E_{m mean}\ll FMK,\qquad
 E_{m cent}\ll N^{3/5+o(1)}K(BMF)^{1/2}.
\]

Thus the **staged** local ridge \(\sigma<1/40\) is restored: the downstream optimizer receives exactly the two estimates it was designed to consume. The staged opposite high-corner minimax is also unchanged:

\[
 \min_\ell\max\!\left({5\over2}-{\ell\over4},{7\over3}+\ellight)
 ={37\over15}\quad	ext{at}\quad \ell={2\over15}.
\]

This high corner is still above the staged fully unconditional target \(12/5=36/15\) by \(1/15\). Consequently the repair does **not** close Problem 3.2.

The separate third-moment sufficient premise remains unproved. From the actual pointwise estimate \(Z(r)\ll N^{2/3+o(1)}\) and AF1 one obtains only

\[
 \sum_{r\asymp N}Z(r)^3\ll N^{43/15+o(1)}.
\]

Nothing here proves either \(N^{29/15-6\sigma+o(1)}\) or \(N^{2-6\sigma+o(1)}\).

---

## 1. Source and dependency boundary

This document is a research proof, not an authoritative insertion.

The connector-visible repository pins used for the audit were

```text
repository: xiangyazi24/Ramanujan_Challenge
canonical main: c5d932b66ce5e4f1657b587215d290ae7a13018b
chatgpt-drop before this write: b012211d9ff90db90e03726213ec40dc840bb6e0
canonical oracleA_result.tex blob on main:
  f0bfccc441bda22c17658c5586fbb2e6c6431238
```

The task's source audit states that the authoritative AF/RC/RR chain in the current canonical working source reaches only roughly RR48, while RR49--RR68 is staged research. The connector-visible `chatgpt-drop` branch is itself a staging/drop branch and does not provide a connected authoritative `problems/3.2/oracleA_result.tex` surface containing RR68c and its consumers. Therefore:

- `problems/3.2/oracleA_result.tex` is **not edited**;
- `problems/3.2/proof.tex` is **not edited**;
- the future insertion block is recorded in Section 10 below, to be used only after RR49--RR68 and their definitions have been promoted into the authoritative source.

### 1.1 Canonical analytic inputs

Only the following current zero-count inputs are used:

\[
 Z(p)\ll p^{2/3+o(1)},
 	ag{A0}
\]

\[
 \sum_{p\asymp N}Z(p)\ll N^{8/5+o(1)},
 	ag{A1}
\]

\[
 \sum_{p\asymp N}|u_p|Z(p)\ll N^{3/5+o(1)},
 	ag{A2}
\]

\[
 \sum_{p\asymp N}|u_p|Z(p)^2\ll N^{6/5+o(1)},
 	ag{A3}
\]

and

\[
 |u_p|\ll N^{-1+o(1)}.
 	ag{A4}
\]

The proof of the packet overlap uses (A1), not a false pointwise estimate \(Z(p)\ll N^{3/5}\).

### 1.2 Staged RR68 inputs

Fix one literal unrefined radial datum

\[
 	heta=(q,\eta,b,c_0,arepsilon),\qquad
 C=b\eta+arepsilon q c_0=bn+ho,\qquad b\asymp_v B.
 	ag{RR68-a}
\]

For every literal RR68c occurrence, the two source primes \(p,p'\) share one genuine integer \(k\) satisfying

\[
 |k|\le C_vB,\qquad (k,b)=(ho,b)=1,
 	ag{RR68-b}
\]

and

\[
 pk\equiv-ho\pmod b,\qquad
 p'k\equiv-ho\pmod b.
 	ag{RR68-c}
\]

The source-prime shell has fixed compact support \(p,p'\asymp_v N\), and the staged physical radial range has \(1\le B\ll_v N\). All source, sign, primitive, actual-zero, chart, radial, and lift masks are retained. When a positive upper bound drops one of these masks, it is dropped only as a \(0\)-\(1\) deletion.

For an exact integer \(k\), let

\[
 \mathcal P_k=\mathcal P_k(	heta)
\]

be the actual restricted source-prime set surviving all those masks. Put

\[
 m_k=|\mathcal P_k|,\qquad
 z_k=\sum_{p\in\mathcal P_k}Z(p).
 	ag{1.1}
\]

Let \(C_k\) be the positive first-carrier submass and

\[
 w_k=|u_q|C_k,\qquad
 W(	heta)=\sum_k w_k.
 	ag{1.2}
\]

The literal staged factorization is

\[
 C_k\le N^{o(1)}
 \left(\sum_{p\in\mathcal P_k}|u_p|Z(p)ight)
 \left(\sum_{p'\in\mathcal P_k}|u_{p'}|ight).
 	ag{RR68-F}
\]

Finally, the staged positive first-moment ledger is

\[
 \sum_	heta W(	heta)\le N^{4/5+o(1)}M,\qquad
 M=\min\!\left(B,{D\over N^{2/5}}ight).
 	ag{RR68-L1}
\]

The host consumer to be audited is

\[
 I_r\le {FZ(r)\over r}\sum_	heta W(	heta)
 +\sqrt{FZ(r)(r-Z(r))}\,B
  \left(\sum_	heta W(	heta)^2ight)^{1/2},
 	ag{RR-host-I}
\]

\[
 E_{m nonper}\le N^{-1+o(1)}K
 \sum_{r\asymp N}|u_r|Z(r)I_r.
 	ag{RR-host-E}
\]

---

## 2. Marked-occurrence multiplicity audit

There are three distinct multiplicities, and conflating them is the main risk.

### 2.1 Zero labels on the first source-prime factor

For a fixed prime \(p\), the first factor in (RR68-F) may carry an actual zero label. Bounding all surviving labels costs at most \(Z(p)\). This is exactly why the first factor is

\[
 \sum_{p\in\mathcal P_k}|u_p|Z(p).
\]

The same zero labels must not be counted again when the exact-\(k\) packets are overlapped. They are already inside \(z_k\).

### 2.2 The second source-prime factor

The second factor in (RR68-F) contains no additional zero-count weight. Fixed signs, charts, primitive branches, and physical lift choices have bounded multiplicity for fixed \(	heta\) and exact \(k\); the staged factorization has already placed this bounded multiplicity in \(N^{o(1)}\). Thus it is bounded by

\[
 N^{o(1)}\sum_{p'\in\mathcal P_k}|u_{p'}|,
\]

not by a second copy of \(z_k\).

### 2.3 Overlap of one prime among different \(k\)'s

Define

\[
 
u_	heta(p)=\#\{k:p\in\mathcal P_k\}.
 	ag{2.1}
\]

The proof below gives

\[
 
u_	heta(p)\ll_v1.
 	ag{2.2}
\]

This is a prime-to-\(k\) overlap. It does not depend on how many actual zero labels the prime has. Consequently

\[
 \sum_k z_k
 =\sum_{p\asymp N}Z(p)
u_	heta(p),
 	ag{2.3}
\]

with one, not two, copies of \(Z(p)\).

---

## 3. The exact-\(k\) carrier bound

### Proposition 3.1

For every exact physical \(k\),

\[
 oxed{w_k\ll N^{-3+o(1)}z_km_k.}
 	ag{3.1}
\]

### Proof

By (A4),

\[
 \sum_{p\in\mathcal P_k}|u_p|Z(p)
 \le N^{-1+o(1)}z_k,
 	ag{3.2}
\]

and

\[
 \sum_{p'\in\mathcal P_k}|u_{p'}|
 \le N^{-1+o(1)}m_k.
 	ag{3.3}
\]

Substituting these two estimates into (RR68-F) gives

\[
 C_k\le N^{-2+o(1)}z_km_k.
 	ag{3.4}
\]

The prime \(q\) is fixed in \(	heta\) and contributes its source weight exactly once. Again by (A4),

\[
 |u_q|\le N^{-1+o(1)}.
\]

Therefore

\[
 w_k=|u_q|C_k\le N^{-3+o(1)}z_km_k.
\]

Every omitted literal mask is a deletion from the two nonnegative factors in (RR68-F). No cancellation and no replacement of an actual-zero condition by a density is used. \(\square\)

---

## 4. The size of one exact-\(k\) prime packet

### Proposition 4.1

Uniformly in every physical \(k\),

\[
 oxed{m_k\ll_v {N^{1+o(1)}\over B}.}
 	ag{4.1}
\]

More robustly, without using \(B\ll N\),

\[
 m_k\ll_v1+{N\over B}.
 	ag{4.2}
\]

### Proof

Since \((k,b)=1\), the congruence

\[
 pk\equiv-ho\pmod b
\]

has the unique solution class

\[
 p\equiv-ho\,\overline{k}\pmod b.
 	ag{4.3}
\]

The fixed compact source-prime shell lies in an interval of length \(O_v(N)\). Any residue class modulo \(b\) has at most

\[
 1+O_v(N/b)
\]

integer representatives in that interval, and the prime representatives form a subset. Since \(b\asymp_vB\), this proves (4.2). In the staged physical range \(B\ll_vN\), the \(1\) is absorbed by \(N/B\), proving (4.1). \(\square\)

### Hostile boundary note

If a future RR version permits \(B>N\), the literal statement must retain \(1+N/B\); the displayed \(N/B\) bound is then false for a nonempty packet. The present RR68 physical range is the reason the shorter form is legal.

---

## 5. Fixed-prime overlap across exact \(k\)'s

### Proposition 5.1

For fixed \(	heta\) and fixed shell prime \(p\),

\[
 \#\{k:p\in\mathcal P_k\}\ll_v1.
 	ag{5.1}
\]

Consequently

\[
 oxed{\sum_k z_k\ll N^{8/5+o(1)}.}
 	ag{5.2}
\]

### Proof

First, \(p\) is invertible modulo \(b\). Indeed, if \(d\mid(p,b)\), then \(d\mid pk\), while RR68c gives \(pk\equiv-ho\pmod b\); hence \(d\midho\). Since \((ho,b)=1\), one gets \(d=1\). Thus

\[
 (p,b)=1.
 	ag{5.3}
\]

The same congruence now fixes \(k\) modulo \(b\):

\[
 k\equiv-ho\,\overline p\pmod b.
 	ag{5.4}
\]

All physical representatives satisfy

\[
 |k|\le C_vB.
\]

This interval has length \(2C_vB\), while \(b\asymp_vB\). Therefore it contains at most

\[
 1+{2C_vB\over b}=O_v(1)
\]

integers in the residue class (5.4). All source, sign, primitive, chart, and lift masks can only delete representatives. This proves (5.1).

Now expand the marked overlap exactly:

\[
 \sum_k z_k
 =\sum_k\sum_{p\in\mathcal P_k}Z(p)
 =\sum_{p\asymp N}Z(p)\#\{k:p\in\mathcal P_k\}.
 	ag{5.5}
\]

By (5.1) and the canonical first moment (A1),

\[
 \sum_kz_k\ll_v\sum_{p\asymp N}Z(p)\ll N^{8/5+o(1)}.
\]

There is no additional zero-label factor in (5.5). A fixed prime's actual zero labels have already been majorized by the single factor \(Z(p)\) in the definition of \(z_k\). \(\square\)

### Corollary 5.2: common residue class for the two source primes

The common integer \(k\) is load-bearing. Since \((k,b)=1\), RR68c also gives

\[
 p\equiv p'\equiv-ho\,\overline k\pmod b.
 	ag{5.6}
\]

If the two source occurrences had unrelated integers \(k,k'\), they would not lie in one common packet \(\mathcal P_k\), and the literal factorization (RR68-F) could not be bounded by \(z_km_k\) as above.

---

## 6. The repaired pointwise radial mass

### Theorem 6.1: RR68 average repair

For every fixed literal unrefined radial datum \(	heta\),

\[
 oxed{W(	heta)\ll {N^{-2/5+o(1)}\over B}.}
 	ag{6.1}
\]

### Proof

All terms are nonnegative. Propositions 3.1, 4.1, and 5.1 give

\[
egin{aligned}
 W(	heta)
 &=\sum_kw_k\\
 &\ll N^{-3+o(1)}\sum_kz_km_k\\
 &\le N^{-3+o(1)}\left(\max_km_kight)\sum_kz_k\\
 &\ll N^{-3+o(1)}\cdot{N\over B}\cdot N^{8/5}\\
 &= {N^{-2/5+o(1)}\over B}.
\end{aligned}
 	ag{6.2}
\]

The exponent identity is

\[
 -3+1+{8\over5}=-{2\over5}.
 	ag{6.3}
\]

This proof uses the first zero moment exactly once. It does not use a pointwise \(3/5\) estimate. \(\square\)

---

## 7. The repaired \(L^2\) radial ledger

### Corollary 7.1

Assuming the staged positive first-moment bound (RR68-L1),

\[
 oxed{\sum_	heta W(	heta)^2
 \ll {N^{2/5+o(1)}M\over B}.}
 	ag{7.1}
\]

### Proof

By nonnegativity,

\[
 \sum_	heta W(	heta)^2
 \le \left(\max_	heta W(	heta)ight)\sum_	heta W(	heta).
 	ag{7.2}
\]

Apply Theorem 6.1 and (RR68-L1):

\[
 \sum_	heta W(	heta)^2
 \ll {N^{-2/5+o(1)}\over B}\cdot N^{4/5+o(1)}M
 ={N^{2/5+o(1)}M\over B}.
\]

The exponent identity is

\[
 -{2\over5}+{4\over5}={2\over5}.
\]

\(\square\)

---

## 8. Literal host summation

Put

\[
 S_1=\sum_	heta W(	heta),\qquad
 S_2=\sum_	heta W(	heta)^2.
\]

Then

\[
 S_1\ll N^{4/5+o(1)}M,\qquad
 S_2\ll {N^{2/5+o(1)}M\over B}.
 	ag{8.1}
\]

We sum over \(r\) before applying any pointwise zero estimate.

### 8.1 Mean term

The mean part of (RR-host-I) is

\[
 I_r^{m mean}\le {FZ(r)\over r}S_1.
\]

Substituting into (RR-host-E) gives

\[
 E_{m mean}\ll N^{-1+o(1)}KF S_1
 \sum_{r\asymp N}{|u_r|Z(r)^2\over r}.
 	ag{8.2}
\]

Since \(r\asymp N\), (A3) yields

\[
 \sum_{r\asymp N}{|u_r|Z(r)^2\over r}
 \ll N^{-1+o(1)}N^{6/5}=N^{1/5+o(1)}.
\]

Together with \(S_1\ll N^{4/5+o(1)}M\),

\[
 E_{m mean}\ll
 N^{-1}KF\cdot N^{4/5}M\cdot N^{1/5+o(1)}
 =oxed{FMK\,N^{o(1)}}.
 	ag{8.3}
\]

Equivalently, before combining the final \(r^{-1}\), the complete exponent check is

\[
 -1-1+{4\over5}+{6\over5}=0.
 	ag{8.4}
\]

The weighted second zero moment, not the weighted first moment, is the correct input here.

### 8.2 Centered term

The centered part of (RR-host-I) is

\[
 I_r^{m cent}\le
 \sqrt{FZ(r)(r-Z(r))}\,B\,S_2^{1/2}.
\]

Hence

\[
 E_{m cent}\ll
 N^{-1+o(1)}K\sqrt F\,B\,S_2^{1/2}
 \sum_{r\asymp N}|u_r|Z(r)^{3/2}\sqrt{r-Z(r)}.
 	ag{8.5}
\]

Since \(0\le Z(r)\le r\asymp N\),

\[
 \sqrt{r-Z(r)}\le N^{1/2+o(1)}.
 	ag{8.6}
\]

Weighted Cauchy--Schwarz gives

\[
egin{aligned}
 \sum_{r\asymp N}|u_r|Z(r)^{3/2}
 &=\sum_rigl(|u_r|Z(r)igr)^{1/2}
          igl(|u_r|Z(r)^2igr)^{1/2}\\
 &\le
 \left(\sum_r|u_r|Z(r)ight)^{1/2}
 \left(\sum_r|u_r|Z(r)^2ight)^{1/2}\\
 &\ll N^{(3/5+6/5)/2+o(1)}\\
 &=N^{9/10+o(1)}.
\end{aligned}
 	ag{8.7}
\]

Furthermore

\[
 S_2^{1/2}\ll N^{1/5+o(1)}\left({M\over B}ight)^{1/2}.
 	ag{8.8}
\]

Substitution into (8.5) yields

\[
egin{aligned}
 E_{m cent}
 &\ll N^{-1+1/2+9/10+1/5+o(1)}
 K\sqrt F\,B\left({M\over B}ight)^{1/2}\\
 &=oxed{N^{3/5+o(1)}K(BMF)^{1/2}}.
\end{aligned}
 	ag{8.9}
\]

The numerical check is

\[
 -1+{1\over2}+{9\over10}+{1\over5}={3\over5}.
 	ag{8.10}
\]

No pointwise \(Z(r)\ll N^{3/5}\) is used. The only pointwise use is the tautological \(\sqrt{r-Z(r)}\le\sqrt r\).

---

## 9. Downstream exponent status

### 9.1 Local ridge

The RR68 average repair is independent of the local \(\sigma\)-optimization. Its entire downstream output is precisely (8.3) and (8.9), with the same \(F,M,K,B\) dependence as the intended staged consumer. Therefore the staged local optimizer receives no exponent loss, and its strict range

\[
 oxed{\sigma<1/40}
\]

is restored.

This is a staged implication, not an authoritative theorem insertion: the definitions and consumer carrying \(\sigma\) are in RR49--RR68 research, not in the authoritative canonical source surface.

### 9.2 Opposite high corner

The staged opposite-corner ledger is the maximum of the two affine exponent lines

\[
 H_1(\ell)={5\over2}-{\ell\over4},\qquad
 H_2(\ell)={7\over3}+\ell.
 	ag{9.1}
\]

Their difference is

\[
 H_1(\ell)-H_2(\ell)={1\over6}-{5\ell\over4}.
\]

They cross at

\[
 \ell={2\over15},
\]

and

\[
 H_1(2/15)=H_2(2/15)={37\over15}.
 	ag{9.2}
\]

Thus

\[
 oxed{\min_\ell\max(H_1(\ell),H_2(\ell))={37\over15}.}
 	ag{9.3}
\]

The staged fully unconditional target is \(12/5=36/15\). The repaired high corner remains above it by

\[
 {37\over15}-{12\over5}={1\over15}.
 	ag{9.4}
\]

The average repair therefore restores the intended local ridge but does not close the opposite high corner.

### 9.3 The third-moment premise is still absent

AF1 and the actual pointwise bound give only

\[
egin{aligned}
 \sum_{r\asymp N}Z(r)^3
 &\le \left(\max_{r\asymp N}Z(r)ight)\sum_{r\asymp N}Z(r)^2\\
 &\ll N^{2/3+o(1)}N^{11/5+o(1)}\\
 &=N^{43/15+o(1)}.
\end{aligned}
 	ag{9.5}
\]

The exponent identity is

\[
 {2\over3}+{11\over5}={43\over15}.
\]

This is far weaker than either staged sufficient premise

\[
 \sum Z(r)^3\ll N^{29/15-6\sigma+o(1)}
 	ag{9.6}
\]

or

\[
 \sum Z(r)^3\ll N^{2-6\sigma+o(1)}.
 	ag{9.7}
\]

Neither (9.6) nor (9.7) is proved here.

### 9.4 Smallest remaining residual

At the current staged endpoint, the smallest explicitly quantified residual is the opposite high-corner excess \(N^{1/15+o(1)}\) over the fully unconditional target. One sufficient route is a genuinely stronger third zero moment of the form (9.6) or (9.7), but such a moment does not follow from AF1 and \(Z(r)\ll r^{2/3+o(1)}\).

Accordingly the final status is:

```text
(i) RR68 common-k average repair:             PROVED
(ii) staged local sigma < 1/40 ridge:          RESTORED
(iii) staged high-corner exponent 37/15:       UNCHANGED
(iv) fully unconditional target 12/5:          NOT REACHED
(v) third-moment sufficient premise:           UNPROVED
(vi) full Problem 3.2:                          OPEN
```

---

## 10. Future authoritative insertion block

Do **not** insert this block until the RR49--RR68 definitions, the exact-\(k\) packet \(\mathcal P_k(	heta)\), \(C_k,w_k,W(	heta)\), and the host consumer have been promoted into `oracleA_result.tex`.

```latex
\begin{lemma}[Average repair for the RR68 common radial integer]
\label{lem:oracleA-rr68-average-repair}
Fix one literal unrefined radial datum
\[
 \theta=(q,\eta,b,c_0,\varepsilon),\qquad
 b\asymp_v B,\qquad 1\le B\ll_vN,
\]
and write
\[
 b\eta+\varepsilon q c_0=bn+\rho.
\]
Assume that every literal RR68c occurrence belonging to an exact
integer $k$ satisfies
\[
 |k|\ll_vB,\qquad (k,b)=(\rho,b)=1,
\]
and, for its two source primes,
\[
 pk\equiv p'k\equiv-\rho\pmod b.
\]
Let $\mathcal P_k$ be the actual restricted source-prime set, and put
\[
 m_k=|\mathcal P_k|,\qquad
 z_k=\sum_{p\in\mathcal P_k}Z(p),\qquad
 w_k=|u_q|C_k.
\]
Suppose the literal positive carrier factorization is
\[
 C_k\ll N^{o(1)}
 \left(\sum_{p\in\mathcal P_k}|u_p|Z(p)\right)
 \left(\sum_{p'\in\mathcal P_k}|u_{p'}|\right).
\]
Then
\[
 w_k\ll N^{-3+o(1)}z_km_k,\qquad
 m_k\ll_v {N^{1+o(1)}\over B},
\]
\[
 \sum_kz_k\ll N^{8/5+o(1)},
\qquad
 W(\theta):=\sum_kw_k
 \ll {N^{-2/5+o(1)}\over B}.
\]
All literal source, sign, primitive, actual-zero, chart and lift masks
may be retained.
\end{lemma}

\begin{proof}
The source weights give
\[
 C_k\ll N^{-2+o(1)}z_km_k,
 \qquad w_k\ll N^{-3+o(1)}z_km_k.
\]
For fixed $k$, the congruence and $(k,b)=1$ put every source prime in
one residue class modulo $b$, whence
$m_k\ll1+N/b\ll N/B$.  Conversely, for fixed $p$, the congruence and
$(\rho,b)=1$ first imply $(p,b)=1$ and then fix
$k\equiv-\rho\bar p\pmod b$.  The physical interval $|k|\ll_vB$ has
length $O_v(B)$ and $b\asymp_vB$, so it contains $O_v(1)$ such lifts.
Therefore
\[
 \sum_kz_k
 \ll_v\sum_{p\asymp N}Z(p)
 \ll N^{8/5+o(1)}.
\]
Combining these estimates proves the result.  The factor $Z(p)$ already
counts all zero labels for $p$ and is not repeated in the $k$-overlap.
\end{proof}

\begin{corollary}[Repaired RR68 radial square]
\label{cor:oracleA-rr68-average-square}
If
\[
 \sum_\theta W(\theta)\ll N^{4/5+o(1)}M,
 \qquad M=\min\!\left(B,{D\over N^{2/5}}\right),
\]
then
\[
 \sum_\theta W(\theta)^2
 \ll {N^{2/5+o(1)}M\over B}.
\]
Consequently the literal RR host inequalities give
\[
 E_{\rm mean}\ll FMK,
 \qquad
 E_{\rm cent}\ll N^{3/5+o(1)}K(BMF)^{1/2}.
\]
\end{corollary}
```

---

## 11. Computation and verification scope

The companion script

```text
problems/3.2/research/scripts/q32_rr68_average_repair_verify.py
```

checks the rational exponent arithmetic and exhaustively tests the elementary residue/lift assertions for small moduli, physical intervals, and prime shells. That computation is **not** the asymptotic proof. The asymptotic proof is Sections 3--8 above and uses only interval lengths, coprimality, and the canonical moment estimates.

The intended verification commands are

```bash
git diff --check
python3 problems/3.2/research/scripts/q32_rr68_average_repair_verify.py
cd problems/3.2
latexmk -pdf -interaction=nonstopmode -halt-on-error proof.tex
```

Because this research write does not edit or include either TeX file, the LaTeX source graph is unchanged. A successful build checks the pre-existing canonical document, not the staged RR68 block.
