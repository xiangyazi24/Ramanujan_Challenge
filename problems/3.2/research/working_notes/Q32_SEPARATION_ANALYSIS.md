# Q3.2 quotient-separation audit: Casoratian and defect-dimension routes

**Date:** 2026-07-29

**Source under audit:** `/tmp/q_separation_sol.txt` and the scripts in the same
temporary scratch directory.

**Goal:** assess whether the proposed CRT quotient-separation mechanism can prove
the fully pointwise estimate `log G_n=o(n)`, with emphasis on Route A
(Casoratian) and Route B (defect-dimension counting).

## 1. Executive verdict

The proposed separation mechanism does **not** currently give a route to the
pointwise theorem.

1. The Casoratian in the source note has the wrong sign:

   \[
   a_n b_{n+1}-a_{n+1}b_n=-\frac{6}{(n+1)^3},
   \]

   with no alternating factor. The generating script prints failures at every
   odd index but continues and later tests only even outer indices, thereby
   masking the error.

2. After the sign is repaired and denominators are treated correctly, the
   Casoratian gives an exact formula for the local quotient `b_n/p`. It does
   not give an independent condition on that quotient. Replacing `b_n` by
   `R_nX` produces a defect which factors **exactly** as a fixed integer times
   `R_n(A_n-X)`. If `X≡A_n (mod R_n^k)`, the cleared defect contains
   `R_n^{k+1}` for the tautological reason that `A_n-X` already contains
   `R_n^k`. There is no compression.

3. The cited Liu 2024 theorem is materially different from the claimed
   endpoint tower. Ji-Cai Liu, arXiv:2404.16636, Theorem 1.1 proves one extra
   Gauss digit (here `\mathbf A` is the generalized Apéry sequence, not the
   quotient defined below):

   \[
   \mathbf A_{np^m}^{(r,s,t)}
   \equiv \mathbf A_{np^{m-1}}^{(r,s,t)}
      +p^{3m}B_{p-3}\,\mathcal A_n^{(r,s,t)}
      \pmod {p^{3m+1}}.
   \]

   It concerns indices `np^m`, adds only `B_{p-3}`, and for `m=1` is a
   congruence modulo `p^4`, not `p^5`. It does not state an arbitrary-order
   shifted endpoint expansion for `n=p+r`, and it does not supply the asserted
   sequence `B_{p-3},B_{p-5},...`.

4. Sun's relevant statements do not supply the missing tower either.
   Conjecture 5.1 of arXiv:2005.02081 is only modulo `p^4` for a few fixed
   multiples. Conjecture 2.4 of arXiv:2409.06544 predicts a modulo-`p^5`
   formula at multiplicative indices, but it is unproved and contains one
   fixed combination of `B_{2p-4}` and `B_{p-3}`, not the claimed sequence
   `B_{p-3},B_{p-5},...`. Thus the unconditional dimensions at
   `p^5,p^7,p^9` are unknown. The temporary dimension script is a speculative
   free-polynomial model, not a computation of an Apéry defect module. It is
   internally inconsistent even as a speculative model.

5. The inverse-limit formulation is backwards. If each finite congruence
   really computes the residue of the fixed integer `A_n` modulo `R_n^k`,
   then the compatible inverse-limit element is **automatically** the diagonal
   image of `A_n`, and its least nonnegative representatives eventually
   stabilize at `A_n`. Proving non-diagonality would contradict the congruence
   input, not prove the desired estimate.

6. A fixed level such as `k=7` could at best improve a linear constant. To
   prove `log R_n=o(n)` by reconstruction, one needs a separation theorem
   uniformly at arbitrarily large fixed levels `k` (or a controlled growing
   level), plus a genuinely independent nonzero low-height certificate.

7. The later small-zero-fiber computation does not change this verdict.
   Even \(|Z_p|\le2\) for every prime would not bound one fixed integer
   column: reflection-invariant two-point rows can all be aligned at the
   same \(n\). The exact scan through two million instead supports a
   Poisson law with slowly unbounded row and column maxima. The claimed
   evenness also has the central exceptions \(p=11,3137\).

8. The clean surviving target is horizontal. It is enough to prove, for
   every fixed quotient \(q\), a power saving for the shell pair energy of
   the affine zero sets \(qp+Z_p\). This is a cross-prime small-CRT theorem,
   not a consequence of a local Casoratian or of defect-space dimension.

9. There is nevertheless one unconditional pointwise pruning. A
   divisor-sensitive Kummer-order sieve shows that all middle primes with
   \[
   \gcd(p-1,n-q)>
   \tau(n-q)\log n\,F(n),\qquad q=\lfloor n/p\rfloor,
   \]
   have total logarithmic weight \(o(n)\) for any fixed subpolynomial choice
   \(F(n)\to\infty\), even before the Apéry zero condition is imposed. Thus
   the unresolved primes select characters of order \(p^{1-o(1)}\). This
   sharpens the target, but does not supply the missing cross-prime
   nonalignment.

10. The safe large-prime depth bound is \(v_p(G_n)\le6\), not 3. The
    claimed bound 3 in `proof.tex` uses \(a_{n-1}\) as an integer Bézout
    coefficient although it is generally rational. Clearing both adjacent
    denominators gives the valid constant 6. This changes no support-to-gcd
    little-oh implication.

11. The entire individual-gap factorial-moment variant of Route A can be
    optimized explicitly. In the direct prefix window it gives
    \[
    \sum_{z_p(H)>B}\log p\ll H^3\log H/B^3.
    \]
    At the available \(B\asymp H^{2/3}\) scale this is only
    \(O(H\log H)\). Higher moments are worse; the missing improvement is a
    joint short-gap compression identity.

12. A sharper top-half reduction is available. A prime-gap collision is
    **degenerate** if either characteristic sees zeros at both endpoint
    indices; the continuant then proves the pointwise bound
    \(M_h^{\rm deg}(n)\ll h\). On every remaining **pure-cross** collision,
    the same continuant is a unit modulo both primes. Consequently the full
    pair-energy theorem can be weakened to one logarithmic-gap statement:
    adjacent pure-cross target pairs of gap \(O(\log n)\) must be
    \(o(n/\log n)\). This is now the narrowest clean top-half obstruction.

13. The CFVZ square-root structure does not yield a hidden low-order
    Casoratian. The reversed convolution \(b_r=\sum_i s_i s_{r-i}\) admits
    no uniformly bounded rational boundary telescoper: the natural
    two-term certificate fails by a translation-by-two divisor imbalance,
    while a general bounded-width collapse would contradict the
    \({\rm SL}_2\) differential Galois group. Its minimal outer telescoper
    is exactly the original symmetric-square Apéry operator.

Thus Routes A and B are closed **as presently formulated**. This is not a
disproof of `log G_n=o(n)`. It identifies the exact missing ingredient:
cross-prime, Apéry-specific Archimedean/p-adic coupling not already contained
in the local residue statement.

## 2. Scope and notation

Let

\[
b_n=\sum_{j=0}^n\binom nj^2\binom{n+j}j^2,
\qquad (b_0,b_1)=(1,5),
\]

and let the rational companion solution have `(a_0,a_1)=(0,6)`. Both satisfy

\[
(n+1)^3u_{n+1}=P(n)u_n-n^3u_{n-1},
\quad
P(n)=34n^3+51n^2+27n+5.
\]

Following the source note, let

\[
R_n=\prod_{\substack{\sqrt n<p\le n\\p\mid b_n}}p,
\qquad A_n=b_n/R_n.
\]

Let \(d_n=\operatorname{lcm}(1,\ldots,n)^3\) and
\(G_n=\gcd(d_na_n,d_nb_n)\). Since
\(d_na_{n-1},d_nb_{n-1}\in\mathbf Z\), the valid integral combination is

\[
 (d_nb_{n-1})(d_na_n)-(d_na_{n-1})(d_nb_n)
 =\frac{6d_n^2}{n^3}.
\]

Consequently, for \(p\ge5\),

\[
 v_p(G_n)\le6\lfloor\log_p n\rfloor-3v_p(n)
 \le6\lfloor\log_p n\rfloor.
\]

In particular \(v_p(G_n)\le6\) for \(p>\sqrt n\); the total contribution of
\(p\le\sqrt n\) is \(O(\sqrt n)\). Hence an \(o(n)\) logarithmic bound for
the relevant large-prime radical is enough for \(\log G_n=o(n)\). The issue
is identifying which factors of the stronger proxy \(R_n\) are actually
relevant to \(G_n\).

There is a scope caveat. For the full interval `sqrt(n)<p<=n`, the block
criterion proved in the repository is

\[
p\mid G_n\quad\Longleftrightarrow\quad D_qb_r\equiv0\pmod p,
\]

where `D_q` is the companion-block constant. The repository's stronger
identification `D_q≡a_q (mod p)` gives the two-channel form

\[
p\mid G_n\quad\Longleftrightarrow\quad a_qb_r\equiv0\pmod p.
\]

The resume note records that the Wronskian argument alone does not prove
this identification at Hasse-zero indices; it points instead to the audited
explicit companion-sum/Kummer proof of the full block congruence. This
foundation issue does not affect the top-half calculation below, where
`q=1` and `D_1=a_1=6`.

The one-channel equivalence with `p|b_r` is valid for `p>n/2`, where `q=1`
and `a_1=6`. The companion `a_q` channel is already controllable with
`O(n^{2/3})` logarithmic weight, so a sufficiently uniform direct-channel
bound would still close the main problem. Nevertheless, `R_n` as defined
above is not literally the middle-prime part of `G_n`. The source note
silently switches between:

- all middle primes `sqrt(n)<p<=n`;
- only top-half primes `n/2<p<=n`;
- the lower-digit support relevant to `G_n`.

The accompanying `crt_tower_digits.py` claims the first definition in its
docstring but implements only the second.

## 3. Reconstruction threshold: the correct calculation

The dominant root of the Apéry recurrence is

\[
\lambda=17+12\sqrt2=(1+\sqrt2)^4,
\qquad c=\log\lambda=3.52549\ldots,
\]

and

\[
\log b_n=cn-\frac32\log n+O(1).
\]

Put `rho_n=(log R_n)/n`. Since `A_n=b_n/R_n`,

\[
\log A_n=(c-\rho_n)n+o(n),
\]

not `cn+o(n)` as stated in `/tmp/q_separation_sol.txt`.

If a level-`k` congruence determines `A_n mod R_n^k`, its least representative
is forced to equal `A_n` once

\[
A_n<R_n^k.
\]

The correct threshold is therefore

\[
cn+o(n)<(k+1)\log R_n,
\qquad\text{or}\qquad
\rho_n>\frac{c}{k+1}+o(1).
\]

The extra `+1` comes from the factor of `R_n` already removed from `b_n`.
This explains the numerical comment in the source note:

\[
c-\frac{7}{2}=0.02549\ldots.
\]

For a top-half radical, the prime number theorem gives the ceiling
`log R_n<=theta(n)-theta(n/2)=n/2+o(n)`. Thus level `k=6` misses reconstruction by
`0.02549n+o(n)`, while level `k=7` is the first level that could improve
the `1/2` constant. The source note's displayed threshold
`cn<k log R_n` is inconsistent with its own statement that `k=7` is first.

More importantly, ruling out

\[
\rho_n>\frac{c}{k+1}+o(1)
\]

for one fixed `k` gives only another positive linear constant. A
little-oh theorem requires the argument for arbitrarily large fixed `k`:
given every `epsilon>0`, choose `k>c/epsilon`, prove the relevant separation
uniformly in the moving `n,p,r` at that level, and only then let `n` grow.
No rate uniform as `k` itself tends to infinity is needed.

## 4. Route A: exact Casoratian audit

### 4.1 Correct Casoratian

Set

\[
W_n=a_nb_{n+1}-a_{n+1}b_n.
\]

Substitution of the recurrence gives

\[
\begin{aligned}
W_n
&=\frac{n^3}{(n+1)^3}
  (a_{n-1}b_n-a_nb_{n-1})\\
&=\frac{n^3}{(n+1)^3}W_{n-1}.
\end{aligned}
\]

Since `W_0=-6`,

\[
\boxed{W_n=-\frac6{(n+1)^3}.}
\]

Equivalently, in the orientation used in the repository,

\[
a_nb_{n-1}-a_{n-1}b_n=\frac6{n^3}.
\]

The temporary script `apery_casoratian.py` makes the sign error exactly
between its lines 47 and 48: line 47 has the expression equal to
`W_{n-1}`, but line 48 inserts an extra minus sign. When executed, it prints

```text
n=1: W_n/W_{n-1} = 1/8, expected -1/8 [FAIL]
n=2: W_n/W_{n-1} = 8/27, expected -8/27 [FAIL]
...
```

and its proposed closed form fails at every odd index. Its later assertions
use only `n=200,300,500`, all even, so the error happens not to trigger an
exception.

### 4.2 Exact local `p`-adic content

Fix a prime

\[
p\in(n/2,n],\qquad p\ge7,\qquad p\mid b_n,
\]

and put `r=n-p`. Then `0<=r<=p-2`: if `r=p-1`, Lucas and reflection give
`b_n≡5b_{p-1}≡5 (mod p)`, contrary to `p|b_n`.

Lucas gives

\[
b_n\equiv5b_r\pmod p,
\]

so `p|b_r`. Two consecutive Apéry values below `p` cannot both vanish,
hence

\[
p\nmid b_{r+1},
\qquad
p\nmid b_{n+1}.
\]

Let

\[
R_n=pS,\qquad
B=b_n=pSA_n,\qquad
C=b_{n+1},
\]

and define the `p`-integral companion values

\[
U=p^3a_n,\qquad V=p^3a_{n+1}.
\]

The block congruence in the audited repository gives

\[
U\equiv6b_r\equiv0\pmod p,
\qquad
V\equiv6b_{r+1}\not\equiv0\pmod p.
\]

Thus `u:=U/p` is `p`-integral and `V,C,S` are `p`-adic units. Multiplying
the correct Casoratian by `p^3` yields

\[
UC-VpSA_n=-\frac{6p^3}{(n+1)^3}.
\]

After division by `p`,

\[
\boxed{
VS A_n
=uC+\frac{6p^2}{(n+1)^3}
\quad\text{in }\mathbf Z_p.
}
\tag{A.1}
\]

Consequently,

\[
\boxed{
A_n
=(VS)^{-1}
\left(uC+\frac{6p^2}{(n+1)^3}\right).
}
\tag{A.2}
\]

Modulo `p`, this reduces to

\[
A_n\equiv (VS)^{-1}uC\pmod p.
\]

This is the exact local content. It expresses the quotient digit in terms of
the normalized companion digit `u=(p^3a_n)/p`; it does not eliminate that
new digit or constrain `A_n` independently.

The formula remains valid when `v_p(b_n)>=2`. In that case `A_n` itself
may still be divisible by `p`, because `R_n` removes only one copy of each
prime.

### 4.3 Exact numerical checks

Independent exact-rational recurrence arithmetic gives:

| `n` | `p` | `S=R_n/p` | `u mod p` | `V mod p` | `C mod p` | `A_n mod p` |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 11 | 5 | 9 | 6 | 5 | 7 |
| 20 | 17 | 1 | 6 | 7 | 3 | 5 |
| 27 | 19 | 11 | 8 | 17 | 11 | 15 |

In each case `(VS)^{-1}uC` gives the last column. For example, at
`n=16,p=11`, the full middle radical is `R_16=5*11`, so

\[
A_{16}\equiv(6\cdot5)^{-1}\cdot9\cdot5\equiv7\pmod {11}.
\]

These checks also give `v_p(U)>=1` and `v_p(V)=0`, exactly as the derivation
requires.

### 4.4 What happens to a fake quotient

Suppose `X` is any integer with

\[
X\equiv A_n\pmod {R_n^k}.
\]

Equation (A.1) gives

\[
uC-VSX+\frac{6p^2}{(n+1)^3}
=VS(A_n-X),
\]

so the left side has `p`-adic valuation at least `k`. In the un-divided
scaled Casoratian, the valuation is at least `k+1`. This is not a new
constraint: it is exactly the original congruence multiplied by a unit.

There is an equally sharp global formulation. Let

\[
L_{n+1}=\operatorname{lcm}(1,\ldots,n+1)^3
\]

and replace only `b_n` by `R_nX` in the Casoratian. The cleared defect is

\[
\begin{aligned}
E_n(X)
&=L_{n+1}\left(
 a_nb_{n+1}-a_{n+1}R_nX+\frac6{(n+1)^3}
\right)\\
&=(L_{n+1}a_{n+1})R_n(A_n-X).
\end{aligned}
\tag{A.3}
\]

All quantities in the last line are integers. If
`A_n-X=R_n^kt`, then

\[
\boxed{
E_n(X)=(L_{n+1}a_{n+1})R_n^{k+1}t.
}
\tag{A.4}
\]

Therefore:

- if `X=A_n`, the defect is exactly zero;
- if `X!=A_n`, then `|E_n(X)|>=R_n^{k+1}`.

The divisibility and the Archimedean lower bound saturate at the same power.
Indeed, in the nonzero case,

\[
\log|E_n(X)|
=(k+1)\log R_n+\log|t|+\log|L_{n+1}a_{n+1}|
\ge (k+1)\log R_n.
\]

The standard denominator and recurrence-growth bounds give
`\log|L_{n+1}a_{n+1}|=O(n)`; they cannot reverse this inequality. The strict
inequality needed for a contradiction,

\[
0<|E_n(X)|<R_n^{k+1},
\]

is impossible from this construction.

### 4.5 Why recurrence rigidity and Padé uniqueness do not repair this

The integer `X` is a proposed value at one index, not a second solution of
the Apéry recurrence.

- Keeping `b_{n-1}` and `b_{n+1}` fixed while replacing `b_n` by `R_nX`
  normally breaks the recurrence.
- Requiring the modified value to satisfy the recurrence forces
  `R_nX=b_n`, hence `X=A_n`; this simply assumes the desired equality.
- Padé uniqueness concerns a whole solution or a polynomial approximant
  with prescribed degree and normalization. A single CRT representative is
  neither.
- Interlacing or real positivity does not constrain the `p`-adic quotient
  digit `(b_n/p) mod p^j`.

The repository already contains the structural versions of this no-go:

- `Q32_CODEX_RESUME_2026-07-23.md`, Section 13: consecutive Casoratian
  certificates are an invertible change of coordinates, with linear
  height;
- Section 15: every multi-shift Casoratian family has rank two;
- Section 143: a comparison section synchronized with a star of bad primes
  carries the full star radical in its determinant;
- Section 149: the rational horizontal endomorphism ring of the canonical
  recurrence is scalar.

Thus “more Casoratians” do not add independent rows. They either remain in
the same rank-two module or encode the target radical in their normalization.

### 4.6 Route A verdict and reopen condition

**Verdict:** the ordinary Casoratian route is exact but non-compressive.
It reconstructs the same local quotient coordinate and gives no
non-coincidence theorem.

A genuine reopening would require a second global observable `C_n` with all
of the following properties:

1. every target prime contributes a power `p^k` to `C_n`;
2. that divisibility is not an algebraic multiple of
   `X-A_n (mod p^k)` and not a rank-two change of coordinates;
3. `C_n` is proved nonzero;
4. `log|C_n|=o(k n)` or, more precisely, its Archimedean height is strictly
   below the accumulated target-prime valuation;
5. the construction is uniform for arbitrarily large fixed `k`.

The exact cubic `P(n)`, the canonical low-height gauge, or the
hypergeometric/Dwork realization would have to supply the independence.
The determinant line alone cannot.

## 5. Route B: literature and dimension audit

### 5.1 What Liu 2024 actually proves

There is a notation collision to avoid: these papers write `A_j` for the
ordinary Apéry number that this note calls `b_j`. In this section,
`\mathbf A_j` denotes the papers' Apéry number; the quotient
`A_n=b_n/R_n` retains its earlier meaning.

The relevant primary source is:

- Ji-Cai Liu, *An extension of Gauss congruences for Apéry numbers*,
  [arXiv:2404.16636](https://arxiv.org/abs/2404.16636), Theorem 1.1.

For `p>=5`, positive `N,ell`, and generalized-Apéry parameters
`r>=2`, `s,t>=0`, it proves

\[
\mathbf A_{Np^\ell}^{(r,s,t)}
\equiv \mathbf A_{Np^{\ell-1}}^{(r,s,t)}
+p^{3\ell}B_{p-3}\mathcal A_N^{(r,s,t)}
\pmod {p^{3\ell+1}}.
\tag{B.1}
\]

The earlier conjectures cited by Liu do not say more. Zhi-Hong Sun,
*Congruences for two types of Apéry-like sequences*,
[arXiv:2005.02081](https://arxiv.org/abs/2005.02081), Conjecture 5.1, gives
only special cases `\mathbf A_p,\mathbf A_{2p},\mathbf A_{3p} (mod p^4)`,
each involving
`B_{p-3}`. Liu's theorem is the uniform generalized-Apéry extension of this
kind of statement.

A later paper of Sun,
[arXiv:2409.06544](https://arxiv.org/abs/2409.06544), both quotes (B.1) and
states the stronger Conjecture 2.4:

\[
\mathbf A_{Mp^\ell}-\mathbf A_{Mp^{\ell-1}}
\equiv
2C_M p^{3\ell}
\left(
\frac{B_{2p-4}}{2p-4}
-2\frac{B_{p-3}}{p-3}
\right)
\pmod {p^{3\ell+2}}.
\tag{B.2}
\]

Here `p>5` and `M,ell` are positive. At `ell=1`, (B.2) is indeed a
conjectural congruence modulo `p^5`. It still concerns multiplicative
indices, not `p+r`; all `M` share the same one displayed prime-local scalar
combination; and it gives neither an arbitrary-order tower nor terms
`B_{p-5},B_{p-7},...`.

The differences from the source note are decisive:

- (B.1) is at multiplicative indices `Np^ell`, not shifted indices `p+r`;
- it is one digit deeper than the order-three Gauss congruence;
- for `ell=1`, its modulus is `p^4`;
- its only Bernoulli coordinate is `B_{p-3}`;
- the correction is additive with an explicit coefficient
  `mathcal A_n`, not the asserted general multiplicative formula
  `b_{mp}=b_m(1+C_m Delta_p) (mod p^5)`;
- it provides no formulas involving `B_{p-5},B_{p-7},...`.

No source for the claimed all-order tower or for “Q2672” was present in the
canonical repository or the temporary research directory.

### 5.2 Honest status of the requested level table

From the cited theorem, the only defensible table is:

| level | unconditional information | conjectural information relevant here |
|---|---|---|
| `p^3` | Gauss congruence `\mathbf A_{Np}≡\mathbf A_N`; 0 new coordinates | none needed |
| `p^4` | Liu (B.1); at most 1 coordinate, `B_{p-3}` | Sun 2020 special cases |
| `p^5` | no formula from Liu; dimension unknown | Sun (B.2), at most one displayed scalar combination |
| `p^7` | no cited formula; dimension unknown | no tower supplied |
| `p^9` | no cited formula; dimension unknown | no tower supplied |

Therefore the requested unconditional dimensions at `p^5,p^7,p^9` have not
been computed. They first require actual congruence theorems defining a
defect module. Even assuming (B.2) would specify only one
`\mathbf Z/p^2\mathbf Z`-valued scalar combination—at most one new
`p`-adic digit beyond (B.1)—not the free Bernoulli algebra used by the
temporary script. Turning that combination into an `F_p`-dimension statement
would itself require its Kummer and higher-digit relations.

### 5.3 Internal errors in `defect_dimension_count.py`

The temporary script assumes, without derivation, that the defect algebra is
the free polynomial algebra on symbols

\[
X_j=B_{p-(2j+1)},\qquad \operatorname{wt}(X_j)=j.
\]

Even under this unproved model, its table is inconsistent.

The number of monomials of **exact** weight `w` is the partition number
`p(w)`. If all weights up to `j` are retained, the nonconstant cumulative
dimension is

\[
D(j)=\sum_{w=1}^j p(w).
\]

Thus the speculative counts would be

| putative level | maximum weight | exact new monomials | cumulative nonconstant monomials |
|---|---:|---:|---:|
| `p^3` | 0 | 0 | 0 |
| `p^5` | 1 | 1 | 1 |
| `p^7` | 2 | 2 | 3 |
| `p^9` | 3 | 3 | 6 |

The script instead:

- says in prose that the `p^3` dimension is zero but prints one;
- calls `p(w)` a cumulative dimension although it counts exact weight;
- obtains “new dimension” by subtracting consecutive partition numbers,
  which has no meaning for either grading;
- counts `m=0` as an equation although Liu's theorem assumes positive `n`,
  and the `n=0` identity would be tautological;
- doubles every row by a “reflected endpoint” without proving that the
  reflected row is independent;
- invokes “Bernoulli independence at different primes”, for which no such
  theorem is supplied or known in the required form.

The output is consequently not evidence for rank growth.

### 5.4 Why equations from varying `m` do not overdetermine a fixed `A_n`

There are four separate quantifier failures.

#### (i) Every interpretation of “varying `m`” changes the sequence index

Varying the exponent `ell` in (B.1) concerns

\[
\mathbf A_{Np},\mathbf A_{Np^2},\mathbf A_{Np^3},\ldots.
\]

The temporary dimension script instead varies the base multiplier
`M=0,1,2,...`. That produces values
`\mathbf A_{Mp^\ell}` at different outer indices (and `M=0` is
tautological). Neither variation supplies several equations for one fixed
Apéry value.

For the moving-prime problem, the outer index is fixed and a top-half bad
prime has

\[
n=p+r.
\]

Changing either the exponent or the base multiplier does not give additional
equations for the same integer `A_n=b_n/R_n`.

#### (ii) Each equation has new left-side data

If the values `\mathbf A_{Np^\ell}` or `\mathbf A_{Mp^\ell}` are treated as
unknowns, every new index adds a new unknown along with the equation. If they
are treated as known Apéry values, the equations are already consistent
identities. More rows can determine or verify `B_{p-3}`; “more equations
than Bernoulli symbols” does not imply a contradiction.

An overdetermined linear system is contradictory only after one proves that
its augmented right-hand side is outside the coefficient column space.
Dimension counting alone says nothing of the kind.

#### (iii) Different primes have different local coordinates

For each `p`, the residues

\[
B_{p-3},B_{p-5},\ldots\pmod p
\]

live in a different residue field and form different local data. For a set
of bad primes `p|R_n`, the natural defect space is a direct product of the
`p`-local spaces. The number of local unknowns therefore grows with the
number of primes at the same rate as the number of local congruence rows.

To turn them into shared variables one needs a proved cross-prime
interpolation, motive, or low-height rational relation. That is precisely
the missing horizontal theorem; it cannot be supplied by counting.

#### (iv) Reflection does not automatically double rank

Modulo `p`, the direct and reflected Apéry values satisfy

\[
b_{p-1-r}\equiv b_r\pmod p,
\]

so they are the same condition. At the next digit the exact repository
calculation is

\[
b_{p-1-r}=b_r-2pW_r\pmod {p^2},
\]

where `W_r` is a new harmonic coordinate. A target condition `p|b_r` does
not force `p|W_r`. Thus reflection either duplicates the first row or
introduces a new local unknown; it does not provide a free independent
equation.

### 5.5 No Bernoulli-independence contradiction

Even a rigorously derived polynomial relation among
`B_{p-3},B_{p-5},... (mod p)` would require a theorem showing that the
relation is nonzero for the relevant primes. Statements such as “false for
almost all primes” are insufficient for a pointwise theorem: a hypothetical
large `R_n` could be supported on the exceptional primes unless their total
logarithmic weight along the moving diagonal is itself proved to be `o(n)`.
That is another form of the original problem.

There is also no general algebraic-independence theorem for these varying
Bernoulli residues across varying finite fields that could be inserted at
this point. Kummer congruences and finite harmonic identities create
dependencies rather than generic independent coordinates.

### 5.6 The inverse-limit no-go

The logical obstruction can be stated independently of Apéry numbers.

**Lemma (diagonal tautology).** Let `R>=2` and `A` be a nonnegative integer.
For every `k>=1`, let `x_k` be the least nonnegative representative of
`A mod R^k`. Then:

1. `x_{k+1}≡x_k (mod R^k)`;
2. `(x_k)_k` is the diagonal image of `A` in
   `lim_k Z/R^kZ`;
3. `x_k=A` for every `k` with `R^k>A`.

The proof is immediate from the definition.

Apply this with `R=R_n` and `A=A_n`. If a proposed endpoint tower is proved
to compute `A_n mod R_n^k` at every level, its compatible inverse-limit point
is diagonal and its canonical lifts stabilize. It is impossible to prove
non-diagonality from the same congruences.

Writing

\[
d_k=\frac{x_{k+1}-x_k}{R^k}
\]

only produces the base-`R` digits of `A`. Showing one digit `d_k` is nonzero
does not contradict eventual stabilization. One would need infinitely many
nonzero digits for the same fixed `A_n`, but an ordinary integer has only
finitely many.

If a symbolic Bernoulli tower instead defines a non-diagonal profinite
integer, the conclusion is that the symbolic tower is not the residue tower
of the actual `A_n`; it is not a contradiction to the size of `R_n`.

### 5.7 Route B verdict and reopen condition

**Verdict:** Route B has no defined defect module at the claimed levels, and
the proposed equation count has the wrong variables and quantifiers. Even a
correct all-order local expansion would reconstruct the actual integer rather
than prove its non-diagonality.

A genuine reopening needs all of the following.

1. A proved shifted expansion for `b_{qp+r}` or `b_{qp+r}/p` to arbitrary
   fixed order, uniformly in the moving residue `r`, not only a Gauss
   congruence at `np^m`.
2. A rigorously defined defect module after all harmonic, Kummer, reflection,
   and integrality relations are imposed.
3. An elimination identity for a **fixed outer `n`** whose rows do not add
   new Apéry-value unknowns as fast as they add equations.
4. A proof that the eliminated integer is nonzero.
5. A cross-prime height bound strong enough that its target-prime divisibility
   exceeds its Archimedean size.

Items 4 and 5 are the actual separation theorem. Rank growth without them is
only local reconstruction.

## 6. The V1 defect identity does not create a growing-rank tower

Let `mathcal P_p(z)` be the truncated Apéry polynomial. The proved V1 identity
over `F_p[z]` is

\[
(z+1)^3\mathcal P_p(z+1)-P(z)\mathcal P_p(z)
 +z^3\mathcal P_p(z-1)
=-16(2z+1)(z^p-z)^2.
\tag{V1}
\]

At a residue `a in F_p`, put `z=a+t`. Then

\[
z^p-z=t^p-t
\]

and the defect is

\[
-16(2a+1+2t)(t^p-t)^2.
\tag{V2}
\]

For Hasse orders below `p+1`, (V2) has only:

\[
[t^2]=-16(2a+1),
\qquad
[t^3]=-32.
\]

Thus:

- the value and first-derivative defects vanish;
- the next two local jets are explicitly forced;
- there is no sequence of new low-order independent Bernoulli coordinates;
- at the central residue `a=-1/2`, the vanishing order rises from two to
  three, but still by only one fixed order.

This explains why V1 is powerful for close-pair transport and the
`|Z_p|<<p^{2/3}` theorem: it gives homogeneous recurrence transport through
first order. It does not give higher `p`-adic divisibility of the normalized
quotient `b_n/p`.

To use V1 modulo `p^k`, one would first need an integral lift of its boundary
term. The factorization by `(z^p-z)^2` is a characteristic-`p` identity; its
higher lifts introduce Wilson quotients, harmonic sums, and normalized
Frobenius jets. Existing exact audits in
`Q32_CODEX_RESUME_2026-07-23.md` show:

- Section 14: the natural `q=1` Lucas defect normally has valuation exactly
  one;
- Section 20/Q376: after the forced reflected factor, the next quotient is
  Gessel's derivative defect, a new unconstrained coordinate;
- Section 118: mod-`p^2` reflection introduces the independent harmonic
  value `W_r`;
- Section 150/Q824: higher Frobenius jets saturate at homogeneous degree;
  one bad equation gives at most one selective prime factor per row degree.

Therefore repeated differentiation or exterior powers increase selective
valuation and algebraic/Archimedean degree together. A new algebraic relation
among the normalized singular jets would be required to break saturation.

## 7. Audit of the remaining temporary scripts

### `crt_tower_digits.py`

This script:

- computes `A_n` exactly first;
- prints the tautological residues `A_n mod R_n^k`;
- never derives those residues from Liu/Sun endpoint formulas;
- uses only top-half primes despite claiming all middle primes;
- correctly observes that no distinct congruent integer can lie within
  distance `<R_n^k` of `A_n`.

It is useful as a base-`R_n` digit display but supplies no separation
evidence.

### `verify_quotient_digit.py`

This script fails when it tries to invert `(r+1)^3` modulo `p^2` at the
singular step `r=p-1`:

```text
ValueError: base is not invertible for the given modulus
```

This is exactly why the scaled companion values and the block theorem are
needed.

### `verify_quotient_digit_v2.py`

The repaired script stays in exact integer arithmetic through the singular
step and then reduces modulo `p^2`. It successfully reproduces
`(b_n/p) mod p`, but this is not an independent endpoint formula: it has
already computed the exact Apéry sequence through the problematic index. It
verifies arithmetic implementation, not the claimed Bernoulli tower.

## 8. What the finite-level countermodel teaches

The source note asks what genuinely global property of the Apéry numbers is
absent from an artificial sequence that matches all finite local data while
prescribing the quotient.

It is **not** merely “more `p`-adic levels”. Any fixed integer already has a
compatible tower at every level, and CRT can impose arbitrarily many finite
local digits. The audited countermodels in the repository show that the
following data are also insufficient by themselves:

- a rank-two recurrence module;
- the determinant/Casoratian law;
- a common initial vector;
- finitely many comparison solutions;
- generic interpolation, Smith, Fitting, or resultant data;
- correct first-order exponential growth.

The feature not preserved cheaply by the strongest countermodels is the
simultaneous combination of:

1. the exact canonical cubic `P(n)`;
2. the canonical bounded-height integral gauge and initial state `(1,5)`;
3. the hypergeometric/G-function or Dwork realization across all indices;
4. a global Archimedean height bound coupled to the same nonzero integer that
   accumulates the local valuations.

The fourth item is essential. The first three are structural sources from
which such an integer might be built, but they do not themselves constitute
a separation theorem. In the fixed-state lattice calculation of repository
Section 143, the unknown radical is exactly the second successive minimum;
changing gauge can realize the prescribed local star only by writing that
radical into the coefficient height.

Thus the true missing theorem is not “the CRT lift is non-diagonal”. It is
closer to:

> Construct an Apéry-specific nonzero integer `C_{n,k}` whose target-prime
> valuation is at least `k log R_n`, while
> `log|C_{n,k}|=o(k n)`, uniformly for arbitrarily large fixed `k`, and prove
> that its divisibility is not a tautological multiple of `A_n-x_{n,k}`.

No such integer is produced by Routes A or B.

## 9. Final assessment

### Route A

**Closed in its present form.** The exact local formula is (A.1), and the
global fake-value defect is (A.3). Both are saturated reformulations of
`X≡A_n`; neither yields non-coincidence, additional `p`-adic depth, or a
sublinear-height cross-prime certificate.

### Route B

**Not yet mathematically instantiated.** The cited unconditional result gives
only one extra Gauss digit involving `B_{p-3}` at multiplicative indices;
Sun's next digit remains conjectural. The claimed high-order endpoint tower
and its partition dimensions are unsupported; the temporary count is
internally incorrect; varying `m` does not preserve the fixed outer index;
and correct residue formulas force a diagonal inverse limit rather than
exclude it.

### Route C

**Impossible as stated.** The true residue tower of the integer `A_n` is
diagonal and eventually stabilizes by definition. A useful replacement would
need a second, independently constructed lift with a strict Archimedean
separation, not the canonical residue lift itself.

### Research implication

These verdicts leave the conjecture open. They narrow the viable frontier:
any further quotient-separation attempt must start with an explicit
Apéry-specific global nonzero certificate and its height ledger. More local
digits, more Casoratians, or symbolic Bernoulli dimension counts do not
address the missing cross-prime Archimedean inequality.

## 10. Follow-up audit: sparse zero fibers do not give pointwise separation

**Date:** 2026-07-29.

This section audits the subsequent computational observation that the
single-prime zero fibers

\[
 Z_p=\{0\le r<p:b_r\equiv0\pmod p\}
\]

are very small, while the top-half target set

\[
 T_n=\{p:n/2<p\le n,\ p\text{ prime},\ n-p\in Z_p\}
\tag{10.1}
\]

has at most three elements in the tested range. The proposed inference was

\[
 |Z_p|\le D\text{ uniformly}\quad\Longrightarrow\quad |T_n|=O(1),
\tag{10.2}
\]

which would make the top-half logarithmic weight \(O(\log n)\).

The data are real and striking, but (10.2) is false. It confuses a bound on
the degree of every **prime row** with a bound on the degree of an
**integer column**. Reflection, nonconsecutivity, the Casoratian, and the
currently known defect identities do not repair that change of quantifiers.
Moreover, the observed Poisson law predicts that both \(\max |Z_p|\) and
\(\max |T_n|\) grow slowly rather than remain absolutely bounded.

The useful outcome is a corrected target: a horizontal, cross-prime
collision estimate for the translated sets \(p+Z_p\). This is exactly the
kind of datum missing from Routes A and B.

### 10.1 Exact data audit

The audit used the independently validated binary scan

`problems/3.2/data_zp_pairs.bin`

and the scanner report

`problems/3.2/bn_bigscan_report.md`.

The file contains every pair \((p,r)\), for \(7\le p\le2,000,000\), with
\(p\) prime, \(0\le r<p\), and \(p\mid b_r\). Its SHA-256, recorded by the
scanner, is

```text
8746d0b400c1b669b001eae955c602908a10c9ee4cb3cac62c6676ea2ddd874d
```

An independent parser reproduced the following.

For \(p\le10,000\), among the 1,226 primes \(p\ge7\):

| \(|Z_p|\) | number of primes |
|---:|---:|
| 0 | 756 |
| 1 | 2 |
| 2 | 369 |
| 4 | 83 |
| 6 | 14 |
| 8 | 2 |

Thus the exact mean is \(1172/1226=0.955954\ldots\). The two omitted rows in
the originally quoted even-valued histogram are

\[
 Z_{11}=\{5\},\qquad Z_{3137}=\{1568\}.
\tag{10.3}
\]

They are both central fixed-point zeros. Consequently, “\(|Z_p|\) is always
even” is already false below \(10^4\); the correct parity statement is
proved in Section 10.2.

Under the scanner's convention \(p\ge7\), the exact top-half histogram for
\(1\le n\le10,000\) is:

| \(|T_n|\) | number of \(n\) |
|---:|---:|
| 0 | 9,204 |
| 1 | 769 |
| 2 | 26 |
| 3 | 1 |

For the literal all-prime definition (10.1), \(Z_5=\{1,3\}\) adds singleton
hits at \(n=6,8\), changing only the first two counts to 9,202 and 771.
Primes 2 and 3 add nothing. All maxima, double hits, triple hits, and
asymptotic statements below are unchanged.

The phrase “27 examples with two primes” is correct only if it means “at
least two”: there are 26 exact double hits and the one triple hit

\[
 n=321,\qquad T_{321}=\{179,193,211\}.
\]

At this \(n\),

\[
 \frac{\log\prod_{p\in T_n}p}{\log n}
 =2.737952929\ldots .
\]

The full scan through two million gives:

| \(|Z_p|\) | number of primes |
|---:|---:|
| 0 | 90,375 |
| 1 | 2 |
| 2 | 45,020 |
| 4 | 11,375 |
| 6 | 1,875 |
| 8 | 257 |
| 10 | 24 |
| 12 | 2 |

There are 148,930 primes and 149,112 zero pairs, so the mean is

\[
 1.001222050628\ldots .
\]

The maximum is now 12, at

\[
 p=159,977,\qquad p=1,823,963.
\]

Thus the earlier maximum 8 was a finite-range value, not a stable ceiling.
The assertion in `Q3.2_density_theorem.tex` that the maximum zero count is
16 through \(80,000\) is also false: the exact full scan gives maximum 8 in
that range. The value 16 occurs in a separate maximum-**fiber** computation
\(\max_a N_p(a)\); it is not the zero-fiber maximum.

For \(1\le n\le2,000,000\), again with the report's \(p\ge7\) convention,
the top-half histogram is:

| \(|T_n|\) | number of \(n\) |
|---:|---:|
| 0 | 1,896,672 |
| 1 | 100,670 |
| 2 | 2,605 |
| 3 | 53 |

Including \(p=5\) changes only the first two counts to 1,896,670 and
100,672.

There is no fourfold hit in this range. The largest logarithmic-radical
exponent is

\[
 2.960803420\ldots
\]

at

\[
 n=1,689,988,\qquad
 T_n=\{1,253,839,\ 1,328,449,\ 1,651,781\}.
\]

This confirms the quoted finite-range phenomenon while extending it by a
factor 200.

### 10.2 The parity is elementary reflection, not a new K3 signal

For \(m<p\), write

\[
 \binom zm\binom{z+m}{m}
 =\frac{1}{(m!)^2}
   \prod_{j=0}^{m-1}(z-j)(z+j+1)
 =\frac{1}{(m!)^2}
   \prod_{j=0}^{m-1}\bigl(z(z+1)-j(j+1)\bigr).
\tag{10.4}
\]

Define

\[
 H_p(X)=
 \sum_{m=0}^{(p-1)/2}
 \frac{\prod_{j=0}^{m-1}(X-j(j+1))^2}{(m!)^4}
 \in\mathbf F_p[X].
\tag{10.5}
\]

For every \(r\in\{0,\ldots,p-1\}\), terms with

\[
 m>\min(r,p-1-r)
\]

vanish at \(X=r(r+1)\). The remaining terms are exactly the defining
Apéry sum, so

\[
 b_r=H_p(r(r+1))\quad\text{in }\mathbf F_p.
\tag{10.6}
\]

Since

\[
 (-1-r)(-r)=r(r+1),
\]

(10.6) immediately gives

\[
 b_{p-1-r}\equiv b_r\pmod p.
\tag{10.7}
\]

Notice that (10.5) has degree \(p-1\) in \(X\), not bounded degree: the
top summand has degree \(2\cdot(p-1)/2=p-1\) and nonzero leading
coefficient.

The involution \(r\mapsto p-1-r\) has the unique fixed point

\[
 c_p=(p-1)/2.
\]

Therefore the exact parity law is

\[
 |Z_p|\equiv
 \mathbf 1_{\{b_{(p-1)/2}\equiv0\pmod p\}}
 \pmod2.
\tag{10.8}
\]

Ahlgren--Ono identify the central value with the \(p\)-th coefficient
\(\gamma(p)\) of

\[
 \eta(2\tau)^4\eta(4\tau)^4
\]

(in fact via a stronger supercongruence). Thus the fixed weight-four form
controls precisely the possible unpaired root. It says nothing about the
off-center reflection pairs which account for essentially all of the data.

There is an additional pointwise reason that the central slice is harmless.
If \(n=qp+r\) and \(r=(p-1)/2\), then

\[
 2n+1=(2q+1)p.
\tag{10.9}
\]

Hence all central bad primes at a fixed \(n\) divide the single integer
\(2n+1\), and their total logarithmic weight is at most

\[
 \log(2n+1)=O(\log n).
\tag{10.10}
\]

For the top-half arc \(q=1\), (10.9) is \(3p=2n+1\), so there is at most one
central candidate. Consequently, a theorem about the nonordinary primes of
the fixed eta-product cannot be the missing pointwise theorem: that entire
fixed-form slice is already negligible. The residual problem consists of
the moving off-center roots.

### 10.3 The row-versus-column counterexample

The logical failure of (10.2) is absolute and does not depend on a weakness
of the known bound for \(Z_p\).

Fix a large integer \(n_0\). For every prime

\[
 n_0/2<p\le n_0,
\]

put

\[
 r_p=n_0-p
\]

and define an artificial reflected zero set

\[
 S_p=\{r_p,\ p-1-r_p\}.
\tag{10.11}
\]

Then:

1. \(|S_p|\le2\) for every \(p\);
2. \(S_p\) is invariant under \(r\mapsto p-1-r\);
3. \(n_0-p=r_p\in S_p\) for every prime in the interval.

Thus the column over \(n_0\) contains every prime in
\((n_0/2,n_0]\), even though every row has size at most two. Its logarithmic
weight is

\[
 \sum_{n_0/2<p\le n_0}\log p
 \sim n_0/2.
\tag{10.12}
\]

The possible coincidences in which the two entries of (10.11) collapse, or
become consecutive, occur for at most \(O(1)\) primes. Omitting those primes
preserves (10.12). Thus the model can also obey the central exception rule
and nonconsecutivity, up to an immaterial bounded deletion.

This countermodel does not claim that the actual Apéry zero sets are
adversarial. It proves the precise logical statement:

> No theorem whose hypotheses see only \(|Z_p|\), reflection, central
> parity, and within-row spacing can bound \(|T_n|\) pointwise.

The reflection partner of an actual edge \((p,n=p+r)\) lies at

\[
 n'=p+(p-1-r)=3p-1-n.
\tag{10.13}
\]

It is a second edge in a generally different integer column. Formula
(10.13) explains why reflection doubles row counts but does not constrain
the number of distinct primes meeting one fixed \(n\).

In bipartite-graph language, a uniform degree bound on the prime side gives
no maximum-degree bound on the integer side. The proposed deduction was the
forbidden transposition of those two degrees.

### 10.4 The Poisson model predicts slow unbounded growth

The finite data do not merely fail to prove an absolute constant. Their
best-fitting model predicts that no such constant exists.

After removing the two central exceptions, put \(Y_p=|Z_p|/2\). The scan is
extremely close to

\[
 Y_p\sim\operatorname{Poisson}(1/2).
\tag{10.14}
\]

For 148,930 primes, (10.14) predicts:

| threshold | observed number | Poisson prediction |
|---:|---:|---:|
| \(|Z_p|\ge8\) | 283 | 260.87 |
| \(|Z_p|\ge10\) | 26 | 25.63 |
| \(|Z_p|\ge12\) | 2 | 2.11 |
| \(|Z_p|\ge14\) | 0 | 0.15 |

The agreement at the extreme tail is especially relevant: the two rows of
size 12 are almost exactly the predicted number. For independent
Poisson variables of fixed mean, the sample maximum grows like

\[
 \frac{\log M}{\log\log M}
\]

over \(M\) samples. In this model, \(\max_{p\le P}|Z_p|\) is therefore
unbounded (with an extra factor two from reflection), although it grows very
slowly.

For the top-half target count, a random sparse-row model gives

\[
 \lambda_n:=\mathbf E|T_n|
 \simeq\sum_{n/2<p\le n}\frac1p
 \simeq\frac{\log2}{\log n}.
\tag{10.15}
\]

Using the literal varying-mean Poisson mixture in (10.15), the expected
numbers of exceptional columns are:

| range | expected \(|T_n|\ge3\) | observed | expected \(|T_n|\ge4\) | observed |
|---:|---:|---:|---:|---:|
| \(n\le10,000\) | 1.251 | 1 | 0.055 | 0 |
| \(n\le2,000,000\) | 45.606 | 53 | 0.636 | 0 |

Under an independence heuristic, seeing no fourfold hit through two million
has probability approximately

\[
 e^{-0.636}=0.529\ldots,
\]

so it is not evidence for a hard ceiling at three. Numerical integration of
the same model makes the expected cumulative number of fourfold hits cross
one near \(N=4,000,000\). This is a heuristic prediction, not a theorem, but
it shows that a modest extension of the finite range could naturally
produce a fourfold hit.

More generally, (10.15) predicts a slowly unbounded top-half maximum,
roughly of order

\[
 \frac{\log N}{2\log\log N},
\]

not \(O(1)\). This growth is still vastly smaller than what Problem 3.2
requires. At the extreme-value scale it predicts

\[
 \max_{n\le N}\sum_{p\in T_n}\log p
 \asymp \frac{(\log N)^2}{2\log\log N},
\]

up to lower-order factors. Thus the finite-range law
\(\log R_n=O(\log n)\) is also not the natural asymptotic prediction, even
though the predicted polylogarithmic growth is still \(o(n)\).

### 10.5 A new full-middle scatter: maximum seven and Poisson mean \(\log2\)

The top-half statistic sees only the quotient \(q=1\). For the actual
middle-prime \(b\)-channel, define

\[
 K_b(n)=
 \#\{\sqrt n<p\le n:p\text{ prime},\ b_{n\bmod p}\equiv0\pmod p\}.
\tag{10.16}
\]

The saved pair file suffices to compute (10.16) exactly through two million.
For each recorded \((p,r)\), scatter it to

\[
 n=qp+r,\qquad q\ge1,\qquad n\le2,000,000,\qquad n<p^2.
\tag{10.17}
\]

The last inequality is exactly \(p>\sqrt n\).

The result is

\[
 \max_{n\le2,000,000}K_b(n)=7,
\tag{10.18}
\]

attained at 14 integers, beginning with

```text
208253, 501092, 513596, 855688, 939041, 985143, 1097287, ...
```

For example, the seven hits at \(n=208,253\) are

```text
(p,r,q) =
(499,170,417), (1381,1103,150), (6949,6732,29),
(22697,3980,9), (47017,20185,4), (54629,44366,3),
(195319,12934,1).
```

The histogram on the final million integers is:

| \(K_b(n)\) | count for \(1,000,001\le n\le2,000,000\) |
|---:|---:|
| 0 | 510,950 |
| 1 | 343,332 |
| 2 | 114,928 |
| 3 | 25,808 |
| 4 | 4,341 |
| 5 | 572 |
| 6 | 61 |
| 7 | 8 |

Its mean, variance, and variance-to-mean ratio are

\[
 0.671258,\qquad
 0.671072697436,\qquad
 0.999723947328.
\tag{10.19}
\]

The falling second moment on the same million-point interval is exactly

\[
 \sum_n(K_b(n))_2=450,402,
\qquad
 \frac1{10^6}\sum_n(K_b(n))_2=0.450402.
\]

The Poisson prediction is
\(0.671258^2=0.450587\ldots\). Thus the specific cross-prime collision
statistic isolated later in (10.23), not just the one-point histogram,
already matches the independence model numerically.

This is the predicted scale:

\[
 \sum_{\sqrt n<p\le n}\frac1p
 =\log\log n-\log\log\sqrt n+o(1)
 =\log2+o(1).
\tag{10.20}
\]

A Poisson variable of mean \(0.671258\), sampled one million times,
predicts 6.79 values at least seven and 0.56 values at least eight. The
observed counts are eight and zero. Thus the full-middle data agree with the
same Poisson picture down to the extreme tail.

This new scatter is a useful warning. The maximum has already grown from
five below \(10^4\), to six below \(10^5\), to seven below \(10^6\). The
natural prediction is

\[
 \max_{n\le N}K_b(n)\asymp
 \frac{\log N}{\log\log N},
\tag{10.21}
\]

not a uniform constant. Of course (10.21), if proved even with a much weaker
subpower upper bound, would solve the \(b\)-channel with enormous room.

The separate companion \(a_q\)-channel in the exact support law is not
included in (10.16). Existing divisor counting already gives that channel
sublinear logarithmic weight. The moving \(b_{n\bmod p}\)-channel remains
the pointwise obstruction.

### 10.6 Why K3/Frobenius does not supply the missing implication

There are three distinct arithmetic objects which must not be conflated.

1. **The fixed central modular form.** Ahlgren--Ono control
   \(b_{(p-1)/2}\). By (10.9), this slice is already pointwise harmless.
2. **The geometric Apéry/K3 family parameter.** Ordinary-reduction theorems
   for a fixed K3 surface vary the prime while holding the surface fixed.
3. **The coefficient-index interpolation (10.5).** Here \(p\) is fixed and
   the character/index \(r\) varies. The polynomial itself depends on \(p\)
   and has degree \(p-1\).

The distinction between the last two objects can be made exact. Put

\[
 \mathcal A_p(t)=\sum_{r=0}^{p-1}b_rt^r\in\mathbf F_p[t].
\]

This truncated generating series is the genuine first Hasse--Witt scalar
for the geometric Apéry K3 pencil (in the standard local trivialization).
For \(1\le r\le p-2\), finite Mellin inversion gives

\[
 b_r=-\sum_{t\in\mathbf F_p^\times}\mathcal A_p(t)t^{-r}.
\]

Thus \(b_r=0\) is the vanishing of one multiplicative Fourier coefficient
of the whole Hasse section. It is not the geometric nonordinary-fiber
condition \(\mathcal A_p(t_0)=0\).

Caruso--Fürnsinn--Vargas-Montoya--Zudilin, Theorem 2, prove the precise
factorization

\[
 \mathcal A_p(t)=
 \begin{cases}
  Q_p(t)^2,&p\equiv1,5,7,11\pmod {24},\\
  (t^2-34t+1)Q_p(t)^2,&p\equiv13,17,19,23\pmod {24}.
 \end{cases}
\]

This is strong structure in the **geometric generating variable** \(t\).
It gives quadratic convolution identities among the coefficients, but no
coefficient-support bound by itself: even a square such as
\((1+t^m)^2\) can have arbitrarily many zero coefficients. In particular,
the theorem neither bounds \(Z_p\) nor identifies (10.5) with a bounded
degree K3 Hasse polynomial in the spectral variable \(X\).

There is a stronger exact coefficient package, but it leads to a
convolution problem rather than a zero bound. Let

\[
 F(t)=\sum_{m\ge0}b_mt^m,\qquad D(t)=t^2-34t+1,
\]

and normalize \(Q_p(0)=1\). Uniqueness of a formal square root in odd
characteristic shows that \(Q_p\) is the appropriate truncation modulo
\(p\) of

\[
 S_+(t)=\sqrt{F(t)}
 \quad\text{or}\quad
 S_-(t)=\sqrt{F(t)/D(t)},
\]

according to the two CFVZ branches. Self-reciprocity of
\(\mathcal A_p\) and \(D\) also gives \(Q_p^*=\pm Q_p\).

Direct symmetric-square reduction of the Apéry differential operator gives
the independently checked equations

\[
\begin{aligned}
 4tD S_+''+4(2t^2-51t+1)S_+'+(t-10)S_+&=0,\\
 4tD S_-''+4(4t^2-85t+1)S_-'
   +3(3t-26)S_-&=0.
\end{aligned}
\tag{K3.1}
\]

Writing \(S_\pm=\sum s_m^\pm t^m\), \(s_{-1}^\pm=0\), these are

\[
\begin{aligned}
 4(m+1)^2s_{m+1}^+
 &-(136m^2+68m+10)s_m^+
 +(2m-1)^2s_{m-1}^+=0,\\
 4(m+1)^2s_{m+1}^-
 &-(136m^2+204m+78)s_m^-
 +(2m+1)^2s_{m-1}^-=0.
\end{aligned}
\tag{K3.2}
\]

Thus an Apéry coefficient is exactly a Cauchy-convolution coefficient

\[
 b_r=\sum_{i=0}^r s_i^+s_{r-i}^+,
\]

or, in the other branch, the fixed three-term \(D\)-filter of the analogous
convolution. This rank-two square-root period is genuine extra structure.
It is not yet transverse structure: the convolution map is triangular,
with \(b_r=2s_r^++P_r(s_1^+,\ldots,s_{r-1}^+)\), and the third-order Apéry
operator is precisely the symmetric square of (K3.1). Substitution therefore
recovers the original recurrence rather than a second hit-specific equation.
A useful escape would require a new anti-cancellation or bounded-boundary
telescoping theorem for this canonical convolution; factorization and
reciprocity alone permit reciprocal squares with linearly many zero
coefficients.

The most natural bounded-boundary escape can in fact be ruled out.  After
the self-adjoint gauges

\[
 s_m^\varepsilon=
 (2m+1)^\varepsilon\frac{\binom{2m}{m}}{4^m}y_m^\varepsilon,
 \qquad \varepsilon\in\{0,1\},
\]

the two square-root recurrences have a genuine same-index
Christoffel--Darboux identity.  It sums
\(y_i(\lambda)y_i(\mu)\) at two spectral parameters.  The Apéry
convolution instead pairs the reversed indices \(s_i s_{r-i}\), so this
identity does not apply.

The obstruction is exact.  Put \(u_i=s_i\), \(v_i=s_{r-i}\) and try the
natural bilinear concomitant

\[
 {\cal W}_i=X_i u_{i+1}v_i-Y_i u_iv_{i+1}.
\]

Cancellation of its two cross terms forces

\[
 \frac{X_{i+1}}{X_{i-1}}
 =
 \left\{
 \frac{(i+2)(2r-2i+2\varepsilon-1)}
 {(2i+2\varepsilon+1)(r-i+1)}
 \right\}^{\!2}.
\tag{K3.3}
\]

No nonzero \(X(i,r)\in\overline{\mathbf Q(r)}(i)\) satisfies (K3.3).
For a rational shift quotient \(X(i+1)/X(i-1)\), the zero-minus-pole
multiplicity must balance on every translation-by-two orbit.  The right
side of (K3.3) has integer zero/pole locations \(-2,r+1\) and
half-integer locations \(r+\varepsilon-\tfrac12,
-\varepsilon-\tfrac12\); exactly one of the two pairs lies on the same
parity orbit.  The other orbit is unbalanced.  Thus even a rationally
weighted reversed Green identity already fails before one asks for unit
weight.

Solving (K3.3) separately on the even and odd subsequences does produce a
hypergeometric weight, but its Pochhammer products have length
\(\asymp r\).  Clearing its moving zeros and poles has degree
\(\asymp r\) and logarithmic height \(O(r\log r)\).  It merely encodes the
whole convolution and is not a bounded boundary certificate.

The obstruction extends beyond this two-term ansatz.  A uniformly
bounded-width rational adjacent-state telescoper would sum to an identity

\[
 Q(r)[t^r]S_\varepsilon(t)^2
 =P_0(r)s_r^\varepsilon+P_1(r)s_{r-1}^\varepsilon
\tag{K3.4}
\]

up to finitely many initial terms (and with the fixed \(D\)-filter in the
minus branch).  Generating functions would turn (K3.4) into a nonzero
quadratic--linear differential relation between
\(S_\varepsilon\) and \(\theta S_\varepsilon\).
The CFVZ algebraic pullback identifies this rank-two equation with the
Gauss equation for
\({}_2F_1(\tfrac13,\tfrac23;1;y)\).  Its irreducible monodromy contains
two nontrivial unipotents with distinct fixed lines, so its connected
differential Galois group is \({\rm SL}_2\).  Hence
\(S_\varepsilon,\theta S_\varepsilon\) are algebraically independent over
\(\overline{\mathbf Q(t)}\), contradicting such a relation.

Accordingly the minimal bounded outer telescoper for the observable
\(F=D^\varepsilon S_\varepsilon^2\) is its irreducible
three-dimensional symmetric-square operator.  In coefficient form it is
exactly the original Apéry recurrence.  The square-root factorization
therefore supplies neither a second zero condition nor a lower-degree
gap certificate.  What remains open is a batch anti-cancellation theorem
for actual convolution returns; it cannot be obtained by collapsing each
coefficient to \(O(1)\) boundary states.

Bogomolov--Zarhin prove density-one ordinary reduction, after a finite
extension, for a fixed K3 surface over a number field. That theorem has
neither of the quantifiers needed here:

- it does not bound the number of \(r\)-values in one field
  \(\mathbf F_p\);
- it does not control the moving specialization \(r=n\bmod p\) while both
  \(p\) and the character vary.

Likewise, generic ordinarity of a family does not make its nonordinary
divisor bounded-degree uniformly in \(p\). Classical Hasse invariants
already provide counterexamples to that inference. In the present problem
the mismatch is sharper: \(H_p\) is an interpolation in the coefficient
index, not the geometric Hasse invariant in the K3 pencil parameter.

Archimedean Frobenius equidistribution and Sato--Tate control fixed
continuous test functions of normalized traces. The exact event
\(b_r\equiv0\pmod p\) is a same-characteristic, \(p\)-adically shrinking
target. Passing from equidistribution to an exact residue-zero count would
require a local limit theorem uniform in the moving character. No cited
Katz/Deligne, ordinary-reduction, or Chebotarev result supplies such a
theorem.

Therefore K3/Frobenius currently explains why modular and hypergeometric
structures are present, but it proves neither:

\[
 |Z_p|=O(1)
\]

nor the genuinely needed horizontal statement

\[
 \#\{p:n\bmod p\in Z_p\}=o(n/\log n)
\quad\text{uniformly in }n.
\]

The reflection parity is already fully explained by (10.4), with no
algebraic geometry.

#### An unconditional divisor-sensitive Kummer-order pruning

The fixed-\(n\) relation does give one useful cross-prime reduction that is
independent of the zero-fiber bound. Write

\[
 n=qp+r,\qquad q=\lfloor n/p\rfloor,\qquad0\le r<p,
\]

and let

\[
 d_{p,n}=\operatorname{ord}(\omega_p^r)
 =\frac{p-1}{\gcd(p-1,r)}
 =\frac{p-1}{\gcd(p-1,n-q)}.
\tag{KO.1}
\]

The last equality follows from \(r\equiv n-q\pmod {p-1}\). For fixed
\(n,q,D\), one has the purely divisor-theoretic bound

\[
 \boxed{
 \#\left\{p:\left\lfloor\frac np\right\rfloor=q,\ d_{p,n}\le D\right\}
 \le
 \tau(n-q)\left(\frac Dq+2\right).
 }
\tag{KO.2}
\]

Indeed, put \(g=\gcd(p-1,r)\), \(p-1=dg\), and \(r=jg\). Then
\(0\le j\le d\), \(\gcd(d,j)=1\), and

\[
 n-q=(qd+j)g.
\]

Writing \(h=qd+j\), one has

\[
 h\mid n-q,\qquad qd\le h\le(q+1)d,\qquad
 p=1+d\frac{n-q}{h}.
\]

For each divisor \(h\mid n-q\), the possible integers \(d\le D\) lie in
\([h/(q+1),h/q]\), an interval containing at most \(D/q+2\) integers.
This proves (KO.2) without using either primality or the Apéry condition.
It sharpens the earlier \(O(QD\,n^{o(1)})\) order split by retaining the
factor \(1/q\).

Let

\[
 M(n)=\max_{m\le n}\tau(m),\qquad
 L_n=4M(n)(\log n)^2=n^{o(1)}.
\tag{KO.3}
\]

For \(p>\sqrt n\), if
\(\gcd(p-1,n-q)>L_n\), then

\[
 d_{p,n}<\frac{n}{qL_n}.
\]

Applying (KO.2) with \(D=n/(qL_n)\), summing over
\(q<\sqrt n\), and bounding every prime by \(\log n\), gives

\[
\begin{aligned}
 &\sum_{\substack{p>\sqrt n\\
   \gcd(p-1,n-\lfloor n/p\rfloor)>L_n}}\log p\\
 &\quad\le
 \sum_{q<\sqrt n}
 \left(\frac{n}{q^2L_n}+2\right)\tau(n-q)\log n\\
 &\quad\ll
 \frac{nM(n)\log n}{L_n}
 +\sqrt n\,M(n)\log n
 =O\!\left(\frac n{\log n}\right)+n^{1/2+o(1)}
 =o(n).
\end{aligned}
\tag{KO.4}
\]

This is an unconditional, all-\(n\) improvement over the older
bounded-order pruning. The full unresolved lower-digit channel may therefore
be restricted to

\[
 \gcd\!\left(p-1,n-\lfloor n/p\rfloor\right)\le L_n,
\qquad
 d_{p,n}\ge\frac{p-1}{L_n}=p^{1-o(1)}.
\]

The uniform threshold (KO.3) is convenient but not optimal. The same proof
allows a quotient-adaptive cutoff. Let \(F(n)\to\infty\) with
\(F(n)=n^{o(1)}\), and put

\[
 L_{n,q}=\tau(n-q)\log n\,F(n).
\tag{KO.5}
\]

For primes with \(\gcd(p-1,n-q)>L_{n,q}\), (KO.2) with
\(D=n/(qL_{n,q})\) gives total logarithmic weight

\[
\begin{aligned}
 &\sum_{q<\sqrt n}
 \left(\frac{n}{q^2L_{n,q}}+2\right)
 \tau(n-q)\log n\\
 &\qquad\le
 \frac n{F(n)}\sum_{q\ge1}\frac1{q^2}
 +2\log n\sum_{q<\sqrt n}\tau(n-q)\\
 &\qquad=O\!\left(\frac n{F(n)}\right)+n^{1/2+o(1)}
 =o(n).
\end{aligned}
\tag{KO.6}
\]

Thus the sharper canonical residual is

\[
 \gcd(p-1,n-q)\le\tau(n-q)\log n\,F(n),
\qquad
 d_{p,n}\ge
 \frac{p-1}{\tau(n-q)\log n\,F(n)}
 =p^{1-o(1)}.
\tag{KO.7}
\]

For the top half this becomes
\(\gcd(p-1,n-1)\le\tau(n-1)\log n\,F(n)\). The reduction is sharp in the
sense relevant here: small kernels appear frequently in the computed target
set, and the sieve gives no estimate for their Apéry-zero condition. It
removes bounded- and moderately growing-order motives; it does not replace
the horizontal collision theorem.

The exact scan confirms that the hard sector is already populated. Among
106,039 top-half hits with \(n\le2,000,000\), 39,804 have
\(\gcd(p-1,n-1)=1\). Of the 53 triple-hit columns, 12 have kernel 1 at all
three primes, and 21 have the same kernel at all three primes. For example,

\[
 n=11576,\qquad p=8893,9319,11437
\]

is a triple with \(\gcd(p-1,n-1)=1\) for every \(p\). Thus the
nearly-primitive restriction is a genuine theorem-level pruning, not an
explanation of the observed low column multiplicity.

### 10.7 Consequences for Routes A and B

The new zero-fiber data do not reopen either audited route.

There is a stronger exact saturation statement. Let
\(\Omega_p=\{z(z+1):z\in Z_p\}\), so
\(q_p=|\Omega_p|\) is the number of reflection orbits. For a fixed target
column \(n\), let \(\mathcal L_d(n)\) be the lattice of integral
polynomials \(F(Y)\) of degree at most \(d\) which vanish modulo each target
prime \(p\) at every point of \(\Omega_p\). Since
\((n-p)(n-p+1)\equiv n(n+1)\pmod p\), every such polynomial satisfies

\[
 \prod_{p\in T_n}p\mid F(n(n+1)).
\]

For \(d\ge\max_p(q_p-1)\), interpolation and CRT give the exact local/global
indices

\[
 [\mathbf Z^{d+1}:\mathcal L_d(n)]
 =\prod_{p\in T_n}p^{q_p},
\qquad
 [(Y-N)\mathbf Z[Y]_{\le d-1}:
   \mathcal L_d(n)\cap(Y-N)\mathbf Z[Y]_{\le d-1}]
 =\prod_{p\in T_n}p^{q_p-1},
\tag{S.1}
\]

where \(N=n(n+1)\). Their transverse quotient is exactly
\(\prod_{p\in T_n}p\). Thus the entire zero fiber contributes \(q_p\)
conditions, but the fake-evaluation subspace \(F(N)=0\) consumes
\(q_p-1\) of them. The observed values \(q_p\le6\) through two million
cancel completely from the only quotient capable of producing a nonzero
global evaluation.

#### Route A: Casoratian

The Casoratian is a same-prime, same-recurrence determinant. It can prove:

- nonvanishing of the local solution pair;
- absence of simultaneous local zeros;
- transport and close-pair certificates inside one \(Z_p\);
- saturation of a locally imposed fake quotient.

More precisely, use integral recurrence rows \(L_m\) which evaluate an
arbitrary initial state at \(m\). If \(p>b\) and \(a,b\in Z_p\), then
\(L_a\) and \(L_b\) are nonzero rows with the same kernel, namely the line
spanned by the Apéry initial state. Hence all rows indexed by \(Z_p\) have
rank exactly one over \(\mathbf F_p\), regardless of whether the fiber has
2, 6, or 12 elements. Every nonzero exterior construction therefore
reduces to a \(2\times2\) row minor; all higher homogeneous exterior powers
vanish.

The forced reflected minor is also saturated. If the two indices are
\(t,t+h\) with \(2t+h+1=p\), its gap continuant has the universal factor
\(2t+h+1\), whether or not \(b_t\) vanishes. Dividing by this factor exposes
the first reflection/derivative defect, not a second zero condition.
Consequently a one-orbit fiber supplies no selective Casoratian divisor
after saturation.

Even a hypothetical proof \(|Z_p|\le2\) obtained from the Casoratian would
not compare the selected root at \(p\) with the selected root at \(q\).
The adversarial section (10.11) obeys the same row bound and reflection.

For a top-half target prime, Lucas gives

\[
 b_n\equiv b_1b_{n-p}=5b_{n-p}\equiv0\pmod p.
\tag{10.22}
\]

Thus the target radical divides the fixed integer \(b_n\). This is real
global coupling, but \(\log b_n=\Theta(n)\), exactly the saturated height
scale. The Casoratian formulas audited in Section 4 reconstruct the same
local quotient once its residue is imposed; they do not construct a
second nonzero integer of \(o(n)\) logarithmic height.

There is one exact narrow reopen condition. A non-reflected simultaneous
zero pair \(a,b\in Z_p\) forces \(p\) into the corresponding gap continuant
after all universal reflection factors have been divided out. Hence a
prime-independent collection of such pairs covering every target prime
would give one positive integer divisible by the target radical, with
logarithmic height

\[
 O\!\left(\log n\sum_{\{a,b\}}|a-b|\right).
\]

To make this height \(o(n)\), however, the total gap budget must already be
\(o(n/\log n)\). Minimal one-orbit fibers provide no non-reflected pair,
and bounded fiber size gives neither a cover nor reuse of one continuant by
many primes. This is a precise collision-cover reformulation, not a height
gain from the Casoratian.

The capacity ledger is sharp. A saturated continuant with gap \(h\) has
degree \(3(h-1)\), less only the one universal reflection factor when
applicable, and

\[
 |\Sigma_n(a,h)|\log(n/2)
 \le \log N_h^\sharp(a)
 \le(3(h-1)+O(1))\log n.
\]

It can therefore cover at most \((3+o(1))(h-1)\) target primes. Any cover
of \(K\) targets already has total gap \(\Omega(K)\); the determinant
product is not an amplification mechanism.

#### Quantitative closure of gap-continuant factorial moments

The continuant obstruction can be made quantitative.  This also records the
strongest new result recovered from the uncommitted `uisai2` notes Q910 and
Q913.

Fix the direct prefix window

\[
 n=3H+1,\qquad
 2H<p\le3H+1,\qquad
 Z_p(H)=Z_p\cap[0,H],\qquad z_p(H)=|Z_p(H)|.
\]

For \(h\ge1\), let

\[
 C_p(h)=\#\{0\le s\le H-h:s,s+h\in Z_p(H)\}.
\]

The cleared recurrence gives a positive gap continuant \(N_h(s)\), and two
zero endpoints imply \(p\mid N_h(s)\).  With

\[
 \Lambda_H=\log118+6\log(H+1),
\]

the elementary tridiagonal recurrence gives

\[
 0<N_h(s)\le
 \exp\bigl((h-1)\Lambda_H\bigr).
\tag{CM.1}
\]

Consequently, for every \(L\le H+1\),

\[
 \begin{aligned}
 W_2(H,L)
 &:=
 \sum_{2H<p\le3H+1}
 \sum_{1\le h<L}C_p(h)\log p\\
 &\le
 \sum_{h=2}^{L-1}\sum_{s=0}^{H-h}\log N_h(s)
 \le \frac12HL^2\Lambda_H .
 \end{aligned}
\tag{CM.2}
\]

This is a weighted factorial-moment inequality, not a heuristic.  If
\(8\le B\le H\), choose

\[
 L=\left\lceil\frac{4(H+1)}B\right\rceil.
\]

Partition \([0,H]\) into at most \(B/2\) intervals of length at most \(L\).
Cauchy--Schwarz on the numbers of zeros in these intervals shows that
\(z_p(H)>B\) creates at least \(B/2\) pairs of zeros at distance \(<L\).
Thus (CM.2) gives the explicit tail bound

\[
 \boxed{
 \sum_{\substack{2H<p\le3H+1\\z_p(H)>B}}\log p
 \le
 \frac{25H(H+1)^2\{\log118+6\log(H+1)\}}{B^3}.
 }
\tag{CM.3}
\]

The same inequality remains valid after imposing the moving target condition
\(p\mid b_{3H+1-p}\), since that only restricts the prime set.  At the
available vertical scale \(B=C H^{2/3}\), however, (CM.3) is only

\[
 O_C(H\log H),
\tag{CM.4}
\]

whereas the desired conclusion is \(o(H)\).  It becomes sublinear only when

\[
 \frac{B^3}{H^2\log H}\longrightarrow\infty,
\tag{CM.5}
\]

which lies beyond the proved \(O(H^{2/3})\) row bound.

Higher moments do not remove the logarithm.  Multiplying the \(r-1\)
consecutive gap certificates over an \(r\)-tuple yields, for fixed \(r\ge2\),

\[
 \sum_{2H<p\le3H+1}(z_p(H))_r\log p
 \ll_r H^{r+1}\log H.
\tag{CM.6}
\]

At \(B=H^{2/3}\), ordinary Markov bounds worsen as \(r\) grows.  Localizing
all \(r\) zeros to short blocks is also minimized by the pair case: the
height pays for each possible gap location at exactly the rate at which
the forced multiplicity grows.  Hence (CM.3) closes the entire method that
multiplies individual gap continuants and estimates them by absolute
height.

The precise reopening statement is a joint compression theorem.  For
\(L_H\asymp H^{1/3}\), let \(c_p(H)\) count simultaneous zero pairs with
gap \(<L_H\).  It would suffice to construct one nonzero integer \(C_H\),
independent of \(p\), such that

\[
 v_p(C_H)\ge c_p(H)
 \quad(2H<p\le3H+1),\qquad
 \log|C_H|=o(H^{5/3}).
\tag{CM.7}
\]

Indeed, \(z_p(H)>H^{2/3}\) forces \(c_p(H)\gg H^{2/3}\), so (CM.7) would
turn the high-fiber tail into \(o(H)\).  The naive product of the individual
continuants has logarithmic height \(O(H^{5/3}\log H)\).  Thus the missing
object is a genuinely Apéry-specific common compression identity, not
another Casoratian or another factorial moment.

#### Cross-prime collisions: the degenerate part is now controlled

There is a sharper pointwise decomposition for the top-half arc. Order two
target primes as \(p<\ell\), and write

\[
 h=\ell-p,\qquad s=n-\ell.
\]

The collision conditions are exactly

\[
 \ell=p+h,\qquad
 p\mid b_{s+h},\qquad
 \ell\mid b_s,\qquad
 0\le s<s+h<p<\ell.
\tag{PC.1}
\]

Call the pair **degenerate** if

\[
 p\mid b_s\quad\hbox{or}\quad \ell\mid b_{s+h};
\tag{PC.2}
\]

otherwise call it **pure cross**. The terminology records the essential
field distinction. In a degenerate pair, one of the two characteristics
sees two zero endpoints. In a pure-cross pair, each characteristic sees
only the endpoint belonging to the other modulus.

Let \(M_h^{\rm deg}(n)\) count the degenerate target pairs with prime gap
\(h\). Then the same-prime continuant gives the new unconditional pointwise
estimate

\[
 \boxed{M_h^{\rm deg}(n)\ll h.}
\tag{PC.3}
\]

Indeed, if \(p\mid b_s,b_{s+h}\), then \(p\mid N_h(s)\). Since
\(s=n-h-p\), this implies

\[
 p\mid N_h(n-h).
\]

If instead \(\ell\mid b_s,b_{s+h}\), then \(s=n-\ell\) gives

\[
 \ell\mid N_h(n).
\]

Both continuants are positive nonzero integers and

\[
 \log\{N_h(n-h)N_h(n)\}\ll h\log n.
\tag{PC.4}
\]

Every target prime exceeds \(n/2\), so (PC.3) follows by counting their
distinct large prime divisors. In particular,

\[
 \sum_{h\le H}M_h^{\rm deg}(n)\ll H^2.
\tag{PC.5}
\]

On the pure-cross locus the same continuant is forced to be a unit in both
fields:

\[
 p\nmid N_h(s),\qquad \ell\nmid N_h(s).
\tag{PC.6}
\]

For \(\ell\), divisibility of \(N_h(s)\), together with
\(\ell\mid b_s\), would propagate the zero to \(b_{s+h}\). For \(p\),
divisibility of \(N_h(s)\), together with \(p\mid b_{s+h}\), would force
\(p\mid N_{h-1}(s+1)\); the continuant Dodgson identity then contradicts
the fact that all intervening indices are \(p\)-units. Thus the ordinary
Casoratian detector does not merely fail to prove divisibility in the hard
case: its relevant factor is provably nonzero modulo both moving primes.

Writing

\[
 b_{s+h}=pu,\qquad b_s=\ell v
\]

turns divisor switching into the exact equation

\[
 u b_s-vb_{s+h}=huv.
\tag{PC.7}
\]

The new co-divisors \(u,v\) can have exponential height. Equivalently, the
transfer identity leaves the two opposite quotient residues

\[
 b_{s+h}/p\pmod\ell,\qquad b_s/\ell\pmod p
\]

uncontrolled. This is the precise different-characteristic obstruction.
A crossed companion/Casoratian expression can be made divisible by
\(p\ell\), but its low-gap term is accompanied by a multiple of
\(b_sb_{s+h}\); its logarithmic height is \(\Theta(s+h)\), not
\(O(h\operatorname{polylog}n)\).

This decomposition yields a sufficient pointwise theorem strictly weaker
than the shell pair-energy hypothesis. Let \(A_H^\times(n)\) be the number
of **adjacent**, in the ordered target-prime set \(T_n\), pure-cross pairs
with gap at most \(H\). It is enough to prove that, for every fixed \(A>0\),

\[
 \boxed{
  \sup_{N<n\le2N}A_{A\log N}^{\times}(n)
  =o_A(N/\log N).
 }
\tag{PC.8}
\]

To see this, a target set of size \(\gg_\varepsilon n/\log n\) has, after
choosing \(A>2/\varepsilon\), all but \(O(n/(A\log n))\) of its adjacent
gaps at most \(A\log n\). Equation (PC.5) makes the degenerate short gaps
only \(O_A(\log^2n)\). The remaining
\(\gg_\varepsilon n/\log n\) adjacent gaps would be pure cross,
contradicting (PC.8).

An algebraic statement far stronger than necessary, but exactly targeted,
would be a pure-cross gap certificate: for every \(h\le A\log n\),
construct one nonzero integer \(D_{n,h}\), independent of the pair location
\(s\), such that every pure-cross gap-\(h\) pair satisfies

\[
 p(p+h)\mid D_{n,h},\qquad
 \log|D_{n,h}|\ll h(\log n)^C.
\tag{PC.9}
\]

The degenerate certificate \(N_h(n-h)N_h(n)\) has the required height but
misses every pure-cross pair by (PC.6). Conversely, all currently known
certificates which see (PC.1) either depend on \(s\), retain the full
exponential Apéry height, or insert the moving product \(p(p+h)\)
universally. Thus (PC.8), rather than the full second moment, is the
narrowest clean top-half horizontal target exposed so far.

There is also a hard uncovered component. If

\[
 Z_p=\{r_p,p-1-r_p\},
\]

then \(p\) has no non-reflected zero pair at all, so no saturated
collision-cover factor is guaranteed to contain it. In the scan through
two million, 45,020 of the 58,555 nonempty rows are such exact doublets,
about 76.9 percent of the active vertical population. This does not prove
the same proportion in a fixed target column, but it shows that the
primitive obstruction is structurally dominant, not exceptional.

Consequently a collision-cover proof must split into two genuinely new
horizontal assertions:

1. the total logarithmic weight of target primes with exact-doublet fibers
   is \(o(n)\);
2. the remaining rich-fiber targets admit a non-reflection cover with
   total gap \(o(n/\log n)\).

Neither assertion follows from the current zero-fiber histogram.

#### Route B: defect dimension

A bound \(|Z_p|\le D\) reduces the number of choices at each prime, but it
does not reduce the number of prime coordinates. The local defect spaces
still form a direct product over \(p\), and a section may choose one of the
\(D\) roots in every row. Varying the local level introduces new local
coordinates just as in Section 5.

Dimension growth becomes relevant only after producing one of:

1. a common rational coordinate in which roots for different primes are
   reductions of bounded-height global points;
2. a nonzero eliminated integer whose target-prime valuation accumulates;
3. a horizontal equidistribution theorem showing that the CRT section has
   no exceptional low-height lift.

None follows from the row dimension. Hence the new computation reinforces,
rather than changes, the earlier verdict: local rank is not global
separation.

The V1 jet tower has the same exact defect. At derivative order \(j\), once
all lower jets are fixed, the highest jet satisfies an inhomogeneous copy of
the original second-order Apéry operator. Its solution space is therefore
an affine torsor under a fresh two-dimensional homogeneous state. Each new
equation arrives with new local coordinates; reflection leaves one free
scalar at every order. Thus higher jets identify successive quotient and
reflection-defect coordinates but do not create growing codimension. This
is the structural reason the speculative dimension table in Route B cannot
be repaired merely by adding more formal defect levels.

#### An exact fixed-\(n\) carrier, but no new height saving

There is a useful deterministic reformulation of the top-half branch. Set

\[
 C_n=\binom n{\lfloor n/2\rfloor}.
\]

For a prime \(p\in(n/2,n]\), the numerator \(n!\) contains exactly one
factor \(p\). Both denominator arguments are below \(p\), except when
\(n=2p-1\), where \(\lceil n/2\rceil=p\) cancels it. Consequently

\[
 v_p(C_n)=
 \begin{cases}
  1,&n\ne2p-1,\\
  0,&n=2p-1.
 \end{cases}
\]

The boundary prime is never a target for \(p\ge7\), since its residue is
\(p-1\) and

\[
 b_{p-1}\equiv1\pmod {p^2}.
\]

The stronger congruence is elementary: for \(1\le k\le p-1\), the factor
\(\binom{p-1+k}{k}\) contains exactly one \(p\), so every nonzero-\(k\)
summand in \(b_{p-1}\) contains \(p^2\), while the \(k=0\) summand is 1.

For \(n>10\), Gessel--Lucas also gives

\[
 b_n\equiv b_1b_{n-p}=5b_{n-p}\pmod p.
\]

It follows in the stronger valuation-by-valuation form that, for every
prime \(p\in(n/2,n]\),

\[
 v_p\!\left(\gcd\!\left(b_n,\binom n{\lfloor n/2\rfloor}\right)\right)
 =
 \mathbf1_{\{p\mid b_{n-p}\}}.
\]

Thus the interval part of the gcd is already squarefree. With
\(\operatorname{rad}_{(n/2,n]}\) denoting the radical restricted to primes
in that interval,

\[
 \boxed{
  \prod_{p\in T_n}p
  =
  \operatorname{rad}_{(n/2,n]}
  \gcd\!\left(b_n,\binom n{\lfloor n/2\rfloor}\right).
 }
\]

This exact identity was already isolated in Section 37 of
`Q32_CODEX_RESUME_2026-07-23.md`; an independent recurrence computation
rechecked it for every \(11\le n\le1000\). It is conceptually cleaner than
the interpolation carriers, but it is saturated:

\[
 \log C_n=(\log2)n+O(\log n),
\]

and every nonboundary candidate prime divides \(C_n\) independently of the
Apéry-zero condition. In fact its interval part is exactly the squarefree
primorial over all nonboundary candidate primes before intersecting with
\(b_n\); the binomial factor has no local selectivity. A subexponential
bound for the displayed gcd would
solve the top-half branch, but no known gcd theorem for fixed-dimensional
\(S\)-unit or constant-coefficient recurrence orbits applies to the pair
consisting of an Apéry \(P\)-recursive term and a factorial ratio. The
carrier is therefore an exact restatement, not the missing second
Casoratian observable.

### 10.8 The corrected sufficient condition: cross-prime collision energy

There is a clean theorem-shaped replacement for (10.2).

On a dyadic shell, define

\[
 K_b(n)=
 \#\{\sqrt n<p\le n:n\bmod p\in Z_p\}
\]

as in (10.16), and its cross-prime falling second moment

\[
 \mathcal C(N)=
 \sum_{N<n\le2N}(K_b(n))_2,
\qquad (x)_2=x(x-1).
\tag{10.23}
\]

If

\[
 \mathcal C(N)=o\!\left(\frac{N^2}{\log^2N}\right),
\tag{10.24}
\]

then for every \(n\in(N,2N]\),

\[
 (K_b(n))_2\le\mathcal C(N),
\]

and therefore

\[
 \max_{N<n\le2N}K_b(n)=o(N/\log N).
\tag{10.25}
\]

Since every middle prime has logarithm \(O(\log N)\), (10.25) gives

\[
 \sum_{\substack{\sqrt n<p\le n\\n\bmod p\in Z_p}}\log p=o(N).
\tag{10.26}
\]

Together with the already controlled small-prime and companion channels,
(10.26) proves Problem 3.2: the Wronskian bound gives
\(v_p(G_n)\le6\) for \(p>\sqrt n\), so passing from the support radical to
the prime-power contribution costs only this absolute factor.

The probabilistic data predict the much stronger natural bound

\[
 \mathcal C(N)\ll N^{1+o(1)}.
\tag{10.27}
\]

Indeed, \(K_b(n)\) has measured mean and variance close to \(\log2\).
Bound (10.27) would give

\[
 K_b(n)\ll N^{1/2+o(1)}
\]

uniformly, already far more than needed.

The best immediate deterministic estimate from the proved vertical bound is
much larger. Dropping the moving endpoint restrictions,

\[
\begin{aligned}
 \sum_{N<n\le2N}K_b(n)
 &\le
 \sum_{\sqrt N<p\le2N}|Z_p|\left(\frac Np+1\right)\\
 &\ll \frac{N^{5/3}}{\log N}.
\end{aligned}
\]

Together with \(\max K_b(n)\ll N/\log N\), this gives only

\[
 \mathcal C(N)
 \ll\frac{N^{8/3}}{\log^2N},
\]

missing (10.24) by the factor \(N^{2/3}\). A black-box additive large sieve
is weaker: if
\(\mu_N=\sum_{\sqrt N<p\le2N}|Z_p|/p\), its natural output is

\[
 \mathcal C(N)\ll N\mu_N^2+N^2\mu_N
 \ll\frac{N^{8/3}}{\log N}.
\]

The \(N^2\) term reflects Farey spacing \(N^{-2}\) while the sampled
\(n\)-interval has length only \(N\). Even a hypothetical
\(|Z_p|=O(1)\) gives only \(\mu_N=O(1)\) and the black-box bound
\(\mathcal C(N)=O(N^2)\), still short of the required logarithmic
little-oh. The needed gain is therefore genuinely Apéry-specific, not a
formal consequence of ordinary additive large-sieve spacing.

Expanding (10.23) shows the exact arithmetic content. It counts solutions
with distinct primes \(p\ne q\) to

\[
 n=ap+r=bq+s,\qquad
 r\in Z_p,\quad s\in Z_q,\quad N<n\le2N,
\tag{10.28}
\]

under the middle-prime restrictions. For the top-half arc, \(a=b=1\), so
the core becomes the fixed-sum collision

\[
 p+r=q+s=n.
\tag{10.29}
\]

Equivalently, one needs additive-energy control for the family of translated
sets

\[
 p+Z_p\subset[p,2p).
\]

A row bound alone permits all these translated sets to pass through one
common \(n_0\), making \(\mathcal C(N)\) as large as the failure scale in
(10.24). Thus (10.24), unlike \(|Z_p|\le D\), explicitly excludes the
adversarial alignment.

For a dyadic prime block \(X<p\le2X\), put

\[
 K_X(m)=\#\{X<p\le2X:m\bmod p\in Z_p\}.
\]

On the full CRT representative range \(0\le m<4X^2\), the corresponding
higher-moment formulation is

\[
 \sum_{m<O(X^2)}(K_X(m))_k
 \ll X^{2+o(1)}\lambda_X^k,
\qquad
 \lambda_X=\sum_{X<p\le2X}\frac{|Z_p|}{p}.
\tag{10.30}
\]

With the proven \(|Z_p|\ll p^{2/3}\), any fixed \(k>6\) at the
independence scale in (10.30) gives a power-saving pointwise bound. If a
uniform constant row bound were also known, a much lower moment would
suffice. But in either case the horizontal moment, not the row bound, is
the decisive new input.

#### A weaker sufficient theorem after quotient localization

The all-middle estimate (10.24) is convenient but stronger than necessary.
For a fixed positive integer \(q\), define the quotient-arc count

\[
 K_q(n)=
 \#\left\{
  p\text{ prime}:\frac n{q+1}<p\le\frac nq,\quad n-qp\in Z_p
 \right\}
\]

There is an exact factorial carrier on every fixed quotient arc, extending
the central-binomial identity above. Define

\[
 C_{n,q}=\binom n{\lfloor n/(q+1)\rfloor}.
\]

If \(n=qp+r<(q+1)p\), \(p>\sqrt n\), and \(r\le p-2\), then
\[
 r<\left\lfloor\frac n{q+1}\right\rfloor<p.
\]
Kummer's theorem therefore gives exactly one base-\(p\) carry and

\[
 v_p(C_{n,q})=1.
\tag{10.31}
\]

At the only excluded boundary \(r=p-1\), there is no carry, but also
\(b_r\equiv1\pmod p\). On the other hand, Gessel--Lucas gives

\[
 b_n\equiv b_qb_r\pmod p.
\]

For fixed \(q\), all sufficiently large primes in the arc are larger than
the fixed integer \(b_q\). Hence, for all sufficiently large \(n\),

\[
 \prod_{\substack{n/(q+1)<p\le n/q\\b_{n-qp}\equiv0\pmod p}}p
 =
 \operatorname{rad}_{(n/(q+1),\,n/q]}
 \gcd(b_n,C_{n,q}).
\tag{10.32}
\]

This was independently checked in 197,070 admissible prime/index cases for
\(1\le q\le7\) and \(n\le2000\), after explicitly excluding
\(p\mid b_q\). The proof, not the computation, establishes (10.32).

The carrier height is

\[
 \log C_{n,q}
 =
 n\,h\!\left(\frac1{q+1}\right)+O(\log n),
\]

where \(h(x)=-x\log x-(1-x)\log(1-x)\). It has positive linear rate for
every fixed \(q\), but rate \(O((\log q)/q)\) as \(q\to\infty\). This is
another exact explanation for the quotient localization: growing-\(q\)
arcs are cheap, while every fixed-\(q\) gcd remains a genuinely new
horizontal theorem.

Now define the shell pair energy

\[
 \mathcal E_q(N)=
 \sum_{N<n\le2N}(K_q(n))_2.
\tag{10.33}
\]

Every collision counted here has the exact form

\[
 n=qp+r=q\ell+s,\qquad
 r\in Z_p,\quad s\in Z_\ell,
\]

or equivalently

\[
 r-s=q(\ell-p).
\tag{10.34}
\]

The following family of estimates suffices for the full lower channel:

\[
 \boxed{
  \text{for every fixed }q,\qquad
  \mathcal E_q(N)
  =o_q\!\left(\frac{N^2}{\log^2N}\right).
 }
\tag{10.35}
\]

Here the quantifiers are important: (10.35) is required on every
sufficiently large dyadic shell, separately for each fixed \(q\).

To prove sufficiency, suppose instead that along an infinite sequence

\[
 \sum_{\substack{\sqrt n<p\le n\\n\bmod p\in Z_p}}\log p
 \ge\varepsilon n.
\]

Choose \(\eta=\eta(\varepsilon)>0\) so small that Chebyshev's bound for
\(\vartheta(\eta n)\) makes all primes \(p\le\eta n\) contribute at most
\(\varepsilon n/2\), even without imposing the zero condition. The
remaining hits have \(q=\lfloor n/p\rfloor<1/\eta\), so they lie in only
finitely many fixed quotient arcs. They include
\(\gg_\varepsilon n/\log n\) distinct primes. By pigeonhole and then an
infinite subsequence, one fixed \(q\) satisfies

\[
 K_q(n)\gg_\varepsilon n/\log n.
\]

The single summand at that \(n\) then forces

\[
 \mathcal E_q(N)\gg_\varepsilon N^2/\log^2N
\]

on its dyadic shell, contradicting (10.35). In particular, the much rougher
bound

\[
 \mathcal E_q(N)\ll_q N^{2-\delta_q}
\]

for any \(\delta_q>0\) would be enough.

This argument does **not** assert
\(\mathcal C(N)=\sum_q\mathcal E_q(N)\); that identity is false because
\(\mathcal C(N)\) also contains cross-\(q\) pairs. Instead it gives a
different \(L^\infty\) route: first discard \(p\le\eta n\), then bound the
finitely many surviving quotient arcs separately. No uniformity in \(q\)
is required, since \(q\) is frozen only after \(\varepsilon\) and \(\eta\)
are fixed.

The exact same-\(q\) equation (10.34) does improve the completely
unrestricted CRT count: for fixed distinct \(p,\ell\), choosing
\(r\in Z_p\) determines \(s\) and \(n\), so there are at most
\(\min(|Z_p|,|Z_\ell|)\) collisions. With only
\(|Z_p|\ll p^{2/3}\), however, this yields merely

\[
 \mathcal E_q(N)\ll_q\frac{N^{8/3}}{\log^2N},
\]

again missing (10.35) by \(N^{2/3}\). Reflection changes only constants.
The saving must decide whether the one marked residue determined modulo
\(\ell\) is actually in \(Z_\ell\); present same-prime continuants do not.

The binary scan gives direct evidence for exactly (10.33), not just for
one-point sparsity. On \(1,000,000<n\le2,000,000\), with \(p\ge7\), the
first four fixed arcs give:

| \(q\) | \(N^{-1}\sum_nK_q(n)\) | \(N^{-1}\mathcal E_q(N)\) | square of the mean | \(\max K_q(n)\) |
|---:|---:|---:|---:|---:|
| 1 | 0.050072 | 0.002396 | 0.002507 | 3 |
| 2 | 0.030336 | 0.000946 | 0.000920 | 3 |
| 3 | 0.022023 | 0.000550 | 0.000485 | 3 |
| 4 | 0.017400 | 0.000328 | 0.000303 | 3 |

The natural prediction for fixed \(q\) is

\[
 \mathbf E K_q(n)
 \sim
 \log\!\left(
  \frac{\log(n/q)}{\log(n/(q+1))}
 \right)
 \sim\frac{\log(1+1/q)}{\log n},
\]

and hence \(\mathcal E_q(N)=N^{1+o(1)}/\log^2N\). Again this is evidence,
not the needed deterministic bound.

#### Moment-order bookkeeping

There are two different higher-moment thresholds. Formula (10.30) uses the
entire CRT representative range \(m<O(X^2)\). With
\(|Z_p|\ll p^{2/3}\), its interval length \(X^2\) explains the condition
\(k>6\).

For the actual dyadic \(n\)-shell, one may localize simultaneously in
\(N<n\le2N\) and \(X<p\le2X\). Put

\[
 K_{N,X}(n)=
 \#\{X<p\le2X:\sqrt n<p\le n,\ n\bmod p\in Z_p\}
\]

and

\[
 A_X=\sum_{X<p\le2X}\frac{|Z_p|}{p}
 \ll\frac{X^{2/3}}{\log X}.
\]

An independence-scale bound

\[
 \sum_{N<n\le2N}(K_{N,X}(n))_k
 \ll N^{1+o(1)}A_X^k
\tag{10.36}
\]

would give

\[
 \max K_X(n)
 \ll N^{1/k+o(1)}\frac{X^{2/3}}{\log X}.
\]

After restoring the logarithmic prime weight, the worst block is \(X\asymp
N\), with exponent \(2/3+1/k\). Thus \(k=4\) is the first integer moment
which closes using only the currently proved vertical exponent:

\[
 \frac23+\frac14=\frac{11}{12}<1.
\]

This fourth-moment route is stronger than the absolute no-spike criterion
(10.35). The hierarchy is:

1. fixed-\(q\) pair energy (10.35), the weakest clean pointwise condition;
2. a sharp first moment plus Poisson-scale pair dispersion;
3. with only \(|Z_p|\ll p^{2/3}\), the relative fourth moment (10.36).

### 10.9 Precise live route and corrected verdict

The strongest genuinely live formulation exposed by the computation is:

> Prove a power-saving upper bound for the number of cross-prime
> fixed-sum coincidences \(p+r=q+s\), with \(r\in Z_p\) and \(s\in Z_q\),
> or more generally prove (10.24)/(10.30) for the actual Apéry zero
> positions.

By (KO.6), this horizontal theorem only has to treat the nearly primitive
sector

\[
 \gcd(p-1,n-q)\le\tau(n-q)\log n\,F(n),
 \qquad q=\lfloor n/p\rfloor,
\]

for any chosen subpolynomial \(F(n)\to\infty\). Thus fixed-order modular
forms and bounded-order Kummer motives can be removed unconditionally before
any collision estimate is attempted.

This requires an Apéry-specific horizontal mechanism, for example:

- a reciprocity law turning \(p\mid b_{n-p}\) into divisibility of a
  bounded-height integer shared across different \(p\);
- a bounded-complexity algebraic parametrization of the off-center root
  position as \(p\) varies;
- a dispersion theorem for the moving character-index Hasse interpolation;
- the global sublinear-height certificate already isolated at the end of
  Section 8.

No such mechanism is presently supplied by the recurrence Casoratian,
finite-level defect dimensions, K3 generic ordinarity, or the central
eta-product.

The revised assessment is therefore:

1. the small zero fibers and target counts are high-quality evidence for
   the conjecture;
2. the claimed evenness has two central exceptions and is completely
   explained by elementary reflection;
3. the exact scan through two million raises the zero-fiber maximum from
   8 to 12 and the full-middle column maximum to 7;
4. the observed Poisson law predicts slow unbounded maxima, not absolute
   constants;
5. even a proved absolute bound \(|Z_p|\le D\) would strengthen averaged
   results but would **not** solve the pointwise problem without an
   additional horizontal collision estimate;
6. Routes A and B remain saturated for exactly the same reason as before:
   neither creates cross-prime low-height separation.
7. The divisor-sensitive order sieve is a genuine unconditional gain: it
   reduces the remaining Mellin problem to characters of order
   \(p^{1-o(1)}\), but leaves the same horizontal alignment obstruction in
   that sector.

This is a negative answer to the proposed shortcut, but a positive
sharpening of the frontier: the missing theorem is now the explicit
scale-zero collision estimate (10.24), not a vertical zero-fiber bound.

### 10.10 Sources checked for the geometry and congruences

- S. Ahlgren and K. Ono, *A Gaussian hypergeometric series evaluation and
  Apéry number congruences*, J. reine angew. Math. 518 (2000), 187--212,
  <https://doi.org/10.1515/crll.2000.004>.
- F. Bogomolov and Y. Zarhin, *Ordinary reduction of K3 surfaces*,
  <https://arxiv.org/abs/0902.1548>.
- A. Malik and A. Straub, *Divisibility properties of sporadic Apéry-like
  numbers*, <https://arxiv.org/abs/1508.00297>.
- X. Caruso, F. Fürnsinn, D. Vargas-Montoya, and W. Zudilin,
  *Galois groups of Apéry-like series modulo primes*, Theorems 1--2,
  <https://arxiv.org/abs/2510.23298>,
  <https://doi.org/10.1017/S0004972725100932>.
- N. Grieve and J. T.-Y. Wang, *Greatest common divisors with moving
  targets and consequences for linear recurrence sequences*, Theorem 1.2,
  <https://arxiv.org/abs/1902.09109>.

These sources support, respectively, the central modular congruence, the
fixed-K3 ordinary-reduction theorem, the Lucas/divisibility framework, the
square/quadratic-times-square factorization of the reduced generating
series, and the scope of moving-target gcd theorems. None states a uniform
bound for the coefficient-index zero fibers or a pointwise horizontal
collision theorem of the form (10.24) or (10.35).

### 10.11 Corrections required in earlier working notes

This audit supersedes several pointwise claims elsewhere in the research
directory.

1. Any table claiming that a bound on \(|Z_p|\) **alone**, even
   \(|Z_p|=O(1)\), implies the all-\(n\) version of Problem 3.2 is invalid.
   Such a bound supports first-moment and density statements, but the
   adversarial construction (10.11) disproves the pointwise implication.
2. In particular, the hierarchy in `Q3.2_paper_draft.md` which lists
   \(|Z_p|=O(p^{1-\delta})\), \(O(\log p)\), or \(O(1)\) as standalone
   pointwise sufficient conditions must be amended by a horizontal
   hypothesis such as (10.24) or (10.30).
3. The maximum-zero claim 16 through 80,000 in
   `Q3.2_density_theorem.tex`/`.md` must be replaced by 8. The number 16
   belongs to \(\max_aN_p(a)\), not \(N_p(0)\).
4. The statement in `Q32_MAIN_THEOREM_writeup.md` that the every-\(n\)
   residual is the density-zero problem for the **fixed** weight-four
   eta-product is not a valid reduction for the off-center roots. The
   later `Q32_CODEX_RESUME_2026-07-23.md` correctly repairs this: the fixed
   form controls only the central residue, whose pointwise weight is already
   \(O(\log n)\) by (10.9).
5. The abstract of `proof.tex` says that the Poisson model predicts
   \(\max_{p\le N}|Z_p|\sim2\log\log N\). This is not the extreme-value law
   for a fixed-mean Poisson sample. Since \(|Z_p|/2\) is modeled by
   \(\operatorname{Poisson}(1/2)\) over \(\pi(N)\) rows, the first-order
   prediction is
   \[
    \max_{p\le N}|Z_p|
    \sim\frac{2\log N}{\log\log N},
   \]
   up to the usual replacement of \(N\) by \(\pi(N)\) inside the
   logarithms.
6. Corollary `cor:divides` in `proof.tex` claims
   \(v_p(G_n)\le3\lfloor\log_p n\rfloor\) from
   \(G_n\mid d_n(a_nb_{n-1}-a_{n-1}b_n)\). The displayed divisibility is
   not justified: in the attempted Bézout combination the coefficient
   \(a_{n-1}\) is generally rational. The denominator-safe combination uses
   \(d_n^2W_{n-1}\) and proves
   \[
    v_p(G_n)\le6\lfloor\log_p n\rfloor-3v_p(n)
   \]
   for \(p\ge5\). Thus the safe middle-prime depth is 6. All density and
   pointwise-support reductions remain valid after multiplying their
   absolute constants by at most two.

The unconditional density theorem and the proved local
\(|Z_p|\ll p^{2/3}\) bound are not affected by these corrections. What is
removed is only the unsupported passage from vertical sparsity to a
pointwise horizontal bound.
