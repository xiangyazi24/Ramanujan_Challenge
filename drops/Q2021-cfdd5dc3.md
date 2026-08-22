ANSWER Q2021 cfdd5dc3

# Hostile audit: determinant-content carrier from Apéry linear forms

## Verdict

**NOT PASS as stated.** The real/Archimedean height argument can be made rigorous, including both parities, the small `epsilon_n` term, the content bound, and the advertised positive constant. The integrality and determinant signs are also correct.

The unsupported step is the strict-top prime assertion after division by the content `C_n`. For a strict-top prime

\[
p>\frac{n+1}{2},\qquad p\le n,\qquad p\ge7,
\]

write

\[
n=p+r,\qquad 0\le r<M=\lfloor n/2\rfloor.
\]

At a target `p | b_r`, put

\[
x=v_p(d_n a_n),\qquad y=v_p(d_n b_n).
\]

Then the exact local calculation is

\[
\boxed{v_p(C_n)=\min(x,y)}
\]

and

\[
\boxed{
p\mid Q_n\iff x<y.
\ }
\tag{V}
\]

Thus the claimed equivalence

\[
p\mid Q_n\iff p\mid b_{n-p}
\]

requires the additional all-target orientation theorem

\[
v_p(d_n a_n)<v_p(d_n b_n).
\tag{O}
\]

That theorem is not supplied by the claimed proof.

In fact the neighboring Casoratian reduces the only possible failure of (O) to one very sharp seam. Since `d_n/p^3` is a `p`-unit on the strict top block,

\[
x=v_p(p^3a_n),\qquad y=3+v_p(b_n).
\]

At a target one has

\[
\boxed{
p\nmid Q_n
\iff
v_p(b_n)=3\ \text{and}\ v_p(p^3a_n)\ge6.
}
\tag{E}
\]

Consequently the strongest theorem justified by the stated ingredients is the corrected theorem in Section 7 below. Proving that the exceptional seam (E) never occurs would upgrade it to the claimed support theorem. I do not infer that nonoccurrence from finite scans.

The height conclusion *does* survive:

\[
\boxed{
\liminf_{n\to\infty}\frac{\log|Q_n|}{n}
\ge
3\log3+\log2+\frac12\log(1-\rho)-3>0,
\qquad
\rho=(\sqrt2-1)^4.
}
\tag{H}
\]

So the theorem is not a total failure: the exponential-height obstruction is valid; only the target-support equivalence is missing one genuinely arithmetic `p`-adic exclusion.

---

# 1. Conventions, signs, and integrality

Let

\[
L_n=\operatorname{lcm}(1,\ldots,n),\qquad d_n=L_n^3,
\]

and write

\[
X_n=d_na_n,\qquad Y_n=d_nb_n.
\]

The standard Apéry denominator theorem gives

\[
L_m^3a_m\in\mathbf Z\qquad(m\ge0),
\]

so for every `j<=n`, both `d_n a_j` and `d_n b_j` are integers. Hence

\[
Z_{n,j}
=\det\bigl((X_n,Y_n),(a_j,b_j)\bigr)
=X_nb_j-Y_na_j
=d_n(a_nb_j-b_na_j)
\in\mathbf Z.
\]

The signs at the first two nodes are exactly

\[
\boxed{Z_{n,0}=X_n,}
\qquad
\boxed{Z_{n,1}=5X_n-6Y_n.}
\tag{1.1}
\]

Indeed `(a_0,b_0)=(0,1)` and `(a_1,b_1)=(6,5)`.

The cardinal weights are the ordinary Lagrange coefficients for the nodes `0,...,M-1` evaluated at `n`:

\[
\lambda_j(n)
=(-1)^{M-1-j}\binom nj\binom{n-j-1}{M-1-j}.
\tag{1.2}
\]

They are integers and satisfy

\[
\sum_{j=0}^{M-1}\lambda_j(n)=1.
\]

Since `C_n` is the gcd of the integer `Z_{n,j}`, it divides their integral linear combination `T_n`; hence

\[
Q_n=T_n/C_n\in\mathbf Z.
\]

No denominator problem occurs here.

---

# 2. Exact interpolation identity for the Apéry error

Put

\[
\epsilon_j=\zeta(3)b_j-a_j.
\]

The standard Beukers positive triple integral may be written

\[
\epsilon_n
=\frac12\iiint_{[0,1]^3}
\frac{[x(1-x)y(1-y)z(1-z)]^n}
{[1-(1-xy)z]^{n+1}}\,dx\,dy\,dz.
\tag{2.1}
\]

Pushing forward the positive measure

\[
\frac{dx\,dy\,dz}{2[1-(1-xy)z]}
\]

under

\[
t=
\frac{x(1-x)y(1-y)z(1-z)}{1-(1-xy)z}
\]

gives

\[
\boxed{
\epsilon_j=\int_0^\rho t^j\,d\mu(t),
\qquad
0\le t\le\rho,
\qquad
\rho=(\sqrt2-1)^4=17-12\sqrt2.
}
\tag{2.2}
\]

Thus, with

\[
\mu_k:=\int_0^\rho(1-t)^k\,d\mu(t)>0,
\]

one has

\[
\Delta^k\epsilon_0
=\int(t-1)^k\,d\mu(t)
=(-1)^k\mu_k.
\tag{2.3}
\]

Let

\[
A_M(n):=\sum_{j<M}\lambda_j(n)\epsilon_j.
\]

Because the cardinal interpolant through the first `M` values is the truncated Newton polynomial,

\[
A_M(n)=\sum_{k=0}^{M-1}(-1)^k\mu_k\binom nk.
\]

Therefore its exact remainder is

\[
R_n:=\epsilon_n-A_M(n)
=\sum_{k=M}^{n}(-1)^k\mu_k\binom nk.
\tag{2.4}
\]

This is the alternating tail that drives the height estimate.

---

# 3. Alternating-tail lower bound: even and odd parity

Set

\[
H_k:=\binom nk\mu_k>0.
\]

Since `0<=1-t<=1`,

\[
\mu_{k+1}\le\mu_k.
\]

Also

\[
\frac{\binom n{k+1}}{\binom nk}=\frac{n-k}{k+1}.
\]

For `k>=M=floor(n/2)`, this ratio is at most `1`; after the central step it is strictly below `1`. Hence `H_k` is nonincreasing on the entire tail. Pairing consecutive terms gives

\[
\boxed{
(-1)^MR_n
\ge H_M-H_{M+1}>0.
}
\tag{3.1}
\]

This finite-tail statement is valid in both parities. What changes is the lower bound for the first difference.

## 3.1 Even case `n=2M`

Here

\[
H_M-H_{M+1}
=\binom{2M}{M}
\left(\mu_M-\frac{M}{M+1}\mu_{M+1}\right).
\]

Using `mu_{M+1}<=mu_M`,

\[
H_M-H_{M+1}
\ge
\frac1{M+1}\binom{2M}{M}\mu_M.
\]

Since `(1-t)^M >= (1-rho)^M`,

\[
\mu_M\ge(1-\rho)^M\mu_0
=(1-\rho)^M\epsilon_0.
\]

Thus

\[
\boxed{
|R_{2M}|
\ge
\frac{\epsilon_0}{M+1}
\binom{2M}{M}(1-\rho)^M.
}
\tag{3.2}
\]

## 3.2 Odd case `n=2M+1`

At the central step the two binomial coefficients are equal, so one must *not* recycle the even estimate. Instead

\[
H_M-H_{M+1}
=\binom{2M+1}{M}(\mu_M-\mu_{M+1}).
\]

But

\[
\mu_M-\mu_{M+1}
=\int_0^\rho t(1-t)^M\,d\mu(t)
\ge(1-\rho)^M\int_0^\rho t\,d\mu(t).
\]

The last integral is

\[
\epsilon_1=5\zeta(3)-6>0.
\]

Hence

\[
\boxed{
|R_{2M+1}|
\ge
\epsilon_1\binom{2M+1}{M}(1-\rho)^M.
}
\tag{3.3}
\]

So the parity issue is real but reparable. Uniformly,

\[
\boxed{
\log|R_n|
\ge
\left(\log2+\frac12\log(1-\rho)\right)n
-O(\log n).
}
\tag{3.4}
\]

The candidate's exponential constant is therefore correct; a proof that used the even `1/(M+1)` estimate unchanged in odd parity would be wrong, but (3.3) supplies the required replacement.

---

# 4. The `epsilon_n B_M(n)` term is negligible

Let

\[
B_M(n):=\sum_{j<M}\lambda_j(n)b_j.
\]

Write

\[
c_k:=\Delta^k b_0.
\]

The current repository file

```text
problems/3.2/research/proofs/Q32_SIGNED_PADE_DEGREE_ONE_NO_GO.md
```

derives the exact recurrence

\[
\begin{aligned}
(k+4)^3c_{k+4}
={}&(2k+7)(15k^2+105k+184)c_{k+3}\\
&+(k+3)(95k^2+570k+864)c_{k+2}\\
&+48(k+2)(k+3)(2k+5)c_{k+1}\\
&+32(k+1)(k+2)(k+3)c_k.
\end{aligned}
\tag{4.1}
\]

With

\[
c_0=1,\quad c_1=4,\quad c_2=64,\quad c_3=1240,
\]

all coefficients on the right are positive for `k>=0`; hence

\[
\boxed{c_k>0\quad(k\ge0).}
\tag{4.2}
\]

The Newton formulas are therefore

\[
b_n=\sum_{k=0}^n c_k\binom nk,
\qquad
B_M(n)=\sum_{k=0}^{M-1}c_k\binom nk,
\]

so

\[
\boxed{0<B_M(n)<b_n.}
\tag{4.3}
\]

From (2.2),

\[
0<\epsilon_n\le\epsilon_0\rho^n.
\tag{4.4}
\]

On the other hand (3.4) gives a lower exponential base

\[
2\sqrt{1-\rho}>1
\]

for `|R_n|`. Thus

\[
\frac{\epsilon_n}{|R_n|}\longrightarrow0.
\tag{4.5}
\]

Now use the determinant sign exactly. Since

\[
a_j=\zeta(3)b_j-\epsilon_j,
\]

we have

\[
Z_{n,j}
=d_n(b_n\epsilon_j-\epsilon_nb_j).
\]

Hence

\[
\frac{T_n}{d_n}
=b_nA_M(n)-\epsilon_nB_M(n).
\tag{4.6}
\]

Because `A_M(n)=epsilon_n-R_n`, one may also write the cleaner identity

\[
\frac{T_n}{d_n}
=-b_nR_n+\epsilon_n\bigl(b_n-B_M(n)\bigr).
\tag{4.7}
\]

By (4.3)--(4.5), the second term in (4.7) is `o(b_n|R_n|)`. Therefore

\[
\boxed{
|T_n|
=d_nb_n|R_n|(1+o(1)).
}
\tag{4.8}
\]

In particular `T_n!=0` for all sufficiently large `n`.

So item (ii) of the requested audit passes after making the Newton-coefficient positivity explicit.

---

# 5. Lower height of `T_n`

The single Apéry-sum term `k=M=floor(n/2)` gives

\[
b_n
\ge
\binom nM^2\binom{n+M}{M}^2.
\]

Stirling's formula, in either parity, yields

\[
\boxed{
\log b_n\ge3n\log3-O(\log n).
}
\tag{5.1}
\]

Combining (3.4), (4.8), and (5.1),

\[
\log|T_n|
\ge
\log d_n
+\left(3\log3+\log2+\frac12\log(1-\rho)\right)n
-O(\log n).
\tag{5.2}
\]

Since

\[
\log d_n
=3\log\operatorname{lcm}(1,\ldots,n)
=3n+o(n),
\tag{5.3}
\]

we obtain

\[
\log|T_n|
\ge
\left(3+3\log3+\log2+\frac12\log(1-\rho)\right)n
+o(n).
\tag{5.4}
\]

---

# 6. Content bound `C_n <= 36 d_n^2/n^3`

From (1.1),

\[
C_n\mid\gcd(X_n,5X_n-6Y_n).
\]

Let `g=gcd(X_n,Y_n)`, write `X_n=gx`, `Y_n=gy`, and `gcd(x,y)=1`. Then

\[
\gcd(X_n,5X_n-6Y_n)
=g\gcd(x,5x-6y)
=g\gcd(x,6y)
\le6g.
\]

Thus

\[
\boxed{C_n\le6\gcd(X_n,Y_n).}
\tag{6.1}
\]

Set

\[
X_{n-1}'=d_na_{n-1},\qquad Y_{n-1}'=d_nb_{n-1},
\]

which are integers. Multiplying the given Casoratian

\[
a_nb_{n-1}-a_{n-1}b_n=\frac6{n^3}
\]

by `d_n^2` gives

\[
X_nY_{n-1}'-X_{n-1}'Y_n
=\frac{6d_n^2}{n^3}.
\]

Therefore

\[
\gcd(X_n,Y_n)\mid\frac{6d_n^2}{n^3},
\]

and hence

\[
\boxed{
C_n\le\frac{36d_n^2}{n^3}.
}
\tag{6.2}
\]

The right side is an integer multiple up to the harmless factor `36`, since `n|L_n`.

Using (5.3),

\[
\log C_n\le6n+o(n).
\tag{6.3}
\]

Subtracting this from (5.4) proves exactly

\[
\boxed{
\liminf_{n\to\infty}\frac{\log|Q_n|}{n}
\ge
3\log3+\log2+\frac12\log(1-\rho)-3.
}
\tag{6.4}
\]

The constant is positive. Indeed `rho=17-12sqrt(2)<3/4`, so `sqrt(1-rho)>1/2`; hence

\[
54\sqrt{1-\rho}>27>e^3,
\]

which is equivalent to positivity of the right side of (6.4).

Thus items (i)--(iii) of the requested audit are valid after the parity repair above.

---

# 7. Strict-top `p`-adic support: exact corrected theorem

Fix a prime

\[
p>\frac{n+1}{2},\qquad p\le n,\qquad p\ge7,
\]

and write

\[
n=p+r.
\]

The strict inequality implies

\[
0\le r<M<p.
\tag{7.1}
\]

Put

\[
D=\frac{d_n}{p^3},
\]

which is a `p`-adic unit.

## 7.1 Two standard block congruences

The Apéry Lucas congruence gives

\[
\boxed{b_{p+s}\equiv5b_s\pmod p}
\qquad(0\le s\le p-2).
\tag{7.2}
\]

For the companion, put `A_m=p^3a_m`. The Casoratian at `m=p` gives

\[
A_pb_{p-1}-(p^3a_{p-1})b_p=6.
\]

Now `a_{p-1}` is `p`-integral and `b_{p-1}\equiv1 (mod p)`, so

\[
A_p\equiv6\pmod p.
\]

Multiplying the Apéry recurrence by `p^3` and reducing successively modulo `p` shows that `A_{p+s}` obeys the same reduced recurrence as `6b_s`. Hence

\[
\boxed{p^3a_{p+s}\equiv6b_s\pmod p}
\qquad(0\le s\le p-2).
\tag{7.3}
\]

Thus

\[
X_n=D(p^3a_n),
\qquad
Y_n=D(p^3b_n).
\tag{7.4}
\]

## 7.2 Cardinal weights modulo `p`

The Lagrange-product form is

\[
\lambda_j(n)
=\prod_{\substack{0\le m<M\\m\ne j}}
\frac{n-m}{j-m}.
\tag{7.5}
\]

All denominators are `p`-units because `M<p`, and `n\equiv r (mod p)`. Therefore

\[
\boxed{
\lambda_r(n)\equiv1\pmod p,
\qquad
v_p(\lambda_j(n))=1\quad(j\ne r).
}
\tag{7.6}
\]

There is exactly one numerator factor `n-r=p` for `j!=r`.

## 7.3 Non-targets

If `p` does not divide `b_r`, then (7.3) makes `X_n` a `p`-unit, whereas `Y_n` is divisible by `p^3`. Hence `Z_{n,0}=X_n` is a unit, so

\[
v_p(C_n)=0.
\]

By (7.6),

\[
T_n\equiv Z_{n,r}\pmod p.
\]

But

\[
Z_{n,r}=X_nb_r-Y_na_r\equiv X_nb_r\not\equiv0\pmod p.
\]

Therefore

\[
\boxed{
p\nmid b_r\Longrightarrow p\nmid Q_n.}
\tag{7.7}
\]

This half of the claimed equivalence is unconditional.

## 7.4 Targets and the exact content valuation

Now assume

\[
p\mid b_r.
\tag{7.8}
\]

Since `p>=7`, neither `r=0` nor `r=1` can be a target; thus `r>=2`.

The Casoratian at `r` gives

\[
a_rb_{r-1}-a_{r-1}b_r=\frac6{r^3}.
\]

Because `r<p`, all quantities are `p`-integral. Reducing modulo `p` under (7.8) shows

\[
\boxed{a_r\in\mathbf Z_p^\times,\qquad b_{r-1}\in\mathbf Z_p^\times.}
\tag{7.9}
\]

Let

\[
x=v_p(X_n),\qquad y=v_p(Y_n).
\]

By (7.2)--(7.3), `x>=1` and `y>=4`.

For every `j<M`, both `a_j,b_j` are `p`-integral, so

\[
v_p(Z_{n,j})\ge\min(x,y).
\]

Conversely the first two determinants are

\[
Z_{n,0}=X_n,
\qquad
Z_{n,1}=5X_n-6Y_n.
\]

The change of coordinates

\[
(X_n,Y_n)\mapsto(Z_{n,0},Z_{n,1})
\]

has determinant `-6`, a `p`-unit. Hence

\[
\boxed{
v_p(C_n)=\min(x,y)=:c.
}
\tag{7.10}
\]

Using (7.6), after division by `p^c` all terms in `T_n` except the `j=r` term vanish modulo `p`. Thus

\[
\boxed{
p\mid Q_n\iff v_p(Z_{n,r})>c.}
\tag{7.11}
\]

But

\[
Z_{n,r}=X_nb_r-Y_na_r.
\]

By (7.8)--(7.9), the two summands have valuations at least `x+1` and exactly `y`, respectively. Therefore:

- if `x<y`, both summands are divisible by `p^{x+1}`, so `v_p(Z_{n,r})>c=x`;
- if `x>=y`, the second summand has exact valuation `y=c` while the first has valuation at least `x+1>c`, so `v_p(Z_{n,r})=c`.

Hence the exact target criterion is

\[
\boxed{
 p\mid Q_n
 \iff
 p\mid b_{n-p}
 \ \text{and}\ 
 v_p(d_na_n)<v_p(d_nb_n).
}
\tag{7.12}
\]

This is the key correction to the candidate theorem.

---

# 8. The missing orientation reduces to one depth-three seam

The condition in (7.12) can be sharpened further.

At a target, (7.9) and the block laws show

\[
p^3a_{n-1}\in\mathbf Z_p^\times,
\qquad
b_{n-1}\in\mathbf Z_p^\times.
\tag{8.1}
\]

Multiply the Casoratian at `n` by `p^3`:

\[
(p^3a_n)b_{n-1}-(p^3a_{n-1})b_n
=\frac{6p^3}{n^3}.
\tag{8.2}
\]

Since `n\equiv r not\equiv0 (mod p)`, every coefficient outside `p^3a_n` and `b_n` in (8.2) is a `p`-unit.

Put

\[
A=p^3a_n,
\qquad
B=b_n,
\qquad
\beta=v_p(B)\ge1.
\]

Equation (8.2) has the form

\[
uA-vB=p^3w,
\qquad u,v,w\in\mathbf Z_p^\times.
\tag{8.3}
\]

Therefore:

\[
\boxed{
\begin{array}{c|c}
\beta & v_p(A)\\ \hline
1,2 & \beta\\
>3 & 3\\
3 & \ge3\ \text{with no upper bound supplied by (8.3)}.
\end{array}}
\tag{8.4}
\]

Because `d_n/p^3` is a `p`-unit,

\[
v_p(d_na_n)=v_p(A),
\qquad
v_p(d_nb_n)=3+\beta.
\]

Thus the orientation in (7.12) is automatic in every case except

\[
\beta=3.
\]

In that case it fails precisely when `v_p(A)>=6`. Consequently

\[
\boxed{
 p\nmid Q_n
 \iff
 p\mid b_{n-p},\quad
 v_p(b_n)=3,\quad
 v_p(p^3a_n)\ge6.
}
\tag{8.5}
\]

This is an exact theorem, not a numerical observation.

The original support claim is therefore equivalent to the additional arithmetic assertion

\[
\boxed{
\text{For every strict-top target }p,
\quad
v_p(b_n)=3\Longrightarrow v_p(p^3a_n)\le5.
}
\tag{8.6}
\]

No part of the proposed height argument, the content estimate, the Casoratian alone, or the first-order Lucas/block congruences proves (8.6). It is a genuinely deeper `p`-adic lifting statement.

---

# 9. Corrected fully quantified theorem

Let `n>=4`, `M=floor(n/2)`, define `lambda_j`, `W_n`, `Z_{n,j}`, `C_n`, `T_n`, and `Q_n` exactly as in the question, and use the standard Apéry denominator theorem and Beukers moment representation.

Then:

1. **Integrality.** `Z_{n,j}`, `C_n`, `T_n`, and `Q_n` are integers (with `C_n>0`).

2. **Content.** For all `n>=4`,
   \[
   C_n\le6\gcd(d_na_n,d_nb_n)\le\frac{36d_n^2}{n^3}.
   \]

3. **Height.** With `rho=(sqrt2-1)^4`,
   \[
   \liminf_{n\to\infty}\frac{\log|Q_n|}{n}
   \ge
   3\log3+\log2+\frac12\log(1-\rho)-3>0.
   \]

4. **Strict-top non-target exclusion.** If
   \[
   p>\frac{n+1}{2},\quad p\le n,\quad p\ge7,
   \]
   and `p` does not divide `b_{n-p}`, then `p` does not divide `Q_n`.

5. **Exact strict-top target criterion.** Under the same conditions on `p`, if `p | b_{n-p}`, then
   \[
   p\mid Q_n
   \iff
   v_p(d_na_n)<v_p(d_nb_n).
   \]
   Equivalently,
   \[
   p\nmid Q_n
   \iff
   v_p(b_n)=3\ \text{and}\ v_p(p^3a_n)\ge6.
   \]

6. Therefore the claimed simpler equivalence
   \[
   p\mid Q_n\iff p\mid b_{n-p}
   \]
   follows **if and only if** the depth-three seam (8.6) is excluded for every strict-top target.

This is the strongest theorem obtained from the claimed argument without inserting an unproved `p`-adic transversality statement.

---

# 10. Audit checklist

- **(i) Alternating finite tail:** PASS after parity split. Even parity uses `mu_{M+1}<=mu_M`; odd parity requires the different identity `mu_M-mu_{M+1}=int t(1-t)^M dmu` and the positive constant `epsilon_1`. The same exponential rate results.

- **(ii) `epsilon_n B_M(n)` negligibility:** PASS. Positive Newton coefficients of `b_n` give `0<B_M(n)<b_n`, while `epsilon_n<=epsilon_0 rho^n` and the alternating remainder grows like `[2 sqrt(1-rho)]^n` up to powers of `n`.

- **(iii) content gcd and integrality:** PASS. `Z_0=X`, `Z_1=5X-6Y`, so `C<=6 gcd(X,Y)`; the Casoratian gives `gcd(X,Y)|6d_n^2/n^3`. All determinants are integral because `d_n a_j` is integral for `j<=n`.

- **(iv) support after division by `C_n`:** **FAIL as claimed.** The exact statement is (7.12)/(8.5). The missing all-target input is the exclusion of `v_p(b_n)=3` together with `v_p(p^3a_n)>=6`.

---

# Final conclusion

The proposed theorem should **not** be published with `p|Q_n iff p|b_{n-p}` as an unconditional clause. The determinant/content normalization does preserve a strict-top target except at one explicitly characterized deep `p`-adic seam, but the claimed proof does not rule that seam out.

The exponential-height lower bound is valid and, in fact, is the strongest part of the construction: even after dividing by the full common determinant content, the quotient remains exponentially large with the advertised rate.

A complete PASS would now require one new lemma only:

\[
\boxed{
 p>\frac{n+1}{2},\ p\ge7,\ p\mid b_{n-p},\ v_p(b_n)=3
 \Longrightarrow
 v_p(p^3a_n)\le5.
}
\]

Without that lemma, the corrected theorem above is the rigorous endpoint.