# Q8407 — five-flat high-prime obstruction: exact normalization and first toric residual

## Verdict

I do **not** have a theorem-level proof or a genuine counterexample to the actual 22-ray five-flat lemma from the connector-visible repository state. The four named local notes are not tracked on visible `main` (current head `3484a0081e8c9ea021b64aa60a1eff6bd0b8eeb4`), and the exact 22-ray Laurent polynomial `Lambda` and the literal coefficient-difference numerator defining `rho_{h,s}` are therefore unavailable to this audit. I will not invent them.

What can be closed rigorously from the normalization in Q8407 is substantial:

1. all half-integer quantities live in `Z[1/2]`; **no odd denominator prime exists**;
2. the five `Y`-values and the five coefficient gaps of `F_h` are related by an explicit unitriangular matrix over `Z`, so their simultaneous vanishing is exactly equivalent over every odd residue field;
3. for `P>2h+2`, shifting to `z=x-(lambda_h-2)` makes the five-flat event exactly
   `q(z) | Y_h(lambda_h-2+z)` in `F_P[z]`, with
   `q(z)=(z+2)(z+1)z(z-1)(z-2)=z^5-5z^3+4z`;
4. the normalization alone cannot rule out five-flatness: for every `h>=3` the abstract `rho`-space has an explicit nonzero **integral** sublattice of rank at least `2h-4` on which all five quantities vanish. Thus the missing input really must use the special 22-ray `Lambda`/coefficient-difference structure, not Newton-basis algebra alone;
5. the first algebraic identity that would close the actual lemma is a single **twisted toric Euler / Griffiths–Dwork certificate**, stated precisely in Section 8. It yields an identity

   `sum_{k=0}^4 A_{h,k} C_{h,k} = h`

   with dyadic coefficients (or, for the high-prime version, coefficients localized only at primes `<=2h+2`), where `C_{h,k}=[u^(h-k)]F_h`. Any common odd prime of the five entries then divides `h`; in the localized version no prime `P>2h+2` can annihilate all five.

This is the earliest residual I can isolate without silently assuming the missing Laurent polynomial. A full five-dimensional recurrence in `h` would be stronger than necessary.

There is also one logical correction to the question's wording. The statement

> no odd prime `P>2h+2` annihilates all five values

is **not by itself equivalent** to

> every odd common prime divisor of the five values divides `h`.

The second statement also rules out primes `P<=2h+2` with `P∤h`. It implies the high-prime statement, but the converse needs an additional small-prime lemma. The unavailable `tmux11_fiveflat_r_minus_one_addendum.md` may contain such a lemma; I cannot certify it from the visible tree. Everything below separates the high-prime theorem from that stronger global-gcd formulation.

---

## 1. Exact ring of definition

Fix `h>=1` and write

`lambda_h = floor(h/2) - 1/2`,

`alpha_h = lambda_h - 4`.

Both are half-integers of the form `m-1/2` with `m in Z`.

### Lemma 1.1 — half-integral binomial coefficients are dyadic

For every integer `m` and every `s>=0`,

`binom(m-1/2,s) in Z[1/2]`.

#### Proof

In `Z[1/2][[T]]`,

`(1+T)^(m-1/2) = (1+T)^m (1+T)^(-1/2)`.

For every integer `m`, the coefficients of `(1+T)^m` are integers (also for negative `m`, using `binom(-a,n)=(-1)^n binom(a+n-1,n)`). Moreover

`(1+T)^(-1/2) = sum_{n>=0} (-1)^n binom(2n,n) T^n / 4^n`,

whose coefficients are dyadic. Their Cauchy product is therefore in `Z[1/2][[T]]`. The coefficient of `T^s` is `binom(m-1/2,s)`. ∎

Consequences:

- every value `Y_h(lambda_h-j)` is in `Z[1/2]` whenever every `rho_{h,s}` is integral;
- every coefficient `[u^n](1-u)^alpha_h R_h(u)` is in `Z[1/2]`;
- reduction modulo any odd prime is intrinsically defined by inverting `2` only.

Thus there is no hidden prime from a half-integer denominator.

For an odd prime `P`, I use `P | a` for `a in Z[1/2]` to mean that the numerator of `a` in lowest dyadic form is divisible by `P`.

---

## 2. Parity and the `r=h+1` indexing

Write either `h=2m` or `h=2m+1`. In both cases

`lambda_h = m - 1/2`,

so the exponent in `F_h` is

`alpha_h = m - 9/2`.

In particular

`lambda_{h+2}=lambda_h+1`,

`alpha_{h+2}=alpha_h+1`.

This explains why any exact coefficient-gap recurrence obtained from the toric model should naturally be organized in parity-preserving steps `h -> h+2`, rather than guessed one step at a time.

If the local notes use `r=h+1`, then:

- `h=2m`, `r=2m+1` odd:
  `lambda_h=(r-2)/2`, `alpha_h=(r-10)/2`;
- `h=2m+1`, `r=2m+2` even:
  `lambda_h=(r-3)/2`, `alpha_h=(r-11)/2`.

Most importantly,

`h = r-1`,

and the high-prime hypothesis becomes

`P > 2h+2 = 2r`.

Hence a prime in the target range is automatically a unit on `h=r-1` and on every nonzero integer whose absolute value is at most `2h+2`.

---

## 3. Coefficient-extraction identity

By definition

`R_h(u)=sum_{t=-h}^h rho_{h,h-t} u^t`.

Set `s=h-t`; then `s=0,...,2h` and `t=h-s`, so

`R_h(u)=sum_{s=0}^{2h} rho_{h,s} u^(h-s)`.

For arbitrary `x`, expand

`(1-u)^x = sum_{s>=0} (-1)^s binom(x,s)u^s`.

The coefficient of `u^h` in `(1-u)^x R_h(u)` is therefore

`sum_{s=0}^{2h} (-1)^s rho_{h,s} binom(x,s)=Y_h(x)`.

Thus the displayed identity in Q8407 is exact:

`Y_h(x)=[u^h](1-u)^x R_h(u)`.

No convergence issue is present: for a fixed coefficient only finitely many terms contribute.

---

## 4. Exact triangular equivalence

Put

`F_h(u)=(1-u)^(lambda_h-4)R_h(u)`

and define

`C_{h,k}:=[u^(h-k)]F_h(u)`, `0<=k<=4`.

Since

`R_h=(1-u)^(4-lambda_h)F_h`,

for `j=0,1,2,3,4` we get

`Y_h(lambda_h-j)=[u^h](1-u)^(4-j)F_h(u)`.

Writing `m=4-j`,

`Y_h(lambda_h-j)=sum_{k=0}^m (-1)^k binom(m,k) C_{h,k}`.     (4.1)

Explicitly, if `A_j=Y_h(lambda_h-j)`,

- `A_4=C_0`;
- `A_3=C_0-C_1`;
- `A_2=C_0-2C_1+C_2`;
- `A_1=C_0-3C_1+3C_2-C_3`;
- `A_0=C_0-4C_1+6C_2-4C_3+C_4`.

The matrix is unitriangular. Binomial inversion gives

`C_k=sum_{m=0}^k (-1)^m binom(k,m) A_{4-m}`.     (4.2)

Therefore over `Z[1/2]`, over `Q`, and after reduction modulo **every odd prime**,

`A_0=...=A_4=0  <=>  C_{h,0}=...=C_{h,4}=0`.     (4.3)

This part requires no bound on `P`.

There is also a direct formula useful for exact checking. Since the exponent of `rho_{h,s}` in `R_h` is `h-s`,

`C_{h,k} = sum_{s=k}^{2h} rho_{h,s} (-1)^(s-k) binom(alpha_h,s-k)`.     (4.4)

The coefficient of `rho_{h,k}` is exactly `1`.

---

## 5. The centered five-root quotient

Let

`a_h=lambda_h-2`,

and define the shifted polynomial

`Ytilde_h(z)=Y_h(a_h+z)`.

The five points `x=lambda_h-j`, `j=0,...,4`, become

`z=2,1,0,-1,-2`.

Define

`q(z)=(z+2)(z+1)z(z-1)(z-2)`

`    =z(z^2-1)(z^2-4)`

`    =z^5-5z^3+4z`.

For `h>=1` and `P>2h+2`, we have `P>=5`, and the five roots are distinct in `F_P`. Also `s<=2h<P`, so the polynomial `binom(x,s)=x(x-1)...(x-s+1)/s!` itself has a valid reduction in `F_P[x]`. Hence

`Y_h(lambda_h-j)=0 (mod P) for all j=0,...,4`

if and only if

`q(z) | Ytilde_h(z)` in `F_P[z]`.     (5.1)

For `P>5`, one may also write

`q(z)=120 binom(z+2,5)`,

but the factored form is the safer normalization because it still makes sense at `P=5`.

This quotient language is useful for Griffiths–Dwork reduction: it says the obstruction is the vanishing of a class in the fixed rank-five algebra

`F_P[z]/(q)`.

It does **not** by itself prove nonvanishing.

---

## 6. A rigorous generic obstruction: the special `Lambda` is indispensable

The normalization alone has far too many degrees of freedom.

Consider the `Z[1/2]`-linear map

`T_h : Z[1/2]^(2h+1) -> Z[1/2]^5`,

`rho=(rho_0,...,rho_{2h}) |-> (C_{h,0},...,C_{h,4})`,

with `C_{h,k}` given by (4.4).

Look at the five columns `rho_0,...,rho_4`. The coefficient of `rho_s` in row `k` is zero for `s<k`, and the diagonal coefficient at `s=k` is `1`. Hence this `5x5` minor is triangular with determinant `1`.

### Proposition 6.1

For every `h>=2`, `T_h` has rank `5`. For every `h>=3`, its kernel has rank

`(2h+1)-5 = 2h-4`.

Moreover the kernel contains nonzero **integral** vectors.

#### Proof

Rank `5` follows from the unit minor. The rank count gives the kernel dimension over `Q`. All matrix coefficients are dyadic by Lemma 1.1, so a nonzero dyadic kernel vector may be multiplied by a power of `2` to become integral. ∎

A particularly explicit construction is: for `h>=3`, set `rho_5=1` and `rho_s=0` for `s>5`; then solve the equations `C_4=0,C_3=0,...,C_0=0` successively for `rho_4,rho_3,...,rho_0`. Each pivot is `1`. Scaling by a power of `2` produces a nonzero integral vector with all five gaps and all five `Y`-values exactly zero.

This is **not** a counterexample to Q8407, because the real `rho_{h,s}` are constrained constant terms coming from the fixed 22-ray `Lambda`. It proves only that no argument which forgets that toric origin can possibly close the lemma.

This is why I do not pursue generic polynomial-root counting, affine Smith form, or fixed-width Turán algebra here.

---

## 7. Exact toric Euler identity available for any Laurent polynomial

Let the exact 22-ray Laurent polynomial, once supplied, lie in

`A = Z[1/2][x_1^±1,...,x_d^±1]`

(or include the staircase variable among the torus variables, as in the repository's literal normalization).

For a torus variable `x_i`, write

`Theta_i = x_i d/dx_i`.

For every Laurent polynomial `G`,

`CT(Theta_i G)=0`.     (7.1)

Apply this to `G Lambda^h`:

`0 = CT Theta_i(G Lambda^h)`

`  = CT((Theta_i G)Lambda^h + h G(Theta_i Lambda)Lambda^(h-1))`.

Equivalently, define the **twisted Euler divergence**

`E_{h,i}(G) := Lambda Theta_i G + h G Theta_i Lambda`.     (7.2)

Then

`CT(Lambda^(h-1) E_{h,i}(G))=0`.     (7.3)

This identity is exact over `Z`, uses no smoothness assumption, and is the algebraic core of toric Griffiths–Dwork integration by parts.

The crucial feature is the literal factor `h` produced by differentiating `Lambda^h`. Since the empirical common-prime statement singles out `h=r-1`, this is the first mechanism in the proposed routes that naturally manufactures exactly the desired factor rather than merely bounding a determinant.

---

## 8. The first precise residual identity that closes the lemma

Because the literal coefficient-difference numerator is unavailable, I state this in a way that can be instantiated mechanically once that numerator is inserted.

By linearity of the constant term, the actual coefficient gaps can be written

`C_{h,k}=CT(Lambda^h Psi_{h,k})`, `0<=k<=4`,     (8.1)

for five explicit finite Laurent expressions `Psi_{h,k}` obtained from the repository's exact coefficient-difference definition and formula (4.4). Equivalently put

`Omega_{h,k}=Lambda Psi_{h,k}`, so

`C_{h,k}=CT(Lambda^(h-1) Omega_{h,k})`.     (8.2)

The exact formulas for `Psi/Omega` are where the missing 22-ray local note is needed.

### Residual [TGD-5] — twisted Griffiths–Dwork five-flat certificate

It is enough to prove the following **finite Laurent-polynomial identity**.

For each `h` (separately for the two parity classes if necessary), exhibit:

- coefficients `A_{h,0},...,A_{h,4}`;
- Laurent polynomials `G_{h,1},...,G_{h,d}`;
- one Laurent polynomial `Omega_h^*`;

such that

`sum_{k=0}^4 A_{h,k} Omega_{h,k} - h Omega_h^*`

`    = sum_{i=1}^d E_{h,i}(G_{h,i})`,     (TGD-5)

and

`CT(Lambda^(h-1) Omega_h^*) = 1`.     (8.3)

If the desired global odd-gcd statement is the target, require

`A_{h,k} in Z[1/2]`.     (8.4-global)

If only the stated high-prime lemma is targeted, it is enough that every denominator of every `A_{h,k}` has prime divisors at most `2h+2`; equivalently

`A_{h,k} in Z[1/(2(2h+2)!)]`.     (8.4-high)

### Why [TGD-5] closes the theorem

Take `CT(Lambda^(h-1) * -)` of (TGD-5). Every twisted Euler term vanishes by (7.3), and (8.2)–(8.3) give

`sum_{k=0}^4 A_{h,k} C_{h,k} = h`.     (8.5)

Now suppose an odd prime `P>2h+2` annihilates all five `C_{h,k}`. All denominators allowed in (8.4-high) are `P`-units, so reducing (8.5) modulo `P` gives

`0 = h (mod P)`.

But `0<h<P`, contradiction.

If the stronger dyadic condition (8.4-global) holds, then any odd prime dividing all five `C_{h,k}` divides the right side `h`, proving the global common-prime statement directly.

Thus [TGD-5] is a complete algebraic bridge from the exact 22-ray torus to the desired arithmetic statement.

### Why this residual is genuinely earlier than a full recurrence

One could instead seek a parity-preserving `5x5` transfer law for

`C_h=(C_{h,0},...,C_{h,4})^t`,

using `lambda_{h+2}=lambda_h+1`, and prove its determinant is a unit for `P>2h+2`. That would allow descent `h -> h-2`. But this asks Griffiths–Dwork reduction to close an entire five-dimensional state under multiplication by `Lambda^2`.

[TGD-5] needs only **one cohomological relation**: a single linear functional on that state whose value is `h`. It is therefore the first residual identity I would try to prove from the 22-ray data.

---

## 9. How to attack [TGD-5] from the exact 22 rays

This is a proof program, not an empirical recurrence search.

### 9.1 Work in the twisted Euler quotient

For fixed `h`, quotient the finite Laurent support needed for the five `Omega_{h,k}` by the submodule generated by the exact twisted divergences

`E_{h,i}(G)=Lambda Theta_iG+hGTheta_iLambda`.

In this quotient [TGD-5] becomes the class identity

`sum_k A_{h,k}[Omega_{h,k}] = h[Omega_h^*]`.     (9.1)

All operations are exact linear algebra over `Z[1/2]` or the high-prime localization. No guessed recurrence is involved.

### 9.2 Use the Newton polytope before solving coefficients

The 22 rays should first be grouped by their exponent vectors and the five obstruction numerators placed in the Newton polytope. Toric Euler relations preserve sharply controlled support. If the right side `Omega_h^*` can be chosen on an exposed face/corner for which the constant term in (8.3) has a unique path, then (8.3) becomes an exact monomial calculation, not a positivity or height argument.

This is the most promising place for the `r=-1` shell mentioned in the local addendum: a singular boundary face often leaves one uncancelled toric residue while interior contributions are Euler-exact.

### 9.3 Singular-residue test

Before attempting full reduction, apply every toric boundary residue functional to (9.1). Twisted Euler divergences have zero residue after the corresponding integration-by-parts normalization. If four of the five obstruction classes die under a boundary residue and the fifth gives a multiple of `h`, that residue itself supplies the coefficients `A_{h,k}`.

Conversely, if a boundary residue gives an obstruction inconsistent with `h[Omega_h^*]`, [TGD-5] in that form is impossible and the residue identifies the exact extra class that a full Griffiths–Dwork reduction must retain.

### 9.4 Parity should be separated at the start

For both `h=2m` and `h=2m+1`, `alpha_h=m-9/2`. Therefore the binomial-convolution part of `Omega_{h,k}` is identical within the pair; only the `Lambda^h` power and the original coefficient-difference term change.

A clean proof should reduce the even and odd families separately, then compare their certificates under multiplication by one copy of `Lambda`. Mixing parity before this step needlessly obscures the half-integer shift.

### 9.5 Denominator audit is part of the theorem

A computer algebra reduction over `Q` is not enough. Every division in Griffiths–Dwork reduction must be recorded. To prove the high-prime theorem it suffices to show the denominator support lies in primes `<=2h+2`. To prove the stronger “common odd primes divide `h`” statement, all odd denominators must cancel, leaving only powers of `2`.

This denominator support is exactly what the verifier for the eventual certificate should check.

---

## 10. The centered quotient `q(z)` and what it does — and does not — buy

The identity

`q(z)=z(z^2-1)(z^2-4)`

is valuable because it packages the five evaluations as one rank-five remainder problem. It also has the Newton form

`q(z)=120 binom(z+2,5)`

away from `2,3,5`.

However, `q | Ytilde_h` is not by itself a contradiction: `deg Y_h` can be as large as `2h`, and Proposition 6.1 shows that abstract coefficient data can satisfy the divisibility with a large kernel.

The useful role of `q` in a toric proof is instead to reduce the obstruction module to five coordinates before applying Euler/GD relations. In a basis `1,z,z^2,z^3,z^4` of the quotient by `q`, [TGD-5] asks for one exact functional of the five remainder coordinates to equal `h`.

---

## 11. Why ordinary height, LLL, or guessed recurrences do not settle this

- A finite scan establishes only finite evidence; Q8407 explicitly asks for all `h` and all large `P`.
- LLL may suggest coefficients `A_{h,k}` for particular rows but does not prove a uniform toric identity or denominator control.
- Generic recurrence guessing in `h` is especially risky because the natural normalization changes by parity and `alpha_{h+2}=alpha_h+1`.
- The known fixed-width Turán and affine Smith obstructions concern what can be extracted from finite windows or ambient linear algebra. Proposition 6.1 already shows the abstract five-flat map has a large kernel; the new content must be the **twisted Euler class of the actual 22-ray constant term**.

Accordingly I do not recycle those routes.

---

## 12. Status classification

### THEOREM

- Half-integer values and `F_h` coefficients are dyadic; no odd denominator prime occurs.
- The five `Y`-values and five `F_h` coefficient gaps are exactly unitriangularly equivalent.
- For `P>2h+2`, the five-flat event is exactly `q | Ytilde_h` in `F_P[z]`.
- The abstract five-flat map has rank `5` and a kernel of rank `2h-4` for `h>=3`; nonzero integral synthetic five-flat vectors exist.
- The twisted Euler integration-by-parts identity (7.3) is exact for every Laurent polynomial.
- [TGD-5], if established with the stated denominator control, proves the desired high-prime lemma; with dyadic coefficients it proves the stronger global odd-common-prime statement.

### NOT PROVED HERE

- The actual 22-ray [TGD-5] certificate, because the literal `Lambda` and coefficient-difference numerator are not connector-visible.
- Any actual counterexample.
- The small-prime implication needed to make the high-prime statement logically equivalent to the global “common odd prime divides `h`” statement.

### FINITE EVIDENCE

None is used as a theorem. I deliberately did not run or extend the finite scan.

---

## 13. Exact verifier

The companion file

`problems/3.2/ORACLE_COMM/chatgpt_q8407_fiveflat_residual_verify.py`

uses only the Python standard library and exact `fractions.Fraction` arithmetic. It verifies, for arbitrary supplied exact `rho` data:

- the `R_h`/`Y_h` coefficient-extraction normalization;
- dyadic denominators for all half-integer binomial coefficients used;
- the five-by-five triangular transform and inverse;
- the centered `q(z)` factorization and polynomial remainder criterion;
- reduction modulo a supplied prime `P>2h+2`;
- the parity and `r=h+1` formulas;
- the explicit synthetic-kernel construction showing that the normalization alone has nonzero integral five-flat vectors.

It does **not** manufacture the missing 22-ray `rho` and therefore does not pretend to prove the actual arithmetic target.

Expected local use once exact `rho` is exported as JSON:

```bash
python3 problems/3.2/ORACLE_COMM/chatgpt_q8407_fiveflat_residual_verify.py --self-test
python3 problems/3.2/ORACLE_COMM/chatgpt_q8407_fiveflat_residual_verify.py \
  --rho-json /path/to/exact_rho.json --prime P --assert-target
```

JSON schema:

```json
{
  "h": 7,
  "rho": ["1", "-3/2", "17", "0", "... exactly 2h+1 entries ..."]
}
```

The verifier was written but **not executed in this response**, in accordance with Q8407's instruction not to use Python/code-interpreter/sandbox.

---

## 14. Exact next proof obligation

Expose the literal 22-ray `Lambda` and the literal coefficient-difference numerator in a tracked/readable form, form the five exact Laurent numerators `Omega_{h,k}` via (4.4), and prove [TGD-5].

The most economical theorem to target is therefore:

> **[FIVEFLAT-TORIC-EULER].** For each `h>=1`, the five obstruction classes in the twisted toric Euler quotient admit coefficients with denominator support at most `2h+2` whose linear combination is `h` times a unit constant-term class. If the coefficients are dyadic, every odd common prime divisor of the five obstruction values divides `h`.

That is the first precise algebraic residual identity whose proof closes the requested lemma without a finite scan, LLL inference, generic recurrence guess, or a return to the known fixed-width/Smith dead ends.
