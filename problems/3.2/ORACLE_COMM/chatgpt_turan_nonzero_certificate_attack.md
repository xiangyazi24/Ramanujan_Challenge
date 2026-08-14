# Turan nonzero-certificate attack for the top-half Apéry radical

## Verdict

**Verdict: smallest exact surviving lemma, not a proof and not a disproof of the requested `O(sqrt(n log n))` radical bound.**

The zero-output branch cannot be eliminated by proving that short kernel relations do not exist.  In fact, for the actual integer Turan row `(E_1,...,E_r)`, short exact kernel relations are **unavoidable at exactly the `exp(O(n/r))` coefficient scale**: a one-line box-pigeonhole/Siegel argument gives a nonzero `c` with

```text
sum_j c_j E_j = 0,
max |c_j| <= exp(O(n/r))
```

whenever `log max|E_j|=O(n)`.  The same remains true after adding any fixed number of extra linear first-jet rows.  Thus the hoped-for route “prove a lower bound on short kernel relations, so Minkowski cannot land in the kernel” is theorem-level impossible.

The actual Apéry/Newton structure does, however, give two useful exact refinements.

1.  With `D_j=F_j-F_{j+1}` and the existing Pascal normalizer `q_j`, define

    ```text
    R_j = (D_{j-1}-D_j)/q_j.
    ```

    Then `R_j` is an integer and

    ```text
    E_j = F_j R_j - D_{j-1} D_j / q_j.                (A)
    ```

    At every common top-half candidate prime `p`, target or non-target,

    ```text
    E_j == b_{n-p} R_j  (mod p).                       (B)
    ```

    This is the exact first-jet content behind the rank-one alias in the current terminal Turan audit.

2.  The most obvious Wronskian that tries to use `(B)` to remove the kernel,

    ```text
    W_j = R_{j+1} E_j - R_j E_{j+1},
    ```

    is **universally candidate-divisible**, not target-selective.  Exact algebra gives

    ```text
    W_j = D_j (D_{j-1} D_{j+1} - D_j^2) / (q_j q_{j+1}).   (C)
    ```

    For the actual Apéry row `n=200,r=6`, the common candidate primorial divides `gcd(W_1,...,W_5)` to exactly the first power at every candidate, including the genuine targets `139,181`.  After removing the universal candidate content, neither target remains.  Thus the first-jet Wronskian does **not** contain a hidden second target digit.

The smallest target-selective object I found that survives this obstruction is instead the **value--Turan second determinantal divisor**.  Put

```text
J_{i,j} = E_i F_j - E_j F_i,
Delta_FE(n,r) = gcd_{1 <= i < j <= r} |J_{i,j}|.
```

For every genuine target prime in the common top-half interval,

```text
p^2 | J_{i,j}  for every i<j,
```

because both `F_i` and `E_i` vanish mod `p`.  In contrast, at a non-target common candidate,

```text
J_{i,j} == b_{n-p}^2 (R_i-R_j)  (mod p),              (D)
```

so there is **no universal candidate-primorial factor**.  This is a genuine Smith/determinantal refinement, not the contaminated `R`-Wronskian.

Exact independent computations give, for six columns,

| `n` | common targets | exact `Delta_FE(n,6)` |
|---:|---|---|
| 200 | `139,181` | `2 * 5^2 * 139^2 * 181^2` |
| 272 | `191,233` | `2 * 191^2 * 233^2` |
| 300 | `191,227` | `37 * 191^2 * 227^2` |
| 321 | `179,193,211` | `179^2 * 193^2 * 211^2` |

These are exact finite identities only; **I do not infer an asymptotic bound from them.**  They identify the sharp surviving theorem:

> **[FE-SMITH]** For some `r=r(n)->infinity`, `r=o(n)`, prove that the `2 x r` matrix with rows `(F_1,...,F_r)` and `(E_1,...,E_r)` has rank two and
> 
> ```text
> 0 < Delta_FE(n,r) <= exp(O(n/r + r log n)).
> ```

Then, without any Bezout coefficient control,

```text
rad_T(n)^2 | Delta_FE(n,r)
```

for the common target set `T`, and therefore

```text
log rad_T(n) = O(n/r + r log n).
```

Choosing `r ~ sqrt(n/log n)` gives `O(sqrt(n log n))` on this top strip.

This **does not prove** the requested small nonzero linear combination `H=sum c_jE_j`; it is a different, Smith-determinantal route.  That distinction matters.  In the same actual row `n=200`,

```text
gcd(E_1,E_2) = 503180 = 2^2 * 5 * 139 * 181,
```

while the coordinatewise-minimal two-term Bezout representation of this gcd has coefficient bit lengths `2092` and `2091`; both `E_1,E_2` themselves have `2112` bits.  A small gcd is therefore emphatically **not** a small-coefficient `H` certificate.

The exact first-jet identities, the actual `n=200` hostile checks, and the four `Delta_FE` records are reproduced by the owned, dependency-free verifier

```text
problems/3.2/ORACLE_COMM/chatgpt_turan_nonzero_certificate_verify.py
```

on `chatgpt-drop`.

---

## 1. Source audit and limitation

The GitHub connector visible in this session resolves `main` at

```text
734a5a84c1e4fd8703a811aadaa2b4c7f532b20e
```

and the five newly named Turan/prefix sources in the question are not present in that visible tree:

```text
problems/3.2/ORACLE_COMM/codex_many_prefix_cancellation.md
problems/3.2/ORACLE_COMM/tmux11_first_jet_turan_audit.md
problems/3.2/research/scripts/q32_many_prefix_cancellation_verify.py
problems/3.2/research/scripts/q32_first_jet_turan_scan.py
problems/3.2/research/scripts/q32_first_jet_bezout_search.py
```

Direct fetches of the named paths returned `404`, repository code search returned no copies, and the visible `chatgpt-drop` tree did not contain them either.  I therefore **did not invent their contents** and do not claim to have re-audited their particular LLL basis/scaling.  I take the definitions and the finite-LLL premise stated in the question as the current interface.

The connector-visible direct ancestors that I did read are:

```text
problems/3.2/research/scripts/q32_terminal_turan_hankel_audit.py
    blob c54d01fd0554201a703a7bf21796682e2be9ff6b

problems/3.2/research/scripts/q32_terminal_family_audit.py
    blob 96d53d33583ac4651f2b178fbe8542f288e690c0

problems/3.2/research/scripts/q32_terminal_bernstein_audit.py
    blob 5dd4f6f22da9fbc0a2fe8084bb5af0442936de90

problems/3.2/research/scripts/q32_cartier_packet_audit.py
    blob e1f802e8f085ae588df52ba7c949974bf4168ca8

problems/3.2/proof.tex
    blob efbede7ea8ac6e040a5d380860ef5009a564fb01
```

The first file contains exactly the `F_j`, `E_j`, Pascal normalization, first rank-one divided identity, adjacent elimination identity, and Hankel checks used below.  The second proves the terminal Newton increment identity.  The fourth supplies the physical scalar shell.  The current `proof.tex` contains the top-half radical reduction.

No shared TeX, `DOCTRINE`, or `RUN_LOG` was edited.

---

## 2. The top-half target and the endpoint

The visible `proof.tex` proves for every prime `p>=7` with

```text
n/2 < p <= n,
```

putting `h=n-p`, that

```text
v_p(G_n)=0  iff  p does not divide b_h.
```

Its Lucas corollary gives, for `p!=5`,

```text
v_p(G_n)>=1  iff  p | b_{n-p}  iff  p | b_n.
```

It then records the exact reduction

```text
log rad_{(n/2,n]}(b_n)=o(n)
```

as the remaining all-`n` top-half radical problem.

For the present Turan construction I enforce the requested strict endpoint exclusion

```text
p < n.
```

If `p=n` is prime then `h=n-p=0` and

```text
b_h=b_0=1,
```

so `p=n` is never a physical target anyway.  The verifier constructs candidate lists with `primes_up_to(n-1)` and separately checks the `b_0=1` endpoint convention.

---

## 3. Exact terminal indexing and Pascal normalizers

Put

```text
M = n-1,
s = floor(M/2)+1,
L0 = M-s.
```

For the scalar shell `C_M(d)`, the existing terminal family uses

```text
f_k = G_{M-k,k}(C_M),
F_j = f_{L0-j}.
```

Thus `F_0,F_1,...` are successive terminal-prefix Newton values with all stencils sharing the upper endpoint `M`.

Set

```text
P_j = binomial(n,L0-j).
```

For `j>=1`, the once-divided Turan carrier is exactly

```text
q_j = gcd(P_{j-1},P_j),
E_j = (F_{j-1}F_{j+1}-F_j^2)/q_j.                    (3.1)
```

The current audit checks that `(3.1)` is integral, not merely rational.

For later use, put

```text
k = L0-j,
g_j = gcd(n-k,k+1).
```

The adjacent-binomial ratio gives the exact primitive decomposition

```text
A_j := P_{j-1}/q_j = (n-k)/g_j,
B_j := P_j/q_j     = (k+1)/g_j,
gcd(A_j,B_j)=1.                                          (3.2)
```

Indeed

```text
P_{j-1}/P_j = (n-k)/(k+1),
```

and cancelling `g_j` leaves a coprime numerator and denominator, so the denominator divides `P_j` and the displayed gcd follows.

The coefficients `A_j,B_j` are at most `n`.  This small primitive part is why it is important not to replace `q_j` by an unnormalized Pascal coefficient.

---

## 4. The exact common-prime interval

To use `E_1,...,E_r` we need `F_0,...,F_{r+1}` to share the physical node `p-1`.  The exact safe condition inherited from the terminal audit is

```text
p > s+r+1.                                               (4.1)
```

Accordingly define

```text
C(n,r) = { prime p : s+r+1 < p < n }.
```

For `r=o(n)`, this is the claimed

```text
n/2 + O(r) < p < n
```

strip, with the exact `+1` bookkeeping fixed.

For `p in C(n,r)`, write `h=n-p`.  Since `p>n/2` and `p<n`,

```text
1 <= h <= p-2.
```

The shell Cartier/Newton selector gives, simultaneously for `0<=j<=r+1`,

```text
F_j == b_h  (mod p).                                     (4.2)
```

The exact audit verifies the same statement by checking that all consecutive increments are divisible by `p`.

Moreover, for every `1<=j<=r`,

```text
v_p(P_{j-1}) = v_p(P_j) = v_p(q_j) = 1.                 (4.3)
```

This is also elementary.  Write `n=p+h<2p`.  Under `(4.1)` the relevant lower binomial indices exceed `h` but are less than `p`.  Hence `n!` contributes one factor of `p`, while the two denominator factorials contribute none.  The current audit checks `(4.3)` candidate by candidate.

Define the actual target set

```text
T(n,r) = { p in C(n,r) : p | b_{n-p} }.
```

Then `(4.2)` says exactly that every `p in T(n,r)` divides every `F_0,...,F_{r+1}`.

---

## 5. Terminal increments and the integer first jet

Let

```text
D_j = F_j-F_{j+1}.                                       (5.1)
```

The exact terminal-family identity in `q32_terminal_family_audit.py` says, in the original `f_k` orientation,

```text
f_k-f_{k-1} = (-1)^k binomial(n,k) B_k,
```

where `B_k` is an explicit integer terminal boundary packet.  In `F_j` orientation this becomes

```text
D_j = P_j X_j,                                           (5.2)
```

with

```text
X_j = (-1)^(L0-j) B_{L0-j} in Z.
```

Thus both adjacent increments entering `E_j` contain the exact Pascal factors appearing in `q_j`.

Now define

```text
R_j = (D_{j-1}-D_j)/q_j.                                 (5.3)
```

This is an **integer**.  Substituting `(5.2)` and `(3.2)` gives the primitive form

```text
R_j = A_j X_{j-1} - B_j X_j.                             (5.4)
```

No division by a target prime has occurred in `(5.3)` or `(5.4)`; all quantities are characteristic-zero integers.

---

## 6. First-jet decomposition of the Turan carrier

### Theorem 6.1 -- exact first-jet identity

For every admissible `j`,

```text
E_j = F_j R_j - D_{j-1}D_j/q_j.                          (6.1)
```

### Proof

From `(5.1)`,

```text
F_{j-1}=F_j+D_{j-1},
F_{j+1}=F_j-D_j.
```

Therefore

```text
F_{j-1}F_{j+1}-F_j^2
 = F_j(D_{j-1}-D_j)-D_{j-1}D_j.
```

Divide by `q_j`.  The first quotient is `R_j` by `(5.3)`.  The second is integral because `q_j` divides each of `D_{j-1}` and `D_j` by `(5.2)` and the definition of `q_j`.  This proves `(6.1)`.

### Corollary 6.2 -- exact candidate alias

For every `p in C(n,r)` and `1<=j<=r`,

```text
E_j == b_{n-p} R_j  (mod p).                              (6.2)
```

### Proof

By `(4.3)`, `p` divides both `D_{j-1}` and `D_j`, while `q_j` contains exactly one factor of `p`.  Hence

```text
D_{j-1}D_j/q_j == 0  (mod p).
```

Reduce `(6.1)` and use `F_j==b_{n-p}` from `(4.2)`.

For a target, `b_{n-p}=0`, so `(6.2)` recovers

```text
p | E_j.
```

For a non-target, `b_{n-p}` is a unit and `(6.2)` identifies the normalized first jet seen by `E_j`.

This is the characteristic-zero form of the rank-one divided identity that the visible Turan audit checks numerically for both target and non-target candidates.

---

## 7. Why short kernel relations cannot be excluded

The question asks whether the zero-output Minkowski branch can be killed by a lower bound on short relations in

```text
c -> sum_j c_j E_j.
```

It cannot.  The obstruction is elementary and unconditional.

### Proposition 7.1 -- box-Siegel zero relation

Let `A` be any `k x r` integer matrix with `1<=k<r` and

```text
max_{i,j}|A_{i,j}| <= B,
B>=1.
```

There exists a nonzero integer vector `c` with `Ac=0` and

```text
||c||_infinity <= Q,
Q <= 2 (3 r B)^(k/(r-k)) + 1.                            (7.1)
```

In particular,

```text
log ||c||_infinity
 <= (k/(r-k))(log B + log(3r)) + O(1).                   (7.2)
```

### Proof

For an integer `Q>=1`, consider all

```text
x in {0,1,...,Q}^r.
```

There are `(Q+1)^r` such points.  Every coordinate of `Ax` has absolute value at most `QrB`, so the number of possible image vectors is at most

```text
(2QrB+1)^k <= (3QrB)^k.
```

If

```text
Q^(r-k) > (3rB)^k,
```

then `(Q+1)^r>(3QrB)^k`, so two distinct points `x,y` have the same image.  Their difference `c=x-y` is nonzero, satisfies `Ac=0`, and has `||c||_infinity<=Q`.  Choosing `Q` just above `(3rB)^(k/(r-k))` gives `(7.1)`.

### Consequence for the actual Turan row

Section 11 below gives

```text
log max_{j<=r}|E_j| = O(n).
```

Applying Proposition 7.1 with `k=1` gives an exact nonzero relation

```text
sum_{j<=r} c_j E_j = 0
```

with

```text
log max|c_j| = O(n/(r-1) + log r/r).                     (7.3)
```

That is the **same exponential scale** as the coefficient box used to seek an `O(n/r)` nonzero combination.

Therefore the zero branch is not an LLL pathology and cannot be ruled out by exploiting “genericity” of the actual Apéry integers.  Integer arithmetic forces short exact relations at that scale for every long row.

More strongly, if one appends any fixed number `k` of extra linear carrier rows, all of characteristic-zero height `exp(O(n))`, Proposition 7.1 still gives a nonzero **common** exact kernel vector with

```text
log ||c||_infinity = O(k n/(r-k)).                        (7.4)
```

For fixed `k`, `(7.4)` is again `O(n/r)`.  Thus the strategy

```text
E row + one more jet row + simultaneous Minkowski
```

cannot by itself remove the zero-output phenomenon.  A finite tower of linear jets still has a short common kernel.

This does **not** disprove the existence of a small nonzero `H`; it proves that a proof cannot proceed by showing that the short kernel is empty.

---

## 8. The first discrete Wronskian and its exact obstruction

The first natural attempt to exploit `(6.2)` is

```text
W_j = R_{j+1}E_j - R_jE_{j+1}.                            (8.1)
```

If the two row vectors `(E_j)` and `(R_j)` were being used only over `Q`, a nonzero `W_j` would certify rank two.  But modulo every common candidate prime the two rows are rank one by `(6.2)`.  Consequently

```text
p | W_j   for every p in C(n,r),                          (8.2)
```

whether or not `p` is a target.

There is also an exact characteristic-zero identity.

### Proposition 8.1 -- adjacent first-jet Wronskian

For every `j`,

```text
W_j
 = D_j (D_{j-1}D_{j+1}-D_j^2)/(q_j q_{j+1}).            (8.3)
```

### Proof

Put

```text
a=D_{j-1}, b=D_j, c=D_{j+1}.
```

Then

```text
R_j=(a-b)/q_j,
R_{j+1}=(b-c)/q_{j+1}.
```

Let

```text
T_j=q_jE_j,
T_{j+1}=q_{j+1}E_{j+1}.
```

The exact adjacent-elimination identity already checked in the current Turan audit is

```text
(b-c)T_j - (a-b)T_{j+1} = b(ac-b^2).
```

Divide by `q_jq_{j+1}` to get `(8.3)`.

Thus the proposed discrete Wronskian is not a mysterious new carrier: it is exactly the next Turan expression in the terminal increments.

### Candidate-primorial height

Put

```text
P_cand(n,r) = product_{p in C(n,r)} p.
```

By `(8.2)`,

```text
P_cand(n,r) | W_j
```

for every `j<r`.

For `r=o(n)`, the prime number theorem gives

```text
log P_cand(n,r)
 = theta(n)-theta(n/2+O(r))
 = (1/2+o(1)) n.                                         (8.4)
```

This is much larger than

```text
n/r+r log n
```

when `r~sqrt(n/log n)`.  Therefore a nonzero universally candidate-divisible fallback cannot have the desired subcritical height.  It must first be deflated by the candidate content, and then one must show that targets retain another `p`-adic digit.

The next subsection shows that this extra digit is false already on actual Apéry data.

---

## 9. Actual-Apéry hostile certificate at `n=200`

I independently reimplemented the exact shell/Newton/Turan arithmetic in a stateless Wolfram kernel and then encoded the same checks in the owned standard-library verifier.

For

```text
n=200,
r=6,
M=199,
s=100,
```

condition `(4.1)` gives the common candidates

```text
109, 113, 127, 131, 137, 139, 149, 151, 157,
163, 167, 173, 179, 181, 191, 193, 197, 199.
```

The actual physical targets are exactly

```text
139, 181.
```

The common candidate primorial is

```text
3039274714336462506043636550560523672423.
```

For the first-jet Wronskians, the exact gcd has one and only one copy of every common candidate prime.  Dividing by the candidate primorial gives

```text
161372143389600
 = 2^5 * 3^2 * 5^2 * 11^2 * 13 * 37 * 59 * 61 * 107.
```

In particular,

```text
139 does not divide 161372143389600,
181 does not divide 161372143389600.
```

So candidate-content removal erases both genuine targets.  The independent Wolfram run also found valuation one for each of the first four individual `W_j` at both targets.

This is an **actual Apéry counterexample** to the tempting strengthening

```text
"a target gives a second p-adic digit in the R-Wronskian".
```

It is not a generic-lattice model.

---

## 10. Small gcd versus small Bezout coefficients

The same `n=200` row sharply separates the two notions required in the question.

The first two Turan carriers have

```text
bit_length(E_1) = bit_length(E_2) = 2112,
```

but

```text
gcd(E_1,E_2)
 = 503180
 = 2^2 * 5 * 139 * 181.                                  (10.1)
```

Thus the gcd is tiny compared with the two inputs and already contains the target product.

However, solve

```text
x E_1 + y E_2 = 503180.                                  (10.2)
```

The centered extended-Euclidean solution has

```text
bit_length(|x|)=2092,
bit_length(|y|)=2091.                                    (10.3)
```

This is not an artifact of a bad choice of Bezout pair.  Every solution of `(10.2)` has the form

```text
x = x0 + t E_2/g,
y = y0 - t E_1/g.
```

The exact computation gives

```text
2|x0| < |E_2|/g,
2|y0| < |E_1|/g.
```

Therefore every `t!=0` makes each shifted coordinate larger in absolute value; the displayed solution is coordinatewise minimal in the two-term family.

So `(10.1)` is a small **gcd carrier**, but `(10.3)` shows it is not a small-coefficient `H` of the requested type.  Any proof that silently turns a small gcd into a small Bezout combination is invalid.

For the first six `E_j` at this row the gcd is even smaller,

```text
251590 = 2 * 5 * 139 * 181,
```

but I do not claim a sharp six-variable Bezout lower bound from this finite computation.

---

## 11. Height ledger

No asymptotic claim below comes from numerical fitting.

### 11.1 Shell height

The Section-48 Laurent polynomial has nonnegative integer coefficients whose sum is `40`.  Therefore the sum of **all** coefficients of its `M`-th power is `40^M`.  The physical shell `C_M(d)` is a sub-sum, so

```text
0 < C_M(d) <= 40^M.                                      (11.1)
```

### 11.2 Newton values

A terminal Newton weight is

```text
(-1)^i binomial(start+i,i) binomial(M+1,L-i).
```

Each binomial factor is at most `2^n`, and there are at most `n` terms.  Hence a convenient crude bound is

```text
|F_j| <= n * 4^n * 40^(n-1) < n * 160^n.                (11.2)
```

Thus

```text
log |F_j| <= n log 160 + O(log n).                       (11.3)
```

### 11.3 Turan, first jet, and determinants

Since `q_j>=1`,

```text
|E_j| <= 2 (n 160^n)^2,                                 (11.4)
```

so

```text
log |E_j| <= 2n log160 + O(log n).                       (11.5)
```

From `(5.3)`,

```text
|R_j| <= 4 n 160^n,                                     (11.6)
```

and therefore both determinant families satisfy crude linear-exponential bounds

```text
log |W_j| = O(n),
log |J_{i,j}| = O(n).                                    (11.7)
```

The precise constants are irrelevant to the exponent calculation; what matters for Proposition 7.1 is that a fixed number of rows has `exp(O(n))` height.

### 11.4 Permissible `r`

We need

```text
2 <= r <= L0-1
```

and a nonempty common strip `(4.1)`.  Any

```text
r -> infinity,
r=o(n)
```

is compatible with the terminal construction.  In particular

```text
r ~ sqrt(n/log n)                                        (11.8)
```

is admissible.

At `(11.8)`,

```text
n/r + r log n = Theta(sqrt(n log n)).                    (11.9)
```

This is the allowed coefficient/clearing budget in the question.  By contrast `(8.4)` is `Theta(n)` and cannot be absorbed.

---

## 12. Why finitely many linear first-jet families do not repair Minkowski

It is tempting to append another target-preserving row---for example the original terminal values `(F_j)`, the first jet `(R_j)` after a suitable target factor, a shifted model, or a higher fixed-order jet---and ask simultaneous Minkowski to return a vector not in the common kernel.

Proposition 7.1 gives the exact obstruction.  If there are `k` fixed rows, each of height `exp(Cn)`, there is always a nonzero common kernel vector with

```text
log ||c||_infinity <= (kC+o(1)) n/(r-k).                 (12.1)
```

For fixed `k` and growing `r`, `(12.1)` is still `O(n/r)`.  Therefore no finite linear enlargement can prove that the short zero branch is absent.

One can of course take `k` growing with `r`, but then the simultaneous Dirichlet/Siegel exponent degrades with `k/(r-k)`.  Taking enough rows to make the common kernel trivial removes the high-dimensional cancellation that produced the `n/r` gain.

So the surviving mechanism must be one of the following:

1. a **nonlinear** invariant/determinantal divisor with target multiplicity, or
2. a genuine quotient-lattice covering theorem that bounds a nonzero image representative despite the short kernel, or
3. an arithmetic theorem directly bounding the gcd/Smith invariant.

The value--Turan determinant below is of type 1.

---

## 13. The value--Turan determinant: a target-square carrier without candidate pollution

For `1<=i<j<=r`, define

```text
J_{i,j}=E_iF_j-E_jF_i.                                   (13.1)
```

Equivalently, these are the `2 x 2` minors of

```text
A_FE(n,r) =
[ F_1 F_2 ... F_r ]
[ E_1 E_2 ... E_r ].                                     (13.2)
```

Let

```text
Delta_FE(n,r)=gcd_{i<j}|J_{i,j}|.                        (13.3)
```

This is the second determinantal divisor of the integer matrix `(13.2)` when its rank is two.

### Theorem 13.1 -- target-square divisibility

For every `p in T(n,r)`,

```text
p^2 | J_{i,j}  for every i<j.                            (13.4)
```

Consequently, if `Delta_FE(n,r)!=0`,

```text
rad(T(n,r))^2 | Delta_FE(n,r).                           (13.5)
```

### Proof

For a target, `(4.2)` gives

```text
F_i == 0 (mod p)
```

for every relevant column, and Corollary 6.2 gives

```text
E_i == 0 (mod p).
```

Each product in `(13.1)` is therefore divisible by `p^2`, proving `(13.4)` and `(13.5)`.

No new denominator is introduced: `E_i` is already an exact integer after the existing Pascal normalization.

### Theorem 13.2 -- non-target first-jet reduction

For every common candidate `p in C(n,r)`,

```text
J_{i,j}
 == b_{n-p}^2 (R_i-R_j)  (mod p).                        (13.6)
```

### Proof

Use

```text
F_i == b_{n-p},
E_i == b_{n-p} R_i
```

from `(4.2)` and `(6.2)` in `(13.1)`.

If `p` is a non-target, `b_{n-p}` is a unit, so

```text
p | Delta_FE(n,r)
```

among common candidates exactly when

```text
R_1 == R_2 == ... == R_r  (mod p).                       (13.7)
```

There is therefore no formal candidate-primorial factor analogous to `(8.2)`.

### Endpoint check

At `p=n`, `b_0=1`, so neither row in `(13.2)` is forced to vanish.  The target-square argument does not include the endpoint, exactly as required.

---

## 14. Independent exact `Delta_FE` records

The following computations were performed with an independent exact Wolfram implementation of the shell/Newton formulas; the owned verifier reproduces them with standard-library Python.

For `r=6`:

```text
n=200
common targets = (139,181)
Delta_FE = 31648764050
         = 2 * 5^2 * 139^2 * 181^2
bit length = 35

n=272
common targets = (191,233)
Delta_FE = 2 * 191^2 * 233^2
bit length = 32

n=300
common targets = (191,227)
Delta_FE = 37 * 191^2 * 227^2
bit length = 37

n=321
common targets = (179,193,211)
Delta_FE = 179^2 * 193^2 * 211^2
bit length = 46
```

For comparison, the individual `E_j` at `n=200` have about `2112` bits.  Thus the second determinantal divisor is spectacularly smaller in these rows.

Again: these are **regression/candidate data, not an asymptotic proof**.  No step below uses the observed factorization pattern as if it were a theorem.

---

## 15. The exact surviving lemma

The strongest theorem I can honestly isolate from the present source state is the following one-step target.

### [FE-SMITH]

There exists a function `r=r(n)` with

```text
r -> infinity,
r=o(n),
```

such that for all sufficiently large `n`, the matrix `A_FE(n,r)` in `(13.2)` has rank two over `Q` and

```text
0 < Delta_FE(n,r)
    <= exp(C (n/r + r log n))                            (15.1)
```

for an absolute constant `C`.

This is **not** the desired radical statement in different notation.  It is an explicit second Smith invariant of a concrete characteristic-zero matrix built from the already-defined terminal value row and once-divided Turan row.  Its local behavior distinguishes targets from non-target candidates by `(13.4)` and `(13.6)`.

If `(15.1)` holds, `(13.5)` gives immediately

```text
2 log rad(T(n,r))
 <= log Delta_FE(n,r)
 = O(n/r+r log n).                                       (15.2)
```

With `r` as in `(11.8)`,

```text
log rad(T(n,r)) = O(sqrt(n log n)).                      (15.3)
```

The omitted lower boundary strip of width `O(r)` is exactly the part that must be handled by the surrounding top-half decomposition in the current many-prefix setup; I do not claim a new estimate for that missing strip from the connector-invisible files.

### Why this is smaller than a generic resultant request

`Delta_FE` is a determinantal divisor of `r` explicit integer columns; it does not introduce a degree-`Theta(n)` polynomial resultant.  It also requires no target-dependent choice and no division by an unknown radical.  The non-target reduction `(13.6)` is first-jet arithmetic, not a tautological statement that targets divide a gcd.

### What remains genuinely hard

A proof of `(15.1)` needs a **characteristic-zero Smith/content theorem across the growing terminal prefix**.  The finite factorizations suggest that the nuisance part may be only `exp(O(r log n))`, but no visible source proves this, and I do not infer it from the table.

One plausible exact formulation is to prove that the primitive second Smith invariant after removing its forced top-target square has only primes arising from a controlled product of `O(r)` Pascal/boundary factors.  That would give the allowed `O(r log n)` nuisance.  The currently visible identities do not establish such a factorization.

---

## 16. Why the original small-`H` theorem remains open

The desired statement asks for

```text
H = sum_{j<=r} c_j E_j != 0
```

with both

```text
log |H| = O(n/r+r log n),
log max|c_j| = O(n/r+r log n).                           (16.1)
```

Nothing above proves `(16.1)`.

Proposition 7.1 shows why the usual “Minkowski gives a short vector” proof is incomplete: the same coefficient box necessarily contains exact zero relations.  A quotient-lattice proof would have to control the **covering radius of the primitive image lattice**, not merely its determinant or gcd.

The actual `n=200` Bezout calculation in Section 10 shows that a small gcd does not automatically provide this covering control.

Thus there are two logically distinct possible closures:

```text
(A) prove the original bounded-coefficient nonzero H theorem;
```

or

```text
(B) prove [FE-SMITH], which closes the radical estimate without producing H.
```

The present attack advances `(B)` and gives a theorem-level reason that the naive kernel-exclusion version of `(A)` cannot work.

---

## 17. Mutation tests

The owned verifier includes several hostile mutations implicitly as assertions.

### Mutation 1: drop the Pascal gcd

If `q_j` is replaced by an arbitrary smaller normalization, exact integrality of

```text
(F_{j-1}F_{j+1}-F_j^2)/q_j
```

or the primitive identities `(5.4)`/`(6.1)` fail.  The verifier checks every division before using it.

### Mutation 2: treat `p=n` as a target

The verifier constructs the candidate interval with `p<n` and checks `b_0=1` at prime endpoints.  A proof that silently includes `p=n` as a target is rejected.

### Mutation 3: assume the `R`-Wronskian has an extra target digit

At `n=200`, after removing the full common candidate primorial from the gcd of the `W_j`, the exact quotient is

```text
161372143389600,
```

which is divisible by neither `139` nor `181`.  This directly rejects the extra-digit claim.

### Mutation 4: replace `Delta_FE` by universal candidate content

For `n=200`, `Delta_FE` is

```text
2 * 5^2 * 139^2 * 181^2
```

and none of the other common candidates divides it.  Thus the target-square determinant is not merely the candidate primorial in disguise.

---

## 18. Reproduction

The owned verifier is

```text
problems/3.2/ORACLE_COMM/chatgpt_turan_nonzero_certificate_verify.py
```

It has no imports from repository code and uses only the Python standard library.  The default run reconstructs `n=200,r=6` and checks:

```text
* terminal Newton values and Pascal q_j integrality;
* the exact first-jet identity (6.1);
* the primitive A_j,B_j formula (3.2)/(5.4);
* the Wronskian identity (8.3);
* every common candidate congruence (4.2), (4.3), (6.2);
* exact target selection with p=n excluded;
* universal candidate divisibility of W;
* loss of 139 and 181 after candidate deflation;
* target-square and non-target formulas for every FE minor;
* Delta_FE(200,6)=2*5^2*139^2*181^2;
* gcd(E_1,E_2)=503180;
* the 2092/2091-bit coordinatewise-minimal two-term Bezout coefficients.
```

Run

```bash
python3 problems/3.2/ORACLE_COMM/chatgpt_turan_nonzero_certificate_verify.py
```

For the four exact `Delta_FE` records, run

```bash
python3 problems/3.2/ORACLE_COMM/chatgpt_turan_nonzero_certificate_verify.py --extended
```

For a dense identity/congruence sweep over `12<=n<=60`, add

```bash
--small-sweep
```

The independent Wolfram run used to cross-check the new formulas produced the exact records stated in Sections 9, 10, and 14.  No floating-point evidence is used in any asserted identity.

---

## 19. Final assessment

The requested small nonzero linear combination is **not proved** from the connector-visible source state, and I found no actual-Apéry family disproving it.

What is proved here is sharper than saying “Minkowski might return zero”:

1. **Short exact zero relations necessarily exist at the same `exp(O(n/r))` scale** for the actual integer row.  Hence a short-kernel lower-bound attack is impossible.
2. The exact first-jet identity is

   ```text
   E_j = F_j R_j - D_{j-1}D_j/q_j,
   E_j == b_{n-p}R_j (mod p).
   ```

3. The corresponding `R`-Wronskian is universally candidate-divisible and, on actual Apéry data at `n=200`, has only one target digit.  Candidate deflation destroys the targets.
4. The **value--Turan second determinantal divisor** escapes that contamination:

   ```text
   p target  =>  p^2 | Delta_FE,
   p non-target candidate  =>
       J_{i,j} == b_{n-p}^2(R_i-R_j) (mod p).
   ```

5. Exact rows `n=200,272,300,321` give target-square `Delta_FE` with only tiny nuisance, but this finite pattern is not promoted to an asymptotic theorem.

Therefore the sharp verdict is:

> **SMALLEST EXACT SURVIVING LEMMA = [FE-SMITH].**  Prove a subcritical height bound for the second determinantal divisor of the concrete `(F,E)` prefix matrix.  This would close the hardest common top strip without any Bezout-coefficient theorem.  The original bounded-coefficient `H` statement remains open, and any proof of it must control a quotient-lattice covering radius rather than try to eliminate short kernel relations.
