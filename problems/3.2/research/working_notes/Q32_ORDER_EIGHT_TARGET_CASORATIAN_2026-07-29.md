# The precision-eight target Casoratian law

Date: 2026-07-29. Owner: Codex.

## 0. Result

The direct/reflected target construction extends one further digit.
The new target-dependent jet is not arbitrary: it is the same
companion value on the two reflected zero rows, and a fixed combination
cancels it.

Let `p>=11`, and put

```text
Delta=b_(p-1)-1,
H=b_p-5+7Delta,

E=b_(2p)-73+824Delta-(752/5)H,
F=b_(2p-1)-5-8Delta-(336/5)H.
```

The proved raw endpoint laws and `H6` identity imply

```text
E,F in p^6 Z_(p),                 769F+103E in p^7 Z_(p).
```

Indeed, if `delta_0=Delta/p^3 mod p`, their first normalized digits are

```text
E/p^6 == (24*769/5)delta_0^2,
F/p^6 ==-(24*103/5)delta_0^2                    (mod p),
```

so the displayed linear combination cancels without dividing by
`769`.

Thus `w=E/p^6` is integral for every `p>=11`.  When `p!=769`,

```text
v=[F+(103/769)E]/p^7
```

is also integral, but this split normalization is not needed in the
theorem.

Suppose `p|b_r`, and set

```text
s=p-1-r,                 x=b_(p+r)/p.
```

Let `U_j(X),V_j(X)` be the shifted fundamental solutions of the
Apéry recurrence, and set

```text
J_j(X)=V_j(X)/(1+X)^3.
```

Use `U` through degree eight and `J` through degree five before
dividing by `p`, and define

```text
D_8=[
 (5-7Delta)U_r(p)-p^3(1+Delta)J_r(p)
]/p                                                   (mod p^8),

Z_8=[
 (5+8Delta)U_s(-2p)
 +8p^3(73-824Delta)J_s(-2p)
]/p                                                   (mod p^8).
```

Then

```text
D_8
 ==x(1-H/5)-(p^2H/5)J_r(0)                           (mod p^8), (0.1)

Z_8
 ==x(1-336H/25-F/5)
   +(166144/25)p^2H J_s(0)                            (mod p^8). (0.2)
```

At the target zero,

```text
J_s(0)==J_r(0)                                        (mod p).   (0.3)
```

Consequently the companion coordinate cancels:

```text
166144D_8+5Z_8
 ==x(166149-33296H-F)                                 (mod p^8).
                                                                  (0.4)
```

Since the correction in parentheses after `166149` has valuation at
least five and `D_8==x (mod p^5)`, it can be moved to the left:

```text
[166144+33296H+F]D_8+5Z_8
 ==166149x                                             (mod p^8). (0.5)
```

This is denominator-free at the exceptional normalization prime
`p=769`.  Clearing the harmless denominator `5` in `F` gives the
ordinary integer-coefficient form

```text
[5b_(2p-1)+830695-40Delta+166144H]D_8+25Z_8
 ==830745x                                             (mod p^8). (0.6)
```

For `p!=769`, substituting

```text
F=-(103/769)p^6w+p^7v
```

recovers the original normalized `w,v` form exactly.

Here

```text
166149=3^2*18461.
```

Thus the right coefficient is a fixed unit outside `p=18461` in the
stated range.  The congruence itself remains valid at `p=18461`; only
its inversion fails there.

Equation `(0.5)` is an unconditional local quotient law, not a global
separation theorem.  Its new endpoint datum is the two-digit residue

```text
F/p^6 mod p^2.
```

Its first digit is fixed by `H6`, while its second digit is new.  These
residues vary with `p` and do not constitute one
characteristic-zero carrier for all target primes.

## 1. Exact shifted decompositions

The recurrence gives

```text
b_(p+r)
 =b_p U_r(p)-p^3b_(p-1)J_r(p),                       (1.1)

b_(2p-1-s)
 =b_(2p-1)U_s(-2p)+8p^3b_(2p)J_s(-2p).               (1.2)
```

The left side of `(1.2)` is `b_(p+r)`.  Put

```text
q =U_r(p)/p,                 q'=U_s(-2p)/p.
```

The target congruences make both quotients `p`-integral.  Dividing
`(1.1)` by `p` and using the definitions of `Delta,H` gives exactly

```text
x=D_8+Hq                                               (mod p^8). (1.3)
```

Modulo `p^3`,

```text
x==5q-p^2J_r(0),
```

because `v_p(Delta)>=3`.  Therefore

```text
q==x/5+(p^2/5)J_r(0)                                  (mod p^3).
```

Substitution in `(1.3)` proves `(0.1)`.

For the reflected row, the definitions of `E,F` are the exact endpoint
expansions

```text
b_(2p)
 =73-824Delta+(752/5)H+E,

b_(2p-1)
 =5+8Delta+(336/5)H+F.                                (1.4)
```

After division of `(1.2)` by `p`,

```text
x==5q'+584p^2J_s(0)                                   (mod p^3),
```

so

```text
q'==x/5-(584/5)p^2J_s(0)                              (mod p^3). (1.5)
```

In the endpoint corrections:

- the `H` coefficient needs `(1.5)` modulo `p^3`;
- the `F` coefficient needs only `q'==x/5 (mod p^2)`;
- the `E` correction in `b_(2p)` is multiplied by `p^2`
  and vanishes modulo `p^8`.

The companion coefficient is

```text
(336*584)/25-(8*752)/5
 =196224/25-30080/25
 =166144/25,
```

which proves `(0.2)`.

The degree-eight truncation of `U` has error `O(p^9)` before
division.  The degree-five truncation of `J` is multiplied by `p^3`,
so its error is also `O(p^9)`.  Both become `O(p^8)` after division,
as required.

## 2. The companion reflection is a Casoratian identity

At `X=0`, `U_j(0)=b_j`, while the companion has initial values

```text
V_0(0)=0,                    V_1(0)=1.
```

Its Casoratian is

```text
U_k(0)V_(k+1)(0)-U_(k+1)(0)V_k(0)=1/(k+1)^3.          (2.1)
```

No target can have `r=0` or `s=0`, because `b_0=1` and
`b_(p-1)==1 (mod p)`.  Hence all denominators below are units.
At `b_r==0 (mod p)`, equation `(2.1)` at `k=r-1` gives

```text
V_r(0)==1/[r^3b_(r-1)]                                (mod p). (2.2)
```

For `s=p-1-r`, reflection gives

```text
b_(s-1)==b_(r+1)                                      (mod p).
```

Using the recurrence at the zero `b_r`,

```text
(r+1)^3b_(r+1)==-r^3b_(r-1)                           (mod p),
```

and using `(2.1)` at the reflected zero yields

```text
V_s(0)
 ==1/[s^3b_(s-1)]
 ==-1/[(r+1)^3b_(r+1)]
 ==1/[r^3b_(r-1)]
 ==V_r(0)                                             (mod p).
```

Since `J_j(0)=V_j(0)`, this proves `(0.3)`.  The central case
`r=s=(p-1)/2` is included and is tautological at the final step.

## 3. Information and height audit

The new coordinate ledger is

```text
endpoint coordinate:       F_p/p^6 mod p^2,
target companion:          p^2H_p J_r(0).
```

Reflection kills the second coordinate with fixed coefficients.  The
first is absorbed by adding `F_p` to the coefficient of the direct
row.  This proves local recovery of `x mod p^8` away from the fixed
right-coefficient exception.  The pair `w_p,v_p` is only a split
normalization of this single residual, and its failure at `p=769` is
artificial.

It does not give a global nonzero integer of sublinear height.  Across
all target primes the laws still recover the actual common integer
`b_n` in its CRT coordinates.  If `R` is the target radical, the
canonical reconstructed carrier is

```text
b_n-(b_n mod R^9)
 =R^9 floor(b_n/R^9).
```

It is zero once `R^9>b_n`, and before that yields only the fixed-order
linear bound.  A global proof still needs a transverse nonzero lift or
a carrier whose height per acquired power tends to zero.

## 4. Exact audit

The independent reproducer is

```text
../scripts/q32_order_eight_target_audit.py
```

For all target rows at primes through `1000`, it reports

```text
all target rows                            163
target rows at p=769                         2
direct/reflected equation checks       163+163
fixed and integer target-law checks    163+163
exact shifted decompositions               326
endpoint integrality checks                 492
universal companion-reflection checks     76110
symbolic eliminations                         3
failures                                      0
```

The only fixed inversion exception relevant to `(0.5)` is

```text
18461.
```

The `w,v` rewrite alone excludes `769`; the denominator-free law does
not.
