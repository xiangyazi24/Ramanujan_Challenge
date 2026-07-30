# RETRACTED: the apparent order-seven scalar barrier

**Correction (2026-07-29).**  The linear-rank conclusion below is
false because `W_p` is not independent:

```text
W_p == (24*769/5)(Delta_p/p^3)^2 mod p.
```

Multiplying the direct row by `Delta_p^2` removes it and proves the
fixed law

```text
(1680+2472Delta_p^2)D_7-25Z_7
 ==1655 b_n/p mod p^7.
```

See `Q32_ORDER_SEVEN_TARGET_FIXED_LAW_2026-07-29.md`.  The remainder
of this file is retained only as the audit trail for the two linear
rows; its claim that `W_p` is a free third coordinate must not be
cited.

---

# Original linear audit

Date: 2026-07-29. Owner: Codex.

## 0. Result

Assume the seventh endpoint theorem from
`Q32_WEIGHT_SEVEN_ENDPOINT_RANK_ONE_2026-07-29.md`.  Put

```text
Delta=b_(p-1)-1,
H=b_p-5+7Delta,

W=[b_(2p)-73+824Delta-(752/5)H]/p^6              (mod p).
```

For `p>=11`, `p!=769`, its four endpoint specializations are

```text
b_(p-1) ==1+Delta,
b_p     ==5-7Delta+H,

b_(2p)  ==73-824Delta+(752/5)H+p^6W,
b_(2p-1)==5+8Delta+(336/5)H-(103/769)p^6W        (mod p^7).
                                                               (0.1)
```

Let `0<=r<p`, assume `p|b_r`, and put

```text
n=p+r,               s=p-1-r,               x=b_n/p.
```

There are explicit direct and reflected residues `D_7,Z_7` such
that

```text
D_7==x(1-H/5),                                      (0.2)

Z_7==x(1-336H/25+(103/(5*769))p^6W)                (mod p^7).
                                                               (0.3)
```

The old weights still eliminate `H`, but now give

```text
336D_7-5Z_7
 ==x(331-(103/769)p^6W)                            (mod p^7).
                                                               (0.4)
```

Thus the endpoint rank-one theorem does **not**, by itself, furnish
a seventh fixed-coefficient quotient digit.  The new scalar `W`
occurs in the reflected leading endpoint and has no direct
counterpart because the direct carrier satisfies `R_1=0`.
Modulo `p^6`, `(0.4)` reduces to the previous `331x` law.  At the
seventh digit, one must either import `W`, obtain a third
controlled-height description of the same target, or prove a
cross-prime constraint on the `W_p`.

If one is allowed to import the local endpoint scalar `W`, then the
factor on the right of `(0.4)` is a unit for
`p notin {331,769}` and `x mod p^7` is recoverable.  This is only a
local reconstruction, not the fixed-height carrier needed for the
global radical bound.

## 1. Shifted recurrence notation

Let `mathcal U_r(X),mathcal V_r(X)` be the shifted fundamental
solutions and put

```text
mathcal J_r(X)=mathcal V_r(X)/(1+X)^3.
```

The exact direct and reflected decompositions are

```text
b_(p+r)
 =b_p mathcal U_r(p)-p^3b_(p-1)mathcal J_r(p),    (1.1)

b_(2p-1-s)
 =b_(2p-1)mathcal U_s(-2p)
  +8p^3b_(2p)mathcal J_s(-2p).                    (1.2)
```

Notice that the first factor in `(1.2)` is `b_(2p-1)`, not
`b_(2p)`.

For order seven after division by `p`, take

```text
T_r^[7](X)=mathcal U_r(X) mod X^8,
J_r^[4](X)=mathcal J_r(X) mod X^5.
```

Define

```text
D_7=
 [(5-7Delta)T_r^[7](p)
  -p^3(1+Delta)J_r^[4](p)]/p                     (mod p^7),
                                                               (1.3)

Z_7=
 [(5+8Delta)T_s^[7](-2p)
  +8p^3(73-824Delta)J_s^[4](-2p)]/p              (mod p^7).
                                                               (1.4)
```

The truncation remainders in `(1.3)` and `(1.4)` are multiples of
`p^8` before division.  The target assumptions ensure that the
displayed numerators are divisible by `p`.

## 2. Direct equation

Set

```text
q=mathcal U_r(p)/p.
```

Divide `(1.1)` by `p` and use `(0.1)`:

```text
x=D_7+Hq                                           (mod p^7). (2.1)
```

Reducing the same divided recurrence modulo `p^2` gives

```text
x==5q                                             (mod p^2), (2.2)
```

because `p^2mathcal J_r(p)` vanishes modulo `p^2` and
`Delta,H` have valuations at least `3,5`.  Since `p^5|H`, `(2.2)`
may be substituted into `(2.1)`, proving `(0.2)`.

## 3. Reflected equation

Put

```text
q'=mathcal U_s(-2p)/p.
```

In `(1.2)`, the `H` term in `b_(2p)` is multiplied by `8p^3` and
therefore has valuation at least `8` before division.  The `W` term
there has still larger valuation.  Only the `H,W` corrections in
`b_(2p-1)` survive.  Dividing by `p` gives

```text
x=Z_7+(336/5)Hq'-(103/769)p^6Wq'                 (mod p^7).
                                                               (3.1)
```

Again, reduction modulo `p^2` gives

```text
x==5q'                                            (mod p^2). (3.2)
```

For the `H` term we need `(3.2)` modulo `p^2`; for the `p^6W`
term only its reduction modulo `p` is needed.  Substitution in
`(3.1)` proves `(0.3)`, and `336*(0.2)-5*(0.3)` proves `(0.4)`.

## 4. Information content

The two normalized equations have three coordinates

```text
x,                  Hx,                  p^6Wx.
```

Their coefficient matrix is

```text
[ 1, -1/5,       0             ],
[ 1, -336/25,    103/(5*769)   ].
```

Hence no fixed linear combination of the two rows can retain `x`
while killing both endpoint scalars.  The lower-order equations are
their reductions modulo `p^6` and add no third row.  This is an
exact rank obstruction, not a numerical failure.

The obstruction identifies the next viable inputs precisely:

```text
1. a third short-gap representation of b_n at the same p;
2. a controlled-height law for W_p across target primes;
3. a same-kernel relation coupling the quotient digits directly;
4. an Archimedean/non-stabilization theorem for the reconstructed
   cofactor tower.
```
