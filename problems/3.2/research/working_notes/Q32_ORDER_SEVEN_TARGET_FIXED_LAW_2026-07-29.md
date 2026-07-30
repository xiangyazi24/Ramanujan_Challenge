# The fixed order-seven target law

Date: 2026-07-29.  Owner: Codex.

This retracts the independence conclusion in
`Q32_ORDER_SEVEN_TARGET_BARRIER_2026-07-29.md`.  The two linear rows
do contain a scalar `W_p`, but H6 makes that scalar a quadratic
function of the already present endpoint coordinate.  A nonlinear
combination removes it.

## 1. Endpoint notation

For a prime `p>=11`, put

```text
Delta=b_(p-1)-1,
H=b_p-5+7Delta,

X=p^2 S_2,
Y=p^4 S_4,
Z=p^5 sum_(k=1)^(p-1) H_k^(2)/k^3.
```

The seventh-grade block calculation gives

```text
Delta == X+5Y-X^2+4Z                            (mod p^7),    (1.1)

E := b_(2p)-73+824Delta-(752/5)H
   == (24/5)(935X^2-830Y-332Z)                  (mod p^7).    (1.2)
```

Define

```text
W=E/p^6 mod p.
```

The now-proved H6 identity is

```text
X^2-5Y-2Z==0                                    (mod p^7).    (1.3)
```

## 2. `W` is a square, not a free digit

Substitute (1.3) in (1.2):

```text
935X^2-830Y-332Z
 ==935(5Y+2Z)-830Y-332Z
 ==3845Y+1538Z
 ==769(5Y+2Z)
 ==769X^2                                       (mod p^7).
```

Consequently

```text
W == (24*769/5)(X/p^3)^2                       (mod p).       (2.1)
```

The valuations are

```text
v_p(X)>=3,       v_p(Y)>=5,       v_p(Z)>=5.
```

Thus (1.1) gives

```text
Delta/p^3 == X/p^3                              (mod p),
```

and (2.1) becomes

```text
W == (24*769/5)(Delta/p^3)^2                   (mod p).       (2.2)
```

Equivalently,

```text
p^6W == (24*769/5)Delta^2                       (mod p^7).     (2.3)
```

This is the point missed by the earlier linear rank audit: `W` is
not an independent third coordinate after the nonlinear endpoint
relations are used.

## 3. Fixed target elimination

Let `p|b_r`, set

```text
n=p+r,        s=p-1-r,        x=b_n/p,
```

and let `D_7,Z_7` be the direct and reflected target residues from
the old note.  Their verified equations are

```text
D_7==x(1-H/5),                                  (mod p^7),    (3.1)

Z_7==x(1-336H/25+(103/(5*769))p^6W)             (mod p^7).    (3.2)
```

The old linear combination gives

```text
336D_7-5Z_7
 ==x(331-(103/769)p^6W)
 ==x(331-(2472/5)Delta^2)                       (mod p^7),    (3.3)
```

where the second equality is (2.3).

Now `v_p(Delta^2)>=6` and `v_p(H)>=5`, so (3.1) also gives

```text
Delta^2D_7==Delta^2x                             (mod p^7).    (3.4)
```

Add `(2472/5)` times (3.4) to (3.3).  The result is the fixed law

```text
(336+(2472/5)Delta^2)D_7-5Z_7
 ==331x                                          (mod p^7),    (3.5)
```

or, with integral coefficients,

```text
(1680+2472Delta^2)D_7-25Z_7
 ==1655 b_n/p                                    (mod p^7).    (3.6)
```

The right coefficient is a unit away from `p=331`.  The existing
endpoint normalization omits `p=769`, and H6 omits `p=7`; all three
are fixed global exceptions with logarithmic cost `O(1)`.

For height accounting, `Delta^2 mod p^7` need not be represented by
the exponentially large endpoint integer.  If

```text
delta_p=(Delta/p^3) mod p,       0<=delta_p<p,
```

then

```text
Delta^2 == p^6 delta_p^2                         (mod p^7),
```

whose representative has polynomial height.  Thus (3.6) does not
import a separate free endpoint digit.

## 4. Consequence and remaining gap

The order-seven local tower does **not** stop at `W_p`: the fixed
quotient carrier persists for one more digit.  This is a genuine
positive correction to the previous barrier note.

It does not yet prove

```text
log R_n^(>=7)=o(n).
```

The global task is now the higher-order version: prove that each new
endpoint scalar is algebraic over the previously exposed normalized
coordinates in a way that permits the same small-height nonlinear
elimination, and combine an unbounded number of grades with a
non-stabilization argument for the common quotient.

## 5. Audit

The updated independent reproducer is

```text
../scripts/q32_order_seven_target_audit.py
```

For all primes through `1000`, excluding the normalization prime
`769`, it reports

```text
endpoint W-square identity             163/163
direct target equation                 161/161
reflected target equation              161/161
old linear combination                 161/161
fixed nonlinear combination            161/161
small coefficient representative       161/161
symbolic W-square reduction               1/1
symbolic nonlinear elimination            1/1
failures                                     0
```
