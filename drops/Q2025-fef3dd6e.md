ANSWER Q2025 fef3dd6e

# P3.2 depth-three seam: exact jet reduction and exhaustive strict-top audit

## Verdict

**Outcome 3.** I do not have an all-prime proof that

\[
(E):\qquad v_p(b_{p+r})=3,\qquad v_p(p^3a_{p+r})\ge 6
\]

is impossible, and I found no counterexample. I can, however, reduce `(E)` rigorously to explicit lower-digit jets and one exact three-digit Wronskian scalar, and I carried out an exhaustive exact scan over every strict-top target for every prime

\[
7\le p\le 10000.
\]

The scan contains `1226` primes and `1172` strict-top target pairs. It finds **no instance of `(E)`**. More importantly, it refutes a tempting but false stronger identity: the first divided companion/Gessel jet can vanish at a target. Exactly five such first-jet cancellations occur in this range:

\[
\boxed{
(p,r,n)=
(41,30,71),\ (97,25,122),\ (151,14,165),\
(1453,1180,2633),\ (6781,3974,10755).
}
\]

Every one of these five lifts exactly to

\[
\boxed{
v_p(b_{p+r})=v_p(p^3a_{p+r})=2,
}
\]

so none even reaches the valuation-three entrance to `(E)`.

The main exact reduction is as follows. Put

\[
A_{p,r}:=p^3a_{p+r}\in\mathbf Z_p,
\qquad B_{p,r}:=b_{p+r}\in\mathbf Z,
\]

and let `G_r,H_r` be the first two translation jets defined below. Then for every prime `p>=7` and every strict digit `0<=r<=p-2`, both shifted Apéry coordinates have the **same quadratic jet**:

\[
\boxed{
B_{p,r}\equiv 5\bigl(b_r+pG_r+p^2H_r\bigr)\pmod{p^3},
}
\tag{Q1}
\]

\[
\boxed{
A_{p,r}\equiv 6\bigl(b_r+pG_r+p^2H_r\bigr)\pmod{p^3}.
}
\tag{Q2}
\]

At a target `p|b_r`, write `c=b_r/p`. The first two exact cancellation scalars are

\[
\theta_1(p,r)=c+G_r\pmod p,
\tag{J1}
\]

and, if `theta_1=0`,

\[
\theta_2(p,r)=\frac{c+G_r}{p}+H_r\pmod p.
\tag{J2}
\]

Here the quotient in `(J2)` is taken in `Z_p`; all denominators of `G_r,H_r` are `p`-units because `r<p-1`.

They give the exact ladder

\[
\theta_1\ne0
\Longrightarrow
v_p(A_{p,r})=v_p(B_{p,r})=1,
\]

\[
\theta_1=0,\ \theta_2\ne0
\Longrightarrow
v_p(A_{p,r})=v_p(B_{p,r})=2,
\]

and only

\[
\theta_1=\theta_2=0
\]

can reach valuation at least `3`.

Once `v_p(B_{p,r})=3`, the exact Wronskian reduces the whole seam to one scalar. Set

\[
c_3:=\frac{B_{p,r}}{p^3}\in\mathbf Z_p^\times,
\quad
A_-:=p^3a_{p+r-1},
\quad
B_-:=b_{p+r-1}.
\]

Then `A_-` and `B_-` are `p`-units and

\[
\boxed{
\Omega_{p,r}:=A_-c_3+\frac6{(p+r)^3}\in\mathbf Z_p
}
\tag{O1}
\]

satisfies

\[
\boxed{
\frac{A_{p,r}}{p^3}=\frac{\Omega_{p,r}}{B_-}.
}
\tag{O2}
\]

Consequently

\[
\boxed{
(E)\iff
v_p(B_{p,r})=3\quad\text{and}\quad
\Omega_{p,r}\equiv0\pmod{p^3}.
}
\tag{O3}
\]

This is not the old reformulation `x<y`: `(Q1)--(J2)` explicitly computes the first two lower digits, and `(O1)` is the first remaining depth-three cancellation scalar. Its first residue is

\[
\Omega_{p,r}\equiv
6b_{r-1}c_3+\frac6{r^3}\pmod p.
\tag{O4}
\]

Thus a putative seam point must first satisfy

\[
\boxed{
c_3\equiv-\frac{a_r}{6}\pmod p,}
\tag{O5}
\]

where `a_r` is reduced modulo `p` using its `p`-unit denominator. The next two digits of `Omega` are the first genuinely unresolved data. I found no identity in the current repository that forces either of them to be nonzero.

---

## 1. Normalization and denominator bookkeeping

The Apéry recurrence is

\[
(m+1)^3y_{m+1}-P(m)y_m+m^3y_{m-1}=0,
\qquad
P(m)=34m^3+51m^2+27m+5,
\tag{1.1}
\]

with

\[
b_0=1,\ b_1=5,
\qquad
a_0=0,\ a_1=6.
\]

The exact Casoratian convention is

\[
\boxed{
a_m b_{m-1}-a_{m-1}b_m=\frac6{m^3}.}
\tag{1.2}
\]

Fix a prime `p>=7` and a strict top digit

\[
n=p+r,
\qquad p>\frac{n+1}{2}.
\]

The strict inequality is equivalent to

\[
0\le r\le p-2.
\tag{1.3}
\]

Apéry's denominator theorem gives

\[
\operatorname{den}(a_m)\mid\operatorname{lcm}(1,\ldots,m)^3.
\]

Since `p+r<2p`, the reduced denominator of `a_{p+r}` contains at most `p^3`; therefore

\[
A_{p,r}:=p^3a_{p+r}\in\mathbf Z_p.
\tag{1.4}
\]

For every lower index `j<p`, `a_j` has a `p`-unit denominator. Hence every reduction of `a_j`, and of the rational jets below whose recurrence denominators use only integers at most `j`, is legitimate in `Z_p`.

The endpoint supercongruences used below are the standard Apéry ones

\[
b_p\equiv5\pmod{p^3},
\qquad
b_{p-1}\equiv1\pmod{p^3}.
\tag{1.5}
\]

They are the same endpoint inputs already used in the repository's first-block/Dwork arguments.

---

## 2. Two explicit translation jets

Write

\[
(\mathcal L f)_r
=(r+1)^3f_{r+1}-P(r)f_r+r^3f_{r-1}.
\]

Define the first jet by

\[
G_0=0,\qquad G_1=12,
\]

and

\[
\boxed{
\mathcal LG_r
=P'(r)b_r
-3(r+1)^2b_{r+1}
-3r^2b_{r-1}.
}
\tag{2.1}
\]

This is Gessel's derivative sequence; equivalently

\[
G_r=
2\sum_{k=0}^r
\binom rk^2\binom{r+k}{k}^2
\bigl(H_{r+k}^{(1)}-H_{r-k}^{(1)}\bigr).
\tag{2.2}
\]

To avoid confusing the second jet with harmonic numbers, denote it here by `H_r` only inside this section. Set

\[
H_0=H_1=0
\]

and

\[
\boxed{
\begin{aligned}
\mathcal LH_r={}&
P'(r)G_r
-3(r+1)^2G_{r+1}
-3r^2G_{r-1}\\
&+\frac{P''(r)}2b_r
-3(r+1)b_{r+1}
-3rb_{r-1}.
\end{aligned}}
\tag{2.3}
\]

Because the recurrence solves each next value by division by `(r+1)^3`, every prime in the reduced denominator of `G_r` or `H_r` is at most `r`. In particular both are `p`-integral for every strict digit `(1.3)`.

Now put

\[
F_r(X)=b_r+XG_r+X^2H_r.
\tag{2.4}
\]

Expanding the shifted recurrence gives the exact polynomial identity modulo `X^3`:

\[
\boxed{
(r+X+1)^3F_{r+1}(X)
-P(r+X)F_r(X)
+(r+X)^3F_{r-1}(X)
\equiv0\pmod{X^3}.
}
\tag{2.5}
\]

Indeed the coefficients of `X^0,X^1,X^2` are respectively the recurrences for `b,G,H` above.

For reference, the first values are

\[
G_0=0,\quad G_1=12,\quad G_2=210,\quad G_3=4438,
\]

\[
H_0=H_1=0,\quad H_2=72,\quad H_3=2160.
\tag{2.6}
\]

A useful consistency correction: if one compares with an older workspace jet note that wrote `J_3=24552` for the first companion jet, that sample is inconsistent with its displayed recurrence. The exact first companion jet is

\[
J_r=6G_r,
\]

so

\[
J_3=6\cdot4438=26628.
\]

No conclusion below uses the erroneous sample.

---

## 3. Shared quadratic block theorem

### Theorem 3.1
For every prime `p>=7` and every `0<=r<=p-2`, in `Z_p`,

\[
B_{p,r}=b_{p+r}
\equiv5F_r(p)\pmod{p^3},
\tag{3.1}
\]

and

\[
A_{p,r}=p^3a_{p+r}
\equiv6F_r(p)\pmod{p^3}.
\tag{3.2}
\]

### Proof
For the `b`-solution, `(1.5)` gives

\[
\frac{b_p}{5}\equiv1\pmod{p^3}.
\]

At the next recurrence step,

\[
(p+1)^3b_{p+1}=P(p)b_p-p^3b_{p-1}.
\]

The last term disappears modulo `p^3`, and

\[
\frac{P(p)}{(1+p)^3}
\equiv5+12p\pmod{p^3};
\tag{3.3}
\]

the coefficient of `p^2` cancels exactly:

\[
(5+27p+51p^2)(1-3p+6p^2)
=5+12p+O(p^3).
\]

Therefore

\[
\frac{b_{p+1}}5\equiv5+12p=F_1(p)\pmod{p^3}.
\]

For the rational companion, multiply `(1.2)` at `m=p` by `p^3`:

\[
(p^3a_p)b_{p-1}-p^3a_{p-1}b_p=6.
\]

The second term is divisible by `p^3` in `Z_p`, while `b_{p-1}\equiv1 (mod p^3)`. Hence

\[
p^3a_p\equiv6\pmod{p^3}.
\tag{3.4}
\]

Multiplying the recurrence at `m=p` by `p^3` similarly yields

\[
\frac{p^3a_{p+1}}6\equiv5+12p\pmod{p^3}.
\tag{3.5}
\]

Now both normalized shifted sequences

\[
\frac{b_{p+r}}5,
\qquad
\frac{p^3a_{p+r}}6
\]

satisfy the exact recurrence with coefficients shifted by `p`. Equation `(2.5)` says `F_r(p)` satisfies that recurrence modulo `p^3`. For `0<=r<=p-2`, every forward coefficient `(p+r+1)^3` is a `p`-unit, so the recurrence has unique propagation modulo `p^3` from its first two values. Equations `(3.1)--(3.2)` follow. ∎

A direct exact-rational regression, independent of the derivation, is given in the verifier below. For `p=7,r=3`,

\[
F_3(7)=1445+7\cdot4438+49\cdot2160=138351,
\]

while exact recurrence gives

\[
a_{10}=\frac{43786938951280269198311}{2667168000},
\qquad
b_{10}=13657436403073.
\]

After multiplying the rational companion by `7^3`, its reduced denominator is a `7`-unit, and exact modular reduction verifies

\[
7^3a_{10}-6F_3(7)\equiv0\pmod{7^3},
\]

\[
b_{10}-5F_3(7)\equiv0\pmod{7^3}.
\]

---

## 4. Exact target ladder through valuation two

Assume now

\[
p\mid b_r.
\]

The cases `r=0,1` are impossible for `p>=7`, so `2<=r<=p-2`. Write

\[
b_r=pc,
\qquad c\in\mathbf Z.
\]

From `(3.1)--(3.2)`,

\[
B_{p,r}\equiv5p(c+G_r+pH_r)\pmod{p^3},
\tag{4.1}
\]

\[
A_{p,r}\equiv6p(c+G_r+pH_r)\pmod{p^3}.
\tag{4.2}
\]

Define

\[
\boxed{\theta_1(p,r)=c+G_r\pmod p.}
\tag{4.3}
\]

Then

\[
\theta_1\ne0
\iff
v_p(B_{p,r})=v_p(A_{p,r})=1.
\tag{4.4}
\]

If `theta_1=0`, then `c+G_r` lies in `p Z_p`; put

\[
\eta_{p,r}=\frac{c+G_r}{p}\in\mathbf Z_p
\]

and define

\[
\boxed{\theta_2(p,r)=\eta_{p,r}+H_r\pmod p.}
\tag{4.5}
\]

Then

\[
\theta_2\ne0
\iff
v_p(B_{p,r})=v_p(A_{p,r})=2.
\tag{4.6}
\]

Only

\[
\boxed{\theta_1=\theta_2=0}
\tag{4.7}
\]

can produce `v_p(B_{p,r})>=3` and therefore enter the seam `(E)`.

The first scalar also has a fixed characteristic-zero packaging. For a positive integer `m`, define its arithmetic derivative

\[
\partial_a(m)=\sum_{\ell\mid m}v_\ell(m)\frac m\ell.
\]

For every prime `p|m`,

\[
\partial_a(m)\equiv m/p\pmod p.
\tag{4.8}
\]

Let `q_r` clear the reduced denominator of `G_r`. For every eligible `p>r+1`, `q_r` is a `p`-unit. Thus

\[
R_r=q_r\bigl(\partial_a(b_r)+G_r\bigr)\in\mathbf Z
\]

satisfies, at every strict target,

\[
\boxed{
\theta_1(p,r)=0\iff p\mid R_r.
}
\tag{4.9}
\]

This does not solve the global height problem, but it shows the first exceptional layer is a genuine fixed-in-`r` arithmetic condition, not a moving-modulus tautology.

---

## 5. Exact Wronskian reduction at valuation three

The following theorem is the clean local endpoint of the present audit.

### Theorem 5.1
At a strict target define

\[
A=A_{p,r}=p^3a_{p+r},\quad B=B_{p,r}=b_{p+r},
\]

\[
A_-=p^3a_{p+r-1},\quad B_-=b_{p+r-1}.
\]

Then `A_-` and `B_-` are `p`-units. If `t=v_p(B)`, then

\[
\boxed{
v_p(A)=
\begin{cases}
t,&t<3,\\
3+v_p(\Omega_{p,r}),&t=3,\\
3,&t>3,
\end{cases}}
\tag{5.1}
\]

where, in the middle case,

\[
\boxed{
\Omega_{p,r}
=A_-\frac{B}{p^3}+\frac6{(p+r)^3}.
}
\tag{5.2}
\]

### Proof
Lucas/first-block propagation gives

\[
B_-\equiv5b_{r-1}\pmod p,
\qquad
A_-\equiv6b_{r-1}\pmod p.
\]

The lower Casoratian at the target says

\[
a_rb_{r-1}\equiv\frac6{r^3}\pmod p,
\]

so `b_{r-1}` and `a_r` are both `p`-units. Hence `A_-,B_-` are units.

Multiplying `(1.2)` at `n=p+r` by `p^3` gives the exact `Z_p` identity

\[
\boxed{
AB_- - A_-B=\frac{6p^3}{(p+r)^3}.
}
\tag{5.3}
\]

If `t<3`, the two terms on the right after solving for `AB_-` have distinct valuations `t` and `3`, so no cancellation is possible and `v_p(A)=t`. If `t>3`, the fixed Wronskian term has the smaller valuation `3`, so `v_p(A)=3`. If `t=3`, divide `(5.3)` by `p^3`; since `B_-` is a unit, `(5.1)--(5.2)` follow. ∎

### Corollary 5.2 — exact seam scalar

\[
\boxed{
(E)\iff
v_p(B)=3
\quad\text{and}\quad
\Omega_{p,r}\in p^3\mathbf Z_p.
}
\tag{5.4}
\]

Its first digit is

\[
\Omega_{p,r}\equiv
6b_{r-1}c_3+\frac6{r^3}\pmod p,
\qquad c_3=B/p^3.
\tag{5.5}
\]

Therefore a necessary first depth-three condition is

\[
\boxed{
c_3\equiv-\frac1{r^3b_{r-1}}
=-\frac{a_r}{6}\pmod p.}
\tag{5.6}
\]

If `(5.6)` holds, two further digits

\[
\frac{\Omega_{p,r}}p\pmod p,
\qquad
\frac{\Omega_{p,r}}{p^2}\pmod p
\]

must both vanish. These are the first genuinely unresolved scalars after the universal quadratic jet. The recurrence and Casoratian used so far do not force them.

---

## 6. Why the order-`p^3` companion difference does not close the seam

The shared quadratic jet implies

\[
D_{p,r}:=5A_{p,r}-6B_{p,r}\in p^3\mathbf Z_p.
\]

Put

\[
\delta_{p,r}=D_{p,r}/p^3\pmod p.
\]

Because `A` and `B` satisfy the same shifted recurrence, `delta_{p,r}` satisfies the ordinary Apéry recurrence in `r` modulo `p`. Let

\[
\kappa_p=\delta_{p,0}.
\]

At the singular first step, writing

\[
D_{p,-1}=5p^3a_{p-1}-6b_{p-1}\equiv-6\pmod p,
\]

gives

\[
\delta_{p,1}=5\kappa_p+6.
\]

Hence uniqueness of the ordinary second-order recurrence yields

\[
\boxed{
\delta_{p,r}\equiv\kappa_pb_r+a_r\pmod p.
}
\tag{6.1}
\]

At a target,

\[
\boxed{
\delta_{p,r}\equiv a_r\ne0\pmod p.
}
\tag{6.2}
\]

Thus `v_p(5A-6B)=3` exactly at every target. This is genuine transversality between the two shifted coordinates, but it does **not** rule out `(E)`. If `(E)` held, then

\[
5A/p^3\equiv0\pmod p,
\qquad B/p^3=c_3,
\]

and `(6.2)` would give exactly

\[
c_3\equiv-a_r/6\pmod p,
\]

which is the same condition `(5.6)` already obtained from the Wronskian. It is not an independent second equation.

This is why merely repeating the first-order block law, the Casoratian, or the apparent `p^3` difference does not settle the depth-three seam.

---

## 7. Exhaustive exact scan through `p<=10000`

I scanned **every** prime

\[
7\le p\le10000
\]

and every strict digit

\[
0\le r\le p-2.
\]

There are `1226` such primes. The scan found `1172` target pairs `p|b_r`.

The scan itself uses only exact modular recurrence arithmetic:

- `b_r` is propagated modulo `p^2`;
- `G_r` is propagated modulo `p` using `(2.1)`;
- for every target, `theta_1=(b_r/p)+G_r (mod p)` is evaluated exactly.

Exactly five target pairs satisfy `theta_1=0`:

| `p` | `r` | `n=p+r` | `b_r/p mod p` | `G_r mod p` |
|---:|---:|---:|---:|---:|
| 41 | 30 | 71 | 35 | 6 |
| 97 | 25 | 122 | 94 | 3 |
| 151 | 14 | 165 | 112 | 39 |
| 1453 | 1180 | 2633 | 858 | 595 |
| 6781 | 3974 | 10755 | 2986 | 3795 |

For each of these five, I then lifted the actual shifted values, not just the first jet:

- `B=b_{p+r}` was generated by the exact integer Apéry recurrence;
- `A=p^3a_{p+r}` was propagated modulo `p^7` from the rational initial data, inverting only denominators known to be `p`-units; the singular step is handled without division by `p` by
  \[
  p^3a_p=P(p-1)a_{p-1}-(p-1)^3a_{p-2}.
  \]

The exact results are:

| `p` | `r` | `B/(5p^2) mod p` | `v_p(B)` | `v_p(A)` |
|---:|---:|---:|---:|---:|
| 41 | 30 | 7 | 2 | 2 |
| 97 | 25 | 70 | 2 | 2 |
| 151 | 14 | 18 | 2 | 2 |
| 1453 | 1180 | 741 | 2 | 2 |
| 6781 | 3974 | 1508 | 2 | 2 |

The third column is exactly `theta_2`; every value is nonzero. Hence no target with `p<=10000` even reaches `v_p(B)>=3`, and therefore `(E)` does not occur in the entire scanned range.

This finite result is not promoted to an all-prime theorem. Its theorem-level value is different: it decisively refutes the stronger universal guess `theta_1!=0` (already `(41,30)` is an exact counterexample), while the five independent cancellations show that first-jet vanishing is a real arithmetic phenomenon rather than one accidental small-prime anomaly.

---

## 8. Standalone exact verifier

The following Python 3 program uses only integer arithmetic and `fractions.Fraction`. It reproduces the full `p<=10000` first-jet scan, checks the five exceptions, verifies their exact valuation-two lift, and includes an exact rational regression of the quadratic block theorem.

```python
from fractions import Fraction
from math import gcd

PMAX = 10_000
EXPECTED_HITS = [
    (41, 30, 35, 6),
    (97, 25, 94, 3),
    (151, 14, 112, 39),
    (1453, 1180, 858, 595),
    (6781, 3974, 2986, 3795),
]
EXPECTED_THETA2 = {
    (41, 30): 7,
    (97, 25): 70,
    (151, 14): 18,
    (1453, 1180): 741,
    (6781, 3974): 1508,
}


def P(n):
    return 34*n**3 + 51*n**2 + 27*n + 5


def P1(n):
    return 102*n*n + 102*n + 27


def primes_upto(N):
    s = bytearray(b"\x01") * (N + 1)
    s[:2] = b"\x00\x00"
    for q in range(2, int(N**0.5) + 1):
        if s[q]:
            s[q*q:N+1:q] = b"\x00" * (((N-q*q)//q) + 1)
    return [q for q in range(7, N + 1) if s[q]]


def inverse_table_below_p(p, mod):
    """inv[k] = k^{-1} mod mod for 1 <= k < p, mod=p^e."""
    inv = [0] * p
    inv[1] = 1
    for k in range(2, p):
        q, r = divmod(mod, k)
        # Since k<p and p is prime, r != 0 for k>1.
        inv[k] = (-q * inv[r]) % mod
        assert (k * inv[k]) % mod == 1
    return inv


def first_jet_scan(PMAX=PMAX):
    targets = 0
    hits = []
    ps = primes_upto(PMAX)

    for p in ps:
        mod = p*p
        inv = inverse_table_below_p(p, mod)

        # b_0,b_1 modulo p^2; G_0,G_1 modulo p.
        b0, b1 = 1, 5
        g0, g1 = 0, 12 % p

        # s=1,...,p-3 computes rr=s+1=2,...,p-2.
        # rr=0,1 can never be targets for p>=7.
        for s in range(1, p-2):
            den = s + 1
            inv3 = pow(inv[den], 3, mod)
            b2 = ((P(s) % mod) * b1 - (s**3 % mod) * b0) * inv3 % mod

            rhs = (
                (P1(s) % p) * (b1 % p)
                - 3 * ((s+1)**2 % p) * (b2 % p)
                - 3 * (s*s % p) * (b0 % p)
            ) % p
            inv3p = pow(inv[den] % p, 3, p)
            g2 = (
                (P(s) % p) * g1
                - (s**3 % p) * g0
                + rhs
            ) * inv3p % p

            rr = s + 1
            if b2 % p == 0:
                targets += 1
                c = (b2 // p) % p  # exact because 0 <= b2 < p^2
                theta1 = (c + g2) % p
                if theta1 == 0:
                    hits.append((p, rr, c, g2))

            b0, b1 = b1, b2
            g0, g1 = g1, g2

    assert len(ps) == 1226
    assert targets == 1172
    assert hits == EXPECTED_HITS
    return ps, targets, hits


def exact_b_sequence(N):
    b = [0] * (N + 1)
    b[0] = 1
    if N:
        b[1] = 5
    for m in range(1, N):
        num = P(m) * b[m] - m**3 * b[m-1]
        den = (m+1)**3
        assert num % den == 0
        b[m+1] = num // den
    return b


def vp_int(x, p):
    assert x != 0
    e = 0
    while x % p == 0:
        x //= p
        e += 1
    return e


def vp_residue(x, p, cap):
    """Exact valuation if it is < cap, for x represented modulo p^cap."""
    x %= p**cap
    if x == 0:
        return cap
    e = 0
    while x % p == 0:
        x //= p
        e += 1
    return e


def companion_A_mod(p, r, cap=7):
    """Return p^3 a_{p+r} mod p^cap, with no illegal p-division."""
    mod = p**cap

    # First propagate rational a_m in Z_p up to m=p-1.
    # Every denominator (m+1)^3 here is a p-unit.
    am2, am1 = 0, 6  # a_0,a_1 modulo p^cap
    for m in range(1, p-1):  # m=1,...,p-2
        den = (m+1)**3
        assert gcd(den, p) == 1
        nxt = (P(m)*am1 - m**3*am2) * pow(den, -1, mod) % mod
        am2, am1 = am1, nxt
    # am2=a_{p-2}, am1=a_{p-1} modulo p^cap.

    # Singular step m=p-1: compute p^3 a_p without dividing by p^3.
    Aminus = (p**3 * am1) % mod                  # p^3 a_{p-1}
    A0 = (P(p-1)*am1 - (p-1)**3*am2) % mod      # p^3 a_p

    if r == 0:
        return A0

    # After the anchor, all forward denominators are again p-units for r<=p-2.
    aprev, acur = Aminus, A0
    for k in range(r):
        den = (p+k+1)**3
        assert gcd(den, p) == 1
        anext = (
            P(p+k)*acur - (p+k)**3*aprev
        ) * pow(den, -1, mod) % mod
        aprev, acur = acur, anext
    return acur


def verify_five_lifts():
    maxn = max(p+r for p, r, _, _ in EXPECTED_HITS)
    b = exact_b_sequence(maxn)

    for p, r, c, g in EXPECTED_HITS:
        assert b[r] % p == 0
        assert (b[r] // p) % p == c

        n = p + r
        B = b[n]
        vB = vp_int(B, p)
        assert vB == 2

        theta2 = (B // (p*p)) * pow(5, -1, p) % p
        assert theta2 == EXPECTED_THETA2[(p, r)]
        assert theta2 != 0

        A = companion_A_mod(p, r, cap=7)
        vA = vp_residue(A, p, cap=7)
        assert vA == 2

        # In particular, none reaches the entry condition v_p(B)>=3.
        assert not (vB == 3 and vA >= 6)


def exact_a_sequence(N):
    a = [Fraction(0) for _ in range(N + 1)]
    a[0] = Fraction(0)
    if N:
        a[1] = Fraction(6)
    for m in range(1, N):
        a[m+1] = (
            P(m)*a[m] - m**3*a[m-1]
        ) / (m+1)**3
    return a


def rat_mod(q, mod):
    q = Fraction(q)
    assert gcd(q.denominator, mod) == 1
    return q.numerator * pow(q.denominator, -1, mod) % mod


def exact_quadratic_regression():
    # b_3=1445, G_3=4438, H_3=2160.
    p, r = 7, 3
    F3 = 1445 + p*4438 + p*p*2160
    assert F3 == 138351

    a = exact_a_sequence(10)
    assert a[10] == Fraction(43786938951280269198311, 2667168000)
    b10 = exact_b_sequence(10)[10]
    assert b10 == 13657436403073

    A10 = p**3 * a[10]
    # After multiplying by p^3, the reduced denominator is a p-unit.
    assert gcd(A10.denominator, p) == 1
    assert rat_mod(A10 - 6*F3, p**3) == 0
    assert (b10 - 5*F3) % (p**3) == 0


if __name__ == "__main__":
    ps, targets, hits = first_jet_scan()
    print("PRIMES", len(ps))
    print("STRICT_TARGETS", targets)
    print("THETA1_ZERO", hits)
    verify_five_lifts()
    exact_quadratic_regression()
    print("ALL_EXACT_CHECKS_OK")
```

Expected terminal summary is

```text
PRIMES 1226
STRICT_TARGETS 1172
THETA1_ZERO [(41, 30, 35, 6), (97, 25, 94, 3), (151, 14, 112, 39), (1453, 1180, 858, 595), (6781, 3974, 2986, 3795)]
ALL_EXACT_CHECKS_OK
```

---

## 9. Final scoped theorem

Combining the exact derivations gives the following fully quantified statement.

### Theorem 9.1
Let `p>=7`, let `n=p+r` with `p>(n+1)/2`, and assume `p|b_r`. Then:

1. `A_{p,r}=p^3a_{p+r}` is `p`-integral, and `(Q1)--(Q2)` hold.
2. If `theta_1!=0`, then
   \[
   v_p(A_{p,r})=v_p(B_{p,r})=1.
   \]
3. If `theta_1=0` but `theta_2!=0`, then
   \[
   v_p(A_{p,r})=v_p(B_{p,r})=2.
   \]
4. If `v_p(B_{p,r})<3`, then in fact
   \[
   v_p(A_{p,r})=v_p(B_{p,r}).
   \]
5. If `v_p(B_{p,r})=3`, then `(E)` is equivalent to the explicit depth-three condition
   \[
   \Omega_{p,r}\in p^3\mathbf Z_p
   \]
   with `Omega` defined by `(5.2)`; its first required digit is `(5.6)`.
6. The apparently new order-`p^3` difference `5A-6B` supplies no independent equation at a target: its first divided residue is `a_r`, and under `(E)` it reduces exactly to `(5.6)`.

Thus the seam is not closed by first-order Lucas, by the Wronskian alone, or by the first nontrivial companion difference. The first two lower digits are explicit; the remaining obstruction is the three-digit scalar `Omega`.

### Exact finite corollary
For every prime `7<=p<=10000` and every strict digit `0<=r<=p-2` with `p|b_r`, `(E)` is false. The exhaustive scan contains `1172` target pairs; exactly five reach valuation two, and none reaches valuation three.

This finite corollary is a certificate, not an asymptotic theorem. An all-prime proof now needs a theorem excluding simultaneous vanishing of `theta_1,theta_2` and the three Wronskian digits in `(O3)`, or a characteristic-zero carrier for those divided jets with a usable height bound.
