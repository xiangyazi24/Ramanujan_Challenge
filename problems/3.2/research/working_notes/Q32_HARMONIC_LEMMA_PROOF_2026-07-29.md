# The weight-three harmonic lemma: H^{(2)}_{p-1} == 2p*Xi_p (mod p^3)
# Owner: Fable (zinan:7), 2026-07-29. Status: PROVED (two standard
# inputs), every step numerically verified 75/75 primes 7<=p<400.

## Statement

For every prime p >= 7, with Xi_p = B_{2p-4}/(2p-4) - 2B_{p-3}/(p-3):

```text
H^{(2)}_{p-1} := sum_{s=1}^{p-1} 1/s^2 == 2p Xi_p   (mod p^3).
```

Combined with the rank-one endpoint theorem
(Q32_ENDPOINT_BERNOULLI_RANK_ONE_2026-07-29.md, (0.2)+(4.10)) this
proves Z.-H. Sun's Conjecture 2.4 (arXiv:2409.06544v2) at r=1, for
all m, in both branches — statement CONFIRMED against the source PDF.

## Proof

### Step 1: Fermat-quotient expansion into three power sums

For each unit s, write s^{p-1} = 1 + p q_s. Then, exactly,
1/s^2 = s^{2p-4} (s^{p-1})^{-2}, and mod p^3

```text
(1+p q_s)^{-2} == 1 - 2p q_s + 3p^2 q_s^2,
```

so, absorbing p q_s = s^{p-1}-1,

```text
H2 == sum s^{2p-4} - 2[sum s^{3p-5} - sum s^{2p-4}]
      + 3[sum s^{4p-6} - 2 sum s^{3p-5} + sum s^{2p-4}]   (mod p^3)
   == 6 S(2p-4) - 8 S(3p-5) + 3 S(4p-6),
S(k) := sum_{s=1}^{p-1} s^k.
```

### Step 2: Faulhaber. For even k with p-3 <= k <= 4p-6 and p >= 7,

```text
S(k) == p B_k   (mod p^3).
```

Faulhaber gives S(k) = (1/(k+1)) sum_j C(k+1,j) B_j p^{k+1-j}. The
j=k term is p B_k; the j=k-1 term vanishes (B_{odd>=3}=0, k even);
the j=k-2 term is (k(k-1)/6) B_{k-2} p^3 with B_{k-2} p-integral by
von Staudt-Clausen ((p-1) | k-2 would force (p-1)|4, i.e. p=5); all
lower j carry p^{>=4} against at most one p in a Bernoulli
denominator. Hence

```text
H2 == p [ 6B_{2p-4} - 8B_{3p-5} + 3B_{4p-6} ]   (mod p^3).   (*)
```

[Verified 75/75.]

### Step 3: strong Kummer congruence along the branch

The exponents p-3, 2p-4, 3p-5, 4p-6 form the arithmetic progression
m_k = (p-3) + k(p-1), k=0..3, all == -2 mod (p-1), never divisible
by p-1. Write beta_k = B_{m_k}/m_k. The strong (second-order) Kummer
congruence — the a=2 case of

```text
sum_{j=0}^{a} (-1)^j C(a,j) (1-p^{m-1+j(p-1)}) B_{m+j(p-1)}/(m+j(p-1))
   == 0 (mod p^a)
```

(equivalently: Mahler coefficients of the branch of the Kubota-
Leopoldt L-function lie in p^j Z_p; see Washington, Introduction to
Cyclotomic Fields, §5.1-5.2, or Z.-H. Sun, Discrete Appl. Math. 105
(2000), Cor. 2.x — exact numbering to be pinned when the tex is
written) — gives, since the Euler factors p^{m-1} here are O(p^{p-4})
= 0 mod p^3:

```text
beta_2 == 2 beta_1 - beta_0,     beta_3 == 3 beta_1 - 2 beta_0
                                              (mod p^2).
```

[Verified 75/75.]

### Step 4: recombination

Substitute B_{m_k} = m_k beta_k into (*):

```text
coefficient of beta_1:
  6(2p-4) - 8(3p-5)*2 + 3(4p-6)*3 = 12p-24 -48p+80 +36p-54 = 2;
coefficient of beta_0:
  8(3p-5) - 6(4p-6)*... : - 8(3p-5)(-1) - 3(4p-6)(2) ->
  +8(3p-5) - 6(4p-6) = 24p-40-24p+36 = -4.
```

Hence H2 == p[2 beta_1 - 4 beta_0] = 2p[B_{2p-4}/(2p-4)
- 2B_{p-3}/(p-3)] = 2p Xi_p (mod p^3).  QED

### Sanity reductions

Mod p^2 this collapses (Kummer beta_1 == beta_0 mod p) to
H2 == -2p beta_0 = -2p B_{p-3}/(p-3) == (2/3) p B_{p-3}, the
classical evaluation — consistent.

## Consequences (all now unconditional at r=1)

With Delta_p = b_{p-1}-1 == p^2 H2 (mod p^5) and the rank-one laws:

```text
b_{mp}  - b_m     == 2 C_m       p^3 Xi_p (mod p^5)   [Sun Conj 2.4, r=1]
b_{mp-1}- b_{m-1} == 2 C_{-m}    p^3 Xi_p (mod p^5)   [reflected 2.4 —
                                            not even conjectured by Sun]
```

C_m = m^3(b_{m-1}-17b_m)/12, C_{-m} = m^3(17b_{m-1}-b_m)/12.

## Verification script

../scripts/q32_harmonic_lemma_verify.py (steps A/B inline
in the r24-r25 transcript; end-to-end 106/106 primes to 600).
