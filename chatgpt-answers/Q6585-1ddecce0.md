ANSWER Q6585 1ddecce0.

# qATR audit: exact third moment, COINC3, and the best current bound

Repository state audited: head e2f2e26beaa0a084b171dc47f695c68c474857f1, with the mounted orbit-energy proof in problems/3.2/FABLE_SECTION_orbit_energy.tex, the corrected nonwrapping triple identity in Q6511.md, and the proved triple-moment interpolation lemma in Q6533-b4b60c0c.md.

## Executive verdict

1. Claims D1-D4 are correct on the stated regular nonwrapping windows. D1 has a short exact proof from the Casoratian and the recurrence.

1. D5 contains a decisive notation error. The sigma-ladder quantity is

```plain text
Q_p(H) = sum_r binom(d_H(r),2)
       = sum_{1 <= h < k <= H} J_p(h,k),
```

whereas D4 controls the first moment

```plain text
S_p(H) = sum_r d_H(r) = sum_{h <= H} C_h.
```

They are not interchangeable. At H=p^(2/3), COINC gives S_p(H)<<p^(5/6), not Q_p(H)<<p. Thus the asserted sigma-3/2 inference is invalid. The numerical conclusion F_p<<p^(4/3) is nevertheless correct by the direct block-energy inequality, without the sigma ladder.

1. D6 understates the generic Weil error. After cross multiplication, the difference curve has bidegree at most (3h,3h') in P^1 x P^1. Its generic arithmetic genus is of order h h', so a generic per-pair error is O(h h' sqrt(p)), not O((h+h')sqrt(p)). Summing gives O(H^4 sqrt(p)), not O(H^3 sqrt(p)). This only reaches the COINC error scale for H<=p^(1/6), where the elementary bound is much better. Component control is also open.

1. With the ordered nonwrapping convention in the question,

```plain text
A_p := sum_{h != k} J_p(h,k)
     = (1/3) sum_v (m_p(v))_3.
```

The factor 1/3 is exact. Each increasing same-fiber triple gives two ordered positive-gap pairs, while (m)_3=6 binom(m,3).

1. The direct trilinear analogue of COINC is a collision estimate for the two-dimensional value vector

```plain text
(Delta_{r,h}, Delta_{r,k}).
```

Its natural Poisson-scale form gives only

```plain text
A_p(H) << H^2/p + H sqrt(p).
```

At full scale this is O(p^(3/2)), not ATR. To get ATR from Cauchy-Parseval alone would require a statistically much stronger, super-Poisson variance saving by a factor p, or else direct signed cancellation on the frequency plane.

1. From the unconditional inputs listed in the question, the best bound I can prove today is

```plain text
sum_{h != k} J_p(h,k) << p^2 log p.
```

This improves the previous p^(13/6) consequence of F_p<<p^(3/2) and max_v m_p(v)<<p^(2/3), but it is not o(p^(3/2)); it is not even a fixed power below p^2.

1. The gap-degree structure gives

```plain text
max_v m_p(v) << p^(2/3),
```

not p^(1/2). The precise optimization is p/H + H^2, whose minimum is at H about p^(1/3).

1. If pair COINC is known up to H=p^beta, then a rigorous two-short-gap/one-long amplification gives

```plain text
T_3 := sum_v (m_p(v))_3
 << p^(3-2 beta) + p^(9/4-3 beta/4) + p^(3/2).
```

Thus a power saving below p^2 begins exactly when beta>1/2. At beta=2/3 one gets T_3<<p^(7/4). Even beta=1 gives only T_3<<p^(3/2) by this mechanism, not o(p^(3/2)) and not ATR.

A further threshold correction is important: by the proved T3 interpolation lemma, any estimate

```plain text
T_3 << p^(2-delta)
```

already yields

```plain text
F_p << p^(3/2-delta/2).
```

Therefore p^(3/2) is not the first publishable threshold for the triple sum. Any fixed saving below p^2 is already a new energy theorem. A bound below p^(3/2) would be much stronger.

# 1. Audit of D1: exact determinant formula

Put

```plain text
v_n = (b_n,c_n),
Delta_{r,h} = det(v_r,v_{r+h}).
```

Fix r and write z_h=Delta_{r,h}. Taking the determinant of the recurrence for v_{r+h} against v_r gives

```plain text
(r+h+1)^3 z_{h+1}
 = P(r+h) z_h - (r+h)^3 z_{h-1}.
```

The initial values are

```plain text
z_0 = 0,
z_1 = b_r c_{r+1} - b_{r+1} c_r
    = -W_{r+1}
    = 1/(r+1)^3.
```

Define

```plain text
D_h(r) = product_{j=1}^h (r+j)^3,
K_h(r) = D_h(r) z_h.
```

Multiplying the recurrence by D_{h+1} gives

```plain text
K_{h+1}(r)
 = P(r+h) K_h(r) - (r+h)^6 K_{h-1}(r).
```

Also

```plain text
K_1=1,
K_2=P(r+1).
```

These are exactly the defining initial values and recurrence for N_h. Hence, on every regular nonwrapping window,

```plain text
Delta_{r,h}
 = N_h(r) / product_{j=1}^h (r+j)^3.
```

So D1 is proved, including its sign. Every denominator is nonzero because 1<=r and r+h<=p-2.

# 2. Audit of D2-D4: the pair-level second moment

Let

```plain text
I_H = {(r,h): 1<=h<=H, 1<=r<=M-h},
L_1(H) = #I_H = H M - H(H+1)/2,
S_p(H) = sum_{h<=H} C_h.
```

For

```plain text
B_t(H) = sum_{(r,h) in I_H} e_p(t Delta_{r,h}),
```

additive orthogonality gives the exact identity

```plain text
S_p(H)
 = L_1(H)/p + (1/p) sum_{t != 0} B_t(H).
```

Define

```plain text
N_coinc(H)
 = #{(a,a') in I_H^2: Delta_a=Delta_a'}.
```

Full Parseval gives

```plain text
sum_{t in F_p} |B_t(H)|^2 = p N_coinc(H).
```

Since B_0(H)=L_1(H), this is exactly

```plain text
sum_{t != 0} |B_t(H)|^2
 = p N_coinc(H) - L_1(H)^2.
```

Therefore the hypothesis

```plain text
N_coinc(H) <= L_1(H)^2/p + K L_1(H)
```

implies

```plain text
S_p(H)
 <= L_1(H)/p + sqrt(((p-1)/p) K L_1(H))
 << H + sqrt(K p H).
```

Thus D2-D4 are correct.

# 3. D5 correction: first moment versus triple moment

The mounted sigma ladder uses

```plain text
d_H(r) = #{1<=h<=H: r+h<=M and pi(r)=pi(r+h)},
Q_p(H) = sum_r binom(d_H(r),2).
```

There are two exact identities:

```plain text
S_p(H) = sum_r d_H(r) = sum_{h<=H} C_h,
Q_p(H) = sum_r binom(d_H(r),2)
       = sum_{1<=h<k<=H} J_p(h,k).
```

D4 controls S_p(H). The sigma hypothesis controls Q_p(H). In particular, at H=p^(2/3), sigma 3/2 would require

```plain text
Q_p(H) << H^(3-3/2) = H^(3/2) = p.
```

But D4 gives only

```plain text
S_p(H) << p^(5/6).
```

There is no implication from this to Q_p(H)<<p. Even using the uniform multiplicity bound below gives only

```plain text
Q_p(H)
 <= (max_r d_H(r)) S_p(H)/2
 <= (max_v m_p(v)) S_p(H)/2
 << p^(3/2).
```

So the sigma-ladder sentence in D5 must be deleted.

The desired energy conclusion survives by a different exact inequality. Partition I_p into consecutive blocks of length at most H. Every within-block collision has gap at most H-1, so

```plain text
F_p
 <= ceil(M/H) (M + 2 S_p(H-1)).
```

Under D4,

```plain text
F_p
 << p^2/H + p^(3/2)/sqrt(H) + p.
```

For H=p^beta, this is

```plain text
F_p << p^(2-beta) + p^(3/2-beta/2) + p.
```

At beta=2/3, the first term dominates and gives

```plain text
F_p << p^(4/3).
```

Thus the 4/3 prize is real, but it is a first-moment block consequence, not a sigma-3/2 consequence.

# 4. D6 correction: the curve error is generically H^4 sqrt(p)

Write

```plain text
Delta_{r,h} = N_h(r)/D_h(r),
D_h(r) = product_{j=1}^h (r+j)^3.
```

The equation

```plain text
Delta_{r,h}=Delta_{r',h'}
```

cross-multiplies to

```plain text
N_h(r) D_{h'}(r') - D_h(r) N_{h'}(r') = 0.
```

Its bidegree in P^1 x P^1 is at most

```plain text
(3h,3h').
```

For a smooth geometrically irreducible curve of bidegree (a,b), the arithmetic genus is (a-1)(b-1). Thus the generic genus here is

```plain text
O(h h'),
```

and the generic Hasse-Weil error is

```plain text
O(h h' sqrt(p)).
```

Summing over h,h'<=H gives

```plain text
sqrt(p) sum_{h,h'<=H} h h'
 << H^4 sqrt(p).
```

The COINC allowance is O(#I_H)=O(pH). Per-pair Weil therefore reaches that scale only if

```plain text
H^4 sqrt(p) <= pH,
```

or

```plain text
H <= p^(1/6).
```

But COINC first improves the elementary S_p(H)<<H^2 bound only when H>p^(1/3). There is no overlap.

There is also a component problem before Weil can be summed:

- for h=h', the fiber product has the diagonal component r=r';

- a bounded number of components would require a uniform monodromy or double-transitivity statement for each rational map Delta_h;

- for h!=h', common compositional factors or functional relations can create further components;

- adjacent coprimality of N_h,N_{h+1} controls only one zero fiber and does not prove irreducibility of these fiber products.

The correct family target is therefore stronger than average cancellation of already-controlled Weil errors: it first needs uniform component classification, then cancellation of genus-sized errors over a discrete growing family.

# 5. Exact third-moment identities

Let

```plain text
T_3 = sum_v (m_p(v))_3,
A_p = sum_{1<=h,k<=M-1, h!=k} J_p(h,k),
S_p = sum_{h=1}^{M-1} C_h.
```

Every unordered same-fiber pair has one positive gap, hence

```plain text
F_p = M + 2 S_p.
```

Every increasing same-fiber triple n_1<n_2<n_3 contributes exactly twice to A_p, through the two ordered gap pairs

```plain text
(h,k)=(n_2-n_1,n_3-n_1)
and
(h,k)=(n_3-n_1,n_2-n_1).
```

Therefore

```plain text
A_p = 2 sum_v binom(m_p(v),3)
    = T_3/3.
```

Equivalently,

```plain text
sum_{1<=h<k<=M-1} J_p(h,k)
 = sum_v binom(m_p(v),3)
 = T_3/6.
```

The raw cubic moment has the exact forms

```plain text
sum_v m_p(v)^3
 = M + 6 S_p + 3 A_p,
```

and

```plain text
sum_v m_p(v)^3
 = T_3 + 3 F_p - 2M.
```

These are the clean identities to put in the paper. The explanation sometimes written as "three choices of base point" is not correct for positive nonwrapping gaps. Only the smallest index can be the base. The factor comes from the two orders of the two positive gaps.

One useful lower relation is

```plain text
sum_v m_p(v)^3 >= F_p^2/M,
```

by Cauchy. Hence

```plain text
T_3 >= F_p^2/M - 3F_p + 2M.
```

This is another way to see why ATR is genuinely stronger than E1.

# 6. Exact trilinear Fourier identity

For 1<=H<=M-1, define the ordered nonwrapping triple domain

```plain text
T_H = {(r,h,k):
       1<=h,k<=H,
       h!=k,
       1<=r<=M-max(h,k)}.
```

Its exact cardinality is

```plain text
L_3(H)
 = 2 sum_{q=2}^H (q-1)(M-q)
 = H(H-1)(M - 2(H+1)/3).
```

At full scale,

```plain text
L_3(M-1) = M(M-1)(M-2)/3.
```

Define

```plain text
B_{t1,t2}^{(3)}(H)
 = sum_{(r,h,k) in T_H}
   e_p(t1 Delta_{r,h} + t2 Delta_{r,k}).
```

Then two-dimensional additive orthogonality gives

```plain text
A_p(H) := sum_{h!=k<=H} J_p(h,k)
 = L_3(H)/p^2
   + (1/p^2)
     sum_{(t1,t2)!=(0,0)} B_{t1,t2}^{(3)}(H).
```

This is the requested exact Fourier formula.

A small but important axis warning: the terms with t2=0 are not the old unweighted pair sums. In fact

```plain text
B_{t,0}^{(3)}(H)
 = sum_{(r,h) in I_H}
   w_H(r,h) e_p(t Delta_{r,h}),
```

where

```plain text
w_H(r,h) = min(H,M-r)-1.
```

Thus pair COINC does not automatically control the frequency axes. A correct trilinear hypothesis must include these weighted marginals or work with the full vector collision count below.

# 7. The correct COINC3 statement and what it yields

Define

```plain text
N_coinc^{(3)}(H)
 = #{(a,a') in T_H^2:
      (Delta_{r,h},Delta_{r,k})
      =(Delta_{r',h'},Delta_{r',k'})}.
```

Two-dimensional Parseval gives the exact identity

```plain text
sum_{(t1,t2)!=(0,0)}
 |B_{t1,t2}^{(3)}(H)|^2
 = p^2 N_coinc^{(3)}(H) - L_3(H)^2.
```

The direct Poisson-scale analogue of pair COINC is

```plain text
[COINC3-0]
N_coinc^{(3)}(H)
 <= L_3(H)^2/p^2 + K L_3(H).
```

Cauchy then yields

```plain text
A_p(H)
 <= L_3(H)/p^2 + sqrt(K L_3(H))
 << H^2/p + H sqrt(p).
```

At full scale H about p, this gives

```plain text
A_p << p^(3/2),
T_3 << p^(3/2).
```

It does not give ATR. It would nevertheless be a major theorem: the proved T3 interpolation lemma would then give

```plain text
F_p << p^(5/4).
```

More generally, an artificial strengthened hierarchy

```plain text
[COINC3-rho]
N_coinc^{(3)}(M-1)
 <= L_3(M-1)^2/p^2
    + K L_3(M-1)/p^rho
```

would give

```plain text
A_p << p^((3-rho)/2) + p.
```

Thus:

- rho=0 gives the natural second-moment scale p^(3/2);

- every rho>0 gives a strict saving below p^(3/2);

- rho=1 gives ATR scale O(p).

But rho=1 is not a random-map prediction. A random map from L samples to p^2 values has collision excess of order L, not L/p. Therefore Cauchy plus a natural global L^2 collision estimate loses exactly the square-root factor that separates p^(3/2) from p.

The realistic ATR-strength target should instead be one of the following signed/local statements:

```plain text
sum_{(t1,t2)!=(0,0)} B_{t1,t2}^{(3)}(M-1)
 << p^3,
```

which is exactly what the Fourier identity needs after division by p^2, or a direct local-value estimate

```plain text
#{a in T_{M-1}:
  (Delta_{r,h},Delta_{r,k})=(0,0)}
 << p.
```

The latter is ATR itself in value-vector language. This is a defining-value local limit problem, not merely an equidistribution-of-all-values problem.

# 8. The best unconditional maximum-fiber bound

First, energy alone gives only

```plain text
max_v m_p(v) <= sqrt(F_p) << p^(3/4).
```

The gap-degree structure improves this to exponent 2/3 for every fiber.

Fix one fiber with positions

```plain text
1<=n_1<...<n_m<=M,
```

and let

```plain text
g_i=n_{i+1}-n_i,
1<=i<=m-1.
```

For an integer H>=1, call a consecutive gap short if g_i<=H. The number of long gaps is at most

```plain text
(M-1)/(H+1),
```

because their total length is at most M-1. Every short consecutive gap is one of the collision pairs counted by C_h, so

```plain text
#{i:g_i<=H}
 <= sum_{h<=H} C_h
 <= 3 sum_{h=2}^H (h-1)
 = 3H(H-1)/2.
```

Consequently

```plain text
m-1
 <= (M-1)/(H+1) + 3H(H-1)/2.
```

Balancing the two terms at

```plain text
H about (M/3)^(1/3)
```

gives

```plain text
m
 <= ((3/2) 3^(1/3) + o(1)) M^(2/3)
 < 3 M^(2/3)
```

for large M; the remaining bounded range is covered trivially. Hence

```plain text
max_v m_p(v) <= 3 p^(2/3).
```

This explains precisely why the suggested p^(1/2) does not follow. The available inequality is

```plain text
m << p/H + H^2,
```

and its optimum is a cube-root balance, not a square-root balance. Counting all pairwise gaps only gives

```plain text
binom(m,2) <= sum_{h<=diameter} C_h << diameter^2,
```

which is essentially tautological and supplies no extra factor.

The known |Z_p|<=3p^(2/3) concerns one distinguished fiber. The argument above is the required uniform extension to every projective fiber; it uses only the collision polynomial degree bound.

# 9. Three immediate unconditional triple bounds

The exact identity A_p=T_3/3 reduces the problem to upper-bounding T_3.

## 9.1 Energy alone

Monotonicity of finite l^q norms gives

```plain text
T_3 <= sum_v m_p(v)^3
    <= F_p^(3/2)
    << p^(9/4).
```

## 9.2 Energy plus the uniform maximum

Since

```plain text
(m)_3 <= (max_v m_p(v)) (m)_2,
```

one has

```plain text
T_3
 <= (max_v m_p(v)) (F_p-M)
 << p^(2/3) p^(3/2)
 = p^(13/6).
```

This is the previously recorded current estimate.

## 9.3 Using the full gap profile: p^2 log p

The degree budget for every individual gap gives a stronger result.

### Theorem

For every prime p>=7,

```plain text
T_3 <= 36 M^2 (1+log M),
```

and therefore

```plain text
sum_{h!=k} J_p(h,k)
 <= 12 M^2 (1+log M)
 << p^2 log p.
```

### Proof

Fix a fiber of size m>=3, with consecutive gaps g_1,...,g_{m-1}. Put

```plain text
q=m-1.
```

Since

```plain text
sum_i g_i <= M-1,
```

at least q/2 of the consecutive gaps satisfy

```plain text
g_i <= 2M/q.
```

Select any such set of at least q/2 short consecutive gaps. For every selected gap h,

```plain text
q <= 2M/h.
```

Because q>=2, also

```plain text
m=q+1 <= 3M/h.
```

Now

```plain text
(m)_3 = m q(q-1) <= m q^2.
```

There are at least q/2 selected gaps, so the fiber's whole factorial third moment can be charged to its selected gaps with charge at most

```plain text
2 m q <= 12 M^2/h^2
```

per selected gap. Summing over all fibers, the number of selected consecutive collision pairs of gap h is at most the total collision count C_h. Therefore

```plain text
T_3
 <= 12 M^2 sum_{h=2}^{M-1} C_h/h^2.
```

Using C_h<=3(h-1) gives

```plain text
T_3
 <= 36 M^2 sum_{h=2}^{M-1} (h-1)/h^2
 <= 36 M^2 (1+log M).
```

This proves the theorem.

The logarithm is exactly the harmonic cost of summing the scale-by-scale degree budget. The energy, maximum-fiber, and gap-degree inputs do not presently remove it.

# 10. What reflection and parity add

The orbit reflection

```plain text
pi(p-1-n)=pi(n)
```

implies:

- every fiber not containing the midpoint has even multiplicity;

- the unique fiber containing (p-1)/2 has odd multiplicity;

- every even gap has one forced mirror collision.

Write schematically

```plain text
C_h = kappa_h + C_h^prim,
kappa_h=1 for h even and 0 for h odd.
```

In the weighted proof above, the entire forced layer contributes only

```plain text
M^2 sum_h kappa_h/h^2 = O(M^2).
```

The primitive part still has the harmonic O(M^2 log M) allowance. Thus reflection improves constants and explains parity, but does not change the exponent or remove the logarithm. The even-gap factorization law likewise removes only a bounded-degree structural factor from each even N_h; the remaining degree is still linear in h.

Therefore the best unconditional conclusion from all four named inputs is

```plain text
sum_{h!=k} J_p(h,k) << p^2 log p,
```

not o(p^(3/2)).

# 11. Pair-COINC amplification: two short gaps plus one long

Assume that for every H<=p^beta,

```plain text
S_p(H)=sum_{h<=H} C_h
 << H + sqrt(pH).
```

There is a clean general amplification lemma.

### Lemma

For every 1<=H<=M, put

```plain text
R=ceil(M/H).
```

Then

```plain text
T_3
 <= 4 M R^2 + 8 (R S_p(H-1))^(3/2).
```

### Proof

Partition I_p into R consecutive blocks of length at most H. For each fiber v, let a_j(v) be its multiplicity in block j, and put

```plain text
P_v = sum_j binom(a_j(v),2).
```

Every pair counted by P_v has gap at most H-1, hence

```plain text
sum_v P_v <= S_p(H-1).
```

Split the fibers into two classes.

For m_v<2R,

```plain text
(m_v)_3 <= m_v^3 <= 4R^2 m_v.
```

Summing gives a contribution at most 4MR^2.

For m_v>=2R, Cauchy across the blocks gives

```plain text
m_v^2
 <= R sum_j a_j(v)^2
 = R(m_v+2P_v).
```

Since R m_v<=m_v^2/2, this implies

```plain text
m_v^2 <= 4R P_v.
```

Therefore

```plain text
sum_{high v} m_v^2
 <= 4R S_p(H-1).
```

Finally, monotonicity of l^q norms gives

```plain text
sum_{high v} (m_v)_3
 <= sum_{high v} m_v^3
 <= (sum_{high v} m_v^2)^(3/2)
 <= 8(R S_p(H-1))^(3/2).
```

Adding the two classes proves the lemma.

### Exponent calculation

Take H about p^beta. Since

```plain text
R << p^(1-beta)
```

and, for beta<=1,

```plain text
S_p(H) << p^((1+beta)/2),
```

the lemma gives

```plain text
T_3
 << p^(3-2beta)
    + p^(9/4-3beta/4)
    + p^(3/2).
```

The dominant exponent is

```plain text
3-2beta                 for beta<=3/5,
9/4-3beta/4             for beta>=3/5.
```

Concrete landmarks are:

```plain text
beta=1/2:  T_3 << p^2,
beta=3/5:  T_3 << p^(9/5),
beta=2/3:  T_3 << p^(7/4),
beta=1:    T_3 << p^(3/2).
```

Thus:

- beta=1/2 removes the unconditional logarithm but gives no power saving below p^2;

- every beta>1/2 gives a fixed power saving below p^2, and hence a new energy theorem through T3 interpolation;

- beta=2/3 gives the strong triple estimate T_3<<p^(7/4);

- no beta<=1 makes this amplification o(p^(3/2));

- ATR requires an additional long-range or higher-order input.

For energy itself, pair COINC should be used directly rather than routed through triples. The direct block bound is

```plain text
F_p
 << p^(2-beta) + p^(3/2-beta/2) + p,
```

which at beta=2/3 gives F_p<<p^(4/3). Passing through T_3<<p^(7/4) and T3 interpolation gives only F_p<<p^(11/8), so the triple amplification is weaker for energy but valuable as progress toward ATR.

# 12. Recommended theorem/ledger structure

The rigorous report should record the following items separately.

- [PROVED-D1] Exact determinant formula:

```plain text
Delta_{r,h}=N_h(r)/product_{j=1}^h(r+j)^3.
```

- [PROVED-D2-D4] Pair orthogonality, Parseval, and COINC implication.

- [CORRECTION-D5] Replace the false sigma-ladder inference by the direct block consequence F_p<<p^(4/3) at H=p^(2/3).

- [CORRECTION-D6] Generic curve error O(H^4 sqrt(p)); component classification remains open.

- [PROVED-T3-IDENTITY]

```plain text
3 sum_{h!=k}J_p(h,k)=sum_v(m_p(v))_3.
```

- [PROVED-MAX-2/3] Uniform max_v m_p(v)<<p^(2/3); no p^(1/2) follows.

- [PROVED-T3-2LOG]

```plain text
sum_{h!=k}J_p(h,k)<<p^2 log p.
```

- [CONDITIONAL-COINC3-0] Natural vector-collision variance gives T_3<<p^(3/2) and hence F_p<<p^(5/4).

- [CONDITIONAL-PAIR-AMPLIFICATION] Pair COINC through p^beta gives the explicit exponent formula above.

- [OPEN-ATR-LOCAL] The remaining ATR target is a local estimate at the single vector value (0,0), or equivalent signed cancellation over the whole two-frequency plane. Natural global L^2 equidistribution alone stops at the p^(3/2) barrier.

## Least-confident point

The algebraic identities, the p^2 log p proof, and the amplification lemma are elementary and I regard them as firm. The least-confident part is not a proof step but an optimality assessment: I have shown that the listed scalar inputs yield p^2 log p by a sharp multiscale charge, but I have not constructed an actual Apery orbit saturating that bound. A further deterministic identity among nonconsecutive gaps could conceivably remove the logarithm. Nothing in the current mounted theorem set supplies such an identity.

## Main blind spot

The shared blind spot is to expect a natural second-moment collision theorem to count one exact value at its random mean. In the trilinear setting the domain has size about p^3 and the value space has size p^2. Natural global L^2 control has square-root error p^(3/2), while the exact zero-vector main term is only p. Bridging that factor requires a local limit theorem, signed frequency-plane cancellation, or a new arithmetic rigidity statement; it is not delivered by ordinary COINC alone.