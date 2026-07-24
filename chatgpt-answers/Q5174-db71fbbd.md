ANSWER Q5174 db71fbbd

# P2.7 in the Conservative Matrix Field framework

## Executive verdict

There are two different objects that should not be conflated.

1. The P2.7 recurrence has a completely explicit, canonical one-dimensional CMF of rank 3: its ordinary companion matrix. This matrix is written below, and its product generates both p_n and q_n exactly from their three initial values.

1. A more informative, genuinely higher-dimensional discovery CMF may have been used by the Ramanujan Machine group to produce P2.7. No such P2.7-specific field, base point, direction, gauge, or boundary vectors are printed in the public CMF paper, the public challenge statement, or the current Ramanujan_Challenge repository.

The public paper On Conservative Matrix Fields: Continuous Asymptotics and Arithmetic develops the general framework and generic constructions. The 2026 challenge paper The Ramanujan Challenge For AI lists P2.7 among the problems whose proof or proof procedure is known to the authors but not yet public. The current project file problems/2.7/proof.tex contains the recurrence, initial data, AESZ #209 fingerprint, and the high-precision adjoint certificate, but no native CMF matrix or CMF integral.

The precise answers are therefore:

- (a) The publicly reconstructible CMF is the rank-3 companion matrix in §2. A polynomial-entry representative is also given. The native higher-dimensional matrix is not publicly specified.

- (b) Conservativity means a discrete flatness/cocycle identity, not an equation such as ∂_xM+∂_yM=0. It gives path independence. It does not imply that the marked error vector is recessive. In one dimension conservativity is automatic, which already proves that it cannot by itself force c₀(e)=0.

- (c) The CMF paper does not give a P2.7-specific integral for ζ(2)+ζ(3). The elementary mixed double integral is known, but no public identity connects its kernel to the P2.7 matrix product.

- (d) The AESZ #209 Hadamard split is visible in companion-matrix language as an exact diagonal coboundary/gauge. No corresponding factorization of the P2.7 companion field has been established.

- (e) The huge initial values are the two boundary covectors of the CMF ratio. They are not intrinsic matrix entries. Their size is gauge- and normalization-dependent. Public sources do not disclose the smaller native vector or gauge, if one exists.

---

# 1. What “conservative” means

Let K be a difference field with commuting shifts σ₁,…,σ_d. A CMF of dimension d and rank r consists of matrices

satisfying the discrete zero-curvature equations

for every pair i,j. Equivalently, there is a multiplicative cocycle M_v(x) such that

Consequently the product between two lattice points depends only on the endpoints, not on the chosen lattice path.

This is the matrix analogue of a closed discrete one-form. It is not the differential condition

For a fixed base point x and direction v, the trajectory product is

A CMF ratio has the form

## The crucial one-dimensional observation

When d=1, there are no nontrivial pairs of coordinate directions in (1.1). Hence every generically invertible rational matrix M(n) defines a one-dimensional CMF by

Thus any finite-order linear recurrence is automatically a CMF trajectory through its companion matrix. This is explicitly part of the general framework: ordinary D-finite/Apéry ratios are special cases of CMF ratios.

This observation is also the simplest reason that conservativity alone cannot prove convergence or recessiveness. An arbitrary recurrence, including one whose marked error has a nonzero dominant component, is already a one-dimensional CMF.

---

# 2. The exact rank-3 P2.7 companion CMF

Write the P2.7 recurrence as

where

For m≥0, define the row state

Set

Then

with the explicit rank-3 matrix

Indeed, multiplying (u_m,u_{m+1},u_{m+2}) by the three columns gives

The determinant is

This is the one-step Casoratian/volume factor. It determines how three-dimensional volume changes; it does not determine the dominant line or the recessive plane.

## 2.1 Exact product formula for p_n and q_n

Let

with

For N≥2, put

with the empty product for N=2. Then

Therefore

This is exactly a CMF ratio with

- CMF dimension d=1;

- matrix rank r=3;

- base point m=0;

- direction +1;

- left boundary covectors \mathbf p^T,\mathbf q^T;

- common right boundary vector e_3.

No 4×4 augmentation is needed because the recurrence is homogeneous and of shift order three. A 4×4 matrix would only be natural after adding an affine forcing term or an additional marked period.

If a column-state convention is preferred, transpose (2.6).

## 2.2 A polynomial-entry representative

The exact transfer matrix (2.6) has rational-function entries, which is allowed in the CMF definition. If one wants polynomial entries, set

and define

Explicitly,

Every entry is a polynomial in m. Moreover,

and the identical scalar multiplies q_N. Thus the projective ratio is unchanged:

The rational matrix is the exact state transfer; the polynomial matrix is its projectively equivalent cleared representative.

## 2.3 Exact verification code

```python
# SageMath
from sage.all import QQ, matrix, vector


def A(n):
    n = QQ(n)
    return (1024*(2*n+5)^4*(2*n+7)^3*(2*n+9)^3
            *(946*n^2+6407*n+10860))


def B(n):
    n = QQ(n)
    return (128*(2*n+7)^3*(2*n+9)^3
            *(104060*n^6+1745370*n^5+12145238*n^4
              +44886481*n^3+92943995*n^2
              +102256019*n+46709052))


def C(n):
    n = QQ(n)
    return (16*(n+3)^4*(2*n+9)^3
            *(3784*n^5+57792*n^4+351019*n^3
              +1059230*n^2+1587211*n+944620))


def D(n):
    n = QQ(n)
    return (n+3)^4*(n+4)^6*(946*n^2+4515*n+5399)


def M(m):
    m = QQ(m)
    return matrix(QQ, [
        [0, 0,  D(m)/A(m)],
        [1, 0, -C(m+1)/A(m+1)],
        [0, 1,  B(m+2)/A(m+2)],
    ])


def term(initial, N):
    if N < 3:
        return QQ(initial[N])
    row = vector(QQ, initial)
    for m in range(N-2):
        row = row * M(m)
    return row[2]


p_init = [
    QQ(-612218384750),
    QQ(-9525021973931919)/18100,
    QQ(-29561828382772029)/65380,
]
q_init = [
    QQ(-215040420000),
    QQ(-167282265043404)/905,
    QQ(-964185327658080)/6071,
]

p = [term(p_init, n) for n in range(10)]
q = [term(q_init, n) for n in range(10)]

assert p[:3] == p_init
assert q[:3] == q_init
for m in range(7):
    assert vector(QQ, q[m:m+3]) * M(m) == vector(QQ, q[m+1:m+4])
    assert vector(QQ, p[m:m+3]) * M(m) == vector(QQ, p[m+1:m+4])

print(p)
print(q)
```

---

# 3. Is this the CMF that discovered P2.7?

Not necessarily.

The companion construction above is canonical once the scalar recurrence and indexing convention are fixed, but it is mathematically tautological: it packages the recurrence into first-order matrix form. It does not explain why these coefficients, these initial values, or the constant ζ(2)+ζ(3) were discovered.

A native discovery CMF would normally require the following additional data:

1. a lattice dimension d≥2;

1. generator matrices M_1(x),…,M_d(x);

1. an exact flatness proof (1.1);

1. a base point x_0 and trajectory direction v;

1. four boundary vectors defining the CMF ratio;

1. a gauge/coboundary that converts the trajectory matrix into the P2.7 companion form;

1. an analytic argument identifying the limiting ratio and its rate.

The 2025 CMF paper proves a general construction from finite-dimensional Ore modules of D-finite functions to CMFs. To apply that construction to P2.7 one must specify the source D-finite function, the Ore-module basis, and the shift representations. None of those P2.7-specific ingredients is printed in that paper.

The public challenge paper gives the recurrence and initials but not the hidden proof procedure. It explicitly states that the problems in Section 2 have proofs or proof procedures known to the authors that are not yet public. Therefore it would be incorrect to manufacture a two-dimensional field and attribute it to Leibtag or to the CMF paper.

There are infinitely many gauge-equivalent CMFs even after one native field is known. Hence the phrase “the CMF matrix” is not unique without a normalization convention.

---

# 4. Why conservativity does not imply c₀(e)=0

Let

Then the P2.7 error is

Suppose the trajectory cocycle has a dominant/recessive splitting. In a suitable normalization its product has an asymptotic rank-one leading term

where r_0 is the dominant direction transported back to the initial fiber. Then

Consequently

and in particular

This is the same condition expressed by the adjoint Lagrange bracket:

Flatness says that products along homotopic lattice paths with the same endpoints agree. It does not say that the specific boundary covector \mathbf p-L\mathbf q annihilates the dominant direction. That is additional marked boundary data.

## A decisive logical test

In dimension one, every recurrence companion matrix is conservative. Choose any initial vector with a nonzero dominant component. Its error will be dominant, despite conservativity. Therefore

The public CMF paper’s broad higher-dimensional asymptotic principles are either conditional on convergence/spectral hypotheses or formulated as conjectures. They are not a theorem that every CMF ratio converges, much less that a prescribed zeta-linear form is recessive.

## What a CMF-native proof would look like

A genuinely CMF-native proof could work if the unreleased field supplies one of the following stronger structures.

### Route 1: two-path evaluation

Use flatness to identify the P2.7 trajectory product with a second path whose factors admit an exact boundary limit. One must then prove

with a rigorous uniform error bound along the second path.

### Route 2: invariant covector or normal function

Construct an exact dual section w of the CMF satisfying a conserved pairing and prove

symbolically, for example by a period, regulator, or boundary-cycle identity.

### Route 3: gauge to a known Barnes system

Exhibit an explicit coboundary and projection carrying a source CMF with a known integral error—such as Zudilin’s simultaneous ζ(2),ζ(3) construction—to the P2.7 marked vectors. The same map must carry both the denominator and numerator, not just the characteristic polynomial.

Any of these would prove (4.6). Discrete flatness alone does not.

---

# 5. Does the CMF paper give the P2.7 integral?

No public P2.7-specific integral is supplied in the CMF paper.

The project already records the elementary identity

This follows by expanding (1-xy)^{-1} and integrating termwise. It identifies the target constant, but it does not identify the P2.7 error.

To close the proof one needs an N-dependent kernel K_N such that

plus a manifest estimate

No such K_N is constructed in the public CMF paper for P2.7.

The paper’s Ore-module construction runs primarily in the direction

It is not a general inverse theorem turning an arbitrary CMF into a definite integral. A native source function for P2.7 might provide such an integral, but that source is not publicly identified.

---

# 6. The AESZ #209 Hadamard split in matrix language

Let

so the AESZ #209 holomorphic period is

The scalar factor satisfies

Suppose the inner sequence has row-state matrix M_S(m), so

Define the diagonal matrix

Since

the outer AESZ companion matrices are related to the inner matrices by

This is exactly a one-dimensional CMF coboundary/gauge transformation. It is the matrix-field shadow of the coefficientwise Hadamard identity

At the level of differential modules, the same operation is a Kummer/Hadamard convolution and can change differential rank.

## Is this split visible in the P2.7 matrix?

Not at present.

The P2.7 companion matrix (2.6) has its own determinant and coboundary/Pochhammer fingerprints. Those are compatible with nontrivial gauges, but they do not provide an equality of matrix fields. The project has already verified that the P2.7 recurrence and the standard AESZ/Zudilin recurrence are different and that their marked solutions do not agree after the obvious scalar gauges.

To claim that (6.7) explains P2.7, one would need an explicit identity such as

where Π is a stated projection/direct-image construction and G(m) is explicit. No such identity is public.

Thus the Hadamard decomposition is clearly visible in a CMF built from the AESZ inner/outer pair, but it is not presently visible as a factorization of the P2.7 trajectory.

---

# 7. Where the huge initial values occur

In the canonical companion CMF, the large numbers appear only in the boundary covectors

They are not entries of M(m). The same trajectory matrix produces both sequences; the different left boundary vectors select different solutions.

## 7.1 Gauge dependence

Let G(m)∈GL_3(\mathbb Q(m)) and change row coordinates by

Then the transfer matrix becomes

The boundary data transform as

Therefore a small native boundary vector can become a vector with huge rational entries after conversion to the scalar companion basis, and conversely. The magnitude of the coordinates has no invariant meaning.

A common nonzero scaling of both p and q also leaves every ratio p_N/q_N unchanged. Large common factors are often chosen to clear denominators or to obtain integral arithmetic data.

## 7.2 What is publicly known

The public challenge statement supplies the initial values but does not explain their derivation. The public CMF paper does not contain these numbers. The current project repository treats them as input data.

Accordingly, the strongest justified statement is:

In the companion CMF, the P2.7 initials are the left boundary covectors. They were presumably obtained from a more natural native construction and then normalized or gauged into the printed recurrence, but the public sources do not identify that native vector or gauge.

This is not a defect of the companion construction; it is precisely the missing discovery metadata.

---

# 8. What must be supplied to turn the CMF origin into a proof

A complete CMF certificate for P2.7 should include the following finite list of data and checks.

## Algebraic data

1. The native dimension d and rank r.

1. Explicit generators M_1(x),…,M_d(x).

1. Exact verification of every flatness identity (1.1).

1. Base point x_0, direction v, and trajectory matrices.

1. Boundary vectors defining the numerator and denominator.

1. An explicit gauge/projection relating the native trajectory to (2.6).

1. Exact matching of the three printed initial values of both p and q.

## Analytic data

1. A theorem identifying the CMF ratio limit with (5.1), a Barnes integral, or a regulator period.

1. A uniform bound showing that the marked error has only the two recessive Poincaré modes.

1. Equivalently, an exact proof of the boundary orthogonality (4.6) or bracket identity (4.7).

Once these are supplied, the CMF framework can indeed furnish the desired native proof. Without items 8–10, the matrix product merely reproduces the recurrence.

---

# 9. Direct answers to the five questions

## (a) Explicit matrix

The exact publicly reconstructible matrix is the 3×3 companion matrix (2.6), or the polynomial representative (2.15). Its product with the two boundary covectors generates p_n,q_n exactly. No P2.7-specific higher-dimensional CMF is printed in the public CMF paper or current project repository.

## (b) Does conservativity imply recessiveness?

No. It implies discrete path independence. In one dimension it is automatic. The missing statement is the marked orthogonality

or the equivalent exact adjoint-bracket identity. A special higher-dimensional second-path argument could prove it, but flatness by itself does not.

## (c) Integral representation

The CMF paper does not give a P2.7-specific integral. The elementary integral (5.1) gives the target constant, not the N-dependent P2.7 error. An exact kernel (5.2) remains necessary.

## (d) AESZ Hadamard split

Yes, for the AESZ inner/outer pair it is exactly the diagonal coboundary (6.7). No explicit factorization or projection identifies this pair with the P2.7 companion matrix.

## (e) Huge initials

They are boundary covectors in (2.8), not matrix entries. Their size depends on the cyclic basis, gauge, and denominator-clearing normalization. Their native derivation is not public.

---

# References

1. S. Weinbaum, E. Leibtag, R. Kalisch, M. Shalyt, and I. Kaminer, On Conservative Matrix Fields: Continuous Asymptotics and Arithmetic, 2025.

1. M. Shalyt et al., The Ramanujan Challenge For AI, 2026; P2.7 is contributed by Elyasheev Leibtag and is listed in the section whose proof procedures are known to the authors but not yet public.

1. Current project state: problems/2.7/proof.tex at commit d702a470….

1. W. Zudilin, Approximations to -, di- and tri-logarithms, for the inner simultaneous ζ(2),ζ(3) construction underlying the AESZ #209 Hadamard period.

1. G. Almkvist, C. van Enckevort, D. van Straten, and W. Zudilin, Tables of Calabi–Yau equations, for AESZ #209.