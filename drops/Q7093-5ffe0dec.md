ANSWER Q7093 5ffe0dec

# Audit of central-interval removal of carrier exceptions

## Claim that can be proved

The endpoint identity gives a precise source of some carrier factors. If
\[
N_h(-s)=(-1)^{s-1}b_{s-1}b_{h-s}((s-1)!(h-s)!)^3,
\]
then for a prime p dividing a carrier value \(b_j\), a forced endpoint zero occurs at the residue
\[
x\equiv -j\pmod p
\]
whenever the corresponding index appears as an endpoint of the gap polynomial. More explicitly, if \(1\le j\le h\) and \(p\mid b_j\), then
\[
N_h(-j-1)\equiv 0\pmod p
\]
from the endpoint formula with \(s=j+1\), and similarly the other endpoint factor gives roots near \(-h+j\).

For a gap quadruple of total length at most H, these forced roots occur in the translated endpoint strips. Hence a purely endpoint-generated common root is confined to residues within O(H) of the wrap boundary (equivalently near 0 or p modulo p).

Therefore restricting to
\[
H\le x\le p-2H
\]
does remove these specific endpoint witnesses.

## What is not proved

The stronger statement in the question is false without an additional hypothesis:

> Every common-root contribution from a carrier prime is forced by an endpoint root.

The endpoint identity only exhibits some roots. It does not classify all roots modulo p. If \(p\mid b_j\) and a resultant
\[
S_{d,r}=\operatorname{Res}(N_d(x),N_r(x+d))
\]
vanishes modulo p, the common root \(\alpha\) may be a different root of both polynomials. The implication
\[
p\mid b_j\Longrightarrow \text{all common roots lie at }-j+O(H)
\]
does not follow.

## Three cases

### (i) Forced boundary root

Here one polynomial vanishes because of the endpoint identity. For example,
\[
\alpha=-j-1,
\]
with
\[
p\mid b_j.
\]
This is a carrier contribution. It is excluded by the central interval condition when the zero start variable is interpreted in the non-wrapping interval.

### (ii) Genuine central common root at the same carrier prime

It is possible that the same prime p divides a carrier value and also satisfies
\[
N_d(\alpha)=0,
\qquad
N_r(\alpha+d)=0
\]
for some central \(\alpha\). The resultant only records existence of a common root:
\[
p\mid S_{d,r}\iff \exists\alpha\in\overline{\mathbf F}_p.
\]
It does not remember whether the root came from an endpoint factor.

Thus a central common root contributes to the primitive/excess resultant even when p is in the carrier set.

### (iii) Multiplicity overlap

If the endpoint root has multiplicity, then the valuation of the resultant can contain both structural and witness contributions. Locally,
\[
v_p(S_{d,r})
\]
is controlled by the intersection multiplicity of the two divisor schemes. A boundary factor and a genuine intersection at the same prime need not separate additively.

## Consequence for primitive resultants

Deleting every prime dividing a carrier \(U_H\) is too aggressive. It removes:

1. pure boundary factors (desired), and
2. possible central primitive intersections at the same prime (undesired).

A valid replacement would require an exact local factorization:
\[
v_p(S_{d,r})=v_p(S_{d,r}^{\rm endpoint})+v_p(S_{d,r}^{\rm excess}),
\]
where the first term is computed from the known endpoint divisors. No such factorization follows from the endpoint formula alone.

## Remaining obstruction

The missing lemma is a local separation theorem:

For every p and every \(d,r\le H\), every common root of
\[
N_d(x),\quad N_r(x+d)
\]
modulo p either

1. is one of the explicit endpoint roots \(x\equiv -j\) or its reflected boundary translate, or
2. contributes to a primitive resultant factor after removing the exact endpoint intersection multiplicity.

The root-strip theorem over characteristic zero does not imply this modulo p. Extension-field roots can appear without any complex endpoint analogue. Therefore the central restriction removes known boundary witnesses but does not by itself remove carrier-prime exceptions.

The quadruple-corank argument still needs a weighted control of carrier primes or an exact local resultant decomposition; the central interval alone is insufficient.