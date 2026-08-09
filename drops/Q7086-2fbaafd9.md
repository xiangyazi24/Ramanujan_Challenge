ANSWER Q7086 2fbaafd9

# Q7086: Generic coprimality of consecutive reduced gap resultants

## Short answer

The proposed root-strip argument gives a useful obstruction, but it does **not by itself prove**
\[
\gcd(S^*_{h_1,h_2},S^*_{h_2,h_3})=1
\]
for generic triples. The reason is that a common prime divisor of two resultants does **not require the same geometric root** to occur in the two collisions.

What can be proved is a stronger structural statement in the "same root" case, and a reduction of the remaining case to an independent arithmetic multiplicity problem.

## 1. Resultant interpretation

For a prime \(\ell\nmid\) structural content,
\[
\ell\mid S^*_{h_1,h_2}
\]
means that over \(\overline{\mathbf F}_\ell\) there exists
\(\alpha\) such that
\[
N_{h_1}(\alpha)=0,
\qquad
N_{h_2}(\alpha+h_1)=0.
\]
Similarly,
\[
\ell\mid S^*_{h_2,h_3}
\]
means that there exists \(\beta\) such that
\[
N_{h_2}(\beta)=0,
\qquad
N_{h_3}(\beta+h_2)=0.
\]

The two witnesses \(\alpha\) and \(\beta\) need not coincide.

Therefore the gcd question splits into:

1. the triple-root case \(\alpha=\beta\);
2. the independent double-collision case \(\alpha\ne\beta\).

Only the first is controlled by strip geometry.

## 2. The same-root case

Assume \(\alpha=\beta\). Then
\[
N_{h_1}(\alpha)=0,
\]
\[
N_{h_2}(\alpha+h_1)=0,
\]
\[
N_{h_3}(\alpha+h_2)=0.
\]

Over characteristic zero the strip theorem says roots of \(N_h(x)\) avoid the boundary regions and lie in the allowed Apéry strip
\[
-h < \operatorname{Re}(x)<-1.
\]

Applying this formally gives:

- from \(N_{h_1}(\alpha)=0\):
\[
-h_1<\operatorname{Re}\alpha<-1;
\]

- from \(N_{h_2}(\alpha+h_1)=0\):
\[
-h_2<\operatorname{Re}(\alpha+h_1)<-1,
\]
hence
\[
-h_1-h_2<\operatorname{Re}\alpha<-h_1-1;
\]

- from \(N_{h_3}(\alpha+h_2)=0\):
\[
-h_3<\operatorname{Re}(\alpha+h_2)<-1,
\]
hence
\[
-h_2-h_3<\operatorname{Re}\alpha<-h_2-1.
\]

The first two strips are already disjoint:
\[
(-h_1,-1)\cap(-h_1-h_2,-h_1-1)=\varnothing.
\]

Thus in characteristic zero there cannot be a genuine triple root. This proves the geometric statement:

> A single complex point cannot witness two consecutive gap collisions.

## 3. Why this does not prove gcd = 1

The obstruction is reduction modulo \(\ell\).

The root-strip theorem is an Archimedean statement. A root \(\alpha\in\overline{\mathbf F}_\ell\) has no canonical real part. A Teichmüller lift only gives a preferred lift for units, and the lift does not preserve the zero-set avoidance statement automatically.

Therefore the implication
\[
\text{three common roots mod }\ell
\Rightarrow
\text{overlapping strips}
\]
is not currently justified.

The correct algebraic replacement would require proving a characteristic-free coprimality statement such as
\[
\gcd(N_{h_1}(x),N_{h_2}(x+h_1),N_{h_3}(x+h_2))=1
\]
over \(\mathbf Z[x]\) after removing structural factors. That is a much stronger theorem.

## 4. What remains for generic gcd

For a common large prime \(\ell\), we only know that there are two (possibly different) collision points:

\[
\alpha:\quad h_1\to h_2,
\]

\[
\beta:\quad h_2\to h_3.
\]

Hence \(\ell\) divides the product of two independent intersection resultants. To obtain
\[
\gcd(S^*_{h_1,h_2},S^*_{h_2,h_3})=1
\]
one needs an additional genericity input, for example:

* irreducibility/transversality of the reduced gap-polynomial family;
* a Hilbert-irreducibility style argument in the height parameters;
* an explicit resultant-of-resultants computation showing that the common discriminant locus is the structural locus.

Without such an input, the root-strip argument only eliminates the diagonal triple-intersection component.

## 5. Consequence for GM-large

The good news is that the GM-large estimate does not actually need exact gcd = 1.

For any pair,
\[
\sum_{\ell>L}v_\ell(\gcd(S_1,S_2))\log\ell
\leq
\log|\gcd(S_1,S_2)|
\leq
\min(\log|S_1|,\log|S_2|).
\]

After structural saturation, the remaining common large-prime contribution is controlled by the reduced gcd height. Thus the "large prime" part is essentially free once the reduced gcd is shown to be small on average.

The missing theorem is therefore not the root-strip argument; it is an arithmetic transversality statement for the family \(N_h\).
