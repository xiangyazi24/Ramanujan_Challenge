ANSWER Q7541 cdee88fb

# Audit verdict

The corrected equivalence is **exactly right** in the stated range.  In fact, if
\[
Q:=a+2h,
\]
then for every prime \(p>Q\)
\[
 p\mid \gcd\bigl(b_a,N_h(a),N_h(a+h)\bigr)
 \quad\Longleftrightarrow\quad
 b_a\equiv b_{a+h}\equiv b_{a+2h}\equiv0\pmod p.
\tag{1}
\]
Equivalently,
\[
 \operatorname{rad}_{>Q}D_{a,h}
 =\operatorname{rad}_{>Q}\gcd(b_a,b_{a+h},b_{a+2h}),
\qquad
D_{a,h}:=\gcd(b_a,N_h(a),N_h(a+h)).
\tag{2}
\]
Thus excluding distinguished triples in the initial quarter is exactly the statement that no prime
\(p>4Q\) divides \(D_{a,h}\).

There **is** a useful nontrivial consequence of the known continuant algebra that should be retained.  Put
\[
 \mathcal R_h:=\operatorname{Res}_x\bigl(N_h(x),N_h(x+h)\bigr).
\tag{3}
\]
The current repository's complex root-strip theorem implies \(\mathcal R_h\ne0\) for \(h\ge2\), and the resultant height estimate gives
\[
 \log|\mathcal R_h|\ll h^2\log(2h).
\tag{4}
\]
Since a Sylvester Bezout identity gives
\[
 U_h(x)N_h(x)+V_h(x)N_h(x+h)=\mathcal R_h,
 \qquad U_h,V_h\in\mathbb Z[x],
\tag{5}
\]
one has the exact divisibility
\[
 \boxed{D_{a,h}\mid \mathcal R_h}
\tag{6}
\]
and hence
\[
 \boxed{\log D_{a,h}\ll h^2\log(2h).}
\tag{7}
\]
So for each **fixed** \(h\), quarter triples are excluded for all sufficiently large \(a\).  More generally, any quarter counterexample must satisfy
\[
 \log(a+2)\ll h^2\log(2h).
\tag{8}
\]
This is genuine information, but it is nowhere near the required uniform linear prime-support bound when \(h\) grows.

The supplied \((p,a,h)=(8431,1601,4)\) double-return witness shows exactly why the transfer-only route stops: it forces
\[
8431\mid \mathcal R_4,
\]
while
\[
8431>4(1601+8)=6436.
\]
Thus no identity involving only the two return polynomials can have a right side supported on primes \(\le4(a+2h)\) at this point.  The fact that the supplied actual Apéry value \(b_{1601}\not\equiv0\pmod{8431}\) is precisely the extra marked-initial-line information that removes this phantom.

The Casoratian does give an especially clean exact saturation identity, but it **only proves (1); it does not give a new scalar obstruction**.  Reflection likewise gives no general extra equation.  The current continuant identities reduce the problem to the separated resultant \(\mathcal R_h\), whose large-prime factors are genuinely new and are not boundary Apéry factors.

My scoped conclusion is therefore:

> **Proved no-go for the stated identity class.**  Recurrence transfer algebra, the Casoratian normalization, reflection, and the known continuant/Dodgson identities cannot by themselves rule out distinguished quarter triples by a controlled-support Bezout identity.  They are basis-covariant/projective identities.  At a separated double return there is a one-dimensional solution line having the three required zeros, and one can choose a companion so that the same Casoratian law holds.  The supplied \(8431\) return gives a concrete local countermodel with a prime larger than the desired support cutoff.
>
> This no-go is deliberately **not** a claim about an identity that uses the actual Apéry initial ratio \(b_0:b_1=1:5\), the binomial sum, the diagonal/constant-term model, or another genuinely Apéry-specific arithmetic theorem.  Such an input is exactly what is still missing.

The smallest exact missing statement is a **marked separated-return avoidance theorem**: on the common-root locus of \(N_h(x)\) and \(N_h(x+h)\), the actual Apéry initial line must avoid the zero line in the quarter.  A simple scalar sufficient strengthening is
\[
 \boxed{
 P^+\!\left(\gcd(b_a,\mathcal R_h)\right)
 \le 4(a+2h).
 }
\tag{9}
\]
Because \(D_{a,h}\mid\gcd(b_a,\mathcal R_h)\), (9) would close this route.  Nothing currently in the Casoratian/reflection/continuant package proves (9), and the repository's present marked-row analysis explains why companion labels alone do not add another equation.

No all-index Apéry smooth-radical theorem is proved here.

---

# 1. Repository state actually used

I audited the connected repository rather than relying on earlier conversational summaries.

Current `main` is

```text
8ee654c3cf80625e935460713579dd2dbc12a7a4
```

and I used in particular:

```text
problems/3.2/proof.tex
problems/3.2/AVENUES.md
problems/3.2/near_bridge_energy.tex
```

The relevant current facts are:

1. `proof.tex` contains the Apéry recurrence and Casoratian
   \[
   a_n b_{n-1}-a_{n-1}b_n=6/n^3,
   \]
   the gap-polynomial recurrence, the endpoint evaluations
   \[
   N_H(-r)=(-1)^{r-1}b_{r-1}b_{H-r}((r-1)!)^3((H-r)!)^3
   \quad(1\le r\le H),
   \tag{10}
   \]
   the adjacent-resultant factorization, and the continuant Bezout identity
   \[
   N_{e-1}(x+1)N_d(x)-N_{d-1}(x+1)N_e(x)
   =S_d(x)^6N_{e-d}(x+d),
   \tag{11}
   \]
   where \(S_d(x)=\prod_{j=2}^d(x+j)\).

2. `near_bridge_energy.tex` proves that the shifted resultants
   \[
   \mathcal R_{s,G}
   =\operatorname{Res}_x(N_s(x),N_G(x+s))
   \]
   are nonzero by disjoint complex root strips, and proves
   \[
   \log|\mathcal R_{s,G}|\ll sG\log(2H)
   \quad(s,G\le H).
   \tag{12}
   \]

3. `AVENUES.md` records the current terminal conclusion for the basis-covariant companion-label refinement.  At a distinguished zero the marked evaluation row is proportional to the fixed linear functional cutting out the actual initial line, and localization in the companion coordinate gives one marked pin plus the ordinary projective-return minors—no extra return-index equation.  This is directly relevant to the present three-return question.

The finite \(8431\) witness in the question is not present in the connected tracked files I found, so I do **not** claim to have independently recomputed the 1600-step Apéry congruence.  The algebraic equivalence and no-go below are independently derived; the numerical congruences supplied in the question are used only as the stated finite witness when specializing the obstruction.

---

# 2. Exact transfer formulation

Put
\[
 B_n=(n!)^3b_n
\]
and
\[
 M(y)=
 \begin{pmatrix}
 P(y)&-y^6\\
 1&0
 \end{pmatrix},
 \qquad
 P(y)=(2y+1)(17y^2+17y+5).
\tag{13}
\]
Then
\[
 \binom{B_{n+1}}{B_n}
 =M(n)\binom{B_n}{B_{n-1}}.
\tag{14}
\]
For \(h\ge1\), define
\[
 T_h(a):=M(a+h)M(a+h-1)\cdots M(a+1).
\tag{15}
\]
Thus
\[
 \binom{B_{a+h+1}}{B_{a+h}}
 =T_h(a)\binom{B_{a+1}}{B_a}.
\tag{16}
\]
The bottom-left entry of \(T_h(a)\) is precisely \(N_h(a)\).  Also
\[
 \det T_h(a)=\prod_{j=1}^h(a+j)^6.
\tag{17}
\]
Hence if \(p>a+h\), the \(h\)-step transfer is invertible modulo \(p\).

A useful full matrix formula is
\[
T_h(a)=
\begin{pmatrix}
N_{h+1}(a)&-(a+1)^6N_h(a+1)\\
N_h(a)&-(a+1)^6N_{h-1}(a+1)
\end{pmatrix},
\tag{18}
\]
with \(N_0=0\).  In particular, modulo a prime not meeting the local determinant factors,
\[
N_h(a)=0
\quad\Longleftrightarrow\quad
T_h(a)\text{ preserves the zero line }L:=\mathbb F_p(1,0)^t.
\tag{19}
\]
This projective interpretation is the core of both the positive equivalence and the no-go.

---

# 3. Independent proof of the corrected equivalence

Assume throughout this section that
\[
p>a+2h.
\tag{20}
\]
All factorials through \((a+2h)!\) are \(p\)-units, so
\[
 b_j\equiv0\pmod p
 \quad\Longleftrightarrow\quad
 B_j\equiv0\pmod p
\tag{21}
\]
for every index used below.

## 3.1 No two consecutive zero states

If \(1\le n<p\) and
\[
B_n\equiv B_{n+1}\equiv0\pmod p,
\]
then the recurrence at the previous step gives
\[
0=B_{n+1}=P(n)B_n-n^6B_{n-1},
\]
so \(B_{n-1}\equiv0\).  Iterating backwards gives \(B_0\equiv0\), contradicting \(B_0=1\).  For \(n=0\) the contradiction is immediate.

Therefore
\[
 B_n\equiv0\pmod p\quad\Longrightarrow\quad B_{n+1}\not\equiv0\pmod p
\tag{22}
\]
for every relevant \(n\).

## 3.2 First return

If \(B_a=0\), (16) gives
\[
 B_{a+h}=N_h(a)B_{a+1}.
\tag{23}
\]
By (22), \(B_{a+1}\ne0\).  Hence
\[
 B_a=B_{a+h}=0
 \quad\Longleftrightarrow\quad
 B_a=0\text{ and }N_h(a)=0.
\tag{24}
\]

## 3.3 Second return

Under these equivalent conditions, \(B_{a+h}=0\), so again by (22),
\[
B_{a+h+1}\ne0.
\]
Applying the same \(h\)-step formula from the new start gives
\[
 B_{a+2h}=N_h(a+h)B_{a+h+1}.
\tag{25}
\]
Therefore
\[
 B_a=B_{a+h}=B_{a+2h}=0
\]
is equivalent to
\[
 B_a=0,\qquad N_h(a)=0,\qquad N_h(a+h)=0.
\tag{26}
\]
Using (21) gives (1).

This proof uses only invertibility of the local recurrence and the initial fact \(B_0=1\).  It does not use the Casoratian.

For \(h=1\), \(N_1=1\), so both sides of the proposed equivalence are impossible; the statement still holds.

---

# 4. What the Casoratian actually gives

The Casoratian can be repackaged in a form that is exceptionally well adapted to the division-free recurrence.

Define
\[
 C_n:=\frac{(n!)^3a_n}{6}.
\tag{27}
\]
Then \(C_0=0,C_1=1\), and \(C_n\) satisfies the same integral recurrence as \(B_n\).  Moreover
\[
 C_n=N_n(0),
 \qquad
 B_n=N_{n+1}(-1).
\tag{28}
\]
The second equality is exactly the current repository evaluation
\(N_H(-1)=((H-1)!)^3b_{H-1}\).

Scaling the Apéry Casoratian gives
\[
 \boxed{
 C_nB_{n-1}-C_{n-1}B_n=((n-1)!)^6.
 }
\tag{29}
\]
There is then an exact endpoint-to-gap determinant identity
\[
 \boxed{
 C_{a+h}B_a-C_aB_{a+h}=(a!)^6N_h(a).
 }
\tag{30}
\]
Likewise,
\[
 \boxed{
 C_{a+2h}B_{a+h}-C_{a+h}B_{a+2h}
 =((a+h)!)^6N_h(a+h).
 }
\tag{31}
\]

These formulas are not a new source of arithmetic.  Formula (30) is literally the current continuant Bezout identity (11) specialized to
\[
x=-1,\qquad d=a+1,\qquad e=a+h+1.
\]
Indeed
\[
N_{a+h}(0)N_{a+1}(-1)-N_a(0)N_{a+h+1}(-1)
=(a!)^6N_h(a),
\]
which is exactly (30) after (28).

So the Casoratian and the boundary continuant Bezout identity are the **same two-dimensional determinant identity in different coordinates**.

## 4.1 It re-proves the equivalence, but nothing stronger

Suppose \(p>a+2h\) and \(p\mid B_a\).  From (29) at index \(a\) (with \(a\ge1\); \(a=0\) cannot be a zero),
\[
 C_aB_{a-1}\equiv((a-1)!)^6\not\equiv0\pmod p,
\]
so \(C_a\) is a unit modulo \(p\).  Reducing (30) then gives
\[
 N_h(a)\equiv0
 \quad\Longleftrightarrow\quad
 B_{a+h}\equiv0.
\tag{32}
\]
Once \(B_{a+h}=0\), (29) makes \(C_{a+h}\) a unit, and (31) gives
\[
 N_h(a+h)\equiv0
 \quad\Longleftrightarrow\quad
 B_{a+2h}\equiv0.
\tag{33}
\]
This is exactly (1) again.

The key point is that (30) has no nonzero scalar remainder.  On the return locus it merely transports the zero from one endpoint to the other.  There is no contradiction to extract.

---

# 5. Reflection: what it does and does not add

Assume the stated exact polynomial reflection
\[
 N_h(-x-h-1)=(-1)^{h-1}N_h(x).
\tag{34}
\]
If
\[
N_h(a)=N_h(a+h)=0\pmod p,
\]
then it also gives roots
\[
 -a-h-1,\qquad -a-2h-1
\tag{35}
\]
modulo \(p\).  Under the quarter hypothesis \(p>4(a+2h)\), these are represented by two points high in \([0,p-1]\), distinct from \(a,a+h\).

So reflection turns a separated double return into four roots of \(N_h\).  This gives one small positive result.

## Proposition — the cases \(h=1,2\)

For distinguished quarter triples, \(h=1\) and \(h=2\) are impossible.

* \(h=1\): \(N_1=1\).
* \(h=2\): \(N_2\) is a cubic.  A distinguished zero has \(a\ge1\); hence
  \(p>4(a+4)\ge20\), so \(p\ge23\) and the leading coefficient \(34\) remains nonzero modulo \(p\).  The four points
  \[
  a,\quad a+2,\quad p-a-3,\quad p-a-5
  \]
  are distinct under the quarter inequality, contradicting that a nonzero cubic has at most three roots.

For \(h\ge3\), however,
\[
\deg N_h=3(h-1)\ge6,
\]
so four reflected roots are perfectly compatible with the degree bound.  The supplied \(h=4\) witness shows that this compatibility is realized by the actual transfer system in the initial quarter.

The reflection of the Apéry solution does not add another condition either.  A lower triple reflects to the upper triple, but the two reflected gap starts are exactly the two roots already supplied by (34).  No new independent return equation appears.

## 5.1 Why the endpoint factorization cannot rescue reflection

The current manuscript has the genuinely Apéry-specific boundary identity (10).  For \(H=h\), however, reflecting a positive root \(a\ge0\) lands at
\[
-a-h-1=-r,
\qquad r=a+h+1>h.
\]
The factorization (10) only applies for
\[
1\le r\le h,
\]
i.e. on the finite boundary strip \(-h\le x\le-1\).  The reflected positive-root locus lies strictly beyond that strip.  Therefore the endpoint factorization that completely explains adjacent-resultant boundary factors gives **no factorization of \(N_h(a)\)** at the separated positive roots relevant here.

That is a precise reason the attractive boundary Apéry factors do not propagate into the present double-return problem.

---

# 6. The strongest unconditional continuant consequence: the separated resultant

For \(h\ge2\), set
\[
\mathcal R_h=\operatorname{Res}_x(N_h(x),N_h(x+h)).
\]
The root-strip result in `near_bridge_energy.tex` applies with \(s=G=h\): every complex root of \(N_h(x)\) has real part in \((-h,-1)\), whereas every complex root of \(N_h(x+h)\) has real part in \((-2h,-h-1)\).  Hence
\[
\boxed{\mathcal R_h\ne0.}
\tag{36}
\]

By the Sylvester adjugate there are \(U_h,V_h\in\mathbb Z[x]\) with (5).  Evaluating at an integer \(a\) gives
\[
 \gcd(N_h(a),N_h(a+h))\mid \mathcal R_h.
\tag{37}
\]
Therefore
\[
 D_{a,h}\mid \mathcal R_h
\]
with no localization and no omitted exceptional primes.

The same repository result gives
\[
\boxed{
 \log|\mathcal R_h|\ll h^2\log(2h).
}
\tag{38}
\]
Consequently any prime \(p\mid D_{a,h}\) satisfies
\[
 p\le P^+(\mathcal R_h)\le|\mathcal R_h|.
\tag{39}
\]
This proves:

* for each fixed \(h\), only finitely many primes can ever occur in any \(D_{a,h}\);
* for each fixed \(h\), the quarter target holds for all sufficiently large \(a\);
* a quarter counterexample must have \(h^2\log(2h)\gg\log(a+2)\).

This is the honest positive output of the known continuant machinery.

It is not the desired all-index result.  The bound (38) is exponential in roughly \(h^2\log h\), and there is no theorem forcing
\[
P^+(\mathcal R_h)=O(h)
\]
or even polynomially bounded prime support.  The supplied \(8431\) example demonstrates the issue already at \(h=4\).

The current manuscript's same-start Bezout identity also makes clear why there is no further universal Euclidean descent.  Taking \(d=h,e=2h\) in (11) reduces the pair \(N_h(x),N_{2h}(x)\), after cut-edge saturation, to the separated pair
\[
N_h(x),\quad N_h(x+h).
\]
The remaining shifted resultant is explicitly called genuinely new in `proof.tex`; the subtraction step does not iterate into boundary factors.

---

# 7. The best scalar coupling with the actual Apéry term

Combining the resultant with \(b_a\) gives the simple but important envelope
\[
\boxed{
 D_{a,h}\mid G_{a,h},
 \qquad
 G_{a,h}:=\gcd(b_a,\mathcal R_h).
}
\tag{40}
\]
This is the first place the actual Apéry initial line enters a scalar object.

There is even a literal three-generator Bezout identity.  Choose integers \(r,s\) with
\[
 r b_a+s\mathcal R_h=G_{a,h},
\]
and substitute (5):
\[
\boxed{
 r b_a+sU_h(a)N_h(a)+sV_h(a)N_h(a+h)=G_{a,h}.
}
\tag{41}
\]
So an exact Bezout coupling exists.

But (41) does **not** control the right side.  Its prime support is the new unresolved object
\[
\operatorname{Supp}\gcd(b_a,\mathcal R_h).
\]
The supplied \(8431\) witness is instructive: the local double return forces \(8431\mid\mathcal R_4\), while the supplied \(b_{1601}\not\equiv0\pmod{8431}\) shows that the marked Apéry factor removes this particular large resultant prime.  What is missing is a theorem saying that this always happens for the root **at the same start \(a\)** in the quarter.

The scalar gcd (40) also loses the common-root label.  A prime can divide \(\mathcal R_h\) because some residue \(x\) is a separated common root, and independently divide \(b_a\) at a different residue \(a\).  Hence (9) is a sufficient strengthening of the exact target, not an equivalent reformulation.  This root-label loss is the same issue highlighted in the current `AVENUES.md` discussion of scalar resultant-gcd formulations.

No known identity in the connected tree controls (40) at the required linear-prime scale.

---

# 8. Rigorous scoped obstruction from solution-space freedom

Here is the clean local no-go promised in the question.

Fix a prime \(p>a+2h\) for which
\[
N_h(a)=N_h(a+h)=0\pmod p.
\tag{42}
\]
Consider the two-dimensional \(\mathbb F_p\)-solution space \(\mathscr S\) of the same division-free recurrence on the index interval in question.  Because all one-step determinants \(j^6\) are nonzero for \(1\le j\le a+2h\), evaluation at time \(a\),
\[
 U\longmapsto (U_{a+1},U_a),
\]
is an isomorphism \(\mathscr S\to\mathbb F_p^2\).

Choose the nonzero solution \(U\) with
\[
(U_{a+1},U_a)=(1,0).
\tag{43}
\]
Then (42) and the definition of \(N_h\) give
\[
U_{a+h}=0,
\]
and because \(U_{a+h+1}\ne0\), the second condition gives
\[
U_{a+2h}=0.
\]
Thus
\[
\boxed{U_a=U_{a+h}=U_{a+2h}=0.}
\tag{44}
\]
So every separated double-return point carries a one-dimensional solution line with the corresponding triple zeros.  This is not a heuristic; it is elementary linear algebra in the exact recurrence.

Now choose any independent solution \(V\).  Its scale may be normalized so that the Casoratian of \((V,U)\) has whichever nonzero normalization one wants at one index.  The determinant recurrence then propagates that normalization to the standard \(\text{constant}/n^3\) Casoratian law.  The polynomial reflection and all continuant/Dodgson identities depend only on the recurrence coefficients, so they remain valid as well.  If one uses the current strong mod-\(p\) reflection theorem for solutions, the new pair satisfies that too.

Therefore the package

```text
same Apéry recurrence
+ nonzero Casoratian normalization
+ reflection
+ continuant/Dodgson identities
```

is perfectly compatible with a triple-zero solution whenever the separated double-return locus is nonempty.

## 8.1 Concrete controlled-support contradiction at the supplied witness

Specialize this to the supplied
\[
(p,a,h)=(8431,1601,4).
\]
Then \(Q=1609\) and
\[
4Q=6436<8431.
\]
Suppose there were a Bezout identity derived **universally from the basis-covariant transfer/Casoratian/continuant package**, valid for every solution \(U\), of the form
\[
 X\,U_a+Y\,N_h(a)+Z\,N_h(a+h)=R,
\tag{45}
\]
where every denominator in \(X,Y,Z\) is supported on primes at most \(4Q\), and where \(R\) is an integer (or localized unit) whose prime support is also at most \(4Q\).

Reduce (45) modulo \(8431\) and evaluate it on the solution (43).  All three terms on the left vanish, while \(R\) is a unit modulo \(8431\).  Contradiction.

Hence:

> **No controlled-support identity of this universal type exists.**

This is stronger than saying “the problem looks hard”: the supplied transfer witness is an explicit local countermodel to the whole basis-covariant identity class.

The scope is essential.  An identity whose coefficients genuinely use the distinguished Apéry initial ratio \(1:5\), or the binomial/constant-term description of \(b_n\), is not universal in \(U\) and is not ruled out by this countermodel.

---

# 9. Why the fixed companion does not secretly add another equation

The current `AVENUES.md` contains an exact marked-row calculation that is directly relevant here.

Take a solution basis consisting of the distinguished Apéry solution and a companion.  Write a projective family as
\[
F_i(q)=B_i-qC_i.
\]
After localizing at a nonzero companion coordinate \(C_0\), the evaluation ideal for several returns satisfies the elementary elimination identity
\[
(F_i(q):0\le i\le m)
=
(C_0q-B_0,\;B_iC_0-B_0C_i:1\le i\le m).
\tag{46}
\]
For the present three-return situation, \(m=2\).

Interpretation:

* the first generator is exactly the **marking pin** selecting the distinguished projective fiber;
* the remaining generators are exactly the ordinary projective-return minors.

When the first marked value is zero, those two minors reduce, up to local units, to the two return conditions
\[
N_h(a)=0,\qquad N_h(a+h)=0.
\]
There is no fourth equation hidden in the companion label or in differentiation with respect to the marked coordinate.

This is the three-return subcase of the current repository's terminal companion-label analysis for four returns.  It is another exact reason the Casoratian should not be expected to manufacture a new right side: it supplies a basis and proves transversality of the marked fiber, but the marked fiber contributes exactly the one condition \(b_a=0\) that is already present in \(D_{a,h}\).

---

# 10. What would actually be enough

There are three increasingly concrete ways to state the missing Apéry-specific input.

## 10.1 Logically minimal: quarter marked-return avoidance

For every prime \(p\) and positive \(a,h\) with
\[
p>4(a+2h),
\]
prove
\[
N_h(a)=N_h(a+h)=0\pmod p
\quad\Longrightarrow\quad
b_a\ne0\pmod p.
\tag{QMRA}
\]
By (1), this is exactly the statement that the lower-quarter Apéry zero set has no three-term arithmetic progression
\[
a,\ a+h,\ a+2h.
\]
This is an exact formulation, not progress by itself.

## 10.2 Smallest simple scalar sufficient theorem

Prove
\[
 P^+\!\left(\gcd(b_a,\mathcal R_h)\right)
 \le4(a+2h).
\tag{MSR}
\]
Because of (40), `(MSR)` implies `(QMRA)`.  It is stronger than necessary because it forgets which root of the separated resultant produced the prime, but it is a clean scalar interaction between the distinguished Apéry solution and the universal separated-return obstruction.

This is, in my view, the smallest genuinely Apéry-specific **scalar** theorem worth testing next.  It is not a consequence of the Casoratian, reflection, or the current continuant identities.

## 10.3 The desired explicit S-unit Bezout certificate

Let
\[
S_Q=\prod_{\ell\le4Q\atop \ell\text{ prime}}\ell.
\]
A direct certificate would be an explicit identity
\[
 \boxed{
 U_{a,h}B_a+V_{a,h}N_h(a)+W_{a,h}N_h(a+h)=R_{a,h},
 }
\tag{47}
\]
with
\[
U_{a,h},V_{a,h},W_{a,h}\in\mathbb Z[1/S_Q]
\]
and
\[
R_{a,h}\in\mathbb Z[1/S_Q]^\times.
\tag{48}
\]
After clearing denominators, this means \(R_{a,h}\) has no prime divisor above \(4Q\).  Such an identity immediately proves the quarter exclusion.

Abstract existence of (47) is just another form of the desired gcd statement, so it becomes useful only if the coefficients and the right side come from an **explicit Apéry-specific construction** whose support can be analyzed without first knowing \(D_{a,h}\).  The binomial sum, the diagonal/constant-term realization, a Frobenius/Cartier identity, or a new modular-form identity are plausible sources.  The two-dimensional transfer algebra is not.

---

# 11. Final answer to the specific proposed ingredients

### Casoratian

**Proved useful, but only for saturation/equivalence.**  It gives the exact determinant identities (30)--(31), which are the boundary specialization of the continuant Bezout law.  At a zero it proves the companion coordinate is a unit and converts \(N_h(a)=0\) exactly into the next Apéry zero.  It yields no nonzero controlled-support remainder.

### Reflection

**No general quarter exclusion.**  It supplies the reflected pair of roots.  It does settle \(h=1,2\), but for \(h\ge3\) the degree permits the four roots, and the supplied \(h=4\) witness shows that the separated double-return geometry really occurs in the quarter.  The endpoint Apéry factorization does not apply at the reflected positive-root locations because they lie outside its boundary strip.

### Known continuant identities

**They give a genuine fixed-\(h\) bound, not the required support theorem.**  The correct terminal scalar is the nonzero separated resultant \(\mathcal R_h\), with
\[
D_{a,h}\mid\mathcal R_h,
\qquad
\log|\mathcal R_h|\ll h^2\log h.
\]
This excludes fixed/sufficiently small \(h\), but the separated resultant has genuinely new large prime factors.  The \(8431\) witness is one such factor at \(h=4\).

### Explicit Bezout coupling

**The best one available from these ingredients is (41)**, whose right side is
\[
\gcd(b_a,\mathcal R_h).
\]
That right side has no known controlled prime support.  Replacing it by an \(S_Q\)-unit would be precisely the new Apéry-specific theorem required.

---

# 12. Recommended route decision

I would **stop trying to extract the quarter theorem by further algebraic recombination of Casoratian + reflection + continuants**.  The current repository's marked-row calculation and the explicit separated-return witness jointly give a rigorous reason: those tools describe the projective return geometry and its basis changes, while the missing fact is that the particular Apéry initial line avoids the phantom common-root locus.

The useful objects to retain are
\[
\mathcal R_h=\operatorname{Res}(N_h(x),N_h(x+h))
\]
and the marked scalar envelope
\[
G_{a,h}=\gcd(b_a,\mathcal R_h).
\]
A new attack should ask for arithmetic information on the **actual Apéry coefficient at a labeled separated-return root**, not for another universal continuant identity.  If an explicit binomial/constant-term/Frobenius identity can turn that marked condition into an \(S_Q\)-unit Bezout remainder, it would be genuinely new and would evade the solution-space countermodel above.

That is the smallest scope in which I see a theorem capable of ruling out the distinguished quarter triples without simply restating the target.