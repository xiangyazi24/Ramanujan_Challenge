ANSWER Q5726 bebaee33

# State-prime descent in the first cell: exact one-step contraction, multiplicity ledger, and the fixed-point obstruction

## Executive verdict

I take the new ray-vanishing theorem in the question as proved and audited.
For a first-cell residue block

\[
I=I(M;R,H)=\{R,R+1,\ldots,R+H-1\},
\qquad R+H-1<M/2,
\]

write, for each prime \(p\),

\[
M=a_p p+s_p,
\qquad 0\le s_p<p.
\]

Define the theorem-generated state-prime set

\[
\mathcal S(M;R,H)=
\left\{
 p:
 \begin{array}{l}
 p\text{ prime},\quad p>R+H-1,\quad 2s_p<R,\\
 p\mid b_{a_p}b_{s_p}
 \end{array}
\right\}.
\tag{0.1}
\]

Then the new theorem gives the exact block divisibility

\[
\boxed{
\prod_{p\in\mathcal S(M;R,H)}p
\mid
\gcd_{r\in I} C_M(M-r).
}
\tag{0.2}
\]

The proposed geometric descent has a real positive part, but it does **not** close:

1. If \(p\mid b_{s_p}\), then \(p\) maps canonically to the ordinary top-half target
   \[
   (N,p,u)=(p+s_p,p,s_p).
   \]
   For quotient \(a_p\ge2\), the outer index contracts uniformly:
   \[
   N\le \frac58 M.
   \]
2. If \(p\mid b_{a_p}\), Apéry--Lucas descent in the base-\(p\) digits of \(a_p\) produces a digit \(u<p\) with \(p\mid b_u\).  The corresponding target row is \(N=p+u\).  For the large-prime branch this is at most \(M/2+O(\sqrt M)\); if the zero occurs in a higher base-\(p\) digit, then \(N\le2\sqrt M\).
3. These are only **one-step** contractions.  At the child row \(N=p+u\), the same prime has quotient one and residue \(u\).  The \(b_u\)-edge either becomes the literal self-loop
   \[
   (N,p,u)\longmapsto(N,p,u)
   \]
   or the theorem has no further applicable edge.  There is no iterated geometric chain for the same prime.
4. The multiplicities can be controlled exactly:
   - for fixed lower residue \(s\), the product of all parent primes divides
     \[
     \gcd(M-s,b_s);
     \]
   - for fixed quotient \(a\), the product of all \(b_a\)-channel primes divides
     \[
     \operatorname{rad}(b_a).
     \]
5. Consequently the complete theorem-generated radical divides the explicit integer
   \[
   \boxed{
   \left(\prod_{a\le A}\operatorname{rad}(b_a)\right)
   \left(\prod_{0\le s\le S}\operatorname{rad}\gcd(M-s,b_s)\right),
   }
   \tag{0.3}
   \]
   where
   \[
   A=\left\lfloor\frac{M}{R+H}\right\rfloor,
   \qquad
   S=\left\lfloor\frac{R-1}{2}\right\rfloor.
   \]
   Since \(b_n\le40^n\), this gives
   \[
   \boxed{
   \log\operatorname{rad}\mathcal S(M;R,H)
   \le
   \frac{A(A+1)}2\log40
   +\sum_{s=0}^{S}\min\{\log M,s\log40\}.
   }
   \tag{0.4}
   \]
   The second term is \(O(R\log M)\), not \(o(H)\).
6. The failure is not merely a resemblance to the original problem.  There is an exact actual-Apéry alignment theorem.  Every direct target prefix
   \[
   \mathcal T_M(S)=
   \prod_{\substack{0\le s\le S\\M-s\text{ prime}\\M-s\mid b_s}}(M-s)
   \]
   divides one state block as soon as \(2S<R\):
   \[
   \boxed{
   \mathcal T_M(S)\mid\gcd_{r=R}^{R+H-1}C_M(M-r).
   }
   \tag{0.5}
   \]
   The descent map for every prime in (0.5) fixes the outer row \(N=M\).  Thus an \(o(H)\) theorem for the state radical would imply the much stronger bound \(\log\mathcal T_M(S)=o(H)\) for an actual direct-target prefix.
7. On the natural near-boundary block with \(H=M^{1/3}\), there is a second exact alignment.  Quotient-two primes in the interval
   \[
   M/2-H<p\le M/2,
   \qquad s=M-2p<2H,
   \]
   all divide the whole block whenever \(p\mid b_2b_s\).  Their universal prime support has logarithmic weight \(O(H)\), and this scale is sharp for prime intervals along subsequences.  The descent sends them to the diagonal family of distinct target rows
   \[
   N=p+s=M-p,
   \]
   spanning another interval of length \(H\).  It supplies no little-\(o\) saving.

Therefore the answer is negative:

\[
\boxed{
\text{the new implication yields a one-step reduction, but not an iterable geometric descent proving }o(H).
}
\]

The narrow missing arithmetic statements are now explicit:

- an upper-quarter direct-prefix theorem for the \(a=1\) fixed-point branch;
- an \(a=2\) diagonal anti-clustering theorem
  \[
  \sum_{M/2-H<p\le M/2\atop p\mid b_{M-2p}}\log p=o(H).
  \]

Without at least these, the theorem-generated actual-state radical can retain the full block scale.

---

## 1. Exact block consequence of the new ray theorem

Let

\[
d_r=M-r,
\qquad r\in I(M;R,H).
\]

For \(p\in\mathcal S(M;R,H)\), one has for every \(r\in I\)

\[
2s_p<R\le r<p,
\qquad d_r=M-r>M/2.
\]

The theorem in the question says that every nonzero ray coefficient

\[
c_M(d_r\kappa),\qquad \kappa\ne0,
\]

vanishes modulo \(p\).  The central coefficient is \(b_M\), and Apéry--Lucas gives

\[
b_M\equiv b_{a_p}b_{s_p}\pmod p.
\]

Hence \(p\mid C_M(d_r)\) for every \(r\in I\), proving (0.2).

This is a block theorem, not merely a pointwise statement.  Once the start threshold \(2s_p+1\) lies to the left of the block and the endpoint lies below \(p\), the same prime persists across the entire block.

### Plain ASCII

```text
I(M;R,H) = {R,...,R+H-1},  R+H-1 < M/2
M = a_p*p + s_p, 0 <= s_p < p

state condition:
    p > R+H-1
    2*s_p < R
    p | b[a_p]*b[s_p]

then:
    p | gcd_{r in I} C_M(M-r)
```

---

## 2. The exact lower-residue descent maps

### 2.1 The \(b_s\)-channel

Assume

\[
p\mid b_s,
\qquad M=ap+s,
\qquad 2s<R<p.
\]

Put

\[
\Phi_s(M,p)=(N,p,u),
\qquad N=p+s,\quad u=s.
\tag{2.1}
\]

Because \(p>2s\), one has \(p>N/2\).  Thus \((N,p,u)\) is an ordinary top-half target pair:

\[
N=p+u,
\qquad p\mid b_u.
\tag{2.2}
\]

The relation with the parent quotient is exact:

\[
\boxed{
N=p+s=\frac{M+(a-1)s}{a}=M-(a-1)p.
}
\tag{2.3}
\]

If \(a\ge2\), then \(s<R/2<M/4\), so

\[
\frac NM
<\frac1a+\frac{a-1}{4a}
=\frac{a+3}{4a}
\le\frac58.
\tag{2.4}
\]

Thus the outer size contracts by a fixed factor on every \(a\ge2\) edge.

However, the local width does not contract.  For fixed \(a\), admissible residues satisfy

\[
s\equiv M\pmod a.
\]

Increasing \(s\) by \(a\) decreases \(p=(M-s)/a\) by one and changes the child row by

\[
N=M-(a-1)p
\quad\Longrightarrow\quad
N\mapsto N+(a-1).
\tag{2.5}
\]

So an interval of parent primes of width \(W\) becomes a child-row progression of span approximately \((a-1)W\).  The outer scale contracts, but the family does not collapse to one child row unless \(a=1\).

### 2.2 The \(b_a\)-channel

Assume \(p\mid b_a\).  Write the base-\(p\) expansion

\[
a=d_0+d_1p+\cdots+d_ep^e,
\qquad0\le d_i<p.
\]

Apéry--Lucas gives

\[
b_a\equiv\prod_{i=0}^e b_{d_i}\pmod p.
\]

Hence there is at least one digit \(u=d_i\) with

\[
p\mid b_u.
\tag{2.6}
\]

Choose the least such digit index to make the map canonical, and put

\[
\Phi_a(M,p)=(p+u,p,u).
\tag{2.7}
\]

There are two cases.

#### Case A: \(a<p\)

Then \(u=a\), and

\[
N_a=p+a\le \frac Ma+a.
\tag{2.8}
\]

For \(a\ge2\), one has \(a<p\) and therefore \(a<\sqrt M\), giving

\[
N_a\le M/2+\sqrt M.
\tag{2.9}
\]

The case \(a=1\) has \(b_1=5\), so its only prime is \(5\); it is absent from every asymptotic block with \(R\ge5\).

#### Case B: \(a\ge p\)

Then \(p^2\le ap\le M\), so \(p\le\sqrt M\), and

\[
N_a=p+u<2p\le2\sqrt M.
\tag{2.10}
\]

Thus the higher-digit branch contracts even more strongly, but it is confined to primes at most \(\sqrt M\).

---

## 3. Why the descent has depth one

The \(b_s\)-child in (2.1) satisfies

\[
N=p+s=1\cdot p+s.
\]

At that child row the same prime has quotient one and the same residue \(s\).  If one chooses another first-cell segment with residue coordinate \(r'\) satisfying

\[
2s<r'<p,
\]

the new theorem sends the prime to

\[
p+s=N
\]

again.  Thus the descent edge is the exact self-loop

\[
\boxed{
(N,p,s)\longmapsto(N,p,s).
}
\tag{3.1}
\]

For a \(b_a\)-child, if \(p>2u\) the same self-loop occurs.  If \(p\le2u\), the theorem does not provide a further edge.  In neither case is there a second guaranteed contraction.

This is the fundamental graph-theoretic obstruction:

```text
parent state prime at M, quotient a>=2
        |
        | one contraction
        v
ordinary target (N=p+u, p | b_u)
        |
        +--> self-loop if p>2u
        +--> no theorem edge otherwise
```

One cannot iterate the same implication down a geometric tower.  Factoring the integer \(b_u\) introduces other primes, but it does not produce a relation that replaces the original prime \(p\) by those primes while preserving the state-radical divisibility being bounded.

---

## 4. Exact multiplicity control

The descent does give clean multiplicity statements.

### 4.1 Fixed lower residue

Fix \(s\).  Every parent prime with residue \(s\) satisfies

\[
p\mid M-s
\]

because \(M-s=ap\).  On the \(b_s\)-channel it also satisfies \(p\mid b_s\).  Therefore

\[
\boxed{
\prod_{\substack{p\in\mathcal S(M;R,H)\\s_p=s\\p\mid b_s}}p
\mid\operatorname{rad}\gcd(M-s,b_s).
}
\tag{4.1}
\]

This controls all quotient values \(a\) at once.  In particular,

\[
\sum_{\substack{p\in\mathcal S\\s_p=s\\p\mid b_s}}\log p
\le\min\{\log(M-s),\log b_s\}.
\tag{4.2}
\]

### 4.2 Fixed quotient

Fix \(a\).  Every prime on the \(b_a\)-channel divides the same integer \(b_a\).  Hence

\[
\boxed{
\prod_{\substack{p\in\mathcal S(M;R,H)\\a_p=a\\p\mid b_a}}p
\mid\operatorname{rad}(b_a).
}
\tag{4.3}
\]

The residues for this quotient obey \(s\equiv M\pmod a\), so among

\[
0\le s\le S=\left\lfloor\frac{R-1}{2}\right\rfloor
\]

there are at most

\[
1+\left\lfloor\frac Sa\right\rfloor
\tag{4.4}
\]

possible parent primes.

### 4.3 Uniqueness

For fixed \(M\) and \(p\), the pair

\[
(a_p,s_p)=\left(\left\lfloor\frac Mp\right\rfloor,M-p\left\lfloor\frac Mp\right\rfloor\right)
\]

is unique.  Thus a state prime is not repeated across quotient classes.  If it belongs to both channels, choosing the \(b_s\)-channel first and the \(b_a\)-channel second gives a canonical single descent edge.

### 4.4 Multiplicity across a short \(M\)-block

Suppose \(M\) itself ranges over an interval of length \(B<p\).

- A fixed \((p,s)\) can occur for at most one \(M\), because the allowed values are \(M=ap+s\), spaced by \(p\).
- A fixed \((p,a)\) can occur for up to \(B\) consecutive residues \(s\), but all of them charge to the single target divisor \(p\mid b_a\).

This explains the observed run

```text
p=73, a=2, M=146,147,148,149,150
s=0,1,2,3,4
zero starts r=1,3,5,7,9
```

without paying five copies of \(\log73\) in a union radical.

---

## 5. The exact height ledger

Put

\[
A=\left\lfloor\frac{M}{R+H}\right\rfloor,
\qquad
S=\left\lfloor\frac{R-1}{2}\right\rfloor.
\tag{5.1}
\]

Every state prime has \(a_p\le A\) and \(s_p\le S\).  Equations (4.1) and (4.3) prove the integer divisibility (0.3):

\[
\prod_{p\in\mathcal S(M;R,H)}p
\mid
\left(\prod_{a=1}^{A}\operatorname{rad}(b_a)\right)
\left(\prod_{s=0}^{S}\operatorname{rad}\gcd(M-s,b_s)\right).
\tag{5.2}
\]

Because \(\Lambda\) has nonnegative coefficients and \(\Lambda(1,1,1)=40\),

\[
b_n\le40^n.
\tag{5.3}
\]

Therefore

\[
\log\prod_{a=1}^{A}\operatorname{rad}(b_a)
\le\frac{A(A+1)}2\log40,
\tag{5.4}
\]

and

\[
\log\operatorname{rad}\gcd(M-s,b_s)
\le\min\{\log M,s\log40\}.
\tag{5.5}
\]

Summing proves (0.4).

For \(R\gg\log M\),

\[
\sum_{s=0}^{S}\min\{\log M,s\log40\}
=O(R\log M).
\tag{5.6}
\]

This is much larger than \(H=M^{1/3}\) on the blocks of interest.

### 5.1 Exact candidate-prime intervals

The conditions

\[
p\ge R+H,
\qquad 0\le s=M-ap\le S
\]

put every quotient-\(a\) state prime in the interval

\[
\boxed{
J_a(M;R,H)=
\left[
\max\left\{R+H,\left\lceil\frac{M-S}{a}\right\rceil\right\},
\left\lfloor\frac Ma\right\rfloor
\right].
}
\tag{5.7}
\]

Let \(\vartheta(x)=\sum_{p\le x}\log p\).  Independently of any Apéry zero condition,

\[
\log\prod_{p\in\mathcal S(M;R,H)}p
\le
\sum_{a=1}^{A}
\left(
\vartheta(\lfloor M/a\rfloor)
-\vartheta(\min J_a-1)
\right).
\tag{5.8}
\]

This interval carrier is often sharper than (5.2), but it still has the wrong scale.

---

## 6. The exact \(a=1\) fixed-point alignment

Define the direct target prefix

\[
\mathcal T_M(S)=
\prod_{\substack{0\le s\le S\\p=M-s\text{ prime}\\p\mid b_s}}p.
\tag{6.1}
\]

### Theorem 6.1

Assume

\[
2S<R,
\qquad R+H-1<M/2.
\]

Then

\[
\boxed{
\mathcal T_M(S)\mid
\gcd_{r=R}^{R+H-1}C_M(M-r).
}
\tag{6.2}
\]

#### Proof

For a factor \(p=M-s\) in (6.1), one has

\[
M=1\cdot p+s,
\qquad p\mid b_1b_s
\]

because \(b_1=5\) and \(p\mid b_s\).  Also

\[
2s\le2S<R.
\]

Since \(S<R/2<M/4\),

\[
p=M-s>3M/4>R+H-1.
\]

The block theorem (0.2) applies.  Multiplying the distinct primes proves (6.2). \(\square\)

This is the promised rigorous adversarial alignment.  It uses the actual Apéry target set, not an arbitrary recurrence model.

The descent map is

\[
N=p+s=M.
\tag{6.3}
\]

Thus every prime in the whole prefix maps back to the same outer row.  There is no decrease at all.

### Near-boundary consequence

Take

\[
H=\lceil M^{1/3}\rceil,
\qquad R=\lfloor M/2\rfloor-2H.
\tag{6.4}
\]

For large \(M\), this is a valid first-cell block.  Equation (6.2) contains all direct targets with

\[
s\le M/4-H+O(1),
\]

equivalently

\[
p=M-s\ge3M/4+H+O(1).
\]

Hence

\[
\boxed{
\prod_{\substack{3M/4+H+O(1)\le p\le M\\p\mid b_{M-p}}}p
\mid
\gcd_{r=R}^{R+H-1}C_M(M-r).
}
\tag{6.5}
\]

Therefore a theorem

\[
\log\gcd_{r=R}^{R+H-1}C_M(M-r)=o(H)
\]

would imply the much stronger actual-Apéry statement

\[
\sum_{\substack{3M/4+H+O(1)\le p\le M\\p\mid b_{M-p}}}\log p=o(M^{1/3}).
\tag{6.6}
\]

This is an exact reduction, not the assertion that the two problems merely look similar.

---

## 7. The quotient-two diagonal obstruction at the block scale

Continue with the near-boundary block (6.4).  Consider primes

\[
M/2-H<p\le M/2
\tag{7.1}
\]

and put

\[
s=M-2p.
\tag{7.2}
\]

Then \(0\le s<2H\) and \(M=2p+s\).  For sufficiently large \(M\),

\[
2s<4H<R,
\qquad p>R+H-1.
\]

Thus every prime satisfying

\[
p\mid b_2b_s
\tag{7.3}
\]

divides the entire block.  Define

\[
\mathcal T_M^{(2)}(H)=
\prod_{\substack{M/2-H<p\le M/2\\p\mid b_2b_{M-2p}}}p.
\tag{7.4}
\]

Then

\[
\boxed{
\mathcal T_M^{(2)}(H)\mid
\gcd_{r=R}^{R+H-1}C_M(M-r).
}
\tag{7.5}
\]

Since \(b_2=73\), all asymptotically large factors in (7.4) are in the \(b_s\)-channel.

The descent map is

\[
N=p+s=M-p.
\tag{7.6}
\]

As \(p\) runs through the interval (7.1), the child rows \(N\) run through an interval of the same length \(H\), in the opposite direction.  Thus the quotient-two descent contracts the outer scale from \(M\) to approximately \(M/2\), but it does **not** compress the block width or combine the targets into one fixed row.

The universal prime support in (7.1) has logarithmic weight \(O(H)\) by Brun--Titchmarsh.  This order is sharp along subsequences at the level of candidate primes: partitioning a dyadic prime interval into subintervals of length comparable to \(x^{1/3}\) and using the prime number theorem shows that some such intervals have logarithmic prime weight \(\gg H\).

Therefore the geometry alone naturally pays the whole \(H\)-budget.  To obtain little-\(o\), one needs the genuinely arithmetic diagonal estimate

\[
\boxed{
\sum_{\substack{M/2-H<p\le M/2\\p\mid b_{M-2p}}}\log p=o(H).
}
\tag{7.7}
\]

The new ray theorem does not imply (7.7); it transforms every prime counted by (7.7) into a block factor.

---

## 8. Why induction on the contracted outer index still fails

Suppose one tries to use an inductive target bound at the child indices.

### 8.1 The \(a=1\) branch

All descendants have the same child row \(N=M\), by (6.3).  The induction parameter does not decrease.

### 8.2 A fixed \(a\ge2\) branch

The child rows are

\[
N=M-(a-1)p
\]

and form an arithmetic progression of step \(a-1\).  The number of admissible residues is at most

\[
1+\left\lfloor\frac Sa\right\rfloor.
\]

Even if every child row of size at most \(5M/8\) satisfied a pointwise bound \(o(N^{1/3})\), summing that bound over \(\asymp S/a\) different rows is not \(o(M^{1/3})\).  In the quotient-two near-boundary case there are potentially \(\asymp H\) distinct child rows in an interval of length \(H\).

A fixed-row target theorem cannot be summed over that diagonal family without a new cross-row anti-clustering result.

### 8.3 Recursion on the residue itself

One may instead try to charge \(p\) to the height of \(b_s\), since \(s<R/2\).  For a single residue this is geometric.  Across all residues, however,

\[
\sum_{s\le S}\log b_s\le \log40\sum_{s\le S}s=\Theta(S^2),
\]

and the sharper fixed-residue divisor bound gives only

\[
\sum_{s\le S}\log\gcd(M-s,b_s)=O(S\log M).
\]

Neither is \(o(H)\) when \(S\) is a positive power or positive fraction of \(M\).

The residue halves on each individual edge, but the number of terminal residues is too large, and the same-prime graph has no second edge.

---

## 9. The verified examples in the descent graph

### 9.1 The \(p=73\) family

For

\[
M=146+s,
\qquad0\le s\le4,
\]

one has

\[
M=2\cdot73+s,
\qquad73\mid b_2.
\]

The theorem gives zero start

\[
r_0=2s+1,
\]

namely

```text
M=146  s=0  start=1
M=147  s=1  start=3
M=148  s=2  start=5
M=149  s=3  start=7
M=150  s=4  start=9
```

All five parent states use the same \(b_a\)-child

\[
N=73+2=75,
\qquad73\mid b_2.
\]

Thus the repetition across \(M\) is genuinely compressible to one target prime, but the target row \(75\) is terminal.

### 9.2 The \(p=61\) example

For

\[
M=126=2\cdot61+4,
\]

one has \(61\mid b_4\).  The zero start is

\[
r_0=2\cdot4+1=9.
\]

This is the \(b_s\)-channel, and the child target is

\[
N=61+4=65.
\]

At row \(65\), the same prime has quotient one and residue four, so the descent is the self-loop \((65,61,4)\mapsto(65,61,4)\).

---

## 10. Precise no-go statement for this descent mechanism

Let a **descent-only proof** mean a proof using only:

1. the new ray-vanishing implication;
2. Apéry--Lucas factorization of \(b_M\) and \(b_a\);
3. the maps \(p\mapsto(p+s,p,s)\) and \(p\mapsto(p+u,p,u)\);
4. pointwise size bounds \(b_n\le40^n\);
5. gcd/radical/lcm operations and prime-interval bounds;
6. induction solely on the child outer index.

### Theorem 10.1

Such a proof cannot establish

\[
\log\prod_{p\in\mathcal S(M;R,H)}p=o(H)
\]

uniformly for all first-cell blocks of length \(H=M^{1/3}\), unless it proves at least one new arithmetic estimate implying both:

\[
\sum_{\substack{3M/4+H+O(1)\le p\le M\\p\mid b_{M-p}}}\log p=o(H),
\tag{10.1}
\]

and

\[
\sum_{\substack{M/2-H<p\le M/2\\p\mid b_{M-2p}}}\log p=o(H).
\tag{10.2}
\]

#### Proof

Equation (6.5) embeds the complete product in (10.1) into one state block, and every descent edge fixes the outer row \(M\).  Equation (7.5) embeds the complete product in (10.2) into the same type of block; its descendants occupy a diagonal family of distinct rows of total width \(H\), and no second descent edge reduces those targets.  The height and interval operations listed above give only \(O(M)\) for (10.1) and \(O(H)\) for (10.2).  Therefore little-\(o\) requires an additional arithmetic theorem that removes a little-\(o\) proportion of the prime weight in these actual Apéry target sets. \(\square\)

This is a demonstrated obstruction for the distinguished shell, not a hostile arbitrary sequence.

---

## 11. What positive content survives

The new theorem is still valuable.

1. **Classification.**  It identifies a large, previously mysterious part of the actual-state gcd as digit-zero inheritance from \(b_a b_s\).
2. **One-step contraction.**  Every \(a\ge2\) \(b_s\)-state prime maps to an ordinary target row at most \(5M/8\); every high \(b_a\)-state prime maps to at most \(M/2+O(\sqrt M)\).
3. **Exact multiplicity compression.**  Repeated state appearances of one \(b_a\)-prime across many consecutive \(M\)-values cost one prime, and all quotient multiplicity at fixed \(s\) is absorbed by \(\gcd(M-s,b_s)\).
4. **Sharp missing statements.**  The obstruction is no longer an unspecified state-content factor.  In the near-boundary first cell it is concentrated in the direct-prefix family (10.1) and the quotient-two diagonal family (10.2).

A successful continuation should therefore attack (10.2) first: it has exactly the desired \(H\)-scale, a one-dimensional affine relation

\[
s=M-2p,
\]

and no broad \(M/4\)-length prefix.  Any power saving or little-\(o\) theorem for this diagonal family would genuinely improve the state-prime ledger.  The direct fixed-point family (10.1) remains the larger obstruction.

---

## 12. Dependency-free audit code

The following standard-library script verifies the descent maps, channel multiplicities, and the displayed \(73\) and \(61\) examples.  It uses Apéry--Lucas recursively, with a direct binomial sum only for one base-\(p\) digit.

```python
#!/usr/bin/env python3
from math import comb, isqrt


def primes_upto(n):
    if n < 2:
        return []
    s = bytearray(b"\x01") * (n + 1)
    s[:2] = b"\x00\x00"
    for p in range(2, isqrt(n) + 1):
        if s[p]:
            s[p*p:n+1:p] = b"\x00" * (((n-p*p)//p) + 1)
    return [p for p in range(2, n + 1) if s[p]]


def apery_digit_mod(n, p):
    assert 0 <= n < p
    return sum((comb(n, k) * comb(n + k, k))**2
               for k in range(n + 1)) % p


def apery_mod(n, p):
    # Apéry--Lucas recursion.
    if n < p:
        return apery_digit_mod(n, p)
    q, r = divmod(n, p)
    return apery_mod(q, p) * apery_digit_mod(r, p) % p


def zero_digit(n, p):
    """Return a base-p digit u with p|b_u, assuming p|b_n."""
    while n:
        n, u = divmod(n, p)
        if apery_digit_mod(u, p) == 0:
            return u
    raise AssertionError("no zero digit")


def state_edges(M, R, H):
    out = []
    for p in primes_upto(M):
        a, s = divmod(M, p)
        if not (p > R + H - 1 and 2*s < R):
            continue
        za = apery_mod(a, p) == 0
        zs = apery_digit_mod(s, p) == 0
        if not (za or zs):
            continue
        if zs:
            u = s
            N = p + s
            channel = "s"
            assert p > 2*u
            assert N == (M + (a - 1)*s) // a
            if a >= 2:
                assert 8*N < 5*M + 8  # integer-safe N <= 5M/8 + O(1)
        else:
            u = zero_digit(a, p)
            N = p + u
            channel = "a"
            assert apery_digit_mod(u, p) == 0
        out.append((p, a, s, channel, u, N))
    return out


def audit_examples():
    p = 73
    assert apery_mod(2, p) == 0
    for M in range(146, 151):
        a, s = divmod(M, p)
        assert a == 2
        assert apery_mod(a, p) * apery_digit_mod(s, p) % p == 0
        assert 2*s + 1 == 2*(M - 146) + 1
        print("p73", "M", M, "a", a, "s", s,
              "zero_start", 2*s + 1, "child", p + a)

    p, M = 61, 126
    a, s = divmod(M, p)
    assert (a, s) == (2, 4)
    assert apery_digit_mod(s, p) == 0
    print("p61", "M", M, "a", a, "s", s,
          "zero_start", 2*s + 1, "child", p + s)


def audit_fixed_s_product(M, R, H):
    by_s = {}
    for p, a, s, channel, u, N in state_edges(M, R, H):
        if channel == "s":
            by_s.setdefault(s, []).append(p)
    for s, ps in by_s.items():
        z = 1
        for p in ps:
            z *= p
        assert (M - s) % z == 0
        assert apery_mod(s, p) == 0
        print("fixed_s", s, "primes", ps, "product_divides", M-s)


if __name__ == "__main__":
    audit_examples()
    for M in (146, 147, 148, 149, 150, 126, 320, 754):
        H = max(1, round(M**(1/3)))
        R = max(1, M//2 - 2*H)
        if R + H - 1 < M/2:
            print("STATE", M, R, H, state_edges(M, R, H))
            audit_fixed_s_product(M, R, H)
```

---

## Final conclusion

The new ray theorem turns the state-prime problem into a highly structured digit-zero problem.  It gives exact one-step contraction and excellent multiplicity bookkeeping.  But the contraction is not iterable: the lower-residue child is an ordinary quotient-one target and is then a fixed point.

More decisively, the actual direct target prefix and the quotient-two diagonal family embed exactly into one first-cell block.  Their natural logarithmic scales are respectively \(O(M)\) and \(O(H)\).  Therefore the theorem by itself cannot prove

\[
\log\operatorname{rad}(\text{state factors})=o(H).
\]

The smallest useful next theorem is the quotient-two diagonal estimate (10.2); after that, the upper-quarter direct fixed-point family (10.1) remains.