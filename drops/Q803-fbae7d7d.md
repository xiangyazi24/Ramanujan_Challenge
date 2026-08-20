ANSWER Q803 fbae7d7d

# Q803 audit: repaired rectangle theorem

## Verdict

The previous Q802 counterexample was invalid. Under the corrected definition of a rectangle, the proposed implication is essentially correct provided the Q8740 aligned-pair estimate is available with the stated uniformity. The old objection confused a missing vertical edge in the fixed-h graph with a missing rectangle unit. A rectangle needs four own hits and six cross-units; it does not require the vertical s-edges to belong to the fixed-h graph.

## 1. Matching/path extraction

Let the fixed-h fully transverse hit-edge graph have R edges. Since each vertex is a residue in the row and edges have the same difference h, every connected component is a path. Choose an alternating matching J. A path graph with R edges has a matching of size

\[
J\ge \lceil R/2\rceil.
\]

Write the matching edges as

\[
(u_i,u_i+h),\qquad i=1,\dots,J.
\]

Their starting points are distinct. They lie in an interval of length at most L in the quotient row coordinate. The ordered gaps between consecutive matching starts have total length at most L (up to the harmless endpoint convention). Hence the pigeonhole step gives many small gaps.

More precisely, among the J-1 gaps, at least (J-1)/2 satisfy

\[
\Delta_i\le \frac{2L}{J-1}.
\]

Because starting residues are congruent in the q-row, possible s values are multiples of q. The number of possible multiples of q in this range is at most

\[
\frac{2L}{q(J-1)}+1.
\]

Keeping the integer term gives a slightly weaker constant than Q802; for R>=6 it can be absorbed, yielding a valid choice with

\[
s\le \frac{6L}{R},\qquad q\mid s,
\]

and s != h after excluding the fixed matching displacement. The exclusion is legitimate because the matching edges are disjoint and a second edge with displacement h would be the same fixed-h graph direction, not the new rectangle direction.

## 2. Raw rectangles from the matching

For every pair of matching edges with common gap s,

\[
(u,u+h),\quad (u+s,u+s+h),
\]

we obtain the four own hits

\[
u,u+h,u+s,u+s+h.
\]

The vertical pairs

\[
(u,u+s),\quad (u+h,u+h+s)
\]

are not required to be fixed-h edges. They are simply additional hit pairs whose gaps must be checked.

Thus the raw count M satisfies

\[
M\ge \frac{q(J-1)^2}{4L}
\]

with the constants coming from the number of small matching gaps and the number of q-spaced starts. Since J>=R/2, this gives the claimed scale

\[
M\gg \frac{qR^2}{L}.
\]

## 3. Bad cross-unit analysis

A raw rectangle fails only if one of the four missing cross-units is not transverse. Each failure is a pair of hit vertices with one of the four gaps

\[
s,\quad h+s,\quad |h-s|.
\]

(The horizontal h-pairs are already transverse by construction.) Every such failed unit is therefore an aligned hit-pair of gap at most

\[
h+s.
\]

Hence the Q8740 estimate applies, provided it is stated for the whole row and all gaps up to h+s:

\[
\#\{(x,y):x,y\text{ hits}, |x-y|\le h+s,\text{ aligned failure}\}
=O((h+s)^2).
\]

This is exactly the needed bridge. No additional divisibility phenomenon escapes Q8740.

## 4. Pair-to-rectangle incidence

A fixed bad aligned pair can occur in only a bounded number of raw rectangles. Let the bad pair have difference d.

It can occupy a side of the rectangle only when d is one of

\[
s,\ h+s,\ |h-s|.
\]

Solving the placement equations gives at most:

* one placement for d=s;
* one placement for d=h+s;
* potentially two placements when d=|h-s|.

The exceptional overlap cases h=2s or s=2h merge two of these possibilities. Therefore the uniform incidence bound is

\[
\le 3.
\]

The factor 3 is sharp in the exceptional ratio cases; otherwise the factor is smaller.

Therefore deletion of all bad rectangles costs

\[
O((h+s)^2).
\]

## 5. Correct theorem statement

For R>=6, there exists a multiple s of q such that

\[
s\ne h,\qquad s\le 6L/R,
\]

and

\[
Q_{q;h,s}(n)
\ge
\frac{qR^2}{36L}-O_q((h+s)^2).
\]

The constants depend only on the harmless endpoint and incidence losses.

## 6. Consequence for LGTQ-2

If

\[
\max_{h\le C_1\log n,\ s\le C_2\log^2 n}
Q_{q;h,s}(n)=o(n/\log^4 n),
\]

then a fixed-h count R cannot reach the LGTQ-2 obstruction size. Indeed the extraction gives a rectangle with s in the allowed logarithmic range and with quadratic amplification in R. Hence the two-edge obstruction collapses, giving the desired all-index implication.

## 7. Remaining caveat

The theorem is not purely combinatorial: its validity depends on the exact quantifier form of Q8740. It must count all aligned hit pairs in the whole row uniformly for every gap up to C_1 log n+C_2 log^2 n. If Q8740 only controls fixed gaps or averages over gaps, the deletion step does not follow.

The earlier E_m construction is not a counterexample: it creates many raw rectangles precisely because two horizontal hit edges already supply the four own hits. The construction fails only by confusing graph edges with rectangle units.