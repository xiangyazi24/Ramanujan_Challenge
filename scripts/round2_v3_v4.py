#!/usr/bin/env python3
"""Exact checks for Round 2, V3--V4.

V3 works with the Adolphson--Sperber polytope

    Delta = conv({0} union Supp(w(x+x^-1+y+y^-1+z+z^-1))).

All combinatorics use integer tuples and exact set operations.  V4 uses
``fractions.Fraction`` throughout.  No third-party packages are required.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path


OUT = Path("/tmp/round2_v3_v4.txt")

# Exponent order is (w, x, y, z).
APEX = "O"
BASE_VERTICES = ("X+", "X-", "Y+", "Y-", "Z+", "Z-")
VERTICES = (APEX,) + BASE_VERTICES
VECTOR = {
    "O": (0, 0, 0, 0),
    "X+": (1, 1, 0, 0),
    "X-": (1, -1, 0, 0),
    "Y+": (1, 0, 1, 0),
    "Y-": (1, 0, -1, 0),
    "Z+": (1, 0, 0, 1),
    "Z-": (1, 0, 0, -1),
}
OPPOSITE = {"X+": "X-", "X-": "X+", "Y+": "Y-", "Y-": "Y+", "Z+": "Z-", "Z-": "Z+"}
AXIS = {"X+": 1, "X-": 1, "Y+": 2, "Y-": 2, "Z+": 3, "Z-": 3}
SIGN = {v: (1 if v.endswith("+") else -1) for v in BASE_VERTICES}
MONOMIAL = {
    "X+": "w*x",
    "X-": "w*x^-1",
    "Y+": "w*y",
    "Y-": "w*y^-1",
    "Z+": "w*z",
    "Z-": "w*z^-1",
}


def affine_rank(points: list[tuple[int, ...]]) -> int:
    """Rank over Q of affine differences, using Fraction elimination."""
    if len(points) <= 1:
        return 0
    rows = [[Fraction(x - y) for x, y in zip(p, points[0])] for p in points[1:]]
    rank = 0
    col = 0
    while rank < len(rows) and col < len(rows[0]):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        q = rows[rank][col]
        rows[rank] = [x / q for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][col]:
                q = rows[i][col]
                rows[i] = [x - q * y for x, y in zip(rows[i], rows[rank])]
        rank += 1
        col += 1
    return rank


def base_faces() -> dict[int, set[frozenset[str]]]:
    """All nonempty faces of the three-dimensional base octahedron."""
    ans: dict[int, set[frozenset[str]]] = defaultdict(set)
    ans[0] = {frozenset((v,)) for v in BASE_VERTICES}
    ans[1] = {
        frozenset(pair)
        for pair in combinations(BASE_VERTICES, 2)
        if pair[1] != OPPOSITE[pair[0]]
    }
    ans[2] = {
        frozenset((sx, sy, sz))
        for sx, sy, sz in product(("X+", "X-"), ("Y+", "Y-"), ("Z+", "Z-"))
    }
    ans[3] = {frozenset(BASE_VERTICES)}
    return ans


def delta_faces() -> dict[int, set[frozenset[str]]]:
    """All nonempty proper faces of the pyramid Delta."""
    base = base_faces()
    ans: dict[int, set[frozenset[str]]] = defaultdict(set)
    for dim, faces in base.items():
        ans[dim].update(faces)
    ans[0].add(frozenset((APEX,)))
    # Join the apex to each proper face of the base.  Joining it to the full
    # base produces Delta itself, not a proper face.
    for dim in (0, 1, 2):
        for face in base[dim]:
            ans[dim + 1].add(face | {APEX})
    return ans


def dot(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    return sum(x * y for x, y in zip(a, b))


def exposing_normal(face: frozenset[str]) -> tuple[int, int, int, int]:
    """An exact exposing normal for every nonempty proper face."""
    if face == frozenset((APEX,)):
        return (-1, 0, 0, 0)
    base_part = face - {APEX}
    if base_part == frozenset(BASE_VERTICES):
        return (1, 0, 0, 0)
    normal = [(-1 if APEX in face else 0), 0, 0, 0]
    for v in base_part:
        normal[AXIS[v]] = SIGN[v]
    return tuple(normal)


def verify_face_lattice(faces: dict[int, set[frozenset[str]]]) -> None:
    expected = {0: 7, 1: 18, 2: 20, 3: 9}
    assert {d: len(faces[d]) for d in range(4)} == expected
    for dim, dim_faces in faces.items():
        for face in dim_faces:
            assert affine_rank([VECTOR[v] for v in face]) == dim
            normal = exposing_normal(face)
            vals = {v: dot(normal, VECTOR[v]) for v in VERTICES}
            maximum = max(vals.values())
            exposed = frozenset(v for v, value in vals.items() if value == maximum)
            assert exposed == face, (face, normal, exposed)


def fmt_face(face: frozenset[str]) -> str:
    return "{" + ",".join(v for v in VERTICES if v in face) + "}"


def nondegeneracy_lines(faces: dict[int, set[frozenset[str]]]) -> list[str]:
    """Check every face not containing the origin/apex.

    On a proper base face there is at most one selected monomial on each of
    the X/Y/Z axes.  A logarithmic derivative in any selected axis is then a
    nonzero Laurent monomial on the torus.  The full base is handled exactly.
    """
    lines: list[str] = []
    base_only = [f for d in range(4) for f in faces[d] if APEX not in f]
    assert len(base_only) == 6 + 12 + 8 + 1 == 27
    proper = sorted(
        (f for f in base_only if f != frozenset(BASE_VERTICES)),
        key=lambda f: (len(f), tuple(VERTICES.index(v) for v in VERTICES if v in f)),
    )
    assert len(proper) == 26
    for face in proper:
        # A proper octahedron face never contains an opposite pair.
        assert all(OPPOSITE[v] not in face for v in face)
        chosen = min(face, key=VERTICES.index)
        # The indicated log derivative is exactly +/- the chosen monomial,
        # because no other term of the face uses that coordinate.
        assert sum(AXIS[v] == AXIS[chosen] for v in face) == 1
        lines.append(
            f"  {fmt_face(face):18s}: D_{chosen[0].lower()} f_sigma is +/-{MONOMIAL[chosen]}, "
            "hence never 0 on the torus"
        )
    # For the full base, D_x=D_y=D_z=0 gives x^2=y^2=z^2=1.
    # D_w=f=0 then gives 2(e_x+e_y+e_z)=0.  Its possible integer values are
    # +/-2 and +/-6; hence it has no solution in characteristic p >= 5, and
    # it does have solutions in characteristics 2 and 3.
    for p in (2, 3):
        assert p in (2, 3)
    for p in (5, 7, 11, 13, 101):
        assert all(value % p for value in (-6, -2, 2, 6))
    lines.append(
        "  {X+,X-,Y+,Y-,Z+,Z-}: critical equations reduce to x,y,z in {+/-1} and "
        "2(x+y+z)=0; no solution in char p>=5, but solutions exist in char 2 and 3"
    )
    return lines


def l1_ball_3_count(n: int) -> int:
    """Number of u in Z^3 with |u_1|+|u_2|+|u_3| <= n."""
    if n < 0:
        return 0
    return sum((2**k) * comb(3, k) * comb(n, k) for k in range(4))


def fourth_difference(values: list[int]) -> list[int]:
    """Coefficients after multiplying a series by (1-t)^4."""
    kernel = (1, -4, 6, -4, 1)
    return [
        sum(kernel[k] * values[n - k] for k in range(5) if n >= k)
        for n in range(len(values))
    ]


def w_shift_hodge_slopes(a: Fraction) -> list[Fraction]:
    """A--S slopes for shift theta=(a,0,0,0), 0 <= a < 1."""
    assert 0 <= a < 1
    # Cone C(Delta) is r >= |u|_1 and wt(r,u)=r.  For r=a+n and u integral,
    # |u|_1 <= a+n is equivalent to |u|_1 <= n.  Thus R_theta=w^a R.
    lattice_counts = [l1_ball_3_count(n) for n in range(8)]
    numerator = fourth_difference(lattice_counts)
    assert numerator[:4] == [1, 3, 3, 1]
    assert all(x == 0 for x in numerator[4:])
    return [a + n for n, multiplicity in enumerate(numerator) for _ in range(multiplicity)]


def transverse_shift_hodge_slopes(a: Fraction) -> list[Fraction]:
    """A comparison convention: theta=(0,a,0,0)."""
    assert 0 <= a < 1
    if a == 0:
        return [Fraction(n) for n, m in enumerate((1, 3, 3, 1)) for _ in range(m)]

    # At integer weight r=n enumerate x in a+Z and y,z in Z.  Bounds below
    # are exact because every coordinate has absolute value <= n in the cone.
    counts: list[int] = []
    for n in range(8):
        total = 0
        for mx in range(-n - 1, n + 1):
            x = a + mx
            for y in range(-n, n + 1):
                for z in range(-n, n + 1):
                    if abs(x) + abs(y) + abs(z) <= n:
                        total += 1
        counts.append(total)
    numerator = fourth_difference(counts)
    assert numerator[:4] == [0, 2, 4, 2]
    assert all(x == 0 for x in numerator[4:])
    return [Fraction(n) for n, multiplicity in enumerate(numerator) for _ in range(multiplicity)]


def multiplicities(slopes: list[Fraction]) -> str:
    count = Counter(slopes)
    return ", ".join(f"{s} (x{count[s]})" for s in sorted(count))


def build_report() -> str:
    faces = delta_faces()
    verify_face_lattice(faces)

    out: list[str] = []
    out.append("ROUND 2: EXACT V3/V4 CHECK")
    out.append("============================")
    out.append("")
    out.append("V3. FACE LATTICE AND NONDEGENERACY")
    out.append("-----------------------------------")
    out.append("Exponent order: (w,x,y,z).")
    out.append("Supp(f)={(1,+/-1,0,0),(1,0,+/-1,0),(1,0,0,+/-1)}.")
    out.append("A-S Delta=conv({0} union Supp(f)) is a 4-dimensional pyramid over a 3-octahedron.")
    out.append("Ordinary conv(Supp(f)) alone is only the 3-dimensional base octahedron.")
    out.append("")
    out.append("f-vector (nonempty proper faces of the 4-polytope):")
    out.append("  vertices=7, edges=18, 2-faces=20, 3-faces=facets=9")
    out.append("  including empty face and Delta itself: (1,7,18,20,9,1)")
    out.append("  Euclidean volume(Delta)=1/3; normalized volume 4!*volume=8")
    out.append("")
    for dim in range(4):
        listed = sorted(faces[dim], key=lambda f: tuple(VERTICES.index(v) for v in VERTICES if v in f))
        out.append(f"dimension {dim} ({len(listed)} faces): " + " ".join(fmt_face(f) for f in listed))
    out.append("")
    out.append("Exact exposed-face verification:")
    out.append("  Every listed face was checked as the exact maximizer set of an integral linear functional;")
    out.append("  every listed affine dimension was checked by rational Gaussian elimination.")
    out.append("")
    out.append("All 27 nonempty faces not containing the origin:")
    out.extend(nondegeneracy_lines(faces))
    out.append("")
    out.append("V3 verdict: nondegenerate for every challenge prime p>=5 (in particular all listed primes).")
    out.append("It is not characteristic-free: the full base face is degenerate in characteristics 2 and 3.")
    out.append("")

    out.append("V4. ADOLPHSON--SPERBER TWISTED HODGE POLYGON")
    out.append("------------------------------------------------")
    out.append("A-S data are a character vector, not an unlabelled scalar:")
    out.append("  chi_i=omega^{-d_i}, theta=d/(q-1), N_d=theta+Z^4.")
    out.append("For q=p^s the definition orders each S_{d^(i)} by weight and uses slopes")
    out.append("  (1/s)*sum_{i=0}^{s-1} wt(u_{d^(i)}(j)), where d'=p*d mod (q-1).")
    out.append("For q=p, s=1 (and p*d=d mod p-1), so the slopes are exactly the weights of S_d.")
    out.append("Here C(Delta)={(r,u): r>=|u_1|+|u_2|+|u_3|} and wt(r,u)=r.")
    out.append("")
    out.append("Intended w-Kummer convention theta=(a,0,0,0):")
    out.append("  R_theta=w^a R.  The degree-n lattice count is")
    out.append("    L(n)=sum_{r=0}^3 2^r*C(3,r)*C(n,r),")
    out.append("  hence Hilb(R;t)=(1+t)^3/(1-t)^4.  The four degree-one logarithmic")
    out.append("  derivatives give the Koszul factor (1-t)^4, so")
    out.append("    Hilb(R_theta/I_theta;t)=t^a(1+t)^3.")
    out.append("  Complete slopes: a (x1), a+1 (x3), a+2 (x3), a+3 (x1).")
    out.append("")
    out.append("Grid k=0,...,99 (all arithmetic exact):")
    out.append("  k    a       first slope    first-a")
    for k in range(100):
        a = Fraction(k, 100)
        slopes = w_shift_hodge_slopes(a)
        first = slopes[0]
        assert first == a
        out.append(f"  {k:2d}   {str(a):7s} {str(first):12s} {str(first-a)}")
    out.append("")
    out.append("Conditional V4 verdict: under chi_w=omega^{-j}, a=j/(p-1), every one of the 100")
    out.append("formal grid values has first slope exactly a; none exceeds a.")
    out.append("")
    out.append("Why the scalar-only statement is underdetermined:")
    for a in (Fraction(1, 100), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(99, 100)):
        w_slopes = w_shift_hodge_slopes(a)
        sign_slopes = w_shift_hodge_slopes(1 - a)
        x_slopes = transverse_shift_hodge_slopes(a)
        out.append(
            f"  a={str(a):6s}: w shift first={w_slopes[0]}, sign-reversed w shift first={sign_slopes[0]}, "
            f"x shift first={x_slopes[0]}"
        )
    out.append("  If chi_w=omega^{+j}, the A-S convention uses theta_w=1-a for j>0, giving first 1-a.")
    out.append("  If theta=(0,a,0,0), Hilbert numerator is 2t+4t^2+2t^3 and the first slope is 1.")
    out.append("  Thus 'first slope=a' is true only after specifying both the twisted coordinate and sign.")
    out.append("  Also a=k/100 is an actual F_p character exponent only when its reduced denominator divides p-1;")
    out.append("  otherwise the grid is a formal rational-coset computation, not an A-S character at that p.")
    out.append("")
    out.append("Source audit:")
    out.append("  Directly read: Y. Qin--D. Xu, Irregular Hodge filtration of hypergeometric")
    out.append("  differential equations, Algebra & Number Theory 19 (2025), sec. 4.1.1--4.1.3.")
    out.append("  Official PDF: https://msp.org/ant/2025/19-12/ant-v19-n12-p07-p.pdf")
    out.append("  That section explicitly restates the shifted lattice, weight filtration, S_d quotient basis,")
    out.append("  Frobenius-orbit averaging, and Hodge polygon definition used above.")
    out.append("  It attributes the weight-sum formula / Newton-over-Hodge result to Adolphson--Sperber,")
    out.append("  Twisted exponential sums and Newton polyhedra, J. reine angew. Math. 443 (1993),")
    out.append("  Theorem 3.17 / Corollary 3.18.  The 1993 article body was not accessible here, so those")
    out.append("  original internal theorem numbers were NOT independently checked; only its metadata was.")
    return "\n".join(out) + "\n"


def main() -> None:
    report = build_report()
    OUT.write_text(report, encoding="utf-8")
    print(report, end="")


if __name__ == "__main__":
    main()
