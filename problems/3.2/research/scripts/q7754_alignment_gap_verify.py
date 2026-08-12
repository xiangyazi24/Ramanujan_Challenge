#!/usr/bin/env python3
"""Q7754 finite audit: X=128 aligned triples and reflected gap certificates.

This is a finite verification only.  It imports the repository's exact HM3
zero-set enumerator, prints the active Z_p, enumerates all aligned triples,
and checks that the gap-polynomial certificate supplied by the mandatory
reflection pair is carried by the forced central factor 2*x+h+1=p.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import hm3_explore as hm3


def P(n: int) -> int:
    return 34*n**3 + 51*n**2 + 27*n + 5


def gap_value(h: int, x: int) -> int:
    """Evaluate N_h(x) by N_0=0,N_1=1 and the exact h-recurrence."""
    if h == 0:
        return 0
    if h == 1:
        return 1
    n0, n1 = 0, 1
    for j in range(1, h):
        n2 = P(x+j)*n1 - (x+j)**6*n0
        n0, n1 = n1, n2
    return n1


def main() -> None:
    X = 128
    all_data, checked = hm3.build_zero_sets(X, True)
    active = [z for z in all_data if z.zeros]
    pairs, triples = hm3.enumerate_pairs_and_triples(active, X*X)

    print(f"X={X} primes={len(all_data)} recurrence_checks={checked}")
    print(f"active={len(active)} sumZ={sum(len(z.zeros) for z in active)}")
    for z in active:
        print(f"Z_{z.p}={list(z.zeros)}")

    # One representative x per reflection orbit of each active zero set.
    orbit_rows = []
    accidental_extra = []
    for z in active:
        p = z.p
        seen = set()
        for r in z.zeros:
            x = min(r, p-1-r)
            if x in seen:
                continue
            seen.add(x)
            h = p-1-2*x
            assert h >= 0 and h % 2 == 0
            if h == 0:
                # central zero: reflection gives no distinct gap pair.
                orbit_rows.append((p,x,h,"central",None))
                continue
            N = gap_value(h,x)
            central_factor = 2*x+h+1
            assert central_factor == p
            assert N % p == 0
            cofactor_mod_p = (N//p) % p
            orbit_rows.append((p,x,h,central_factor,cofactor_mod_p))
            if cofactor_mod_p == 0:
                accidental_extra.append((p,x,h))

    print(f"reflection_orbits={len(orbit_rows)}")
    print("REFLECTION_GAP_ROWS p x h central_factor cofactor_mod_p")
    for row in orbit_rows:
        print(*row)
    print(f"accidental_p2_central_gap_count={len(accidental_extra)}")
    print(f"accidental_p2_central_gap_rows={accidental_extra}")

    print(f"canonical_aligned_triples={len(triples)}")
    for t in triples:
        hits = [(t.p,t.rp,t.m//t.p),(t.q,t.rq,t.m//t.q),(t.s,t.rs,t.m//t.s)]
        certs = []
        for p,r,q in hits:
            x = min(r,p-1-r)
            h = p-1-2*x
            if h == 0:
                certs.append((p,r,q,x,h,"central"))
            else:
                N = gap_value(h,x)
                assert 2*x+h+1 == p and N % p == 0
                certs.append((p,r,q,x,h,(N//p)%p))
        print(f"TRIPLE m={t.m} hits={hits} reflected_gap_certs={certs}")

    # HM3 cross-check from the exact scatter.
    k = hm3.scatter_k(active, X*X)
    S3 = sum(v*(v-1)*(v-2) for v in k)
    assert S3 == 6*len(triples)
    lam = sum(len(z.zeros)/z.p for z in all_data)
    expected = X*X*lam**3
    print(f"S3={S3}")
    print(f"lambda={lam:.12f}")
    print(f"X2_lambda3={expected:.12f}")
    print(f"R3={S3/expected:.12f}")
    print("Q7754_ALIGNMENT_GAP_VERIFY PASS")


if __name__ == '__main__':
    main()
