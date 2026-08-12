#!/usr/bin/env sage -python
"""Q7734: exact horizontal gap-carrier audit for Apéry coefficient zeros.

Run locally with Sage, e.g.

  sage -python problems/3.2/research/scripts/q7734_horizontal_gap_carrier.sage \
      --X 64 --mmax 4096

The script is deliberately finite and exact.  It checks:

1. the Apéry three-term recurrence modulo every prime p in (X,2X];
2. coefficient reflection b_{p-1-r}=b_r;
3. the canonical gap polynomials

       N_{h+1}(x)=P(x+h)N_h(x)-(x+h)^6 N_{h-1}(x),
       N_1=1, N_2=P(x+1),

   including degree 3(h-1), reflection
       N_h(-x-h-1)=(-1)^(h-1)N_h(x),
   and the even-h center factor;
4. the Q7734 short-nonreflection-gap carrier theorem on all scanned hits:
   if s=min(m mod p,p-1-(m mod p)) is an Apéry zero and t is another
   first-half zero at distance h=|t-s|, then p divides

       C_{m,h}=N_h(m) * (N_h(m-h) if m>=h else 1).

   The script additionally checks the orientation-specific stronger statement:
   the prime divides exactly the appropriate one of N_h(m), N_h(m-h).
5. a scoped recurrence-only countermodel.  For each p in (X,2X] it chooses
   the genuine Apéry recurrence with p-dependent initial data so that the same
   prescribed interior residue r=m0 is zero; direct recurrence propagation
   verifies that the reflected residue p-1-r is then also zero.  This preserves
   the actual recurrence coefficients but intentionally does NOT preserve the
   Apéry initial condition (1,5), Hasse normalization, or Cartier/Lucas data.

No probabilistic model is used.
"""

import argparse
from sage.all import GF, PolynomialRing, ZZ, binomial, prime_range


def P_int(n):
    return 34*n**3 + 51*n**2 + 27*n + 5


def apery_int(n):
    return ZZ(sum(binomial(n, k)**2 * binomial(n + k, k)**2
                  for k in range(n + 1)))


def apery_row_mod(p):
    """b_0,...,b_{p-1} from the actual recurrence over F_p."""
    K = GF(p)
    b = [K(0)] * p
    b[0] = K(1)
    if p > 1:
        b[1] = K(5)
    for n in range(1, p - 1):
        b[n + 1] = (K(P_int(n))*b[n] - K(n**3)*b[n - 1]) / K((n + 1)**3)
    return b


def make_gap_polynomials(hmax):
    R = PolynomialRing(ZZ, 'x')
    x = R.gen()

    def P(y):
        return 34*y**3 + 51*y**2 + 27*y + 5

    N = {1: R(1)}
    if hmax >= 2:
        N[2] = P(x + 1)
    for h in range(2, hmax):
        N[h + 1] = (P(x + h)*N[h] - (x + h)**6*N[h - 1]).expand()
    return R, x, N


def audit_gap_polynomials(x, N, hmax):
    for h in range(1, hmax + 1):
        Nh = N[h]
        assert Nh.degree() == 3*(h - 1)
        reflected = Nh(-x - h - 1)
        assert reflected == ((-1)**(h - 1))*Nh
        if h % 2 == 0:
            q, rem = Nh.quo_rem(2*x + h + 1)
            assert rem == 0
            assert q != 0


def first_half_zeros(b, p):
    half = (p - 1)//2
    return [r for r in range(half + 1) if b[r] == 0]


def fold_residue(r, p):
    s = min(r, p - 1 - r)
    branch = 'left' if r == s else 'right'
    return s, branch


def carrier_tag(branch, s, t):
    """0 means N_h(m); 1 means N_h(m-h)."""
    if branch == 'left':
        return 0 if t > s else 1
    return 1 if t > s else 0


def audit_actual_hits(X, mmax, N):
    primes = list(prime_range(X + 1, 2*X + 1))
    rows = {p: apery_row_mod(p) for p in primes}

    for p, b in rows.items():
        assert all(b[p - 1 - r] == b[r] for r in range(p))
        # Spot-check recurrence construction against the closed Apéry sum.
        for r in range(min(p, 12)):
            assert b[r] == GF(p)(apery_int(r))

    hits = 0
    nonreflection_gap_checks = 0
    carrier0 = 0
    carrier1 = 0

    for m in range(mmax):
        for p in primes:
            b = rows[p]
            r = m % p
            if b[r] != 0:
                continue
            hits += 1
            s, branch = fold_residue(r, p)
            assert b[s] == 0
            fh = first_half_zeros(b, p)
            assert s in fh

            for t in fh:
                if t == s:
                    continue
                h = abs(t - s)
                x0 = min(s, t)
                assert h >= 2  # no consecutive Apéry zeros
                assert h in N

                # Same-row recurrence constraint.
                assert N[h](x0) % p == 0

                tag = carrier_tag(branch, s, t)
                if tag == 0:
                    # s or its reflected congruence transports x0 to m.
                    assert N[h](m) % p == 0
                    carrier0 += 1
                else:
                    # This orientation can occur only when m>=h.
                    assert m >= h
                    assert N[h](m - h) % p == 0
                    carrier1 += 1

                C = ZZ(N[h](m))
                if m >= h:
                    C *= ZZ(N[h](m - h))
                assert C != 0
                assert C % p == 0
                nonreflection_gap_checks += 1

    return {
        'X': X,
        'primes': len(primes),
        'mmax': mmax,
        'hits': hits,
        'nonreflection_gap_checks': nonreflection_gap_checks,
        'carrier_Nh_m': carrier0,
        'carrier_Nh_m_minus_h': carrier1,
    }


def recurrence_solution_with_zero(p, r):
    """Full nonzero F_p solution with u_r=0,u_{r+1}=1.

    This uses the exact Apéry recurrence coefficients, but the initial pair is
    allowed to depend on p.  All divisions are by nonzero residues because the
    recurrence is propagated only through indices 1,...,p-2.
    """
    assert 1 <= r <= p - 2
    K = GF(p)
    u = [None] * p
    u[r] = K(0)
    u[r + 1] = K(1)

    # Forward from r+1 through p-1.
    for n in range(r + 1, p - 1):
        u[n + 1] = (K(P_int(n))*u[n] - K(n**3)*u[n - 1]) / K((n + 1)**3)

    # Backward from r through 0.
    for n in range(r, 0, -1):
        u[n - 1] = (K(P_int(n))*u[n] - K((n + 1)**3)*u[n + 1]) / K(n**3)

    assert all(v is not None for v in u)
    for n in range(1, p - 1):
        assert K((n + 1)**3)*u[n + 1] == K(P_int(n))*u[n] - K(n**3)*u[n - 1]
    return u


def audit_recurrence_countermodel(X):
    primes = list(prime_range(X + 1, 2*X + 1))
    if not primes:
        return {'countermodel_primes': 0}

    # Same fixed m0 for every prime; p>X>4*m0-4 guarantees m0<p/2.
    m0 = max(1, X // 4)
    rows = []
    for p in primes:
        assert 2*m0 < p
        u = recurrence_solution_with_zero(p, m0)
        reflected = p - 1 - m0
        assert u[m0] == 0
        assert u[reflected] == 0
        rows.append((p, int(u[0]), int(u[1]), reflected))

    return {
        'countermodel_primes': len(rows),
        'aligned_m': m0,
        'sample_initial_pairs': rows[:8],
    }


def audit_height_bound(X, N, hmax):
    """Check the explicit evaluation majorant used in the proof note.

    For n>=0 and R=n+h+1, the claimed bound is
      0 < N_h(n) <= 118^(h-1) R^(3(h-1)).
    We test it on the finite rectangle n<=min(X^2,256), h<=hmax.
    The mathematical proof is by induction from the N_h recurrence.
    """
    nmax = min(X*X, 256)
    checks = 0
    for h in range(1, hmax + 1):
        for n in range(nmax + 1):
            val = ZZ(N[h](n))
            bound = ZZ(118)**(h - 1) * ZZ(n + h + 1)**(3*(h - 1))
            assert 0 < val <= bound
            checks += 1
    return checks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--X', type=int, default=64)
    ap.add_argument('--mmax', type=int, default=None,
                    help='scan 0<=m<mmax; default X^2')
    args = ap.parse_args()

    X = args.X
    assert X >= 8
    mmax = X*X if args.mmax is None else min(args.mmax, X*X)

    # Any two first-half zeros have gap < p/2 < X for p<2X.
    hmax = X - 1
    R, x, N = make_gap_polynomials(hmax)
    audit_gap_polynomials(x, N, hmax)
    height_checks = audit_height_bound(X, N, min(hmax, 24))
    actual = audit_actual_hits(X, mmax, N)
    counter = audit_recurrence_countermodel(X)

    print('Q7734 GAP POLYNOMIALS', hmax)
    print('Q7734 HEIGHT CHECKS', height_checks)
    print('Q7734 ACTUAL', actual)
    print('Q7734 RECURRENCE COUNTERMODEL', counter)
    print('Q7734_HORIZONTAL_GAP_CARRIER PASS')


if __name__ == '__main__':
    main()
