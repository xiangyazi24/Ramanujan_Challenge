from sage.all import ZZ, binomial, gcd, prime_range, prod


def T(n, k):
    if k < 0 or k > n:
        return ZZ(0)
    return ZZ(binomial(n, k)**2 * binomial(n + k, k)**2)


def apery(n):
    return sum(T(n, k) for k in range(n + 1))


def slices(n):
    J = (n - 1) // 3
    L = sum(T(n, k) for k in range(max(J + 1, 0)))
    H = sum(T(n, k) for k in range((n + 1)//2, n + 1))
    return ZZ(L), ZZ(H)


def top_primes(n):
    return [ZZ(p) for p in prime_range(n//2 + 1, n + 1)]


def folded_index(n, p):
    r = n - p
    return min(r, p - 1 - r)


def verify_congruences(N=500):
    for n in range(1, N + 1):
        A = apery(n)
        L, H = slices(n)
        for p in top_primes(n):
            j = folded_index(n, p)
            Aj = apery(j)
            assert (L - Aj) % p == 0
            assert (H - 4*Aj) % p == 0
            if p >= 7:
                assert ((A % p) == 0) == ((Aj % p) == 0)
                assert ((A % p) == 0) == ((L % p) == 0)
                assert ((A % p) == 0) == ((H % p) == 0)
    print("Lucas truncation congruences: PASS")


def verify_square_primorial(N=500):
    for n in range(1, N + 1):
        L, H = slices(n)
        P = prod(top_primes(n), ZZ(1))
        assert (H - 4*L) % (P*P) == 0
    print("(top-half primorial)^2 divisibility: PASS")


def audit_one(n):
    A = apery(n)
    L, H = slices(n)
    ps = top_primes(n)
    P = prod(ps, ZZ(1))
    R = prod([p for p in ps if A % p == 0], ZZ(1))
    g = gcd(L, H)
    print("n =", n)
    print("top-half hits =", [p for p in ps if A % p == 0])
    print("R_n =", R)
    print("gcd(L,H) =", g)
    print("top support of gcd =", gcd(g, P))
    print("v_p(H-4L) at hits =",
          [(p, (H-4*L).valuation(p)) for p in ps if A % p == 0])
    assert gcd(g, P) == R
    assert (H - 4*L) % (P*P) == 0


verify_congruences(300)
verify_square_primorial(300)
audit_one(321)
