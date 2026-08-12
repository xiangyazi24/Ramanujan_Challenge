#!/usr/bin/env python3
"""Q7721: finite structural countermodels for the moving Apéry ray problem.

This is a finite obstruction audit, not an asymptotic theorem.

It checks two complementary models at a fixed anchor N.

(A) CRT-Lucas model, uniform over every prime p<=N at this scale.
For target primes p in (3N/4,N], set r=N-p and choose a digit map f_p
with zeros exactly at {r,p-1-r}; for every other prime p<=N use the
constant digit map 1.  Extend every f_p multiplicatively over base-p digits
and CRT-lift all primewise residues to one integer c_m for each m<=N.
Then the resulting finite integer sequence satisfies the p-Lucas law and
reflection for every prime p<=N, has no zero digits at unselected primes,
and every selected p divides c_N.  Its height is at most the full primorial
modulus exp(theta(N))=exp(O(N)).

(B) One common exact-Apery-recurrence model.
For each target p,r, choose a nonzero initial state modulo p for the *same*
Apéry order-two recurrence whose solution vanishes at r.  CRT-lift all these
local initial states to one characteristic-zero integer pair (A,B), run the
exact rational Apéry recurrence once, and clear denominators on the whole ray
prefix by a p-unit integer.  This yields one integral finite-prefix sequence
of exponential height satisfying the exact Apéry recurrence and realizing all
selected ray zeros.  Modulo each target p the full homogeneous solution is
reflection-fixed, hence p-1-r is also a zero.

The two models expose the remaining seam: model (A) has exact Lucas +
reflection + integrality uniformly over all primes at the scale, but not the
Apéry recurrence; model (B) has one common initial state, exact Apéry recurrence
+ reflection + integrality on the ray prefix, but not the Apéry Lucas law of
the distinguished initial vector (1,5).
"""

from fractions import Fraction
from math import isqrt, lcm, log


def primes_upto(n):
    sieve = bytearray(b"\x01") * (n + 1)
    if n >= 0:
        sieve[0] = 0
    if n >= 1:
        sieve[1] = 0
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            sieve[p*p:n+1:p] = b"\x00" * (((n-p*p)//p)+1)
    return [p for p in range(2, n+1) if sieve[p]]


def P(n):
    return 34*n**3 + 51*n**2 + 27*n + 5


def homogeneous_mod(p, y0, y1):
    y = [0]*p
    y[0], y[1] = y0 % p, y1 % p
    for n in range(1, p-1):
        y[n+1] = (P(n)*y[n] - n**3*y[n-1]) * pow((n+1)**3, -1, p) % p
    return y


def homogeneous_Q(limit, y0, y1):
    y = [Fraction(y0), Fraction(y1)]
    for n in range(1, limit):
        y.append((P(n)*y[n] - n**3*y[n-1]) / (n+1)**3)
    return y


def digit_value(n, p, zero_pair):
    if n == 0:
        return 1
    while n:
        d = n % p
        if d in zero_pair:
            return 0
        n //= p
    return 1


def crt_pair(a, m, b, p):
    # x=a mod m, x=b mod p, gcd(m,p)=1
    t = (b-a) * pow(m, -1, p) % p
    return a + m*t, m*p


def symmetric_rep(x, m):
    x %= m
    return x-m if x > m//2 else x


def main():
    N = 500
    primes = primes_upto(N)
    targets = []
    zero_pairs = {p: set() for p in primes}
    for p in primes:
        if not (3*N < 4*p <= 4*N):
            continue
        r = N-p
        if 2 <= r < (p-1)//2:
            targets.append((p,r))
            zero_pairs[p] = {r,p-1-r}
    assert targets

    # (A) Build CRT-lifted integers c_m for m<=N using *all* primes p<=N.
    all_modulus = 1
    for p in primes:
        all_modulus *= p
    c = []
    for n in range(N+1):
        x, mod = 0, 1
        for p in primes:
            residue = digit_value(n,p,zero_pairs[p])
            x, mod = crt_pair(x,mod,residue,p)
        assert mod == all_modulus
        c.append(symmetric_rep(x,all_modulus))

    target_dict = dict(targets)
    for p in primes:
        z = zero_pairs[p]
        # Exact reflected first-block zero set: target pair or empty.
        assert all((c[j] % p == 0) == (j in z) for j in range(p))
        # For all m<=N, verify p-Lucas against the actual CRT lift.
        for m in range(N+1):
            q, rr = divmod(m,p)
            assert c[m] % p == (c[q] * c[rr]) % p
        assert all(c[j] % p == c[p-1-j] % p for j in range(p))
        if p in target_dict:
            r = target_dict[p]
            assert c[r] % p == 0 and c[p-1-r] % p == 0
            assert c[N] % p == 0
        else:
            assert c[N] % p != 0

    # (B1) Determine one desired initial state modulo each target p.
    local_initial = []
    target_modulus = 1
    for p,r in targets:
        target_modulus *= p
        b = homogeneous_mod(p,1,5)
        v = homogeneous_mod(p,0,1)
        assert all(b[p-1-j] == b[j] for j in range(p))
        assert all(v[p-1-j] == v[j] for j in range(p))
        # Evaluation at r cannot annihilate both basis vectors: transfer is invertible.
        assert (b[r],v[r]) != (0,0)
        if v[r] != 0:
            lam = (-b[r] * pow(v[r],-1,p)) % p
            a0, a1 = 1, (5+lam) % p
        else:
            a0, a1 = 0, 1
        u = homogeneous_mod(p,a0,a1)
        assert u[r] == 0 and u[p-1-r] == 0 and any(u)
        assert all(u[p-1-j] == u[j] for j in range(p))
        local_initial.append((p,r,a0,a1))

    # (B2) CRT-lift those local states to one characteristic-zero initial pair.
    A, modA = 0, 1
    B, modB = 0, 1
    for p,r,a0,a1 in local_initial:
        A, modA = crt_pair(A,modA,a0,p)
        B, modB = crt_pair(B,modB,a1,p)
    assert modA == target_modulus and modB == target_modulus
    A, B = symmetric_rep(A,target_modulus), symmetric_rep(B,target_modulus)

    # The relevant ray indices r are all < N/4. Run one exact Q-solution
    # through the largest ray index and clear denominators by 6*lcm(1..R)^3.
    R = max(r for _,r in targets)
    uQ = homogeneous_Q(R,A,B)
    L = 1
    for k in range(1,R+1):
        L = lcm(L,k)
    clear = 6*L**3
    assert all(clear % p for p,_ in targets)
    C = []
    for value in uQ:
        z = value*clear
        assert z.denominator == 1
        C.append(z.numerator)

    for p,r,a0,a1 in local_initial:
        assert A % p == a0 and B % p == a1
        assert C[r] % p == 0
        full = homogeneous_mod(p,A,B)
        assert full[r] == 0 and full[p-1-r] == 0
        assert all(full[p-1-j] == full[j] for j in range(p))
        for n in range(1,p-1):
            assert ((n+1)**3*full[n+1] - P(n)*full[n] + n**3*full[n-1]) % p == 0

    # Exact recurrence over Z on the denominator-cleared ray prefix.
    for n in range(1,R):
        assert (n+1)**3*C[n+1] == P(n)*C[n] - n**3*C[n-1]

    log_mass = sum(log(p) for p,_ in targets)
    finite_height = max(log(max(1,abs(x))) for x in C)
    lucas_height = max(log(max(1,abs(x))) for x in c)
    print(f"N={N}")
    print(f"all_prime_count={len(primes)}")
    print(f"target_count={len(targets)}")
    print("targets=" + repr(targets))
    print(f"log_target_product={log_mass:.12f}")
    print(f"log_target_product_over_N={log_mass/N:.12f}")
    print(f"log_all_prime_CRT_modulus={log(all_modulus):.12f}")
    print(f"lucas_model_log_height={lucas_height:.12f}")
    print(f"lucas_model_log_height_over_N={lucas_height/N:.12f}")
    print(f"common_initial_A_digits={len(str(abs(A)))}")
    print(f"common_initial_B_digits={len(str(abs(B)))}")
    print(f"ray_prefix_R={R}")
    print(f"cleared_recurrence_log_height={finite_height:.12f}")
    print(f"cleared_recurrence_log_height_over_N={finite_height/N:.12f}")
    print("ALL_PRIMES_CRT_LUCAS_REFLECTION_INTEGRAL_HEIGHT_MODEL=PASS")
    print("COMMON_INITIAL_EXACT_APERY_RECURRENCE_REFLECTION_MODEL=PASS")
    print("finite_countermodel_only=True")


if __name__ == "__main__":
    main()
