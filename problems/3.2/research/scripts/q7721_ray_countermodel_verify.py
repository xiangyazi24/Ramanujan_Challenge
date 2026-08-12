#!/usr/bin/env python3
"""Q7721: finite structural countermodels for the moving Apéry ray problem.

This is a finite obstruction audit, not an asymptotic theorem.

It checks two complementary models at a fixed anchor N.

(A) CRT-Lucas model.
For every target prime p in (3N/4,N], set r=N-p and choose a digit map
f_p on [0,p) with zeros exactly at {r,p-1-r} and value 1 elsewhere.
Extend multiplicatively over base-p digits and CRT-lift the simultaneous
residues to integers c_m.  Then c_{ap+r} == c_a*c_r (mod p), reflection
holds, the zero set has size two, and every selected p divides c_N.  The
integer height is at most the primorial modulus, hence exp(O(N)).

(B) One common exact-Apery-recurrence model.
For each target p,r, choose a nonzero initial state modulo p for the *same*
Apéry order-two recurrence whose solution vanishes at r.  CRT-lift all these
local initial states to one characteristic-zero integer pair (A,B), run the
exact rational Apéry recurrence once, and clear denominators on the whole ray
prefix by a p-unit integer.  This yields one integral finite-prefix sequence
of exponential height satisfying the exact Apéry recurrence and realizing all
selected ray zeros.  Modulo each target p the full homogeneous solution is
reflection-fixed, hence p-1-r is also a zero.

The two models deliberately expose the remaining seam: model (A) has exact
Lucas+reflection+integrality but not the Apéry recurrence; model (B) has one
common initial state, exact Apéry recurrence+reflection+integrality on the ray
prefix, but not the Apéry Lucas law of the distinguished initial vector (1,5).
"""

from fractions import Fraction
from math import gcd, isqrt, lcm, log


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
    out = 1
    while n:
        d = n % p
        if d in zero_pair:
            return 0
        n //= p
    return out


def crt_pair(a, m, b, p):
    # x=a mod m, x=b mod p, gcd(m,p)=1
    t = (b-a) * pow(m, -1, p) % p
    return a + m*t, m*p


def symmetric_rep(x, m):
    x %= m
    return x-m if x > m//2 else x


def main():
    N = 500
    targets = []
    for p in primes_upto(N):
        if not (3*N < 4*p <= 4*N):
            continue
        r = N-p
        if 2 <= r < (p-1)//2:
            targets.append((p,r))
    assert targets

    # (A) Build CRT-lifted integers c_m for m<=N from multiplicative digit data.
    modulus = 1
    for p,_ in targets:
        modulus *= p
    c = []
    for n in range(N+1):
        x, mod = 0, 1
        for p,r in targets:
            z = {r, p-1-r}
            residue = digit_value(n,p,z)
            x, mod = crt_pair(x,mod,residue,p)
        assert mod == modulus
        c.append(symmetric_rep(x,modulus))

    for p,r in targets:
        z = {r,p-1-r}
        assert c[r] % p == 0
        assert c[p-1-r] % p == 0
        assert all((c[j] % p == 0) == (j in z) for j in range(p))
        assert c[N] % p == 0
        # For all m<=N, verify the Lucas recursion against the actual CRT lift.
        for m in range(N+1):
            q, rr = divmod(m,p)
            assert c[m] % p == (c[q] * c[rr]) % p
        assert all(c[j] % p == c[p-1-j] % p for j in range(p))

    # (B1) Determine one desired initial state modulo each target p.
    local_initial = []
    for p,r in targets:
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
    assert modA == modulus and modB == modulus
    A, B = symmetric_rep(A,modulus), symmetric_rep(B,modulus)

    # The relevant reflected ray indices r are all < N/4.  Run one exact Q-solution
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
        # The common initial state reduces to the chosen local state.
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
    print(f"N={N}")
    print(f"target_count={len(targets)}")
    print("targets=" + repr(targets))
    print(f"log_target_product={log_mass:.12f}")
    print(f"log_target_product_over_N={log_mass/N:.12f}")
    print(f"log_CRT_modulus={log(modulus):.12f}")
    print(f"common_initial_A_digits={len(str(abs(A)))}")
    print(f"common_initial_B_digits={len(str(abs(B)))}")
    print(f"ray_prefix_R={R}")
    print(f"cleared_recurrence_log_height={finite_height:.12f}")
    print(f"cleared_recurrence_log_height_over_N={finite_height/N:.12f}")
    print("CRT_LUCAS_REFLECTION_INTEGRAL_HEIGHT_MODEL=PASS")
    print("COMMON_INITIAL_EXACT_APERY_RECURRENCE_REFLECTION_MODEL=PASS")
    print("finite_countermodel_only=True")


if __name__ == "__main__":
    main()
