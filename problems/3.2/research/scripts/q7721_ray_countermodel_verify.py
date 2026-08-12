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

(B) Local exact-Apery-recurrence model.
For the same target p,r, choose a nonzero solution of the *same* Apéry
second-order recurrence modulo p that vanishes at r.  The full homogeneous
solution space is checked to be reflection-fixed, so p-1-r vanishes too.
Thus recurrence order/coefficient data plus reflection do not constrain the
cross-prime ray unless the actual Apéry initial vector is used.
"""

from math import gcd, isqrt, log


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
        if x > modulus//2:
            x -= modulus
        c.append(x)

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

    # (B) Same Apéry recurrence, p-dependent local state, forced ray zero.
    for p,r in targets:
        b = homogeneous_mod(p,1,5)
        v = homogeneous_mod(p,0,1)
        assert all(b[p-1-j] == b[j] for j in range(p))
        assert all(v[p-1-j] == v[j] for j in range(p))
        # Evaluation at r cannot annihilate both basis vectors: transfer is invertible.
        assert (b[r],v[r]) != (0,0)
        if v[r] != 0:
            lam = (-b[r] * pow(v[r],-1,p)) % p
            u = [(b[j] + lam*v[j]) % p for j in range(p)]
        else:
            u = v[:]  # then v[r]=0 and v is nonzero
        assert u[r] == 0
        assert u[p-1-r] == 0
        assert any(u)
        for n in range(1,p-1):
            assert ((n+1)**3*u[n+1] - P(n)*u[n] + n**3*u[n-1]) % p == 0

    log_mass = sum(log(p) for p,_ in targets)
    print(f"N={N}")
    print(f"target_count={len(targets)}")
    print("targets=" + repr(targets))
    print(f"log_target_product={log_mass:.12f}")
    print(f"log_target_product_over_N={log_mass/N:.12f}")
    print(f"log_CRT_modulus={log(modulus):.12f}")
    print("CRT_LUCAS_REFLECTION_INTEGRAL_HEIGHT_MODEL=PASS")
    print("LOCAL_EXACT_APERY_RECURRENCE_REFLECTION_MODEL=PASS")
    print("finite_countermodel_only=True")


if __name__ == "__main__":
    main()
