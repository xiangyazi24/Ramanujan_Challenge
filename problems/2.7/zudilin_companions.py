"""
Compute Zudilin's exact companion sequences p^{(2)}_n and p^{(3)}_n
from the rational function R_n(t) = ((t-1)...(t-n))^3 / ((n!)^2 * t(t+1)...(t+n)).

r_n^{(2)} = -sum_{nu>=1} R_n'(nu) = b_n * zeta(2) - p_n^{(2)}
r_n^{(3)} = (1/2) sum_{nu>=1} R_n''(nu) = b_n * zeta(3) - p_n^{(3)}

p^{(2)}_n and p^{(3)}_n are rational, and b_n * L - (p^{(2)}_n + p^{(3)}_n) -> 0
where L = zeta(2) + zeta(3).
"""
from fractions import Fraction as Q
from math import comb, factorial

def compute_b(n):
    return sum(comb(n,k)**2 * comb(n+k,n) * comb(n+2*k,n) for k in range(n+1))

# R_n(t) = prod_{j=1}^n (t-j)^3 / (n!^2 * prod_{j=0}^n (t+j))
# At positive integer nu > n, this is:
# R_n(nu) = prod_{j=1}^n (nu-j)^3 / (n!^2 * prod_{j=0}^n (nu+j))
# For nu = 1,...,n: R_n(nu) = 0 (triple zero from numerator)

# For R_n'(nu), R_n''(nu) at nu=1,...,n: these are also 0 (order 3 zero)

# For nu > n:
# log R_n(t) = 3 sum_{j=1}^n log(t-j) - 2 log(n!) - sum_{j=0}^n log(t+j)
# R_n'(t)/R_n(t) = 3 sum_{j=1}^n 1/(t-j) - sum_{j=0}^n 1/(t+j)
# R_n'(t) = R_n(t) * [3 sum 1/(t-j) - sum 1/(t+j)]

# For the SECOND derivative:
# (R'/R)' = -3 sum 1/(t-j)^2 + sum 1/(t+j)^2
# R''/R = (R'/R)' + (R'/R)^2

# Since R_n(nu) = 0 for nu <= n, the sums are Σ_{nu=n+1}^∞

def R_n_val(n_val, t):
    """Compute R_n(t) as exact fraction."""
    numer = Q(1)
    for j in range(1, n_val+1):
        numer *= (t - Q(j))**3
    denom = Q(factorial(n_val))**2
    for j in range(0, n_val+1):
        denom *= (t + Q(j))
    return numer / denom

def R_n_logderiv_parts(n_val, nu):
    """Compute R_n'/R_n at t=nu, and R_n''/R_n - (R_n'/R_n)^2 at t=nu."""
    # d/dt log R = 3 * sum_{j=1}^n 1/(t-j) - sum_{j=0}^n 1/(t+j)
    S1 = Q(0)
    for j in range(1, n_val+1):
        S1 += Q(1) / Q(nu - j)
    S2 = Q(0)
    for j in range(0, n_val+1):
        S2 += Q(1) / Q(nu + j)
    psi = Q(3) * S1 - S2

    # d^2/dt^2 log R = -3 * sum 1/(t-j)^2 + sum 1/(t+j)^2
    S1sq = Q(0)
    for j in range(1, n_val+1):
        S1sq += Q(1) / Q(nu - j)**2
    S2sq = Q(0)
    for j in range(0, n_val+1):
        S2sq += Q(1) / Q(nu + j)**2
    phi = Q(-3) * S1sq + S2sq

    return psi, phi

def compute_zudilin_companions(n_val, max_nu=None):
    """
    Compute p^{(2)}_n and p^{(3)}_n.

    r_n^{(2)} = -sum_{nu>=1} R_n'(nu) = b_n*zeta(2) - p^{(2)}_n
    r_n^{(3)} = (1/2)*sum_{nu>=1} R_n''(nu) = b_n*zeta(3) - p^{(3)}_n

    R_n'(nu) = R_n(nu) * psi(nu)
    R_n''(nu) = R_n(nu) * (phi(nu) + psi(nu)^2)

    For nu=1,...,n: all zero.
    For nu>n: R_n(nu) is nonzero, sum converges since R_n(nu) ~ 1/nu^3.

    p^{(2)}_n is the rational part of b_n*zeta(2) + sum R_n'(nu).
    Since r_n^{(2)} = b_n*zeta(2) - p^{(2)}_n -> 0 exponentially,
    p^{(2)}_n = b_n*zeta(2) - r_n^{(2)} = b_n*zeta(2) + sum R_n'(nu).

    But p^{(2)}_n is supposed to be RATIONAL. The sum has no zeta values?
    Actually: R_n'(nu) at integer nu involves only rationals.
    R_n(nu) = product of rational numbers.
    psi(nu) = sum of 1/(integer) terms.
    So R_n'(nu) = R_n(nu)*psi(nu) is rational.
    And r_n^{(2)} = -sum R_n'(nu) is a sum of rationals = rational.
    Then p^{(2)}_n = b_n*zeta(2) - r_n^{(2)}... but this contains zeta(2)!

    Wait, this can't be right. p^{(2)}_n must be rational. So:
    b_n*zeta(2) - p^{(2)}_n = r_n^{(2)} = -sum R_n'(nu)
    => p^{(2)}_n = b_n*zeta(2) + sum R_n'(nu)

    But the left side is rational, so b_n*zeta(2) + sum R_n'(nu) must
    be rational. But b_n*zeta(2) is irrational...

    UNLESS the sum has a zeta(2) component that cancels.

    Let me reconsider: R_n'(nu) = R_n(nu) * psi(nu). The sum_{nu>=n+1} is finite
    in terms of zeta? No, it's an infinite sum of rationals.

    Actually, the resolution is that the partial fractions of R_n(t) at
    poles t = 0, -1, ..., -n contain zeta values. The sum over positive
    integers involves these. Let me think again...

    Actually, I think the correct formulation involves evaluating R_n at
    specific points that produce harmonic numbers, and the zeta values come
    from extending the sums. Let me just compute directly.

    For nu > n+1, R_n(nu) decreases as ~1/nu^3, so the tail is small.
    p^{(2)}_n can be extracted from:
      sum_{nu=1}^M R_n'(nu) = -b_n*zeta(2) + p^{(2)}_n + O(tail)
    But we need to handle the zeta part carefully.

    Actually, the standard approach: expand 1/(t+j) = sum_m 1/(nu+j) etc.
    Let me use the EXPLICIT formula via harmonic numbers.

    R_n(nu) for nu > n:
    prod_{j=1}^n (nu-j) = (nu-1)!/(nu-n-1)! = C(nu-1, n) * n!
    So (nu-j)^3 product: ((nu-1)!/(nu-n-1)!)^3 * ... no wait.

    prod_{j=1}^n (nu-j)^3 = [prod_{j=1}^n (nu-j)]^3 = [(nu-1)!/(nu-n-1)!]^3

    prod_{j=0}^n (nu+j) = (nu+n)!/((nu-1)!)

    R_n(nu) = [(nu-1)!/(nu-n-1)!]^3 / [(n!)^2 * (nu+n)!/(nu-1)!]
            = [(nu-1)!]^3 / [(nu-n-1)!]^3 * (nu-1)! / [(n!)^2 * (nu+n)!]
            = [(nu-1)!]^4 / [(nu-n-1)!]^3 * 1/[(n!)^2 * (nu+n)!]

    This can be written in terms of binomials:
    = C(nu-1,n)^3 * (n!)^3 / (n!)^2 * (nu-1)! / (nu+n)!
    = C(nu-1,n)^3 * n! * (nu-1)! / (nu+n)!
    = C(nu-1,n)^3 * n! / [(nu)(nu+1)...(nu+n)]
    = C(nu-1,n)^3 / C(nu+n, n+1) ... hmm getting complicated.

    Let me just compute numerically for small n.
    """
    if max_nu is None:
        max_nu = 3*n_val + 100  # Enough for convergence

    sum_Rprime = Q(0)
    sum_Rdblprime = Q(0)

    for nu in range(n_val + 1, max_nu + 1):
        Rval = R_n_val(n_val, Q(nu))
        psi, phi = R_n_logderiv_parts(n_val, nu)

        Rprime = Rval * psi
        Rdblprime = Rval * (phi + psi**2)

        sum_Rprime += Rprime
        sum_Rdblprime += Rdblprime

    # r_n^{(2)} = -sum R_n'(nu)
    # r_n^{(3)} = (1/2) sum R_n''(nu)
    r2 = -sum_Rprime
    r3 = sum_Rdblprime / Q(2)

    return r2, r3

N = 8
b_vals = [Q(compute_b(n)) for n in range(N)]

print("Computing Zudilin companion sequences...")
for n in range(N):
    r2, r3 = compute_zudilin_companions(n, max_nu=n + 200)
    print(f"\nn={n}:")
    print(f"  b_n = {b_vals[n]}")
    print(f"  r_n^(2) = {float(r2):.15e}")
    print(f"  r_n^(3) = {float(r3):.15e}")
    print(f"  r_n^(2) + r_n^(3) = {float(r2 + r3):.15e}")
    # p^{(2)}_n = b_n*zeta(2) - r_n^{(2)} -- contains zeta!
    # Instead check: b_n*(zeta(2)+zeta(3)) - (r_n^{(2)} + r_n^{(3)}) should give p^{(2)}_n + p^{(3)}_n
    # r2 + r3 is rational and small
    # Hmm, but r2 = b_n*zeta(2) - p^{(2)}_n, so r2 is NOT rational -- it's the error!
    # p^{(2)}_n IS rational. So r2 = b_n*zeta(2) - p^{(2)}_n ~ 0 but irrational.

    # Wait, the sum_{nu} R_n'(nu) IS rational (sum of rationals at integers).
    # But r_n^{(2)} = -sum R_n'(nu) is supposed to equal b_n*zeta(2) - p^{(2)}_n.
    # That means -sum R_n'(nu) = b_n*zeta(2) - p^{(2)}_n.
    # So p^{(2)}_n = b_n*zeta(2) + sum R_n'(nu).
    # But this contains zeta(2)!

    # Resolution: The sum is INFINITE (nu from 1 to infinity).
    # sum_{nu=1}^infty 1/nu^2 = zeta(2), so the infinite sum of R_n'(nu) DOES
    # contain a zeta(2) component. The decomposition is:
    # -sum R_n'(nu) = b_n*zeta(2) - p^{(2)}_n
    # where the zeta(2) part of -sum R_n'(nu) is exactly b_n*zeta(2), and
    # the remaining rational part is -p^{(2)}_n.

    # So I need to SEPARATE the zeta(2) and rational parts of the sum.
    # This requires partial fraction decomposition + harmonic number identities.

    # For now, just check the COMBINED: r2 + r3 should be close to 0 for large n.
    # Nah, r2 and r3 are the ERRORS, they should go to 0.

    # Actually let me check: r2 approaches 0? It equals b_n*zeta(2) - p^{(2)}_n.
    # Since this is the error in the diophantine approximation, yes it -> 0.
    # But for n=0, r_0^{(2)} should be nonzero.

    # Hmm, but I computed r2 as -sum R_n'(nu) which is RATIONAL (a finite truncation).
    # The FULL sum -sum_{nu=1}^infty R_n'(nu) contains zeta values.
    # My max_nu truncation misses the tail which contains the zeta values.

    # I think the correct formula uses a different approach: partial fractions.
    # Let me try a different method.

print("\n\nNote: The above approach has convergence issues because the infinite")
print("sum of R_n'(nu) involves zeta values, not just rationals.")
print("Need partial fraction decomposition to extract p^{(2)}_n exactly.")
