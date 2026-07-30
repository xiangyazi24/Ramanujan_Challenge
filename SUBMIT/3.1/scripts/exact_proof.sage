## Exact verification of the torsion mechanism for BOTH endpoints.
##
## Step 1: the palindromic minimal polynomial factors through w = a + 1/a;
##         if that degree-n/2 polynomial is TOTALLY REAL and we count how many
##         of its roots lie in [-2,2], we know exactly how many roots of the
##         original lie on the unit circle.
## Step 2: |a| = 1  ==>  u is real   (pure algebra, shown symbolically below)
## Step 3: u real   ==>  T real
## Step 4: |a| = 1 and u real  ==>  W = 1/(1 - conj V)
## Step 5: Bloch-Wigner: D(1/(1-z)) = D(z), D(conj z) = -D(z)
##  ==>  sum_j D(z_j) = 0 exactly at every complex embedding.

R.<z> = PolynomialRing(QQ)
falpha = (z^12 - 3*z^11 + 4*z^10 - 5*z^9 + 6*z^8 - 7*z^7
          + 7*z^6 - 7*z^5 + 6*z^4 - 5*z^3 + 4*z^2 - 3*z + 1)
fbeta  = (z^16 - 7*z^15 + 22*z^14 - 48*z^13 + 87*z^12 - 133*z^11 + 178*z^10
          - 211*z^9 + 223*z^8 - 211*z^7 + 178*z^6 - 133*z^5 + 87*z^4
          - 48*z^3 + 22*z^2 - 7*z + 1)

def analyse(f, name):
    print("="*66)
    print(name, " degree", f.degree())
    c = f.coefficients(sparse=False)
    print("  palindromic:", c == list(reversed(c)))
    print("  irreducible:", f.is_irreducible())
    # depress:  f(a) = a^{n/2} * g(a + 1/a)
    n = f.degree()
    S.<w> = PolynomialRing(QQ)
    # build g by Chebyshev-like reduction
    a = R.gen()
    # a^k + a^{-k} = p_k(w)
    p = [S(2), S(w)]
    for k in range(2, n//2 + 1):
        p.append(w*p[k-1] - p[k-2])
    g = S(0)
    for k in range(n//2 + 1):
        coef = c[n//2 + k] if k > 0 else c[n//2]
        g += coef * (p[k] if k > 0 else S(1))
    print("  g(w) =", g)
    print("  g degree:", g.degree())
    rr = g.roots(RR, multiplicities=False)
    print("  g real roots:", len(rr), "of", g.degree(), " => totally real:", len(rr) == g.degree())
    inside  = [r for r in rr if -2 <= r <= 2]
    outside = [r for r in rr if not (-2 <= r <= 2)]
    print("  roots in [-2,2] :", len(inside), " -> that many UNIT-CIRCLE pairs of f")
    print("  roots outside   :", len(outside), " -> that many REAL pairs of f")
    print("  so f has %d roots on |z|=1 and %d real roots"
          % (2*len(inside), 2*len(outside)))
    # cross-check numerically
    allr = f.roots(CC, multiplicities=False)
    onunit = sum(1 for t in allr if abs(abs(t) - 1) < 1e-25)
    realr  = sum(1 for t in allr if abs(t.imag()) < 1e-25)
    print("  direct count: %d on the unit circle, %d real" % (onunit, realr))
    print("  CONSISTENT:", onunit == 2*len(inside) and realr == 2*len(outside))
    return g

ga = analyse(falpha, "f_alpha")
gb = analyse(fbeta,  "f_beta")

print()
print("="*66)
print("Step 2 symbolically:  u is real whenever |a| = 1")
print("="*66)
F.<A> = FractionField(PolynomialRing(QQ, 'A'))
# u(a) = (L + X^3)/(X(L+X)) with L = a, X = a^4   (alpha chart)
u_alpha = (A + A^12) / (A^4 * (A + A^4))
print("alpha chart:  u(a) =", u_alpha.factor())
# conj on |a|=1 is a -> 1/a ; check u(1/a) == u(a)
u_inv = u_alpha.subs(A=1/A)
print("u(1/a) - u(a) =", (u_inv - u_alpha).simplify_full() if hasattr(u_inv-u_alpha,'simplify_full') else (u_inv - u_alpha))
print("u(1/a) == u(a) :", bool(u_inv == u_alpha))

# beta chart: L = a, X = a^2
u_beta = (A + A^6) / (A^2 * (A + A^2))
print()
print("beta chart:   u(a) =", u_beta.factor())
print("u(1/a) == u(a) :", bool(u_beta.subs(A=1/A) == u_beta))
