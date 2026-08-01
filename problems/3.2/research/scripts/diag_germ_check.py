#!/usr/bin/env python3
"""Verify: #{(r,h<=H): sigma_h(r)=sigma_{h+1}(r), type-I} equals
#{(x,r): 2 <= x-r <= H+1, Y(x,r)=1} where Y(x,r) = y_r for the germ with state (y_x, y_{x-1}) = (1,1).
p=101, H=10."""
p, H = 101, 12
def AB(d, rr):
    Ap_, Bp_ = 1, 0; Ac = Bc = None
    for j in range(d):
        x = (rr + j) % p
        den = pow((x+1) % p, 3, p)
        if den == 0: return None, None
        di = pow(den, p-2, p)
        a = (2*x+1)*(17*x*x+17*x+5) % p * di % p
        be = (-(x**3)) % p * di % p
        if j == 0: Ac, Bc = a, be
        else:
            Ac, Ap_ = (a*Ac + be*Ap_) % p, Ac
            Bc, Bp_ = (a*Bc + be*Bp_) % p, Bc
    return Ac, Bc

# (a) direct count
direct = set()
for r in range(1, p-2):
    for h in range(1, H+1):
        if r + h + 1 > p-2: break
        A1, B1 = AB(h, r); A2, B2 = AB(h+1, r)
        if A1 is None or A2 is None or not B1 or not B2: continue
        s1 = (1-A1)*pow(B1, p-2, p) % p
        s2 = (1-A2)*pow(B2, p-2, p) % p
        if s1 == s2: direct.add((r, h))
# (b) diagonal-germ count: germ with (y_x, y_{x-1}) = (1,1); propagate BACKWARD to r
diag = set()
for x in range(2, p-1):
    # backward: y_{n-1} = (P(n) y_n - (n+1)^3 y_{n+1}) / n^3
    y = {x: 1, x-1: 1}
    n = x-1
    while n >= 1:
        # recurrence at m=n: (n+1)^3 y_{n+1} = P(n) y_n - n^3 y_{n-1}
        yn1 = ( (2*n+1)*(17*n*n+17*n+5) % p * y[n] - pow(n+1,3,p)*y[n+1] ) * pow(pow(n,3,p), p-2, p) % p
        y[n-1] = yn1
        n -= 1
        if x - (n) > H + 3: break
    for r in range(max(1, x-H-1), x-1):
        if y.get(r) == 1:
            h = x - r - 1
            if 1 <= h <= H: diag.add((r, h))
print(f"direct sigma_h=sigma_(h+1) events: {len(direct)}")
print(f"diagonal-germ Y=1 events:          {len(diag)}")
print(f"match: {direct == diag}; sym-diff sample: {list(direct ^ diag)[:6]}")
