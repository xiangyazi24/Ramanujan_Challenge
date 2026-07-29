#!/usr/bin/env python3
"""Numerical pre-screen for Wirtinger-to-Seifert word map (P3.1).

Searches for words w_x, w_y in {s, s^-1, t, t^-1} such that:
- X = ρ(w_x) has tr(X) = 0  (i.e. X² = -I)
- Y = ρ(w_y) has tr(Y) = 1  (i.e. Y³ = -I)
- (XY)^17 = I
- tr(XY) = -2cos(π/17)

Uses floating-point arithmetic for speed; candidates verified in Sage.
"""
import mpmath as mp
import itertools

mp.mp.dps = 30

# Beta endpoint: s ≈ 0.4068
# Polynomial: s^16 - 7s^15 + ... + 1 = 0
coeffs = [1, -7, 22, -48, 87, -133, 178, -211, 223,
          -211, 178, -133, 87, -48, 22, -7, 1]
roots = mp.polyroots(coeffs, maxsteps=200)
real_roots = [r for r in roots if abs(r.imag) < 1e-20 and 0 < r.real < 1]
real_roots.sort(key=lambda r: abs(r.real - 0.4068))
s_beta = real_roots[0].real
print(f"s_beta = {mp.nstr(s_beta, 20)}")

xi = (1 - s_beta)**2 / (s_beta**2 - s_beta + 1)
M = s_beta
Mi = 1/s_beta

S = mp.matrix([[M, 1], [0, Mi]])
Si = mp.matrix([[Mi, -1], [0, M]])
T = mp.matrix([[M, 0], [xi, Mi]])
Ti = mp.matrix([[Mi, 0], [-xi, M]])

# Verify knot relation: S W T^-1 W^-1 = I, W = C^{-3}, C = T S^-1 T^-1 S
C = T * Si * Ti * S
W = mp.inverse(C * C * C)
rel = S * W * Ti * mp.inverse(W)
print(f"Knot relation check: |SWT⁻¹W⁻¹ - I| = {mp.nstr(mp.norm(rel - mp.eye(2)), 5)}")

gen = {'s': S, 'S': Si, 't': T, 'T': Ti}
inv_letter = {'s': 'S', 'S': 's', 't': 'T', 'T': 't'}

def word_matrix(word):
    A = mp.eye(2)
    for c in word:
        A = A * gen[c]
    return A

def tr(M):
    return M[0,0] + M[1,1]

target_trXY = -2 * mp.cos(mp.pi / 17)
print(f"Target tr(XY) = {mp.nstr(target_trXY, 15)}")

# Generate reduced words
def gen_words(maxlen):
    alphabet = 'sStT'
    yield ''
    frontier = list(alphabet)
    for w in frontier:
        yield w
    for length in range(2, maxlen + 1):
        new_frontier = []
        for w in frontier:
            for a in alphabet:
                if inv_letter[w[-1]] == a:
                    continue
                v = w + a
                yield v
                new_frontier.append(v)
        frontier = new_frontier

# Phase 1: find all words with tr = 0 or tr = 1 up to given length
MAXLEN = 12
print(f"\nSearching words up to length {MAXLEN}...")
print(f"(4 × 3^(L-1) words per length L)")

x_candidates = []
y_candidates = []
count = 0
for word in gen_words(MAXLEN):
    if not word:
        continue
    count += 1
    A = word_matrix(word)
    t = float(tr(A).real)
    if abs(t) < 1e-10:
        x_candidates.append((word, A))
    if abs(t - 1) < 1e-10:
        y_candidates.append((word, A))

print(f"Total words checked: {count}")
print(f"X candidates (tr=0): {len(x_candidates)}")
print(f"Y candidates (tr=1): {len(y_candidates)}")

# Phase 2: check pairs
print(f"\nChecking {len(x_candidates)} × {len(y_candidates)} = {len(x_candidates)*len(y_candidates)} pairs...")
found = 0
for wx, X in x_candidates:
    for wy, Y in y_candidates:
        Q = X * Y
        trQ = float(tr(Q).real)
        if abs(trQ - float(target_trXY)) > 0.001:
            continue
        # Check (XY)^17 = I
        Qpow = mp.eye(2)
        for _ in range(17):
            Qpow = Qpow * Q
        if mp.norm(Qpow - mp.eye(2)) > 1e-5:
            continue
        found += 1
        print(f"\n*** CANDIDATE FOUND ***")
        print(f"  w_x = {wx} (length {len(wx)})")
        print(f"  w_y = {wy} (length {len(wy)})")
        print(f"  tr(X) = {mp.nstr(tr(X), 15)}")
        print(f"  tr(Y) = {mp.nstr(tr(Y), 15)}")
        print(f"  tr(XY) = {mp.nstr(tr(Q), 15)}")
        print(f"  |(XY)^17 - I| = {mp.nstr(mp.norm(Qpow - mp.eye(2)), 10)}")

if found == 0:
    print(f"\nNo certificate found through length {MAXLEN}.")
    print("Try MAXLEN = 14 or 16.")
else:
    print(f"\nFound {found} candidate pair(s). Verify in Sage for exact proof.")
