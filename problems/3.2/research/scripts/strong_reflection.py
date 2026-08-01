#!/usr/bin/env python3
"""STRONG REFLECTION THEOREM test: EVERY solution germ of the Apery recurrence mod p
satisfies y_{p-1-n} = y_n on [0, p-1]. Test: random seeds (y_0, y_1), propagate, check."""
import random
rng = random.Random(3)
for p in (101, 199, 1009, 4001):
    okall = True
    for trial in range(20):
        y = [rng.randrange(p), rng.randrange(p)]
        ok = True
        for n in range(1, p-1):
            den = pow((n+1) % p, 3, p)
            yn1 = ((2*n+1)*(17*n*n+17*n+5) % p * y[n] - pow(n,3,p)*y[n-1]) * pow(den, p-2, p) % p
            y.append(yn1)
        # y defined on [0, p-1]
        sym = all(y[p-1-n] == y[n] for n in range(p))
        if not sym:
            okall = False
            bad = [n for n in range(p) if y[p-1-n] != y[n]][:5]
            print(f"p={p} trial{trial}: NOT symmetric, first bad n: {bad}")
            break
    if okall: print(f"p={p}: ALL 20 random germs reflection-symmetric  ✓")
