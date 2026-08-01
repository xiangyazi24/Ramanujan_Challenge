# Big-p family root statistics via projective coordinates + numpy vectorization.
# State (X,Y) mod p; map: (X,Y) -> ((u+1)^6 * Y, P(u)*Y - X); hit iff Y==0 (then state=(1,0)->next gives (0? ...)):
# careful: after hit (Y=0, X!=0), next step: X' = (u+1)^6*0 = 0? NO: map uses Y: X'=(u+1)^6*Y=0, Y'=P*Y-X=-X -> state (0,-X) ~ x=0 ✓ restart automatic.
import numpy as np, math, sys
def run(p, H):
    s = np.arange(p - H, dtype=np.int64)  # wrap-free windows
    X = np.zeros(p - H, dtype=np.int64)
    Y = np.ones(p - H, dtype=np.int64)
    total = 0; per_block = {}
    mx = 0
    for j in range(H):
        u = (s + j) % p
        Pu = (((34*u + 51) % p * u + 27) % p * u + 5) % p
        c6 = pow6 = np.mod((u + 1)**2 % p * (u + 1) % p, p); c6 = (c6 * c6) % p  # (u+1)^6 = ((u+1)^3)^2
        newX = (c6 * Y) % p
        newY = (Pu * Y - X) % p
        X, Y = newX, newY
        hits = int(np.count_nonzero(Y == 0))
        total += hits
        mx = max(mx, hits)
        blk = 1 << (j.bit_length())
        per_block[blk] = per_block.get(blk, [0, 0]); per_block[blk][0] += hits; per_block[blk][1] += 1
    R = total
    zb = 1 + R + p // H
    print(f"p={p} H={H} R={R} R/H={R/H:.3f} max_m={mx} implied|Z_p|<={zb} (=p^{math.log(zb)/math.log(p):.4f})", flush=True)
    print("  dyadic means:", {k: round(v[0]/v[1], 2) for k, v in sorted(per_block.items())}, flush=True)
for p, H in [(200003, 447), (500009, 707), (1000003, 1000)]:
    run(p, H)
