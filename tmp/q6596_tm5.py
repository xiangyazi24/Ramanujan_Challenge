"""TM5 stream module: the Q5221 D33 test (i = 3, j = 5, m = 32) as adjacent-block address comparison.
Substrate (ASSUMED static and flushed, never written here): ph_i = Lemma-P phase mod 8 of the label window starting at site i
(None if the mu-coordinate is not a Thue-Morse factor there), computed by the oracle `phases`.
Stream state per site (l3, l4, cy, v3, v4, r): 64 states. Synchronous, one-way, radius one.
 l3 = address bit 3 of the latest standard 3-block (launched at block offset b_3 = 3), l4 likewise (offset b_4 = 6),
 cy = the previous block's bit 3 (the carry into bit 4), v3 / v4 = lane validity, r = reject.
Mutants: M_NOPOS (no phase-continuity / validity check), M_NOLAUNCH (lane never rewritten from the label)."""
import sys, itertools
def tail_codes(count=64):
    codes, inv, value, index = {}, {}, 3, 3
    while index < count + 3:
        b = format(value, "b")
        if b.startswith("11") and "11" not in b[1:]: codes[index] = value; inv[value] = index; index += 1
        value += 1
    return codes, inv
CODES, INV = tail_codes()
def alpha(k):
    mu = bin(k).count("1") & 1; t = 0; level = None
    while (k >> t) >= 3:
        if ((k >> t) & 3) == 3: level = INV[k & ((1 << (t + 2)) - 1)]; break
        t += 1
    d = (k >> level) & 1 if level is not None and level <= max(1, k.bit_length()) else 2
    return 2 * d + mu
def mu_phase(bits, h):
    if h == 0: return 0
    pos = next((p for p in range(min(5, len(bits) - 1)) if bits[p] == bits[p + 1]), None)
    if pos is None: return None
    par = (1 - pos) & 1; q = mu_phase(bits[1 if par else 0::2], h - 1)
    return None if q is None else (2 * q - par) % (1 << h)
def residue(word, i=3, j=5):
    r = 1 << i; p = mu_phase([c & 1 for c in word], i)
    if p is None: return None
    t = (-p) % r; addr = 0
    for lv in range(i, j):
        s = word[t + CODES[lv]] // 2; addr += (0 if s == 2 else s) << lv
    return (addr - t) % (1 << j)
def reference_bad(c, m=32):
    n = len(c); R = [residue([c[(s + t) % n] for t in range(m)]) for s in range(n)]
    return [s for s in range(n) if R[s] is None or R[(s + 1) % n] is None or R[(s + 1) % n] != (R[s] + 1) % 32], R
def phases(c, W=20):
    n = len(c); return [mu_phase([c[(s + t) % n] & 1 for t in range(W)], 3) for s in range(n)]
def abit(s): return 0 if s >= 4 else s // 2
def g(L, lph, ph, lab, mut=None):
    l3, l4, cy, v3, v4, _ = L; a = abit(lab)
    ok = ph is not None and lph is not None and ph == (lph + 1) % 8
    if mut == "M_NOPOS": ok = ph is not None; 
    if not ok:
        v3 = v4 = 0
        if ph is None: return (0, 0, 0, 0, 0, 0)
    r = 0
    n3, n4, ncy, nv3, nv4 = l3, l4, cy, v3, v4
    if ph == 3:
        if v3 and a == l3: r = 1
        ncy = l3
        if mut != "M_NOLAUNCH": n3 = a; nv3 = 1
    if ph == 6:
        if v3 and v4 and a != (l4 ^ cy): r = 1
        if mut != "M_NOLAUNCH": n4 = a; nv4 = 1 if (v3 or mut == "M_NOPOS") else 0
    if mut == "M_NOPOS": nv3 = nv4 = 1
    return (n3, n4, ncy, nv3, nv4, r)
def step(X, c, PH, mut=None):
    n = len(c); return tuple(g(X[i-1], PH[i-1], PH[i], c[i], mut) for i in range(n))
ALL = list(itertools.product((0,1), repeat=6))
def cartesian(c, mut=None, T=200):
    n = len(c); PH = phases(c); S = [set(ALL) for _ in range(n)]
    for t in range(1, T + 1):
        S2 = [{g(x, PH[i-1], PH[i], c[i], mut) for x in S[i-1]} for i in range(n)]
        if S2 == S: return t - 1, S
        S = S2
    return None, S
def report(name, c, mut=None):
    bad, R = reference_bad(c); f, S = cartesian(c, mut)
    single = all(len(s) == 1 for s in S)
    rej = sorted(i for i in range(len(c)) if any(x[5] for x in S[i]))
    print(name, "n", len(c), "mutant", mut, "| reference bad starts", bad, "| flush", f, "unique fixed point", single,
          "| reject sites", rej, flush=True)
if __name__ == "__main__":
    A16 = [int(ch) for ch in "4554544154434512"]
    assert [alpha(k) for k in range(16)] == A16
    for mut in (None, "M_NOPOS", "M_NOLAUNCH"):
        report("alpha16-cycle", A16, mut)
        report("alpha64-cycle", [alpha(k) for k in range(64)], mut)
        report("alpha[100:164]+alpha[300:364] splice", [alpha(k) for k in range(100, 164)] + [alpha(k) for k in range(300, 364)], mut)
    c = [alpha(k) for k in range(1000, 1512)]
    f, S = cartesian(c); print("alpha[1000:1512] ring: flush", f, "reject sites", sorted(i for i in range(512) if any(x[5] for x in S[i])), flush=True)
