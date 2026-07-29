#!/usr/bin/env python3
"""Verify the corrected Seifert presentation and trace certificate for P3.1.

Corrected presentation (Q4833/Q4844): M(-1; (2,1),(3,1),(17,3))
  x²h = 1,  y³h = 1,  z¹⁷h³ = 1,  xyz = h⁻¹
  ρ(h) = -I

Trace certificate: tr(X) = 0, tr(Y) = +1, tr(XY) = -2cos(π/17)

Explicit matrices (Q4844 eq 1.2):
  X₀ = [[0, -1], [1, 0]]
  Y₀ = [[1/2, d-c], [d+c, 1/2]]  where c = cos(π/17), d = √(c²-3/4)
"""
import mpmath as mp

mp.mp.dps = 80

c = mp.cos(mp.pi / 17)
d = mp.sqrt(c**2 - mp.mpf(3)/4)
print(f"c = cos(π/17) = {mp.nstr(c, 30)}")
print(f"d = √(c²-3/4) = {mp.nstr(d, 30)}")

def tr(M):
    return M[0,0] + M[1,1]

X = mp.matrix([[0, -1], [1, 0]])
Y = mp.matrix([[mp.mpf(1)/2, d - c], [d + c, mp.mpf(1)/2]])

print(f"\nX = [[{mp.nstr(X[0,0],5)}, {mp.nstr(X[0,1],5)}], [{mp.nstr(X[1,0],5)}, {mp.nstr(X[1,1],5)}]]")
print(f"Y = [[{mp.nstr(Y[0,0],5)}, {mp.nstr(Y[0,1],10)}], [{mp.nstr(Y[1,0],10)}, {mp.nstr(Y[1,1],5)}]]")

print(f"\ndet(X) = {mp.nstr(X[0,0]*X[1,1] - X[0,1]*X[1,0], 10)}")
print(f"det(Y) = {mp.nstr(Y[0,0]*Y[1,1] - Y[0,1]*Y[1,0], 10)}")

print(f"\ntr(X) = {mp.nstr(tr(X), 10)}")
print(f"tr(Y) = {mp.nstr(tr(Y), 10)}")

XY = X * Y
print(f"tr(XY) = {mp.nstr(tr(XY), 20)}")
print(f"-2cos(π/17) = {mp.nstr(-2*c, 20)}")
print(f"|tr(XY) + 2cos(π/17)| = {mp.nstr(abs(tr(XY) + 2*c), 10)}")

# Check X² = -I
X2 = X * X
print(f"\nX² = [[{mp.nstr(X2[0,0],5)}, {mp.nstr(X2[0,1],5)}], [{mp.nstr(X2[1,0],5)}, {mp.nstr(X2[1,1],5)}]]")
ok = abs(X2[0,0] + 1) < 1e-50 and abs(X2[1,1] + 1) < 1e-50
print(f"X² = -I? {ok}")

# Check Y³ = -I
Y2 = Y * Y
Y3 = Y2 * Y
print(f"\nY³ = [[{mp.nstr(Y3[0,0],10)}, {mp.nstr(Y3[0,1],10)}], [{mp.nstr(Y3[1,0],10)}, {mp.nstr(Y3[1,1],10)}]]")
ok = abs(Y3[0,0] + 1) < 1e-50 and abs(Y3[1,1] + 1) < 1e-50 and abs(Y3[0,1]) < 1e-50 and abs(Y3[1,0]) < 1e-50
print(f"Y³ = -I? {ok}")

# Z = -(XY)⁻¹
Z = -mp.inverse(XY)
print(f"\ntr(Z) = {mp.nstr(tr(Z), 20)}")
print(f"2cos(π/17) = {mp.nstr(2*c, 20)}")
print(f"|tr(Z) - 2cos(π/17)| = {mp.nstr(abs(tr(Z) - 2*c), 10)}")

# Check Z¹⁷ = -I
Zpow = mp.eye(2)
for _ in range(17):
    Zpow = Zpow * Z
print(f"\nZ¹⁷ = [[{mp.nstr(Zpow[0,0],10)}, {mp.nstr(Zpow[0,1],10)}], [{mp.nstr(Zpow[1,0],10)}, {mp.nstr(Zpow[1,1],10)}]]")
ok = abs(Zpow[0,0] + 1) < 1e-50 and abs(Zpow[1,1] + 1) < 1e-50
print(f"Z¹⁷ = -I? {ok}")

# Check (XY)¹⁷ = I
XYpow = mp.eye(2)
for _ in range(17):
    XYpow = XYpow * XY
print(f"\n(XY)¹⁷ = [[{mp.nstr(XYpow[0,0],10)}, {mp.nstr(XYpow[0,1],10)}], [{mp.nstr(XYpow[1,0],10)}, {mp.nstr(XYpow[1,1],10)}]]")
ok = abs(XYpow[0,0] - 1) < 1e-50 and abs(XYpow[1,1] - 1) < 1e-50
print(f"(XY)¹⁷ = I? {ok}")

# Check XYZ = -I
XYZ = X * Y * Z
print(f"\nXYZ = [[{mp.nstr(XYZ[0,0],10)}, {mp.nstr(XYZ[0,1],10)}], [{mp.nstr(XYZ[1,0],10)}, {mp.nstr(XYZ[1,1],10)}]]")
ok = abs(XYZ[0,0] + 1) < 1e-50 and abs(XYZ[1,1] + 1) < 1e-50
print(f"XYZ = -I? {ok}")

# Commutator trace (non-elementary check)
comm = X * Y * mp.inverse(X) * mp.inverse(Y)
tr_comm = tr(comm)
print(f"\ntr[X,Y] = {mp.nstr(tr_comm, 20)}")
print(f"4c² - 1 = {mp.nstr(4*c**2 - 1, 20)}")
print(f"tr[X,Y] > 2? {float(tr_comm) > 2}")

# Euler number and GV
print("\n=== Topological invariants ===")
e = mp.mpf(-1) + mp.mpf(1)/2 + mp.mpf(1)/3 + mp.mpf(3)/17
chi = mp.mpf(-1) + mp.mpf(1)/2 + mp.mpf(1)/3 + mp.mpf(1)/17
print(f"e = {mp.nstr(e, 20)} = {mp.nstr(e*102, 5)}/102")
print(f"χ_orb = {mp.nstr(chi, 20)} = {mp.nstr(chi*102, 5)}/102")
GV = 4 * mp.pi**2 * chi**2 / e
target = 242 * mp.pi**2 / 51
print(f"GV(Fuchsian) = {mp.nstr(GV, 20)}")
print(f"242π²/51 = {mp.nstr(target, 20)}")
print(f"|GV - 242π²/51| = {mp.nstr(abs(GV - target), 10)}")

# Target identity check
integral_target = 4 * mp.pi**2 / 85
print(f"\n4π²/85 = {mp.nstr(integral_target, 20)}")
gv_alpha = GV - integral_target
print(f"GV(α) = GV(β) - 4π²/85 = {mp.nstr(gv_alpha, 20)}")
print(f"1198π²/255 = {mp.nstr(1198 * mp.pi**2 / 255, 20)}")
print(f"|GV(α) - 1198π²/255| = {mp.nstr(abs(gv_alpha - 1198*mp.pi**2/255), 10)}")
