from sage.all import *

# Exploratory normalization probe for the four d=1,2,3,6 Eisenstein-Eichler pieces.
N = 24
PS = PowerSeriesRing(QQ, 'q', default_prec=N+4)
q = PS.gen()

# A=q/t = product_{3 not | a} (1+q^a)^12.
A = PS(1)
for a in range(1, N+4):
    if a % 3 != 0:
        A *= (1 + q^a)^12
        A = A.add_bigoh(N+4)
tq = (q/A).add_bigoh(N+4)

# E = eta(2)^7 eta(3)^7 / (eta(1)^5 eta(6)^5), with q-powers cancelled.
E = PS(1)
for n in range(1, N+4):
    E *= (1-q^(2*n))^7 * (1-q^(3*n))^7
    E /= (1-q^n)^5 * (1-q^(6*n))^5
    E = E.add_bigoh(N+4)

# compositional inverse q=q(t) by Newton in a t-series ring
TS = PowerSeriesRing(QQ, 't', default_prec=N+2)
t = TS.gen()
tqT = TS([tq[i] for i in range(N+2)])
qt = t
for _ in range(8):
    err = tqT(qt) - t
    der = tqT.derivative()(qt)
    qt = (qt - err/der).add_bigoh(N+2)
assert (tqT(qt)-t).add_bigoh(N+1) == 0

Et = TS([E[i] for i in range(N+2)])(qt).add_bigoh(N+1)
b = [ZZ(Et[n]) for n in range(N+1)]
print('b', b[:8])
assert b[:6] == [1,5,73,1445,33001,819005]

# S_d(q)=sum sigma_3(m)/m^3 q^(d*m); K_d(t)=E(q) S_d(q), expressed in t.
kaps = {}
for d in [1,2,3,6]:
    Sd = PS(0)
    for m in range(1, (N+2)//d + 1):
        Sd += QQ(sigma(m,3), m^3) * q^(d*m)
    Kq = (E*Sd).add_bigoh(N+2)
    Kt = TS([Kq[i] for i in range(N+2)])(qt).add_bigoh(N+1)
    kaps[d] = [QQ(Kt[n]) for n in range(N+1)]
    print('d',d,'kap0..8',kaps[d][:9])

# W_r^(d)=r^3 (b_{r-1} kappa_r - b_r kappa_{r-1}).
Ws = {d:[QQ(0)]*(N+1) for d in kaps}
for d in kaps:
    for r in range(1,N+1):
        Ws[d][r] = QQ(r^3)*(QQ(b[r-1])*kaps[d][r] - QQ(b[r])*kaps[d][r-1])
    print('d',d,'W1..10',Ws[d][1:11])

for r in range(1,N+1):
    rel = Ws[1][r] - 28*Ws[2][r] + 63*Ws[3][r] - 36*Ws[6][r]
    print('REL',r,rel)

print('DONE')
