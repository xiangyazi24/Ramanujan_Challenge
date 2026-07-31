from fractions import Fraction as F

# ---- challenge (index m = n+3, so Q_0..Q_2 are the given initial values)
def ch_coeffs(n):
    return (-8*n**3-51*n**2-105*n-68,
            24*n**5+337*n**4+1833*n**3+4818*n**2+6092*n+2928,
            -(n+2)*(n+3)*(24*n**5+273*n**4+1150*n**3+2154*n**2+1635*n+268),
            (n+1)*(n+2)**4*(n+3)*(8*n**3+75*n**2+231*n+232))
def ch_solve(init, N):
    u={-3:F(init[0]),-2:F(init[1]),-1:F(init[2])}
    for n in range(0,N+1):
        c0,c1,c2,c3=ch_coeffs(n)
        u[n]=-(c1*u[n-1]+c2*u[n-2]+c3*u[n-3])/c0
    return [u[n] for n in range(-3,N+1)]     # indexed by m=n+3 -> position

# ---- Aptekarev (3): (16n-15)y_{n+1} = (128n^3+40n^2-82n-45)y_n
#                      - n^2(256n^3-240n^2+64n-7)y_{n-1} + n^2(n-1)^2(16n+1)y_{n-2}
def apt_solve(init, N):
    y={0:F(init[0]),1:F(init[1]),2:F(init[2])}
    for n in range(2,N):
        lead=16*n-15
        y[n+1]=(F((128*n**3+40*n**2-82*n-45))*y[n]
                - F(n**2*(256*n**3-240*n**2+64*n-7))*y[n-1]
                + F(n**2*(n-1)**2*(16*n+1))*y[n-2])/lead
    return [y[i] for i in range(N+1)]

# ---- Rivoal: (n+3)^2(8n+11)(8n+19) y_{n+3} = (n+3)(8n+11)(24n^2+145n+215) y_{n+2}
#              - (8n+27)(24n^3+105n^2+124n+25) y_{n+1} + (n+2)^2(8n+19)(8n+27) y_n
def riv_solve(init, N):
    y={0:F(init[0]),1:F(init[1]),2:F(init[2])}
    for n in range(0,N-2):
        lead=F((n+3)**2*(8*n+11)*(8*n+19))
        y[n+3]=(F((n+3)*(8*n+11)*(24*n**2+145*n+215))*y[n+2]
                - F((8*n+27)*(24*n**3+105*n**2+124*n+25))*y[n+1]
                + F((n+2)**2*(8*n+19)*(8*n+27))*y[n])/lead
    return [y[i] for i in range(N+1)]

N=10
chQ=ch_solve([1,12,306],N); chP=ch_solve([0,7,179],N)
aptQ=apt_solve([1,3,50],N); aptP=apt_solve([0,2,31],N)
rivQ=riv_solve([1,7,F(65,2)],N); rivP=riv_solve([-1,4,F(77,4)],N)

print("challenge Q:", [str(v) for v in chQ[:8]])
print("Aptekarev Q:", [str(v) for v in aptQ[:8]])
print("Rivoal    Q:", [str(v) for v in rivQ[:8]])
print()
print("challenge P:", [str(v) for v in chP[:8]])
print("Aptekarev P:", [str(v) for v in aptP[:8]])
print("Rivoal    P:", [str(v) for v in rivP[:8]])
print()
print("ratio challenge/Rivoal (Q):", [str(chQ[i]/rivQ[i]) for i in range(8)])
print("ratio challenge/Aptekarev (Q):", [str(chQ[i]/aptQ[i]) for i in range(8)])
