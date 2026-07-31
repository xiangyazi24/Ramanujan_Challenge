from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 120

def coeffs(n):
    c0 = -8*n**3 - 51*n**2 - 105*n - 68
    c1 = 24*n**5 + 337*n**4 + 1833*n**3 + 4818*n**2 + 6092*n + 2928
    c2 = -(n+2)*(n+3)*(24*n**5 + 273*n**4 + 1150*n**3 + 2154*n**2 + 1635*n + 268)
    c3 = (n+1)*(n+2)**4*(n+3)*(8*n**3 + 75*n**2 + 231*n + 232)
    return c0,c1,c2,c3

def solve(init, N):
    u = {-3:F(init[0]), -2:F(init[1]), -1:F(init[2])}
    for n in range(0, N+1):
        c0,c1,c2,c3 = coeffs(n)
        assert c0 != 0, n
        u[n] = -(c1*u[n-1] + c2*u[n-2] + c3*u[n-3])/c0
    return u

N=60
p = solve([0,7,179], N); q = solve([1,12,306], N)
print("integral p:", all(v.denominator==1 for v in p.values()))
print("integral q:", all(v.denominator==1 for v in q.values()))
print("first values q:", [int(q[n]) for n in range(-3,4)])
print("first values p:", [int(p[n]) for n in range(-3,4)])
g = mp.euler
for n in [5,10,20,30,40,50,60]:
    r = mp.mpf(int(p[n]))/mp.mpf(int(q[n]))
    print("n=%2d  p/q-gamma = %s   digits=%.1f" % (n, mp.nstr(r-g,6), -mp.log10(abs(r-g)) if r!=g else 999))
# growth: q_n/(n!)^2
print("\nq_n / ((n+3)!)^2 :")
for n in [5,10,20,30]:
    print("  n=%2d  %s" % (n, mp.nstr(mp.mpf(int(q[n]))/mp.factorial(n+3)**2, 12)))
