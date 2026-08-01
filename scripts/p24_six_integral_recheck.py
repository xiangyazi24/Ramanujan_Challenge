from mpmath import mp, mpf, log, polylog, zeta, quad, pi, exp, nstr
mp.dps = 40
Z2 = pi**2/6; L = log(2); A4 = polylog(4, mpf(1)/2); Z3 = zeta(3)
def W0(t): return Z2 - 2*polylog(2, t/2) - log(t/2)**2
def H1(t): return -log(1-t)
def H2(t): return -log(1-t/2)
# substitutions kill the endpoint singularities: t=e^-u near 0, 1-t=e^-u near 1
def I(num, den):
    f = lambda t: W0(t)*num(t)/den(t)
    a = quad(lambda u: f(exp(-u))*exp(-u), [mp.log(2), 5, 40])       # t in (0,1/2]
    b = quad(lambda u: f(1-exp(-u))*exp(-u), [mp.log(2), 5, 40])     # t in [1/2,1)
    return a + b
defs = {'I10':(H1,lambda t:t), 'I11':(H1,lambda t:1-t), 'I12':(H1,lambda t:2-t),
        'I20':(H2,lambda t:t), 'I21':(H2,lambda t:1-t), 'I22':(H2,lambda t:2-t)}
rows = {'I10':(-2,mpf(-1)/12,-1,mpf(-7)/4,mpf(1)/10),
        'I11':(0,0,0,mpf(-7)/2,mpf(3)/4),
        'I12':(-6,mpf(-1)/4,3,mpf(-21)/4,mpf(9)/5),
        'I20':(-2,mpf(-1)/12,mpf(1)/2,mpf(-7)/4,mpf(1)/4),
        'I21':(-6,mpf(-1)/4,mpf(-3)/2,mpf(-7)/2,mpf(51)/20),
        'I22':(-6,mpf(-1)/4,mpf(3)/2,mpf(-21)/4,mpf(23)/10)}
basis = [A4, L**4, L**2*Z2, L*Z3, Z2**2]
c = {'I10':-2,'I11':-2,'I12':2,'I20':4,'I21':6,'I22':-5}
tot_n = mpf(0); tot_c = mpf(0)
print(f"{'row':<5}{'numeric':>26}{'closed form':>26}   |diff|")
for k in ['I10','I11','I12','I20','I21','I22']:
    n = I(*defs[k]); cf = sum(a*b for a,b in zip(rows[k], basis))
    tot_n += c[k]*n; tot_c += c[k]*cf
    print(f"{k:<5}{nstr(n,22):>26}{nstr(cf,22):>26}   {nstr(abs(n-cf),3)}")
tgt = -22*A4 - mpf(11)/12*L**4 - mpf(13)/2*L**2*Z2 - mpf(7)/4*L*Z3 + mpf(67)/10*Z2**2
print(f"\ncombination (numeric integrals) = {nstr(tot_n,25)}")
print(f"combination (closed forms)      = {nstr(tot_c,25)}")
print(f"repo alternatingQuadEulerValue24= {nstr(tgt,25)}")
print(f"|numeric - repo target|         = {nstr(abs(tot_n-tgt),3)}")
