# Machine check of Q6372 section 1-2: operator annihilates tau=sqrt(F); indicial data
from sympy import symbols, sqrt, series, Rational, expand, simplify, Poly, solve
from sympy import O as bigO
t = symbols('t')
N = 40
# b_n exact rationals via recurrence
b = [Rational(1), Rational(5)]
for n in range(1, N+2):
    b.append(Rational((34*n**3+51*n**2+27*n+5)*b[n] - n**3*b[n-1], (n+1)**3))
F = sum(b[n]*t**n for n in range(N+2))
# tau = sqrt(F) series
tau = [Rational(1)]
for k in range(1, N+1):
    acc = sum(tau[i]*tau[k-i] for i in range(1, k))
    tau.append((b[k]-acc)/2)
# (1) check half-index recurrence (1.1)
ok1 = all(4*(j+2)**2*tau[j+2] == 2*(68*j**2+170*j+107)*tau[j+1] - (2*j+1)**2*tau[j] for j in range(N-2))
# (2) check ODE 4t^2 q y'' + 4t(1-51t+2t^2) y' + t(t-10) y = 0 as series
y = sum(tau[n]*t**n for n in range(N+1))
q = 1 - 34*t + t**2
expr = expand(4*t**2*q*y.diff(t,2) + 4*t*(1-51*t+2*t**2)*y.diff(t) + t*(t-10)*y)
pe = Poly(expr, t)
ok2 = all(pe.coeff_monomial(t**k) == 0 for k in range(N-1))
# (3) indicial at root a of q: A1(a)/A2'(a) = 1/2 given q(a)=0
a = symbols('a')
A2p = (4*t**2*q).diff(t).subs(t, a)
A1a = (4*t*(1-51*t+2*t**2)).subs(t, a)
rel = simplify((A1a/A2p).subs(a**2, 34*a-1))  # reduce mod q(a)=0
# careful: substitute a^2 = 34a - 1 repeatedly
from sympy import cancel, together
r = cancel(together(A1a/A2p))
num, den = r.as_numer_denom()
red = lambda p: Poly(p, a).rem(Poly(a**2-34*a+1, a))
val = simplify(red(num).as_expr()/red(den).as_expr())
# evaluate numerically at a = 17+12*sqrt(2)
import mpmath
av = 17+12*mpmath.sqrt(2)
A2pv = 8*av**3 - 3*4*34*av**2/ (1)  # just numeric direct
f_A2 = lambda x: 8*x**3 - 408*x**2 + 8*x  # d/dt[4t^2(1-34t+t^2)] = 8t -408t^2+16t^3 -> recompute
# d/dt [4t^2 -136 t^3 + 4 t^4] = 8t - 408 t^2 + 16 t^3
A2pv = 8*av - 408*av**2 + 16*av**3
A1v = 4*av*(1-51*av+2*av**2)
print("half-index recurrence:", ok1)
print("ODE annihilates tau (coeffs to t^%d):"%(N-2), ok2)
print("A1(a)/A2'(a) =", A1v/A2pv, "(claim: 0.5)")
# (4) sigma recurrence (1.6)
g = [Rational(1), Rational(34)]
for k in range(2, N+1): g.append(34*g[k-1]-g[k-2])
c = [sum(b[i]*g[k-i] for i in range(k+1)) for k in range(N+1)]
sig = [Rational(1)]
for k in range(1, N+1):
    acc = sum(sig[i]*sig[k-i] for i in range(1, k))
    sig.append((c[k]-acc)/2)
ok4 = all(4*(j+2)**2*sig[j+2] == 2*(68*j**2+238*j+209)*sig[j+1] - (2*j+3)**2*sig[j] for j in range(N-2))
print("sigma half-index recurrence (1.6):", ok4)
