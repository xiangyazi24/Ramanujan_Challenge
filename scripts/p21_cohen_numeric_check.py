import mpmath as mp
mp.mp.dps = 120

def alpha(n): return 220*n**3 - 176*n**2 - 7*n + 5
def beta(n):  return 4*n**2*(2*n+1)**2*(5*n-4)*(5*n+6)
def a(n): return -220*n**3 - 484*n**2 - 301*n - 42
def b(n): return 4*n**2*(2*n+1)**2*(5*n-4)*(5*n+6)

# (0) the claimed index-shift identities
ok = all(a(n) == -alpha(n+1) and b(n) == beta(n) for n in range(0, 60))
print("index-shift identities a_n = -alpha(n+1), b_n = beta(n):", ok)

# (1) Cohen entry 5.3.22:  pi = 3 + 6/(alpha(1) + beta(1)/(alpha(2) + ...))
def cohen(N):
    t = mp.mpf(alpha(N))
    for n in range(N-1, 0, -1):
        t = alpha(n) + mp.mpf(beta(n))/t
    return 3 + mp.mpf(6)/t

# (2) the challenge PCF: a_0 + b_1/(a_1 + b_2/(a_2 + ...))
def challenge(N):
    t = mp.mpf(a(N))
    for n in range(N-1, 0, -1):
        t = a(n) + mp.mpf(b(n+1))/t     # careful: b_{n+1} sits above a_n
    return a(0) + mp.mpf(b(1))/t

print("\nCohen 5.3.22 vs pi:")
for N in (10, 20, 40, 80):
    v = cohen(N)
    print("   N=%3d  value=%s  err=%s" % (N, mp.nstr(v, 30), mp.nstr(v-mp.pi, 6)))

target = 6/(3-mp.pi)
print("\ntarget 6/(3-pi) =", mp.nstr(target, 30))
print("challenge PCF:")
for N in (10, 20, 40, 80):
    v = challenge(N)
    print("   N=%3d  value=%s  err=%s" % (N, mp.nstr(v, 30), mp.nstr(v-target, 6)))

# (3) golden-ratio rate check
phi = (1+mp.sqrt(5))/2
print("\nphi^-10 =", mp.nstr(phi**-10, 8), " -> digits/term =", mp.nstr(10*mp.log10(phi), 6))
