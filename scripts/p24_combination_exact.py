from fractions import Fraction as F
rows = {  # (A4, L^4, L^2*Z2, L*Z3, Z2^2)  -- from DOCTRINE_P24_LAYERD.md:71-76
 'I10': (F(-2), F(-1,12), F(-1),    F(-7,4),  F(1,10)),
 'I11': (F(0),  F(0),     F(0),     F(-7,2),  F(3,4)),
 'I12': (F(-6), F(-1,4),  F(3),     F(-21,4), F(9,5)),
 'I20': (F(-2), F(-1,12), F(1,2),   F(-7,4),  F(1,4)),
 'I21': (F(-6), F(-1,4),  F(-3,2),  F(-7,2),  F(51,20)),
 'I22': (F(-6), F(-1,4),  F(3,2),   F(-21,4), F(23,10)),
}
c = {'I10':-2,'I11':-2,'I12':2,'I20':4,'I21':6,'I22':-5}
comb = [sum(c[k]*rows[k][i] for k in rows) for i in range(5)]
target = (F(-22), F(-11,12), F(-13,2), F(-7,4), F(67,10))  # Problem24.lean:2215
names = ['Li4(1/2)','log^4 2','log^2 2*Z2','log2*Z3','Z2^2']
print(f"{'basis':<12}{'combination':>12}{'target':>12}   match")
for n,a,b in zip(names,comb,target):
    print(f"{n:<12}{str(a):>12}{str(b):>12}   {'OK' if a==b else 'MISMATCH'}")
print("\nALL MATCH" if list(comb)==list(target) else "\nFAILED")
