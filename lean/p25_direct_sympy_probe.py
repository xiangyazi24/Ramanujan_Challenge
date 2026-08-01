#!/usr/bin/env python3
"""Plain-SymPy construction and low-degree gauge scan for direct Catalan 3F2."""

import sympy as s

n=s.symbols('n')

def theta(pos):
    x0,x1,x2,y0,y1=pos
    c2=((y0-1)+(y1-1)+x0+x1+x2)/2
    c1=((y0-1)*(y1-1)+x0*x1+x0*x2+x1*x2)/2
    c0=x0*x1*x2/2
    return s.Matrix([[0,0,-c0],[1,0,-c1],[0,1,-c2]])

def xpos(pos,i): return s.eye(3)+theta(pos)/pos[i]

def xneg(pos,i):
    q=list(pos);q[i]-=1
    return xpos(q,i).inv()

def diag(pos):
    q=list(pos)
    out=xpos(q,1);q[1]+=1
    out=out*xneg(q,0);q[0]-=1
    return s.simplify(out),q

pos=[s.Rational(1,2)-2*n,s.Rational(1,2)+2*n,1,s.Rational(3,2),s.Rational(3,2)]
b1,pos1=diag(pos)
b2,pos2=diag(pos1)
B=s.simplify(b1*b2)
print('B')
for row in B.tolist(): print([s.factor(v) for v in row])
print('detB',s.factor(B.det()))
