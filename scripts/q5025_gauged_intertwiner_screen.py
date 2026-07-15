#!/usr/bin/env python3
"""Screen polynomial intertwiners after the necessary 64^{-n} gauge.

If a_n is annihilated by LZ and f_n = a_n/64^n, then f is annihilated by
    LZ64 = sum_j 64^j z_j(n+2) S^j.
This has characteristic polynomial P27(64*r), matching L27.
"""
from __future__ import annotations
from typing import List, Optional, Sequence, Tuple

from q5025_poly_intertwiner_screen import (
    l27_at, lz_at, rank_mod
)


def lz64_at(n: int, p: int) -> List[int]:
    z = lz_at(n,p)
    return [(z[j] * pow(64,j,p)) % p for j in range(4)]


def basis_remainder_scaled(i: int, k: int, n0: int, p: int) -> Optional[Tuple[int,int,int]]:
    l27 = l27_at(n0,p)
    coeff = [0]*(i+4)
    for j in range(4):
        coeff[i+j] = l27[j]*pow((n0+j)%p,k,p)%p
    for m in range(i+3,2,-1):
        shift=m-3
        z=lz64_at(n0+shift,p)
        if z[3]==0:
            return None
        q=coeff[m]*pow(z[3],p-2,p)%p
        if q:
            for j in range(4):
                coeff[shift+j]=(coeff[shift+j]-q*z[j])%p
    return coeff[0],coeff[1],coeff[2]


def matrix_for(d: int,D: int,p: int,extra_blocks: int=12):
    cols=[(i,k) for i in range(d+1) for k in range(D+1)]
    need=(len(cols)+2)//3+extra_blocks
    rows=[]; ns=[]; n0=0
    while len(ns)<need and n0<30*need+100:
        rr=[]; ok=True
        for i,k in cols:
            v=basis_remainder_scaled(i,k,n0,p)
            if v is None:
                ok=False; break
            rr.append(v)
        if ok:
            for c in range(3):
                rows.append([v[c] for v in rr])
            ns.append(n0)
        n0+=1
    return rows,ns


def screen(max_d=6,max_D=40):
    primes=[1000003,1000033,1000037]
    print('Q5025 corrected LZ64 polynomial-intertwiner screen')
    print('LZ64 coefficients are 64^j*z_j(n+2); f_n=a_n/64^n')
    survivors=[]
    for d in range(max_d+1):
        print(f'=== order d={d} ===')
        for DD in range(max_D+1):
            u=(d+1)*(DD+1)
            ranks=[]; full=None
            for p in primes:
                rows,ns=matrix_for(d,DD,p)
                rk=rank_mod(rows,p); ranks.append(rk)
                if rk==u:
                    full=p; break
            if full:
                print(f'd={d:2d} D={DD:2d} unknowns={u:3d}: FULL rank mod {full}')
            else:
                print(f'd={d:2d} D={DD:2d} unknowns={u:3d}: DEFICIENT ranks={ranks}')
                survivors.append((d,DD,u,ranks))
        print()
    print('SURVIVORS',survivors)

if __name__=='__main__':
    screen()
