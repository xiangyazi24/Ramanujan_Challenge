#!/usr/bin/env python3
"""
P2.5: AZ certificate search at eps != 0 to capture the Catalan direction.

The eps=0 certificate (verified) gives only the trivial Delannoy-square solution.
For G to appear, we need the t-connection pole at t=0, which enters at eps != 0.

The t-connection D_t = d_t + eps/t - 2t/(1+t^2) has a simple pole at t=0
when eps != 0. This means:
  - C_n certificates may need NEGATIVE t powers
  - The integral picks up a residue at t=0, bringing in arctan(1) = pi/4 or G

Search plan:
  1. eps=1 with polynomial P support [0,2]^3, varying n-degrees
  2. eps=1 with extended t support [-1,2] or [-2,2]
  3. eps=generic (random prime-field element)
"""

import os
import sys
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p25_az_correct import solve_config, PRIME1, PRIME2, describe_P

def search(dP, dC, eps_mode, Pu=(0,2), Pv=(0,2), Pt=(0,2),
           Ct=None, label=""):
    print(f"\n{'='*70}")
    print(f"Search: dP={dP}, dC={dC}, eps={eps_mode}, "
          f"Pu={Pu}, Pv={Pv}, Pt={Pt}, Ct_override={Ct}  [{label}]")
    print(f"{'='*70}", flush=True)

    kwargs = dict(Pu=Pu, Pv=Pv, Pt=Pt)
    if Ct is not None:
        kwargs['Ct'] = Ct

    res = solve_config(dP, dC, eps_mode, PRIME1, seed=42, **kwargs)

    if res['Pdim'] > 0 and res['verify_ok']:
        print(f"  >>> P-dim = {res['Pdim']} with verify=OK!")
        res2 = solve_config(dP, dC, eps_mode, PRIME2, seed=42, **kwargs)
        if res2['Pdim'] == res['Pdim']:
            print(f"  >>> Two-prime agreement: P-dim = {res['Pdim']}")
            for line in describe_P(res, res['Pbasis'][0]):
                print(line)
        else:
            print(f"  WARNING: p1 P-dim={res['Pdim']}, p2 P-dim={res2['Pdim']}")
        return True
    elif res['Pdim'] > 0 and not res['verify_ok']:
        print(f"  P-dim={res['Pdim']} but verify=FAIL (underdetermined)")
    else:
        print(f"  P-dim=0 (no certificate)")
    return False


print("=" * 70)
print("P2.5 AZ certificate search: eps != 0")
print("Looking for Catalan-direction certificates")
print("=" * 70)

# Phase 1: eps=1 with standard polynomial support, varying degrees
for dP in range(3):
    dC = dP
    if search(dP, dC, '1', label=f"eps=1 poly, degree {dP}"):
        break

# Phase 2: eps=1 with extended t-support for P (allow t^{-1})
for dP in range(3):
    dC = dP
    if search(dP, dC, '1', Pt=(-1, 2), Ct=(-1, 3),
              label=f"eps=1 extended-t, degree {dP}"):
        break

# Phase 3: eps=1 with Laurent t-support [-2,2]
for dP in range(3):
    dC = dP
    if search(dP, dC, '1', Pt=(-2, 2), Ct=(-2, 3),
              label=f"eps=1 full-Laurent-t, degree {dP}"):
        break

# Phase 4: eps=2 (another specific value)
for dP in range(3):
    dC = dP
    if search(dP, dC, '2', label=f"eps=2 poly, degree {dP}"):
        break

# Phase 5: eps=generic
for dP in range(3):
    dC = dP
    if search(dP, dC, 'generic', label=f"eps=generic poly, degree {dP}"):
        break

# Phase 6: eps=1, wider u,v support [-1,2] x [-1,2]
for dP in range(3):
    dC = dP
    if search(dP, dC, '1', Pu=(-1, 2), Pv=(-1, 2), Pt=(-1, 2), Ct=(-1, 3),
              label=f"eps=1 wider-uv-t, degree {dP}"):
        break

# Phase 7: eps=free (eps as free variable)
for dP in range(3):
    dC = dP
    if search(dP, dC, 'free', Pt=(-1, 2), Ct=(-1, 3),
              label=f"eps=free extended-t, degree {dP}"):
        break

print("\n\nDone.")
