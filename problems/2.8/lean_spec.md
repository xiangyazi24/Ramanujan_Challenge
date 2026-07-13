# Lean Proof Spec for Problem 2.8

## Goal
Prove in Lean 4 (within the Ripple project) that the 4×4 CMF matrix product
from Problem 2.8 converges to √10005/π.

## What Ripple already has (0 sorry)
- `Ripple.Number.Chudnovsky1989.a` — coefficient sequence a_k = (6k)!/((3k)!(k!)³)
- `a_recurrence` — the 2-term recurrence for a_k
- `a_eq_3F2_coeff` — identification a_k = ₃F₂ coeff × 1728^k
- `chudnovskyF` — the ₃F₂(1/6,1/2,5/6;1,1;z) function
- `chudnovskyX` — evaluation point -1728/640320³
- `chudnovskyLinear3F2Series` — the full A+Bk series
- `chudnovskySeries` — the real-valued series
- `chudnovsky_one_over_pi` — final theorem (conditional on CM evaluation)
- All summability lemmas

## What needs to be added

### Step 1: Define the CMF matrix M(n) and initial matrix A
Define the 4×4 matrix M(n) with the polynomial entries from Problem 2.8.
Define the 2×4 initial matrix A with the large integer entries.
Define the product M_N = M(0)·...·M(N-1) and P_{N,j}/Q_{N,j} ratios.

### Step 2: Extract the scalar recurrence
Show that the (1,1) entry of A·M_N satisfies a scalar recurrence
equivalent (after gauge transformation) to the Chudnovsky recurrence
`a_recurrence`.

### Step 3: Match initial conditions
Show the initial conditions from A select the Chudnovsky partial sums.

### Step 4: Conclude
Apply `chudnovsky_one_over_pi` (conditional on CM) to get
P_{N,j}/Q_{N,j} → √10005/π.

## Key identity
R = 151931373056001 = 1 + 640320³/1728
√10005/π = (1/426880) × chudnovskySeries (after normalization)
The connection: 640320 = 64 × 10005, so √640320 = 8√10005.

## Priority
Step 2 is the core — showing the matrix recurrence = Chudnovsky recurrence.
This is a polynomial identity that can be verified by `ring` or `norm_num`
once both sides are expanded.
