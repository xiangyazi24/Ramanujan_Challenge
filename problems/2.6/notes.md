# Problem 2.6: A Series for ζ(2) + ζ(3)

## Statement

Let (u_n)_{n≥1} with:
- u_1 = -93/4480
- u_2 = -117/14000

and for n ≥ 3, the recurrence:
```
0 = -2(n+3)³(2n+5)(3n+5) u_n
  + (n+2)²(15n³ + 85n² + 155n + 93) u_{n-1}
  - (n+1)³(n+2)(3n+8) u_{n-2}
```

Prove:
```
2077/720 + ∑_{j=1}^∞ u_j = ζ(2) + ζ(3)
```

## Contributor
Hila Barkan (Ramanujan Machine Group, Technion)

## Approach Ideas
- 3-term recurrence, so the u_n satisfy a holonomic equation
- The target ζ(2) + ζ(3) = π²/6 + ζ(3) ≈ 2.8469...
- This is in the Apéry tradition — find the integral representation
- Connection to HolonomicCRN project: holonomic sequences converging to periods
- Check if the recurrence factors or relates to known Apéry-like recurrences
- Creative telescoping with Zeilberger's algorithm
- Possible connection to Beukers-type integrals
