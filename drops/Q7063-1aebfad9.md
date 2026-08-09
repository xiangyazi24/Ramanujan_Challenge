ANSWER Q7063 1aebfad9

# Transversality of the Apéry polynomial family and CRT repulsion

## Short answer

The proposed geometric intuition is correct at the level of philosophy, but pairwise coprimality of the one-parameter family `b_r(λ)` is not by itself strong enough to imply the desired bound on primes hitting a fixed integer `m`. The reason is that the obstruction is not a common component of the hypersurface in the `(r,λ)` direction; it is the arithmetic freedom coming from choosing different residue characteristics `p_i` and different points `(r_i,1)`. A genuine unlikely-intersection theorem would need an additional height statement controlling the arithmetic complexity of the collection of points.

The right framework is not simply the variety `V`, but an incidence variety with the CRT parameter `m` and with the primes as varying residue characteristics.

## 1. Reformulating the event geometrically

For a fixed integer `m` and primes `p_i` in `(X,2X]`, write

```
r_i = m mod p_i.
```

The bad-prime condition is

```
b_{r_i}(1) = 0 mod p_i.
```

Equivalently, we have points

```
(r_i, λ=1) ∈ A^1 × F_{p_i}
```

lying on the reduction of the corresponding polynomial divisor. The CRT condition says

```
m ≡ r_i (mod p_i).
```

Thus the tuple `(p_i,r_i)` is constrained by two independent-looking equations:

1. the local geometric condition `b_{r_i}(1)=0 mod p_i`;
2. the compatibility condition that all `r_i` arise from one integer `m`.

The desired conclusion would be that many simultaneous local intersections force an impossible global intersection.

## 2. What pairwise coprimality really gives

The hypothesis

```
Res_λ(b_r,b_s) ≠ 0
```

means that for fixed characteristic zero there is no algebraic number `λ` satisfying infinitely many equations

```
b_r(λ)=0.
```

It rules out a horizontal component of the universal zero locus. In geometric language, the family is transverse in the `λ`-direction.

However, our problem specializes to

```
λ = 1
```

and then reduces modulo different primes. Pairwise coprimality over `Q` does not prevent infinitely many unrelated reductions from vanishing at `λ=1`.

A simple analogy: a sequence of coprime integer polynomials can still have many different prime divisors at a fixed integer argument. Coprimality controls common algebraic zeros, not the distribution of residue characteristics.

Therefore the missing ingredient is not transversality alone but arithmetic transversality.

## 3. What an unlikely-intersection theorem would need to say

A Bombieri--Masser--Zannier/Zannier style argument would require an ambient algebraic group or a controlled family where anomalous intersections have smaller dimension. Here the prime coordinate is not an algebraic parameter in a fixed variety; it is a varying residue field.

The natural object would be something like

```
W = { (m,p,r) : m-r = 0 mod p, b_r(1)=0 mod p }.
```

A tuple of `T` bad primes corresponds to an intersection of `T` Frobenius-type conditions on `W`.

To obtain a repulsion theorem one would need a statement of the form:

```
height(m) large + many independent reductions on V
        => arithmetic degree/complexity contradiction.
```

This is closer to an arithmetic unlikely-intersection theorem than to the usual geometric one.

## 4. Relation with Chabauty-Coleman

Chabauty works because one has a single global curve `C/Q` and asks for rational points satisfying local constraints. The key input is a global differential annihilating all rational points.

Here the roles are reversed:

- the index `r` varies;
- the prime varies;
- the residue condition changes with the characteristic.

There is no single curve over `Q` whose rational points encode all bad primes.

A possible Chabauty analogue would require embedding the Apéry sequence into a fixed algebraic dynamical object (for example the projective orbit of the Apéry recurrence) and proving that the zero condition defines a sparse set of reductions. The existing transfer-matrix viewpoint is closer to this direction.

## 5. A more realistic geometric route

The useful object is probably not the zero variety itself, but the correlation varieties:

```
b_r(1)=0 mod p,
b_{r+h}(1)=0 mod p.
```

or across different primes:

```
b_{m mod p}(1)=0 mod p,
b_{m mod q}(1)=0 mod q.
```

The existing project results already exploit this philosophy:

- no-consecutive-zero structure;
- reflection symmetry;
- gap polynomial bounds;
- codegree bounds for common bad primes.

These are arithmetic intersection estimates. They show that intersections of the zero locus with shifted copies have unexpectedly small size.

A possible next step would be to prove a higher-order incidence bound:

```
#{(p1,...,pT): p_i ~ X,
 b_{m mod p_i}=0 mod p_i}
```

using a determinant/resultant construction in several variables. This would be the closest analogue of an unlikely-intersection theorem.

## 6. Conclusion

The transversality heuristic is pointing at the correct phenomenon:

> many primes forcing the same integer `m` to land on Apéry zero fibers should be arithmetically unlikely.

But the statement

```
Res_λ(b_r,b_s) ≠ 0
```

only removes global components and does not by itself create quantitative repulsion.

The missing theorem would be an arithmetic unlikely-intersection statement for the family of reductions of the Apéry zero variety, probably involving:

1. a higher-dimensional incidence variety including `m` and the prime reductions;
2. height bounds for simultaneous intersections;
3. a quantitative estimate on the codimension gain from each additional prime.

In the language of the current Apéry GCD program, this would be a replacement for the unavailable square-root bound on `|Z_p|`: instead of proving each fiber is small, prove that a single CRT point cannot lie in too many fibers simultaneously.