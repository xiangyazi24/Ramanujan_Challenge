ANSWER Q6543 3b0405d8

# P3.2 Orbit rigidity survey

## 1. Literature audit

The short verdict is: existing orbit-intersection theorems are mostly autonomous. They do not directly cover the Apery cocycle.

### Ostafe-Shparlinski

Relevant works:

- A. Ostafe and I. E. Shparlinski, "On the length of critical orbits of rational functions modulo p" (2013).
- A. Ostafe and I. E. Shparlinski, "Orbits of polynomial dynamical systems modulo primes" (Math. Ann. 2014).

The main results concern an autonomous map f in F_p[X]. Typical statement shape:

    f^n(x) = a

has few solutions when the orbit length is large and f is not dynamically degenerate.

The hypotheses are:

    x_{n+1}=f(x_n)

with one fixed polynomial/rational function f. The Apery recurrence has

    x_{n+1}=M_n(x_n)

where M_n changes with n. Therefore the theorem does not apply.

### Chang

M.-C. Chang, works on intersections of polynomial orbits and algebraic sets, proves bounds for autonomous polynomial dynamics. The mechanism is algebraic independence of iterates

    f(x), f^2(x), ...

and degree growth.

Again the iteration identity is lost for a nonautonomous cocycle.

### Kurlberg-Rudnick

Kurlberg-Rudnick results on cat maps and quantum chaos concern value distribution of autonomous linear dynamics on finite tori. The relevant mechanism is spectral expansion of a fixed operator.

Not applicable: the Apery transfer matrices are time dependent.

### Gomez-Nicolas-Ostafe-Sadornil

Their results on polynomial/rational dynamics modulo p again assume a fixed dynamical system. They give exceptional-set bounds and orbit statistics under non-degeneracy assumptions.

No theorem found giving the required nonautonomous Mobius cocycle rigidity.

## 2. What the Apery Casoratian gives

The cocycle is

    pi(n+1)=M_n(pi(n))

with

    M_n = [[P(n)/(n+1)^3, -n^3/(n+1)^3], [1,0]].

The identity

    b_{n+1} c_n - b_n c_{n+1} = -1/n^3

is a very strong non-degeneracy statement.

In autonomous terminology it plays the role of:

1. separability: distinct nearby states cannot collapse;
2. non-isotriviality: the dynamics does not become a fixed map after conjugation;
3. absence of invariant algebraic subvarieties.

However, it is not enough by itself. Autonomous proofs use repeated composition

    f^m

and degree growth. Here the composition

    M_{n+m-1} ... M_n

has arithmetic variation with n. A new theorem would need to exploit this variation.

Potential advantage: the explicit Wronskian may make the cocycle easier than a generic nonautonomous system because every two-step collision gives an explicit rational equation.

## 3. Three target rigidity statements

Let

    Z_v(T) = #{0<=n<p-1 : y^(v)(n)=0}

and m_v=Z_v(p-1).

### Level A (minimum needed)

For some delta>0:

    #{v : m_v >= p^(1/3)} <= p^(1-delta).

This only removes a positive-density exceptional family of bad fibers.

It is enough for a weak fiber-tail estimate:

    N(T) << p^(1-delta)

at T=p^(1/3).

### Level B (incidence useful)

For every T>=p^(1/3):

    N(T) <= p/T^2 * p^(epsilon)

for any epsilon>0.

Equivalently the second moment satisfies:

    sum_v m_v^2 <= p^(1+epsilon).

This is the natural threshold because

    sum_v m_v = p

and Cauchy converts second moment control into fiber tails.

### Level C (full rigidity)

Uniformly for all v:

    m_v <= (log p)^C.

This gives near-random fibers and strongest collision estimates, but is much stronger than required.

## 4. Why nonautonomous is the real barrier

A tempting transfer-operator approach would define a graph on states P^1(F_p):

    state -> M_n(state).

For a fixed graph one studies the second eigenvalue. Here the graph changes every time step.

A naive automaton has p states and p time labels. The state count itself is not fatal for second moments, but it kills the usual spectral gap argument because there is no stationary transition operator.

The right object would be a skew-product transfer operator:

    (n,x) -> (n+1,M_n(x)).

The clock variable has size p, so the spectrum is on a p*p space. Any useful theorem must exploit the explicit Apery structure instead of generic expansion.

## 5. Practical computational test

For a prime p:

1. Compute b_n,c_n modulo p for 0<=n<p-1.
2. For every v=(alpha,beta) in F_p^2 up to scalar, compute

       y_v(n)=beta*b_n-alpha*c_n.

3. Count zeros m_v.
4. Compare

       sum_v m_v^2

against the random prediction p^2.

Pseudo-code:

```python
for p in [13,17]:
    b,c=aper_y_sequences_mod_p(p)
    fibers=[]
    for a in range(p):
        for d in range(p):
            if (a,d)!=(0,0):
                m=sum((d*b[n]-a*c[n])%p==0 for n in range(p-1))
                fibers.append(m)
    print(max(fibers))
    print(sum(x*x for x in fibers))
```

The key diagnostic is not max(m_v), but the second moment. A random model predicts

    E[m_v]=1

and

    sum_v m_v^2 approximately 2 p^2

up to projective normalization constants.

## Final answers

(1) Existing Ostafe-Shparlinski/Chang/Gomez-Nicolas-Ostafe-Sadornil orbit theorems are essentially autonomous. They do not provide a direct theorem for the Apery nonautonomous Mobius cocycle.

(2) The Casoratian identity is exactly the missing non-degeneracy input: it gives separability/non-collapse analogous to autonomous non-isotriviality. It may make a new cocycle theorem possible, but it does not follow from current autonomous results.
