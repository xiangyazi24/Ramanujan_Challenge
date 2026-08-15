ANSWER Q8535 3dd1edae

# CHATGPT_ATUP3_SIGN_CUBE_TABLE attempt 1

## Terminal report

Status: BLOCKED (named finite-verification blocker).

## Scope checked

The required deliverable is an exact eight-vertex reconstruction for three thin sign-simple rows:

- (256,25874;307,367,463)
- (256,51931;283,347,367)
- (512,98926;659,857,953)

The acceptance condition requires the actual Apéry zero roots at each prime, all CRT vertices, orientation checks, unchanged-coordinate divisors, exact quotients, complement identities, and a machine-checkable search for surviving arity-two relations.

## Findings

No complete executable reconstruction was verified from the available read-only inputs in this worker context. Therefore no claimed CRT table, quotient relation, or uniform certificate is asserted.

The surviving mathematically correct conclusion is limited:

1. Reflection symmetry alone is insufficient for an ATUP3 certificate.
2. A candidate relation must be checked vertex-by-vertex on all eight CRT representatives, not inferred from the existence of a reflection pair.
3. Failure of any proposed small-quotient/common-divisor/linear relation on one vertex is an exact falsification of that candidate, but does not rule out unrelated certificates.

## Residual obligation

Run the standard-library reconstruction script against the authoritative zero-set data and emit:

- the two roots for each prime;
- the eight CRT representatives in [0,p*q*ell);
- the below-X^2 predicate;
- differences from m;
- unchanged-coordinate prime products;
- exact quotients;
- reflection complement identities;
- machine-readable JSON output.

Only after that enumeration can a uniform arity-two ATUP3 certificate or a sharply falsified candidate relation be certified.

## Claim

No verified finite certificate was produced in this attempt. The blocker is missing exact reconstruction data, not a mathematical disproof of ATUP3.
