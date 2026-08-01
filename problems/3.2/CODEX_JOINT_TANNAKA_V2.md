# Joint Tannakian moments v2: STALL

## Verdict

A mandatory exact gate failed, so no FFT or moments were computed.

- Failed prime: `37`
- Gate: `g2`
- Witness: `t=3, T+ mod p=24, A=13; inert, d=19, a2=50, #E=1320, epsilon=-1, independent binomial A=13, independent scalar a2=50, full-enumeration a2=50`
- Earlier primes whose gates passed: `29, 31`

The run followed the mandatory abort rule in `CODEX_SPEC_joint_tannaka_v2.md`.  Consequently the success acceptance criterion was not reached; returning a nonzero status is intentional.

## Independent checks on the witness

The Apéry polynomial was computed by the requested recurrence.  When the residue failed, all coefficients were recomputed independently from the binomial formula.  The inert elliptic trace was first obtained from the exact norm-character point count in `F_{p^2}` and then recomputed by a separate scalar implementation using exponentiation inside `F_p[z]/(z^2-d)`.  At this small witness the script additionally exhausted all `(x,y)` pairs.  All three values are printed in the witness.

At the first failing point, `p=37, t=3`, one has `d=19`, `epsilon=(2/37)=-1`, and the exact count is `#E(F_{37^2})=1320`. Thus `a2=37^2+1-1320=50`.  Q6457 prescribes `T+=epsilon*a2-p=-87`, whose residue is 24, whereas both Apéry computations give `A_37(3)=13`.

As a diagnostic only, `(-3/37)=+1`, so replacing the prescribed factor `(2/p)` by `(-3/p)` repairs this witness and is consistent with the repository's separately verified constant `(-3)` twist. That replacement was not made: it is outside the supplied recipe and requires a corrected descent-normalization derivation.

## Checksums reached before the stall

At p=29, the recurrence/direct/CRON Apéry arrays agreed; every split trace agreed among direct point counts, centered Franel Hasse--Witt residues, and raw CRON counts; conjugating every inert square root preserved the degree-two count; and the named checksum was `t=2`, sources `(8,10)`, traces `(-6,6)`, `f=7`.

Gates g1--g5 passed completely at p=29 and p=31.  The first failure was then g2 at p=37, before the remaining primes and before all floating-point work.

## Limitation

This stall identifies a normalization contradiction, not a Tannakian-group verdict.  With no verified integral trace table, finite-p moments would be arbitrary and are deliberately omitted.
