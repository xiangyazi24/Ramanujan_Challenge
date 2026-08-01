# Ramanujan Challenge — Problem 2.2

**Euler's constant gamma as an Apery limit**

Submitter: Xiang Huang (University of Illinois Springfield)
Contact: xhuan5@uis.edu

## What is in this package

- solution.pdf / solution.tex — human-readable proof
- lean/ — Lean 4 formalization (0 sorry)
- README.md — this file

## Summary

The three-term recurrence from the challenge is identified as a first-order
Ore transform of Rivoal's factorial-scaled construction. Positivity of the
denominator is proved directly from the recurrence via a positive cone.
The limit P_n/Q_n -> gamma is proved unconditionally via harmonic concentration.

## Lean verification

4 files, 0 sorry. Main theorem: problem22_solved : Problem22Claim.
Axioms: {propext, Classical.choice, Quot.sound}.
