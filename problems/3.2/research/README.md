# Problem 3.2 research archive

This directory is the authoritative archive for the continuing attack on
Ramanujan Challenge Problem 3.2.

## Layout

- `working_notes/`: chronological doctrine, audits, failed routes, and the
  current frontier.  Start with
  `working_notes/Q32_CODEX_RESUME_2026-07-23.md`.
- `proofs/`: proof drafts for exact intermediate theorems.  A file in this
  directory is not automatically a completed proof of Problem 3.2; consult
  the audit classification in the current frontier note.
- `scripts/`: exact arithmetic experiments and reproducible counterexamples.
- `../chatgpt-answers/`: raw ChatGPT answer drops and extracted code
  attachments.

## Provenance and consolidation

On 2026-07-23 the active research files were found in two non-authoritative
locations:

1. a local scratch `working_notes/` directory;
2. a local `Q-series-and-Chan-s-work/` gitlink, stale and without an active
   repository checkout.

They were copied here without deleting the source copies.  The canonical
GitHub repository is now `xiangyazi24/Ramanujan_Challenge`; new Problem 3.2
work must be run from that checkout and written here.

## Proof status

The fully unconditional pointwise target remains open:

```text
gcd(d_n a_n, d_n b_n) = exp(o(n)).
```

The current exact bottleneck is the moving-zero sum

```text
L(n) = sum_{sqrt(n)<p<=n, p | b_(n mod p)} log p = o(n).
```

Do not promote numerical evidence, a generic-density statement, or an
unverified ChatGPT argument to a theorem.  The frontier log records explicit
counterexamples and the conditions under which failed routes may be reopened.
