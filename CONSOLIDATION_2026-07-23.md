# Repository consolidation — 2026-07-23

`/home/xhuan5/Ramanujan_Challenge` is the sole authoritative checkout for
this project.

The consolidation imported:

- source files for every challenge problem that existed only in the
  non-Git snapshot `/home/xhuan5/repos/Ramanujan_Challenge`;
- the Mini-synchronized Lean project already arriving on `origin/main`,
  excluding `.lake` and generated artifacts;
- archived ChatGPT answers;
- the full 2026-07-23 Problem 3.2 research campaign, including failed routes
  and executable audits.

The two old locations were deliberately not deleted.  They remain temporary
recovery sources until this commit is verified on GitHub.

## Lean provenance

While consolidation was in progress, seven Mini-generated Lean commits
landed on `origin/main`.  Those remote files were kept as authoritative in
every overlap; the stale non-Git snapshot was not allowed to overwrite them.
Lean is retained here for provenance and future formalization, but the
current priority for Problem 3.2 is the unconditional paper proof.

`Problem32/Main.lean` still contains three explicit `sorry`
declarations.  It must not be described as a formal proof of the open
pointwise conjecture.
