#!/usr/bin/env bash
# bank.sh — record a Lean result in the ledger, and REFUSE to record anything
# this script did not just verify itself.
#
#   usage: scripts/bank.sh <LeanProjectDir> <Module.Name> [ledger.jsonl]
#   e.g.   scripts/bank.sh lean RamanujanChallenge.Problem32.Main
#
# Design (from Danus's fact graph, adapted so the verifier is the Lean kernel
# instead of an LLM):
#
#   1. NOBODY CAN CLAIM. There is no "mark as done" input. The only way an entry
#      exists is for this script to build the module and read `#print axioms`
#      itself. A subagent reporting success changes nothing.
#   2. `#print axioms` is the gate, not grep. It is transitive: a `sorryAx`
#      three layers down still surfaces. A source grep misses subdirectories,
#      re-exports, and dependencies.
#   3. CONTENT-ADDRESSED. id = sha256(module + name + normalized statement)[:16].
#      Edit the statement and the id changes, so every dependent entry that
#      cites the old id is mechanically orphaned — you cannot silently move a
#      lemma out from under its consumers.
#   4. THE STATEMENT IS STORED IN FULL. `#print axioms` is SILENT on carried
#      hypotheses: a theorem can be axiom-clean and still be conditional, or
#      vacuous. The ledger records the whole type so an audit can see every
#      `(h : ...)` that was carried. Machine-clean is necessary, never sufficient.
set -uo pipefail

PROJ="${1:?usage: bank.sh <LeanProjectDir> <Module.Name> [ledger.jsonl]}"
MOD="${2:?missing module, e.g. RamanujanChallenge.Problem32.Main}"
LEDGER="${3:-$(cd "$(dirname "$PROJ")" && pwd)/ledger.jsonl}"

cd "$PROJ" || { echo "no such project dir: $PROJ" >&2; exit 2; }
SRC="${MOD//.//}.lean"
[ -f "$SRC" ] || { echo "no source for module $MOD (looked for $SRC)" >&2; exit 2; }

ts=$(date -u +%Y-%m-%dT%H:%M:%SZ)
probe=$(mktemp -t bank).lean
trap 'rm -f "$probe" "${probe%.lean}"' EXIT

# --- discover the declarations this module states itself -----------------
names=$(grep -oE '^(private[[:space:]]+)?(protected[[:space:]]+)?(theorem|lemma)[[:space:]]+[A-Za-z_][A-Za-z0-9_'"'"'.]*' "$SRC" \
        | sed -E 's/.*(theorem|lemma)[[:space:]]+//' | sort -u)
[ -n "$names" ] || { echo "no theorems found in $SRC" >&2; exit 2; }

# --- step 1: build. a module that does not build banks nothing. ----------
echo "[1/3] lake build $MOD"
if ! lake build "$MOD" >/tmp/bank_build.log 2>&1; then
  echo "REFUSED: $MOD does not build. nothing written." >&2
  tail -20 /tmp/bank_build.log >&2
  exit 1
fi

# --- step 2: ask the kernel, do not ask the author -----------------------
{
  echo "import $MOD"
  echo "open Real in"
  while IFS= read -r n; do
    echo "#print axioms $n"
    echo "#check @$n"
  done <<< "$names"
} > "$probe"

echo "[2/3] #print axioms via lake env lean"
out=$(lake env lean "$probe" 2>&1)

# --- step 3: one ledger line per declaration -----------------------------
echo "[3/3] writing ledger -> $LEDGER"
banked=0; refused=0
while IFS= read -r n; do
  ax=$(echo "$out" | grep -F "'$n' depends on axioms:" | head -1 | sed "s/.*depends on axioms: //")
  [ -z "$ax" ] && ax=$(echo "$out" | grep -F "'$n' does not depend on any axioms" | head -1)
  stmt=$(echo "$out" | grep -A2 -F "@$n :" | head -3 | tr '\n' ' ' | sed 's/  */ /g')
  [ -z "$stmt" ] && stmt="(type not captured)"

  if [ -z "$ax" ]; then
    echo "  ?? $n — no axiom line captured, REFUSED"; refused=$((refused+1)); continue
  fi
  if echo "$ax" | grep -q "sorryAx"; then
    echo "  ✗ $n — sorryAx, REFUSED"; refused=$((refused+1)); continue
  fi

  norm=$(echo "$stmt" | tr -s ' ')
  id=$(printf '%s|%s|%s' "$MOD" "$n" "$norm" | shasum -a 256 | cut -c1-16)
  # carried hypotheses: a crude but useful count of binders in the type
  hyps=$(echo "$norm" | grep -o '∀\|→' | wc -l | tr -d ' ')

  python3 - "$id" "$MOD" "$n" "$ax" "$norm" "$ts" "$hyps" "$LEDGER" <<'PY'
import json,sys
id_,mod,name,ax,stmt,ts,hyps,ledger = sys.argv[1:9]
rec = {"fact_id":id_,"module":mod,"name":name,"axioms":ax.strip(),
       "statement":stmt.strip(),"binder_count":int(hyps),"banked_utc":ts,
       "gate":"lake build + #print axioms"}
with open(ledger,"a") as f: f.write(json.dumps(rec,ensure_ascii=False)+"\n")
PY
  echo "  ✓ $n  [$id]  axioms: $(echo "$ax" | cut -c1-60)"
  banked=$((banked+1))
done <<< "$names"

echo
echo "banked=$banked refused=$refused  ledger=$LEDGER"
echo "NOTE: axiom-clean != correct. The ledger stores the full statement so a"
echo "human can check what hypotheses were carried; #print axioms cannot."
[ "$banked" -gt 0 ] || exit 1
