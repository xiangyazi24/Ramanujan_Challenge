#!/bin/bash
# Build the submission archive for the Ramanujan Challenge.
#
# Usage:  bash make_zip.sh 2.8            -> SUBMIT/dist/ramanujan-2.8-huang.zip
#         bash make_zip.sh 2.8 3.1        -> a single archive containing both
#
# The challenge accepts only .zip, max 50 MB, and each solution must contain a
# solution.tex or solution.pdf.  This script checks both conditions.

set -e
cd "$(dirname "$0")"
OUT=dist
mkdir -p "$OUT"

if [ $# -eq 0 ]; then
  echo "usage: bash make_zip.sh <problem> [<problem> ...]"; exit 1
fi

STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT

for p in "$@"; do
  if [ ! -d "$p" ]; then echo "no such problem dir: $p"; exit 1; fi
  if [ ! -f "$p/solution.pdf" ] && [ ! -f "$p/solution.tex" ]; then
    echo "REFUSING: $p has neither solution.pdf nor solution.tex (required artifact)"
    exit 1
  fi
  mkdir -p "$STAGE/$p"
  # copy everything except build residue
  rsync -a --exclude '*.aux' --exclude '*.log' --exclude '*.out' \
        --exclude '.lake' --exclude '__pycache__' --exclude '*.olean' \
        "$p/" "$STAGE/$p/"
done

if [ $# -eq 1 ]; then
  NAME="ramanujan-$1-huang.zip"
else
  NAME="ramanujan-huang.zip"
fi

( cd "$STAGE" && zip -qr "$OLDPWD/$OUT/$NAME" . )

SZ=$(du -m "$OUT/$NAME" | cut -f1)
echo "wrote $OUT/$NAME  (${SZ} MB)"
if [ "$SZ" -ge 50 ]; then
  echo "*** OVER THE 50 MB LIMIT ***"; exit 1
fi

echo
echo "--- archive contents (top level) ---"
unzip -l "$OUT/$NAME" | head -25
echo
echo "--- required artifact check ---"
for p in "$@"; do
  if unzip -l "$OUT/$NAME" | grep -q "$p/solution.pdf"; then
    echo "  $p/solution.pdf  present"
  else
    echo "  $p/solution.pdf  MISSING"
  fi
done
