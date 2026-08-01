#!/usr/bin/env bash
# deps.sh — 从 Lean 内核里扒出真实依赖边, 补进账本。
#
#   usage: scripts/deps.sh <LeanProjectDir> <Module.Name> [ledger.jsonl]
#
# 和 bank.sh 里那个 "← 卡在: X" 不同 —— 那个是文本 grep 猜的。
# 这里问的是 Lean 环境本身: 一条定理的证明项里真正引用了哪些常量
# (thmInfo.value.getUsedConstants), 再筛掉 Mathlib、只留本项目的。
#
# 有了真依赖边才能做两件 grep 做不到的事:
#   1. 改了上游 -> 哈希变 -> 机械列出所有必须重验的下游
#   2. 证明链里的循环依赖, 机械查得出来
set -uo pipefail

PROJ="${1:?usage: deps.sh <LeanProjectDir> <Module.Name> [ledger.jsonl]}"
MOD="${2:?missing module}"
LEDGER="${3:-/tmp/deps_ledger.jsonl}"
ROOTNS="${MOD%%.*}"          # RamanujanChallenge.Problem32.Main -> RamanujanChallenge

cd "$PROJ" || exit 2
SRC="${MOD//.//}.lean"
[ -f "$SRC" ] || { echo "no source for $MOD" >&2; exit 2; }

names=$(grep -oE '^(private[[:space:]]+)?(protected[[:space:]]+)?(theorem|lemma)[[:space:]]+[A-Za-z_][A-Za-z0-9_'"'"'.]*' "$SRC" \
        | sed -E 's/.*(theorem|lemma)[[:space:]]+//' | sort -u)
nss=$(grep -oE '^namespace[[:space:]]+[A-Za-z_][A-Za-z0-9_.]*' "$SRC" | awk '{print $2}')
[ -n "$names" ] || { echo "no theorems in $SRC" >&2; exit 2; }

probe=$(mktemp -t deps).lean
trap 'rm -f "$probe" "${probe%.lean}"' EXIT

{
  echo "import Lean"
  echo "import $MOD"
  echo "open Lean"
  echo ""
  echo "run_cmd do"
  echo "  let env ← getEnv"
  # 名字用 open 过的短名解析不到 run_cmd 里的 Name 字面量, 所以两种都试
  printf '  let targets : List Name := ['
  sep=""
  while IFS= read -r n; do
    [ -z "$n" ] && continue
    printf '%s`%s' "$sep" "$n"; sep=", "
    for ns in $nss; do printf ', `%s.%s' "$ns" "$n"; done
  done <<< "$names"
  echo "]"
  cat <<'LEAN'
  for n in targets do
    match env.find? n with
    | some (.thmInfo tv) =>
        let cs := tv.value.getUsedConstants
        let ours := cs.filter fun c =>
          match env.getModuleFor? c with
LEAN
  echo "          | some m => (\`$ROOTNS).isPrefixOf m"
  cat <<'LEAN'
          | none => false
        if ours.size > 0 then
          let joined := String.intercalate "," (ours.toList.map toString)
          IO.println s!"DEP {n} -> {joined}"
    | _ => pure ()
LEAN
} > "$probe"

echo "[deps] 抽取 $MOD 的依赖边..."
out=$(lake env lean "$probe" 2>&1)
if echo "$out" | grep -q "error:"; then
  echo "探针编译失败:" >&2; echo "$out" | grep "error:" | head -5 >&2; exit 1
fi

echo "$out" | grep '^DEP ' | while IFS= read -r line; do
  n=$(echo "$line" | awk '{print $2}')
  d=$(echo "$line" | sed 's/.*-> //')
  short=${n##*.}
  cnt=$(echo "$d" | tr ',' '\n' | grep -c . )
  printf "  %-38s 依赖本项目 %2s 条\n" "$short" "$cnt"
  python3 - "$short" "$MOD" "$d" "$LEDGER" <<'PY'
import json,sys,os
name,mod,deps,ledger = sys.argv[1:5]
import re
# 滤掉编译器生成的辅助常量: _proof_1_3 / match_1_1 / foo.eq_def / _sunfold ...
NOISE = re.compile(r'^_|^match_\d|_proof_\d|^.*\.(eq_def|eq_\d+|_sunfold|_unfold|induct|below|brecOn)$')
seen = sorted({d.split(".")[-1] for d in deps.split(",") if d})
clean = [d for d in seen if not NOISE.match(d)]
rec = {"name":name,"module":mod,"depends_on":clean}
with open(ledger,"a") as f: f.write(json.dumps(rec,ensure_ascii=False)+"\n")
PY
done

echo "[deps] 写入 $LEDGER"
