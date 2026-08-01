#!/usr/bin/env bash
# truth.sh — "这个项目到底证完了什么" 的人话报告。
#
#   usage: scripts/truth.sh <LeanProjectDir> [模块前缀]
#   e.g.   scripts/truth.sh lean RamanujanChallenge.Problem32
#
# 你不需要知道 Lean 内部怎么回事。这个脚本问的是 Lean 内核本身，
# 不是问写代码的人、也不是读文档里的 ✅。它只回答一件事:
#
#   哪些结论是真的成立的, 哪些只是看起来成立。
#
# 底下用的是 #print axioms —— 它是传递的: 一个定理只要在任意深度上
# 依赖了一个未完成的证明, 这里就会显示出来。文档、注释、✅ 标记
# 都不参与判断。
set -uo pipefail

PROJ="${1:?usage: truth.sh <LeanProjectDir> [模块前缀]}"
PREFIX="${2:-}"
BANK="$(cd "$(dirname "$0")" && pwd)/bank.sh"   # 必须在 cd 之前定住绝对路径
cd "$PROJ" || exit 2

mods=$(find . -name "*.lean" -not -path "./.lake/*" -not -path "./scratch/*" 2>/dev/null \
       | sed 's|^\./||; s|/|.|g; s|\.lean$||' | sort)
[ -n "$PREFIX" ] && mods=$(echo "$mods" | grep "^$PREFIX")
[ -n "$mods" ] || { echo "找不到模块"; exit 2; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
real=0; fake=0; unknown=0

for m in $mods; do
  src="${m//.//}.lean"
  [ -f "$src" ] || continue
  bash "$BANK" . "$m" "$tmp/ledger.jsonl" > "$tmp/out_$m" 2>&1
  grep -qE '^\s+[✓✗?]' "$tmp/out_$m" || continue

  echo "━━━ ${m#*.} ━━━"
  while IFS= read -r line; do
    case "$line" in
      *"✓"*)
        n=$(echo "$line" | awk '{print $2}')
        st=$(python3 -c "
import json,sys
for l in open('$tmp/ledger.jsonl'):
    r=json.loads(l)
    if r['name']=='$n': print(r['statement'][:88]); break
" 2>/dev/null)
        printf "  \033[32m真的证完了\033[0m  %-32s %s\n" "$n" "$st"; real=$((real+1)) ;;
      *"sorryAx"*)
        n=$(echo "$line" | awk '{print $2}')
        # 它到底卡在谁身上: 在证明正文里找其它未完成的引理名
        dep=$(awk "/^(theorem|lemma|private lemma).* $n[ (:]/,/^\$/" "$src" 2>/dev/null \
              | grep -oE '[a-zA-Z_][a-zA-Z0-9_]*' | sort -u \
              | grep -Fxf <(grep -oE '^(private )?(theorem|lemma) [a-zA-Z_][a-zA-Z0-9_]*' "$src" | awk '{print $NF}') 2>/dev/null \
              | grep -v "^$n$" | head -3 | tr '\n' ' ')
        printf "  \033[31m没证完\033[0m      %-32s %s\n" "$n" "${dep:+← 卡在: $dep}"; fake=$((fake+1)) ;;
      *"??"*)
        n=$(echo "$line" | awk '{print $2}')
        printf "  \033[33m验不了\033[0m      %-32s %s\n" "$n" "(private 或特殊名字, 需要人看)"; unknown=$((unknown+1)) ;;
    esac
  done < <(grep -E '^\s+[✓✗?]' "$tmp/out_$m")
  echo
done

echo "════════════════════════════════════════"
echo "  真的证完了: $real 条"
echo "  没证完:     $fake 条   ← 文档里的 ✅ 不算数, 以这里为准"
echo "  验不了:     $unknown 条"
echo "════════════════════════════════════════"
echo
echo "注意: '真的证完了' 只保证逻辑上无漏洞, 不保证这条定理"
echo "      说的是你论文里想说的那件事。陈述本身对不对, machine 判不了。"
