#!/usr/bin/env bash
# impact.sh — "我动了这条, 谁必须重验?"
#
#   usage: scripts/impact.sh <定理名> [deps_ledger.jsonl]
#
# 这是依赖图存在的理由。改一条上游引理的陈述, 它的哈希就变了,
# 所有(传递地)依赖它的下游结论都不再算已验证 —— 不是"建议看一下",
# 是逻辑上真的失效了。这个脚本把那份名单机械地列出来。
#
# 反过来也用得上: 一条结论卡住时, 看它压在谁身上 (--why)。
set -uo pipefail

TARGET="${1:?usage: impact.sh <定理名> [ledger] | impact.sh --why <定理名> [ledger]}"
if [ "$TARGET" = "--why" ]; then
  MODE=why; TARGET="${2:?missing name}"; LEDGER="${3:-/tmp/deps_ledger.jsonl}"
else
  MODE=impact; LEDGER="${2:-/tmp/deps_ledger.jsonl}"
fi
[ -f "$LEDGER" ] || { echo "找不到依赖账本 $LEDGER, 先跑 deps.sh" >&2; exit 2; }

python3 - "$MODE" "$TARGET" "$LEDGER" <<'PY'
import json,sys
from collections import defaultdict, deque

mode, target, ledger = sys.argv[1:4]
fwd = {}                      # name -> 它依赖谁
rev = defaultdict(set)        # name -> 谁依赖它
for line in open(ledger):
    r = json.loads(line)
    fwd[r["name"]] = r["depends_on"]
    for d in r["depends_on"]:
        rev[d].add(r["name"])

def walk(start, graph):
    seen, order, q = {start}, [], deque([(start,0)])
    while q:
        n, d = q.popleft()
        for m in sorted(graph.get(n, [])):
            if m not in seen:
                seen.add(m); order.append((m, d+1)); q.append((m, d+1))
    return order

if mode == "impact":
    hits = walk(target, rev)
    if not hits:
        print(f"没有任何已记录的结论依赖 {target} —— 改它不影响别人。")
        print("(注意: 只覆盖跑过 deps.sh 的模块; 没跑过的不在图里)")
    else:
        print(f"改动 {target} 之后, 以下 {len(hits)} 条必须重验:\n")
        for n, d in hits:
            print(f"  {'  '*(d-1)}└─ {n}   (第 {d} 层)")
        print(f"\n这不是建议 —— 它们的正确性是建立在 {target} 的当前陈述上的。")
else:
    deps = walk(target, fwd)
    if not deps:
        print(f"{target} 不依赖本项目任何其它结论 (只用了 Mathlib)。")
    else:
        print(f"{target} 压在以下 {len(deps)} 条上:\n")
        for n, d in deps:
            leaf = "" if n in fwd else "   ← 不在账本里, 未经验证"
            print(f"  {'  '*(d-1)}└─ {n}{leaf}")
PY
