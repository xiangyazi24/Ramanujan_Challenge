#!/usr/bin/env python3
"""Archive exact Problem 3.2 ChatGPT prompts from the local bridge.

Historical Notion answer drops were saved as /tmp/gpt_Q<N>.md, while the
exact submitted questions remained in the bridge database.  The global
ledger maps Q numbers to run numbers, and runs.log maps run numbers to task
IDs.  This script joins those records, fetches /api/result/<task>, writes the
exact question text next to the already archived answers, and rebuilds a
compact index.

Only Q numbers with an existing answer in problems/3.2/chatgpt-answers are
processed.  This prevents unrelated projects from leaking into this archive.
"""

from __future__ import annotations

import argparse
import json
import re
import urllib.request
from pathlib import Path


LEDGER_ROW = re.compile(
    r"^\| Q(?P<q>\d+) \| (?P<time>[^|]*)\| (?P<channel>[^|]*)\| "
    r"(?P<topic>[^|]*)\| (?P<status>[^|]*)\| (?P<run>[^|]*)\|$"
)
RUN_ROW = re.compile(r"^RUN#(?P<run>\d+).* task=(?P<task>\S+)")


def parse_ledger(path: Path) -> dict[int, dict[str, str]]:
    rows: dict[int, dict[str, str]] = {}
    for line in path.read_text().splitlines():
        match = LEDGER_ROW.match(line)
        if not match:
            continue
        row = {key: value.strip() for key, value in match.groupdict().items()}
        rows[int(row.pop("q"))] = row
    return rows


def parse_runs(path: Path) -> dict[str, str]:
    rows: dict[str, str] = {}
    for line in path.read_text().splitlines():
        match = RUN_ROW.match(line)
        if match:
            rows[match.group("run")] = match.group("task")
    return rows


def fetch_task(base_url: str, task: str) -> dict[str, object]:
    url = f"{base_url.rstrip('/')}/api/result/{task}"
    with urllib.request.urlopen(url, timeout=15) as response:
        return json.load(response)


def first_heading(path: Path) -> str:
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("#"):
            return line.lstrip("#").strip()
    return ""


def write_question(
    path: Path,
    q: int,
    row: dict[str, str],
    task: str,
    payload: dict[str, object],
) -> None:
    question = str(payload.get("question") or "")
    # Delivery routing contains connector IDs and is not part of the research
    # question.  Keep the exact user/model prompt while excluding that transport
    # footer from the public archive.
    question = question.split("\n\n----\nDELIVERY INSTRUCTIONS", 1)[0].rstrip()
    metadata = [
        f"# Q{q} question",
        "",
        f"- Channel: `{row.get('channel', '')}`",
        f"- Ledger topic: `{row.get('topic', '')}`",
        f"- Bridge task: `{task}`",
        f"- Bridge status when archived: `{payload.get('status', '')}`",
        "",
        "## Exact submitted text",
        "",
        question,
        "",
    ]
    path.write_text("\n".join(metadata))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bridge",
        default="http://127.0.0.1:8801",
        help="Local ChatGPT bridge base URL",
    )
    parser.add_argument(
        "--ledger",
        type=Path,
        default=Path.home() / ".chatgpt-bridge" / "ASK_LEDGER.md",
    )
    parser.add_argument(
        "--runs",
        type=Path,
        default=Path.home() / ".chatgpt-bridge" / "runs.log",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    archive = root / "problems" / "3.2" / "chatgpt-answers"
    questions = archive / "questions"
    questions.mkdir(parents=True, exist_ok=True)

    ledger = parse_ledger(args.ledger)
    runs = parse_runs(args.runs)
    answers = sorted(
        (
            (int(match.group(1)), path)
            for path in archive.glob("Q*.md")
            if (match := re.fullmatch(r"Q(\d+)\.md", path.name))
        ),
        key=lambda pair: pair[0],
    )

    fetched = 0
    missing_task = 0
    fetch_failures = 0
    index_rows: list[str] = []
    for q, answer_path in answers:
        row = ledger.get(q, {})
        run_match = re.search(r"RUN#(\d+)", row.get("run", ""))
        task = runs.get(run_match.group(1), "") if run_match else ""
        question_path = questions / f"Q{q}.md"

        if task:
            try:
                payload = fetch_task(args.bridge, task)
                if payload.get("question"):
                    write_question(question_path, q, row, task, payload)
                    fetched += 1
                else:
                    fetch_failures += 1
            except Exception:
                fetch_failures += 1
        elif not task:
            missing_task += 1

        topic = row.get("topic", "").replace("|", "/")
        heading = first_heading(answer_path).replace("|", "/")
        question_link = f"[question](questions/Q{q}.md)" if question_path.exists() else "—"
        index_rows.append(
            f"| Q{q} | {topic} | {heading} | {question_link} | "
            f"[answer](Q{q}.md) |"
        )

    index = [
        "# Problem 3.2 ChatGPT question/answer index",
        "",
        "| Q | Ledger topic | Answer heading | Question | Answer |",
        "|---:|---|---|---|---|",
        *index_rows,
        "",
        (
            f"Archive run: fetched {fetched} new exact questions; "
            f"{missing_task} entries lacked a run-to-task mapping; "
            f"{fetch_failures} bridge fetches failed."
        ),
        "",
    ]
    (archive / "INDEX.md").write_text("\n".join(index))
    print(index[-2])


if __name__ == "__main__":
    main()
