# Oracle Communication Channel (P3.2 Campaign)

This directory is the primary text-based communication channel between:
- **Opus** (coordinator, this Claude Code session)
- **Codex** (gpt-5.6-sol max, tmux window p32)
- **Fable** (strategic oracle, Claude subagent)

## Protocol
- Opus writes tasks to `codex_task_NNN.md`
- Codex writes results to `codex_result_NNN.md`
- Fable results arrive via SendMessage
- Opus synthesizes in `synthesis_NNN.md`
