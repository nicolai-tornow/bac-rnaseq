# Codex CLI tool mapping

- "run a shell command" → shell tool
- plugin asset paths → resolve relative to the skill's own directory (Codex has
  no `${CLAUDE_PLUGIN_ROOT}`; use the skill dir as the anchor)
- "ask the user to confirm" → ask inline in the conversation
