---
name: report-feedback
description: Use to report a bug or suggest an improvement for the bac-rnaseq plugin. Files a GitHub issue on the plugin's own repository using the user's GitHub account (no separate feedback repo).
---

# Report Feedback

Files feedback as a GitHub issue on **this plugin's own repo**
(`nicolai-tornow/bac-rnaseq`) — there is no separate feedback repository. Anyone
with repo access (i.e. org members) can file. Tool names per
`skills/_shared/references/<harness>-tools.md`.

1. Ask the user what kind of feedback (**bug** or **idea**), a one-line summary,
   and any details or steps to reproduce.
2. Offer to attach environment info from `${CLAUDE_PLUGIN_ROOT}/bin/bac-rnaseq doctor`
   (cores, tool presence) — helpful for bugs. Never include secrets, tokens, or
   private data in the issue.
3. File it with the user's own GitHub auth:
   ```bash
   gh issue create --repo nicolai-tornow/bac-rnaseq \
     --title "<bug|idea>: <summary>" \
     --body "<details, repro, optional environment info>"
   ```
   (Requires `gh auth login`; org members with repo access can open issues.)
4. Report the issue URL back to the user.
