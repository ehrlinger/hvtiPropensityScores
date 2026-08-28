# Claude Code specifics

@AGENTS.md

[`AGENTS.md`](https://ehrlinger.github.io/hvtiRpropensity/AGENTS.md),
imported above, is the operational contract and applies in full. It is
written to be tool neutral so that Codex and other agents read the same
rules. Only the Claude Code affordances live here.

## Before you touch code

`AGENTS.md` says to orient before editing. In Claude Code the way to do
that is the codemap: it lives in the Obsidian vault under
`Claude/repomaps/` and is read via the `read-codemap` skill
(`/codemap hvtiRpropensity`). If the codemap looks stale, say so and
offer to refresh it (`/regenerate-codemap`) rather than working from a
guess.

If the vault is not available, say so rather than staying quiet about
it, then orient from the repo itself — `NAMESPACE`, `R/ps-data.R` for
the object contract every function returns, and the README — before
editing.

## Prose

`AGENTS.md` points at the composed house style. In Claude Code, apply
the `ehrlinger-writing` skill instead: it carries the same voice, reader
persona and project context, kept in sync from the vault sources
`.claude/house-style.md` is composed from. For documentation
*structure*, the `r-package-style` skill is the companion.
