# AGENTS.md

Instructions for AI coding agents working in this repository.

## Project Context

VDS-Jointcall is a set of scripts used to perform joint variant calling via Hail VDS format.
We use it to perform joint calling:
* On a huge datasets that can't be handled by classin joint callers (like GATK)
* On datasets called by Deepvariant, that do not have genotype annotations required gor GATK joint calling

The current development direction is iterative improvement of coding style,
maintainability, and reliability.
Changes should be small, scoped, and easy to review.

Follow the plan-implement-validate approach:
* Plan the changes and ask for user confirmation
* When the user confirms the plan, proceed with the implementation
* Finally, run validation to check that the changes are safe and do not break anything.

**Always ask to confirm the proposed changes plan before proceeding with the implementation.**

## Non-Negotiable Scope Rules

- Only modify files that are clearly inside the approved task scope.
- Do not perform opportunistic cleanup, broad formatting, import sorting, or
  refactoring in files that are not part of the requested change.
- If a required fix appears to need changes outside the approved scope,
  stop and ask the user for confirmation.
- If task instructions conflict with the actual codebase, stop and ask for
  further instructions.
- If the README, docs, tests, type hints, comments, or function names disagree
  with executable behavior, inspect the code and ask user before making assumptions that affect behavior.

## Environment Assumptions

- Assume the required configuration, tools, credentials, permissions, and data access have been set up correctly for the assigned task.
- NEVER work around missing tooling, unavailable data, permission errors, cloud access issues,
  or broken local configuration by changing project code or widening scope.
- If tooling, permissions, credentials, environment variables, data access,
  or local configuration appear to be missing or incorrect, stop and ask the user to fix the environment before continuing.

## Working Style

- Keep changes minimal and behavior-preserving unless the task explicitly asks
  for behavior changes.
- Add comments only where they clarify algorithm, non-obvious Hail/Spark behavior,
  historical constraints, or complex transformations.

## Python and Hail Conventions

- The project targets Python 3.12.

## Documentation

- Keep documentation changes scoped to the behavior or workflow being changed.
- Do not rewrite broad README or docs sections unless requested.
- If code behavior and documentation disagree, stop and ask before changing one
  to match the other, unless the task explicitly identifies which source should
  be corrected.

## Final Response Expectations

When finishing a task, summarize:

- What changed.
- Which files were touched.
- Which tests or checks were run.
- Any tests or checks that were not run and why.
- Any out-of-scope issue discovered that still needs user approval.
