# Copilot Instructions (Repo-wide)

This repo contains multiple independent project folders. Treat each issue as scoped to a single project folder named in the issue.

## Scope safety
- Never modify other top-level project folders.
- You may only edit:
  - the target <PROJECT_FOLDER> from the issue
  - history/<PROJECT_FOLDER>/**
  - report.md, naming.txt, suggestion.txt, suggestions_done.txt

## Output requirements
- Portfolio readiness: fix absolute paths, remove assignment traces, update README, make runnable only if broken.
- Git historian outputs go under history/<PROJECT_FOLDER>/:
  - github_steps.md
  - steps/step_01..step_N (full snapshots)

## Logging
- All proposed changes go to suggestion.txt (with locator + snippet).
- All applied changes go to suggestions_done.txt (before/after + locator).
- report.md must contain a section per project.

## Do not stop early
Maintain a checklist and do not finish until all acceptance criteria are satisfied.
