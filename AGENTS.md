# AGENTS.md — vaxrank

Guide for coding agents working in this repo. Read this before touching code.

---

## Golden Rules

1. **Never commit to `main`.** Always `git checkout -b <feature-branch>` before editing. Land via PR.
2. **Every PR bumps the version.** Even doc-only PRs — at minimum a patch bump in `vaxrank/version.py`. `deploy.sh <version>` does the bump + commit + push for you.
3. **Always deploy when you merge a PR.** Merging without deploying is never done — every merge is followed by `./deploy.sh` from a clean `main`. No "wait for the next PR," no "deploy later." The exact sequence after `gh pr merge`:
   ```
   git checkout main && git pull
   ./deploy.sh                  # uses current version in vaxrank/version.py
   ```
   `./deploy.sh` runs lint + tests + build + `twine upload` + tag + push. If lint or tests fail on a clean `main`, fix the root cause; never bypass. Verify the new version is live on PyPI (`pip index versions vaxrank`) before declaring the task complete.
4. **File problems as issues, don't silently work around them.** If you hit a bug — here or in a dependency (openvax/pirl-unc: varcode, isovar, mhctools, topiary, mhcgnomes, pepdata, etc.) — open a GitHub issue on the correct repo. Link it from the PR.
5. **After a PR ships, look for the next block of work.** Read open issues across the relevant openvax repos, group by dependency + urgency. Prefer *foundational* changes that unblock multiple downstream improvements. If nothing foundational is pending, pick the smallest independent improvements and do them in a row.

---

## Before Completing Any Task

Before telling the user a change is "complete":

1. **`./lint.sh`** — must pass (runs `ruff check vaxrank tests`)
2. **`./test.sh`** — must pass (pytest + coverage)
3. For a PR: **CI must be green on GitHub**, then merge, then **`./deploy.sh`**.

`deploy.sh` itself gates on lint + test, refuses to run off `main`/`master`, and refuses a dirty tree — don't work around these. If deploy fails, fix the root cause; fallback recipe only when truly needed:

```
rm -rf dist build *.egg-info
python3 -m build --no-isolation
python3 -m twine check dist/*
python3 -m twine upload dist/*
git tag "v${VERSION}" && git push --tags
```

## Scripts

- `./develop.sh` — editable install (uses `uv` if available, else pip)
- `./lint.sh` — `ruff check` on `vaxrank/` and `tests/`
- `./test.sh` — pytest with coverage (handles macOS arm64 WeasyPrint DYLD)
- `./deploy.sh [version]` — lint → test → optional version bump → build → twine upload → tag → push
- `./run-vaxrank-b16-test-data.sh` — end-to-end smoke on B16 test data

**macOS arm64 note:** `test.sh`/`deploy.sh` re-export `DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib` so WeasyPrint can find Pango/Cairo. Invoke the scripts, not `pytest` directly.

## Project Layout

- `vaxrank/` — library + CLI (`vaxrank.cli.entry_point:main`)
- `vaxrank/templates/`, `vaxrank/data/` — Jinja report templates, static data
- `tests/` — pytest suite
- `setup.py` + `requirements.txt` — packaging (no `pyproject.toml` yet)
- `vaxrank/version.py` — single source of truth for version

## Code Style

- Python 3.9+, Linux + macOS (no Windows guarantee)
- Lint: `ruff check` (config minimal; no `ruff format` step in this repo)
- Docstrings: numpy style
- Bugfixes include a regression test where feasible

---

## Workflow Orchestration

### 1. Upfront Planning
- For any non-trivial task (3+ steps or architectural): write a short spec first. If something goes sideways, STOP and re-plan — don't keep pushing.

### 2. Verification Before Done
- Never claim complete without proof: tests green, CI green, PyPI version live.
- For CLI/report changes, run `./run-vaxrank-b16-test-data.sh` or equivalent smoke.

### 3. Autonomous Bug Fixing
- Given a bug report: just fix it. Point at logs/errors/failing tests and resolve them without hand-holding.

### 4. Demand Elegance (Balanced)
- For non-trivial changes pause and ask "is there a more elegant way?" — skip for trivial fixes.
- Treat workarounds as bugs, not new abstractions. Rip out legacy paths decisively rather than accumulating special cases.

### 5. Issue Triage After Each Ship
- Close superseded/outdated issues as you notice them.
- New problems encountered mid-task → file as issues (on the right repo, even if it's not vaxrank), don't bury.

---

## Domain Rules (vaxrank specifics)

- **Allele parsing**: use `mhcgnomes`. Never `startswith("HLA-")` or other string hacks — alleles aren't always human.
- **Epitope prediction**: go through `mhctools` / `topiary`, not MHCflurry directly.
- **Filtering + ranking**: use the topiary DSL (integrated in 2.1.0). Don't add ad-hoc filters inside vaxrank that duplicate topiary predicates.
- **Variants → peptides**: `varcode` for annotation (>=4.0.0 — uses the `fast` annotator as the registry default), `isovar` for RNA-backed variant peptides. There's a DNA-only fallback (2.3.0) and a prediction cache — use them, don't bypass.
- **Provenance**: carry gene name, gene ID, transcript, source species wherever peptides flow. Flanking context is not unique — a peptide can map to multiple source genes.
- **Reports**: Jinja templates → WeasyPrint PDF / ASCII / Excel. Keep report logic in templates + small Python glue, not sprawling formatters.

## Core Principles

- **Simplicity first.** Minimal diffs, minimal abstractions.
- **No laziness.** Find root causes; no temporary fixes, no empty-category fudges.
- **Minimal blast radius.** Touch only what the task requires.

## Scientific Domain Knowledge

- If a change touches immunology/genomics semantics, check primary sources (papers, UniProt, GenBank) before edits.
- If the code expresses a scientific model at odds with your understanding, flag it — don't silently "fix" it into something wrong.
