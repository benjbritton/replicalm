# Replicalm publication workflow — status

**Written 2026-10-02. Paused awaiting human action.** Resume by reading this
file rather than restarting the workflow.

---

## Where it stopped, and why

`gh` is not installed on this machine (`gh: command not found`), so the GitHub
repository cannot be created, the remote cannot be authenticated, and the
release cannot be published from here. Every local step is complete. The exact
commands to run are below.

---

## Done

| step | state |
|---|---|
| Pre-publication audit | passed — see below |
| `LICENSE` (MIT, Benjamin Jay Britton, 2026) | already present, unchanged |
| `README.md` | already present and project-specific; Zenodo badge placeholder and Citation section added |
| `CITATION.cff` | created |
| `.gitignore` | already present and sufficient, unchanged |
| Branch renamed `master` → `main` | done |
| Commits | done, see below |
| Tag `v1.0.0` | created at `a41e3af` |
| Remote `origin` | **not added** — pending |
| Push | **not done** — pending |
| GitHub release + installer asset | **not done** — pending |

## Repository state

- Branch: `main`
- Latest commit: `75b9d8c` — "Replicalm v1.0.0 repository setup"
- Preceding commit: `ededa0e` — batch scheduler sized by measurement
- Tag: `v1.0.0` → `a41e3af` ("Remove a stray backspace byte from the README
  build command", 2026-09-24 12:37)
- Remotes: none
- History: 34 commits, all authored by Benjamin Jay Britton
  <benjaminbritton@yahoo.com>, no tool trailers
- Tracked content: 251 files, 54 MB

## The installer

`C:\Replicalm\packaging\dist\Replicalm-1.0.0-setup.exe`
1,211,198,530 bytes (1.21 GB), built 2026-09-24 12:52.

Not staged, not tracked, not modified, and excluded by `.gitignore`
(`packaging/dist/`). It is the authoritative v1.0.0 release installer and must
not be rebuilt or replaced. It is within GitHub's 2 GB per-asset limit.

## Pre-publication audit

- **Credentials**: none found in tracked files (scanned for API keys, secrets,
  passwords, private key headers, `ghp_`, `AKIA`).
- **Authorship**: single author across the whole history, no injected trailers.
- **.gitignore coverage**: holds back the packed conda environment
  (`packaging/stage/`), the 1.2 GB installer (`packaging/dist/`), the test and
  render trees, and all `.las`, `.laz`, `.tif`. Tracked total 54 MB.
- **Editor litter**: none tracked.
- **Machine paths**: `C:\Users\benja` appears as a fallback default in four
  benchmark scripts — `cleanup_sweep.py`, `clear_fulltile.py`,
  `fit_cleanup_rule.py`, `nugget_sweep.py`. Not secret; `benchmarks/paths.py`
  exists to override them by environment variable. Left as they are; worth
  tidying before a paper cites the harness.

---

## A decision for Ben before the push

The prompt asked for a fresh repository with one commit, "Replicalm v1.0.0
repository setup". `C:\Replicalm` was **already a git repository with 32
commits** of real development history, so that instruction was written without
that being known.

History was **preserved**, not squashed, because the alternative is
irreversible once pushed and because a methods-replication project's value
rests partly on its development record being inspectable. The v1.0.0 commit
message was used for the metadata commit at the tip.

Consequence to be aware of: pushing `main` publishes the full history,
including a week of post-1.0.0 work (the Lamanai ladder, the PMF correction,
the threshold experiment). The `v1.0.0` **tag** points at `a41e3af`, the exact
source the installer was built from, so the release is clean even though the
branch tip is ahead.

If a v1.0.0-only repository is wanted instead, say so **before the push** — it
means publishing an orphan branch from `a41e3af` and is easy now, impossible to
undo later.

---

## Next action required of Ben

1. Install the GitHub CLI, or create the repository through the web interface:
   - winget: `winget install --id GitHub.cli`
   - or create `benjbritton/replicalm`, public, **without** a README, license
     or .gitignore, at https://github.com/new
2. If using `gh`, authenticate: `gh auth login`
3. Confirm the history decision above.

## Commands to run once that is done

```bash
cd /c/Replicalm

# create the repository (skip if created through the web interface)
gh repo create benjbritton/replicalm --public --source=. --remote=origin --push

# or, if it was created on the web:
git remote add origin git@github.com:benjbritton/replicalm.git
git push -u origin main
git push origin v1.0.0
```

Then report: `git status`, current branch, latest commit hash, `git remote -v`,
and the repository URL.

## After the push

- **Step B** — Ben logs into Zenodo and toggles `benjbritton/replicalm` to ON.
  Nothing further happens here until that is confirmed, because Zenodo only
  mints a DOI for releases created *after* the switch is enabled.
- **Step C** — once confirmed:
  ```bash
  gh release create v1.0.0 \
    "packaging/dist/Replicalm-1.0.0-setup.exe" \
    --title "Replicalm 1.0.0" \
    --notes-file RELEASE_NOTES_v1.0.0.md
  ```
  The 1.21 GB upload will take a while.
- **Step D** — Ben copies the Zenodo DOI badge into `README.md` (replacing the
  two `PENDING` markers, one in the badge near the top and one in the Citation
  section), and uploads the installer to Scholar@UC as `brittobj`.
