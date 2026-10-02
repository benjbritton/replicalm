# Replicalm publication status

**v1.0.0 is published, archived and citable as of 2026-10-02.**

| | |
|---|---|
| Repository | https://github.com/benjbritton/replicalm |
| Release | https://github.com/benjbritton/replicalm/releases/tag/v1.0.0 |
| DOI | [10.5281/zenodo.23110738](https://doi.org/10.5281/zenodo.23110738) |
| Branch | `main`, tracking `origin/main` |
| Tag | `v1.0.0` → `e1f15b2` |
| Installer asset | `Replicalm-1.0.0-setup.exe`, 1,211,198,530 bytes, uploaded |
| SHA256 | `9590b7cfc064b99b6210c557c2e2bf6f2418340ca71624fea938889d5db325bb` |

## What was published

241 files, 39.1 MB: the software under `src/replicalm`, the benchmark harness
and its recorded results, `docs/open_observations.md`, `docs/processing_report.md`,
`docs/replicalm_overview.md`, two posts, the packaging recipe, the installer
screenshots, and the current overview and report.

Held back deliberately, on disk and gitignored: the superseded overview and
report drafts, `docs/paper/` (the ISPRS draft is unfinished and parts of it were
found wrong on 2026-10-02), `HANDOFF.md`, the two `.docx` post drafts whose
finished versions are published as markdown, the packed conda environment, the
installer itself, and all point clouds and rasters.

## The email rewrite, and why the hashes moved

GitHub refused the first push: the account blocks command-line pushes that
expose the author's email, and all commits carried a personal address. Rather
than disable that protection, every commit was rewritten to the account's
GitHub noreply address, `317455538+benjbritton@users.noreply.github.com`.

**Every commit hash changed.** `v1.0.0` points at `e1f15b2`, which was
`a41e3af` before the rewrite -- same tree, same message, same date, same author
name. That commit is still the last one before the installer was compiled at
12:52 on 2026-09-24, so the provenance claim is unchanged.

A trap worth recording: rewriting commits does **not** rewrite an annotated
tag's own tagger identity. The tag had to be deleted and recreated, and GitHub
rejected it separately until that was done.

The pre-rewrite history survives locally on the `backup-pre-email-rewrite`
branch and under `refs/original/`. Neither was pushed, and both still carry the
personal address, so neither should be.

## What Zenodo holds

Zenodo archives a zip of the repository source at the tag -- about 39 MB. It
does **not** hold the 1.21 GB installer, which exists only as a GitHub release
asset. If the installer should be in the citable record as well, it has to be
uploaded to the Zenodo record by hand; the default record limit is 50 GB, so
size is not the obstacle.

**Undecided.** Worth settling, because the point of archiving is that the
artefact outlives the host.

## Still outstanding

- **Scholar@UC deposit** of the installer, as `brittobj`. Not started.
- **Work not yet published.** The slope-aware cleanup built on 2026-10-02 is
  stashed, not committed: `git stash list` shows it as "today's slope-aware
  cleanup work, held during the email rewrite". It replaces the cleanup's
  percentile floor with a fitted plane, which stops the rule removing ground for
  being on a slope. Validated on synthetic ground and on l0s395; the fit radius
  is not settled, and a second tile was being processed to test whether
  1.00-1.20 m generalises.
- **The tonal gap.** Every configuration renders about 18 levels darker than the
  archive, and nothing tried so far touches it. The unexamined lead is what
  normalisation the published imagery actually used, in ArcGIS Pro or Surfer.
