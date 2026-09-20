# Archive

Superseded renders and experiments. Nothing here is tracked by git — these are
derived products measured in gigabytes, and every one of them can be
regenerated from the tracked code and the recorded parameters.

They are kept on disk because they are the evidence behind the claims in
`docs/posts/` and `README.md`, and because the comparison figures were cut from
them.

| directory | what it is |
|---|---|
| `render_r15` | 1.5 m search radius. Superseded: the apparent improvement was a scoring artifact, since cells that could not be solved were dropped rather than counted. |
| `../render` | the flat-calibrated configuration (CSF rigidness 2, one pass, filters off, 14.4 m radius). The middle panel in every comparison figure. |
| `../render_smrf` | SMRF single pass at a 5 m radius. The intermediate step. |
| `../render_ncalm` | **the accepted product** — locked baseline. Not superseded; left in place. |

Three of those directories are still at the repository root rather than here
because they were held open by another process at the time of the move. They
are gitignored either way.

If any of this needs to travel, git-lfs is the option; at roughly 3.4 GB it is
worth deciding deliberately rather than by default.
