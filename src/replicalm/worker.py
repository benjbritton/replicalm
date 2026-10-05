"""Run a job in a child process and report progress as JSON lines.

WHY A SUBPROCESS AND NOT A THREAD
---------------------------------
The window used to process on a worker thread. Long numpy and PDAL operations
hold the interpreter lock for stretches of several seconds, so the main thread
could not service its own timer; repaints and clicks queued up and Windows
appended "(Not Responding)" to the title. That was not a false alarm -- the
window really had stopped pumping messages. Suppressing the symptom is possible
(`user32.DisableProcessWindowsGhosting`) and would only hide a frozen window.

In a child process nothing competes for the parent's interpreter lock, so the
window stays live for the whole run.

It also makes stopping possible. A Python thread cannot be killed: there is no
safe way to interrupt one part-way through a kriging. A child process can be
terminated, which is what lets the window offer Cancel and return to the state
it was in before Run, with every selection still filled in.

THE PROTOCOL
------------
One JSON object per line on stdout, flushed as it is written, so the parent can
read it with a blocking `readline` on a thread that does nothing else:

    {"t": "start",    "clouds": N}
    {"t": "file",     "index": i, "of": N, "path": "..."}
    {"t": "progress", "stage": "...", "message": "...", "fraction": 0.0-1.0}
    {"t": "result",   ...everything pipeline.process returned...}
    {"t": "error",    "path": "...", "message": "..."}
    {"t": "done",     "ok": N, "failed": N}

Anything the child writes that is not valid JSON is a bug or a library writing
to stdout behind our back; the parent shows such lines verbatim rather than
discarding them, because a message nobody sees is how the image step stayed
broken for a whole release.

The job arrives as a file rather than on the command line: it carries paths that
may contain quotes and spaces, and a Windows command line has a length limit a
few hundred tiles would exceed.
"""
import json
import os
import sys


def emit(obj):
    """One JSON object, one line, flushed. The parent blocks on these."""
    sys.stdout.write(json.dumps(obj) + "\n")
    sys.stdout.flush()


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        emit({"t": "error", "path": "", "message":
              "usage: python -m replicalm.worker <job.json>"})
        return 2
    with open(argv[0], encoding="utf-8") as fh:
        job = json.load(fh)

    from .config import ReplicalmConfig
    from . import pipeline

    cfg = ReplicalmConfig.from_json(job["config"])
    clouds = job["clouds"]
    out_dir = job["out"]
    total = len(clouds)
    emit({"t": "start", "clouds": total})

    ok = failed = 0
    for i, cloud in enumerate(clouds, start=1):
        emit({"t": "file", "index": i, "of": total, "path": cloud})
        # Each file owns its slice of the bar, so progress runs once across the
        # whole job instead of resetting per tile.
        base = (i - 1) / total
        span = 1.0 / total

        def progress(stage, message, fraction=None, _b=base, _s=span):
            emit({"t": "progress", "stage": stage, "message": message,
                  "fraction": None if fraction is None else _b + _s * fraction})

        try:
            result = pipeline.process(
                cloud, out_dir, cfg=cfg, cell_m=job.get("cell_m"),
                make_g1=job.get("make_g1", False),
                keep_ground=job.get("keep_ground", False),
                ground_stage=job.get("ground_stage", "after"),
                progress=progress)
        except Exception as exc:                      # one bad tile, not a job
            failed += 1
            emit({"t": "error", "path": cloud, "message": "%s" % exc})
            continue
        ok += 1
        # Paths are the point of the result; everything else is for the log.
        emit(dict({"t": "result"}, **{k: v for k, v in result.items()
                                      if _sendable(v)}))

    emit({"t": "done", "ok": ok, "failed": failed})
    return 0 if failed == 0 else 1


def _sendable(value):
    """Is this JSON-serialisable without special handling?"""
    return isinstance(value, (str, int, float, bool, type(None), list, dict))


if __name__ == "__main__":
    raise SystemExit(main())
