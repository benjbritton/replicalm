"""Replicalm — a small desktop front end.

Deliberately thin. Every decision lives in `replicalm.pipeline` and
`replicalm.batch`; this window chooses inputs, shows progress, and reports what
was written. A front end that also decides things is one whose behaviour cannot
be reproduced from a script, and the whole point of the locked baseline is that
a result can be traced to a configuration.

tkinter rather than PySide: it ships with Python, adds nothing to the installer,
and this needs a dozen widgets and a log pane.

Processing runs in a CHILD PROCESS, not a thread. It used to be a thread, and
long numpy and PDAL operations hold the interpreter lock for stretches of
several seconds, so the main thread could not service its own timer: the window
stopped repainting and Windows appended "(Not Responding)" to the title. That
was accurate rather than a false alarm. In a child process nothing competes for
this process's lock, so the window stays live for the whole run -- and a child
can be terminated, which is what makes Cancel possible at all. A Python thread
cannot be stopped part-way through a kriging.

The protocol is one JSON object per line on the child's stdout; see
`replicalm.worker`. The only thread here reads those lines and does nothing
else, so it cannot starve the window either.

WHAT CHANGED SINCE 1.0
----------------------
Four methods became two. `ncalm` and `deep` remain preset keys, so an old
`_config.json` can still be identified by `profile_of`, but they are no longer
offered: `deep` was `clear` with a finer grid, which the resolution control
already does, and `ncalm` differed only in using the superseded cleanup rule.
Four near-identical choices made the window harder to use without making the
output better.

Resolution is in pixels per square metre rather than metres per cell. The
relationship is quadratic, so the cell size each choice means is shown beside
it -- 4 px/m2 is a 0.50 m cell, not a 0.25 m one -- and only perfect squares
are offered, since those are the values that land on round cell sizes.

Input accepts files and folders, expanded through `batch.collect_inputs` so the
window finds exactly the tiles a script would. The single-file picker was the
biggest practical limit in 1.0: the batch machinery existed and nothing in the
window could reach it.

The G1 checkbox now produces an image. It used to call a script that is not in
the installer, and silently did nothing on every machine but the one it was
written on.
"""
import os
import queue
import sys
import threading

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

HERE = os.path.dirname(os.path.abspath(__file__))
for candidate in (os.path.join(HERE, "src"),
                  os.path.join(os.path.dirname(HERE), "src")):
    if os.path.isdir(candidate) and candidate not in sys.path:
        sys.path.insert(0, candidate)

APP = "Replicalm"

# The two methods, as the window names them: preset key and the one line shown
# beside the selector. The descriptions are the whole explanation most users
# will read, so they say what the choice does to the output, not how.
PROFILES = {
    "Clear": ("clear", "removes low vegetation the ground filter keeps"),
    "Base":  ("baseline", "the published method, nothing extra removed"),
}

# Resolution, as pixels per square metre against the cell size it means.
RESOLUTIONS = [
    ("auto — matches the point density", None),
    ("1 px/m² — 1.00 m cells", 1),
    ("4 px/m² — 0.50 m cells", 4),
    ("9 px/m² — 0.33 m cells", 9),
    ("16 px/m² — 0.25 m cells", 16),
    ("25 px/m² — 0.20 m cells", 25),
]

# The interpolation neighbourhood: two measured choices rather than four free
# numbers, because the pair moves together -- 64 points only condition well
# inside a radius wide enough to hold them.
#
# 64 points within 20 m is first and is the default. It is the source method's
# own figure, and a smaller neighbourhood cuts spurious narrow grooves into
# every slope: on one tile the 99th percentile of local slope was 32.9 degrees
# with 16 points against the reference's 22, and 23.9 with 64 inside 20 m.
#
# 16 was briefly the default on 2026-10-05, chosen on tile-average measures that
# cannot see a fault occupying 1% of cells arranged in thin lines. It stays on
# the menu because it interpolates seven to nine times faster, and the label
# says what that costs rather than calling it "standard".
NEIGHBOURHOODS = [
    ("published — 64 points within 20 m", "published"),
    ("fast — 16 points, quicker but rougher on slopes", "fast"),
]

POINT_CLOUD_TYPES = [("Point clouds", "*.las *.laz"), ("All files", "*.*")]


def cell_from_pixels(px_per_m2):
    """Pixels per square metre -> metres per cell. 4 px/m2 is a 0.50 m cell."""
    return (1.0 / float(px_per_m2)) ** 0.5


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP)
        self.minsize(780, 500)
        self.messages = queue.Queue()
        self.worker = None          # the reader thread, not the processing
        self.proc = None            # the child doing the work
        self.cancelled = False
        self.inputs = []            # files and folders, in the order chosen
        self._build()
        # The ground-stage selector is meaningless until saving is asked for.
        self._toggle_ground_stage()
        self.after(120, self._drain)
        self.protocol("WM_DELETE_WINDOW", self._close)

    # ---------------------------------------------------------------- layout
    def _build(self):
        pad = dict(padx=10, pady=5)
        frm = ttk.Frame(self)
        frm.pack(fill="x", **pad)

        self.out = tk.StringVar()
        self.profile = tk.StringVar(value="Clear")
        self.resolution = tk.StringVar(value=RESOLUTIONS[0][0])
        self.neighbourhood = tk.StringVar(value=NEIGHBOURHOODS[0][0])
        self.threshold = tk.StringVar(value="auto")
        self.reason = tk.StringVar()
        # The G1 composite is the point of the tool for most users, so
        # it is on unless they turn it off. It adds a couple of minutes
        # to a tile; the elevation model is written either way.
        self.g1 = tk.BooleanVar(value=True)
        self.keep = tk.BooleanVar(value=False)
        # Which ground to save. "after" is the surface the DEM was built from;
        # "before" is the filter's own output, which is where re-deriving the
        # cleanup at another threshold starts without re-classifying.
        self.ground_stage = tk.StringVar(value="after cleanup")

        # --- inputs: any number of files, any number of folders -------------
        ttk.Label(frm, text="Input").grid(row=0, column=0, sticky="nw")
        box = ttk.Frame(frm)
        box.grid(row=0, column=1, sticky="ew")
        self.input_list = tk.Listbox(box, height=4, activestyle="none")
        self.input_list.pack(side="left", fill="both", expand=True)
        scroll = ttk.Scrollbar(box, orient="vertical",
                               command=self.input_list.yview)
        scroll.pack(side="left", fill="y")
        self.input_list.configure(yscrollcommand=scroll.set)

        picks = ttk.Frame(frm)
        picks.grid(row=0, column=2, sticky="nw", padx=(6, 0))
        ttk.Button(picks, text="Add file(s)…", width=13,
                   command=self._add_files).pack(fill="x")
        ttk.Button(picks, text="Add folder…", width=13,
                   command=self._add_folder).pack(fill="x", pady=(3, 0))
        ttk.Button(picks, text="Remove", width=13,
                   command=self._remove_selected).pack(fill="x", pady=(3, 0))

        ttk.Label(frm, text="Output folder").grid(row=1, column=0, sticky="w",
                                                  pady=(8, 0))
        ttk.Entry(frm, textvariable=self.out, width=62).grid(
            row=1, column=1, sticky="ew", pady=(8, 0))
        ttk.Button(frm, text="Browse…", command=self._pick_out).grid(
            row=1, column=2, sticky="ew", padx=(6, 0), pady=(8, 0))

        # Clear is the default because it is the measured improvement. Base
        # stays reachable for terrain where the cleanup threshold has not been
        # checked, and is what the published comparison figures are measured
        # against.
        meth = ttk.Frame(frm)
        meth.grid(row=2, column=1, sticky="w", pady=(10, 0))
        ttk.Label(meth, text="Method").pack(side="left")
        sel = ttk.Combobox(meth, textvariable=self.profile, width=8,
                           state="readonly", values=tuple(PROFILES))
        sel.pack(side="left", padx=(6, 8))
        sel.bind("<<ComboboxSelected>>", self._describe_profile)
        self.profile_note = ttk.Label(meth, text=PROFILES["Clear"][1],
                                      foreground="#444")
        self.profile_note.pack(side="left")

        res = ttk.Frame(frm)
        res.grid(row=3, column=1, sticky="w", pady=(8, 0))
        ttk.Label(res, text="Resolution").pack(side="left")
        ttk.Combobox(res, textvariable=self.resolution, width=30,
                     state="readonly",
                     values=tuple(label for label, _ in RESOLUTIONS)).pack(
                         side="left", padx=(6, 18))
        ttk.Checkbutton(res, text="Build the G1 image",
                        variable=self.g1).pack(side="left")

        nb = ttk.Frame(frm)
        nb.grid(row=4, column=1, sticky="w", pady=(8, 0))
        ttk.Label(nb, text="Interpolation").pack(side="left")
        ttk.Combobox(nb, textvariable=self.neighbourhood, width=34,
                     state="readonly",
                     values=tuple(label for label, _ in NEIGHBOURHOODS)).pack(
                         side="left", padx=(6, 0))

        opts = ttk.Frame(frm)
        opts.grid(row=5, column=1, sticky="w", pady=(8, 0))
        ttk.Label(opts, text="Ground threshold (m)").pack(side="left")
        ttk.Entry(opts, textvariable=self.threshold, width=7).pack(
            side="left", padx=(6, 18))
        # "Keep ground points" said nothing about what it kept or where it went,
        # and the after/before selector sat beside it at equal weight meaning
        # nothing to anyone who did not already know what the cleanup was. The
        # selector now only appears once the box is ticked.
        ttk.Checkbutton(opts, text="Save ground-filtered point cloud (LAZ)",
                        variable=self.keep,
                        command=self._toggle_ground_stage).pack(side="left")
        self.stage_box = ttk.Combobox(opts, textvariable=self.ground_stage,
                                      width=18, state="readonly",
                                      values=("after cleanup",
                                              "before cleanup"))
        self.stage_box.pack(side="left", padx=(6, 0))
        frm.columnconfigure(1, weight=1)

        # Whatever is typed here is stored in the settings file beside the
        # output, so a result can carry the operator's own account of why it was
        # produced that way. It used to be called "Reason for changing" and a
        # modal appeared if it was left empty while anything differed from the
        # defaults -- which interrupted the person who chose those defaults to
        # imply they had done something questionable. The field earns its place;
        # the interrogation did not.
        rs = ttk.Frame(frm)
        rs.grid(row=6, column=1, sticky="ew", pady=(6, 0))
        ttk.Label(rs, text="Processing notes (optional)").pack(side="left")
        ttk.Entry(rs, textvariable=self.reason).pack(side="left", fill="x",
                                                     expand=True, padx=(6, 0))

        runbar = ttk.Frame(self)
        runbar.pack(fill="x", padx=10)
        self.run_btn = ttk.Button(runbar, text="Run", command=self._run)
        self.run_btn.pack(side="left")
        self.cancel_btn = ttk.Button(runbar, text="Cancel",
                                     command=self._cancel, state="disabled")
        self.cancel_btn.pack(side="left", padx=(6, 0))
        self.progress = ttk.Progressbar(runbar, mode="determinate",
                                        maximum=1.0)
        self.progress.pack(side="left", fill="x", expand=True, padx=10)

        self.log = tk.Text(self, height=16, wrap="word", state="disabled",
                           font=("Consolas", 9))
        self.log.pack(fill="both", expand=True, padx=10, pady=(8, 4))

        self.status = ttk.Label(self, anchor="w", text="Ready")
        self.status.pack(fill="x", padx=10, pady=(0, 8))
        self._say("%s — bare-earth processing at the locked baseline." % APP)
        self._say("Add point cloud files or folders, choose an output folder, "
                  "then Run.")

    # ---------------------------------------------------------------- helpers
    def _describe_profile(self, _event=None):
        self.profile_note.configure(text=PROFILES[self.profile.get()][1])

    def _toggle_ground_stage(self):
        """Only show which ground to save once saving it has been asked for."""
        self.stage_box.configure(
            state="readonly" if self.keep.get() else "disabled")

    def _close(self):
        """Closing the window must not leave a child processing in the dark."""
        if self.proc is not None and self.proc.poll() is None:
            if not messagebox.askokcancel(
                    APP, "Processing is still running. Close anyway and stop "
                         "it?"):
                return
            self.cancelled = True
            try:
                self.proc.terminate()
            except Exception:
                pass
        self.destroy()

    def _add_inputs(self, paths):
        added = 0
        for p in paths:
            p = os.path.abspath(p)
            if p not in self.inputs:
                self.inputs.append(p)
                self.input_list.insert(
                    "end", ("[folder]  " if os.path.isdir(p) else "") + p)
                added += 1
        if added and not self.out.get():
            first = self.inputs[0]
            root = first if os.path.isdir(first) else os.path.dirname(first)
            self.out.set(os.path.join(root, "replicalm_out"))
        return added

    def _add_files(self):
        paths = filedialog.askopenfilenames(title="Choose point clouds",
                                            filetypes=POINT_CLOUD_TYPES)
        if paths:
            n = self._add_inputs(paths)
            self._say("added %d file%s" % (n, "" if n == 1 else "s"))

    def _add_folder(self):
        p = filedialog.askdirectory(
            title="Choose a folder of point clouds (searched recursively)")
        if p and self._add_inputs([p]):
            n = self._count_clouds(p)
            self._say("added folder: %s — %d point cloud%s found"
                      % (p, n, "" if n == 1 else "s"))

    @staticmethod
    def _count_clouds(folder):
        n = 0
        for _root, _dirs, files in os.walk(folder):
            n += sum(1 for f in files
                     if os.path.splitext(f)[1].lower() in (".las", ".laz"))
        return n

    def _remove_selected(self):
        for i in sorted(self.input_list.curselection(), reverse=True):
            self.input_list.delete(i)
            del self.inputs[i]

    def _pick_out(self):
        p = filedialog.askdirectory(title="Choose an output folder")
        if p:
            self.out.set(p)

    def _say(self, text):
        self.log.configure(state="normal")
        self.log.insert("end", text.rstrip() + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def _drain(self):
        while True:
            try:
                kind, payload = self.messages.get_nowait()
            except queue.Empty:
                break
            if kind == "log":
                self._say(payload)
            elif kind == "status":
                self.status.configure(text=payload)
            elif kind == "fraction":
                self.progress["value"] = payload
            elif kind == "done":
                # Back to exactly the pre-Run state, with every selection the
                # operator made still filled in. A cancelled run (payload None)
                # leaves the bar where it stopped rather than claiming either
                # completion or failure.
                self.run_btn.configure(state="normal")
                self.cancel_btn.configure(state="disabled")
                if payload is not None:
                    self.progress["value"] = 1.0 if payload else 0.0
        self.after(120, self._drain)

    # ------------------------------------------------------------------- run
    def _run(self):
        if not self.inputs:
            messagebox.showerror(APP, "Add at least one point cloud file or "
                                      "folder first.")
            return
        out = self.out.get().strip()
        if not out:
            messagebox.showerror(APP, "Choose an output folder first.")
            return

        px = dict(RESOLUTIONS)[self.resolution.get()]
        cell = None if px is None else cell_from_pixels(px)

        # The ground filter's elevation tolerance. 'auto' leaves it at the
        # selected method's default, which is the measured value.
        raw = self.threshold.get().strip().lower()
        if raw in ("", "auto"):
            threshold = None
        else:
            try:
                threshold = float(raw)
                if threshold <= 0:
                    raise ValueError
            except ValueError:
                messagebox.showerror(
                    APP, "Ground threshold must be a positive number of "
                         "metres, or 'auto' to use the method's default.")
                return

        neighbourhood = dict(NEIGHBOURHOODS)[self.neighbourhood.get()]

        self._say("")
        self._start(list(self.inputs), out, cell, self.g1.get(),
                    self.keep.get(), self.profile.get(), threshold,
                    neighbourhood, self.reason.get().strip(),
                    self.ground_stage.get())

    # ------------------------------------------------------------------ work
    @staticmethod
    def _configure(profile, threshold, neighbourhood, reason):
        """The config for a run, and the lines that explain it to the log."""
        from dataclasses import replace
        from replicalm.config import BASELINE, PRESETS
        preset, note = PROFILES[profile]
        cfg = PRESETS[preset]
        if threshold is not None:
            cfg = replace(cfg, passes=[replace(cfg.passes[0],
                                               threshold_m=threshold)]
                               + list(cfg.passes[1:]))
        if neighbourhood == "fast":
            cfg = replace(cfg, max_points=16, min_points=3,
                          search_radius_mode="density", search_radius_m=None)
        if reason:
            cfg = replace(cfg, note=reason)
        fixed = (cfg.search_radius_mode == "fixed" and cfg.search_radius_m)
        lines = [
            "baseline locked %s" % BASELINE["locked"],
            "method: %s — %s" % (profile, note),
            "interpolation: %d points, %s"
            % (cfg.max_points, "%.0f m fixed" % cfg.search_radius_m if fixed
               else "radius from density"),
            "ground threshold: %.2f m%s"
            % (cfg.passes[0].threshold_m,
               "" if threshold is not None else " (method default)"),
        ]
        if reason:
            lines.append("reason recorded: %s" % reason)
        return cfg, lines

    def _child_python(self):
        """python.exe beside this interpreter.

        The launcher runs the window under pythonw.exe so no console appears
        behind it. A child started from pythonw has no usable stdout, and the
        whole protocol is stdout, so the console build is used for the child.
        """
        exe = sys.executable or ""
        for a, b in (("pythonw.exe", "python.exe"), ("pythonw", "python")):
            if exe.endswith(a):
                candidate = exe[:-len(a)] + b
                if os.path.exists(candidate):
                    return candidate
        return exe

    def _start(self, inputs, out, cell, make_g1, keep, profile, threshold,
               neighbourhood, reason, ground_stage_label):
        """Write the job, launch the child, and read it on a reader thread."""
        import json
        import subprocess
        import tempfile

        ground_stage = ("before" if str(ground_stage_label).startswith("before")
                        else "after")
        try:
            from replicalm import batch
            cfg, lines = self._configure(profile, threshold, neighbourhood,
                                         reason)
            clouds = batch.collect_inputs(inputs)
        except Exception as exc:
            self._say("FAILED before starting: %s" % exc)
            return
        if not clouds:
            self._say("no .las or .laz files found in what was added")
            self.status.configure(text="Nothing to do")
            return
        for line in lines:
            self._say(line)
        self._say("%d point cloud%s to process"
                  % (len(clouds), "" if len(clouds) == 1 else "s"))
        self._say("")

        self._jobdir = tempfile.mkdtemp(prefix="replicalm_job_")
        cfg_path = os.path.join(self._jobdir, "config.json")
        cfg.to_json(cfg_path)
        job_path = os.path.join(self._jobdir, "job.json")
        with open(job_path, "w", encoding="utf-8") as fh:
            json.dump({"config": cfg_path, "out": out, "clouds": clouds,
                       "cell_m": cell, "make_g1": bool(make_g1),
                       "keep_ground": bool(keep),
                       "ground_stage": ground_stage}, fh)

        src = os.path.join(HERE, "src")
        env = dict(os.environ)
        if os.path.isdir(src):
            env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")
        # CREATE_NO_WINDOW so the child does not flash a console window. The
        # flag only exists on Windows.
        flags = getattr(subprocess, "CREATE_NO_WINDOW", 0)
        try:
            self.proc = subprocess.Popen(
                [self._child_python(), "-u", "-m", "replicalm.worker",
                 job_path],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding="utf-8", errors="replace",
                env=env, creationflags=flags)
        except Exception as exc:
            self._say("FAILED to start the processing job: %s" % exc)
            return

        self.cancelled = False
        self.run_btn.configure(state="disabled")
        self.cancel_btn.configure(state="normal")
        self.progress["value"] = 0.0
        self.worker = threading.Thread(target=self._read_child, daemon=True)
        self.worker.start()

    def _read_child(self):
        """Blocking readline on the child, forwarded to the window's queue.

        This thread does nothing but read, so it cannot hold the interpreter
        lock long enough to stop the window repainting -- which is the whole
        reason the processing moved out of process.
        """
        import json
        q = self.messages
        proc = self.proc
        try:
            for line in proc.stdout:
                line = line.strip()
                if not line:
                    continue
                try:
                    msg = json.loads(line)
                except ValueError:
                    # Not ours. Show it rather than swallow it: a message nobody
                    # sees is how the image step stayed broken for a release.
                    q.put(("log", line))
                    continue
                self._handle(q, msg)
        except Exception as exc:
            q.put(("log", "lost contact with the processing job: %s" % exc))
        code = proc.wait()
        if self.cancelled:
            q.put(("log", ""))
            q.put(("log", "Cancelled. Whatever had already finished is on "
                          "disk; nothing further was written."))
            q.put(("status", "Cancelled"))
            q.put(("done", None))
        elif code not in (0, 1):
            q.put(("log", "the processing job stopped unexpectedly (exit %s)"
                          % code))
            q.put(("status", "Failed"))
            q.put(("done", False))

    def _handle(self, q, msg):
        """One message from the child onto the window's queue."""
        kind = msg.get("t")
        if kind == "progress":
            q.put(("log", "[%-11s] %s" % (msg.get("stage", ""),
                                          msg.get("message", ""))))
            q.put(("status", msg.get("message", "")))
            if msg.get("fraction") is not None:
                q.put(("fraction", msg["fraction"]))
        elif kind == "file":
            q.put(("log", "--- %d of %d: %s"
                          % (msg["index"], msg["of"],
                             os.path.basename(msg["path"]))))
        elif kind == "result":
            self._report(q, msg)
        elif kind == "error":
            q.put(("log", "FAILED: %s" % msg.get("message", "")))
        elif kind == "done":
            q.put(("log", ""))
            total = msg.get("ok", 0) + msg.get("failed", 0)
            q.put(("status", "Finished %d of %d" % (msg.get("ok", 0), total)))
            q.put(("done", not msg.get("failed")))

    def _cancel(self):
        """Stop processing and put the window back as it was before Run.

        Terminating a child is the whole reason the work is in one. A thread
        could not be stopped part-way through a kriging at all.
        """
        proc = getattr(self, "proc", None)
        if proc is None or proc.poll() is not None:
            return
        self.cancelled = True
        self.cancel_btn.configure(state="disabled")
        self.status.configure(text="Stopping...")
        self._say("")
        self._say("stopping the processing job...")
        try:
            proc.terminate()
        except Exception as exc:
            self._say("could not stop it cleanly: %s" % exc)

    @staticmethod
    def _report(q, result):
        q.put(("log", "DEM written: %s" % result["dem"]))
        if result.get("ground_las"):
            q.put(("log", "Ground-filtered cloud (%s cleanup): %s"
                          % (result.get("ground_stage"), result["ground_las"])))
        q.put(("log", "%d ground points, %.2f per m\u00b2, cell %.2f m, "
                      "radius %.2f m"
                      % (result["ground_points"], result["density"],
                         result["cell_m"], result["search_radius_m"])))
        q.put(("log", "%d cells kept after trimming %d from the edge"
                      % (result["cells_after_trim"], result["erode_cells"])))
        if result.get("filled_cells"):
            q.put(("log", "%d cells interpolated across holes (%.3f%%)"
                          % (result["filled_cells"],
                             100 * result["filled_fraction"])))
        g1 = result.get("g1")
        if isinstance(g1, dict) and g1.get("g1"):
            q.put(("log", "G1 image: %s" % g1["g1"]))
            if g1.get("layers"):
                q.put(("log", "G1 source layers: %s"
                              % os.path.dirname(
                                  next(iter(g1["layers"].values())))))
        elif isinstance(g1, dict):
            q.put(("log", "image step not completed: %s"
                          % (g1.get("skipped") or g1.get("failed"))))
        q.put(("log", "finished in %.0f s" % result["seconds"]))



if __name__ == "__main__":
    App().mainloop()
