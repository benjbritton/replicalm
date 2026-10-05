"""Replicalm — a small desktop front end.

Deliberately thin. Every decision lives in `replicalm.pipeline` and
`replicalm.batch`; this window chooses inputs, shows progress, and reports what
was written. A front end that also decides things is one whose behaviour cannot
be reproduced from a script, and the whole point of the locked baseline is that
a result can be traced to a configuration.

tkinter rather than PySide: it ships with Python, adds nothing to the installer,
and this needs a dozen widgets and a log pane.

Processing runs on a worker thread. Classification and kriging on a full tile
take minutes to hours, and a frozen window that may or may not be working is
worse than a slow one that says what it is doing.

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
import traceback

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
NEIGHBOURHOODS = [
    ("published — 64 points within 20 m", "published"),
    ("fast — 16 points, radius from density", "fast"),
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
        self.worker = None
        self.inputs = []            # files and folders, in the order chosen
        self._build()
        self.after(120, self._drain)

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
        self.g1 = tk.BooleanVar(value=False)
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
        ttk.Checkbutton(res, text="Also build the G1 image",
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
        ttk.Checkbutton(opts, text="Keep ground points",
                        variable=self.keep).pack(side="left")
        ttk.Combobox(opts, textvariable=self.ground_stage, width=18,
                     state="readonly",
                     values=("after cleanup", "before cleanup")).pack(
                         side="left", padx=(6, 0))
        frm.columnconfigure(1, weight=1)

        # A departure from the defaults is a claim about this survey, so the
        # window asks for it in words and stores the answer beside the output.
        rs = ttk.Frame(frm)
        rs.grid(row=6, column=1, sticky="ew", pady=(6, 0))
        ttk.Label(rs, text="Reason for changing").pack(side="left")
        ttk.Entry(rs, textvariable=self.reason).pack(side="left", fill="x",
                                                     expand=True, padx=(6, 0))

        runbar = ttk.Frame(self)
        runbar.pack(fill="x", padx=10)
        self.run_btn = ttk.Button(runbar, text="Run", command=self._run)
        self.run_btn.pack(side="left")
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
                self.run_btn.configure(state="normal")
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

        # Not enforced -- a required box teaches people to type a full stop --
        # but asked for once, and recorded either way.
        changed = [n for n, v in (
            ("method", self.profile.get() != "Clear"),
            ("resolution", cell is not None),
            ("ground threshold", threshold is not None),
            ("interpolation", neighbourhood != NEIGHBOURHOODS[0][1]))
            if v]
        if changed and not self.reason.get().strip():
            gap = "\n\n"
            if not messagebox.askokcancel(
                    APP,
                    "You have changed %s from the default."
                    % " and ".join(changed) + gap
                    + "The defaults are what the measurements behind the "
                      "locked baseline were taken at. If this survey needs "
                      "something else, saying why in the Reason box records it "
                      "beside the output." + gap + "Run anyway?"):
                return

        self.run_btn.configure(state="disabled")
        self.progress["value"] = 0.0
        self._say("")
        self.worker = threading.Thread(
            target=self._work,
            args=(list(self.inputs), out, cell, self.g1.get(), self.keep.get(),
                  self.profile.get(), threshold, neighbourhood,
                  self.reason.get().strip(), self.ground_stage.get()),
            daemon=True)
        self.worker.start()

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

    def _work(self, inputs, out, cell, make_g1, keep, profile, threshold,
              neighbourhood, reason, ground_stage_label="after cleanup"):
        q = self.messages
        # the control reads in words; the pipeline takes a keyword
        ground_stage = ("before" if str(ground_stage_label).startswith("before")
                        else "after")

        def progress(stage, message, fraction=None):
            q.put(("log", "[%-11s] %s" % (stage, message)))
            q.put(("status", message))
            if fraction is not None:
                q.put(("fraction", fraction))

        try:
            from replicalm import batch, pipeline
            cfg, lines = self._configure(profile, threshold, neighbourhood,
                                         reason)
            for line in lines:
                q.put(("log", line))

            clouds = batch.collect_inputs(inputs)
            if not clouds:
                q.put(("log", "no .las or .laz files found in what was added"))
                q.put(("status", "Nothing to do"))
                q.put(("done", False))
                return
            q.put(("log", "%d point cloud%s to process"
                          % (len(clouds), "" if len(clouds) == 1 else "s")))
            q.put(("log", ""))

            failed = 0
            for i, cloud in enumerate(clouds, start=1):
                q.put(("log", "--- %d of %d: %s"
                              % (i, len(clouds), os.path.basename(cloud))))
                base = (i - 1) / len(clouds)
                span = 1.0 / len(clouds)

                def scaled(stage, message, fraction=None, _b=base, _s=span):
                    progress(stage, message,
                             None if fraction is None else _b + _s * fraction)

                try:
                    result = pipeline.process(cloud, out, cfg=cfg,
                                              cell_m=cell, make_g1=make_g1,
                                              keep_ground=keep,
                                              ground_stage=ground_stage,
                                              progress=scaled)
                except Exception as exc:
                    failed += 1
                    q.put(("log", "FAILED: %s" % exc))
                    continue
                self._report(q, result)

            q.put(("log", ""))
            q.put(("status", "Finished %d of %d"
                             % (len(clouds) - failed, len(clouds))))
            q.put(("done", failed == 0))
        except Exception as exc:
            q.put(("log", ""))
            q.put(("log", "FAILED: %s" % exc))
            q.put(("log", traceback.format_exc(limit=3)))
            q.put(("status", "Failed"))
            q.put(("done", False))

    @staticmethod
    def _report(q, result):
        q.put(("log", "DEM written: %s" % result["dem"]))
        if result.get("ground_las"):
            q.put(("log", "Ground (%s cleanup): %s"
                          % (result["ground_stage"], result["ground_las"])))
        q.put(("log", "%d ground points, %.2f per m², cell %.2f m, "
                      "radius %.2f m"
                      % (result["ground_points"], result["density"],
                         result["cell_m"], result["search_radius_m"])))
        q.put(("log", "%d cells kept after trimming %d from the edge"
                      % (result["cells_after_trim"], result["erode_cells"])))
        if result["filled_cells"]:
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
