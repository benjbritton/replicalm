"""Replicalm — a small desktop front end.

Deliberately thin. Every decision lives in `replicalm.pipeline`; this window
chooses files, shows progress, and reports what was written. A front end that
also decides things is one whose behaviour cannot be reproduced from a script,
and the whole point of the locked baseline is that a result can be traced to a
configuration.

tkinter rather than PySide: it ships with Python, adds nothing to the installer,
and this needs four widgets and a log pane.

Processing runs on a worker thread. Classification and kriging on a full tile
take minutes, and a frozen window that may or may not be working is worse than
a slow one that says what it is doing.
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

# The three methods, as the window names them: preset key and the one line
# shown beside the selector. The descriptions are the whole explanation most
# users will read, so they say what the choice does to the output, not how.
PROFILES = {
    "Clear":    ("clear",    "removes low vegetation the ground filter keeps"),
    "Baseline": ("baseline", "the published method, nothing extra removed"),
    "Deep":     ("deep",     "Clear on a finer grid; slower, for viewing"),
}


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP)
        self.minsize(720, 430)
        self.messages = queue.Queue()
        self.worker = None
        self._build()
        self.after(120, self._drain)

    # ---------------------------------------------------------------- layout
    def _build(self):
        pad = dict(padx=10, pady=5)
        frm = ttk.Frame(self)
        frm.pack(fill="x", **pad)

        self.las = tk.StringVar()
        self.out = tk.StringVar()
        self.cell = tk.StringVar(value="auto")
        self.profile = tk.StringVar(value="Clear")
        self.threshold = tk.StringVar(value="auto")
        self.maxpts = tk.StringVar(value="auto")
        self.minpts = tk.StringVar(value="auto")
        self.g1 = tk.BooleanVar(value=False)
        self.keep = tk.BooleanVar(value=False)

        ttk.Label(frm, text="Point cloud").grid(row=0, column=0, sticky="w")
        ttk.Entry(frm, textvariable=self.las, width=62).grid(row=0, column=1,
                                                             sticky="ew")
        ttk.Button(frm, text="Browse…", command=self._pick_las).grid(row=0,
                                                                     column=2)

        ttk.Label(frm, text="Output folder").grid(row=1, column=0, sticky="w")
        ttk.Entry(frm, textvariable=self.out, width=62).grid(row=1, column=1,
                                                             sticky="ew")
        ttk.Button(frm, text="Browse…", command=self._pick_out).grid(row=1,
                                                                     column=2)

        # Clear is the default because it is the measured improvement, but the
        # cleanup threshold was fitted on one survey. Baseline stays reachable
        # for terrain where that threshold has not been checked.
        meth = ttk.Frame(frm)
        meth.grid(row=2, column=1, sticky="w", pady=(10, 0))
        ttk.Label(meth, text="Method").pack(side="left")
        box = ttk.Combobox(meth, textvariable=self.profile, width=10,
                           state="readonly",
                           values=("Clear", "Baseline", "Deep"))
        box.pack(side="left", padx=(6, 8))
        box.bind("<<ComboboxSelected>>", self._describe_profile)
        self.profile_note = ttk.Label(meth, text=PROFILES["Clear"][1],
                                      foreground="#444")
        self.profile_note.pack(side="left")

        opts = ttk.Frame(frm)
        opts.grid(row=3, column=1, sticky="w", pady=(8, 0))
        ttk.Label(opts, text="Cell size (m)").pack(side="left")
        ttk.Entry(opts, textvariable=self.cell, width=7).pack(side="left",
                                                              padx=(6, 4))
        ttk.Label(opts, text="'auto' matches the point density").pack(
            side="left", padx=(0, 18))
        ttk.Checkbutton(opts, text="Also build the G1 image",
                        variable=self.g1).pack(side="left", padx=(0, 18))

        ttk.Label(opts, text="Ground threshold (m)").pack(side="left",
                                                          padx=(0, 0))
        ttk.Entry(opts, textvariable=self.threshold, width=7).pack(
            side="left", padx=(6, 18))
        ttk.Checkbutton(opts, text="Keep classified points",
                        variable=self.keep).pack(side="left")
        frm.columnconfigure(1, weight=1)

        # The search neighbourhood. Left alone on most surveys: past the
        # variogram's correlation range extra neighbours carry almost no weight,
        # so the count only matters where that range is long relative to the
        # spacing between returns.
        nb = ttk.Frame(frm)
        nb.grid(row=4, column=1, sticky="w", pady=(6, 0))
        ttk.Label(nb, text="Search neighbours  max").pack(side="left")
        ttk.Entry(nb, textvariable=self.maxpts, width=7).pack(side="left",
                                                              padx=(6, 10))
        ttk.Label(nb, text="min").pack(side="left")
        ttk.Entry(nb, textvariable=self.minpts, width=7).pack(side="left",
                                                              padx=(6, 10))
        ttk.Label(nb, text="'auto' uses the method's own values",
                  foreground="#444").pack(side="left")

        bar = ttk.Frame(self)
        bar.pack(fill="x", padx=10)
        self.run_btn = ttk.Button(bar, text="Run", command=self._run)
        self.run_btn.pack(side="left")
        self.progress = ttk.Progressbar(bar, mode="determinate", maximum=1.0)
        self.progress.pack(side="left", fill="x", expand=True, padx=10)

        self.log = tk.Text(self, height=16, wrap="word", state="disabled",
                           font=("Consolas", 9))
        self.log.pack(fill="both", expand=True, padx=10, pady=(8, 4))

        self.status = ttk.Label(self, anchor="w", text="Ready")
        self.status.pack(fill="x", padx=10, pady=(0, 8))
        self._say("%s — bare-earth processing at the locked baseline." % APP)
        self._say("Choose a LAS or LAZ file and an output folder, then Run.")

    # ---------------------------------------------------------------- helpers
    def _describe_profile(self, _event=None):
        name = self.profile.get()
        self.profile_note.configure(text=PROFILES[name][1])
        # Deep chooses its own grid; a typed cell size would contradict it
        if name == "Deep" and self.cell.get().strip().lower() not in ("", "auto"):
            self.cell.set("auto")
            self._say("Deep sets its own cell size; reverted to auto.")

    def _pick_las(self):
        p = filedialog.askopenfilename(
            title="Choose a point cloud",
            filetypes=[("Point clouds", "*.las *.laz"), ("All files", "*.*")])
        if p:
            self.las.set(p)
            if not self.out.get():
                self.out.set(os.path.join(os.path.dirname(p), "replicalm_out"))

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
        las, out = self.las.get().strip(), self.out.get().strip()
        if not las or not os.path.isfile(las):
            messagebox.showerror(APP, "Choose a point cloud file first.")
            return
        if not out:
            messagebox.showerror(APP, "Choose an output folder first.")
            return
        # "auto" lets the pipeline derive the cell from measured ground-return
        # density -- 1/sqrt(density), which is where the residual curve turns.
        # A fixed default is only right for one density.
        raw = self.cell.get().strip().lower()
        if raw in ("", "auto"):
            cell = None
        else:
            try:
                cell = float(raw)
                if cell <= 0:
                    raise ValueError
            except ValueError:
                messagebox.showerror(
                    APP, "Cell size must be a positive number, or 'auto' to "
                         "match it to the point density.")
                return

        # The ground filter's elevation tolerance. 'auto' leaves it at the
        # selected method's default -- 0.25 m for Clear and Deep, 0.50 m for
        # Baseline, which is the published figure.
        raw_t = self.threshold.get().strip().lower()
        if raw_t in ("", "auto"):
            threshold = None
        else:
            try:
                threshold = float(raw_t)
                if threshold <= 0:
                    raise ValueError
            except ValueError:
                messagebox.showerror(
                    APP, "Ground threshold must be a positive number of "
                         "metres, or 'auto' to use the method's default.")
                return

        def whole(var, label):
            raw = var.get().strip().lower()
            if raw in ("", "auto"):
                return None
            try:
                v = int(raw)
                if v < 1:
                    raise ValueError
                return v
            except ValueError:
                messagebox.showerror(
                    APP, "%s must be a whole number of at least 1, or 'auto' "
                         "to use the method's own value." % label)
                raise ValueError
        try:
            maxpts = whole(self.maxpts, "Maximum search neighbours")
            minpts = whole(self.minpts, "Minimum search neighbours")
        except ValueError:
            return
        if maxpts is not None and minpts is not None and minpts > maxpts:
            messagebox.showerror(APP, "Minimum search neighbours cannot exceed "
                                      "the maximum.")
            return

        self.run_btn.configure(state="disabled")
        self.progress["value"] = 0.0
        self._say("")
        self.worker = threading.Thread(
            target=self._work, args=(las, out, cell, self.g1.get(),
                                     self.keep.get(), self.profile.get(),
                                     threshold, maxpts, minpts),
            daemon=True)
        self.worker.start()

    def _work(self, las, out, cell, make_g1, keep, profile, threshold,
              maxpts, minpts):
        q = self.messages

        def progress(stage, message, fraction=None):
            q.put(("log", "[%-11s] %s" % (stage, message)))
            q.put(("status", message))
            if fraction is not None:
                q.put(("fraction", fraction))

        try:
            from dataclasses import replace
            from replicalm import pipeline
            from replicalm.config import BASELINE, PRESETS
            preset, note = PROFILES[profile]
            cfg = PRESETS[preset]
            if threshold is not None:
                cfg = replace(cfg, passes=[replace(cfg.passes[0],
                                                   threshold_m=threshold)]
                                          + list(cfg.passes[1:]))
            if maxpts is not None or minpts is not None:
                cfg = replace(cfg,
                              max_points=maxpts if maxpts is not None
                              else cfg.max_points,
                              min_points=minpts if minpts is not None
                              else cfg.min_points)
            q.put(("log", "baseline locked %s" % BASELINE["locked"]))
            q.put(("log", "method: %s -- %s" % (profile, note)))
            q.put(("log", "search neighbours: %d max, %d min"
                          % (cfg.max_points, cfg.min_points)))
            q.put(("log", "ground threshold: %s m"
                          % ("%.2f" % threshold if threshold is not None
                             else "%.2f (method default)" % cfg.passes[0].threshold_m)))
            result = pipeline.process(las, out, cfg=cfg,
                                      cell_m=cell, make_g1=make_g1,
                                      keep_ground=keep, progress=progress)
            q.put(("log", ""))
            q.put(("log", "DEM written: %s" % result["dem"]))
            q.put(("log", "%d ground points, %.2f per m2, cell %.2f m, "
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
            if isinstance(g1, dict) and g1.get("rvt_root"):
                q.put(("log", "G1 and RVT layers: %s" % g1["rvt_root"]))
            elif isinstance(g1, dict):
                q.put(("log", "image step skipped: %s"
                              % (g1.get("skipped") or g1.get("stderr"))))
            q.put(("status", "Finished in %.0f s" % result["seconds"]))
            q.put(("done", True))
        except Exception as exc:
            q.put(("log", ""))
            q.put(("log", "FAILED: %s" % exc))
            q.put(("log", traceback.format_exc(limit=3)))
            q.put(("status", "Failed"))
            q.put(("done", False))


if __name__ == "__main__":
    App().mainloop()
