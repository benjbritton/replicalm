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

        opts = ttk.Frame(frm)
        opts.grid(row=2, column=1, sticky="w", pady=(8, 0))
        ttk.Label(opts, text="Cell size (m)").pack(side="left")
        ttk.Entry(opts, textvariable=self.cell, width=7).pack(side="left",
                                                              padx=(6, 4))
        ttk.Label(opts, text="'auto' matches the point density").pack(
            side="left", padx=(0, 18))
        ttk.Checkbutton(opts, text="Also build the G1 image",
                        variable=self.g1).pack(side="left", padx=(0, 18))
        ttk.Checkbutton(opts, text="Keep classified points",
                        variable=self.keep).pack(side="left")
        frm.columnconfigure(1, weight=1)

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

        self.run_btn.configure(state="disabled")
        self.progress["value"] = 0.0
        self._say("")
        self.worker = threading.Thread(
            target=self._work, args=(las, out, cell, self.g1.get(),
                                     self.keep.get()), daemon=True)
        self.worker.start()

    def _work(self, las, out, cell, make_g1, keep):
        q = self.messages

        def progress(stage, message, fraction=None):
            q.put(("log", "[%-11s] %s" % (stage, message)))
            q.put(("status", message))
            if fraction is not None:
                q.put(("fraction", fraction))

        try:
            from replicalm import pipeline
            from replicalm.config import BASELINE
            q.put(("log", "baseline locked %s" % BASELINE["locked"]))
            result = pipeline.process(las, out, cell_m=cell, make_g1=make_g1,
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
