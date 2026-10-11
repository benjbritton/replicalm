"""Is the G1 on by default in all three entry points, and still optional?"""
import re
import sys

sys.path.insert(0, r"C:\Replicalm\src")

print("WINDOW")
src = open(r"C:\Replicalm\packaging\gui.py", encoding="utf-8").read()
m = re.search(r"self\.g1 = tk\.BooleanVar\(value=(\w+)\)", src)
print("  checkbox default: %s" % m.group(1))
lbl = re.search(r'Checkbutton\(res, text="([^"]+)"', src)
print("  label: %s" % lbl.group(1))
assert m.group(1) == "True"

print("")
print("SINGLE-FILE COMMAND LINE")
from replicalm import pipeline
ap = None
for name in ("_argparser", "_parser", "main"):
    if hasattr(pipeline, name):
        break
# drive argparse the way the CLI does
import importlib
import argparse


def probe(mod_path, argv):
    import subprocess
    code = (
        "import sys;sys.path.insert(0,r'C:\\Replicalm\\src');"
        "import %s as M;"
        "p=M._build_parser() if hasattr(M,'_build_parser') else None;"
        "print('no parser factory' if p is None else p.parse_args(%r).g1)"
        % (mod_path, argv))
    return subprocess.run([sys.executable, "-c", code],
                          capture_output=True, text=True).stdout.strip()


for mod in ("replicalm.pipeline", "replicalm.batch"):
    out_default = probe(mod, ["x.las", "-o", "out"])
    out_off = probe(mod, ["x.las", "-o", "out", "--no-g1"])
    print("  %-20s default %-24s --no-g1 %s"
          % (mod, out_default or "(n/a)", out_off or "(n/a)"))

print("")
print("BATCH run() SIGNATURE")
import inspect
from replicalm import batch
sig = inspect.signature(batch.run)
print("  make_g1 default: %s" % sig.parameters["make_g1"].default)
assert sig.parameters["make_g1"].default is True

print("")
print("LIBRARY ENTRY POINT (unchanged on purpose)")
sigp = inspect.signature(pipeline.process)
print("  pipeline.process make_g1 default: %s"
      % sigp.parameters["make_g1"].default)
print("  callers pass it explicitly; the front ends decide the default")
