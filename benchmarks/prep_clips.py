r"""Cut a 400 m window from each calibration tile, cloud and reference alike."""
import json, sys, time
sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import clip

with open(r"C:\Replicalm\tests\calibration_tiles.json", encoding="utf-8") as fh:
    sel = json.load(fh)

OUT = r"C:\Replicalm\tests\clips"
recs = []
for t in sel["tiles"]:
    t0 = time.time()
    try:
        r = clip.prepare(t, OUT, win_m=400.0)
        r["seconds"] = round(time.time() - t0, 1)
        recs.append(r)
        print("  -> %d points, %.0f s\n" % (r["points"], time.time() - t0))
    except Exception as e:
        print("  FAILED %s: %s\n" % (t["tile"], e))

with open(r"C:\Replicalm\tests\clips\index.json", "w", encoding="utf-8") as fh:
    json.dump(recs, fh, indent=1)
print("prepared %d of %d tiles" % (len(recs), len(sel["tiles"])))
