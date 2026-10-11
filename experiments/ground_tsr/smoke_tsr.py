import sys
sys.path.insert(0, r"C:\Replicalm\src")
from replicalm.config import PRESETS
from replicalm import ground_tsr as T

c = PRESETS["clear"]
print("ground_engine default:", c.ground_engine)
print("overrides default:", c.ground_tsr_overrides)
print("baseline cell_m:", c.passes[0].cell_m)
p = T.GroundTSRParams()
print("gate above %.2f, below %.2f, slope allowance %.2f, prune %s"
      % (p.gate_k_above, p.gate_k_below, p.slope_allowance,
         p.prune_confident))
print("density_ref %.1f" % p.density_ref)
for d in (40.0, 20.0, 10.0, 5.0, 1.0):
    q = p.scale_for_density(d)
    print("  %6.1f/m2 -> smrf_cell %.2f  surface_cell %.2f  mcc %s"
          % (d, q.smrf_cell_m, q.surface_cell_m,
             [round(v, 2) for v in q.mcc_scales_m]))

# the override path the pipeline uses
import dataclasses
bad = {"nonsense_setting": 1}
ok = {"gate_k_above": 0.5}
print("")
print("override validation: known %s, unknown %s"
      % (all(k in p.__dataclass_fields__ for k in ok),
         [k for k in bad if k not in p.__dataclass_fields__]))
