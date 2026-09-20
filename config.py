"""Replicalm configuration and baseline enforcement.

This is the front door. The values below are the locked baseline for
bare-earth processing, every one of them measured rather than chosen, and the
implementation lives in `src/replicalm/config.py`. Importing from here gives the
same objects; it exists so that the configuration is the first thing visible in
the repository rather than four directories down.

USAGE
-----
    import config
    cfg = config.load()              # the locked baseline
    config.verify_baseline(cfg)      # raises if anything has drifted

    cfg = config.load(dem_cell_m=0.5)    # override what you mean to override
    config.describe()                    # print the baseline and its evidence

WHAT "LOCKED" MEANS
-------------------
The baseline was fixed on 2026-09-20 after a recalibration on terrain that can
discriminate between settings. It is not a suggestion and it is not a set of
defaults that drifted into place: each parameter carries the measurement that
chose it, recorded in `src/replicalm/config.py` and in
`docs/open_observations.md`.

`verify_baseline` compares a configuration against the locked values and raises
listing every field that differs. Call it before a production run so that an
edited setting stops the run rather than quietly producing something else.

Changing a locked value is a deliberate act. It should come with a measurement
taken on terrain that can show the difference -- which, for this data, means
ground with relief. Three of the four parameters this baseline replaced were
chosen on sample windows containing between 0.00% and 0.20% of cells steeper
than twenty degrees, and all three were wrong.
"""
import os
import sys
from dataclasses import replace as _replace

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from replicalm.config import (           # noqa: E402
    BASELINE,
    GroundPass,
    NCALM_TERRASCAN,
    PRESETS,
    ReplicalmConfig,
    verify_baseline,
)

__all__ = ["BASELINE", "GroundPass", "NCALM_TERRASCAN", "PRESETS",
           "ReplicalmConfig", "verify_baseline", "load", "describe"]


def load(preset="ncalm", **overrides):
    """The locked baseline, with any overrides applied.

    Overrides are not validated against the baseline -- that is what
    `verify_baseline` is for, and separating the two keeps an experiment easy
    to run while keeping a production run honest about what it changed.
    """
    if preset not in PRESETS:
        raise KeyError("no preset %r; have %s" % (preset, ", ".join(PRESETS)))
    cfg = PRESETS[preset]
    return _replace(cfg, **overrides) if overrides else cfg


def describe():
    """Print the locked baseline and the source method it translates."""
    print("Replicalm baseline, locked %s" % BASELINE["locked"])
    print("-" * 52)
    for k, v in BASELINE.items():
        if k != "locked":
            print("  %-26s %s" % (k, v))
    print("\nTranslated from: %s" % NCALM_TERRASCAN["source"])
    print("  source ground classification: two TerraScan passes, "
          "70 m window, 88 deg max terrain angle,")
    print("  9 deg then 12 deg iteration angle, 3.0 m iteration distance")
    print("\nSee docs/processing_report.md for what is the same, similar and "
          "different,\nand docs/open_observations.md for what is not settled.")


if __name__ == "__main__":
    describe()
    print("\nbaseline verifies:", verify_baseline())
