# Replicalm Pipeline Handoff Note

## Current State
- **Core Pipeline:** Relational/relief calibration active. Two-pass SMRF default locked in (9° and 12°, active noise filters, fixed 20m search radius).
- **Recent Additions:** evaluate.py (archive-independent metrics like buffered LOO, returns-below-surface, sharpness), inalise.py, and tokenized benchmarks under enchmarks/results/deep/.
- **Active Tracks:** "Clear" rules are validated but not forced into the default baseline yet. "Deep" track established that the resolution knee sits at mean point spacing.
- **Immediate Next Steps:** Reviewing benchmark reports and deciding whether to formally promote "Clear" or operationalize derived cell sizing.
