# Registration record for the Rescue.Co re-execution

Tracked by git and pushed before the outcome store is built (decision 4 of 2026-10-23, spec 2026-10-20, section 9), so that a dated copy of every hash exists off this machine before any outcome is scored. Patient-level files stay in `rescueCo/data/`; only hashes, approvals and aggregate tables come here.

| File | Written by | When |
|---|---|---|
| `build_fingerprint.csv` | `scripts/22c_record_build.R` | Once, before anything is installed into `rescueCo/lib` |
| `confirmatory.csv`, `lmtp_cohort.csv` | `scripts/23_registration_bundle.R pick` | When L1 to L5 of the final run all have design reports |
| `registration_locks.csv`, `population_gate.csv`, `registration_bundle.md`, `approvals_final_template.csv` | `scripts/23_registration_bundle.R final` | When every sensitivity and 180-day lock has a design report, or at the cut-off |
| `approvals_final.csv` | Andrew, from the template | After reviewing each design report; scripts 24 and 25 check it |
| `estimates_hashes.csv` | `scripts/25_estimate.R run` | After estimation, before the unblinded analysis is opened for the comparison |
