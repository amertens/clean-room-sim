# cleanroomGov

Governance and estimand-reporting templates for the staged, outcome-blind
clean-room workflow. These helpers **record and format** the analytic
specification a clean-room review requires; they do not estimate any causal
quantity and have no dependency on the estimation engine.

This package was split out of [`cleanTMLE`](https://github.com/amertens/cleanTMLE)
so that the estimation package carries only estimators, diagnostics, and the
enforcing pre-outcome gate, while this package holds the structured note-taking
layer. See the audit action items (section Z of the clean-room-sim `TODO.md`)
for the rationale and the planned second phase.

## What's here (Phase 1)

Pure, base-R formatters and declarations:

- `clean_event_process_table()` — classify post-baseline events/processes by
  their role in the estimand, aligned with ICH E9(R1) intercurrent-event
  strategies.
- `clean_check_event_processes()` — coherence check on event-of-interest,
  competing-event, and composite cumulative incidences.
- `clean_target_population()` — explicit target-population declaration.
- `clean_missing_data_plan()` — per-process missing-data plan, separating
  missingness from censoring.
- `clean_risk_report_table()` — cumulative-risk reporting table (formats
  estimator output; does not estimate).

## Status

Phase 1 of the estimation-vs-governance split. Not yet published to its own
repository; once it is, `cleanTMLE` will reference it via `Suggests`, and the
driver-coupled governance functions (`attach_estimand`, the decision-log
family, `record_stage`) can move here in Phase 2.

## Installation

```r
# install.packages("devtools")
devtools::install_local("cleanroomGov")
```

## Licence

MIT.
