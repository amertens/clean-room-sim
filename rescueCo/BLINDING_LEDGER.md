# Blinding ledger: Rescue.Co re-execution

Started 2026-10-18. What Andrew and the co-authors have seen of the outcomes, and when; what each design-stage input reveals about the outcomes; and every approval. Entries are added, never edited; a correction is a new entry.

## 1. What the team has seen of the outcomes

The re-execution is not an analysis by analysts blind to the outcomes. Every member of the team had seen effect estimates from the parallel unblinded analysis (`rescueCo_analysis`) before the re-execution's design stage was rerun.

| Date | Who | What was seen | Source |
|---|---|---|---|
| 2026-02-26 | Analysis team | First effect estimates on every outcome (matched analysis: GOSE, severe disability, death in the emergency department, death by six months, time to arrival) | `rescueCo_analysis` 54f0ef7, `analysis_outputs/matched_effects_summary.csv` |
| 2026-03-13 | Analysis team | TMLE estimates from the `clean_room` fork | `rescueCo_analysis` efc9d7a, `clean_room/results/binary_outcome_comparison.csv` |
| 2026-05-27 to 28 | Andrew and co-authors | Co-author brief and meeting, with unblinded estimates listed | `reports/archive/coauthor_meeting_brief_2026-05-27.md`, `meeting_decisions_2026-05-28.md` |
| 2026-09 | Andrew | Clean-room multi-arm estimates on the 0.2.0 development build; reconciliation against the unblinded analysis | `clean-room-sim/rescueCo/results/reconciliation_2026-09.md`, `reconciliation_table.csv` |
| 2026-10-15 | Andrew and co-authors | Results call: primary and wider-comparison estimates | inventory rows 4 and 7 |
| 2026-10-16 | Andrew and co-authors | Short results report, 2026-10-16a | `rescueCo_analysis/reports/rescueco_results_summary_2026-10-16a.docx` |

To add: the dates co-authors received each earlier report, if Andrew wants them recorded.

## 2. What each design-stage input reveals

| Input | Reveals | Approved |
|---|---|---|
| Covariates, treatment, timing, eligibility and cohort flags | Nothing about outcome values | Ground rule |
| Emergency-department disposition recorded (indicator) | Whether disposition was recorded; nothing about its value | D3, 2026-10-18 |
| Six-month status recorded (indicator) | A patient lost at six months is known to have left the emergency department alive, so this indicator is a lower bound on survival to disposition among the lost | D3, 2026-10-18 |
| Death at emergency-department disposition | Outcome: stays in the outcome store; enters only at `unblind()` | D3, 2026-10-18 |

No pooled outcome quantity has been approved for the design stage.

## 3. Approvals

| Date | Approval | By | Record |
|---|---|---|---|
| 2026-10-18 | Inventory of departures | Andrew | `rescueCo/deviations_inventory.md` |
| 2026-10-18 | Decisions D1 to D11 | Andrew | `docs/superpowers/specs/2026-10-18-outcome-blind-design-stage.md` |
