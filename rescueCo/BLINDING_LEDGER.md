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

## 4. Entries

| Date | Entry |
|---|---|
| 2026-10-20 | The steward script (20) linked the follow-up file to the registry (7,784 linked, as in the unblinded analysis), scored the outcomes in memory, and wrote only two indicators per patient: disposition recorded, and six-month status recorded (the primary outcome is known, counting a death at disposition as known). Registry-wide totals were logged: disposition recorded 8,281 of 8,326; six-month status recorded 6,519 of 8,326. No outcome value was written, printed or logged. |
| 2026-10-20 | What the six-month indicator reveals: a patient with status not recorded did not die at emergency-department disposition (a death there is always recorded). It does not separate deaths from survivors among patients with status recorded. |
| 2026-10-20 | The analyst's design build (21) read registry columns from a whitelist and never opened the follow-up file; it joined the two indicators by patient id. |
| 2026-10-22 | Three automated agents in the read-only audit of 22 October (`audit/2026-10-22/`) read parts of `rescueCo/logs/pipeline.log`. The audit of the package's estimation code printed the first 40 lines while looking for run times; by its own note, those lines come from the earlier unblinded analysis and include outcome counts (`audit/2026-10-22/02-wrapper.md:264`). The audit of the simulation study cited two timing lines (`pipeline.log:1204-1205`, in `03-simulation.md:162` and `:350`), and the check of the audit's recommendations cited one progress line (`pipeline.log:1226`, in `06-check-guarantees.md:17`). No outcome value appears in any audit file (checked with every digit masked), no person read the log through the audit, and no recommendation about the re-execution's design rests on it; the timings were replaced by `rescueCo/logs/design_stage2_stderr.txt:3`. Whether the lines cited by the second and third agents hold outcome information was not checked, because the log was not opened again. Cause: the audit's instructions listed the outcome-bearing folders but not this log. The log is added to section 5. Recorded 2026-10-23. |
| 2026-10-23 | Before any outcome store exists, Andrew decided how the steps after the design stage run (spec 2026-10-20, section 9). Script 25 attempts every lock, estimates each lock that `unblind()` admits without an override with the candidate the simulation selected, and reports each refusal as a result. The candidate the simulation did not select is not estimated. Script 24 runs only against an approvals file listing each approved design report's hash, committed and pushed to GitHub. The confirmatory comparison's lock for death within 180 days is built by script 22b from that comparison's lock. No comparator is run. A sensitivity lock or the 180-day lock that has no design report 7 days after the last of L1 to L5 has one, or that finishes later, is reported as not run. |
| 2026-10-23 | The design stage stopped at about 21:11 on 22 October: the supervisor and its R process were ended from outside during the first four months' L4, after one repetition; the final cohort's L4 had stopped the same morning for lack of memory. Both L4 locks, built on 22 October under the earlier plan and with no design report, are set aside in `abandoned_<date>` folders and not used. Before any further lock was built, and without any outcome, Andrew decided how to finish the design stage in about a day (spec 2026-10-20, section 10): L1, L2, L4, S2 and S3 use main-terms outcome regression only and 100 repetitions without extension; L5 keeps both candidates without extension; L3, S1, S4 and the 180-day lock are unchanged; the first four months are limited to L3 and L1; the run restarts as a scheduled task with up to 14 workers. |

## 5. Files not to open

Started 2026-10-23. Files and folders that no one working on the re-execution's design stage, including automated agents, opens before unblinding, because they hold outcome values or outcome counts. The data steward's scripts (20 and 24) read the outcome sources through their own column lists and are not covered by this list.

| Path | Why | Added |
|---|---|---|
| `rescueCo/logs/pipeline.log` | Its first lines come from the earlier unblinded analysis and include outcome counts (entry of 2026-10-22) | 2026-10-23 |
| `rescueCo/data/outcomes/` | The outcome store, once script 24 builds it | 2026-10-23 |
| `rescueCo/data/Final Data Cut 03.27.2026/6 Month Follow-up/` | The six-month follow-up file | 2026-10-23 |
| Files directly in `rescueCo/results/` (not `reexec_050/`) | Outputs of earlier analyses on the outcomes, such as `binary_outcome_comparison.csv` | 2026-10-23 |
| `rescueCo/reports/` and `rescueCo/README_legacy.md` | Treated as outcome-bearing by the 22 October audit; contents not checked | 2026-10-23 |
| `../rescueCo_analysis`, `../rescueCo_analysis-audit`, `../rescueCo_analysis-severity` | The unblinded analysis and its offshoots (section 1) | 2026-10-23 |
