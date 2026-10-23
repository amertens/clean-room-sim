# Departures from the protocol in the unblinded Rescue.Co analysis

Approved by Andrew on 18 October 2026. The table is `deviations_inventory.csv`, one row per departure.

## What this is

The unblinded analysis (`rescueCo_analysis`) departed from its protocol (version 1, 9 May 2023) 34 times. The inventory starts from section 4 of `reports/cleantmle_rescueco_review_2026-10-18.md`, which listed 31. Each row was checked against the protocol text, Supplementary Table S1 (b), the departures document and the git history of `rescueCo_analysis`. The check found three more departures (rows 32 to 34). It also found two identification problems (I1 and I2) that the protocol shares rather than departs from, but which a design stage should raise. Effect estimates, P values and outcome rates are deliberately left out of the table.

Every item was found after unblinding. The table says whether each was *detectable before outcome access*, never that it was detected.

## How the rows are coded

The date is that of the first change, found with `git log -S` on the relevant code. Where Supplementary Table S1 (b) gives a date, it is shown beside the first-change date. S1 (b) usually dates a change's current form. For rows 1, 2, 3, 13, 14, 15, 18 and 20 its date is weeks or months later than the first change.

The code says what prompted the change:

- **D**: design data only.
- **S**: how the outcome is defined or observed, or arm-blind checks a data steward could run.
- **O**: estimates were available when the change was made, or are cited as its reason.
- **X**: outside information.

O records a fact about the workflow, not a judgment of motive. Effect estimates on every outcome existed from 26 February 2026 (54f0ef7, `analysis_outputs/matched_effects_summary.csv`). So O applies to every change after that date. The column `estimates_at_first_change` separates the two senses:

- *available*: estimates existed;
- *cited*: estimates appear in the commit, review or stated reason for the change.

## What the inventory shows

- **O almost everywhere.** 31 of the 34 departures were made while estimates existed. The exceptions are rows 10, 31 and 34.
  - In 14 of the 31, estimates appear in the record of the change.
  - Several stated reasons changed after the change was made. Examples are row 14, where the header of script 78 was rewritten 3.5 hours after the switch, and row 22, where the reason changed three times in a week.
  - In an unblinded workflow, none of these changes can be shown to be independent of the estimates, whatever its reason. A workflow that makes its changes before any estimate exists can show that for every one of them.
- **D or S departures.** 22 departures carry D or S.
  - 18 could have been raised by an outcome-free output of a design stage.
  - 4 only in part, because outside information was also needed:
    - row 5: the gap is visible in design data, but only Rescue.co's dispatch records could settle it;
    - row 7: overlap by recorded mode shows the problem, but whether a walk-in could have used an ambulance is clinical knowledge;
    - row 12: a covariate's provenance can be fixed in advance, but which covariates matter came partly from co-authors and pooled outcome associations;
    - row 28: the bias analysis for deaths before arrival can be locked, but its plausible range needs outside data.
- **The checks behind those D or S rows.**
  - Only three draw on cleanTMLE 0.4.1 as it stands:
    - the locked outcome (row 34);
    - the simulation verdict (row 2, in part);
    - the trim band (row 21, in part).
  - The rest need the proposed items:
    - P1: time zero, eligibility timing and a variable dictionary (rows 3, 6, 10, 11, 12, 13, 14, 15, 16, 19, 30);
    - P2: the observation structure (rows 14, 17);
    - P3: several comparisons, restrictions and who is unsupported (rows 1, 2, 3, 5, 7, 21, 27);
    - P5: selection after time zero (row 28);
    - P7: comparators and libraries as declared candidates (rows 2, 24).
  - Rows 6, 18, 19 and 30 need arm-blind checks run by a data steward script at build time, not the package.
- **Not caught.** Row 31, calibration of the control arm against the population before Rescue.Co, needs control-arm outcomes. A design stage could only have asked for the pre-period data.

## Corrections to the review

- **Row 1, the comparator.** The change was first made on 1 April 2026 in `clean_room` (e5cce28), labelled "pre-specified". The 5 May date (2029052) is when the main pipeline followed. The review says the estimate moved "from a benefit to null". That holds for the GOSE outcomes only, and 2029052 changed scoring, caliper and covariates at the same time.
- **Row 2.** `clean_room`'s outcome-blind simulation selected propensity-score-matched regression on 12 March and 4 April 2026. The report of 4 April called TMLE primary. The 2024 checks passed on the protocol comparator: 80% of Rescue.Co patients matched through June 2024, 97% through October. The retention problem arises only after the comparator was narrowed (row 1).
- **Row 8.** The ordinal GOSE dates from 5 May, not 27 May.
- **Row 28.** Script 61's grid adds deaths to either arm. Only its tipping point is one-sided.
- **Row 29.** The bias curve is coded O, not D. It is a function of the estimates, and it was dropped without mention on 31 August.
- **Rows 2, 3, 5, 6, 7, 14, 18, 24, 27 and 30** gain O.

## Approval

Andrew approved the inventory on 18 October 2026, including the split of O into estimates available and estimates cited, and keeping I1 and I2 as rows. Three points that rested on records outside both repositories were dropped at his request, because they cannot be backed up: the source of the prior-care exclusion (row 3), the amendment cited for the ordinal GOSE (row 8), and how cross-validated risks were pooled to choose the boosted-tree settings (row 26).

The column `blinded_workflow_decision` stays empty until Phase 5 fills it with what the blinded workflow decided and when.
