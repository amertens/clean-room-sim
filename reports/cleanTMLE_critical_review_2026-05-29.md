# Critical review: cleanTMLE manuscript, vignette, and package documentation

Date: 2026-05-29. Reviewer pass over `reports/manuscript_outcome_blind_dq.qmd`
(3,200 lines), `cleanTMLE/vignettes/cleanTMLE-staged-analysis.qmd`,
`cleanTMLE/vignettes/cleanTMLE-functions.qmd`, `cleanTMLE/DESCRIPTION`, and
`cleanTMLE/README.md`. Direct edits already applied are listed at the end under
"Changes already made." Everything else is a recommendation, with a
ready-to-paste prompt in the final section.

The materials are in good shape. The manuscript is carefully argued, the scope
is stated honestly, and the prose already follows most of the project's style
rules. The recommendations below are concentrated, not a wholesale rewrite.

---

## 1. Package design and description

### 1.1 Reconcile the version number the paper documents

The DESCRIPTION declares `Version: 0.1.5`. The manuscript documents the case
study as run on `cleanTMLE 0.1.1` and refers in places to `0.1.0` and `0.2.1`
(line 2826). A software paper must name one released version and use it
consistently, because reviewers and future users will install exactly that
version to reproduce the results. Recommendation: pin the paper to a single
tagged release, archive that tag (for example with a Zenodo DOI), and either
re-run the case study under that version or add one sentence explaining why the
applied example was executed under 0.1.1 and what changed between 0.1.1 and the
documented release. Remove stray references to 0.2.1 from the body except in the
explicit roadmap.

### 1.2 Tighten the DESCRIPTION Description field

The Description field is a single long paragraph that mixes the tested v0.1
scope with the experimental time-to-event surface and the DQ stress test
internals. CRAN and JSS reviewers read this field first. Recommendation: cut it
to roughly six to eight sentences covering what the package does, the tested
scope (binary point treatment, binary outcome, marginal risk difference), the
two methodological contributions (outcome-blind plasmode candidate selection and
the data-quality stress test), and a one-line pointer to experimental
extensions. Move the per-mechanism detail to the README and vignette where it
already appears.

### 1.3 The experimental surface area is a reviewer risk

A large number of exported functions (the `identify_*` DSL, `estimate_ipwrisk`,
`estimate_gcomprisk`, `estimate_aipwrisk`, `estimate_ipwhr`,
`estimate_surv_tmle`, `estimate_lmtp`) are documented as experimental and
outside the tested scope, yet they are exported on equal footing with the tested
functions. Exporting untested functions invites a reviewer to test them and file
issues against behaviour the paper disclaims. Recommendation: pick one of two
postures and apply it uniformly. Either (a) keep them exported but add a
`@section Experimental:` block to each `.Rd` and emit a one-time
`packageStartupMessage` or per-call `rlang::warn(..., .frequency = "once")`
identifying the function as experimental, or (b) move them to internal
(`@keywords internal`, drop from NAMESPACE) until the 0.2 survival release. The
current README disclaimer is good but does not bind the functions themselves.

### 1.4 Rename the `override_clean_room` argument

The README already concedes the name is "historical" and "the check records the
outcome-blind state inside the package." The argument appears 34 times across
`R/`. A guard whose name does not describe what it does is a usability defect in
a package whose entire selling point is auditable governance. Recommendation:
introduce `allow_outcome_access` (or `confirm_unblinded`) as the documented
argument, keep `override_clean_room` as a deprecated alias that warns via
`lifecycle::deprecate_warn()` (or a manual check), and update the vignette and
manuscript walkthroughs. Do this before any 1.0 tag, since renaming after a
stable release is costlier.

### 1.5 Fix the "GO/FLAG/STOP" usage against the project's own rule

The project style guide says to avoid "GO/FLAG/STOP" in running prose and to use
"the staged decision rule" or "the prespecified threshold" instead. The
DESCRIPTION uses the literal label twice and the README twelve times in prose,
and the manuscript uses it in eight prose sentences. The rule is sensible:
GO/FLAG/STOP should appear only when naming the literal returned value (for
example "the rule returns STOP"), not as a generic noun for the mechanism.
Recommendation: sweep all three documents. Keep the labels where a literal
verdict is reported; replace generic prose uses with "the staged decision rule."

### 1.6 Align claimed artefacts with implemented functions

The README dossier list and the manuscript dossier table reference a checkpoint
dashboard. The README is honest that `checkpoint_dashboard()` is "planned" while
`gate_all()` is the current composite. Make sure no table presents a planned
artefact as if it ships. A reviewer who calls a named function and finds it
absent will distrust the rest of the inventory. Recommendation: in any
user-facing inventory, mark planned items explicitly or omit them.

### 1.7 Document the test suite and the reproduction pipeline in the paper

The package carries a real test suite (17 `testthat` files, including
`test-variance-methods.R` and `test-plasmode_dq.R`), which is a strength the
paper does not currently advertise. JSS reviewers look for this. Separately, the
repository root holds several scratch-style scripts
(`run_simulation.R`, `run_simulation_full.R`, `_bootstrap_variance.R`,
`_smoke_bootstrap.R`, `_diag_c.R`) and results land in both `results/` and
`results_new/`. Recommendation: add a short "Reproducibility" subsection stating
the test coverage and providing a single documented entry point (a `Makefile` or
one `reproduce.R`) that maps each script to the figure or table it generates,
and consolidate the results directories. Remove or relocate the underscore-prefixed
diagnostic scripts so the released tree is clean.

### 1.8 Naming and vocabulary consistency

The codebase mixes "clean room", "cleanroom", and "staged analysis"
(`run_clean_tmle`, `cleanroom_workflow`, `summarize_cleanroom_results`,
`tmle_clean_room_wrapper`). Code identifiers can keep their names, but the
documentation should consistently frame the method as the staged-analysis
framework of Muntner et al. (2024) and define "clean room" once as that paper's
manufacturing analogy. The manuscript already does this at first use; the README
and DESCRIPTION should match. A one-paragraph glossary near the top of the README
would resolve it.

---

## 2. Simulation example and outcome-blind workflow clarity

### 2.1 The strongest device is buried: surface the workflow contrast early

The "Traditional pipeline vs outcome-blind staged workflow" table (manuscript
lines 1861-1867) is the clearest demonstration of what the package buys an
analyst: identical fits, different decisions, with the Scenario C STOP blocking a
biased estimate the traditional pipeline would publish. It currently sits deep in
the simulation section. Recommendation: reference it in the introduction or the
walkthrough as the one-line value claim, and consider promoting the table itself
nearer the front of the simulation section so the reader meets the payoff before
the estimator-by-estimator detail.

### 2.2 Add a degradation-gradient figure

The degradation gradient is the conceptual centre of the DQ stress test, and the
text describes it well (the "shallow gradient" versus "cliff" contrast), but no
figure shows it. The simulation figures are boxplots, a coverage bar chart, an
SE-calibration scatter, and a forest plot, none of which plots performance
against threat severity. Recommendation: add one figure with severity on the x
axis (for example MCAR fraction, or outcome-misclassification severity), bias or
RMSE ratio on the y axis, one line per candidate, and a horizontal line at the
locked threshold. This single visual would make the screen's mechanism concrete
and is the most persuasive missing piece.

### 2.3 The case-study DQ result at n_sims = 3 undercuts the headline contribution

The data-quality stress test is the paper's primary methodological novelty, yet
the one real-data demonstration runs it at three replicates, flagged repeatedly
as a smoke test, with a companion `rescueco_dq_full.yml` at 200 replicates
mentioned but not reported. A reviewer will read this as the central method not
being exercised on real data. Recommendation: run the 200-replicate
configuration, report those numbers and the degradation gradient as the
case-study DQ result, and retain `n_sims = 3` only as the fast-reproducibility
config in the bundle. This also lets the four bracketed `[NOTE: ...]` caveats be
deleted (see 3.2).

### 2.4 Resolve the Scenario A/B/C label collision

The main simulation defines Scenario A (good overlap), B (marginal overlap), and
C (unmeasured confounding). The variance sub-study (lines 1689-1734) reuses
"Scenario A/B/C" but redefines C as "very good overlap, minimal PS spread,"
which directly contradicts the main Scenario C. This is a real source of reader
confusion in the most technical part of the paper. Recommendation: rename the
variance-study regimes (for example "overlap regime 1/2/3" or "V-good / good /
marginal") so the two label sets do not collide.

### 2.5 Elevate the "STOP is not detected bias" guardrail

The distinction that a STOP is a failure of prespecified thresholds under the
chosen DGP, not detection of real-data bias (lines 1636-1648), is the single
most important interpretive caveat for the whole method. It is currently a
mid-section paragraph. Recommendation: state it once, prominently, at the start
of the simulation interpretation (or as a set-off definition), and cross-reference
rather than restate it elsewhere.

### 2.6 Report Monte Carlo error uniformly

The main simulation reports 200 replicates with the Monte Carlo standard error
of coverage, which is good practice. The variance sub-study reports 100
replicates without an explicit MC SE, and the gate operating-characteristics
section should state its replicate count and uncertainty too. Recommendation:
report MC SE for every coverage and bias figure the paper quotes, including the
sub-studies, so the reader can separate signal from noise consistently.

---

## 3. Writing: AI-sounding language, motivation, and usage clarity

The manuscript prose is already clean. After the direct edits below it contains
no em-dashes, no banned filler vocabulary, and only one defensible "not ... but"
construction. The remaining issues are concentrated.

### 3.1 Consolidate the repeated disclaimer

The formula "cleanTMLE does not replace target-trial design, fit-for-purpose
review, validation, governance, or post-outcome sensitivity analysis" recurs in
the abstract, the introduction, the framework section, the DQ section, and the
discussion. Each instance is individually reasonable, but the repetition reads as
defensive and is one of the clearest machine-generated tells in the document.
Recommendation: state the boundary once, authoritatively (the dossier or scope
subsection is the natural home), and replace later instances with a short
cross-reference.

### 3.2 Remove the bracketed `[NOTE: ...]` editorial markers

The manuscript contains inline `[NOTE: ...]` drafting notes (around lines 2383,
2396, 2434, and 2568) that explain the n_sims = 3 caveat. These are
work-in-progress artifacts and must not appear in a submitted manuscript.
Recommendation: once the 200-replicate case-study DQ run is in (2.3), delete
them; if any caveat must survive, fold it into ordinary prose or a footnote.

### 3.3 Shorten over-long captions

Several figure and table captions carry three to six sentences of methodological
detail (the case-study effect-estimate table caption at line 2503 is the worst
case). JSS captions should be one to two sentences; substantive detail belongs in
the body. Recommendation: move the IPCW fallback explanation and the
`method_used` note out of the caption and into the surrounding paragraph.

### 3.4 Sweep em-dashes from the README and vignette

The manuscript is now clean, but the README still contains roughly 28 inline
em-dashes (` --- `) and the staged-analysis vignette about 15, in violation of
the project's no-em-dash rule. These render as em-dashes under Pandoc.
Recommendation: replace each with a comma, colon, parentheses, or a restructured
sentence. This is mechanical but should be reviewed per instance because many sit
in definitional list items ("**Estimand-first design** --- declare ...").

### 3.5 Lead with the usage payoff for a JSS audience

JSS values "what the software does and how to use it" over methodological
advocacy. The introduction is strong on governance motivation but the concrete
"why install this" answer (the workflow-contrast result and the four-step TMLE
that mirrors the stages) arrives late. Recommendation: add one or two sentences
early that state the operational payoff plainly: the package turns prespecified
targeted-learning decisions, outcome-blind estimator selection, and data-quality
stress testing into a single reproducible record, and it refuses to run the
primary analysis until the prespecified checkpoints are recorded.

### 3.6 Minor: stabilise the abstract's computed values

The abstract embeds several `r round(...)` inline expressions, so it cannot be
read without rendering and will shift if the cached results change.
Recommendation: once the reported numbers are final, consider hard-coding the
headline coverage figures in the abstract for stability, leaving the dynamic
values in the body.

---

## Changes already made (direct edits)

Applied to `reports/manuscript_outcome_blind_dq.qmd`:

- Replaced all seven ` -- ` en-dash constructions in prose with commas, a
  semicolon, or parentheses (abstract line 113; variance discussion lines
  1700-1701, 1707-1708, 1725, 1763).
- Removed the editorializing "Notably," at line 1721.

Applied to `cleanTMLE/README.md`:

- Corrected "Muntner et al. (2020)" to "Muntner et al. (2024)" (line 122), which
  was inconsistent with every other citation of that work in the manuscript and
  README.

---

## Ready-to-paste prompt for the remaining work

The following prompt covers the bulk and judgment-dependent items. It is written
to be pasted into a coding session with the repository open.

```
Work in the clean-room-sim repository on the cleanTMLE software paper. Follow
the project's writing rules: flowing academic prose, no em-dashes, no "not X but
Y", no editorializing adverbs, no coined terms, and do not use "GO/FLAG/STOP" in
running prose (use "the staged decision rule" or "the prespecified threshold";
keep the literal labels only when reporting a returned verdict such as "returns
STOP"). Make the following changes and show me a diff before saving large edits.

1. Version reconciliation. Decide on one released version for the paper to
   document (default: the current DESCRIPTION version). In
   reports/manuscript_outcome_blind_dq.qmd, make every version reference
   consistent with it, except the explicit 0.2 roadmap. If the case study was
   run on 0.1.1, add one sentence explaining the gap; otherwise update the
   references.

2. DESCRIPTION. Cut the Description field to 6-8 sentences: what the package
   does, the tested v0.1 scope, the two methodological contributions
   (outcome-blind plasmode candidate selection and the data-quality stress
   test), and a one-line pointer to experimental extensions. Remove
   "GO/FLAG/STOP" from the field.

3. Experimental functions. For each experimental function (the identify_* DSL
   and the time-to-event estimators), add a roxygen "@section Experimental:"
   block stating it is outside the tested v0.1 scope, and add a once-per-session
   rlang::warn(.frequency = "once") at the top of each. Re-run roxygen.

4. Rename the outcome guard. Introduce allow_outcome_access as the documented
   argument across R/, keeping override_clean_room as a deprecated alias that
   warns once when used. Update the vignettes, README, and manuscript
   walkthroughs to the new name.

5. Em-dash sweep. In cleanTMLE/README.md and
   cleanTMLE/vignettes/cleanTMLE-staged-analysis.qmd and
   cleanTMLE/vignettes/cleanTMLE-functions.qmd, replace every inline " --- " and
   " -- " used as a dash with a comma, colon, parentheses, or a restructured
   sentence. Review each definitional list item individually.

6. GO/FLAG/STOP prose sweep. In the manuscript, README, and DESCRIPTION, replace
   generic prose uses of "GO / FLAG / STOP" with "the staged decision rule",
   keeping the literal labels only where a returned verdict is reported.

7. Remove the [NOTE: ...] markers. In the manuscript, delete the bracketed
   editorial notes around lines 2383, 2396, 2434, and 2568, folding any
   surviving caveat into prose or a footnote.

8. Consolidate the disclaimer. Find the recurring "does not replace target-trial
   design / fit-for-purpose review / validation / governance / post-outcome
   sensitivity" statement (abstract, intro, framework, DQ, discussion). Keep one
   authoritative instance and replace the others with a short cross-reference.

9. Caption length. Shorten figure and table captions to one or two sentences,
   moving methodological detail (especially the case-study effect-estimate table
   caption near line 2503) into the body.

10. Scenario label collision. In the variance sub-study (manuscript lines
    1689-1734), rename the overlap regimes so they no longer reuse the labels A,
    B, C, which collide with the main simulation's Scenario C (unmeasured
    confounding). Use "overlap regime 1/2/3" or similar and update the table and
    surrounding prose.

11. Case-study DQ at full replication. Run cleanTMLE's
    inst/extdata/rescueco_dq_full.yml at n_sims = 200, save the results, and
    update the case-study DQ subsection to report those numbers. Retain the
    n_sims = 3 config only as a fast-reproducibility example.

12. Degradation-gradient figure. Add one figure to the simulation results
    section plotting threat severity (x) against bias or RMSE ratio (y), one line
    per candidate, with a horizontal line at the locked threshold. Caption it in
    one sentence.

13. Reproducibility. Add a short Reproducibility subsection to the manuscript
    stating the testthat coverage and providing a single entry point that maps
    each script to the figure or table it produces. Consolidate results/ and
    results_new/ and relocate the underscore-prefixed diagnostic scripts at the
    repo root out of the released tree.

After each edit, re-render the affected document and confirm it knits without
error. Report which figures or tables changed.
```
