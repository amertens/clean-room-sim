# clean-room-sim

Companion repository for the manuscript

> **cleanTMLE: An R Package for Staged, Outcome-Blind Targeted Learning
> with Plasmode Data-Quality Stress Testing**

The repository holds the manuscript and talk, a tutorial, the Rescue.Co
case study, the simulation study behind the manuscript's results, the
revision protocol under which all of it is being revised, and a git
submodule with the [`cleanTMLE`](cleanTMLE/) R package itself.

## What is here

| Path | Purpose |
|------|---------|
| [`cleanTMLE/`](cleanTMLE/) | The R package (version 0.3.0), a submodule pointing at `https://github.com/amertens/cleanTMLE`. Its workflow surface is sixteen verbs that run in stage order against one analysis-lock object. |
| [`reports/`](reports/) | The manuscript (`manuscript_outcome_blind_dq.qmd`), the talk (`cleanTMLE_presentation.qmd`), the shared bibliography, and two companion analyses. See [`reports/README.md`](reports/README.md). |
| [`tutorials/`](tutorials/) | `clean_room_targeted_learning_tutorial.qmd`: from the clean-room design to the causal roadmap, worked on a simulated cohort with a known truth. |
| [`rescueCo/`](rescueCo/) | The Rescue.Co Kenya Trauma Registry case study: scripts, configuration, results, and the reconciliation against the parallel unblinded analysis. Patient-level data are access-restricted and not tracked. |
| [`revision/`](revision/) | The September 2026 revision protocol: the Stage 0 audit, findings, decisions, the claim-level source map, and the tooling (`style_check.R`, `build_source_map.R`, `render_all.R`). `CHANGELOG.md` is the running log. |
| [`sandbox/`](sandbox/) | Focused side studies the manuscript cites: candidate divergence, estimator validation against `tmle::tmle()`, plasmode fidelity, and others. |
| [`docs/revision/`](docs/revision/) | Design notes and literature memos written during the revision. |
| [`plasmode_selection_paper/`](plasmode_selection_paper/) | A companion pilot simulation on plasmode-based candidate selection. |
| [`archive/`](archive/) | Retired material, kept verbatim for provenance: the 0.1.x simulation drivers, the May to July 2026 result files, and the short-lived cleanroomGov package. See [`archive/README.md`](archive/README.md). |
| [`external docs/`](external%20docs/) | Reference PDFs (targeted-learning and real-world-evidence roadmap papers, SPIFD2, Muntner 2024, FIORD 2026). Local copies for convenience; not redistributed. |

## What the project does

The package and the studies together make two methodological
contributions, both computed in the blind phase of a staged analysis in
the sense of Muntner et al. (2024).

1. **An estimability answer before outcome access.** `cleanTMLE` grades
   the overlap on the fitted propensity score, reads the same fit per
   estimand (the effect among the treated and the overlap-weighted
   effect can stay estimable where the average treatment effect is
   not), maps feasibility with an outcome-blind generate-treatment
   plasmode simulation, and makes any switch to a declared fallback a
   logged design decision through a pre-registered estimand ladder.
2. **A plasmode data-quality stress test.** `stress_test()` extends the
   standard outcome-blind plasmode loop with five prespecified threat
   families (covariate missingness under MCAR, MAR, and MNAR; treatment
   misclassification; outcome misclassification; near-positivity;
   unmeasured confounding) and reports, for every candidate estimator,
   a degradation gradient of bias, RMSE, and coverage against
   thresholds declared and fingerprinted on the lock.

Around both, the analysis lock fingerprints the specification and the
design data, keeps the primary outcome in a separate store that only
the estimation stage joins back, and records every declaration, switch,
and approval in a design log that `export_design_log()` releases in the
decision-log format of the staged-analysis literature.

## Status of the simulation results

The manuscript's simulation numbers were produced by the 0.1.x drivers
now in `archive/scripts_0.1x/`, and their result files live in
`archive/results_new_0.1x/`. Under the revision protocol (decision D4)
the drivers are being rebuilt on the sixteen-verb surface and every
scenario rerun under cleanTMLE 0.3.0 into a fresh `results_new/`
directory; the manuscript reads from that path and will not render its
simulation sections until the reruns exist. `RUN_MANIFEST.md` lists the
runs and their approximate cost.

## Reproducing what can be reproduced today

Install the package from the submodule:

```r
devtools::install("cleanTMLE")
```

Render the tutorial (the first render runs its stress test, about ten
minutes; the chunk is cached for later renders):

```bash
quarto render tutorials/clean_room_targeted_learning_tutorial.qmd
```

Render the talk:

```bash
quarto render reports/cleanTMLE_presentation.qmd
```

The case-study scripts are `rescueCo/scripts/10_multiarm_build.R`
through `14_multiarm_report.R`; they require the registry extract and
are documented in `rescueCo/README.md`. Before committing a document,
run the style gate:

```bash
Rscript revision/style_check.R --strict <files>
```

## Citation

If you use `cleanTMLE` or the staged, outcome-blind workflow described
here, please cite:

- Mertens, A. *cleanTMLE: An R Package for Staged, Outcome-Blind
  Targeted Learning with Plasmode Data-Quality Stress Testing.* (in
  preparation).
- Muntner P., Hernandez R. K., Kent S. T., et al. *Staging and
  Clean Room: Constructs Designed to Facilitate Transparency and
  Reduce Bias in Comparative Analyses of Real-World Data.*
  Pharmacoepidemiology and Drug Safety, 33(3):e5770, 2024.
- Gatto N. M., Vititoe S. E., Rubinstein E., Reynolds R. F.,
  Campbell U. B. *A Structured Process to Identify Fit-for-Purpose
  Study Design and Data to Generate Valid and Transparent
  Real-World Evidence for Regulatory Uses.* Clinical Pharmacology
  and Therapeutics, 113(6):1235-1239, 2023.

## Licence

MIT. See [LICENSE](LICENSE).
