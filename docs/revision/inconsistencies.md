# Stage 1: Inconsistency audit (report only, no edits)

Date: 2026-06-08. Source of truth for edits is the `.qmd` and the R sources, not
the rendered `.docx`/`.html`. Every claim below is cited to file and line. No
numbers, methods, or files were changed in producing this report.

Status legend:
- **OPEN** = a real inconsistency that needs a decision or edit.
- **RESOLVED** = the prompt's premise does not hold against current source; at most a wording tightening remains.
- **VERIFIED** = a concrete bug confirmed in the source, ready to fix in a later gated stage.

---

## 0. Source-of-truth: is the `.docx` built from the `.qmd`?  -- RESOLVED (reassuring case)

This was the load-bearing question. Answer: **the `.docx` is a render target of the
`.qmd`, not a separately maintained file.**

- `reports/manuscript_outcome_blind_dq.qmd:7-20` declares both `format: html` and
  `format: docx` in the YAML front matter. Running `quarto render
  manuscript_outcome_blind_dq.qmd` produces both the `.html` and the `.docx`.
- No external render script or Makefile references the file; it is rendered
  directly by Quarto.
- Modification times: `.qmd` = 2026-06-08 04:28 (today); `.docx` and `.html` =
  2026-06-04 06:56 (four days older). The committed `.docx`/`.html` are **stale
  renders**, four days behind the `.qmd`.

Consequences:
1. Editing `reports/manuscript_outcome_blind_dq.qmd` is the correct target. All
   line-referenced edits in Stages 3 and 4 land in the file the `.docx` is built
   from.
2. **A re-render is a required, explicit step**, not an assumption. After Stage 3
   and Stage 4 edits, the manuscript must be re-rendered (`quarto render
   reports/manuscript_outcome_blind_dq.qmd`) so the changes reach the `.docx` a
   reviewer reads. This will be a named acceptance step, gated, and run only with
   `results_new/` present so figures resolve (see item 6c) and no numbers change.
3. The apparent `.qmd`-vs-`.docx` divergence you flagged (broken 6.3.1 header,
   dangling cross-refs "in the .docx but not the .qmd") was an **audit miss on my
   first pass, not a real divergence**: those artifacts are in the `.qmd` too
   (`:1091-1092`, `:3379`) and in the stale `.html` render
   (`...html:3925`, `:6565`). The `.docx` is a faithful but old render of the same
   `.qmd`. Section numbering in the `.docx` may differ cosmetically (the
   cross-ref table sits in a late appendix that `.docx` may number as Section 13
   or 14), but the content source is the one `.qmd`. There is no separate
   document to keep in sync.

---

## 1. Scope contradiction: time-to-event / causalRisk-style surface  -- OPEN

**The prompt's premise about a hidden causalRisk dependency is false; the scope conflict itself is real.**

Dependency status (`cleanTMLE/DESCRIPTION`, lines 38-61): no `causalRisk` in
Depends, Imports, Suggests, or Remotes. `Suggests` lists `survtmle` and `lmtp`.
The time-to-event surface is an **original reimplementation inside cleanTMLE**,
not a wrapper and not a dependency.

Exported and implemented (all in `NAMESPACE`, all with tests under `tests/testthat/`):

| Function | R file | approx lines |
|---|---|---|
| `estimate_ipwrisk` | `R/estimate_ipwrisk.R` | 364 |
| `estimate_gcomprisk` | `R/estimate_gcomprisk.R` | 234 |
| `estimate_aipwrisk` | `R/estimate_aipwrisk.R` | 222 |
| `estimate_ipwhr` | `R/estimate_ipwhr.R` | 188 |
| `specify_models` | `R/specify_models.R` | ~47 |
| `identify_outcome/censoring/treatment/competing_risk/interval/missing/subject` | `R/identify.R` | 32-192 |
| `make_table1`, `make_table2`, `make_wt_summary_table`, `extreme_weights` | `R/tables.R` | ~325 total |

Each estimator file carries a header warning that it is "outside the tested v0.1
scope of cleanTMLE (binary point exposure / binary outcome only)." Tests:
`test-estimate_ipwrisk.R`, `test-estimate_ipwhr.R`, `test-gcomp_aipw.R`,
`test-specify_models.R`, `test-tables.R`.

Manuscript says these estimands are out of scope:
- `reports/manuscript_outcome_blind_dq.qmd:321-329` "Hazard ratios, time-varying treatments, competing risks, and multinomial exposures are out of scope for the v0.1 release ... Analysts who need those estimands can compose cleanTMLE's prespecification and audit functions with `survtmle` or `lmtp`."
- Also `:1551`, `:2783-2792`, `:3177-3193`.

Vignette nonetheless documents the full surface as usable:
- `cleanTMLE/vignettes/cleanTMLE-functions.qmd:451-527` (Model Specification DSL, Time-to-Event Estimation, Tables and Weight Diagnostics).

**The contradiction a reviewer will see:** the manuscript directs users away from
estimands that the package exports and the function vignette demonstrates.

**This is a deletion decision, not a bug fix.** Removing the surface deletes
roughly 1000 tested lines and is irreversible. Two clean resolutions:
- (a) Quarantine: keep the code, mark it `@keywords internal` / move its vignette
  coverage into an explicit "experimental, outside tested v0.1 scope" appendix,
  and add one manuscript sentence acknowledging an experimental survival surface
  exists but is unvalidated. Lower risk, preserves work.
- (b) Remove from exports and the vignette, defer to survtmle/lmtp as the
  manuscript already says. Matches the prompt's default, larger blast radius.

Decision deferred to you (Stage 2 `surface_plan.md`).

---

## 2. MAR / MNAR missingness  -- RESOLVED (premise does not hold)

The prompt expects a live contradiction (manuscript says MAR/MNAR unsupported;
vignettes pass `_mar`/`_mnar`). It is not live.

The DQ stress test fully implements both:
- `R/plasmode_dq.R` `.degrade_missingness_mar()` lines 347-366 (treatment-dependent, `treatment_OR`).
- `R/plasmode_dq.R` `.degrade_missingness_mnar()` lines 385-405 (value-dependent, `strength`).
- `default_dq_scenarios()` lines 39-75 exposes `covariate_missingness_mar` and `covariate_missingness_mnar`.

The manuscript's "not yet supported" refers to native missingness *handling*
(IPW-of-missingness models), not stress-test *degradation*, and it already draws
that distinction:
- `manuscript_outcome_blind_dq.qmd:3114-3119` "MAR and MNAR mechanisms require an explicit missingness model and are not yet supported. The principled alternative is inverse-probability-of-missingness weighting ...".
- `:325` "the data-quality stress test now degrades covariates under MAR and MNAR to test the default median-imputation handling against them."
- `:3233-3237` "the data-quality stress test now degrades covariates under MAR and MNAR, but the only handling it exercises is median imputation".

Both vignettes pass `_mar`/`_mnar` correctly (`cleanTMLE-functions.qmd:332-333`,
`cleanTMLE-staged-analysis.qmd:865-870`).

**No edit required.** Optional Stage-3 wording: make `:3114-3119` say "native
handling models" explicitly so a fast reader does not misread it as "the stress
test cannot degrade under MAR/MNAR."

---

## 3. Variance stage  -- MOSTLY RESOLVED (one stale vignette sentence)  -- VERIFIED

Both functions are implemented, non-stub:
- `R/variance_methods.R` `select_variance_method()` lines 216-291.
- `R/variance_methods.R` `bootstrap_rd_variance()` lines 160-186, helper `.rd_point_estimate()` 63-112.

Manuscript describes them accurately at `:137-140`, `:1067-1070`, `:3094-3099`.

The only conflict is a vignette sentence:
- `cleanTMLE/vignettes/cleanTMLE-functions.qmd:676` "Variance-method selection on the locked point estimator (FIORD Stage 2 proper) is still pending; this rule currently selects the point-estimator candidate."

Both statements are true but read as contradictory: the standalone
`select_variance_method()` exists, yet it is **not wired into** the
`fiord_two_stage` rule of `select_tmle_candidate()`. The fix is a one-line vignette
clarification (state that the standalone selector exists and that what is pending
is its integration into the `fiord_two_stage` rule), not a code change.

---

## 4. Version drift  -- OPEN

| Version | Where | Line |
|---|---|---|
| 0.1.5 | DESCRIPTION; "This paper documents version 0.1.5" | `DESCRIPTION:Version`; `manuscript:330` |
| 0.1.1 | case study estimates produced under 0.1.1; "All five stages were run on cleanTMLE 0.1.1" | `manuscript:331`, `:2514`, `:2571`, `:2830`, `:2885` |
| 0.1.0 | mentioned in passing | `manuscript:2907` |
| 0.2 | "The 0.2 release adds eight axes to the plasmode candidate grid" | `cleanTMLE-functions.qmd:654`; `manuscript:2903` |
| 0.2.1 | "cleanTMLE 0.2.1 adds a `q0_library` argument" | `manuscript:3128` |
| 0.3 | "on the 0.3 roadmap" | `cleanTMLE-functions.qmd:721` |

Recommendation: canonical = **0.1.5** (matches DESCRIPTION and the paper's own
claim at `:330`). Keep the honest note that the case study ran on 0.1.1 with a
backward-compatibility statement (already at `:331-332`). Move 0.2 / 0.2.1 / 0.3
material into a single clearly-labelled roadmap note so the paper does not appear
to document unreleased versions. Decision: confirm 0.1.5 as canonical.

---

## 5. G-value  -- RESOLVED / VERIFIED (legitimate term, attribution can be sharpened)

"G-value" is a published term, not a coinage by analogy to the E-value. The cited
source defines it verbatim (`Evaluating and improving RWE with TL.pdf`, body text):

> "The recently proposed G-value calls attention to the gap size that would be
> needed to negate the finding from the current study ... For a 95% CI
> G-value = min(|theta_hat - 1.96 sigma - null|, |theta_hat + 1.96 sigma - null|)
> ... [12]. The G-value takes both bias and variance into account ..."

Two points for the record:
- Gruber calls it "recently proposed" and attributes it to reference **[12]**, so
  the originator is [12], not Gruber. The manuscript's "(Gruber 2023 RWE p.5)"
  (`manuscript:1730`) cites where it is *discussed*, which is acceptable, but
  citing the primary source [12] alongside would be more precise.
- The package implementation matches the source formula:
  `R/plasmode_extensions.R` `compute_G_value()` lines 310-332 returns
  `min(|ci_lower - null|, |ci_upper - null|)`, which equals the source's
  min-distance-from-CI-bound-to-null. Correct.

**No code change.** Optional Stage-3 citation tightening only.

---

## Additional verified manuscript bugs (these are real; queued for Stage 3)

The first Explore pass missed these in the `.qmd` because they were searched as
rendered `?@` strings; checking the `.qmd` source and the rendered `.html`
confirms all three are live.

### 6a. Broken Section 6.3.1 header  -- VERIFIED
`manuscript_outcome_blind_dq.qmd:1091-1092`:
```
### When STOP motivates estimand revision, not just estimator
swapping
```
The word "swapping" wrapped onto line 1092, so markdown renders the heading as
"...not just estimator" and drops "swapping" into body text. Confirmed in the
render: `manuscript_outcome_blind_dq.html:3925`. Fix: join the two lines.

### 6b. Duplicate "Table 1"  -- VERIFIED
- `manuscript_outcome_blind_dq.qmd:1035-1037` hard-codes
  `caption = "Table 1: Mapping data-quality concerns to plasmode-parameter inputs ..."`
  in a bare `knitr::kable()` with no Quarto label.
- `manuscript_outcome_blind_dq.qmd:542-543` uses `#| label: tbl-functions`, which
  Quarto auto-numbers as "Table 1".
Result: two "Table 1"s in the render. Fix: give the mapping table a real
`#| label: tbl-dq-mapping` and a caption without the literal "Table 1:" prefix,
then reference it as `@tbl-dq-mapping`.

### 6c. Unresolved cross-references  -- VERIFIED, and partly entangled with results
Rendered render shows them dangling: `manuscript_outcome_blind_dq.html:6565`
(`?@tbl-mc-metrics`, `?@fig-mc-boxplots`, `?@fig-coverage-bar`,
`?@fig-se-calibration`, `?@fig-mc-forest`, `?@fig-dq-degradation`). They are all
referenced from the run-manifest table at `manuscript_outcome_blind_dq.qmd:3379`.

Two distinct causes:
- `tbl-mc-metrics` is **never defined**. The Monte Carlo metrics chunk is
  `mc-table` at `:1714` and carries no `tbl-` label. So `@tbl-mc-metrics` dangles
  unconditionally. Fix: add `#| label: tbl-mc-metrics` + a `tbl-cap` to that
  chunk, or rename the reference.
- The five `fig-*` chunks are defined (`:1938`, `:1957`, `:1977`, `:1991`,
  `:2117`) but gated on `eval = has_mc` / `eval = has_dq`. `has_mc` is set by
  whether `results_new/simulation_results.rds` exists at render
  (`:1700-1708`); `has_dq` similarly at `:2041`. The committed render had these
  FALSE (the `results_new/interim_*` files were deleted/modified per `git
  status`), so the figures were not produced and the refs dangled.

**Guardrail flag:** fixing 6c fully means re-rendering against present results,
which touches the "reported numbers are sacred" rule. The label fix for
`tbl-mc-metrics` is safe and number-neutral. The `fig-*` refs will resolve on
their own when the manuscript is rendered with `results_new/` populated; they do
not require any prose or number change. I will not regenerate or alter any
result. Recommend: in Stage 3, fix the missing `tbl-mc-metrics` label only, and
note in `structure_notes.md` that a clean render requires `results_new/` present.

---

## Jargon and AI-language inventory (for Stage 4; reported, not yet edited)

- `operationalis*` ~12: `manuscript:267, 440, 454-458 (x3), 463, 748, 995, 1118, 1153, 3506`. British spelling throughout.
- "quantitative bridge" x2: `manuscript:234, 262`. Also a section title `:1198` "A quantitative supplement to fit-for-purpose review" (different phrase, keep or rename consciously).
- "degradation gradient" 8+: section header `:1047`; also `:349, 1050-1052, 2076, 2086, 2117 (fig cap), 2289`. Decision: justify as a distinct object or rename to "sensitivity curve over threat severity" / "estimator-fragility curve".
- "paired estimator-and-estimand object" x2: `manuscript:1081-1082, 2638`.
- "robust"/"robustness" ~18: legitimate "doubly robust" at `:1653, 2847`; the rest are general-virtue uses (`:217, 347-349, 870-871, 1055-1057, 1132, 1157, 1202, 1675, 2331, 2344`, plus several truncation/candidate labels like `:1679, 2222, 2294`). Note: candidate *labels* that contain "robust" may be data-bound strings; check before editing those so as not to change a factor level that a result table prints.
- Selection-rule name used three ways: "minimax" (`:2216, 2233, 2324, 2685, 2698`), "min_max_rmse" (`:151, 518, 1268, 2236, 2310, 2590`), "worst-case" (`:1269, 2100, 2326, 2331`). Standardize on one. Note `min_max_rmse` is also the literal `rule =` argument value, so the code-facing string must stay; standardize only the prose name.
- Em-dashes: none found in the manuscript or either vignette prose.

---

## Coupling notes for later stages (from your feedback)

These are recorded here so the gates do not treat coupled work as independent.

- **TTE quarantine must reach the manuscript text, not only the vignette.** The
  contradiction lives in the scope sentences (`manuscript:321-329` around Sec 8.3,
  and `:3177-3193` around Sec 12.5), which say those estimands are out of scope
  while the package exports them. If we keep the estimators (quarantine path),
  Stage 3 must add an explicit sentence to those scope passages stating that an
  experimental, out-of-tested-scope TTE surface ships and is not validated.
  Otherwise a reviewer reads "out of scope" and then finds `estimate_ipwrisk`
  exported. This is named as a required Stage 3 edit, contingent on the Stage 2
  TTE decision.

- **Lock/audit surface reduction (Stage 2) is coupled to manuscript positioning
  (Stage 3).** How aggressively we deprecate the ~35 lock/audit functions must
  follow from how the paper positions the lock. Coherent pairings: keep the
  surface and add the "records prespecification, does not provide integrity or
  access control" sentence; or deprecate heavily and stop leaning on the lock as
  governance in the prose. The incoherent case is deprecating the code while the
  paper still presents the lock as governance. Therefore the surface decision and
  the framing sentence will be decided **together at the `surface_plan.md`
  gate**, not sequentially.

- **Stage 5 novelty framing depends on the separating-grid task.** Whether
  Section 9.4 should be reframed as a generic property or an existence
  demonstration is an empirical question answered by the separate
  separating-grid study (`sandbox/separating_grid/`). If that study runs first
  and yields a divergence at anchored severities, Stage 5 folds it in; if this
  revision pass runs first, Stage 5 explicitly defers the empirical question and
  does not pre-judge the outcome. The two tasks must not emit conflicting
  recommendations.

---

## Decisions needed from you before Stage 2

1. **TTE surface (item 1):** quarantine-as-experimental (recommended, low risk) or remove from exports + vignette (prompt default, deletes ~1000 tested lines)?
2. **Version (item 4):** confirm 0.1.5 as the canonical paper version; move 0.2/0.2.1/0.3 to a roadmap note.
3. **G-value (item 5):** add the primary source [12] alongside the Gruber citation, or leave as is?
4. **Cross-refs (item 6c):** agree to fix only the safe `tbl-mc-metrics` label now and leave the `fig-*` refs to resolve at render time, with no touching of results?

Stage 2 will not begin until these are answered.
