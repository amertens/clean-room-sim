# Design-log notes (multi-arm rebuild)

Entries recorded by the design team before outcome estimation results
were read. Times are local (America/Los_Angeles).

## 2026-09-13 18:05 — Support-map reading rule

Decision: the support maps from `simulate_support()` (200
generate-treatment replicates per contrast; surface grid = synthetic
confounding strength {0, 1, 2} x effect modification {off, on}, linear
complexity; bias tolerance 0.01, coverage floor 0.90) are read over the
moderate-confounding region (synthetic confounding strength <= 1) for
the support/feasibility question. The whole-grid pass summary is
reported alongside as a conservative bound, not as the feasibility
rule.

Rationale: the harshest rung (strength 2) stresses the synthetic
outcome's dependence on the separating direction to a degree that
behaves like strong unmeasured confounding for every weighting-based
estimand at once, including on the PRIMARY contrast whose static
support verdict is PASS with maximum ATE weight 8. A rule that requires
a whole-grid pass therefore conflates support failure with
confounding-strength stress and would keep only the matched ATT on
every contrast. The diagnostic signature used instead: failure at the
mildest surface (before any synthetic confounding) indicates a support
failure (observed for the ATE on C1, coverage 0.87 at strength 0);
failure that appears only as the stress strengthens and hits all
unmatched estimands together indicates stress, not support (observed on
PRIMARY at strength 2).

Status at recording time: stage 13 (outcome estimation) has not run;
no treatment-outcome estimate for any contrast has been computed or
seen. Support maps completed for C1 and PRIMARY; C2-C4 maps queued.

Recorded by: design stage (A. Mertens / automated pipeline assistant).
