cleanTMLE separates the design of an observational analysis from its outcome. You declare the plan, lock the design data without the outcome, learn from an outcome-free simulation which estimands the design can deliver and which TMLE candidate to use for each, record that decision in a hashed dossier, and only then admit the outcome.

This is the first public release.

**New in 0.5.0**

- The plan declares time zero, where follow-up starts, and when each eligibility criterion was assessed.
- The lock takes a variable dictionary: for each column, its source, when it was recorded and when the attribute it describes was fixed, and its role. The lock refuses a covariate fixed after time zero, a covariate from the follow-up interview, a node after treatment in the adjustment set, and a treatment not fixed at time zero, and it asks for a written response to each remaining flag.
- The dossier has an identification table listing what the plan declared and each flag with its response. It lists declarations; it does not grade identification.
- Each flag is also a row in the design log.

Locks made by earlier versions must be re-made with a variable dictionary. `time_zero` is now the third argument of `analysis_plan()`, so a call that passed later arguments by position must name them.

Planned additions are listed in `WISHLIST.md`. The full list of changes is in `NEWS.md`.

Install with `remotes::install_github("amertens/cleanTMLE@v0.5.0")`.
