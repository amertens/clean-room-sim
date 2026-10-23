#!/bin/bash
# Logged searches for the 2026-10-22 literature audit. Run from the project root.
S="C:/Users/andre/OneDrive/Documents/scientific-writing-skills/checks/lit_search.py"
O=audit/2026-10-22/lit
SRC=pubmed,europepmc,openalex,semanticscholar
export PYTHONIOENCODING=utf-8

run() {  # name, then block args
  name="$1"; shift
  python "$S" search "$@" --sources $SRC --max 30 --date-from 2019-01-01 --run-label "$name" --out-prefix "$O/$name" --export ids > "$O/$name.stdout.txt" 2>&1
  sleep 3
}

run b_blind --block blind="\"blind analysis\"|\"blinded analysis\"|\"blind data analysis\"|\"outcome-blind\"|\"outcome blind\"|\"blinded to the outcome\"|\"outcome-free\"" --block ctx="observational|real-world|causal|epidemiolog*|propensity|registry|\"social science\"|physics"
run b_design --block design="\"design trumps analysis\"|\"design stage\"|\"design phase\"|\"outcome-free design\"|\"before examining outcomes\"|\"without access to outcome*\"" --block ctx="\"observational stud*\"|\"causal inference\"|\"propensity score\"|\"real-world evidence\""
run c_prereg --block prereg="\"protocol template\"|preregistration|pre-registration|\"registered protocol\"|\"pre-specified\"|prespecified|\"protocol registration\"" --block rwe="\"real-world evidence\"|\"real-world data\"|pharmacoepidemiolog*|\"non-interventional\"|\"observational stud*\""
run c_target --block tte="\"target trial emulation\"|\"target trial\"" --block rep="protocol|\"reporting guideline\"|TARGET|preregist*|transparency|registration"
run d_duplicate --block dup="\"RCT-DUPLICATE\"|\"RCT DUPLICATE\"|\"emulat* randomized\"|\"emulation of randomized\"|\"trial emulation\" benchmark" --block ctx="claims|\"real-world\"|database|agreement|concordance"
run e_ohdsi --block ohdsi="OHDSI|\"Observational Health Data Sciences\"|LEGEND|\"empirical calibration\"|\"study diagnostics\"" --block ctx="\"negative control*\"|diagnostic*|\"large-scale\"|reproducib*|unblind*|calibrat*"
run e_sentinel --block sent="Sentinel|\"FDA Sentinel\"" --block tl="\"targeted learning\"|TMLE|\"targeted maximum likelihood\"|\"targeted minimum loss\"|\"causal roadmap\"|\"super learner\""
run e_software --block soft="\"R package\"|\"Python package\"|software|toolkit" --block plan="\"analysis plan\"|preregist*|\"pre-specified\"|prespecified|\"blind analysis\"|\"unblind*\"|\"design stage\"" --block caus="\"real-world\"|observational|causal|propensity"
run f_plasmode --block pl="plasmode|\"synthetic data\"|\"outcome simulation\"|\"semi-synthetic\"|\"semisynthetic\"" --block est="\"estimator selection\"|\"select* an estimator\"|\"method selection\"|\"choose* the estimator\"|\"selecting among estimators\"|\"data-adaptive\"|\"prespecified estimator\"" --block caus="causal|\"treatment effect\"|propensity|TMLE"
run f_sl --block sl="\"super learner\"|SuperLearner|\"Super Learner\"" --block spec="\"pre-specif*\"|prespecif*|\"specifying a\"|\"practical considerations\"|\"library of learners\"|\"statistical analysis plan\""
run f_roadmap --block rm="\"causal roadmap\"|\"Causal Roadmap\"|\"targeted learning\"" --block ctx="\"real-world\"|\"statistical analysis plan\"|pre-specif*|prespecif*|simulation|registry"
run g_compete --block ob="\"outcome-blind\"|\"outcome blind\"|\"blinded to outcomes\"|\"before the outcome\"|\"without the outcome\"|\"outcome-free\"" --block sim="simulation|plasmode|\"synthetic outcome*\"" --block est="estimand|estimator|\"treatment effect\"|TMLE|\"causal inference\""
run g_estimand_choice --block est="\"estimand selection\"|\"choice of estimand\"|\"choosing the estimand\"|\"selecting the estimand\"|\"estimand choice\"" --block pos="positivity|overlap|\"design-based\"|\"feasib*\"" --block caus="causal|observational|\"treatment effect\""
echo done > "$O/run_searches.done"
