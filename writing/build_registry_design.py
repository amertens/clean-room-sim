"""Add the Rescue.Co design-stage numbers to writing/registry.json.

Reads the outcome-free tables written by the design stage and by
writing/export_results.R and writes one entry per number, each pointing to its
source cell. Entries this script does not generate are kept. Rerun after
writing/export_results.R whenever a design-stage run finishes:
    python writing/build_registry_design.py
"""
import csv
import math
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REG = ROOT / "writing" / "registry.json"
EXP = "writing/results-export"
RES = "rescueCo/results/reexec_050"

RUNS = {"final": "final", "first4m": "first_four_months"}
LOCKS = {
    "L1": "protocol_all_arrivals", "L2": "protocol_direct", "L3": "operators",
    "L4": "rescueco_vs_other_means", "L5": "ambulance_vs_other_means",
    "S1": "operators_arrival_state", "S2": "rescueco_vs_other_means_no_walkin",
    "S3": "ambulance_vs_other_means_no_walkin", "S4": "operators_no_probable_match",
}
EST = {"trimmed_ATE": "trimmed_ate", "ATT": "att", "ATO": "ato", "ATE": "ate"}


def rows(rel):
    with open(ROOT / rel, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fmt_int(x):
    return f"{round(float(x)):,}"


def fmt_pct(x):
    return f"{100 * float(x):.1f}%"


def fmt(x, d):
    return f"{float(x):.{d}f}"


def sig3(x):
    """Three significant digits, keeping trailing zeros (0.0551, 0.119)."""
    x = float(x)
    d = max(0, 2 - math.floor(math.log10(abs(x)))) if x else 3
    return f"{x:.{d}f}"


entries = {}


def add(key, text, value, source):
    entries[key] = {"text": text, "value": float(value), "source": source}


# Cohort counts (registry-wide, and the first four months).
for r in rows(f"{EXP}/cohort_counts.csv"):
    key = f"{RUNS[r['run']]}_count_{r['item']}"
    add(key, fmt_int(r["n"]), r["n"], f"{EXP}/cohort_counts.csv#where=run:{r['run']}&where=item:{r['item']}&col=n")

# Overlap and cohort size by comparison.
for run, rs in RUNS.items():
    src = f"{RES}/{run}/design/design_summary.csv"
    for r in rows(src):
        lk = LOCKS[r["lock"]]
        loc = f"{src}#where=lock:{r['lock']}&col="
        for col in ("n", "n_treated", "n_comparator", "trimmed_kept_treated", "trimmed_kept_comparator"):
            add(f"{rs}_{lk}_{col}", fmt_int(r[col]), r[col], loc + col)
        add(f"{rs}_{lk}_share_outside_band", fmt_pct(r["share_outside_band"]), r["share_outside_band"],
            loc + "share_outside_band")
        add(f"{rs}_{lk}_c_statistic", fmt(r["c_statistic"], 2), r["c_statistic"], loc + "c_statistic")
        w = float(r["ate_max_weight"])
        add(f"{rs}_{lk}_ate_max_weight", fmt(w, 1) if w < 100 else fmt_int(w), w, loc + "ate_max_weight")
        add(f"{rs}_{lk}_att_ess_comparator", fmt_int(r["att_ess_comparator"]), r["att_ess_comparator"],
            loc + "att_ess_comparator")

# The population gate, by comparison and estimand below the ATE.
for r in rows(f"{EXP}/population_gate.csv"):
    rs, lk, e = RUNS[r["run"]], LOCKS[r["lock"]], EST[r["estimand"]]
    loc = f"{EXP}/population_gate.csv#where=run:{r['run']}&where=lock:{r['lock']}&where=estimand:{r['estimand']}&col="
    if r["max_abs_shift"] != "NA":   # NA: fewer than five patients removed, shift not shown
        add(f"{rs}_{lk}_{e}_max_shift", fmt(r["max_abs_shift"], 2), r["max_abs_shift"], loc + "max_abs_shift")
    if r["estimand"] == "trimmed_ATE":
        add(f"{rs}_{lk}_trimmed_share_kept_treated", fmt_pct(r["kept_treated"]), r["kept_treated"], loc + "kept_treated")
        add(f"{rs}_{lk}_trimmed_share_kept_comparator", fmt_pct(r["kept_comparator"]), r["kept_comparator"],
            loc + "kept_comparator")

# Follow-up: patients before the complete-case restriction and the share of each arm lost.
for lock in LOCKS:
    src = f"{RES}/final/loss_by_arm_{lock}.csv"
    if not (ROOT / src).exists():
        continue
    for r in rows(src):
        arm = r["arm"]
        loc = f"{src}#where=arm:{arm}&col="
        add(f"final_{LOCKS[lock]}_{arm}_before_complete_case", fmt_int(r["patients"]), r["patients"], loc + "patients")
        add(f"final_{LOCKS[lock]}_{arm}_six_month_not_recorded", fmt_pct(r["six_month_not_recorded"]),
            r["six_month_not_recorded"], loc + "six_month_not_recorded")
        add(f"final_{LOCKS[lock]}_{arm}_disposition_not_recorded", fmt_pct(r["disposition_not_recorded"]),
            r["disposition_not_recorded"], loc + "disposition_not_recorded")

# Numerators: patients outside the band, and patients lost by arm.
for r in rows(f"{EXP}/overlap_counts.csv"):
    rs, lk = RUNS[r["run"]], LOCKS[r["lock"]]
    loc = f"{EXP}/overlap_counts.csv#where=run:{r['run']}&where=lock:{r['lock']}&col="
    for col in ("outside_treated", "outside_comparator", "outside"):
        add(f"{rs}_{lk}_n_{col}", fmt_int(r[col]), r[col], loc + col)
for r in rows(f"{EXP}/loss_counts.csv"):
    lk = LOCKS[r["lock"]]
    loc = f"{EXP}/loss_counts.csv#where=lock:{r['lock']}&where=arm:{r['arm']}&col="
    for col in ("six_month_not_recorded", "disposition_not_recorded"):
        add(f"final_{lk}_{r['arm']}_n_{col}", fmt_int(r[col]), r[col], loc + col)

# The outcome-blind simulation: ranges by learner library, overall and by estimand.
for name, by_est in (("simulation_by_library.csv", False), ("simulation_by_library_estimand.csv", True)):
    for r in rows(f"{EXP}/{name}"):
        rs, lk = RUNS[r["run"]], LOCKS[r["lock"]]
        lib = {"glm": "main_terms", "sl": "super_learner"}[r["library"]]
        stem = f"{rs}_{lk}_sim_{lib}" + (f"_{EST[r['estimand']]}" if by_est else "")
        loc = f"{EXP}/{name}#where=run:{r['run']}&where=lock:{r['lock']}&where=library:{r['library']}" + \
              (f"&where=estimand:{r['estimand']}" if by_est else "") + "&col="
        for col in ("worst_bias_min", "worst_bias_max", "se_mean_min", "se_mean_max", "emp_se_min", "emp_se_max"):
            add(f"{stem}_{col}", fmt(r[col], 3), r[col], loc + col)
        for col in ("halfwidth_min", "halfwidth_max"):   # compared with the tolerance, so three significant digits
            add(f"{stem}_{col}", sig3(r[col]), r[col], loc + col)
        for col in ("coverage_min", "coverage_max", "coverage_lower_min", "coverage_upper_max"):
            add(f"{stem}_{col}", fmt(r[col], 2), r[col], loc + col)
        add(f"{stem}_bias_mcse_max", fmt(r["bias_mcse_max"], 3), r["bias_mcse_max"], loc + "bias_mcse_max")
        for col in ("cells", "unresolved_cells"):
            add(f"{stem}_{col}", fmt_int(r[col]), r[col], loc + col)

# Simulated datasets per outcome model (by estimand, after any extension), and the number of outcome models.
for r in rows(f"{EXP}/simulation_summary.csv"):
    if r["candidate"].startswith("sl") and r["candidate"].endswith("0.01"):
        rs, lk, e = RUNS[r["run"]], LOCKS[r["lock"]], EST[r["estimand"]]
        loc = f"{EXP}/simulation_summary.csv#where=run:{r['run']}&where=lock:{r['lock']}&where=candidate:{r['candidate']}&where=estimand:{r['estimand']}&col="
        add(f"{rs}_{lk}_sim_datasets_{e}", fmt_int(r["reps_min"]), r["reps_min"], loc + "reps_min")
        if r["estimand"] == "ATE":
            add(f"{rs}_{lk}_sim_datasets_per_model", fmt_int(r["reps_min"]), r["reps_min"], loc + "reps_min")
            add(f"{rs}_{lk}_sim_outcome_models", fmt_int(r["models"]), r["models"], loc + "models")

# The selected candidate of each feasible estimand: its truncation level and expected CI half-width.
for r in rows(f"{EXP}/simulation_verdicts.csv"):
    if r["status"] != "feasible":
        continue
    rs, lk, e = RUNS[r["run"]], LOCKS[r["lock"]], EST[r["estimand"]]
    loc = f"{EXP}/simulation_verdicts.csv#where=run:{r['run']}&where=lock:{r['lock']}&where=estimand:{r['estimand']}&col="
    add(f"{rs}_{lk}_{e}_selected_halfwidth", sig3(r["ci_halfwidth"]), r["ci_halfwidth"], loc + "ci_halfwidth")
    add(f"{rs}_{lk}_{e}_selected_truncation", fmt(r["truncation"], 2), r["truncation"], loc + "truncation")

# Expected events per arm, before the outcome.
for r in rows(f"{EXP}/expected_events_range.csv"):
    rs, lk = RUNS[r["run"]], LOCKS[r["lock"]]
    loc = f"{EXP}/expected_events_range.csv#where=run:{r['run']}&where=lock:{r['lock']}&col="
    for col in ("treated_min", "treated_max", "control_min", "control_max"):
        add(f"{rs}_{lk}_expected_events_{col}", fmt_int(r[col]), r[col], loc + col)

# Negative controls: the null band.
for f in sorted((ROOT / EXP).glob("nc_verdict_*.csv")):
    r = rows(f"{EXP}/{f.name}")[0]
    rs, lk = RUNS[r["run"]], LOCKS[r["lock"]]
    loc = f"{EXP}/{f.name}#row=1&col="
    add(f"{rs}_{lk}_nc_null_band_upper", fmt(r["null_band_upper"], 2), r["null_band_upper"], loc + "null_band_upper")

# Negative controls.
for r in rows(f"{EXP}/negative_controls_summary.csv"):
    if r["rung"] != "full cohort" or r["estimand"] != "ATE":
        continue
    rs, lk = RUNS[r["run"]], LOCKS[r["lock"]]
    loc = f"{EXP}/negative_controls_summary.csv#where=run:{r['run']}&where=lock:{r['lock']}&where=rung:full cohort&where=estimand:ATE&col="
    for col in ("controls", "estimated", "inestimable", "inside", "outside", "inconclusive"):
        add(f"{rs}_{lk}_nc_{col}", fmt_int(r[col]), r[col], loc + col)
for r in rows(f"{EXP}/negative_controls.csv"):
    if r["status"] != "estimated" or r["rung"] != "full cohort":
        continue
    rs, lk, e = RUNS[r["run"]], LOCKS[r["lock"]], EST[r["estimand"]]
    nc = r["negative_control"].removeprefix("nc_")
    src = (f"{EXP}/negative_controls.csv#where=run:{r['run']}&where=lock:{r['lock']}&where=rung:full cohort"
           f"&where=estimand:{r['estimand']}&where=negative_control:{r['negative_control']}"
           f"&col=estimate&lower_col=ci_lower&upper_col=ci_upper")
    entries[f"{rs}_{lk}_nc_{nc}_{e}"] = {
        "text": f"{float(r['estimate']):.3f} (95% CI: {float(r['ci_lower']):.3f}, {float(r['ci_upper']):.3f})",
        "value": float(r["estimate"]), "lower": float(r["ci_lower"]), "upper": float(r["ci_upper"]),
        "source": src}

reg = json.loads(REG.read_text(encoding="utf-8-sig") or "{}")
# Drop entries this script generated earlier, so a number that disappears from the outputs leaves the registry too.
reg = {k: v for k, v in reg.items() if not k.startswith(tuple(f"{r}_" for r in RUNS.values()))}
reg.update(entries)
REG.write_text(json.dumps(dict(sorted(reg.items())), indent=1) + "\n", encoding="utf-8")
print(f"{len(entries)} entries written; registry now holds {len(reg)}")
