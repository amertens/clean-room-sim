# Running the Rescue.Co design stage on the desktop (brief, 23 October 2026)

For a Claude Code session on Andrew's desktop, in a git clone of clean-room-sim (`amertens/clean-room-sim`, branch `main`). The task: run what is left of the Rescue.Co design stage there, through the scheduled task `RescueCo_design_stage`, and report how it starts. Andrew has asked for this run. Nothing else.

## Rules that hold on this machine too

- Outcome-blind. Never open `rescueCo/data/outcomes/`, `rescueCo/logs/pipeline.log`, `rescueCo/data/Final Data Cut 03.27.2026/6 Month Follow-up/`, files directly in `rescueCo/results/` (the `reexec_050/` subfolder is fine), `rescueCo/reports/`, `rescueCo/README_legacy.md`, or `../rescueCo_analysis*`. The design stage reads no outcome.
- Do not edit `rescueCo/lib` or `rescueCo/data`; do not install or update R packages without asking Andrew.
- The folders Andrew copies in hold patient-level registry data. They are ignored by git and must never be committed or pushed; the repository is public.
- One R job at a time on this machine. Do not commit, merge or push.
- The run must not also run on the laptop.

## Check before starting (report every mismatch to Andrew and stop)

1. The clone is on `main` at commit `ac41c94` or later, with `rescueCo/scripts/design_supervisor.ps1` finding the project folder from its own location (`$PSScriptRoot`), not from a fixed path. If it still has a fixed `C:\Users\andre\OneDrive\...` path, pull `main` again.
2. Andrew has copied these folders from the laptop into the same places in the clone, none of which git carries:
   - `rescueCo/data/design/`: `design_data_final.rds`, `design_data_first4m.rds`, `observation_indicators.rds`, and `locks/final/` and `locks/first4m/`, each with `L3_lock.rds`, `L3_design.rds`, `L3_simulation.rds`, `L3_nc.rds` and `L3_dossier.rds` (an `L4_lock.rds` and `L4_design.rds` left by the interrupted run are expected; the supervisor sets them aside);
   - `rescueCo/lib/cleanTMLE` (version 0.5.0, from its DESCRIPTION file);
   - `rescueCo/results/reexec_050/` (aggregate tables and the L3 design reports of both runs);
   - `rescueCo/logs/reexec_050.log`, the run's log so far (optional; the supervisor creates `rescueCo/logs/` if it is missing).
3. R 4.4.2 at `C:\Program Files\R\R-4.4.2`, the path the supervisor uses.
4. The packages the build fingerprint records match the laptop's versions: tmle 2.1.1, glmnet 4.1-8, xgboost 1.7.8.1, earth 5.3.4, nnet 7.3-19, SuperLearner 2.0-29, lmtp 1.5.4, generics 0.1.4, digest 0.6.39, future 1.75.0, future.apply 1.20.2, and cleanTMLE 0.5.0 from `rescueCo/lib`. Read the versions from each package's DESCRIPTION file. The finished L3 locks were built on the laptop with these versions, so a different version would mix two builds in one design stage.
5. `rescueCo/registration/build_fingerprint.csv` does not exist yet. The supervisor's first step records this machine's build, and scripts 22b and 25 stop unless the build matches it, so the outcome steward and the estimation after approval (scripts 24 and 25) must run on this machine as well. Tell Andrew this.
6. Free memory and processor: no other R job running. The supervisor sets workers from free memory (`$headroomGB = 3.0`, about 0.25 + 0.0003 n GB per worker) up to `$maxWorkers = 14`, which was set for the laptop's 22 logical processors. If this machine has more, propose a new cap to Andrew (logical processors minus about 6); do not change it without his yes.

## Start

From PowerShell in the clone's root folder:

    powershell -NoProfile -ExecutionPolicy Bypass -File rescueCo\scripts\register_design_task.ps1
    Start-ScheduledTask -TaskName RescueCo_design_stage

Then read `rescueCo\logs\supervisor.txt` (last 20 lines) and `rescueCo\logs\reexec_050.log` after a few minutes, and report: the fingerprint step, whether the hash check reproduced the L3 lock in both runs ("VERIFY final L3: hash reproduced", "VERIFY first4m L3: hash reproduced"), and how many workers the first lock got. The steps, in order: fingerprint; hash check of both runs; final L1, L2, L4, L5; the confirmatory pick and its 180-day lock; S1, S4, S2, S3; the first four months' L1; the registration bundle. Expect about 25 to 35 hours on 14 workers. The task resumes at logon after a restart.

If the hash check fails, stop the task (`Stop-ScheduledTask -TaskName RescueCo_design_stage`, then end its Rscript process as described at the top of `register_design_task.ps1`) and report: the lock or its plan does not rebuild the same way on this machine.

## After the run

The new locks, design reports and the registration bundle are written only in this clone. Andrew will need `rescueCo/results/reexec_050/`, `rescueCo/registration/` and the logs back on the laptop for the registry and the manuscript; `rescueCo/registration/` is tracked and is pushed by Andrew after his approval.
