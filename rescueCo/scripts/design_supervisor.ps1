# Runs what is left of the Rescue.Co design stage unattended, as decided on
# 23 October 2026 (spec 2026-10-20, section 10). Started by the scheduled task
# RescueCo_design_stage (register_design_task.ps1), at logon and on demand, so
# it survives a closed terminal and resumes after a restart: every step is
# checkpointed and skipped once done, and the files of a lock interrupted
# before its design report are moved to an abandoned_<date> folder, never
# overwritten. One R process at a time:
#   1. the build fingerprint (script 22c), if none is recorded;
#   2. the hash check of every finished lock in both runs (script 22, verify);
#   3. the final cohort's primary comparisons L1, L2, L4, then L5;
#   4. the confirmatory comparison (script 23, pick) and its 180-day lock
#      (script 22b, check then run);
#   5. the sensitivity locks S1, S4, S2, S3, then the first four months' L1;
#   6. the registration bundle and approvals template (script 23, final).
# Workers per lock are set from the memory free at launch, up to $maxWorkers.
# An out-of-memory failure is retried with fewer workers; any other failure
# stops the supervisor. The computer is kept awake (not the display) while it
# runs. It reads no outcome.

$root = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path   # the project folder this copy sits in
$rs = @('C:\Program Files\R\R-4.4.2\bin\Rscript.exe',
        (Join-Path $env:LOCALAPPDATA 'Programs\R\R-4.4.2\bin\Rscript.exe')) |
  Where-Object { Test-Path $_ } | Select-Object -First 1   # R 4.4.2, installed for all users or for this user
# On the desktop, the laptop's versions of glmnet, xgboost, earth and future sit
# in their own library, ahead of the user library, so other R work keeps its own.
$pinnedLib = Join-Path $env:LOCALAPPDATA 'R\rescueCo-design-lib-4.4.2'
if (Test-Path $pinnedLib) { $env:R_LIBS = $pinnedLib }
$logs = Join-Path $root 'rescueCo\logs'
New-Item -ItemType Directory -Force -Path $logs | Out-Null      # git does not carry it
$log = Join-Path $logs 'supervisor.txt'
$locks = Join-Path $root 'rescueCo\data\design\locks'
$reg = Join-Path $root 'rescueCo\registration'
$headroomGB = 3.0     # left free for other programs and the system
$mainGB = 0.6         # the R process that coordinates the workers
$maxWorkers = 14      # of 22 logical processors
$reps = 100
function Say($msg) { Add-Content -Path $log -Value ("{0}  {1}" -f (Get-Date -Format s), $msg) }
function FreeGB { (Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1MB }

Add-Type -Namespace Win32 -Name Power -MemberDefinition `
  '[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);'
[void][Win32.Power]::SetThreadExecutionState([uint32]2147483649)   # ES_CONTINUOUS | ES_SYSTEM_REQUIRED

# R scripts started by something other than this supervisor; they share the CPU.
function OtherR {
  @(Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" |
    Where-Object { $_.CommandLine -match '\.R\b' -and $_.CommandLine -notmatch 'rescueCo/scripts/' }).Count
}

function RunR([string[]]$argList, [string]$tag) {
  $out = Join-Path $logs ($tag + '_stdout.txt')
  $err = Join-Path $logs ($tag + '_stderr.txt')
  $p = Start-Process -FilePath $rs -WorkingDirectory $root -WindowStyle Hidden -PassThru -Wait `
    -ArgumentList $argList -RedirectStandardOutput $out -RedirectStandardError $err
  return $p.ExitCode
}

function Stop-Supervisor($msg) {
  Say ($msg + '; supervisor stopped')
  [void][Win32.Power]::SetThreadExecutionState([uint32]2147483648)
  exit 1
}

# Moves the files of a lock that has no design report, left by an interrupted
# run, so that the new run neither overwrites nor reuses them.
function SetAside([string]$run, [string]$lock) {
  $dir = Join-Path $locks $run
  $files = @(Get-ChildItem -Path $dir -File -ErrorAction SilentlyContinue |
    Where-Object { $_.Name -match ('^{0}_(lock|design|simulation|nc)\.rds$' -f [regex]::Escape($lock)) })
  if ($files.Count -eq 0) { return }
  $dest = Join-Path $dir ('abandoned_{0}' -f (Get-Date -Format 'yyyy-MM-dd_HHmm'))
  New-Item -ItemType Directory -Force -Path $dest | Out-Null
  foreach ($f in $files) { Move-Item -Path $f.FullName -Destination (Join-Path $dest $f.Name) }
  Say ('{0} {1}: set aside {2} file(s) of an interrupted run in {3}' -f $run, $lock, $files.Count, $dest)
}

# One lock: skipped if its design report exists; otherwise launched with the
# workers the free memory allows ({W} in the arguments).
function RunLock([string]$run, [string]$lock, [int]$n, [string[]]$template) {
  $dossier = Join-Path $locks ('{0}\{1}_dossier.rds' -f $run, $lock)
  if (Test-Path $dossier) { Say ('{0} {1}: design report exists, skipped' -f $run, $lock); return }
  $perWorker = 0.25 + 0.0003 * $n   # GB per worker, scaled from the operators lock
  $cap = $maxWorkers
  $tries = 0
  while ($true) {
    SetAside $run $lock
    $free = FreeGB
    while ($free -lt ($headroomGB + $mainGB + $perWorker)) {
      Say ('{0} {1}: only {2:N1} GB free; waiting' -f $run, $lock, $free)
      Start-Sleep -Seconds 600
      $free = FreeGB
    }
    $w = [int][Math]::Max(1, [Math]::Min($cap, [Math]::Floor(($free - $headroomGB - $mainGB) / $perWorker)))
    $other = OtherR
    if ($other -gt 0) { Say ('note: {0} other R job(s) running; they share the CPU' -f $other) }
    Say ('{0} {1}: launching with {2} worker(s), {3:N1} GB free' -f $run, $lock, $w, $free)
    $argList = $template | ForEach-Object { $_ -replace '\{W\}', "$w" }
    $tag = 'ds_{0}_{1}' -f $run, $lock
    $code = RunR $argList $tag
    if ($code -eq 0) { Say ('{0} {1}: finished' -f $run, $lock); return }
    $tries++
    $err = Join-Path $logs ($tag + '_stderr.txt')
    $text = Get-Content -Path $err -Raw -ErrorAction SilentlyContinue
    if ($text -match 'Out of memory|bad_alloc|cannot allocate' -and $tries -le 3) {
      $cap = [Math]::Max(1, $w - 2)
      Copy-Item $err ($err -replace '_stderr.txt$', ('_stderr_failed{0}.txt' -f $tries))
      Say ('{0} {1}: out of memory with {2} worker(s); retry {3} in 10 minutes with at most {4}' -f $run, $lock, $w, $tries, $cap)
      Start-Sleep -Seconds 600
      continue
    }
    Stop-Supervisor ('{0} {1}: failed with exit code {2} (see {3})' -f $run, $lock, $code, $err)
  }
}

Say 'supervisor started (decisions of 2026-10-23)'

# 1. The build fingerprint, before anything else runs on rescueCo/lib.
if (-not (Test-Path (Join-Path $reg 'build_fingerprint.csv'))) {
  if ((RunR @('rescueCo/scripts/22c_record_build.R') 'ds_22c') -ne 0) { Stop-Supervisor 'script 22c failed' }
  Say 'build fingerprint recorded; commit rescueCo/registration/build_fingerprint.csv'
}

# 2. Every finished lock must rebuild with the same hash under the current script.
foreach ($run in @('final', 'first4m')) {
  if ((RunR @('rescueCo/scripts/22_design_stage.R', $run, 'verify', "$reps") ('ds_verify_' + $run)) -ne 0) {
    Stop-Supervisor ('{0}: a finished lock did not reproduce its hash (see ds_verify_{0}_stderr.txt)' -f $run)
  }
  Say ('{0}: finished locks reproduce their hashes' -f $run)
}
# The first four months' L4 was interrupted on 22 October and is no longer run
# (the four-month run is limited to L3 and L1); its files are set aside.
SetAside 'first4m' 'L4'

# 3. The final cohort's primary comparisons (L3 is done).
foreach ($q in @(@{lock = 'L1'; n = 6494}, @{lock = 'L2'; n = 4207},
                 @{lock = 'L4'; n = 3739}, @{lock = 'L5'; n = 4207})) {
  RunLock 'final' $q.lock $q.n @('rescueCo/scripts/22_design_stage.R', 'final', 'run', "$reps", '{W}', $q.lock)
}

# 4. The confirmatory comparison by the rule in the plan, and its 180-day lock.
$conf = Join-Path $reg 'confirmatory.csv'
if (-not (Test-Path $conf)) {
  if ((RunR @('rescueCo/scripts/23_registration_bundle.R', 'pick') 'ds_23_pick') -ne 0) { Stop-Supervisor 'script 23 (pick) failed' }
}
$pick = (Import-Csv $conf).lock
Say ('confirmatory comparison: {0}' -f $pick)
if ($pick -and $pick -ne 'NA') {
  $sec = $pick + '_death180'
  if (-not (Test-Path (Join-Path $locks ('final\{0}_dossier.rds' -f $sec)))) {
    if ((RunR @('rescueCo/scripts/22b_secondary_lock.R', 'check') 'ds_22b_check') -ne 0) { Stop-Supervisor 'script 22b (check) failed' }
  }
  RunLock 'final' $sec 1204 @('rescueCo/scripts/22b_secondary_lock.R', 'run', '{W}')
}

# 5. The sensitivity locks, then the first four months' protocol comparison.
foreach ($q in @(@{lock = 'S1'; n = 1204}, @{lock = 'S4'; n = 1195},
                 @{lock = 'S2'; n = 3478}, @{lock = 'S3'; n = 3946})) {
  RunLock 'final' $q.lock $q.n @('rescueCo/scripts/22_design_stage.R', 'final', 'run', "$reps", '{W}', $q.lock)
}
RunLock 'first4m' 'L1' 1661 @('rescueCo/scripts/22_design_stage.R', 'first4m', 'run', "$reps", '{W}', 'L1')

# 6. The bundle with every hash, and the approvals template.
if ((RunR @('rescueCo/scripts/23_registration_bundle.R', 'final') 'ds_23_final') -ne 0) {
  Stop-Supervisor 'script 23 (final) did not finish (see ds_23_final_stderr.txt)'
}
Say 'all steps done: review the design reports and the bundle, then rehearse script 25'
[void][Win32.Power]::SetThreadExecutionState([uint32]2147483648)   # ES_CONTINUOUS: sleep allowed again
