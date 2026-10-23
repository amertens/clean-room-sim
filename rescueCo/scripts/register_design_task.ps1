# Registers the scheduled task RescueCo_design_stage, which runs
# rescueCo/scripts/design_supervisor.ps1 at logon (so an interrupted run
# resumes) and on demand. Settings follow SUSTAIN6_sim_supervisor: no time
# limit, runs on battery, one instance at a time, your own account without
# elevation. Decision of 2026-10-23 (spec 2026-10-20, section 10).
#
# Run once, from your own PowerShell window in the project folder:
#   powershell -NoProfile -ExecutionPolicy Bypass -File rescueCo\scripts\register_design_task.ps1
# Start now:  Start-ScheduledTask -TaskName RescueCo_design_stage
# Progress:   Get-Content rescueCo\logs\supervisor.txt -Tail 20
# Stop:       Stop-ScheduledTask -TaskName RescueCo_design_stage, then end the R
#             process it started:
#             Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" |
#               Where-Object { $_.CommandLine -match 'rescueCo/scripts/2' } |
#               ForEach-Object { Stop-Process -Id $_.ProcessId }
# Remove:     Unregister-ScheduledTask -TaskName RescueCo_design_stage -Confirm:$false

$root = 'C:\Users\andre\OneDrive\Documents\clean-room-sim'
$script = Join-Path $root 'rescueCo\scripts\design_supervisor.ps1'
if (-not (Test-Path $script)) { throw "No supervisor at $script" }
$user = '{0}\{1}' -f $env:USERDOMAIN, $env:USERNAME
$action = New-ScheduledTaskAction -Execute 'powershell.exe' -WorkingDirectory $root `
  -Argument ('-NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File "{0}"' -f $script)
$trigger = New-ScheduledTaskTrigger -AtLogOn -User $user
$settings = New-ScheduledTaskSettingsSet -ExecutionTimeLimit ([TimeSpan]::Zero) -AllowStartIfOnBatteries `
  -DontStopIfGoingOnBatteries -StartWhenAvailable -MultipleInstances IgnoreNew
$principal = New-ScheduledTaskPrincipal -UserId $user -LogonType Interactive -RunLevel Limited
Register-ScheduledTask -TaskName 'RescueCo_design_stage' -Action $action -Trigger $trigger -Settings $settings `
  -Principal $principal -Description 'Rescue.Co design stage (clean-room-sim): rescueCo/scripts/design_supervisor.ps1'
Write-Output 'Registered RescueCo_design_stage. Start it with: Start-ScheduledTask -TaskName RescueCo_design_stage'
