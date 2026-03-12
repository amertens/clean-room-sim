test_that("Stage 1 produces report and checkpoint", {
  # Use temp directory for outputs
  tmp_dir <- tempfile("stage1_test")
  dir.create(tmp_dir, recursive = TRUE)
  out_dir <- file.path(tmp_dir, "stage1")

  cfg <- load_config()
  cfg$output$base_dir <- tmp_dir
  cfg$output$stage1_dir <- out_dir

  d <- generate_hcv_data(N = 1000, seed = 55, complexity = FALSE,
                         np_hazard = FALSE, dep_censor = FALSE,
                         switch_on = FALSE)
  result <- stage1_build_cohort(data = d, cfg = cfg, output_dir = out_dir)

  expect_true("cohort" %in% names(result))
  expect_true("report" %in% names(result))
  expect_true(result$checkpoint %in% c("PASS", "FAIL"))
  expect_true(file.exists(file.path(out_dir, "stage1_report.json")))
  expect_true(file.exists(file.path(tmp_dir, "checkpoint_1.json")))

  unlink(tmp_dir, recursive = TRUE)
})

test_that("Stage 3 cannot run if checkpoint_2 is FAIL", {
  tmp_dir <- tempfile("gate_test")
  dir.create(tmp_dir, recursive = TRUE)

  # Write a FAIL checkpoint
  write_checkpoint("checkpoint_2", "FAIL",
                   list(test = list(pass = FALSE)),
                   output_dir = tmp_dir)

  cfg <- load_config()
  cfg$output$base_dir <- tmp_dir
  cfg$output$stage3_dir <- file.path(tmp_dir, "stage3")

  d <- generate_hcv_data(N = 500, seed = 1, complexity = FALSE,
                         np_hazard = FALSE, dep_censor = FALSE,
                         switch_on = FALSE)

  expect_error(
    stage3_estimation(cohort = d, cfg = cfg,
                      output_dir = file.path(tmp_dir, "stage3")),
    "FAIL"
  )

  unlink(tmp_dir, recursive = TRUE)
})

test_that("Checkpoint read/write works correctly", {
  tmp_dir <- tempfile("checkpoint_test")
  dir.create(tmp_dir, recursive = TRUE)

  write_checkpoint("checkpoint_test", "PASS",
                   list(metric = list(value = 0.5, threshold = 0.3,
                                      pass = TRUE)),
                   output_dir = tmp_dir)

  cp <- read_checkpoint("checkpoint_test", output_dir = tmp_dir)
  expect_equal(cp$status, "PASS")
  expect_equal(cp$stage, "checkpoint_test")

  # Should not error
  expect_true(require_checkpoint_pass("checkpoint_test", output_dir = tmp_dir))

  # FAIL checkpoint should error
  write_checkpoint("checkpoint_fail", "FAIL",
                   list(metric = list(pass = FALSE)),
                   output_dir = tmp_dir)
  expect_error(require_checkpoint_pass("checkpoint_fail", output_dir = tmp_dir),
               "FAIL")

  unlink(tmp_dir, recursive = TRUE)
})

test_that("Decision log round-trips correctly", {
  tmp_dir <- tempfile("log_test")
  dir.create(tmp_dir, recursive = TRUE)
  log_path <- file.path(tmp_dir, "decision_log.csv")

  mtg <- start_meeting("stage1", approver = "test_user")
  mtg <- log_decision(mtg, "PS", "Use GLM", "Simple model for testing")
  mtg <- log_decision(mtg, "Q", "Use SL.glm", "Baseline model",
                      outcome_blind_confirm = TRUE)
  close_meeting(mtg, log_path = log_path)

  log <- utils::read.csv(log_path, stringsAsFactors = FALSE)
  expect_equal(nrow(log), 2)
  expect_true(all(c("meeting_id", "decision", "justification") %in%
                    names(log)))

  unlink(tmp_dir, recursive = TRUE)
})

test_that("Outcome blind guard blocks if Stage 3 exists", {
  tmp_dir <- tempfile("guard_test")
  s3_dir <- file.path(tmp_dir, "stage3")
  dir.create(s3_dir, recursive = TRUE)

  # Create a stage3 file
  writeLines("test", file.path(s3_dir, "stage3_estimates.csv"))

  # Create a decision log with stage3 at protocol_version = 1
  log_path <- file.path(tmp_dir, "decision_log.csv")
  log <- data.frame(
    meeting_id = "MTG-test", date_time = Sys.time(), stage = "stage3",
    model_component = "test", decision = "test", justification = "test",
    triggered_by = "test", outcome_blind_confirm = TRUE,
    approver = "test", protocol_version = 1,
    stringsAsFactors = FALSE
  )
  utils::write.csv(log, log_path, row.names = FALSE)

  # Should error with protocol_version = 1
  expect_error(
    guard_outcome_blind(s3_dir, current_protocol_version = 1L,
                        log_path = log_path),
    "OUTCOME BLINDNESS"
  )

  # Should warn (but not error) with protocol_version = 2
  expect_warning(
    guard_outcome_blind(s3_dir, current_protocol_version = 2L,
                        log_path = log_path)
  )

  unlink(tmp_dir, recursive = TRUE)
})
