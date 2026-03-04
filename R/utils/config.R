#' @title Configuration Management
#' @description Load and validate the project configuration from YAML.
#' @name config
NULL

#' Load Configuration
#'
#' Reads the YAML configuration file and returns a nested list.
#' Falls back to config/default.yml if no path is specified.
#'
#' @param config_path Character path to YAML config file.
#' @return Named list of configuration values.
#' @export
load_config <- function(config_path = NULL) {
  if (is.null(config_path)) {
    candidates <- c(
      file.path("config", "default.yml"),
      file.path("..", "config", "default.yml"),
      system.file("config", "default.yml", package = "cleanroom.sim")
    )
    config_path <- Find(file.exists, candidates)
    if (is.null(config_path)) {
      stop("No config file found. Provide config_path or place config/default.yml ",
           "in the project root.")
    }
  }
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }
  cfg <- yaml::yaml.load_file(config_path)
  validate_config(cfg)
  cfg
}

#' Validate Configuration
#'
#' Checks that required fields exist in the config.
#'
#' @param cfg Named list from \code{load_config}.
#' @return Invisibly returns TRUE on success; errors otherwise.
#' @keywords internal
validate_config <- function(cfg) {
  required_top <- c("dgp", "scenarios", "simulation", "tmle", "stage1",
                    "stage2", "output")
  missing <- setdiff(required_top, names(cfg))
  if (length(missing) > 0) {
    stop("Config missing required sections: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

#' Get Scenario Parameters
#'
#' Extract DGP toggle parameters for a named scenario, merged with base DGP
#' parameters. Only returns parameters accepted by \code{generate_hcv_data}.
#'
#' @param cfg Config list from \code{load_config}.
#' @param scenario_name Character name of the scenario.
#' @return Named list of DGP parameters suitable for \code{generate_hcv_data}.
#' @export
get_scenario_params <- function(cfg, scenario_name) {
  if (!scenario_name %in% names(cfg$scenarios)) {
    stop("Unknown scenario: ", scenario_name,
         ". Available: ", paste(names(cfg$scenarios), collapse = ", "))
  }
  base <- cfg$dgp
  toggles <- cfg$scenarios[[scenario_name]]
  # Merge toggles into base DGP params
  params <- modifyList(base, toggles)

  # Only keep arguments that generate_hcv_data() accepts
  dgp_args <- names(formals(generate_hcv_data))
  params[intersect(names(params), dgp_args)]
}

#' Ensure Output Directory Exists
#'
#' Creates the output directory if it does not exist.
#'
#' @param dir_path Character path to directory.
#' @return The directory path, invisibly.
#' @export
ensure_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(dir_path)
}
