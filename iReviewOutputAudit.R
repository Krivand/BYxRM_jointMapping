########################################
########################################
## GENETICS REVIEW OUTPUT AUDIT       ##
## CLEAN STALE REVIEWER ARTIFACTS     ##
########################################
########################################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
review_imed = paste0(imed, "review/")
ensure_directory(review_asys)
ensure_directory(review_imed)

clean_outputs = Sys.getenv("REVIEW_OUTPUT_CLEAN", "1") == "1"

########################################
# CURRENT REVIEW OUTPUT CONTRACT       #
########################################

# Each current review analysis owns one subdirectory. Keeping the contract
# explicit makes stale files from superseded iterations visible and removable.
asys_allowed = list(
  go = c(
    "GO_significant_terms.tsv",
    "GO_key_terms.tsv",
    "GO_discordant_vs_other_interactions.tsv"
  ),
  genes = c(
    "supported_gene_characteristics_summary.tsv"
  ),
  missingness = c(
    "missingness_filter_summary.tsv"
  ),
  interaction_dependence = c(
    "clone_clustered_locus_results.tsv",
    "clone_clustered_summary.tsv",
    "clone_clustered_power.tsv",
    "clone_clustered_prevalence.tsv",
    "strict_qtl_clone_sensitivity.tsv"
  ),
  interaction_power = c(
    "interaction_design_868.tsv",
    "lrt_wald_calibration.tsv",
    "interaction_thresholds.tsv",
    "interaction_power_by_magnitude.tsv",
    "interaction_power_by_magnitude_prevalence.tsv",
    "observed_interaction_summary.tsv"
  ),
  power_effect_curve = c(
    "critical_r_power_table.tsv",
    "albert_reference_crossstudy.tsv",
    "albert_2018_power_benchmark.tsv",
    "effect_curve_validation_summary.tsv",
    "figure3_calibration"
  )
)

figure3_allowed = c(
  "albert_figure3_observed_vs_power.tsv",
  "albert_figure3_power_fit_summary.tsv",
  "albert_figure3_observed_vs_power_lod4_5.pdf",
  "albert_figure3_observed_vs_power_lod4_5.png"
)

imed_allowed = list(
  interaction_dependence = c(
    "review_interaction_dependence.rds"
  ),
  interaction_power = c(
    "review_interaction_power.rds"
  )
)

legacy_directories = c(
  paste0(review_asys, "power_pair_conditional/"),
  paste0(review_asys, "power_reviewer_summary/"),
  paste0(review_imed, "power_pair_conditional/"),
  paste0(review_imed, "power_reviewer_summary/"),
  paste0(review_imed, "power_effect_curve/"),
  paste0(review_imed, "power/")
)

####################
# CLEANING HELPERS #
####################

removed = character(0)

remove_path = function(path) {
  if (file.exists(path) || dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
    removed <<- c(removed, path)
  }
}

clean_exact_directory = function(path, allowed) {

  if (!dir.exists(path)) return(invisible(NULL))

  entries = list.files(
    path,
    full.names = TRUE,
    recursive = FALSE,
    all.files = FALSE
  )
  stale = entries[!(basename(entries) %in% allowed)]

  if (length(stale) > 0L) {
    for (item in stale) remove_path(item)
  }

  invisible(NULL)
}

########################
# REMOVE KNOWN STALE   #
########################

if (clean_outputs) {

  for (path in legacy_directories) remove_path(path)

  for (directory in names(asys_allowed)) {
    clean_exact_directory(
      paste0(review_asys, directory, "/"),
      asys_allowed[[directory]]
    )
  }

  # Figure 3 calibration is nested under the current individual-QTL power
  # analysis and is recreated by iReviewPowerFigure3Compare.R on every run.
  clean_exact_directory(
    paste0(review_asys, "power_effect_curve/figure3_calibration/"),
    figure3_allowed
  )

  # review/power is retained only for the historical permutation-threshold
  # audit. Superseded individual-QTL power outputs in this directory are stale.
  threshold_dir = paste0(review_asys, "power/")
  if (dir.exists(threshold_dir)) {
    threshold_entries = list.files(
      threshold_dir,
      full.names = TRUE,
      recursive = FALSE,
      all.files = FALSE
    )
    keep_threshold = grepl(
      "^threshold_.*\\.tsv$",
      basename(threshold_entries)
    )
    stale_threshold = threshold_entries[!keep_threshold]
    if (length(stale_threshold) > 0L) {
      for (item in stale_threshold) remove_path(item)
    }
  }

  for (directory in names(imed_allowed)) {
    clean_exact_directory(
      paste0(review_imed, directory, "/"),
      imed_allowed[[directory]]
    )
  }
}

########################
# BUILD OUTPUT MANIFEST #
########################

manifest_path = paste0(
  review_asys,
  "review_output_manifest.tsv"
)

# Recreate the manifest itself on every audit rather than listing yesterday's
# manifest as though it were an analysis result.
if (file.exists(manifest_path)) unlink(manifest_path, force = TRUE)

review_tables = list.files(
  review_asys,
  pattern = "\\.(tsv|csv)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

manifest = if (length(review_tables) > 0L) {

  info = file.info(review_tables)
  rel_path = substring(review_tables, nchar(review_asys) + 1L)
  top_directory = sub("/.*$", "", rel_path)

  data.frame(
    path = rel_path,
    top_directory = top_directory,
    size_bytes = info$size,
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )

} else {

  data.frame(
    path = character(0),
    top_directory = character(0),
    size_bytes = numeric(0),
    modified = character(0),
    stringsAsFactors = FALSE
  )
}

manifest = manifest[order(manifest$path), , drop = FALSE]

write.table(
  manifest,
  file = manifest_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

########################
# FINAL AUDIT           #
########################

cat("\nREVIEW OUTPUT AUDIT\n")
cat("Review table root: ", review_asys, "\n", sep = "")
cat("Tables found: ", nrow(manifest), "\n", sep = "")
cat("Known stale paths removed: ", length(removed), "\n", sep = "")

if (length(removed) > 0L) {
  cat("\nRemoved:\n")
  cat(paste0("  ", removed, collapse = "\n"), "\n")
}

cat("\nCurrent review tables:\n")
if (nrow(manifest) > 0L) {
  print(manifest, row.names = FALSE)
} else {
  cat("  none\n")
}

cat("Manifest: ", manifest_path, "\n", sep = "")
