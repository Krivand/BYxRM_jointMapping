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

# This is a whitelist, not a recency heuristic. Only outputs supporting the
# current reviewer response remain under review/. Superseded exploratory
# analyses are removed even if their files have recent modification times.
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
    "segregant_panel_overlap.tsv",
    "segregant_dependence_locus_results.tsv",
    "segregant_dependence_summary.tsv",
    "segregant_dependence_validation_summary.tsv"
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

# No current reviewer analysis requires an intermediate under imed/review/.
# The former interaction_power and clone-dependence RDS files are superseded.
imed_allowed = list()

asys_allowed_directories = c(names(asys_allowed), "power")
imed_allowed_directories = names(imed_allowed)

####################
# CLEANING HELPERS #
####################

removed_state = new.env(parent = emptyenv())
removed_state$paths = character(0)

remove_path = function(path) {
  if (file.exists(path) || dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
    removed_state$paths = c(removed_state$paths, path)
  }
  invisible(NULL)
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

clean_top_level = function(root, allowed_directories) {

  entries = list.files(
    root,
    full.names = TRUE,
    recursive = FALSE,
    all.files = FALSE
  )
  if (length(entries) == 0L) return(invisible(NULL))

  is_directory = file.info(entries)$isdir

  # Current review analyses write no direct analysis files at the review root.
  # The manifest is recreated after cleanup, so all direct files here are stale.
  direct_files = entries[!is_directory]
  if (length(direct_files) > 0L) {
    for (item in direct_files) remove_path(item)
  }

  directories = entries[is_directory]
  stale_directories = directories[
    !(basename(directories) %in% allowed_directories)
  ]
  if (length(stale_directories) > 0L) {
    for (item in stale_directories) remove_path(item)
  }

  invisible(NULL)
}

########################
# REMOVE KNOWN STALE   #
########################

if (clean_outputs) {

  clean_top_level(
    review_asys,
    asys_allowed_directories
  )
  clean_top_level(
    review_imed,
    imed_allowed_directories
  )

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
  # audit used to justify the LOD sensitivity range. Everything else there is
  # superseded.
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

removed = removed_state$paths

########################
# BUILD OUTPUT MANIFEST #
########################

manifest_path = paste0(
  review_asys,
  "review_output_manifest.tsv"
)

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
  rel_path = sub("^/+", "", rel_path)
  top_directory = sub("/.*$", "", rel_path)

  data.frame(
    path = rel_path,
    top_directory = top_directory,
    size_bytes = info$size,
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    status = "current",
    stringsAsFactors = FALSE
  )

} else {

  data.frame(
    path = character(0),
    top_directory = character(0),
    size_bytes = numeric(0),
    modified = character(0),
    status = character(0),
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
cat("Current tables found: ", nrow(manifest), "\n", sep = "")
cat("Stale paths removed: ", length(removed), "\n", sep = "")

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
