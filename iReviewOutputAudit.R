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
# CANONICAL REVIEW OUTPUT LAYOUT       #
########################################

canonical_directories = c(
  "go",
  "genes",
  "missingness",
  "interaction_dependence",
  "interaction_power",
  "power_effect_curve",
  "power"
)

# The power_effect_curve directory is the current individual-QTL power
# analysis. The older power, pair-conditional and reviewer-summary analyses
# were superseded. The remaining review/power directory is retained only for
# the historical permutation-threshold audit (threshold_*.tsv).
power_effect_allowed = c(
  "critical_r_power_table.tsv",
  "albert_reference_crossstudy.tsv",
  "albert_2018_power_benchmark.tsv",
  "effect_curve_validation_summary.tsv"
)

legacy_directories = c(
  paste0(review_asys, "power_pair_conditional/"),
  paste0(review_asys, "power_reviewer_summary/"),
  paste0(review_imed, "power_pair_conditional/"),
  paste0(review_imed, "power_reviewer_summary/"),
  paste0(review_imed, "power_effect_curve/"),
  paste0(review_imed, "power/")
)

########################
# REMOVE KNOWN STALE   #
########################

removed = character(0)

if (clean_outputs) {

  for (path in legacy_directories) {
    if (dir.exists(path)) {
      unlink(path, recursive = TRUE, force = TRUE)
      removed = c(removed, path)
    }
  }

  # review/power now belongs only to the threshold audit. Remove outputs from
  # superseded power implementations while preserving threshold audit tables.
  threshold_dir = paste0(review_asys, "power/")
  if (dir.exists(threshold_dir)) {
    threshold_files = list.files(
      threshold_dir,
      full.names = TRUE,
      recursive = FALSE,
      all.files = FALSE
    )
    threshold_keep = grepl(
      "^threshold_.*\\.tsv$",
      basename(threshold_files)
    )
    stale_threshold_dir = threshold_files[!threshold_keep]
    if (length(stale_threshold_dir) > 0L) {
      unlink(stale_threshold_dir, recursive = TRUE, force = TRUE)
      removed = c(removed, stale_threshold_dir)
    }
  }

  # The current effect-curve directory has four direct audit/analysis tables
  # plus the dedicated Figure 3 calibration subdirectory. Remove any direct
  # files left by earlier iterations.
  effect_dir = paste0(review_asys, "power_effect_curve/")
  if (dir.exists(effect_dir)) {
    effect_files = list.files(
      effect_dir,
      full.names = TRUE,
      recursive = FALSE,
      all.files = FALSE
    )
    effect_is_dir = file.info(effect_files)$isdir
    effect_keep =
      basename(effect_files) %in% power_effect_allowed |
      (effect_is_dir & basename(effect_files) == "figure3_calibration")
    stale_effect = effect_files[!effect_keep]
    if (length(stale_effect) > 0L) {
      unlink(stale_effect, recursive = TRUE, force = TRUE)
      removed = c(removed, stale_effect)
    }
  }
}

########################
# BUILD OUTPUT MANIFEST #
########################

manifest_path = paste0(
  review_asys,
  "review_output_manifest.tsv"
)

# Recreate the manifest itself on each audit rather than listing yesterday's
# manifest as though it were an analysis output.
if (file.exists(manifest_path)) {
  unlink(manifest_path, force = TRUE)
}

review_tables = list.files(
  review_asys,
  pattern = "\\.(tsv|csv)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

manifest = if (length(review_tables) > 0L) {

  info = file.info(review_tables)
  rel_path = substring(
    review_tables,
    nchar(review_asys) + 1L
  )
  top_directory = sub("/.*$", "", rel_path)

  data.frame(
    path = rel_path,
    top_directory = top_directory,
    size_bytes = info$size,
    modified = format(
      info$mtime,
      "%Y-%m-%d %H:%M:%S"
    ),
    canonical_review_directory =
      top_directory %in% canonical_directories,
    stringsAsFactors = FALSE
  )

} else {

  data.frame(
    path = character(0),
    top_directory = character(0),
    size_bytes = numeric(0),
    modified = character(0),
    canonical_review_directory = logical(0),
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

unexpected = manifest[
  !manifest$canonical_review_directory,
  ,
  drop = FALSE
]

cat("\nREVIEW OUTPUT AUDIT\n")
cat("Review table root: ", review_asys, "\n", sep = "")
cat("Tables found: ", nrow(manifest), "\n", sep = "")
cat("Known stale paths removed: ", length(removed), "\n", sep = "")

if (length(removed) > 0L) {
  cat("\nRemoved:\n")
  cat(paste0("  ", removed, collapse = "\n"), "\n")
}

if (nrow(unexpected) > 0L) {
  cat("\nTables outside canonical review subdirectories:\n")
  print(unexpected, row.names = FALSE)
} else {
  cat("All discovered review tables are under canonical review subdirectories.\n")
}

cat("Manifest: ", manifest_path, "\n", sep = "")
