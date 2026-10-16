########################################
## FRANK REVIEW DEPENDENCE CLARITY     ##
########################################

# Summarize existing refits; no models are fitted.
# Usage: Rscript --vanilla iReviewFrankDependence.R P05_ROOT NEW_OUTPUT_DIRECTORY
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply p05 root and new output.")
root = normalizePath(args[[1]], mustWork = TRUE)
out = args[[2]]
if (file.exists(out)) stop("Output exists; refusing to overwrite.")
d = read.delim(file.path(root, "asys/review/interaction_dependence",
                         "segregant_dependence_locus_results.tsv"),
               stringsAsFactors = FALSE)
required = c("gene", "marker", "in_joint_868", "in_strict_64",
             "segregant_random_full_singular", "segregant_random_variance",
             "segregant_random_study_variance", "segregant_random_fit_success")
stopifnot(all(required %in% names(d)),
          !anyDuplicated(paste(d$gene, d$marker, sep = "::")),
          !anyNA(d$in_joint_868), !anyNA(d$in_strict_64),
          all(d$in_joint_868 | d$in_strict_64),
          sum(d$in_joint_868) == 868L, sum(d$in_strict_64) == 64L)
counts = data.frame(
  joint = sum(d$in_joint_868), strict = sum(d$in_strict_64),
  intersection = sum(d$in_joint_868 & d$in_strict_64),
  union = nrow(d), fit_success = sum(d$segregant_random_fit_success %in% TRUE),
  full_singular = sum(d$segregant_random_full_singular %in% TRUE),
  missing_singularity = sum(is.na(d$segregant_random_full_singular))
)
variance = do.call(rbind, lapply(
  c("segregant_random_variance", "segregant_random_study_variance"),
  function(column) {
    x = d[[column]]
    data.frame(component = column, n_finite = sum(is.finite(x)),
      n_missing = sum(!is.finite(x)), n_exact_zero = sum(x == 0, na.rm = TRUE),
      minimum = min(x, na.rm = TRUE), median = median(x, na.rm = TRUE),
      maximum = max(x, na.rm = TRUE))
  }
))
# Exact zero counts describe reported variances; they do not reproduce
# lme4's relative-scale isSingular tolerance or identify all boundary fits.
if (!dir.create(out, recursive = TRUE)) stop("Could not create new output.")
write.table(counts, file.path(out, "targeted_fit_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(variance, file.path(out, "random_intercept_variances.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(counts, row.names = FALSE)
print(variance, row.names = FALSE)
