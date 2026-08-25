########################################
########################################
# VALIDATE SEGREGANT DEPENDENCE       #
########################################
########################################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
dependence_asys = paste0(review_asys, "interaction_dependence/")
validation_path = paste0(
  dependence_asys,
  "segregant_dependence_validation_summary.tsv"
)

read_tsv = function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

load_single_rda = function(path, object_name) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  env = new.env(parent = emptyenv())
  loaded = load(path, envir = env)
  if (!(object_name %in% loaded)) {
    stop("Missing object ", object_name, " in ", path)
  }
  env[[object_name]]
}

canonical_clone = function(x) {
  x = trimws(as.character(x))
  x = gsub("_", ".", x, fixed = TRUE)
  gsub("\\s+", "", x)
}

check_state = new.env(parent = emptyenv())
check_state$rows = list()
record_check = function(name, passed, detail) {
  i = length(check_state$rows) + 1L
  check_state$rows[[i]] = data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

panel = read_tsv(
  paste0(dependence_asys, "segregant_panel_overlap.tsv")
)
locus = read_tsv(
  paste0(dependence_asys, "segregant_dependence_locus_results.tsv")
)
summary_table = read_tsv(
  paste0(dependence_asys, "segregant_dependence_summary.tsv")
)

jBH_analyte_geno = load_single_rda(
  paste0(imed, "jBH_analyte_geno.rda"),
  "jBH_analyte_geno"
)
jBH_analyte_ixn = load_single_rda(
  paste0(imed, "jBH_analyte_ixn.rda"),
  "jBH_analyte_ixn"
)

joint_keys = paste(
  jBH_analyte_geno$gene,
  jBH_analyte_geno$marker,
  sep = "::"
)
strict_keys = paste(
  jBH_analyte_ixn$gene,
  jBH_analyte_ixn$marker,
  sep = "::"
)
target_keys = unique(c(joint_keys, strict_keys))

primary_outputs = sort(
  setdiff(
    list.files(
      dependence_asys,
      pattern = "^segregant_.*\\.tsv$",
      full.names = FALSE
    ),
    "segregant_dependence_validation_summary.tsv"
  )
)
expected_outputs = sort(
  c(
    "segregant_panel_overlap.tsv",
    "segregant_dependence_locus_results.tsv",
    "segregant_dependence_summary.tsv"
  )
)
record_check(
  "analysis_emits_three_primary_tables",
  identical(primary_outputs, expected_outputs),
  paste(primary_outputs, collapse = ", ")
)

record_check(
  "targeted_locus_keys_are_exact",
  nrow(locus) == 928L &&
    !anyDuplicated(locus$key) &&
    setequal(locus$key, target_keys),
  paste0("rows=", nrow(locus), "; expected=928")
)

record_check(
  "joint_and_strict_membership_counts_are_exact",
  sum(locus$in_joint_868) == 868L &&
    sum(locus$in_strict_64) == 64L,
  paste0(
    "joint=", sum(locus$in_joint_868),
    "; strict=", sum(locus$in_strict_64)
  )
)

joint = locus[locus$in_joint_868, , drop = FALSE]
strict = locus[locus$in_strict_64, , drop = FALSE]

record_check(
  "original_joint_posthoc_count_is_reproduced",
  sum(joint$original_interaction_p < 0.05) == 149L,
  paste0("nominal joint interactions=", sum(joint$original_interaction_p < 0.05))
)

joint_error = max(
  abs(joint$original_interaction_p - joint$stored_joint_interaction_p),
  na.rm = TRUE
)
strict_error = max(
  abs(strict$original_interaction_p - strict$stored_strict_interaction_p),
  na.rm = TRUE
)
record_check(
  "original_targeted_refits_reproduce_stored_models",
  joint_error < 1e-8 && strict_error < 1e-8,
  paste0(
    "joint max p error=", joint_error,
    "; strict max p error=", strict_error
  )
)

alignment_error = max(
  locus$maximum_genotype_alignment_difference,
  na.rm = TRUE
)
record_check(
  "segregant_reconstruction_preserves_model_rows",
  is.finite(alignment_error) && alignment_error <= 1e-8,
  paste0("max genotype difference=", alignment_error)
)

robust_ok =
  is.finite(locus$robust_interaction_effect) &
  is.finite(locus$robust_interaction_se) &
  locus$robust_interaction_se > 0 &
  is.finite(locus$robust_interaction_p)
record_check(
  "cluster_robust_inference_is_finite_for_all_loci",
  all(robust_ok),
  paste0(sum(robust_ok), "/", nrow(locus), " loci finite")
)

record_check(
  "strict_interactions_are_stable_to_cluster_robust_inference",
  sum(strict$robust_interaction_p < 0.05) == 64L &&
    all(
      sign(strict$robust_interaction_effect) ==
      sign(strict$original_interaction_effect)
    ),
  paste0(
    "robust nominal strict interactions=",
    sum(strict$robust_interaction_p < 0.05),
    "; sign concordance=",
    mean(
      sign(strict$robust_interaction_effect) ==
      sign(strict$original_interaction_effect)
    )
  )
)

random_success = locus$clone_random_fit_success %in% TRUE
record_check(
  "clone_random_refits_are_fully_accounted_for",
  all(random_success) &&
    all(is.finite(locus$clone_random_interaction_p)) &&
    all(is.finite(locus$clone_random_interaction_effect)),
  paste0(sum(random_success), "/", nrow(locus), " successful")
)

record_check(
  "strict_interactions_are_stable_to_clone_random_refits",
  sum(strict$clone_random_interaction_p < 0.05) == 64L &&
    all(
      sign(strict$clone_random_interaction_effect) ==
      sign(strict$original_interaction_effect)
    ),
  paste0(
    "clone-random nominal strict interactions=",
    sum(strict$clone_random_interaction_p < 0.05),
    "; singular full fits=",
    sum(strict$clone_random_full_singular %in% TRUE)
  )
)

metadata_analyte = scaMetadata
studies_analyte = lapply(
  metadata_analyte$study,
  function(study) {
    load_single_rda(
      paste0(imed, study, ".rda"),
      study
    )
  }
)
names(studies_analyte) = metadata_analyte$study

clone_ids = lapply(
  studies_analyte,
  function(study) canonical_clone(rownames(study$gdata))
)
pairs = utils::combn(names(clone_ids), 2, simplify = FALSE)
panel_recomputed = do.call(
  rbind,
  lapply(
    pairs,
    function(pair) {
      data.frame(
        study_1 = pair[[1]],
        study_2 = pair[[2]],
        n_shared = length(
          intersect(
            clone_ids[[pair[[1]]]],
            clone_ids[[pair[[2]]]]
          )
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

panel_key = paste(panel$study_1, panel$study_2, sep = "::")
recomputed_key = paste(
  panel_recomputed$study_1,
  panel_recomputed$study_2,
  sep = "::"
)
panel_recomputed = panel_recomputed[
  match(panel_key, recomputed_key),
  ,
  drop = FALSE
]
panel_error = max(abs(panel$n_shared - panel_recomputed$n_shared))
record_check(
  "segregant_panel_overlap_is_recomputed_exactly",
  nrow(panel) == 21L &&
    panel_error == 0 &&
    sum(panel$n_shared > 0) == 21L,
  paste0(
    "pairs=", nrow(panel),
    "; max shared-count error=", panel_error,
    "; pairs with reuse=", sum(panel$n_shared > 0)
  )
)

joint_summary = summary_table[
  summary_table$set == "868_joint_QTLs",
  ,
  drop = FALSE
]
strict_summary = summary_table[
  summary_table$set == "64_strict_interaction_QTLs",
  ,
  drop = FALSE
]
record_check(
  "reviewer_summary_matches_locus_table",
  nrow(joint_summary) == 1L &&
    nrow(strict_summary) == 1L &&
    joint_summary$original_nominal_p_lt_0_05 ==
      sum(joint$original_interaction_p < 0.05) &&
    joint_summary$robust_nominal_p_lt_0_05 ==
      sum(joint$robust_interaction_p < 0.05) &&
    strict_summary$clone_random_nominal_p_lt_0_05 ==
      sum(strict$clone_random_interaction_p < 0.05),
  "Reviewer-facing summary agrees with locus-level results."
)

analysis_text = paste(
  readLines("iReviewInteractionDependence.R", warn = FALSE),
  collapse = "\n"
)
validator_text = paste(
  readLines("iReviewInteractionDependenceValidate.R", warn = FALSE),
  collapse = "\n"
)
record_check(
  "sensitivity_models_preserve_intended_structure",
  grepl(
    "pheno ~ study + geno + geno:analyte",
    analysis_text,
    fixed = TRUE
  ) &&
    grepl("(1 | clone)", analysis_text, fixed = TRUE) &&
    grepl("REML = FALSE", analysis_text, fixed = TRUE),
  "Robust model uses study fixed intercepts; mixed model adds clone intercept."
)

assignment_token = paste0("<", "-")
assignment_pattern = paste0("(?<!<)", assignment_token)
has_standard_left_assignment = function(text) {
  grepl(assignment_pattern, text, perl = TRUE)
}
record_check(
  "review_dependence_scripts_follow_equals_style",
  !has_standard_left_assignment(analysis_text) &&
    !has_standard_left_assignment(validator_text),
  "No standard left-arrow assignments detected."
)

validation_summary = do.call(rbind, check_state$rows)
write.table(
  validation_summary,
  file = validation_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(validation_summary, row.names = FALSE)

if (!all(validation_summary$passed)) {
  stop(
    "Shared-segregant validation failed: ",
    paste(
      validation_summary$check[!validation_summary$passed],
      collapse = ", "
    )
  )
}

message(
  "Shared-segregant validation passed: ",
  sum(validation_summary$passed),
  "/",
  nrow(validation_summary),
  " checks."
)
