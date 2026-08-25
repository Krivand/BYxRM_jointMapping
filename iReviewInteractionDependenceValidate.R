########################################
########################################
# VALIDATE SEGREGANT DEPENDENCE       #
########################################
########################################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
dependence_asys = paste0(review_asys, "interaction_dependence/")
ensure_directory(dependence_asys)

analysis_path = "iReviewInteractionDependence.R"
validator_path = "iReviewInteractionDependenceValidate.R"

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

load_single_rda = function(path, object_name = NULL) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  env = new.env(parent = emptyenv())
  loaded = load(path, envir = env)
  if (!is.null(object_name)) {
    if (!(object_name %in% loaded)) stop("Missing object ", object_name, " in ", path)
    return(env[[object_name]])
  }
  if (length(loaded) != 1L) stop("Expected one object in ", path)
  env[[loaded[[1]]]]
}

canonical_clone = function(x) {
  x = trimws(as.character(x))
  x = gsub("_", ".", x, fixed = TRUE)
  gsub("\\s+", "", x)
}

interaction_term = function(beta_names) {
  term = grep(
    "(^geno:analyte|^analyte.*:geno$)",
    beta_names,
    value = TRUE
  )
  if (length(term) != 1L) stop("Could not identify interaction term.")
  term[[1]]
}

cluster_vcov_cr1 = function(fit, cluster) {
  X = stats::model.matrix(fit)
  residual = stats::residuals(fit)
  cluster = factor(cluster)
  N = nrow(X)
  K = ncol(X)
  G = nlevels(cluster)
  bread = solve(crossprod(X))
  score_cluster = rowsum(X * residual, group = cluster, reorder = FALSE)
  correction = (G / (G - 1)) * ((N - 1) / (N - K))
  list(
    vcov = correction * bread %*% crossprod(score_cluster) %*% bread,
    df = G - 1L
  )
}

fit_clone_random_lrt = function(model_frame, clone) {
  dat = model_frame
  dat$clone = factor(clone)
  dat$study = factor(dat$study)
  dat$analyte = factor(dat$analyte)
  ctrl = lme4::lmerControl(
    check.conv.singular = "ignore",
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
  m0 = suppressWarnings(
    lme4::lmer(
      pheno ~ (1 | study) + analyte + geno + (1 | clone),
      data = dat,
      control = ctrl,
      REML = FALSE
    )
  )
  m1 = suppressWarnings(
    lme4::lmer(
      pheno ~ (1 | study) + analyte + geno + geno:analyte + (1 | clone),
      data = dat,
      control = ctrl,
      REML = FALSE
    )
  )
  term = interaction_term(names(lme4::fixef(m1)))
  lrt = stats::anova(m0, m1)
  vc = as.data.frame(lme4::VarCorr(m1))
  clone_variance = vc$vcov[vc$grp == "clone" & vc$var1 == "(Intercept)"]
  data.frame(
    effect = unname(lme4::fixef(m1)[term]),
    p = as.numeric(lrt[2, "Pr(>Chisq)"]),
    lod = (
      as.numeric(stats::logLik(m1)) -
      as.numeric(stats::logLik(m0))
    ) / log(10),
    clone_variance = if (length(clone_variance) == 1L) clone_variance else NA_real_,
    full_singular = lme4::isSingular(m1, tol = 1e-5),
    stringsAsFactors = FALSE
  )
}

checks = list()
check_i = 1L
add_check = function(name, passed, detail) {
  data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}
record_check = function(name, passed, detail) {
  checks[[check_i]] <<- add_check(name, passed, detail)
  check_i <<- check_i + 1L
}

################
# READ OUTPUTS #
################

panel = read_tsv(paste0(dependence_asys, "segregant_panel_overlap.tsv"))
locus = read_tsv(paste0(dependence_asys, "segregant_dependence_locus_results.tsv"))
summary_table = read_tsv(paste0(dependence_asys, "segregant_dependence_summary.tsv"))

jBH_analyte_geno = load_single_rda(
  paste0(imed, "jBH_analyte_geno.rda"),
  "jBH_analyte_geno"
)
jBH_analyte_ixn = load_single_rda(
  paste0(imed, "jBH_analyte_ixn.rda"),
  "jBH_analyte_ixn"
)

joint_keys = paste(jBH_analyte_geno$gene, jBH_analyte_geno$marker, sep = "::")
strict_keys = paste(jBH_analyte_ixn$gene, jBH_analyte_ixn$marker, sep = "::")
target_keys = unique(c(joint_keys, strict_keys))

analysis_outputs = sort(
  list.files(
    dependence_asys,
    pattern = "^segregant_.*\\.tsv$",
    full.names = FALSE
  )
)
expected_analysis_outputs = sort(
  c(
    "segregant_panel_overlap.tsv",
    "segregant_dependence_locus_results.tsv",
    "segregant_dependence_summary.tsv"
  )
)
record_check(
  "analysis_emits_three_primary_tables",
  identical(analysis_outputs, expected_analysis_outputs),
  paste(analysis_outputs, collapse = ", ")
)

record_check(
  "targeted_locus_keys_are_exact",
  nrow(locus) == length(target_keys) &&
    !anyDuplicated(locus$key) &&
    setequal(locus$key, target_keys),
  paste0("rows=", nrow(locus), "; expected=", length(target_keys))
)

record_check(
  "joint_and_strict_membership_counts_are_exact",
  sum(locus$in_joint_868) == 868L && sum(locus$in_strict_64) == 64L,
  paste0(
    "joint=", sum(locus$in_joint_868),
    "; strict=", sum(locus$in_strict_64)
  )
)

record_check(
  "original_joint_posthoc_count_is_reproduced",
  sum(locus$in_joint_868 & locus$original_interaction_p < 0.05) == 149L,
  paste0(
    "nominal joint interactions=",
    sum(locus$in_joint_868 & locus$original_interaction_p < 0.05)
  )
)

joint_error = max(
  abs(
    locus$original_interaction_p[locus$in_joint_868] -
    locus$stored_joint_interaction_p[locus$in_joint_868]
  ),
  na.rm = TRUE
)
strict_error = max(
  abs(
    locus$original_interaction_p[locus$in_strict_64] -
    locus$stored_strict_interaction_p[locus$in_strict_64]
  ),
  na.rm = TRUE
)
record_check(
  "original_targeted_refits_reproduce_stored_models",
  joint_error < 1e-8 && strict_error < 1e-8,
  paste0("joint max p error=", joint_error, "; strict=", strict_error)
)

alignment_error = max(locus$maximum_genotype_alignment_difference, na.rm = TRUE)
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

random_success = sum(locus$clone_random_fit_success %in% TRUE)
record_check(
  "clone_random_effect_refits_are_attempted_and_mostly_estimable",
  random_success >= floor(0.90 * nrow(locus)),
  paste0(random_success, "/", nrow(locus), " successful")
)

########################
# PANEL OVERLAP CHECK   #
########################

metadata_analyte = scaMetadata
studies_analyte = lapply(
  metadata_analyte$study,
  function(study) load_single_rda(paste0(imed, study, ".rda"), study)
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
      x = clone_ids[[pair[[1]]]]
      y = clone_ids[[pair[[2]]]]
      data.frame(
        study_1 = pair[[1]],
        study_2 = pair[[2]],
        n_study_1 = length(x),
        n_study_2 = length(y),
        n_shared = length(intersect(x, y)),
        stringsAsFactors = FALSE
      )
    }
  )
)

panel_key = paste(panel$study_1, panel$study_2, sep = "::")
recomputed_key = paste(panel_recomputed$study_1, panel_recomputed$study_2, sep = "::")
panel_recomputed = panel_recomputed[match(panel_key, recomputed_key), , drop = FALSE]
panel_error = max(abs(panel$n_shared - panel_recomputed$n_shared))
record_check(
  "segregant_panel_overlap_is_recomputed_exactly",
  nrow(panel) == length(pairs) && panel_error == 0 && all(panel$n_shared > 0),
  paste0("pairs=", nrow(panel), "; max shared-count error=", panel_error)
)

########################################
# INDEPENDENT LOCUS SUBSET RECOMPUTE  #
########################################

metadata_analyte$genes = lapply(
  metadata_analyte$study,
  function(study) colnames(studies_analyte[[study]]$pheno)
)
metadata_analyte$segs = lapply(
  metadata_analyte$study,
  function(study) rownames(studies_analyte[[study]]$pheno)
)
supported_analyte = suppTwo(metadata_analyte)
unionBG_analyte = load_single_rda(
  paste0(imed, "unionAnalyteBG.rda"),
  "unionAnalyteBG"
)

strict_candidates = locus[
  locus$in_strict_64 & locus$clone_random_fit_success,
  ,
  drop = FALSE
]
joint_sig_candidates = locus[
  locus$in_joint_868 & !locus$in_strict_64 &
    locus$original_interaction_p < 0.05 &
    locus$clone_random_fit_success,
  ,
  drop = FALSE
]
joint_null_candidates = locus[
  locus$in_joint_868 & !locus$in_strict_64 &
    locus$original_interaction_p >= 0.05 &
    locus$clone_random_fit_success,
  ,
  drop = FALSE
]

choose_three = function(x) {
  if (nrow(x) == 0L) return(x)
  x = x[order(x$key), , drop = FALSE]
  x[seq_len(min(3L, nrow(x))), , drop = FALSE]
}

subset_rows = unique(
  rbind(
    choose_three(strict_candidates),
    choose_three(joint_sig_candidates),
    choose_three(joint_null_candidates)
  )
)
subset_pairs = subset_rows[, c("gene", "marker"), drop = FALSE]

subset_models = drillQTL(
  studies = studies_analyte,
  dfQTL = subset_pairs,
  unionBg = unionBG_analyte,
  smetadata = metadata_analyte,
  suppGenes = supported_analyte,
  fixedEffect = "analyte",
  correctBG = TRUE,
  correctGRM = FALSE,
  cores = min(8L, nrow(subset_pairs)),
  chromosomes = NULL
)

reconstruct_clone = function(gene, marker, model_frame) {
  support = metadata_analyte$study[
    vapply(metadata_analyte$genes, function(x) gene %in% x, logical(1))
  ]
  clone = unlist(
    lapply(
      support,
      function(study) canonical_clone(rownames(studies_analyte[[study]]$gdata))
    ),
    use.names = FALSE
  )
  study_expected = unlist(
    lapply(
      support,
      function(study) rep(study, nrow(studies_analyte[[study]]$gdata))
    ),
    use.names = FALSE
  )
  if (!all(as.character(model_frame$study) == study_expected)) {
    stop("Validator study alignment failed.")
  }
  clone
}

subset_recomputed = do.call(
  rbind,
  lapply(
    subset_models,
    function(x) {
      mf = stats::model.frame(x$interaction)
      clone = reconstruct_clone(as.character(x$gene), as.character(x$marker), mf)

      robust_fit = stats::lm(
        pheno ~ study + geno + geno:analyte,
        data = mf
      )
      term = interaction_term(names(stats::coef(robust_fit)))
      robust = cluster_vcov_cr1(robust_fit, clone)
      robust_se = sqrt(robust$vcov[term, term])
      robust_effect = unname(stats::coef(robust_fit)[term])
      robust_p = 2 * stats::pt(
        -abs(robust_effect / robust_se),
        df = robust$df
      )

      random = fit_clone_random_lrt(mf, clone)

      data.frame(
        key = paste(x$gene, x$marker, sep = "::"),
        robust_effect = robust_effect,
        robust_se = robust_se,
        robust_p = robust_p,
        random_effect = random$effect,
        random_p = random$p,
        random_lod = random$lod,
        random_variance = random$clone_variance,
        stringsAsFactors = FALSE
      )
    }
  )
)

subset_output = locus[match(subset_recomputed$key, locus$key), , drop = FALSE]
robust_subset_error = max(
  abs(
    c(
      subset_recomputed$robust_effect - subset_output$robust_interaction_effect,
      subset_recomputed$robust_se - subset_output$robust_interaction_se,
      subset_recomputed$robust_p - subset_output$robust_interaction_p
    )
  ),
  na.rm = TRUE
)
record_check(
  "independent_subset_reproduces_cluster_robust_results",
  robust_subset_error < 1e-8,
  paste0("n=", nrow(subset_recomputed), "; max error=", robust_subset_error)
)

random_subset_error = max(
  abs(
    c(
      subset_recomputed$random_effect - subset_output$clone_random_interaction_effect,
      subset_recomputed$random_p - subset_output$clone_random_interaction_p,
      subset_recomputed$random_lod - subset_output$clone_random_lod,
      subset_recomputed$random_variance - subset_output$clone_random_variance
    )
  ),
  na.rm = TRUE
)
record_check(
  "independent_subset_reproduces_clone_random_refits",
  random_subset_error < 1e-7,
  paste0("n=", nrow(subset_recomputed), "; max error=", random_subset_error)
)

########################
# SUMMARY CONSISTENCY  #
########################

joint_summary = summary_table[summary_table$set == "868_joint_QTLs", , drop = FALSE]
strict_summary = summary_table[
  summary_table$set == "64_strict_interaction_QTLs",
  ,
  drop = FALSE
]
summary_consistent =
  nrow(summary_table) == 3L &&
  nrow(joint_summary) == 1L &&
  nrow(strict_summary) == 1L &&
  joint_summary$n_loci == 868L &&
  strict_summary$n_loci == 64L &&
  joint_summary$original_nominal_p_lt_0_05 == 149L &&
  joint_summary$robust_nominal_p_lt_0_05 ==
    sum(locus$in_joint_868 & locus$robust_interaction_p < 0.05) &&
  strict_summary$clone_random_fit_success ==
    sum(locus$in_strict_64 & locus$clone_random_fit_success)
record_check(
  "reviewer_summary_matches_locus_table",
  summary_consistent,
  "Counts in the three-row summary agree with the locus-level output."
)

########################
# MODEL AND STYLE AUDIT #
########################

analysis_text = paste(readLines(analysis_path, warn = FALSE), collapse = "\n")
validator_text = paste(readLines(validator_path, warn = FALSE), collapse = "\n")
record_check(
  "sensitivity_models_preserve_intended_structure",
  grepl("pheno ~ study + geno + geno:analyte", analysis_text, fixed = TRUE) &&
    grepl("(1 | clone)", analysis_text, fixed = TRUE) &&
    grepl("REML = FALSE", analysis_text, fixed = TRUE),
  "Robust model uses dataset fixed intercepts; mixed sensitivity adds a clone random intercept."
)

assignment_token = paste0("<", "-")
record_check(
  "review_dependence_scripts_follow_equals_style",
  !grepl(assignment_token, analysis_text, fixed = TRUE) &&
    !grepl(assignment_token, validator_text, fixed = TRUE),
  "No left-arrow assignments detected."
)

########################
# WRITE VALIDATION     #
########################

validation_summary = do.call(rbind, checks)
write.table(
  validation_summary,
  file = paste0(dependence_asys, "segregant_dependence_validation_summary.tsv"),
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
