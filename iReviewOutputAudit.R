########################################
########################################
## GENETICS REVIEW FINAL OUTPUT BUNDLE ##
## SIX TABLES USED IN REVIEWER TEXT    ##
########################################
########################################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
review_imed = paste0(imed, "review/")
final_asys = paste0(review_asys, "final/")
build_asys = paste0(review_asys, ".final_build/")

ensure_directory(review_asys)

final_files = c(
  "replicate_filtering.tsv",
  "supported_gene_characteristics.tsv",
  "interaction_GO.tsv",
  "individual_power_and_rediscovery.tsv",
  "power_validation_and_sensitivity.tsv",
  "segregant_dependence.tsv"
)

legacy_final_files = c(
  "replicate_filtering.tsv",
  "supported_gene_characteristics.tsv",
  "interaction_GO.tsv",
  "individual_power_and_rediscovery.tsv",
  "segregant_dependence.tsv"
)

read_tsv = function(path) {
  if (!file.exists(path)) stop("Missing required review output: ", path)
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_build = function(x, filename) {
  write.table(
    x,
    file = paste0(build_asys, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

bundle_matches = function(path, files) {
  if (!dir.exists(path)) return(FALSE)
  present = sort(list.files(path, pattern = "\\.tsv$"))
  identical(present, sort(files))
}

bundle_is_exact = function(path) {
  bundle_matches(path, final_files)
}

make_clean_power_table = function(
  study,
  n,
  pve50,
  pve80,
  observed,
  predicted
) {

  out = data.frame(
    study = as.character(study),
    n = as.numeric(n),
    pve_50pct_power_lod3_03_percent = as.numeric(pve50),
    pve_80pct_power_lod3_03_percent = as.numeric(pve80),
    figure3_observed_albert_query_percent = as.numeric(observed),
    power_predicted_same_direction_lod4_5_percent = as.numeric(predicted),
    observed_minus_predicted_lod4_5_pp = as.numeric(observed) - as.numeric(predicted),
    stringsAsFactors = FALSE
  )

  target =
    out$study != "Albert_2018" &
    is.finite(out$figure3_observed_albert_query_percent) &
    is.finite(out$power_predicted_same_direction_lod4_5_percent)

  if (sum(target) != 9L) {
    stop("Expected nine finite Albert-to-target rediscovery comparisons.")
  }

  mean_row = data.frame(
    study = "MEAN_9_ALBERT_TARGETS",
    n = NA_real_,
    pve_50pct_power_lod3_03_percent = NA_real_,
    pve_80pct_power_lod3_03_percent = NA_real_,
    figure3_observed_albert_query_percent = mean(
      out$figure3_observed_albert_query_percent[target]
    ),
    power_predicted_same_direction_lod4_5_percent = mean(
      out$power_predicted_same_direction_lod4_5_percent[target]
    ),
    observed_minus_predicted_lod4_5_pp = mean(
      out$observed_minus_predicted_lod4_5_pp[target]
    ),
    stringsAsFactors = FALSE
  )

  rbind(out, mean_row)
}

make_power_support_table = function(
  observed,
  predicted,
  scale08,
  scale06,
  albert_public,
  albert_ours,
  max_engine_difference,
  max_simulation_difference
) {

  keep =
    is.finite(observed) &
    is.finite(predicted) &
    is.finite(scale08) &
    is.finite(scale06)

  if (sum(keep) != 9L) {
    stop("Expected nine finite target rows for power sensitivity summary.")
  }

  data.frame(
    metric = c(
      "mean_observed_figure3_rediscovery",
      "mean_predicted_same_direction_lod4_5",
      "mean_observed_minus_predicted_lod4_5",
      "mean_predicted_lod4_5_effect_scale_0_8",
      "mean_predicted_lod4_5_effect_scale_0_6",
      "observed_vs_predicted_pearson_r",
      "observed_vs_predicted_spearman_rho",
      "albert_public_power_at_2_5pct_pve",
      "our_power_at_2_5pct_pve",
      "max_analytic_engine_difference_0_5_to_10pct_pve",
      "max_direct_simulation_difference"
    ),
    value = c(
      mean(observed[keep]),
      mean(predicted[keep]),
      mean(observed[keep] - predicted[keep]),
      mean(scale08[keep]),
      mean(scale06[keep]),
      stats::cor(observed[keep], predicted[keep], method = "pearson"),
      stats::cor(observed[keep], predicted[keep], method = "spearman"),
      albert_public,
      albert_ours,
      max_engine_difference,
      max_simulation_difference
    ),
    unit = c(
      "percent",
      "percent",
      "percentage_points",
      "percent",
      "percent",
      "correlation",
      "correlation",
      "percent",
      "percent",
      "percentage_points",
      "percentage_points"
    ),
    stringsAsFactors = FALSE
  )
}

source_paths = c(
  missingness = paste0(
    review_asys,
    "missingness/missingness_filter_summary.tsv"
  ),
  genes = paste0(
    review_asys,
    "genes/supported_gene_characteristics_summary.tsv"
  ),
  go_key = paste0(review_asys, "go/GO_key_terms.tsv"),
  go_contrast = paste0(
    review_asys,
    "go/GO_discordant_vs_other_interactions.tsv"
  ),
  critical = paste0(
    review_asys,
    "power_effect_curve/critical_r_power_table.tsv"
  ),
  crossstudy = paste0(
    review_asys,
    "power_effect_curve/albert_reference_crossstudy.tsv"
  ),
  benchmark = paste0(
    review_asys,
    "power_effect_curve/albert_2018_power_benchmark.tsv"
  ),
  power_validation = paste0(
    review_asys,
    "power_effect_curve/effect_curve_validation_summary.tsv"
  ),
  figure3 = paste0(
    review_asys,
    "power_effect_curve/figure3_calibration/",
    "albert_figure3_observed_vs_power.tsv"
  ),
  dependence = paste0(
    review_asys,
    "interaction_dependence/segregant_dependence_summary.tsv"
  ),
  panel = paste0(
    review_asys,
    "interaction_dependence/segregant_panel_overlap.tsv"
  )
)

source_ready = all(file.exists(source_paths))
final_ready = bundle_is_exact(final_asys)
legacy_ready = bundle_matches(final_asys, legacy_final_files)

########################################
# BUILD INTO A TEMPORARY DIRECTORY     #
########################################

if (source_ready) {

  if (dir.exists(build_asys)) {
    unlink(build_asys, recursive = TRUE, force = TRUE)
  }
  ensure_directory(build_asys)

  ########################################
  # 1. REPLICATE / MISSINGNESS RESPONSE  #
  ########################################

  replicate_filtering = read_tsv(source_paths[["missingness"]])
  write_build(replicate_filtering, "replicate_filtering.tsv")

  ########################################
  # 2. REPRESENTATIVENESS RESPONSE       #
  ########################################

  supported_gene_characteristics = read_tsv(source_paths[["genes"]])
  write_build(
    supported_gene_characteristics,
    "supported_gene_characteristics.tsv"
  )

  ########################################
  # 3. GO RESPONSE                       #
  ########################################

  go_key = read_tsv(source_paths[["go_key"]])
  go_contrast = read_tsv(source_paths[["go_contrast"]])

  go_key$`Discordant vs other interaction QTLs` = NA_character_
  go_key$`Class contrast FDR` = NA_character_

  contrast_row = go_contrast[
    go_contrast$GO_term == "ribosomal subunit",
    ,
    drop = FALSE
  ]
  key_row = which(go_key$`GO term` == "ribosomal subunit")

  if (nrow(contrast_row) != 1L || length(key_row) != 1L) {
    stop("Could not identify the ribosomal-subunit GO contrast.")
  }

  go_key$`Discordant vs other interaction QTLs`[[key_row]] = paste(
    contrast_row$discordant[[1]],
    "vs",
    contrast_row$other_interactions[[1]]
  )
  go_key$`Class contrast FDR`[[key_row]] = formatC(
    contrast_row$fisher_q[[1]],
    format = "e",
    digits = 2
  )

  write_build(go_key, "interaction_GO.tsv")

  ########################################
  # 4-5. POWER RESPONSE + SUPPORT        #
  ########################################

  critical = read_tsv(source_paths[["critical"]])
  crossstudy = read_tsv(source_paths[["crossstudy"]])
  benchmark = read_tsv(source_paths[["benchmark"]])
  power_validation = read_tsv(source_paths[["power_validation"]])
  figure3 = read_tsv(source_paths[["figure3"]])

  critical_303 = critical[
    abs(critical$lod - 3.03) < 1e-10,
    ,
    drop = FALSE
  ]

  if (nrow(critical_303) != 10L || nrow(figure3) != 9L) {
    stop("Unexpected study counts in final power summary inputs.")
  }

  observed = rep(NA_real_, nrow(critical_303))
  predicted = rep(NA_real_, nrow(critical_303))
  figure_index = match(figure3$target, critical_303$study)

  if (any(is.na(figure_index))) {
    stop("Could not align Figure 3 targets to the study power table.")
  }

  observed[figure_index] = 100 * figure3$observed_figure3_rediscovery
  predicted[figure_index] = 100 * figure3$predicted_same_direction_lod4_5

  power_summary = make_clean_power_table(
    study = critical_303$study,
    n = critical_303$n,
    pve50 = 100 * critical_303$pve_50pct_power,
    pve80 = 100 * critical_303$pve_80pct_power,
    observed = observed,
    predicted = predicted
  )
  write_build(power_summary, "individual_power_and_rediscovery.tsv")

  cross_index = match(figure3$target, crossstudy$target)
  if (any(is.na(cross_index))) {
    stop("Could not align Figure 3 targets to Albert-reference sensitivities.")
  }

  simulation_row = power_validation[
    power_validation$check ==
      "direct_phenotype_simulation_supports_power_engine",
    ,
    drop = FALSE
  ]

  if (nrow(simulation_row) != 1L || !simulation_row$passed[[1]]) {
    stop("Direct phenotype-simulation validation is unavailable or failed.")
  }

  simulation_error = suppressWarnings(
    as.numeric(
      sub(
        ".*maximum absolute error =\\s*([^; ]+).*",
        "\\1",
        simulation_row$detail[[1]]
      )
    )
  )

  if (!is.finite(simulation_error) || nrow(benchmark) != 1L) {
    stop("Could not recover validated power benchmark quantities.")
  }

  power_support = make_power_support_table(
    observed = 100 * figure3$observed_figure3_rediscovery,
    predicted = 100 * figure3$predicted_same_direction_lod4_5,
    scale08 = 100 * crossstudy$expected_detection_lod4_5_scale0_8[cross_index],
    scale06 = 100 * crossstudy$expected_detection_lod4_5_scale0_6[cross_index],
    albert_public = 100 * benchmark$public_code_power_t_test[[1]],
    albert_ours = 100 * benchmark$current_engine_any_detection[[1]],
    max_engine_difference = 100 * benchmark$engine_grid_max_abs_difference[[1]],
    max_simulation_difference = 100 * simulation_error
  )
  write_build(
    power_support,
    "power_validation_and_sensitivity.tsv"
  )

  ########################################
  # 6. SHARED-SEGREGANT RESPONSE         #
  ########################################

  dependence = read_tsv(source_paths[["dependence"]])
  panel = read_tsv(source_paths[["panel"]])

  dependence = dependence[
    dependence$set %in% c(
      "868_joint_QTLs",
      "64_strict_interaction_QTLs"
    ),
    ,
    drop = FALSE
  ]

  dependence_final = dependence[
    ,
    c(
      "set",
      "n_loci",
      "original_nominal_p_lt_0_05",
      "robust_nominal_p_lt_0_05",
      "original_and_robust_nominal_p_lt_0_05",
      "median_robust_to_original_se_ratio",
      "robust_effect_pearson_r",
      "robust_effect_spearman_rho",
      "robust_sign_concordance",
      "clone_random_fit_success",
      "clone_random_full_singular",
      "clone_random_nominal_p_lt_0_05",
      "original_and_clone_random_nominal_p_lt_0_05",
      "clone_random_effect_pearson_r",
      "clone_random_effect_spearman_rho",
      "clone_random_sign_concordance"
    ),
    drop = FALSE
  ]

  dependence_final$n_dataset_pairs = nrow(panel)
  dependence_final$n_dataset_pairs_with_shared_segregants =
    sum(panel$n_shared > 0)

  write_build(dependence_final, "segregant_dependence.tsv")
}

########################################
# MIGRATE THE PREVIOUS FIVE-TABLE BUNDLE
########################################

if (!source_ready && legacy_ready) {

  if (dir.exists(build_asys)) {
    unlink(build_asys, recursive = TRUE, force = TRUE)
  }
  ensure_directory(build_asys)

  for (filename in setdiff(
    legacy_final_files,
    "individual_power_and_rediscovery.tsv"
  )) {
    write_build(
      read_tsv(paste0(final_asys, filename)),
      filename
    )
  }

  old_power = read_tsv(
    paste0(final_asys, "individual_power_and_rediscovery.tsv")
  )

  required_old = c(
    "study",
    "n",
    "pve_50pct_power_lod3_03_percent",
    "pve_80pct_power_lod3_03_percent",
    "figure3_observed_albert_query_percent",
    "power_predicted_same_direction_lod4_5_percent",
    "predicted_lod4_5_scale0_8_percent",
    "predicted_lod4_5_scale0_6_percent",
    "albert_public_power_2_5pct_pve_percent",
    "our_power_2_5pct_pve_percent",
    "max_analytic_engine_difference_0_5_to_10pct_pve_pp",
    "max_direct_simulation_difference_pp"
  )

  if (!all(required_old %in% colnames(old_power))) {
    stop("Legacy final power table lacks columns required for clean migration.")
  }

  study_rows = old_power[old_power$study != "MEAN_9_ALBERT_TARGETS", , drop = FALSE]

  clean_power = make_clean_power_table(
    study = study_rows$study,
    n = study_rows$n,
    pve50 = study_rows$pve_50pct_power_lod3_03_percent,
    pve80 = study_rows$pve_80pct_power_lod3_03_percent,
    observed = study_rows$figure3_observed_albert_query_percent,
    predicted = study_rows$power_predicted_same_direction_lod4_5_percent
  )
  write_build(clean_power, "individual_power_and_rediscovery.tsv")

  target =
    study_rows$study != "Albert_2018" &
    is.finite(study_rows$figure3_observed_albert_query_percent) &
    is.finite(study_rows$power_predicted_same_direction_lod4_5_percent)
  albert = study_rows[study_rows$study == "Albert_2018", , drop = FALSE]

  if (sum(target) != 9L || nrow(albert) != 1L) {
    stop("Legacy power table does not contain expected target and Albert rows.")
  }

  power_support = make_power_support_table(
    observed = study_rows$figure3_observed_albert_query_percent[target],
    predicted = study_rows$power_predicted_same_direction_lod4_5_percent[target],
    scale08 = study_rows$predicted_lod4_5_scale0_8_percent[target],
    scale06 = study_rows$predicted_lod4_5_scale0_6_percent[target],
    albert_public = albert$albert_public_power_2_5pct_pve_percent[[1]],
    albert_ours = albert$our_power_2_5pct_pve_percent[[1]],
    max_engine_difference =
      albert$max_analytic_engine_difference_0_5_to_10pct_pve_pp[[1]],
    max_simulation_difference =
      albert$max_direct_simulation_difference_pp[[1]]
  )
  write_build(
    power_support,
    "power_validation_and_sensitivity.tsv"
  )
}

########################################
# VALIDATE BEFORE DELETING ANY SOURCE  #
########################################

if (dir.exists(build_asys)) {

  if (!bundle_is_exact(build_asys)) {
    stop(
      "Temporary final reviewer bundle is not exactly six tables: ",
      paste(list.files(build_asys), collapse = ", ")
    )
  }

  for (filename in final_files) {
    tab = read_tsv(paste0(build_asys, filename))
    if (nrow(tab) < 1L) {
      stop("Final reviewer table is empty: ", filename)
    }
  }

  clean_power_check = read_tsv(
    paste0(build_asys, "individual_power_and_rediscovery.tsv")
  )
  if (
    nrow(clean_power_check) != 11L ||
    ncol(clean_power_check) != 7L ||
    sum(clean_power_check$study == "Albert_2018") != 1L ||
    sum(clean_power_check$study == "MEAN_9_ALBERT_TARGETS") != 1L
  ) {
    stop("Clean reviewer power table has unexpected shape or study rows.")
  }

  support_check = read_tsv(
    paste0(build_asys, "power_validation_and_sensitivity.tsv")
  )
  if (nrow(support_check) != 11L || ncol(support_check) != 3L) {
    stop("Power validation/sensitivity table has unexpected shape.")
  }

  if (dir.exists(final_asys)) {
    unlink(final_asys, recursive = TRUE, force = TRUE)
  }

  if (!file.rename(build_asys, final_asys)) {
    stop("Could not promote validated temporary reviewer bundle to final/.")
  }

  final_ready = TRUE
}

if (!final_ready || !bundle_is_exact(final_asys)) {
  stop(
    "A validated six-table final bundle is unavailable. Recreate the ",
    "reviewer analysis outputs, then rerun this finalizer."
  )
}

########################################
# REMOVE ALL NONFINAL REVIEW ARTIFACTS #
########################################

review_entries = list.files(
  review_asys,
  full.names = TRUE,
  recursive = FALSE,
  all.files = FALSE
)

keep_name = basename(sub("/$", "", final_asys))
stale = review_entries[basename(review_entries) != keep_name]

if (length(stale) > 0L) {
  unlink(stale, recursive = TRUE, force = TRUE)
}

if (dir.exists(review_imed)) {
  unlink(review_imed, recursive = TRUE, force = TRUE)
}

########################################
# FINAL MANIFEST                       #
########################################

manifest = data.frame(
  table = final_files,
  purpose = c(
    "Replicate averaging and missingness filtering",
    "Representativeness of the 385 supported genes",
    "GO enrichment of interaction QTLs",
    "Study-level power and Figure 3 rediscovery",
    "Power validation, calibration, and attenuation sensitivity",
    "Dependence from shared segregants"
  ),
  stringsAsFactors = FALSE
)

write.table(
  manifest,
  file = paste0(review_asys, "review_output_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

if (!bundle_is_exact(final_asys)) {
  stop("Final six-table bundle changed unexpectedly during cleanup.")
}

cat("\nFINAL REVIEWER OUTPUT BUNDLE\n")
cat("Substantive tables: 6\n")
print(manifest, row.names = FALSE)
cat("\nDirectory: ", final_asys, "\n", sep = "")
cat("Power results are split into a clean study table and a compact support table.\n")
cat("All nonfinal review outputs and intermediates removed.\n")
