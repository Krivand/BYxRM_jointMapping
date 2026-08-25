########################################
########################################
## GENETICS REVIEW FINAL OUTPUT BUNDLE ##
## FIVE TABLES USED IN REVIEWER TEXT   ##
########################################
########################################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
review_imed = paste0(imed, "review/")
final_asys = paste0(review_asys, "final/")

ensure_directory(review_asys)
ensure_directory(final_asys)

final_files = c(
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

write_final = function(x, filename) {
  write.table(
    x,
    file = paste0(final_asys, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

source_paths = c(
  missingness = paste0(review_asys, "missingness/missingness_filter_summary.tsv"),
  genes = paste0(review_asys, "genes/supported_gene_characteristics_summary.tsv"),
  go_key = paste0(review_asys, "go/GO_key_terms.tsv"),
  go_contrast = paste0(review_asys, "go/GO_discordant_vs_other_interactions.tsv"),
  critical = paste0(review_asys, "power_effect_curve/critical_r_power_table.tsv"),
  crossstudy = paste0(review_asys, "power_effect_curve/albert_reference_crossstudy.tsv"),
  benchmark = paste0(review_asys, "power_effect_curve/albert_2018_power_benchmark.tsv"),
  power_validation = paste0(review_asys, "power_effect_curve/effect_curve_validation_summary.tsv"),
  figure3 = paste0(
    review_asys,
    "power_effect_curve/figure3_calibration/albert_figure3_observed_vs_power.tsv"
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
final_ready = all(file.exists(paste0(final_asys, final_files)))

if (source_ready) {

  ########################################
  # 1. REPLICATE / MISSINGNESS RESPONSE  #
  ########################################

  replicate_filtering = read_tsv(source_paths[["missingness"]])
  write_final(replicate_filtering, "replicate_filtering.tsv")

  ########################################
  # 2. REPRESENTATIVENESS RESPONSE       #
  ########################################

  supported_gene_characteristics = read_tsv(source_paths[["genes"]])
  write_final(
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

  write_final(go_key, "interaction_GO.tsv")

  ########################################
  # 4. INDIVIDUAL-QTL POWER RESPONSE     #
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

  power_summary = data.frame(
    study = critical_303$study,
    analyte = critical_303$analyte,
    medium = critical_303$medium,
    n = critical_303$n,
    pve_50pct_power_lod3_03_percent = 100 * critical_303$pve_50pct_power,
    pve_80pct_power_lod3_03_percent = 100 * critical_303$pve_80pct_power,
    figure3_observed_albert_query_percent = NA_real_,
    power_predicted_same_direction_lod4_5_percent = NA_real_,
    observed_minus_predicted_lod4_5_pp = NA_real_,
    predicted_lod4_5_scale0_8_percent = NA_real_,
    predicted_lod4_5_scale0_6_percent = NA_real_,
    albert_public_power_2_5pct_pve_percent = NA_real_,
    our_power_2_5pct_pve_percent = NA_real_,
    max_analytic_engine_difference_0_5_to_10pct_pve_pp = NA_real_,
    max_direct_simulation_difference_pp = NA_real_,
    stringsAsFactors = FALSE
  )

  figure_index = match(figure3$target, power_summary$study)
  if (any(is.na(figure_index))) {
    stop("Could not align Figure 3 targets to the study power table.")
  }

  power_summary$figure3_observed_albert_query_percent[figure_index] =
    100 * figure3$observed_figure3_rediscovery
  power_summary$power_predicted_same_direction_lod4_5_percent[figure_index] =
    100 * figure3$predicted_same_direction_lod4_5
  power_summary$observed_minus_predicted_lod4_5_pp[figure_index] =
    figure3$difference_lod4_5_percentage_points

  cross_index = match(crossstudy$target, power_summary$study)
  if (any(is.na(cross_index))) {
    stop("Could not align Albert-reference targets to the study power table.")
  }
  power_summary$predicted_lod4_5_scale0_8_percent[cross_index] =
    100 * crossstudy$expected_detection_lod4_5_scale0_8
  power_summary$predicted_lod4_5_scale0_6_percent[cross_index] =
    100 * crossstudy$expected_detection_lod4_5_scale0_6

  albert_i = match("Albert_2018", power_summary$study)
  if (is.na(albert_i) || nrow(benchmark) != 1L) {
    stop("Albert benchmark row is not uniquely defined.")
  }

  power_summary$albert_public_power_2_5pct_pve_percent[[albert_i]] =
    100 * benchmark$public_code_power_t_test[[1]]
  power_summary$our_power_2_5pct_pve_percent[[albert_i]] =
    100 * benchmark$current_engine_any_detection[[1]]
  power_summary$max_analytic_engine_difference_0_5_to_10pct_pve_pp[[albert_i]] =
    100 * benchmark$engine_grid_max_abs_difference[[1]]

  simulation_row = power_validation[
    power_validation$check == "direct_phenotype_simulation_supports_power_engine",
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
  if (!is.finite(simulation_error)) {
    stop("Could not parse the direct-simulation validation error.")
  }
  power_summary$max_direct_simulation_difference_pp[[albert_i]] =
    100 * simulation_error

  mean_row = data.frame(
    study = "MEAN_9_ALBERT_TARGETS",
    analyte = "",
    medium = "",
    n = NA_real_,
    pve_50pct_power_lod3_03_percent = NA_real_,
    pve_80pct_power_lod3_03_percent = NA_real_,
    figure3_observed_albert_query_percent =
      mean(100 * figure3$observed_figure3_rediscovery),
    power_predicted_same_direction_lod4_5_percent =
      mean(100 * figure3$predicted_same_direction_lod4_5),
    observed_minus_predicted_lod4_5_pp =
      mean(figure3$difference_lod4_5_percentage_points),
    predicted_lod4_5_scale0_8_percent =
      mean(100 * crossstudy$expected_detection_lod4_5_scale0_8),
    predicted_lod4_5_scale0_6_percent =
      mean(100 * crossstudy$expected_detection_lod4_5_scale0_6),
    albert_public_power_2_5pct_pve_percent = NA_real_,
    our_power_2_5pct_pve_percent = NA_real_,
    max_analytic_engine_difference_0_5_to_10pct_pve_pp = NA_real_,
    max_direct_simulation_difference_pp = NA_real_,
    stringsAsFactors = FALSE
  )

  power_summary = rbind(power_summary, mean_row)
  write_final(
    power_summary,
    "individual_power_and_rediscovery.tsv"
  )

  ########################################
  # 5. SHARED-SEGREGANT RESPONSE         #
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

  write_final(dependence_final, "segregant_dependence.tsv")

  final_ready = TRUE
}

if (!final_ready) {
  stop(
    "Neither the complete source-output set nor an existing five-table final ",
    "bundle is available. Recreate the relevant reviewer analyses first."
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

keep_paths = c(
  sub("/$", "", final_asys),
  paste0(review_asys, "review_output_manifest.tsv")
)

stale = review_entries[!(review_entries %in% keep_paths)]
if (length(stale) > 0L) {
  unlink(stale, recursive = TRUE, force = TRUE)
}

if (dir.exists(review_imed)) {
  unlink(review_imed, recursive = TRUE, force = TRUE)
}

########################################
# VALIDATE THE FINAL FIVE-TABLE BUNDLE #
########################################

present = sort(list.files(final_asys, pattern = "\\.tsv$"))
if (!identical(present, sort(final_files))) {
  stop(
    "Final reviewer bundle is not exactly five tables: ",
    paste(present, collapse = ", ")
  )
}

manifest = data.frame(
  table = final_files,
  purpose = c(
    "Reviewer response on replicate averaging and missingness filtering",
    "Reviewer response on representativeness of the 385 supported genes",
    "Reviewer response on GO enrichment of interaction QTLs",
    "Reviewer response on individual-study power and Figure 3 rediscovery",
    "Reviewer response on dependence from shared segregants"
  ),
  stringsAsFactors = FALSE
)

write.table(
  manifest,
  file = paste0(review_asys, "review_output_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nFINAL REVIEWER OUTPUT BUNDLE\n")
cat("Substantive tables: 5\n")
print(manifest, row.names = FALSE)
cat("\nDirectory: ", final_asys, "\n", sep = "")
cat("All superseded, validation-only, and intermediate review outputs removed.\n")
