################################################
################################################
## GENETICS REVIEW POWER ANALYSIS VALIDATOR   ##
## ALL ORDERED DATASET PAIRS                  ##
################################################
################################################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
power_asys = paste0(review_asys, "power_effect_curve/")
allpair_asys = paste0(power_asys, "all_pairs/")
validation_path = paste0(allpair_asys, "all_pairs_validation_summary.tsv")

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

add_check = function(name, passed, detail) {
  data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

power = read_tsv(paste0(allpair_asys, "pairwise_power_for_frank.tsv"))
comparison = read_tsv(paste0(allpair_asys, "pairwise_observed_vs_power.tsv"))
fit = read_tsv(paste0(allpair_asys, "pairwise_power_fit_summary.tsv"))
legacy_albert = read_tsv(paste0(power_asys, "albert_reference_crossstudy.tsv"))
legacy_figure3 = read_tsv(
  paste0(
    power_asys,
    "figure3_calibration/albert_figure3_observed_vs_power.tsv"
  )
)

frank_name_map = c(
  Brem_2002 = "Brem_2002",
  Brem_2005 = "Brem_2005",
  SmithEth_2008 = "Smith_2008_Ethanol",
  SmithGlc_2008 = "Smith_2008_Glucose",
  Khan_2009 = "Khan_2009",
  Foss_2011 = "Foss_2011",
  Picotti_2013 = "Picotti_2013",
  Albert_2018 = "Albert_2018",
  GrossbachRNA_2022 = "Grossbach_2022_mRNA",
  GrossbachPRT_2022 = "Grossbach_2022_protein"
)
expected_names = unname(frank_name_map)

checks = list()
check_i = 1L

checks[[check_i]] = add_check(
  "exact_90_directed_pairs",
  nrow(power) == 90L &&
    nrow(comparison) == 90L &&
    !anyDuplicated(power[, c("fromData", "toData")]) &&
    !anyDuplicated(comparison[, c("fromData", "toData")]) &&
    all(power$fromData != power$toData),
  paste(
    "power rows =", nrow(power),
    "; comparison rows =", nrow(comparison)
  )
)
check_i = check_i + 1L

checks[[check_i]] = add_check(
  "uses_exact_frank_dataset_names",
  setequal(unique(power$fromData), expected_names) &&
    setequal(unique(power$toData), expected_names),
  paste("datasets =", paste(sort(unique(power$fromData)), collapse = ", "))
)
check_i = check_i + 1L

required_join_columns = c(
  "fromData",
  "toData",
  "fromData_analyte",
  "toData_analyte",
  "fromData_medium",
  "toData_medium",
  "analyteCombo",
  "mediumSame",
  "fromData_n",
  "toData_n",
  "nSourceQTLsUsedForPower",
  "meanAbsSourceEffect",
  "medianAbsSourceEffect",
  "expectedRediscoveryLOD2",
  "expectedRediscoveryLOD3_03",
  "expectedRediscoveryLOD4_5",
  "expectedRediscoveryLOD7"
)
checks[[check_i]] = add_check(
  "merge_ready_columns_present",
  all(required_join_columns %in% names(power)),
  paste(
    "required columns present =",
    sum(required_join_columns %in% names(power)),
    "/", length(required_join_columns)
  )
)
check_i = check_i + 1L

prediction_columns = c(
  "expectedRediscoveryLOD2",
  "expectedRediscoveryLOD3_03",
  "expectedRediscoveryLOD4_5",
  "expectedRediscoveryLOD7"
)
prediction_matrix = as.matrix(power[, prediction_columns])
checks[[check_i]] = add_check(
  "predictions_are_finite_probabilities",
  all(is.finite(prediction_matrix)) &&
    all(prediction_matrix >= 0) &&
    all(prediction_matrix <= 1),
  paste(
    "range =",
    paste(signif(range(prediction_matrix), 6), collapse = " to ")
  )
)
check_i = check_i + 1L

monotone =
  power$expectedRediscoveryLOD2 >= power$expectedRediscoveryLOD3_03 &
  power$expectedRediscoveryLOD3_03 >= power$expectedRediscoveryLOD4_5 &
  power$expectedRediscoveryLOD4_5 >= power$expectedRediscoveryLOD7
checks[[check_i]] = add_check(
  "higher_lod_never_increases_predicted_rediscovery",
  all(monotone),
  paste("passing pairs =", sum(monotone), "/", length(monotone))
)
check_i = check_i + 1L

checks[[check_i]] = add_check(
  "observed_figure3_counts_are_valid",
  all(comparison$numberOfQTLs > 0L) &&
    all(comparison$overlapsWithStrand >= 0L) &&
    all(comparison$overlapsWithStrand <= comparison$numberOfQTLs) &&
    max(abs(
      comparison$fractionOverlapsWithStrand -
        comparison$overlapsWithStrand / comparison$numberOfQTLs
    )) < 1e-14,
  paste(
    "minimum denominator =", min(comparison$numberOfQTLs),
    "; maximum denominator =", max(comparison$numberOfQTLs)
  )
)
check_i = check_i + 1L

albert_power = power[power$fromData == "Albert_2018", , drop = FALSE]
albert_target_internal = names(frank_name_map)[
  match(albert_power$toData, unname(frank_name_map))
]
legacy_power_index = match(albert_target_internal, legacy_albert$target)
power_error = max(abs(
  albert_power$expectedRediscoveryLOD4_5 -
    legacy_albert$expected_same_direction_lod4_5[legacy_power_index]
))
checks[[check_i]] = add_check(
  "albert_lod4_5_matches_existing_analysis",
  nrow(albert_power) == 9L &&
    !any(is.na(legacy_power_index)) &&
    power_error < 1e-12,
  paste("maximum difference =", signif(power_error, 7))
)
check_i = check_i + 1L

albert_observed = comparison[
  comparison$fromData == "Albert_2018",
  ,
  drop = FALSE
]
legacy_observed_index = match(
  names(frank_name_map)[
    match(albert_observed$toData, unname(frank_name_map))
  ],
  legacy_figure3$target
)
observed_error = max(abs(
  albert_observed$fractionOverlapsWithStrand -
    legacy_figure3$observed_figure3_rediscovery[legacy_observed_index]
))
count_same = all(
  albert_observed$numberOfQTLs ==
    legacy_figure3$n_figure3_query_eligible[legacy_observed_index]
) && all(
  albert_observed$overlapsWithStrand ==
    legacy_figure3$n_figure3_matches[legacy_observed_index]
)
checks[[check_i]] = add_check(
  "albert_observed_rows_match_existing_figure3_reconstruction",
  nrow(albert_observed) == 9L &&
    !any(is.na(legacy_observed_index)) &&
    count_same &&
    observed_error < 1e-14,
  paste(
    "counts identical =", count_same,
    "; maximum fraction difference =", signif(observed_error, 7)
  )
)
check_i = check_i + 1L

checks[[check_i]] = add_check(
  "fit_summary_has_four_lod_rows",
  nrow(fit) == 4L &&
    setequal(fit$lod, c(2, 3.03, 4.5, 7)) &&
    all(fit$n_pairs == 90L) &&
    all(is.finite(fit$pearson_r)) &&
    all(is.finite(fit$spearman_rho)),
  paste("LOD rows =", paste(fit$lod, collapse = ", "))
)
check_i = check_i + 1L

validation_summary = do.call(rbind, checks)
rownames(validation_summary) = NULL

write.table(
  validation_summary,
  file = validation_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (!all(validation_summary$passed)) {
  failed = validation_summary$check[!validation_summary$passed]
  stop(
    "All-pair power validation failed: ",
    paste(failed, collapse = ", ")
  )
}

cat("\nALL-PAIR POWER VALIDATION PASSED\n")
cat("Checks: ", nrow(validation_summary), "\n", sep = "")
cat("Directed pairs: 90\n")
cat("Merge keys: fromData + toData\n")
cat("Albert LOD 4.5 reproduces existing analysis\n")
