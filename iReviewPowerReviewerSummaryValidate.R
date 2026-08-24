############################################
############################################
## GENETICS REVIEW POWER SUMMARY VALIDATOR##
## REVIEWER-FACING CALIBRATION TABLES     ##
############################################
############################################

#################
# PRELIMINARIES #
#################

wd = "/PHShome/kv1026/p05/"
out_dir = paste0(
  wd,
  "asys/review/power_pair_conditional/"
)
validation_out = paste0(
  out_dir,
  "reviewer_summary_validation.tsv"
)

required_files = c(
  "study_design.tsv",
  "observed_individual_rediscovery.tsv",
  "pair_conditional_power_sweep.tsv",
  "pair_conditional_validation_summary.tsv",
  "analytic_power_validation.tsv",
  "chance_overlap_null.tsv",
  "reviewer_power_summary.tsv",
  "reviewer_power_table.tsv",
  "reviewer_reciprocal_asymmetry.tsv",
  "reviewer_reciprocal_asymmetry_summary.tsv",
  "reviewer_reciprocal_asymmetry_key.tsv",
  "reviewer_supporting_validation.tsv"
)
missing_files = required_files[
  !file.exists(paste0(out_dir, required_files))
]
if (length(missing_files) > 0L) {
  stop(
    "Missing reviewer-summary outputs: ",
    paste(missing_files, collapse = ", ")
  )
}

####################
# HELPER FUNCTIONS #
####################

read_tsv = function(path) {
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

mean_expected_at_lod = function(
  sweep,
  pairs,
  lod
) {

  key = paste(pairs$query, pairs$target, sep = "|")
  dat = sweep[
    paste(sweep$query, sweep$target, sep = "|") %in% key &
    abs(sweep$lod - lod) < 1e-12,
    ,
    drop = FALSE
  ]
  mean(dat$expected_rediscovery, na.rm = TRUE)
}

################
# READ OUTPUTS #
################

study = read_tsv(paste0(out_dir, "study_design.tsv"))
observed = read_tsv(
  paste0(out_dir, "observed_individual_rediscovery.tsv")
)
sweep = read_tsv(
  paste0(out_dir, "pair_conditional_power_sweep.tsv")
)
pair_validation = read_tsv(
  paste0(out_dir, "pair_conditional_validation_summary.tsv")
)
analytic_validation = read_tsv(
  paste0(out_dir, "analytic_power_validation.tsv")
)
chance = read_tsv(
  paste0(out_dir, "chance_overlap_null.tsv")
)
summary_table = read_tsv(
  paste0(out_dir, "reviewer_power_summary.tsv")
)
presentation_table = read_tsv(
  paste0(out_dir, "reviewer_power_table.tsv")
)
asymmetry = read_tsv(
  paste0(out_dir, "reviewer_reciprocal_asymmetry.tsv")
)
asymmetry_summary = read_tsv(
  paste0(out_dir, "reviewer_reciprocal_asymmetry_summary.tsv")
)
asymmetry_key = read_tsv(
  paste0(out_dir, "reviewer_reciprocal_asymmetry_key.tsv")
)
supporting_validation = read_tsv(
  paste0(out_dir, "reviewer_supporting_validation.tsv")
)

checks = list()
check_i = 1L

########################################
# UPSTREAM PAIR-CONDITIONAL VALIDATION #
########################################

checks[[check_i]] = add_check(
  "upstream_pair_conditional_validator_passed",
  nrow(pair_validation) > 0L && all(pair_validation$passed),
  paste(
    "passed =",
    sum(pair_validation$passed),
    "/",
    nrow(pair_validation)
  )
)
check_i = check_i + 1L

############################
# CHANCE-NULL CONSISTENCY  #
############################

chance_key = paste(chance$query, chance$target, sep = "|")
observed_key = paste(observed$query, observed$target, sep = "|")

checks[[check_i]] = add_check(
  "chance_null_has_all_90_ordered_pairs",
  nrow(chance) == 90L &&
    length(unique(chance_key)) == 90L &&
    setequal(chance_key, observed_key),
  paste(
    "rows =", nrow(chance),
    "; unique pairs =", length(unique(chance_key))
  )
)
check_i = check_i + 1L

chance_observed = merge(
  chance[, c(
    "query",
    "target",
    "n_query_eligible",
    "n_matches",
    "observed_rediscovery"
  )],
  observed[, c(
    "query",
    "target",
    "n_query_eligible",
    "n_matches",
    "observed_rediscovery"
  )],
  by = c("query", "target"),
  suffixes = c("_chance", "_observed"),
  all = TRUE
)

count_same = all(
  chance_observed$n_query_eligible_chance ==
  chance_observed$n_query_eligible_observed
)
match_same = all(
  chance_observed$n_matches_chance ==
  chance_observed$n_matches_observed
)
obs_diff = max(
  abs(
    chance_observed$observed_rediscovery_chance -
    chance_observed$observed_rediscovery_observed
  ),
  na.rm = TRUE
)

checks[[check_i]] = add_check(
  "chance_null_uses_exact_figure3_denominators",
  count_same && match_same && obs_diff < 1e-14,
  paste(
    "denominators same =", count_same,
    "; matches same =", match_same,
    "; max observed diff =", obs_diff
  )
)
check_i = check_i + 1L

chance_values = chance$chance_rediscovery[
  is.finite(chance$chance_rediscovery)
]
expected_match_values = chance$expected_chance_matches[
  is.finite(chance$expected_chance_matches)
]
checks[[check_i]] = add_check(
  "chance_probabilities_are_valid",
  all(chance_values >= 0 & chance_values <= 1) &&
    all(expected_match_values >= 0) &&
    all(
      chance$expected_chance_matches <=
      chance$n_query_eligible + 1e-12
    ),
  paste(
    "chance range =",
    paste(signif(range(chance_values), 6), collapse = " to ")
  )
)
check_i = check_i + 1L

prior_chance_path = paste0(
  wd,
  "asys/review/power/chance_overlap_null.tsv"
)
if (file.exists(prior_chance_path)) {

  prior_chance = read_tsv(prior_chance_path)
  compare = merge(
    chance[, c(
      "query",
      "target",
      "expected_chance_matches",
      "chance_rediscovery"
    )],
    prior_chance[, c(
      "query",
      "target",
      "expected_chance_matches",
      "chance_rediscovery"
    )],
    by = c("query", "target"),
    suffixes = c("_new", "_prior"),
    all = TRUE
  )

  expected_diff = max(
    abs(
      compare$expected_chance_matches_new -
      compare$expected_chance_matches_prior
    ),
    na.rm = TRUE
  )
  fraction_diff = max(
    abs(
      compare$chance_rediscovery_new -
      compare$chance_rediscovery_prior
    ),
    na.rm = TRUE
  )

  checks[[check_i]] = add_check(
    "chance_null_matches_prior_implementation",
    nrow(compare) == 90L &&
      expected_diff < 1e-12 &&
      fraction_diff < 1e-12,
    paste(
      "max expected-match diff =", expected_diff,
      "; max fraction diff =", fraction_diff
    )
  )
  check_i = check_i + 1L
}

########################################
# REVIEWER TABLE GROUP DEFINITIONS     #
########################################

expected_groups = c(
  "all_ordered_pairs",
  "Albert_2018_query",
  "other_studies_query",
  "mRNA_to_mRNA",
  "mRNA_to_protein",
  "protein_to_mRNA",
  "protein_to_protein"
)
expected_counts = c(90L, 9L, 81L, 30L, 24L, 24L, 12L)

checks[[check_i]] = add_check(
  "reviewer_table_has_only_critical_seven_rows",
  nrow(summary_table) == length(expected_groups) &&
    identical(summary_table$group, expected_groups) &&
    identical(as.integer(summary_table$n_pairs), expected_counts) &&
    nrow(presentation_table) == length(expected_groups),
  paste(
    "groups =", paste(summary_table$group, collapse = ", ")
  )
)
check_i = check_i + 1L

qmeta = study[, c("study", "analyte"), drop = FALSE]
names(qmeta) = c("query", "query_analyte")
tmeta = study[, c("study", "analyte"), drop = FALSE]
names(tmeta) = c("target", "target_analyte")

pair_meta = merge(
  observed[, c("query", "target"), drop = FALSE],
  qmeta,
  by = "query",
  all.x = TRUE,
  sort = FALSE
)
pair_meta = merge(
  pair_meta,
  tmeta,
  by = "target",
  all.x = TRUE,
  sort = FALSE
)
pair_meta$analyte_pair = paste0(
  pair_meta$query_analyte,
  "_to_",
  pair_meta$target_analyte
)

pair_sets = list(
  all_ordered_pairs = pair_meta,
  Albert_2018_query = pair_meta[
    pair_meta$query == "Albert_2018",
    ,
    drop = FALSE
  ],
  other_studies_query = pair_meta[
    pair_meta$query != "Albert_2018",
    ,
    drop = FALSE
  ],
  mRNA_to_mRNA = pair_meta[
    pair_meta$analyte_pair == "mRNA_to_mRNA",
    ,
    drop = FALSE
  ],
  mRNA_to_protein = pair_meta[
    pair_meta$analyte_pair == "mRNA_to_protein",
    ,
    drop = FALSE
  ],
  protein_to_mRNA = pair_meta[
    pair_meta$analyte_pair == "protein_to_mRNA",
    ,
    drop = FALSE
  ],
  protein_to_protein = pair_meta[
    pair_meta$analyte_pair == "protein_to_protein",
    ,
    drop = FALSE
  ]
)

summary_recomputed = do.call(
  rbind,
  lapply(
    expected_groups,
    function(group) {

      pairs = pair_sets[[group]]
      key = paste(pairs$query, pairs$target, sep = "|")
      chance_dat = chance[
        paste(chance$query, chance$target, sep = "|") %in% key,
        ,
        drop = FALSE
      ]

      data.frame(
        group = group,
        n_pairs = nrow(pairs),
        mean_chance = mean(
          chance_dat$chance_rediscovery,
          na.rm = TRUE
        ),
        mean_observed = mean(
          chance_dat$observed_rediscovery,
          na.rm = TRUE
        ),
        expected_lod2 = mean_expected_at_lod(
          sweep,
          pairs,
          2
        ),
        expected_lod4_5 = mean_expected_at_lod(
          sweep,
          pairs,
          4.5
        ),
        expected_lod7 = mean_expected_at_lod(
          sweep,
          pairs,
          7
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)
summary_recomputed$residual_lod4_5 =
  summary_recomputed$mean_observed -
  summary_recomputed$expected_lod4_5
rownames(summary_recomputed) = NULL

numeric_columns = c(
  "mean_chance",
  "mean_observed",
  "expected_lod2",
  "expected_lod4_5",
  "expected_lod7",
  "residual_lod4_5"
)
summary_numeric_diff = max(
  abs(
    as.matrix(summary_table[, numeric_columns]) -
    as.matrix(summary_recomputed[, numeric_columns])
  ),
  na.rm = TRUE
)

checks[[check_i]] = add_check(
  "reviewer_power_summary_recomputes_from_primary_outputs",
  summary_numeric_diff < 1e-12,
  paste("max absolute difference =", summary_numeric_diff)
)
check_i = check_i + 1L

monotone_summary =
  summary_table$expected_lod2 >=
  summary_table$expected_lod4_5 - 1e-12 &
  summary_table$expected_lod4_5 >=
  summary_table$expected_lod7 - 1e-12
checks[[check_i]] = add_check(
  "reviewer_expected_values_decline_with_target_lod",
  all(monotone_summary),
  paste("rows monotone =", sum(monotone_summary), "/", nrow(summary_table))
)
check_i = check_i + 1L

########################################
# RECIPROCAL ASYMMETRY RECOMPUTATION  #
########################################

reverse = data.frame(
  query = sweep$target,
  target = sweep$query,
  lod = sweep$lod,
  observed_reverse = sweep$observed_rediscovery,
  expected_reverse = sweep$expected_rediscovery,
  stringsAsFactors = FALSE
)
independent_asym = merge(
  sweep[, c(
    "query",
    "target",
    "lod",
    "observed_rediscovery",
    "expected_rediscovery"
  )],
  reverse,
  by = c("query", "target", "lod"),
  all = FALSE,
  sort = FALSE
)
independent_asym = independent_asym[
  independent_asym$query < independent_asym$target,
  ,
  drop = FALSE
]
independent_asym$observed_asymmetry =
  independent_asym$observed_rediscovery -
  independent_asym$observed_reverse
independent_asym$expected_asymmetry =
  independent_asym$expected_rediscovery -
  independent_asym$expected_reverse

independent_asym_summary = do.call(
  rbind,
  lapply(
    sort(unique(independent_asym$lod)),
    function(lod) {
      dat = independent_asym[
        abs(independent_asym$lod - lod) < 1e-12,
        ,
        drop = FALSE
      ]
      data.frame(
        lod = lod,
        n_unordered_pairs = nrow(dat),
        spearman = stats::cor(
          dat$observed_asymmetry,
          dat$expected_asymmetry,
          method = "spearman"
        ),
        pearson = stats::cor(
          dat$observed_asymmetry,
          dat$expected_asymmetry,
          method = "pearson"
        ),
        direction_agreement = mean(
          sign(dat$observed_asymmetry) ==
          sign(dat$expected_asymmetry)
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

asym_cols = c(
  "spearman",
  "pearson",
  "direction_agreement"
)
asym_diff = max(
  abs(
    as.matrix(asymmetry_summary[, asym_cols]) -
    as.matrix(independent_asym_summary[, asym_cols])
  ),
  na.rm = TRUE
)

checks[[check_i]] = add_check(
  "reciprocal_asymmetry_has_45_pairs_at_each_lod",
  nrow(asymmetry_summary) == 21L &&
    all(asymmetry_summary$n_unordered_pairs == 45L),
  paste(
    "LOD rows =", nrow(asymmetry_summary),
    "; pair-count range =",
    paste(range(asymmetry_summary$n_unordered_pairs), collapse = " to ")
  )
)
check_i = check_i + 1L

checks[[check_i]] = add_check(
  "reciprocal_asymmetry_summary_recomputes",
  asym_diff < 1e-12,
  paste("max absolute difference =", asym_diff)
)
check_i = check_i + 1L

key_recomputed = data.frame(
  metric = c(
    "spearman_reciprocal_asymmetry",
    "pearson_reciprocal_asymmetry",
    "reciprocal_direction_agreement"
  ),
  minimum = c(
    min(independent_asym_summary$spearman),
    min(independent_asym_summary$pearson),
    min(independent_asym_summary$direction_agreement)
  ),
  maximum = c(
    max(independent_asym_summary$spearman),
    max(independent_asym_summary$pearson),
    max(independent_asym_summary$direction_agreement)
  ),
  stringsAsFactors = FALSE
)
key_compare = merge(
  asymmetry_key[, c("metric", "minimum", "maximum")],
  key_recomputed,
  by = "metric",
  suffixes = c("_output", "_recomputed"),
  all = TRUE
)
key_diff = max(
  abs(
    c(
      key_compare$minimum_output - key_compare$minimum_recomputed,
      key_compare$maximum_output - key_compare$maximum_recomputed
    )
  ),
  na.rm = TRUE
)
checks[[check_i]] = add_check(
  "reciprocal_asymmetry_key_range_recomputes",
  nrow(key_compare) == 3L && key_diff < 1e-12,
  paste("max absolute difference =", key_diff)
)
check_i = check_i + 1L

########################################
# MINIMAL SUPPORTING VALIDATION OUTPUT #
########################################

max_validation_error = max(
  abs(
    analytic_validation$simulated_power -
    analytic_validation$analytic_power
  ),
  abs(
    analytic_validation$simulated_same_direction -
    analytic_validation$analytic_same_direction
  ),
  na.rm = TRUE
)

checks[[check_i]] = add_check(
  "direct_simulation_error_remains_small",
  max_validation_error < 0.02,
  paste("maximum absolute error =", signif(max_validation_error, 6))
)
check_i = check_i + 1L

checks[[check_i]] = add_check(
  "supporting_validation_table_is_minimal",
  nrow(supporting_validation) == 4L,
  paste("rows =", nrow(supporting_validation))
)
check_i = check_i + 1L

##################
# WRITE / REPORT #
##################

validation_summary = do.call(rbind, checks)
write.table(
  validation_summary,
  file = validation_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nREVIEWER POWER SUMMARY VALIDATION\n")
print(validation_summary, row.names = FALSE)
cat("\n")

failed = validation_summary$check[!validation_summary$passed]
if (length(failed) > 0L) {
  stop(
    "Reviewer-summary validation failed: ",
    paste(failed, collapse = ", ")
  )
}

message("All reviewer-summary checks passed.")
message("Validation summary: ", validation_out)
