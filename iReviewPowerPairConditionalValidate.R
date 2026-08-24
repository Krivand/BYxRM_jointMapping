############################################
############################################
## GENETICS REVIEW POWER VALIDATOR        ##
## PAIR-CONDITIONAL REDISCOVERY           ##
############################################
############################################

#################
# PRELIMINARIES #
#################

wd = "/PHShome/kv1026/p05/"
out_dir = paste0(wd, "asys/review/power_pair_conditional/")
analysis_script = paste0(wd, "iReviewPowerPairConditional.R")
validation_out = paste0(
  out_dir,
  "pair_conditional_validation_summary.tsv"
)

required_files = c(
  "study_design.tsv",
  "observed_individual_rediscovery.tsv",
  "pair_eligible_query_qtls.tsv.gz",
  "pair_conditional_power_sweep.tsv",
  "pair_conditional_power_summary.tsv",
  "pair_conditional_power_envelope.tsv",
  "analytic_power_validation.tsv",
  "model_definition.tsv"
)

missing_files = required_files[
  !file.exists(paste0(out_dir, required_files))
]
if (length(missing_files) > 0L) {
  stop(
    "Missing pair-conditional outputs: ",
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

lod_to_rcrit_independent = function(lod, n) {
  sqrt(pmax(0, 1 - 10 ^ (-2 * lod / n)))
}

same_direction_prob_independent = function(rho, n, lod) {

  rcrit = lod_to_rcrit_independent(lod = lod, n = n)
  df = n - 2
  tcrit = rcrit * sqrt(df / (1 - rcrit ^ 2))
  ncp = rho * sqrt((n - 1) / (1 - rho ^ 2))

  p_positive = stats::pt(
    tcrit,
    df = df,
    ncp = ncp,
    lower.tail = FALSE
  )
  p_negative = stats::pt(
    -tcrit,
    df = df,
    ncp = ncp,
    lower.tail = TRUE
  )

  ifelse(rho > 0, p_positive, p_negative)
}

add_check = function(name, passed, detail) {
  data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

################
# READ OUTPUTS #
################

study_design = read_tsv(paste0(out_dir, "study_design.tsv"))
observed = read_tsv(
  paste0(out_dir, "observed_individual_rediscovery.tsv")
)
eligible = read_tsv(
  gzfile(paste0(out_dir, "pair_eligible_query_qtls.tsv.gz"))
)
sweep = read_tsv(
  paste0(out_dir, "pair_conditional_power_sweep.tsv")
)
summary_table = read_tsv(
  paste0(out_dir, "pair_conditional_power_summary.tsv")
)
envelope = read_tsv(
  paste0(out_dir, "pair_conditional_power_envelope.tsv")
)
analytic_validation = read_tsv(
  paste0(out_dir, "analytic_power_validation.tsv")
)
model_definition = read_tsv(
  paste0(out_dir, "model_definition.tsv")
)

checks = list()
check_i = 1L

############################
# BASIC DIMENSIONAL CHECKS #
############################

checks[[check_i]] = add_check(
  "ten_studies",
  nrow(study_design) == 10L && length(unique(study_design$study)) == 10L,
  paste("n_studies =", nrow(study_design))
)
check_i = check_i + 1L

checks[[check_i]] = add_check(
  "ninety_ordered_pairs",
  nrow(observed) == 90L,
  paste("n_pairs =", nrow(observed))
)
check_i = check_i + 1L

lod_values = sort(unique(sweep$lod))
expected_lod = seq(2, 7, by = 0.25)
checks[[check_i]] = add_check(
  "lod_grid_exact",
  length(lod_values) == length(expected_lod) &&
    max(abs(lod_values - expected_lod)) < 1e-12,
  paste(
    "LOD grid =",
    paste(signif(lod_values, 4), collapse = ",")
  )
)
check_i = check_i + 1L

checks[[check_i]] = add_check(
  "sweep_row_count",
  nrow(sweep) == 90L * length(expected_lod),
  paste("n_rows =", nrow(sweep))
)
check_i = check_i + 1L

###############################
# FIXED-DENOMINATOR VALIDATION #
###############################

eligible_counts = aggregate(
  query_qtl_id ~ query + target,
  data = eligible,
  FUN = length
)
names(eligible_counts)[3] = "eligible_rows"

count_check = merge(
  observed[, c("query", "target", "n_query_eligible")],
  eligible_counts,
  by = c("query", "target"),
  all.x = TRUE,
  sort = FALSE
)
count_check$eligible_rows[is.na(count_check$eligible_rows)] = 0L

checks[[check_i]] = add_check(
  "eligible_denominator_matches_export",
  all(count_check$n_query_eligible == count_check$eligible_rows),
  paste(
    "max count difference =",
    max(abs(count_check$n_query_eligible - count_check$eligible_rows))
  )
)
check_i = check_i + 1L

observed_ratio = observed$n_matches / observed$n_query_eligible
observed_ratio[observed$n_query_eligible == 0L] = NA_real_
ratio_diff = abs(observed_ratio - observed$observed_rediscovery)
checks[[check_i]] = add_check(
  "observed_fraction_recomputes",
  max(ratio_diff, na.rm = TRUE) < 1e-14,
  paste("max absolute difference =", max(ratio_diff, na.rm = TRUE))
)
check_i = check_i + 1L

########################################
# INDEPENDENT POWER RECOMPUTATION      #
########################################

recomputed = vector("list", nrow(sweep))

for (i in seq_len(nrow(sweep))) {

  query = sweep$query[[i]]
  target = sweep$target[[i]]
  lod = sweep$lod[[i]]
  target_n = study_design$n[match(target, study_design$study)]

  dat = eligible[
    eligible$query == query & eligible$target == target,
    ,
    drop = FALSE
  ]

  if (nrow(dat) == 0L) {
    expected = NA_real_
  } else {
    expected = mean(
      same_direction_prob_independent(
        rho = dat$cor_coef,
        n = target_n,
        lod = lod
      )
    )
  }

  recomputed[[i]] = expected
}

recomputed = unlist(recomputed)
power_diff = abs(recomputed - sweep$expected_rediscovery)

checks[[check_i]] = add_check(
  "independent_probability_recomputation",
  max(power_diff, na.rm = TRUE) < 1e-12,
  paste("max absolute difference =", max(power_diff, na.rm = TRUE))
)
check_i = check_i + 1L

########################################
# MONOTONICITY UNDER FIXED QUERY SET   #
########################################

pair_key = paste(sweep$query, sweep$target, sep = "|")
monotone_ok = TRUE
worst_increase = 0

for (key in unique(pair_key)) {
  dat = sweep[pair_key == key, , drop = FALSE]
  dat = dat[order(dat$lod), , drop = FALSE]
  d = diff(dat$expected_rediscovery)
  d = d[is.finite(d)]
  if (length(d) > 0L) {
    worst_increase = max(worst_increase, max(d))
    if (any(d > 1e-12)) monotone_ok = FALSE
  }
}

checks[[check_i]] = add_check(
  "expected_rediscovery_nonincreasing_with_lod",
  monotone_ok,
  paste("largest upward step =", signif(worst_increase, 6))
)
check_i = check_i + 1L

summary_ordered = summary_table[order(summary_table$lod), , drop = FALSE]
mean_d = diff(summary_ordered$mean_expected)
checks[[check_i]] = add_check(
  "mean_expected_nonincreasing_with_lod",
  all(mean_d <= 1e-12),
  paste("largest upward step =", signif(max(mean_d), 6))
)
check_i = check_i + 1L

###########################
# PROBABILITY RANGE CHECKS #
###########################

prob_columns = c(
  "expected_rediscovery",
  "expected_any_detection",
  "expected_wrong_direction"
)
prob_values = unlist(sweep[, prob_columns], use.names = FALSE)
prob_values = prob_values[is.finite(prob_values)]

checks[[check_i]] = add_check(
  "all_probabilities_in_unit_interval",
  all(prob_values >= -1e-14 & prob_values <= 1 + 1e-14),
  paste(
    "range =",
    paste(signif(range(prob_values), 6), collapse = " to ")
  )
)
check_i = check_i + 1L

##################################
# DIRECT-SIMULATION VALIDATION   #
##################################

mc_tolerance = function(p, n_rep) {
  5 * sqrt(pmax(p * (1 - p), 1 / n_rep) / n_rep) + 1 / n_rep
}

power_tol = mc_tolerance(
  analytic_validation$analytic_power,
  analytic_validation$n_rep
)
dir_tol = mc_tolerance(
  analytic_validation$analytic_same_direction,
  analytic_validation$n_rep
)

power_dev = abs(
  analytic_validation$simulated_power -
  analytic_validation$analytic_power
)
dir_dev = abs(
  analytic_validation$simulated_same_direction -
  analytic_validation$analytic_same_direction
)

checks[[check_i]] = add_check(
  "direct_simulation_power_within_mc_tolerance",
  all(power_dev <= power_tol),
  paste(
    "max dev =", signif(max(power_dev), 6),
    "; max tol =", signif(max(power_tol), 6)
  )
)
check_i = check_i + 1L

checks[[check_i]] = add_check(
  "direct_simulation_direction_within_mc_tolerance",
  all(dir_dev <= dir_tol),
  paste(
    "max dev =", signif(max(dir_dev), 6),
    "; max tol =", signif(max(dir_tol), 6)
  )
)
check_i = check_i + 1L

#####################################
# MODEL-SPECIFICATION SOURCE AUDIT  #
#####################################

script_text = if (file.exists(analysis_script)) {
  paste(readLines(analysis_script, warn = FALSE), collapse = "\n")
} else {
  ""
}

forbidden = c(
  "reference_study = \"Albert_2018\"",
  "expected_pairwise_shared(",
  "query_detection_weight_sum",
  "query_ascertainment_summary("
)
forbidden_present = vapply(
  forbidden,
  function(x) grepl(x, script_text, fixed = TRUE),
  logical(1)
)

checks[[check_i]] = add_check(
  "no_albert_reference_or_query_reascertainment_code",
  !any(forbidden_present),
  if (any(forbidden_present)) {
    paste("found:", paste(forbidden[forbidden_present], collapse = "; "))
  } else {
    "forbidden legacy constructs absent"
  }
)
check_i = check_i + 1L

required_model_text = c(
  "actual Figure 3 eligible query QTLs for each ordered pair",
  "Albert_2018 is treated like every other query study",
  "the observed query denominator is fixed"
)
model_blob = paste(model_definition$definition, collapse = "\n")
model_present = vapply(
  required_model_text,
  function(x) grepl(x, model_blob, fixed = TRUE),
  logical(1)
)

checks[[check_i]] = add_check(
  "model_definition_states_fixed_pair_conditional_design",
  all(model_present),
  paste("required statements present =", sum(model_present), "/", length(model_present))
)
check_i = check_i + 1L

############################################
# OPTIONAL COMPARISON TO PRIOR OBSERVED RUN #
############################################

prior_observed_path = paste0(
  wd,
  "asys/review/power/observed_individual_rediscovery.tsv"
)

if (file.exists(prior_observed_path)) {

  prior = read_tsv(prior_observed_path)
  old_new = merge(
    prior[, c(
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
    suffixes = c("_prior", "_new"),
    all = TRUE
  )

  denominator_same = all(
    old_new$n_query_eligible_prior == old_new$n_query_eligible_new
  )
  matches_same = all(old_new$n_matches_prior == old_new$n_matches_new)
  rediscovery_diff = max(
    abs(
      old_new$observed_rediscovery_prior -
      old_new$observed_rediscovery_new
    ),
    na.rm = TRUE
  )

  checks[[check_i]] = add_check(
    "observed_figure3_logic_matches_prior_run",
    denominator_same && matches_same && rediscovery_diff < 1e-14,
    paste(
      "denominators same =", denominator_same,
      "; matches same =", matches_same,
      "; max fraction diff =", rediscovery_diff
    )
  )
  check_i = check_i + 1L
}

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

cat("\nPAIR-CONDITIONAL POWER VALIDATION\n")
print(validation_summary, row.names = FALSE)
cat("\n")

failed = validation_summary$check[!validation_summary$passed]
if (length(failed) > 0L) {
  stop(
    "Validation failed: ",
    paste(failed, collapse = ", ")
  )
}

message("All pair-conditional power checks passed.")
message("Validation summary: ", validation_out)
