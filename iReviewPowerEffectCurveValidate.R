################################################
################################################
## GENETICS REVIEW POWER ANALYSIS VALIDATOR   ##
## EFFECT-SIZE POWER CURVES                   ##
################################################
################################################

#################
# PRELIMINARIES #
#################

wd = "/PHShome/kv1026/p05/"
out_dir = paste0(wd, "asys/review/power_effect_curve/")
validation_out = paste0(out_dir, "effect_curve_validation_summary.tsv")

required_files = c(
  "study_design.tsv",
  "frank_qtl_count_check.tsv",
  "effect_size_power_surface.tsv.gz",
  "effect_size_power_thresholds.tsv",
  "albert_reference_effects.tsv.gz",
  "albert_effect_summary.tsv",
  "albert_reference_expected_detection.tsv",
  "analytic_power_validation.tsv",
  "model_definition.tsv"
)
missing_files = required_files[
  !file.exists(paste0(out_dir, required_files))
]
if (length(missing_files) > 0L) {
  stop(
    "Missing effect-curve outputs: ",
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

read_tsv_gz = function(path) {
  con = gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(
    con,
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

independent_probability = function(rho, n, lod) {

  rcrit = sqrt(pmax(0, 1 - 10 ^ (-2 * lod / n)))
  df = n - 2
  tcrit = rcrit * sqrt(df / (1 - rcrit ^ 2))
  ncp = rho * sqrt((n - 1) / (1 - rho ^ 2))

  positive = stats::pt(
    tcrit,
    df = df,
    ncp = ncp,
    lower.tail = FALSE
  )
  negative = stats::pt(
    -tcrit,
    df = df,
    ncp = ncp,
    lower.tail = TRUE
  )

  c(
    same = positive,
    wrong = negative,
    any = positive + negative
  )
}

is_nondecreasing = function(x, tolerance = 1e-12) {
  if (length(x) <= 1L) return(TRUE)
  all(diff(x) >= -tolerance)
}

is_nonincreasing = function(x, tolerance = 1e-12) {
  if (length(x) <= 1L) return(TRUE)
  all(diff(x) <= tolerance)
}

################
# READ OUTPUTS #
################

study = read_tsv(paste0(out_dir, "study_design.tsv"))
frank_check = read_tsv(paste0(out_dir, "frank_qtl_count_check.tsv"))
power_surface = read_tsv_gz(
  paste0(out_dir, "effect_size_power_surface.tsv.gz")
)
power_thresholds = read_tsv(
  paste0(out_dir, "effect_size_power_thresholds.tsv")
)
albert_effects = read_tsv_gz(
  paste0(out_dir, "albert_reference_effects.tsv.gz")
)
albert_summary = read_tsv(
  paste0(out_dir, "albert_effect_summary.tsv")
)
albert_detection = read_tsv(
  paste0(out_dir, "albert_reference_expected_detection.tsv")
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
# STUDY / SOURCE INTEGRITY #
############################

checks[[check_i]] = add_check(
  "ten_independent_individual_studies",
  nrow(study) == 10L &&
    length(unique(study$study)) == 10L &&
    !("AlbertAlt_2018" %in% study$study),
  paste(
    "rows =", nrow(study),
    "; studies =", paste(study$study, collapse = ", ")
  )
)
check_i = check_i + 1L

checks[[check_i]] = add_check(
  "frank_qtl_counts_match_current_individual_qtls",
  nrow(frank_check) == 10L && all(frank_check$match),
  paste(
    "matched =", sum(frank_check$match),
    "/", nrow(frank_check)
  )
)
check_i = check_i + 1L

########################################
# POWER SURFACE STRUCTURE AND FORMULA  #
########################################

expected_rho = seq(0, 0.80, by = 0.005)
expected_lod = seq(2, 7, by = 0.25)
expected_surface_rows =
  nrow(study) * length(expected_rho) * length(expected_lod)

checks[[check_i]] = add_check(
  "power_surface_has_exact_requested_grid",
  nrow(power_surface) == expected_surface_rows &&
    isTRUE(all.equal(
      sort(unique(power_surface$abs_r)),
      expected_rho,
      tolerance = 1e-14
    )) &&
    isTRUE(all.equal(
      sort(unique(power_surface$lod)),
      expected_lod,
      tolerance = 1e-14
    )),
  paste(
    "rows =", nrow(power_surface),
    "; expected =", expected_surface_rows,
    "; rho points =", length(unique(power_surface$abs_r)),
    "; LOD points =", length(unique(power_surface$lod))
  )
)
check_i = check_i + 1L

prob_cols = c(
  "p_same_direction",
  "p_wrong_direction",
  "p_any_detection"
)
prob_values = unlist(power_surface[, prob_cols], use.names = FALSE)
checks[[check_i]] = add_check(
  "power_surface_probabilities_are_valid",
  all(is.finite(prob_values)) &&
    all(prob_values >= -1e-14 & prob_values <= 1 + 1e-14) &&
    max(
      abs(
        power_surface$p_any_detection -
        power_surface$p_same_direction -
        power_surface$p_wrong_direction
      )
    ) < 1e-12,
  paste(
    "probability range =",
    paste(signif(range(prob_values), 7), collapse = " to ")
  )
)
check_i = check_i + 1L

independent_rows = do.call(
  rbind,
  lapply(
    seq_len(nrow(power_surface)),
    function(i) {
      p = independent_probability(
        rho = power_surface$abs_r[[i]],
        n = power_surface$n[[i]],
        lod = power_surface$lod[[i]]
      )
      data.frame(
        same = p[["same"]],
        wrong = p[["wrong"]],
        any = p[["any"]],
        stringsAsFactors = FALSE
      )
    }
  )
)
formula_diff = max(
  abs(power_surface$p_same_direction - independent_rows$same),
  abs(power_surface$p_wrong_direction - independent_rows$wrong),
  abs(power_surface$p_any_detection - independent_rows$any),
  na.rm = TRUE
)
checks[[check_i]] = add_check(
  "independent_noncentral_t_recomputation_matches",
  formula_diff < 1e-12,
  paste("maximum absolute difference =", formula_diff)
)
check_i = check_i + 1L

########################################
# MONOTONICITY OF THE POWER ENGINE     #
########################################

rho_groups = split(
  power_surface,
  paste(power_surface$study, power_surface$lod, sep = "|")
)
rho_monotone = vapply(
  rho_groups,
  function(dat) {
    dat = dat[order(dat$abs_r), , drop = FALSE]
    is_nondecreasing(dat$p_same_direction)
  },
  logical(1)
)
checks[[check_i]] = add_check(
  "power_increases_with_effect_magnitude",
  all(rho_monotone),
  paste(
    "monotone study-LOD curves =",
    sum(rho_monotone), "/", length(rho_monotone)
  )
)
check_i = check_i + 1L

lod_groups = split(
  power_surface,
  paste(power_surface$study, power_surface$abs_r, sep = "|")
)
lod_monotone = vapply(
  lod_groups,
  function(dat) {
    dat = dat[order(dat$lod), , drop = FALSE]
    is_nonincreasing(dat$p_same_direction)
  },
  logical(1)
)
checks[[check_i]] = add_check(
  "power_decreases_with_more_stringent_lod",
  all(lod_monotone),
  paste(
    "monotone study-effect curves =",
    sum(lod_monotone), "/", length(lod_monotone)
  )
)
check_i = check_i + 1L

sample_power = aggregate(
  p_same_direction ~ n + lod + abs_r,
  data = power_surface[power_surface$abs_r >= 0.05, , drop = FALSE],
  FUN = mean
)
n_groups = split(
  sample_power,
  paste(sample_power$lod, sample_power$abs_r, sep = "|")
)
n_monotone = vapply(
  n_groups,
  function(dat) {
    dat = dat[order(dat$n), , drop = FALSE]
    is_nondecreasing(dat$p_same_direction, tolerance = 1e-10)
  },
  logical(1)
)
checks[[check_i]] = add_check(
  "power_increases_with_sample_size",
  all(n_monotone),
  paste(
    "monotone LOD-effect curves =",
    sum(n_monotone), "/", length(n_monotone)
  )
)
check_i = check_i + 1L

########################################
# 20/50/80% POWER INVERSION            #
########################################

expected_threshold_rows = nrow(study) * length(expected_lod) * 3L
threshold_structure_ok =
  nrow(power_thresholds) == expected_threshold_rows &&
  setequal(unique(power_thresholds$target_power), c(0.2, 0.5, 0.8)) &&
  all(is.finite(power_thresholds$abs_r_required)) &&
  all(
    power_thresholds$abs_r_required >= 0 &
    power_thresholds$abs_r_required < 1
  )

threshold_recomputed = vapply(
  seq_len(nrow(power_thresholds)),
  function(i) {
    independent_probability(
      rho = power_thresholds$abs_r_required[[i]],
      n = power_thresholds$n[[i]],
      lod = power_thresholds$lod[[i]]
    )[["same"]]
  },
  numeric(1)
)
threshold_error = max(
  abs(threshold_recomputed - power_thresholds$target_power),
  na.rm = TRUE
)

checks[[check_i]] = add_check(
  "power_threshold_inversion_is_accurate",
  threshold_structure_ok && threshold_error < 1e-8,
  paste(
    "rows =", nrow(power_thresholds),
    "; maximum probability error =", threshold_error
  )
)
check_i = check_i + 1L

threshold_groups = split(
  power_thresholds,
  paste(power_thresholds$study, power_thresholds$lod, sep = "|")
)
threshold_order = vapply(
  threshold_groups,
  function(dat) {
    dat = dat[order(dat$target_power), , drop = FALSE]
    is_nondecreasing(dat$abs_r_required)
  },
  logical(1)
)
checks[[check_i]] = add_check(
  "higher_power_requires_larger_effect",
  all(threshold_order),
  paste(
    "ordered study-LOD sets =",
    sum(threshold_order), "/", length(threshold_order)
  )
)
check_i = check_i + 1L

########################################
# ALBERT EMPIRICAL REFERENCE           #
########################################

albert_n_qtl = study$n_qtl[study$study == "Albert_2018"]
albert_effect_integrity =
  nrow(albert_summary) == 1L &&
  nrow(albert_effects) == albert_n_qtl &&
  albert_summary$n_qtl[[1]] == albert_n_qtl &&
  all(is.finite(albert_effects$abs_r)) &&
  all(albert_effects$abs_r > 0 & albert_effects$abs_r < 1) &&
  max(abs(albert_effects$abs_r - abs(albert_effects$cor_coef))) < 1e-14

checks[[check_i]] = add_check(
  "albert_reference_is_exact_discovered_qtl_effect_set",
  albert_effect_integrity,
  paste(
    "Albert current QTL count =", albert_n_qtl,
    "; reference effects =", nrow(albert_effects)
  )
)
check_i = check_i + 1L

expected_targets = setdiff(study$study, "Albert_2018")
expected_scopes = c("all_albert_qtls", "target_measured_genes")
expected_scales = c(0.6, 0.8, 1.0)
expected_albert_rows =
  length(expected_targets) *
  length(expected_scopes) *
  length(expected_scales) *
  length(expected_lod)

checks[[check_i]] = add_check(
  "albert_calibration_has_exact_target_scope_scale_lod_grid",
  nrow(albert_detection) == expected_albert_rows &&
    setequal(unique(albert_detection$target), expected_targets) &&
    setequal(unique(albert_detection$reference_scope), expected_scopes) &&
    setequal(unique(albert_detection$effect_scale), expected_scales) &&
    isTRUE(all.equal(
      sort(unique(albert_detection$lod)),
      expected_lod,
      tolerance = 1e-14
    )),
  paste(
    "rows =", nrow(albert_detection),
    "; expected =", expected_albert_rows
  )
)
check_i = check_i + 1L

scale_groups = split(
  albert_detection,
  paste(
    albert_detection$target,
    albert_detection$reference_scope,
    albert_detection$lod,
    sep = "|"
  )
)
scale_monotone = vapply(
  scale_groups,
  function(dat) {
    dat = dat[order(dat$effect_scale), , drop = FALSE]
    is_nondecreasing(dat$mean_same_direction_detection)
  },
  logical(1)
)
checks[[check_i]] = add_check(
  "albert_expected_detection_increases_with_effect_scale",
  all(scale_monotone),
  paste(
    "monotone target-scope-LOD sets =",
    sum(scale_monotone), "/", length(scale_monotone)
  )
)
check_i = check_i + 1L

albert_lod_groups = split(
  albert_detection,
  paste(
    albert_detection$target,
    albert_detection$reference_scope,
    albert_detection$effect_scale,
    sep = "|"
  )
)
albert_lod_monotone = vapply(
  albert_lod_groups,
  function(dat) {
    dat = dat[order(dat$lod), , drop = FALSE]
    is_nonincreasing(dat$mean_same_direction_detection)
  },
  logical(1)
)
checks[[check_i]] = add_check(
  "albert_expected_detection_decreases_with_lod",
  all(albert_lod_monotone),
  paste(
    "monotone target-scope-scale sets =",
    sum(albert_lod_monotone), "/", length(albert_lod_monotone)
  )
)
check_i = check_i + 1L

scope_counts = reshape(
  unique(albert_detection[, c(
    "target",
    "reference_scope",
    "n_reference_qtl"
  )]),
  idvar = "target",
  timevar = "reference_scope",
  direction = "wide"
)
all_count_col = "n_reference_qtl.all_albert_qtls"
measured_count_col = "n_reference_qtl.target_measured_genes"
scope_count_ok =
  all(scope_counts[[all_count_col]] == nrow(albert_effects)) &&
  all(scope_counts[[measured_count_col]] > 0L) &&
  all(scope_counts[[measured_count_col]] <= scope_counts[[all_count_col]])

checks[[check_i]] = add_check(
  "target_measured_reference_is_only_a_measurability_subset",
  scope_count_ok,
  paste(
    "measured-effect count range =",
    paste(
      range(scope_counts[[measured_count_col]]),
      collapse = " to "
    ),
    "; total Albert effects =", nrow(albert_effects)
  )
)
check_i = check_i + 1L

########################################
# DIRECT PHENOTYPE SIMULATION          #
########################################

analytic_validation$any_error = abs(
  analytic_validation$analytic_any_detection -
  analytic_validation$simulated_any_detection
)
analytic_validation$same_error = abs(
  analytic_validation$analytic_same_direction -
  analytic_validation$simulated_same_direction
)

any_se = sqrt(
  pmax(
    0,
    analytic_validation$analytic_any_detection *
    (1 - analytic_validation$analytic_any_detection) /
    analytic_validation$n_rep
  )
)
same_se = sqrt(
  pmax(
    0,
    analytic_validation$analytic_same_direction *
    (1 - analytic_validation$analytic_same_direction) /
    analytic_validation$n_rep
  )
)
analytic_validation$any_tolerance = 0.005 + 4 * any_se
analytic_validation$same_tolerance = 0.005 + 4 * same_se

simulation_ok = all(
  analytic_validation$any_error <= analytic_validation$any_tolerance &
  analytic_validation$same_error <= analytic_validation$same_tolerance
)
max_sim_error = max(
  analytic_validation$any_error,
  analytic_validation$same_error,
  na.rm = TRUE
)
max_sim_tolerance = max(
  analytic_validation$any_tolerance,
  analytic_validation$same_tolerance,
  na.rm = TRUE
)

checks[[check_i]] = add_check(
  "direct_phenotype_simulation_supports_analytic_engine",
  simulation_ok,
  paste(
    "maximum absolute error =", signif(max_sim_error, 7),
    "; maximum tolerance =", signif(max_sim_tolerance, 7),
    "; settings =", nrow(analytic_validation)
  )
)
check_i = check_i + 1L

########################################
# MODEL-SEPARATION ASSERTIONS          #
########################################

model_map = setNames(
  model_definition$definition,
  model_definition$item
)
required_model_items = c(
  "primary_effect_grid",
  "albert_reference_role",
  "albert_ascertainment_caveat",
  "query_reascertainment",
  "target_qtl_calls_in_power_model",
  "positional_matching_in_power_model",
  "chance_overlap"
)
model_items_present = all(required_model_items %in% names(model_map))
model_separation_ok =
  model_items_present &&
  identical(model_map[["query_reascertainment"]], "None.") &&
  identical(model_map[["target_qtl_calls_in_power_model"]], "None.") &&
  grepl(
    "None; no confidence intervals",
    model_map[["positional_matching_in_power_model"]],
    fixed = TRUE
  ) &&
  grepl(
    "significance-selected",
    model_map[["albert_ascertainment_caveat"]],
    fixed = TRUE
  ) &&
  grepl(
    "not an assumed probability distribution",
    model_map[["primary_effect_grid"]],
    fixed = TRUE
  )

checks[[check_i]] = add_check(
  "model_cleanly_separates_power_from_qtl_matching",
  model_separation_ok,
  paste(
    "required model-definition items present =",
    model_items_present
  )
)
check_i = check_i + 1L

####################
# WRITE / FAIL FAST #
####################

validation_summary = do.call(rbind, checks)
write.table(
  validation_summary,
  file = validation_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nEFFECT-SIZE POWER ANALYSIS VALIDATION\n")
print(validation_summary, row.names = FALSE)
cat("\n")

failed = validation_summary$check[!validation_summary$passed]
if (length(failed) > 0L) {
  stop(
    "Effect-size power validation failed: ",
    paste(failed, collapse = ", ")
  )
}

message("All effect-size power checks passed.")
message("Validation summary: ", validation_out)
