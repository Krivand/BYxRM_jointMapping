################################################
################################################
## GENETICS REVIEW POWER ANALYSIS VALIDATOR   ##
## EFFECT-SIZE POWER CALIBRATION              ##
################################################
################################################

#################
# PRELIMINARIES #
#################

source("/PHShome/kv1026/p05/remap_functions.R")

wd = "/PHShome/kv1026/p05/"
out_dir = paste0(asys, "review/power_effect_curve/")
validation_out = paste0(out_dir, "effect_curve_validation_summary.tsv")
analysis_script = paste0(wd, "iReviewPowerEffectCurve.R")
validator_script = paste0(wd, "iReviewPowerEffectCurveValidate.R")

analysis_files = c(
  "critical_r_power_table.tsv",
  "albert_reference_crossstudy.tsv",
  "albert_2018_power_benchmark.tsv"
)

missing_files = analysis_files[
  !file.exists(paste0(out_dir, analysis_files))
]
if (length(missing_files) > 0L) {
  stop(
    "Missing reviewer-power outputs: ",
    paste(missing_files, collapse = ", ")
  )
}

validation_seed = suppressWarnings(
  as.integer(Sys.getenv("REVIEW_POWER_SEED", "20260818"))
)
if (is.na(validation_seed)) validation_seed = 20260818L

n_validation = suppressWarnings(
  as.integer(Sys.getenv("REVIEW_POWER_NVALID", "5000"))
)
if (is.na(n_validation) || n_validation < 1000L) {
  stop("REVIEW_POWER_NVALID must be at least 1000.")
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

load_single_rda = function(path, object_name = NULL) {

  if (!file.exists(path)) {
    stop("Missing required file: ", path)
  }

  env = new.env(parent = emptyenv())
  loaded = load(path, envir = env)

  if (!is.null(object_name)) {
    if (!(object_name %in% loaded)) {
      stop("Expected object '", object_name, "' in ", path)
    }
    return(env[[object_name]])
  }

  if (length(loaded) != 1L) {
    stop(
      "Expected one object in ", path,
      "; found: ", paste(loaded, collapse = ", ")
    )
  }

  env[[loaded[[1]]]]
}

add_check = function(name, passed, detail) {
  data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

independent_rcrit = function(lod, n) {
  sqrt(pmax(0, 1 - 10 ^ (-2 * lod / n)))
}

independent_probability = function(rho, n, lod) {

  rcrit = independent_rcrit(lod = lod, n = n)
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

validate_direct_simulation = function(
  n,
  rho,
  lod,
  n_rep,
  seed
) {

  set.seed(seed)

  x = rep(c(-1, 1), length.out = n)
  x = as.numeric(scale(x))
  residual_sd = sqrt(1 - rho ^ 2)

  noise = matrix(
    stats::rnorm(n * n_rep),
    nrow = n,
    ncol = n_rep
  )
  y = rho * x + residual_sd * noise
  y = sweep(y, 2, colMeans(y), "-")

  y_sd = sqrt(colSums(y ^ 2) / (n - 1))
  r_hat = drop(crossprod(x, y)) / ((n - 1) * y_sd)
  r_hat = pmax(-0.999999, pmin(0.999999, r_hat))
  lod_hat = -n * log10(1 - r_hat ^ 2) / 2

  analytic = independent_probability(
    rho = rho,
    n = n,
    lod = lod
  )

  data.frame(
    n = n,
    rho = rho,
    lod = lod,
    n_rep = n_rep,
    analytic_any = analytic[["any"]],
    simulated_any = mean(lod_hat > lod),
    analytic_same = analytic[["same"]],
    simulated_same = mean(lod_hat > lod & r_hat > 0),
    stringsAsFactors = FALSE
  )
}

################
# READ OUTPUTS #
################

critical = read_tsv(
  paste0(out_dir, "critical_r_power_table.tsv")
)
crossstudy = read_tsv(
  paste0(out_dir, "albert_reference_crossstudy.tsv")
)
benchmark = read_tsv(
  paste0(out_dir, "albert_2018_power_benchmark.tsv")
)

checks = list()
check_i = 1L

########################################
# OUTPUT SURFACE IS INTENTIONALLY SMALL #
########################################

existing_analysis_tables = setdiff(
  list.files(out_dir, pattern = "\\.tsv$"),
  "effect_curve_validation_summary.tsv"
)

checks[[check_i]] = add_check(
  "only_three_analysis_tables_are_emitted",
  setequal(existing_analysis_tables, analysis_files) &&
    length(existing_analysis_tables) == length(analysis_files),
  paste(
    "analysis tables =",
    paste(sort(existing_analysis_tables), collapse = ", ")
  )
)
check_i = check_i + 1L

########################################
# CRITICAL-R TABLE                     #
########################################

expected_lod = c(2, 3.03, 4.5, 7)
expected_studies = setdiff(classicMetadata$study, "AlbertAlt_2018")

checks[[check_i]] = add_check(
  "critical_table_has_exact_study_lod_grid",
  nrow(critical) == length(expected_studies) * length(expected_lod) &&
    setequal(unique(critical$study), expected_studies) &&
    isTRUE(all.equal(
      sort(unique(critical$lod)),
      sort(expected_lod),
      tolerance = 1e-14
    )),
  paste(
    "rows =", nrow(critical),
    "; studies =", length(unique(critical$study)),
    "; LOD values =", paste(sort(unique(critical$lod)), collapse = ", ")
  )
)
check_i = check_i + 1L

rcrit_expected = mapply(
  independent_rcrit,
  lod = critical$lod,
  n = critical$n
)
rcrit_error = max(abs(critical$critical_abs_r - rcrit_expected))
pve_error = max(
  abs(critical$critical_pve - critical$critical_abs_r ^ 2)
)

checks[[check_i]] = add_check(
  "critical_r_and_pve_are_exact",
  rcrit_error < 1e-12 && pve_error < 1e-12,
  paste(
    "maximum r error =", signif(rcrit_error, 7),
    "; maximum PVE error =", signif(pve_error, 7)
  )
)
check_i = check_i + 1L

power_specs = list(
  c("abs_r_20pct_power", "pve_20pct_power", 0.20),
  c("abs_r_50pct_power", "pve_50pct_power", 0.50),
  c("abs_r_80pct_power", "pve_80pct_power", 0.80)
)

power_errors = numeric(0)
pve_power_errors = numeric(0)
for (spec in power_specs) {
  r_col = spec[[1]]
  pve_col = spec[[2]]
  target_power = as.numeric(spec[[3]])

  recomputed = vapply(
    seq_len(nrow(critical)),
    function(i) {
      independent_probability(
        rho = critical[[r_col]][[i]],
        n = critical$n[[i]],
        lod = critical$lod[[i]]
      )[["any"]]
    },
    numeric(1)
  )

  power_errors = c(
    power_errors,
    max(abs(recomputed - target_power))
  )
  pve_power_errors = c(
    pve_power_errors,
    max(abs(critical[[pve_col]] - critical[[r_col]] ^ 2))
  )
}

checks[[check_i]] = add_check(
  "20_50_80pct_power_inversions_are_exact",
  max(power_errors) < 1e-8 && max(pve_power_errors) < 1e-12,
  paste(
    "maximum power error =", signif(max(power_errors), 7),
    "; maximum PVE error =", signif(max(pve_power_errors), 7)
  )
)
check_i = check_i + 1L

########################################
# INTERNAL POWER-SURFACE BEHAVIOR      #
########################################

rho_grid = seq(0, 0.8, by = 0.005)
lod_grid = seq(2, 7, by = 0.25)
unique_n = sort(unique(critical$n))

rho_monotone = logical(0)
lod_monotone = logical(0)
n_monotone = logical(0)

for (n in unique_n) {
  for (lod in lod_grid) {
    power = vapply(
      rho_grid,
      function(rho) independent_probability(rho, n, lod)[["any"]],
      numeric(1)
    )
    rho_monotone = c(rho_monotone, is_nondecreasing(power))
  }
}

for (n in unique_n) {
  for (rho in rho_grid[rho_grid >= 0.05]) {
    power = vapply(
      lod_grid,
      function(lod) independent_probability(rho, n, lod)[["any"]],
      numeric(1)
    )
    lod_monotone = c(lod_monotone, is_nonincreasing(power))
  }
}

for (lod in lod_grid) {
  for (rho in rho_grid[rho_grid >= 0.05]) {
    power = vapply(
      unique_n,
      function(n) independent_probability(rho, n, lod)[["any"]],
      numeric(1)
    )
    n_monotone = c(
      n_monotone,
      is_nondecreasing(power, tolerance = 1e-10)
    )
  }
}

checks[[check_i]] = add_check(
  "power_engine_has_expected_monotonicity",
  all(rho_monotone) && all(lod_monotone) && all(n_monotone),
  paste(
    "effect curves =", sum(rho_monotone), "/", length(rho_monotone),
    "; LOD curves =", sum(lod_monotone), "/", length(lod_monotone),
    "; n curves =", sum(n_monotone), "/", length(n_monotone)
  )
)
check_i = check_i + 1L

########################################
# ALBERT 2018 PUBLISHED BENCHMARK      #
########################################

ab_lod = 3.03
ab_total_n = 1012L
ab_n_per_group = 506L
ab_pve = 0.025
ab_r = sqrt(ab_pve)
ab_delta = 2 * sqrt(ab_pve / (1 - ab_pve))
ab_alpha = stats::pchisq(
  ab_lod * 2 * log(10),
  df = 1,
  lower.tail = FALSE
)
ab_public_power = stats::power.t.test(
  n = ab_n_per_group,
  delta = ab_delta,
  sd = 1,
  sig.level = ab_alpha,
  type = "two.sample",
  alternative = "two.sided"
)$power
ab_current = independent_probability(
  rho = ab_r,
  n = ab_total_n,
  lod = ab_lod
)

benchmark_error = max(
  abs(benchmark$total_n - ab_total_n),
  abs(benchmark$n_per_genotype_group_in_public_code - ab_n_per_group),
  abs(benchmark$lod_threshold_in_public_code - ab_lod),
  abs(benchmark$pve - ab_pve),
  abs(benchmark$equivalent_abs_r - ab_r),
  abs(benchmark$equivalent_two_group_delta - ab_delta),
  abs(benchmark$genomewide_sig_level_from_lod - ab_alpha),
  abs(benchmark$public_code_power_t_test - ab_public_power),
  abs(benchmark$current_engine_any_detection - ab_current[["any"]]),
  abs(benchmark$current_engine_same_direction - ab_current[["same"]])
)

checks[[check_i]] = add_check(
  "albert_2018_public_power_calculation_is_reproduced",
  nrow(benchmark) == 1L &&
    benchmark_error < 1e-12 &&
    ab_public_power > 0.90 &&
    ab_current[["any"]] > 0.90,
  paste(
    "public-code power =", signif(ab_public_power, 8),
    "; current engine =", signif(ab_current[["any"]], 8),
    "; maximum table error =", signif(benchmark_error, 7)
  )
)
check_i = check_i + 1L

ab_method_difference = abs(ab_current[["any"]] - ab_public_power)
checks[[check_i]] = add_check(
  "albert_2018_and_current_power_engines_agree",
  ab_method_difference < 0.002,
  paste(
    "absolute power difference at PVE 2.5% =",
    signif(ab_method_difference, 7)
  )
)
check_i = check_i + 1L

########################################
# ALBERT CROSS-STUDY TABLE             #
########################################

expected_targets = setdiff(expected_studies, "Albert_2018")
expected_scales = c(0.6, 0.8, 1.0)
expected_cross_rows =
  length(expected_targets) * length(expected_scales) * length(expected_lod)

checks[[check_i]] = add_check(
  "crossstudy_table_has_exact_target_scale_lod_grid",
  nrow(crossstudy) == expected_cross_rows &&
    setequal(unique(crossstudy$target), expected_targets) &&
    setequal(unique(crossstudy$effect_scale), expected_scales) &&
    isTRUE(all.equal(
      sort(unique(crossstudy$lod)),
      sort(expected_lod),
      tolerance = 1e-14
    )),
  paste(
    "rows =", nrow(crossstudy),
    "; expected =", expected_cross_rows
  )
)
check_i = check_i + 1L

studies = load_single_rda(
  paste0(imed, "studies.rda"),
  object_name = "studies"
)
frank_qtls_list = load_single_rda(
  paste0(imed, "Frank/qtls_list_260204.rda"),
  object_name = "qtls_list"
)

albert = frank_qtls_list[["Albert_2018"]]
if (is.null(albert)) {
  stop("Albert_2018 missing from Frank QTL object.")
}
albert_abs_r = abs(albert$cor_coef)

reference_summary_error = numeric(0)
expected_detection_error = numeric(0)

for (target in expected_targets) {

  target_genes = colnames(studies[[target]]$pheno)
  keep = albert$gene_ID %in% target_genes
  ref_r = albert_abs_r[keep]

  rows = crossstudy[crossstudy$target == target, , drop = FALSE]

  reference_summary_error = c(
    reference_summary_error,
    abs(rows$n_albert_qtl_on_measured_genes - length(ref_r)),
    abs(rows$reference_mean_abs_r - mean(ref_r)),
    abs(rows$reference_median_abs_r - stats::median(ref_r)),
    abs(rows$reference_mean_pve - mean(ref_r ^ 2)),
    abs(rows$reference_median_pve - stats::median(ref_r ^ 2))
  )

  for (i in seq_len(nrow(rows))) {
    scaled_r = pmin(0.999, ref_r * rows$effect_scale[[i]])
    probs = vapply(
      scaled_r,
      function(rho) {
        independent_probability(
          rho = rho,
          n = rows$target_n[[i]],
          lod = rows$lod[[i]]
        )[["any"]]
      },
      numeric(1)
    )
    probs_same = vapply(
      scaled_r,
      function(rho) {
        independent_probability(
          rho = rho,
          n = rows$target_n[[i]],
          lod = rows$lod[[i]]
        )[["same"]]
      },
      numeric(1)
    )

    expected_detection_error = c(
      expected_detection_error,
      abs(rows$expected_any_detection[[i]] - mean(probs)),
      abs(rows$expected_same_direction_detection[[i]] - mean(probs_same))
    )
  }
}

checks[[check_i]] = add_check(
  "crossstudy_reference_and_power_recompute_exactly",
  max(reference_summary_error) < 1e-12 &&
    max(expected_detection_error) < 1e-12,
  paste(
    "maximum reference-summary error =",
    signif(max(reference_summary_error), 7),
    "; maximum power error =",
    signif(max(expected_detection_error), 7)
  )
)
check_i = check_i + 1L

scale_groups = split(
  crossstudy,
  paste(crossstudy$target, crossstudy$lod, sep = "|")
)
scale_monotone = vapply(
  scale_groups,
  function(dat) {
    dat = dat[order(dat$effect_scale), , drop = FALSE]
    is_nondecreasing(dat$expected_any_detection)
  },
  logical(1)
)

lod_groups = split(
  crossstudy,
  paste(crossstudy$target, crossstudy$effect_scale, sep = "|")
)
lod_monotone = vapply(
  lod_groups,
  function(dat) {
    dat = dat[order(dat$lod), , drop = FALSE]
    is_nonincreasing(dat$expected_any_detection)
  },
  logical(1)
)

checks[[check_i]] = add_check(
  "crossstudy_sensitivities_are_monotone",
  all(scale_monotone) && all(lod_monotone),
  paste(
    "effect-scale sets =", sum(scale_monotone), "/", length(scale_monotone),
    "; LOD sets =", sum(lod_monotone), "/", length(lod_monotone)
  )
)
check_i = check_i + 1L

########################################
# DIRECT PHENOTYPE SIMULATION          #
########################################

validation_n = unique(c(
  min(critical$n),
  critical$n[which.min(abs(critical$n - stats::median(critical$n)))],
  max(critical$n)
))
validation_r = c(0.12, 0.16, 0.32)
validation_lod = c(2, 4.5, 7)
validation_grid = expand.grid(
  n = validation_n,
  rho = validation_r,
  lod = validation_lod,
  stringsAsFactors = FALSE
)

simulation = do.call(
  rbind,
  lapply(
    seq_len(nrow(validation_grid)),
    function(i) {
      validate_direct_simulation(
        n = validation_grid$n[[i]],
        rho = validation_grid$rho[[i]],
        lod = validation_grid$lod[[i]],
        n_rep = n_validation,
        seed = validation_seed + i
      )
    }
  )
)

simulation$any_error = abs(simulation$analytic_any - simulation$simulated_any)
simulation$same_error = abs(
  simulation$analytic_same - simulation$simulated_same
)
any_se = sqrt(
  pmax(
    0,
    simulation$analytic_any * (1 - simulation$analytic_any) /
      simulation$n_rep
  )
)
same_se = sqrt(
  pmax(
    0,
    simulation$analytic_same * (1 - simulation$analytic_same) /
      simulation$n_rep
  )
)
simulation$any_tolerance = 0.005 + 4 * any_se
simulation$same_tolerance = 0.005 + 4 * same_se

simulation_ok = all(
  simulation$any_error <= simulation$any_tolerance &
    simulation$same_error <= simulation$same_tolerance
)
max_sim_error = max(simulation$any_error, simulation$same_error)

checks[[check_i]] = add_check(
  "direct_phenotype_simulation_supports_power_engine",
  simulation_ok,
  paste(
    "maximum absolute error =", signif(max_sim_error, 7),
    "; settings =", nrow(simulation)
  )
)
check_i = check_i + 1L

########################################
# MODEL SEPARATION / CODE HYGIENE      #
########################################

analysis_text = paste(readLines(analysis_script, warn = FALSE), collapse = "\n")
validator_text = paste(readLines(validator_script, warn = FALSE), collapse = "\n")
forbidden_analysis_terms = c(
  "countOverlaps",
  "CI_min_base",
  "CI_max_base",
  "interval_padding",
  "observed_rediscovery",
  "query_qtl"
)
forbidden_present = vapply(
  forbidden_analysis_terms,
  function(term) grepl(term, analysis_text, fixed = TRUE),
  logical(1)
)

checks[[check_i]] = add_check(
  "power_analysis_contains_no_qtl_matching_machinery",
  !any(forbidden_present),
  paste(
    "forbidden terms present =",
    paste(forbidden_analysis_terms[forbidden_present], collapse = ", ")
  )
)
check_i = check_i + 1L

checks[[check_i]] = add_check(
  "review_power_scripts_follow_equals_style",
  !grepl("<-", analysis_text, fixed = TRUE) &&
    !grepl("<-", validator_text, fixed = TRUE),
  "No '<-' assignments detected in primary or validator script."
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

cat("\nSTREAMLINED EFFECT-SIZE POWER VALIDATION\n")
print(validation_summary, row.names = FALSE)
cat("\n")

failed = validation_summary$check[!validation_summary$passed]
if (length(failed) > 0L) {
  stop(
    "Streamlined reviewer-power validation failed: ",
    paste(failed, collapse = ", ")
  )
}

message("All streamlined reviewer-power checks passed.")
message("Validation summary: ", validation_out)
