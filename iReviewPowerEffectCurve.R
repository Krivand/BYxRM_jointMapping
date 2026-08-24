############################################
############################################
## GENETICS REVIEW POWER ANALYSIS         ##
## EFFECT-SIZE POWER CALIBRATION          ##
############################################
############################################

#################
# PRELIMINARIES #
#################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
power_asys = paste0(review_asys, "power_effect_curve/")

ensure_directory(review_asys)
ensure_directory(power_asys)

# These are the only reviewer-power tables produced by this script.
output_files = c(
  "critical_r_power_table.tsv",
  "albert_reference_crossstudy.tsv",
  "albert_2018_power_benchmark.tsv"
)

# Remove superseded outputs from earlier iterations in this same directory.
deprecated_outputs = c(
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
unlink(paste0(power_asys, deprecated_outputs), force = TRUE)

lod_report = c(2, 3.03, 4.5, 7)
power_targets = c(0.20, 0.50, 0.80)
albert_scales = c(1.00, 0.80, 0.60)
engine_grid_pve = seq(0.005, 0.10, by = 0.001)

power_cores = suppressWarnings(
  as.integer(Sys.getenv("REVIEW_POWER_CORES", "1"))
)
if (is.na(power_cores) || power_cores < 1L) power_cores = 1L

message(
  "Reviewer power calibration: LOD values ",
  paste(lod_report, collapse = ", ")
)
message(
  "Albert-reference effect scales: ",
  paste(albert_scales, collapse = ", ")
)
message("Parallel worker limit: ", power_cores)

####################
# HELPER FUNCTIONS #
####################

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

write_table = function(x, filename) {
  write.table(
    x,
    file = paste0(power_asys, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

lod_to_rcrit = function(lod, n) {

  if (any(!is.finite(n)) || any(n <= 2)) {
    stop("Sample sizes must be finite and exceed 2.")
  }
  if (any(!is.finite(lod)) || any(lod <= 0)) {
    stop("LOD values must be finite and positive.")
  }

  sqrt(pmax(0, 1 - 10 ^ (-2 * lod / n)))
}

correlation_detection_prob = function(rho, n, lod) {

  if (length(n) != 1L || !is.finite(n) || n <= 2) {
    stop("n must be one finite scalar greater than 2.")
  }
  if (length(lod) != 1L || !is.finite(lod) || lod <= 0) {
    stop("lod must be one finite positive scalar.")
  }
  if (any(!is.finite(rho))) {
    stop("rho contains non-finite values.")
  }
  if (any(rho < 0 | rho >= 1)) {
    stop("This effect-magnitude engine requires 0 <= rho < 1.")
  }

  rcrit = lod_to_rcrit(lod = lod, n = n)
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

  p_positive[!is.finite(p_positive)] = 0
  p_negative[!is.finite(p_negative)] = 0
  p_positive = pmax(0, pmin(1, p_positive))
  p_negative = pmax(0, pmin(1, p_negative))

  data.frame(
    abs_r = rho,
    p_same_direction = p_positive,
    p_wrong_direction = p_negative,
    p_any_detection = p_positive + p_negative,
    stringsAsFactors = FALSE
  )
}

parallel_lapply = function(X, FUN) {
  if (
    .Platform$OS.type == "unix" &&
    power_cores > 1L &&
    length(X) > 1L
  ) {
    return(parallel::mclapply(
      X,
      FUN,
      mc.cores = min(power_cores, length(X)),
      mc.preschedule = TRUE
    ))
  }
  lapply(X, FUN)
}

solve_r_for_power = function(n, lod, target_power) {

  if (!is.finite(target_power) || target_power <= 0 || target_power >= 1) {
    stop("target_power must lie strictly between 0 and 1.")
  }

  objective = function(rho) {
    correlation_detection_prob(
      rho = rho,
      n = n,
      lod = lod
    )$p_any_detection - target_power
  }

  lower = 0
  upper = 0.999

  if (objective(lower) >= 0) return(lower)
  if (objective(upper) <= 0) return(NA_real_)

  stats::uniroot(
    objective,
    interval = c(lower, upper),
    tol = 1e-11
  )$root
}

########################
# LOAD INDIVIDUAL DATA #
########################

studies = load_single_rda(
  paste0(imed, "studies.rda"),
  object_name = "studies"
)

analysis_studies = setdiff(
  classicMetadata$study,
  "AlbertAlt_2018"
)

missing_studies = setdiff(analysis_studies, names(studies))
if (length(missing_studies) > 0L) {
  stop(
    "Missing processed individual studies: ",
    paste(missing_studies, collapse = ", ")
  )
}

study_design = classicMetadata[
  match(analysis_studies, classicMetadata$study),
  c("study", "analyte", "medium"),
  drop = FALSE
]
study_design$n = vapply(
  study_design$study,
  function(study) nrow(studies[[study]]$pheno),
  integer(1)
)
study_design$n_genes = vapply(
  study_design$study,
  function(study) ncol(studies[[study]]$pheno),
  integer(1)
)

frank_qtls_list = load_single_rda(
  paste0(imed, "Frank/qtls_list_260204.rda"),
  object_name = "qtls_list"
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

missing_frank = setdiff(unname(frank_name_map), names(frank_qtls_list))
if (length(missing_frank) > 0L) {
  stop(
    "Missing individual datasets from Frank QTL object: ",
    paste(missing_frank, collapse = ", ")
  )
}

frank_individual_qtls = lapply(
  frank_name_map,
  function(name) frank_qtls_list[[name]]
)
names(frank_individual_qtls) = names(frank_name_map)

########################################
# TABLE 1: CRITICAL R AND POWER        #
########################################

message("Building critical-r and detectable-effect table...")

critical_r_power_table = do.call(
  rbind,
  lapply(
    seq_len(nrow(study_design)),
    function(i) {

      study = study_design$study[[i]]
      n = study_design$n[[i]]

      do.call(
        rbind,
        lapply(
          lod_report,
          function(lod) {

            r_required = vapply(
              power_targets,
              function(target_power) {
                solve_r_for_power(
                  n = n,
                  lod = lod,
                  target_power = target_power
                )
              },
              numeric(1)
            )

            data.frame(
              study = study,
              analyte = study_design$analyte[[i]],
              medium = study_design$medium[[i]],
              n = n,
              n_genes = study_design$n_genes[[i]],
              lod = lod,
              critical_abs_r = lod_to_rcrit(lod = lod, n = n),
              critical_pve = lod_to_rcrit(lod = lod, n = n) ^ 2,
              abs_r_20pct_power = r_required[[1]],
              pve_20pct_power = r_required[[1]] ^ 2,
              abs_r_50pct_power = r_required[[2]],
              pve_50pct_power = r_required[[2]] ^ 2,
              abs_r_80pct_power = r_required[[3]],
              pve_80pct_power = r_required[[3]] ^ 2,
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)
rownames(critical_r_power_table) = NULL

################################################
# TABLE 2: ALBERT EFFECTS THROUGH EACH TARGET #
################################################

message("Building Albert-reference cross-study power table...")

albert_qtls = frank_individual_qtls[["Albert_2018"]]
required_albert = c("gene_ID", "cor_coef")
missing_albert = setdiff(required_albert, colnames(albert_qtls))
if (length(missing_albert) > 0L) {
  stop(
    "Albert QTL table missing required columns: ",
    paste(missing_albert, collapse = ", ")
  )
}

invalid_albert =
  !is.finite(albert_qtls$cor_coef) |
  albert_qtls$cor_coef == 0 |
  abs(albert_qtls$cor_coef) >= 1
if (any(invalid_albert)) {
  stop(
    "Albert QTL table contains ", sum(invalid_albert),
    " invalid effect estimates."
  )
}

albert_reference = data.frame(
  gene_ID = albert_qtls$gene_ID,
  abs_r = abs(albert_qtls$cor_coef),
  stringsAsFactors = FALSE
)

albert_targets = setdiff(study_design$study, "Albert_2018")

albert_reference_crossstudy = do.call(
  rbind,
  parallel_lapply(
    albert_targets,
    function(target) {

      target_i = match(target, study_design$study)
      target_n = study_design$n[[target_i]]
      target_genes = colnames(studies[[target]]$pheno)

      reference = albert_reference[
        albert_reference$gene_ID %in% target_genes,
        ,
        drop = FALSE
      ]

      if (nrow(reference) == 0L) {
        stop("No Albert QTL effects remain for target ", target)
      }

      reference_mean_abs_r = mean(reference$abs_r)
      reference_median_abs_r = stats::median(reference$abs_r)
      reference_mean_pve = mean(reference$abs_r ^ 2)
      reference_median_pve = stats::median(reference$abs_r ^ 2)

      do.call(
        rbind,
        lapply(
          albert_scales,
          function(effect_scale) {

            scaled_r = pmin(0.999, reference$abs_r * effect_scale)

            do.call(
              rbind,
              lapply(
                lod_report,
                function(lod) {

                  probs = correlation_detection_prob(
                    rho = scaled_r,
                    n = target_n,
                    lod = lod
                  )

                  data.frame(
                    target = target,
                    target_analyte = study_design$analyte[[target_i]],
                    target_medium = study_design$medium[[target_i]],
                    target_n = target_n,
                    target_n_genes = study_design$n_genes[[target_i]],
                    n_albert_qtl_on_measured_genes = nrow(reference),
                    reference_mean_abs_r = reference_mean_abs_r,
                    reference_median_abs_r = reference_median_abs_r,
                    reference_mean_pve = reference_mean_pve,
                    reference_median_pve = reference_median_pve,
                    effect_scale = effect_scale,
                    lod = lod,
                    expected_same_direction_detection = mean(
                      probs$p_same_direction
                    ),
                    expected_any_detection = mean(
                      probs$p_any_detection
                    ),
                    stringsAsFactors = FALSE
                  )
                }
              )
            )
          }
        )
      )
    }
  )
)
rownames(albert_reference_crossstudy) = NULL

# Keep the full target x scale x LOD grid in memory for validation, but emit
# only the compact table needed for interpretation and reviewer-facing use.
extract_crossstudy_power = function(effect_scale, lod, column) {
  rows = albert_reference_crossstudy[
    albert_reference_crossstudy$effect_scale == effect_scale &
      albert_reference_crossstudy$lod == lod,
    ,
    drop = FALSE
  ]
  rows = rows[match(albert_targets, rows$target), , drop = FALSE]
  rows[[column]]
}

base_crossstudy = albert_reference_crossstudy[
  albert_reference_crossstudy$effect_scale == 1 &
    albert_reference_crossstudy$lod == 3.03,
  c(
    "target",
    "target_analyte",
    "target_medium",
    "target_n",
    "target_n_genes",
    "n_albert_qtl_on_measured_genes",
    "reference_mean_abs_r",
    "reference_median_abs_r",
    "reference_mean_pve",
    "reference_median_pve"
  ),
  drop = FALSE
]
base_crossstudy = base_crossstudy[
  match(albert_targets, base_crossstudy$target),
  ,
  drop = FALSE
]

albert_reference_crossstudy_table = base_crossstudy
albert_reference_crossstudy_table$expected_detection_lod2 =
  extract_crossstudy_power(1, 2, "expected_any_detection")
albert_reference_crossstudy_table$expected_detection_lod3_03 =
  extract_crossstudy_power(1, 3.03, "expected_any_detection")
albert_reference_crossstudy_table$expected_detection_lod4_5 =
  extract_crossstudy_power(1, 4.5, "expected_any_detection")
albert_reference_crossstudy_table$expected_detection_lod7 =
  extract_crossstudy_power(1, 7, "expected_any_detection")
albert_reference_crossstudy_table$expected_detection_lod4_5_scale0_8 =
  extract_crossstudy_power(0.8, 4.5, "expected_any_detection")
albert_reference_crossstudy_table$expected_detection_lod4_5_scale0_6 =
  extract_crossstudy_power(0.6, 4.5, "expected_any_detection")
albert_reference_crossstudy_table$expected_same_direction_lod4_5 =
  extract_crossstudy_power(1, 4.5, "expected_same_direction_detection")
rownames(albert_reference_crossstudy_table) = NULL

################################################
# TABLE 3: ALBERT ET AL. 2018 POWER BENCHMARK #
################################################

message("Reproducing the Albert et al. 2018 Figure 1C power benchmark...")

# Albert et al. 2018 public code used:
#   LOD threshold = 3.03
#   sig.level = pchisq(3.03 * 2 * log(10), 1, lower.tail = FALSE)
#   power.t.test(n = 506, delta = ..., sig.level = sig.level)
# and reported >90% power for an eQTL explaining >=2.5% of variance.
# For a balanced RI cross with residual variance 1,
#   PVE = a^2 / (a^2 + 1),
# where a is half the homozygote mean difference. Therefore the two-group
# difference passed to power.t.test is delta = 2 * sqrt(PVE / (1 - PVE)).

ab2018_lod = 3.03
ab2018_total_n = 1012L
ab2018_n_per_group = 506L
ab2018_pve = 0.025
ab2018_abs_r = sqrt(ab2018_pve)
ab2018_delta = 2 * sqrt(ab2018_pve / (1 - ab2018_pve))
ab2018_sig_level = stats::pchisq(
  ab2018_lod * 2 * log(10),
  df = 1,
  lower.tail = FALSE
)

ab2018_power_t_test = stats::power.t.test(
  n = ab2018_n_per_group,
  delta = ab2018_delta,
  sd = 1,
  sig.level = ab2018_sig_level,
  type = "two.sample",
  alternative = "two.sided"
)$power

our_ab2018_power = correlation_detection_prob(
  rho = ab2018_abs_r,
  n = ab2018_total_n,
  lod = ab2018_lod
)

# Compare the two analytic engines across a range of effect sizes, not only
# at the published 2.5%-PVE benchmark. This is deterministic; stochastic
# phenotype simulation remains an independent validator-only check.
ab2018_public_grid_power = vapply(
  engine_grid_pve,
  function(pve) {
    delta = 2 * sqrt(pve / (1 - pve))
    stats::power.t.test(
      n = ab2018_n_per_group,
      delta = delta,
      sd = 1,
      sig.level = ab2018_sig_level,
      type = "two.sample",
      alternative = "two.sided"
    )$power
  },
  numeric(1)
)
ab2018_current_grid_power = vapply(
  sqrt(engine_grid_pve),
  function(rho) {
    correlation_detection_prob(
      rho = rho,
      n = ab2018_total_n,
      lod = ab2018_lod
    )$p_any_detection
  },
  numeric(1)
)
ab2018_grid_difference =
  ab2018_current_grid_power - ab2018_public_grid_power

albert_2018_power_benchmark = data.frame(
  source = "Albert_et_al_2018_Figure_1C",
  published_statement = ">90% power for PVE >= 0.025",
  total_n = ab2018_total_n,
  n_per_genotype_group_in_public_code = ab2018_n_per_group,
  lod_threshold_in_public_code = ab2018_lod,
  pve = ab2018_pve,
  equivalent_abs_r = ab2018_abs_r,
  equivalent_two_group_delta = ab2018_delta,
  genomewide_sig_level_from_lod = ab2018_sig_level,
  public_code_power_t_test = ab2018_power_t_test,
  current_engine_any_detection = our_ab2018_power$p_any_detection,
  current_engine_same_direction = our_ab2018_power$p_same_direction,
  current_minus_public_code =
    our_ab2018_power$p_any_detection - ab2018_power_t_test,
  engine_grid_pve_min = min(engine_grid_pve),
  engine_grid_pve_max = max(engine_grid_pve),
  engine_grid_pve_step = unique(diff(engine_grid_pve))[[1]],
  engine_grid_n = length(engine_grid_pve),
  engine_grid_max_abs_difference = max(abs(ab2018_grid_difference)),
  engine_grid_mean_abs_difference = mean(abs(ab2018_grid_difference)),
  engine_grid_rmse = sqrt(mean(ab2018_grid_difference ^ 2)),
  agrees_with_published_gt90 =
    ab2018_power_t_test > 0.90 &&
    our_ab2018_power$p_any_detection > 0.90,
  stringsAsFactors = FALSE
)

#################
# WRITE OUTPUTS #
#################

write_table(
  critical_r_power_table,
  "critical_r_power_table.tsv"
)
write_table(
  albert_reference_crossstudy_table,
  "albert_reference_crossstudy.tsv"
)
write_table(
  albert_2018_power_benchmark,
  "albert_2018_power_benchmark.tsv"
)

message("")
message("Reviewer power calibration complete.")
message("Produced exactly three analysis tables:")
for (filename in output_files) {
  message("  ", paste0(power_asys, filename))
}
