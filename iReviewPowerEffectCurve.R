############################################
############################################
## GENETICS REVIEW POWER ANALYSIS         ##
## EFFECT-SIZE POWER CURVES               ##
############################################
############################################

#################
# PRELIMINARIES #
#################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
review_imed = paste0(imed, "review/")
power_asys = paste0(review_asys, "power_effect_curve/")
power_imed = paste0(review_imed, "power_effect_curve/")

ensure_directory(review_asys)
ensure_directory(review_imed)
ensure_directory(power_asys)
ensure_directory(power_imed)

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

rho_grid = seq(0, 0.80, by = 0.005)
lod_sweep = seq(2, 7, by = 0.25)
lod_display = 4.5
power_targets = c(0.20, 0.50, 0.80)
albert_scales = c(1.00, 0.80, 0.60)

message(
  "Effect-size power analysis: |rho| ",
  min(rho_grid), " to ", max(rho_grid),
  " by ", unique(diff(rho_grid))[1]
)
message(
  "LOD sensitivity: ",
  min(lod_sweep), " to ", max(lod_sweep),
  " by ", unique(diff(lod_sweep))[1]
)
message(
  "Direct-validation replicates per setting: ",
  format(n_validation, big.mark = ",")
)

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

write_table_gz = function(x, filename) {
  con = gzfile(paste0(power_asys, filename), open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(
    x,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

lod_to_rcrit = function(lod, n) {

  if (any(n <= 2)) {
    stop("Sample sizes must exceed 2.")
  }

  sqrt(
    pmax(
      0,
      1 - 10 ^ (-2 * lod / n)
    )
  )
}

correlation_detection_prob = function(rho, n, lod) {

  if (length(n) != 1L || !is.finite(n) || n <= 2) {
    stop("n must be one finite scalar greater than 2.")
  }
  if (length(lod) != 1L || !is.finite(lod)) {
    stop("lod must be one finite scalar.")
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
    rho = rho,
    p_same_direction = p_positive,
    p_wrong_direction = p_negative,
    p_any_detection = p_positive + p_negative,
    stringsAsFactors = FALSE
  )
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
    )$p_same_direction - target_power
  }

  lower = 0
  upper = 0.999
  f_lower = objective(lower)
  f_upper = objective(upper)

  if (f_lower >= 0) return(lower)
  if (f_upper <= 0) return(NA_real_)

  stats::uniroot(
    objective,
    interval = c(lower, upper),
    tol = 1e-11
  )$root
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

  analytic = correlation_detection_prob(
    rho = rho,
    n = n,
    lod = lod
  )

  data.frame(
    n = n,
    rho = rho,
    lod = lod,
    n_rep = n_rep,
    analytic_any_detection = analytic$p_any_detection,
    simulated_any_detection = mean(lod_hat > lod),
    analytic_same_direction = analytic$p_same_direction,
    simulated_same_direction = mean(lod_hat > lod & r_hat > 0),
    stringsAsFactors = FALSE
  )
}

########################
# LOAD INDIVIDUAL DATA #
########################

studies = load_single_rda(
  paste0(imed, "studies.rda"),
  object_name = "studies"
)
indivFullQTLs = load_single_rda(
  paste0(imed, "indivFullQTLs.rda"),
  object_name = "indivFullQTLs"
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
qtl_counts = table(indivFullQTLs$study)
study_design$n_qtl = as.integer(qtl_counts[study_design$study])
study_design$n_qtl[is.na(study_design$n_qtl)] = 0L

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

frank_qtl_count_check = data.frame(
  study = study_design$study,
  indivFullQTLs = study_design$n_qtl,
  frank_qtls_list = vapply(
    study_design$study,
    function(study) nrow(frank_individual_qtls[[study]]),
    integer(1)
  ),
  stringsAsFactors = FALSE
)
frank_qtl_count_check$match =
  frank_qtl_count_check$indivFullQTLs ==
  frank_qtl_count_check$frank_qtls_list

if (!all(frank_qtl_count_check$match)) {
  stop(
    "Frank/current QTL-count mismatch for: ",
    paste(
      frank_qtl_count_check$study[!frank_qtl_count_check$match],
      collapse = ", "
    )
  )
}

############################################
# PURE EFFECT-SIZE x SAMPLE-SIZE POWER GRID #
############################################

message("Calculating power surface without QTL matching or positional rules...")

power_surface = do.call(
  rbind,
  lapply(
    seq_len(nrow(study_design)),
    function(i) {
      study = study_design$study[[i]]
      n = study_design$n[[i]]

      do.call(
        rbind,
        lapply(
          lod_sweep,
          function(lod) {
            probs = correlation_detection_prob(
              rho = rho_grid,
              n = n,
              lod = lod
            )
            data.frame(
              study = study,
              analyte = study_design$analyte[[i]],
              medium = study_design$medium[[i]],
              n = n,
              lod = lod,
              abs_r = rho_grid,
              p_same_direction = probs$p_same_direction,
              p_wrong_direction = probs$p_wrong_direction,
              p_any_detection = probs$p_any_detection,
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)
rownames(power_surface) = NULL

############################################
# EFFECT MAGNITUDE FOR 20/50/80% POWER     #
############################################

power_thresholds = do.call(
  rbind,
  lapply(
    seq_len(nrow(study_design)),
    function(i) {
      study = study_design$study[[i]]
      n = study_design$n[[i]]

      do.call(
        rbind,
        lapply(
          lod_sweep,
          function(lod) {
            data.frame(
              study = study,
              analyte = study_design$analyte[[i]],
              medium = study_design$medium[[i]],
              n = n,
              lod = lod,
              target_power = power_targets,
              abs_r_required = vapply(
                power_targets,
                function(target_power) {
                  solve_r_for_power(
                    n = n,
                    lod = lod,
                    target_power = target_power
                  )
                },
                numeric(1)
              ),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)
rownames(power_thresholds) = NULL

############################################
# ALBERT EFFECT DISTRIBUTION AS REFERENCE  #
############################################

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

albert_reference_effects = data.frame(
  qtl_id = seq_len(nrow(albert_qtls)),
  gene_ID = albert_qtls$gene_ID,
  cor_coef = albert_qtls$cor_coef,
  abs_r = abs(albert_qtls$cor_coef),
  stringsAsFactors = FALSE
)

albert_quantiles = stats::quantile(
  albert_reference_effects$abs_r,
  probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
  names = FALSE,
  na.rm = TRUE
)
albert_effect_summary = data.frame(
  n_qtl = nrow(albert_reference_effects),
  mean_abs_r = mean(albert_reference_effects$abs_r),
  q10_abs_r = albert_quantiles[[1]],
  q25_abs_r = albert_quantiles[[2]],
  median_abs_r = albert_quantiles[[3]],
  q75_abs_r = albert_quantiles[[4]],
  q90_abs_r = albert_quantiles[[5]],
  max_abs_r = max(albert_reference_effects$abs_r),
  stringsAsFactors = FALSE
)

########################################################
# EXPECTED DETECTION OF ALBERT EFFECT MAGNITUDES       #
########################################################

# This is a separate empirical calibration layered on top of the pure power
# curves. Albert QTLs are explicitly an ascertained significant set, not an
# estimate of the unselected genome-wide effect distribution.
#
# Two scopes are exported:
#   all_albert_qtls       = the same empirical effect distribution for every
#                           target, answering a pure power question;
#   target_measured_genes = only Albert QTLs whose affected gene was measured
#                           in the target, providing a bridge to rediscovery
#                           without using target QTL calls or positional matching.

albert_targets = setdiff(study_design$study, "Albert_2018")

albert_reference_expected_detection = do.call(
  rbind,
  lapply(
    albert_targets,
    function(target) {

      target_n = study_design$n[match(target, study_design$study)]
      target_genes = colnames(studies[[target]]$pheno)
      target_measured = albert_reference_effects[
        albert_reference_effects$gene_ID %in% target_genes,
        ,
        drop = FALSE
      ]

      scopes = list(
        all_albert_qtls = albert_reference_effects,
        target_measured_genes = target_measured
      )

      do.call(
        rbind,
        lapply(
          names(scopes),
          function(reference_scope) {

            reference = scopes[[reference_scope]]
            if (nrow(reference) == 0L) {
              stop(
                "No Albert reference effects available for target ",
                target, " under scope ", reference_scope
              )
            }

            do.call(
              rbind,
              lapply(
                albert_scales,
                function(effect_scale) {

                  scaled_r = pmin(
                    0.999,
                    reference$abs_r * effect_scale
                  )

                  do.call(
                    rbind,
                    lapply(
                      lod_sweep,
                      function(lod) {
                        probs = correlation_detection_prob(
                          rho = scaled_r,
                          n = target_n,
                          lod = lod
                        )
                        data.frame(
                          target = target,
                          target_analyte = study_design$analyte[
                            match(target, study_design$study)
                          ],
                          target_n = target_n,
                          reference_scope = reference_scope,
                          n_reference_qtl = nrow(reference),
                          effect_scale = effect_scale,
                          lod = lod,
                          mean_abs_r = mean(scaled_r),
                          mean_same_direction_detection = mean(
                            probs$p_same_direction
                          ),
                          mean_any_detection = mean(
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
    }
  )
)
rownames(albert_reference_expected_detection) = NULL

############################################
# DIRECT VALIDATION OF ANALYTIC ENGINE     #
############################################

n_values = sort(unique(study_design$n))
median_n = stats::median(n_values)
validation_n_values = unique(c(
  min(n_values),
  n_values[which.min(abs(n_values - median_n))],
  max(n_values)
))

validation_r_values = as.numeric(
  stats::quantile(
    albert_reference_effects$abs_r,
    probs = c(0.10, 0.50, 0.90),
    names = FALSE,
    na.rm = TRUE
  )
)
validation_lod_values = c(2, lod_display, 7)

validation_grid = expand.grid(
  n = validation_n_values,
  rho = validation_r_values,
  lod = validation_lod_values,
  stringsAsFactors = FALSE
)

analytic_power_validation = do.call(
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
rownames(analytic_power_validation) = NULL

############################################
# MODEL DEFINITION / PROVENANCE            #
############################################

model_definition = data.frame(
  item = c(
    "primary_question",
    "primary_effect_grid",
    "primary_probability",
    "lod_sensitivity",
    "albert_reference_role",
    "albert_ascertainment_caveat",
    "albert_scale_sensitivity",
    "target_measured_scope",
    "query_reascertainment",
    "target_qtl_calls_in_power_model",
    "positional_matching_in_power_model",
    "chance_overlap"
  ),
  definition = c(
    paste0(
      "For a stipulated true absolute genotype-phenotype correlation and a ",
      "study sample size, calculate the probability that the study exceeds ",
      "a candidate LOD threshold in the correct direction."
    ),
    paste0(
      "Equally spaced absolute correlations from ", min(rho_grid),
      " to ", max(rho_grid), " by ", unique(diff(rho_grid))[1],
      "; this grid is not an assumed probability distribution."
    ),
    "Noncentral-t approximation corresponding to the mapper LOD transformation.",
    paste0(
      "Common candidate target LOD from ", min(lod_sweep),
      " to ", max(lod_sweep), " by ", unique(diff(lod_sweep))[1], "."
    ),
    paste0(
      "Observed significant Albert_2018 QTL effect magnitudes are used only ",
      "as an empirical reference distribution layered onto the primary power curves."
    ),
    paste0(
      "Albert_2018 QTL effects are significance-selected and are not treated ",
      "as an unbiased estimate of the latent genome-wide effect distribution."
    ),
    paste0(
      "Reference magnitudes are also multiplied by ",
      paste(albert_scales, collapse = ", "),
      " as a transparent sensitivity analysis for effect inflation."
    ),
    paste0(
      "A secondary target-specific summary restricts the Albert reference to ",
      "genes measured in each target; it does not use target QTL calls."
    ),
    "None.",
    "None.",
    "None; no confidence intervals, genomic overlap, or padding enter the power calculation.",
    "Positional chance overlap remains a separate analysis and is not part of this power model."
  ),
  stringsAsFactors = FALSE
)

#################
# WRITE OUTPUTS #
#################

write_table(study_design, "study_design.tsv")
write_table(frank_qtl_count_check, "frank_qtl_count_check.tsv")
write_table_gz(power_surface, "effect_size_power_surface.tsv.gz")
write_table(power_thresholds, "effect_size_power_thresholds.tsv")
write_table_gz(
  albert_reference_effects,
  "albert_reference_effects.tsv.gz"
)
write_table(albert_effect_summary, "albert_effect_summary.tsv")
write_table(
  albert_reference_expected_detection,
  "albert_reference_expected_detection.tsv"
)
write_table(
  analytic_power_validation,
  "analytic_power_validation.tsv"
)
write_table(model_definition, "model_definition.tsv")

review_power_effect_curve = list(
  settings = list(
    validation_seed = validation_seed,
    n_validation = n_validation,
    rho_grid = rho_grid,
    lod_sweep = lod_sweep,
    lod_display = lod_display,
    power_targets = power_targets,
    albert_scales = albert_scales
  ),
  study_design = study_design,
  frank_qtl_count_check = frank_qtl_count_check,
  power_surface = power_surface,
  power_thresholds = power_thresholds,
  albert_reference_effects = albert_reference_effects,
  albert_effect_summary = albert_effect_summary,
  albert_reference_expected_detection =
    albert_reference_expected_detection,
  analytic_power_validation = analytic_power_validation,
  model_definition = model_definition
)

saveRDS(
  review_power_effect_curve,
  paste0(power_imed, "review_power_effect_curve.rds")
)

message("")
message("Effect-size power analysis complete.")
message(
  "Power surface: ",
  paste0(power_asys, "effect_size_power_surface.tsv.gz")
)
message(
  "20/50/80% thresholds: ",
  paste0(power_asys, "effect_size_power_thresholds.tsv")
)
message(
  "Albert calibration: ",
  paste0(power_asys, "albert_reference_expected_detection.tsv")
)
message(
  "Validation grid: ",
  paste0(power_asys, "analytic_power_validation.tsv")
)
