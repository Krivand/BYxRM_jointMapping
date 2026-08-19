########################################
########################################
## GENETICS REVIEW POWER SIMULATION   ##
## INDIVIDUAL QTL REDISCOVERY         ##
########################################
########################################

#################
# PRELIMINARIES #
#################

source("/PHShome/kv1026/p05/remap_functions.R")

review_imed = paste0(
  imed,
  "review/"
)

review_asys = paste0(
  asys,
  "review/"
)

power_imed = paste0(
  review_imed,
  "power/"
)

power_asys = paste0(
  review_asys,
  "power/"
)

ensure_directory(review_imed)
ensure_directory(review_asys)
ensure_directory(power_imed)
ensure_directory(power_asys)

seed = suppressWarnings(
  as.integer(
    Sys.getenv(
      "REVIEW_POWER_SEED",
      "20260818"
    )
  )
)

if (is.na(seed)) {
  seed = 20260818L
}

n_synthetic = suppressWarnings(
  as.integer(
    Sys.getenv(
      "REVIEW_POWER_NQTL",
      "250000"
    )
  )
)

if (is.na(n_synthetic) || n_synthetic < 10000L) {
  stop(
    "REVIEW_POWER_NQTL must be at least 10000."
  )
}

n_validation = suppressWarnings(
  as.integer(
    Sys.getenv(
      "REVIEW_POWER_NVALID",
      "5000"
    )
  )
)

if (is.na(n_validation) || n_validation < 1000L) {
  stop(
    "REVIEW_POWER_NVALID must be at least 1000."
  )
}

message(
  "Synthetic QTL draws: ",
  format(
    n_synthetic,
    big.mark = ","
  )
)

message(
  "Direct-validation replicates per setting: ",
  format(
    n_validation,
    big.mark = ","
  )
)

####################
# HELPER FUNCTIONS #
####################

load_single_rda = function(
  path,
  object_name = NULL
) {

  if (!file.exists(path)) {
    stop(
      "Missing required file: ",
      path
    )
  }

  env = new.env(
    parent = emptyenv()
  )

  loaded = load(
    path,
    envir = env
  )

  if (!is.null(object_name)) {
    if (!(object_name %in% loaded)) {
      stop(
        "Expected object '",
        object_name,
        "' in ",
        path
      )
    }

    return(
      env[[object_name]]
    )
  }

  if (length(loaded) != 1L) {
    stop(
      "Expected one object in ",
      path,
      "; found: ",
      paste(
        loaded,
        collapse = ", "
      )
    )
  }

  env[[loaded[[1]]]]
}

lod_to_rcrit = function(
  lod,
  n
) {

  if (any(n <= 2)) {
    stop(
      "Sample sizes must exceed 2."
    )
  }

  sqrt(
    pmax(
      0,
      1 -
        10 ^ (
          -2 * lod / n
        )
    )
  )
}

correlation_call_prob = function(
  rho,
  n,
  lod
) {

  if (length(n) != 1L || n <= 2) {
    stop(
      "n must be a scalar greater than 2."
    )
  }

  if (length(lod) != 1L || !is.finite(lod)) {
    stop(
      "lod must be one finite scalar."
    )
  }

  rho = pmax(
    -0.999999,
    pmin(
      0.999999,
      rho
    )
  )

  rcrit = lod_to_rcrit(
    lod = lod,
    n = n
  )

  df = n - 2

  tcrit = rcrit *
    sqrt(
      df /
        (1 - rcrit ^ 2)
    )

  ncp = rho *
    sqrt(
      (n - 1) /
        (1 - rho ^ 2)
    )

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

  p_positive = pmax(
    0,
    pmin(
      1,
      p_positive
    )
  )

  p_negative = pmax(
    0,
    pmin(
      1,
      p_negative
    )
  )

  data.frame(
    rho = rho,
    p_positive = p_positive,
    p_negative = p_negative,
    power = p_positive + p_negative,
    stringsAsFactors = FALSE
  )
}

parse_fdr_thresholds = function(
  path
) {

  if (!file.exists(path)) {
    return(
      numeric(0)
    )
  }

  lines = readLines(
    path,
    warn = FALSE
  )

  hits = grep(
    "FDR Thresh 5%",
    lines,
    value = TRUE,
    fixed = TRUE
  )

  if (length(hits) == 0L) {
    return(
      numeric(0)
    )
  }

  values = suppressWarnings(
    as.numeric(
      sub(
        ".*FDR Thresh 5%[[:space:]]+([-+0-9.eE]+).*$",
        "\\1",
        hits
      )
    )
  )

  values[
    is.finite(values)
  ]
}

make_probability_lookup = function(
  study_design,
  lod,
  rho_grid = seq(
    -0.95,
    0.95,
    by = 0.0005
  )
) {

  out = lapply(
    seq_len(
      nrow(study_design)
    ),
    function(i) {

      probs = correlation_call_prob(
        rho = rho_grid,
        n = study_design$n[[i]],
        lod = lod
      )

      probs$study = study_design$study[[i]]
      probs
    }
  )

  names(out) = study_design$study
  out
}

probability_matrix = function(
  rho,
  lookup,
  component
) {

  if (!(component %in% c(
    "p_positive",
    "p_negative",
    "power"
  ))) {
    stop(
      "Unknown probability component: ",
      component
    )
  }

  out = vapply(
    names(lookup),
    function(study) {

      stats::approx(
        x = lookup[[study]]$rho,
        y = lookup[[study]][[component]],
        xout = rho,
        rule = 2,
        ties = "ordered"
      )$y
    },
    numeric(
      length(rho)
    )
  )

  colnames(out) = names(lookup)
  out
}

pair_table_from_components = function(
  numerator_same,
  numerator_any,
  denominator_query,
  study_design,
  scenario,
  effect_scale,
  lod,
  heterogeneity_mode = "none",
  heterogeneity_fraction = 0
) {

  study_names = study_design$study

  pairs = expand.grid(
    query = study_names,
    target = study_names,
    stringsAsFactors = FALSE
  )

  pairs = pairs[
    pairs$query != pairs$target,
    ,
    drop = FALSE
  ]

  query_index = match(
    pairs$query,
    study_names
  )

  target_index = match(
    pairs$target,
    study_names
  )

  pairs$expected_query_calls =
    denominator_query[
      query_index
    ]

  pairs$expected_both_detected =
    numerator_any[
      cbind(
        query_index,
        target_index
      )
    ]

  pairs$expected_both_same_direction =
    numerator_same[
      cbind(
        query_index,
        target_index
      )
    ]

  pairs$rediscovery_any_direction =
    pairs$expected_both_detected /
    pairs$expected_query_calls

  pairs$rediscovery_same_direction =
    pairs$expected_both_same_direction /
    pairs$expected_query_calls

  pairs$scenario = scenario
  pairs$effect_scale = effect_scale
  pairs$lod = lod
  pairs$heterogeneity_mode = heterogeneity_mode
  pairs$heterogeneity_fraction = heterogeneity_fraction

  metadata_index = match(
    pairs$query,
    study_design$study
  )

  pairs$query_analyte =
    study_design$analyte[
      metadata_index
    ]

  pairs$query_medium =
    study_design$medium[
      metadata_index
    ]

  pairs$query_n =
    study_design$n[
      metadata_index
    ]

  metadata_index = match(
    pairs$target,
    study_design$study
  )

  pairs$target_analyte =
    study_design$analyte[
      metadata_index
    ]

  pairs$target_medium =
    study_design$medium[
      metadata_index
    ]

  pairs$target_n =
    study_design$n[
      metadata_index
    ]

  pairs$comparison = paste0(
    pairs$query_analyte,
    "_to_",
    pairs$target_analyte
  )

  pairs$same_medium =
    pairs$query_medium ==
    pairs$target_medium

  pairs
}

expected_pairwise = function(
  rho_query,
  rho_target,
  lookup,
  study_design,
  scenario,
  effect_scale,
  lod,
  heterogeneity_mode = "none",
  heterogeneity_fraction = 0
) {

  query_positive = probability_matrix(
    rho = rho_query,
    lookup = lookup,
    component = "p_positive"
  )

  query_negative = probability_matrix(
    rho = rho_query,
    lookup = lookup,
    component = "p_negative"
  )

  target_positive = probability_matrix(
    rho = rho_target,
    lookup = lookup,
    component = "p_positive"
  )

  target_negative = probability_matrix(
    rho = rho_target,
    lookup = lookup,
    component = "p_negative"
  )

  query_detect =
    query_positive +
    query_negative

  target_detect =
    target_positive +
    target_negative

  numerator_same =
    crossprod(
      query_positive,
      target_positive
    ) +
    crossprod(
      query_negative,
      target_negative
    )

  numerator_any =
    crossprod(
      query_detect,
      target_detect
    )

  denominator_query =
    colSums(
      query_detect
    )

  pair_table_from_components(
    numerator_same = numerator_same,
    numerator_any = numerator_any,
    denominator_query = denominator_query,
    study_design = study_design,
    scenario = scenario,
    effect_scale = effect_scale,
    lod = lod,
    heterogeneity_mode = heterogeneity_mode,
    heterogeneity_fraction = heterogeneity_fraction
  )
}

validate_direct_simulation = function(
  n,
  rho,
  lod,
  n_rep,
  seed
) {

  set.seed(seed)

  x = rep(
    c(
      -1,
      1
    ),
    length.out = n
  )

  x = as.numeric(
    scale(x)
  )

  residual_sd = sqrt(
    1 - rho ^ 2
  )

  noise = matrix(
    stats::rnorm(
      n * n_rep
    ),
    nrow = n,
    ncol = n_rep
  )

  y = rho * x +
    residual_sd * noise

  y = sweep(
    y,
    2,
    colMeans(y),
    "-"
  )

  y_sd = sqrt(
    colSums(
      y ^ 2
    ) /
    (n - 1)
  )

  r_hat = drop(
    crossprod(
      x,
      y
    )
  ) /
    (
      (n - 1) *
      y_sd
    )

  r_hat = pmax(
    -0.999999,
    pmin(
      0.999999,
      r_hat
    )
  )

  lod_hat =
    -n *
    log(
      1 - r_hat ^ 2
    ) /
    (
      2 * log(10)
    )

  analytic = correlation_call_prob(
    rho = rho,
    n = n,
    lod = lod
  )

  data.frame(
    n = n,
    rho = rho,
    lod = lod,
    n_rep = n_rep,
    analytic_power = analytic$power,
    simulated_power = mean(
      lod_hat > lod
    ),
    analytic_same_direction = analytic$p_positive,
    simulated_same_direction = mean(
      lod_hat > lod &
      r_hat > 0
    ),
    stringsAsFactors = FALSE
  )
}

########################
# LOAD INDIVIDUAL DATA #
########################

studies = load_single_rda(
  paste0(
    imed,
    "studies.rda"
  ),
  object_name = "studies"
)

indivFullQTLs = load_single_rda(
  paste0(
    imed,
    "indivFullQTLs.rda"
  ),
  object_name = "indivFullQTLs"
)

# AlbertAlt_2018 is the same Albert & Bloom 2018 phenotype data
# ported into Grossbach marker coordinates for downstream joint analyses.
# The individual-remapping power analysis therefore uses Albert_2018 once.
analysis_studies = setdiff(
  classicMetadata$study,
  "AlbertAlt_2018"
)

missing_studies = setdiff(
  analysis_studies,
  names(studies)
)

if (length(missing_studies) > 0L) {
  stop(
    "Missing processed individual studies: ",
    paste(
      missing_studies,
      collapse = ", "
    )
  )
}

study_design = classicMetadata[
  match(
    analysis_studies,
    classicMetadata$study
  ),
  c(
    "study",
    "analyte",
    "medium"
  ),
  drop = FALSE
]

study_design$n = vapply(
  study_design$study,
  function(study) {
    nrow(
      studies[[study]]$pheno
    )
  },
  integer(1)
)

study_design$n_genes = vapply(
  study_design$study,
  function(study) {
    ncol(
      studies[[study]]$pheno
    )
  },
  integer(1)
)

qtl_counts = table(
  indivFullQTLs$study
)

study_design$n_qtl = as.integer(
  qtl_counts[
    study_design$study
  ]
)

study_design$n_qtl[
  is.na(study_design$n_qtl)
] = 0L

study_design$qtl_per_gene =
  study_design$n_qtl /
  study_design$n_genes

############################
# EMPIRICAL LOD CALIBRATION #
############################

threshold_log = paste0(
  logs,
  "iRemapAll.out"
)

empirical_thresholds = parse_fdr_thresholds(
  threshold_log
)

if (length(empirical_thresholds) == 0L) {
  warning(
    "No saved empirical FDR thresholds were recovered; using LOD 4.5 as the reference."
  )

  lod_reference = 4.5

  threshold_summary = data.frame(
    statistic = "fallback",
    lod = lod_reference,
    stringsAsFactors = FALSE
  )
} else {

  threshold_quantiles = stats::quantile(
    empirical_thresholds,
    probs = c(
      0,
      0.10,
      0.25,
      0.50,
      0.75,
      0.90,
      1
    ),
    na.rm = TRUE,
    names = FALSE
  )

  threshold_summary = data.frame(
    statistic = c(
      "minimum",
      "q10",
      "q25",
      "median",
      "q75",
      "q90",
      "maximum"
    ),
    lod = threshold_quantiles,
    stringsAsFactors = FALSE
  )

  lod_reference = stats::median(
    empirical_thresholds,
    na.rm = TRUE
  )
}

lod_override = suppressWarnings(
  as.numeric(
    Sys.getenv(
      "REVIEW_POWER_LOD",
      ""
    )
  )
)

if (is.finite(lod_override)) {
  lod_reference = lod_override
  message(
    "Using REVIEW_POWER_LOD override: ",
    lod_reference
  )
} else {
  message(
    "Reference LOD threshold: ",
    signif(
      lod_reference,
      4
    )
  )
}

study_design$reference_lod = lod_reference
study_design$reference_rcrit = lod_to_rcrit(
  lod = lod_reference,
  n = study_design$n
)

############################
# EMPIRICAL EFFECT-SIZE FIT #
############################

reference_study = "Albert_2018"

reference_qtls = indivFullQTLs[
  indivFullQTLs$study == reference_study,
  ,
  drop = FALSE
]

reference_abs_r = abs(
  reference_qtls$r
)

reference_abs_r = reference_abs_r[
  is.finite(reference_abs_r) &
  reference_abs_r > 0 &
  reference_abs_r < 0.999999
]

if (length(reference_abs_r) < 100L) {
  stop(
    "Too few Albert_2018 QTL effects available for empirical effect-size calibration."
  )
}

reference_fisher_z = atanh(
  reference_abs_r
)

log_fisher_z = log(
  reference_fisher_z
)

fit_meanlog = mean(
  log_fisher_z
)

fit_sdlog = stats::sd(
  log_fisher_z
)

if (!is.finite(fit_meanlog) ||
    !is.finite(fit_sdlog) ||
    fit_sdlog <= 0) {
  stop(
    "Failed to fit the Albert_2018 effect-size distribution."
  )
}

effect_fit = data.frame(
  reference_study = reference_study,
  n_reference_qtl = length(reference_abs_r),
  fit_family = "lognormal_on_absolute_fisher_z",
  meanlog = fit_meanlog,
  sdlog = fit_sdlog,
  stringsAsFactors = FALSE
)

set.seed(seed)

synthetic_fisher_z = stats::rlnorm(
  n_synthetic,
  meanlog = fit_meanlog,
  sdlog = fit_sdlog
)

synthetic_abs_r = tanh(
  synthetic_fisher_z
)

synthetic_abs_r = pmin(
  synthetic_abs_r,
  0.95
)

fit_quantiles = data.frame(
  quantile = c(
    0.01,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99
  ),
  observed_abs_r = as.numeric(
    stats::quantile(
      reference_abs_r,
      probs = c(
        0.01,
        0.05,
        0.10,
        0.25,
        0.50,
        0.75,
        0.90,
        0.95,
        0.99
      ),
      names = FALSE
    )
  ),
  fitted_abs_r = as.numeric(
    stats::quantile(
      synthetic_abs_r,
      probs = c(
        0.01,
        0.05,
        0.10,
        0.25,
        0.50,
        0.75,
        0.90,
        0.95,
        0.99
      ),
      names = FALSE
    )
  ),
  stringsAsFactors = FALSE
)

#########################
# GENERAL POWER CURVES  #
#########################

rho_grid = seq(
  0,
  0.70,
  by = 0.005
)

power_curves = do.call(
  rbind,
  lapply(
    seq_len(
      nrow(study_design)
    ),
    function(i) {

      probs = correlation_call_prob(
        rho = rho_grid,
        n = study_design$n[[i]],
        lod = lod_reference
      )

      data.frame(
        study = study_design$study[[i]],
        analyte = study_design$analyte[[i]],
        medium = study_design$medium[[i]],
        n = study_design$n[[i]],
        lod = lod_reference,
        rho = probs$rho,
        power = probs$power,
        p_correct_direction = probs$p_positive,
        p_wrong_direction = probs$p_negative,
        stringsAsFactors = FALSE
      )
    }
  )
)

if (length(empirical_thresholds) > 0L) {
  threshold_sensitivity_lod = as.numeric(
    stats::quantile(
      empirical_thresholds,
      probs = c(
        0.10,
        0.50,
        0.90
      ),
      names = FALSE
    )
  )

  threshold_sensitivity_label = c(
    "empirical_q10",
    "empirical_median",
    "empirical_q90"
  )
} else {
  threshold_sensitivity_lod = c(
    3.5,
    4.5,
    6
  )

  threshold_sensitivity_label = c(
    "lod_3.5",
    "lod_4.5",
    "lod_6"
  )
}

threshold_sensitivity = do.call(
  rbind,
  lapply(
    seq_along(
      threshold_sensitivity_lod
    ),
    function(k) {

      do.call(
        rbind,
        lapply(
          seq_len(
            nrow(study_design)
          ),
          function(i) {

            rho_values = seq(
              0.05,
              0.70,
              by = 0.05
            )

            probs = correlation_call_prob(
              rho = rho_values,
              n = study_design$n[[i]],
              lod = threshold_sensitivity_lod[[k]]
            )

            data.frame(
              threshold = threshold_sensitivity_label[[k]],
              lod = threshold_sensitivity_lod[[k]],
              study = study_design$study[[i]],
              n = study_design$n[[i]],
              rho = rho_values,
              power = probs$power,
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

##############################
# DIRECT PHENOTYPE VALIDATION #
##############################

validation_n = unique(
  c(
    min(
      study_design$n
    ),
    study_design$n[
      which.min(
        abs(
          study_design$n -
          stats::median(
            study_design$n
          )
        )
      )
    ],
    max(
      study_design$n
    )
  )
)

validation_rho = c(
  0.15,
  0.30,
  0.50
)

validation_grid = expand.grid(
  n = validation_n,
  rho = validation_rho,
  stringsAsFactors = FALSE
)

validation = do.call(
  rbind,
  lapply(
    seq_len(
      nrow(validation_grid)
    ),
    function(i) {

      validate_direct_simulation(
        n = validation_grid$n[[i]],
        rho = validation_grid$rho[[i]],
        lod = lod_reference,
        n_rep = n_validation,
        seed = seed + i
      )
    }
  )
)

validation$power_difference =
  validation$simulated_power -
  validation$analytic_power

validation$direction_difference =
  validation$simulated_same_direction -
  validation$analytic_same_direction

##################################
# SHARED-EFFECT REDISCOVERY      #
##################################

message(
  "Building probability lookup at reference LOD..."
)

reference_lookup = make_probability_lookup(
  study_design = study_design,
  lod = lod_reference
)

effect_scales = c(
  0.25,
  0.50,
  0.75,
  1.00
)

shared_rediscovery = do.call(
  rbind,
  lapply(
    effect_scales,
    function(effect_scale) {

      message(
        "Shared-effect rediscovery: effect scale ",
        effect_scale
      )

      rho = pmin(
        0.95,
        synthetic_abs_r *
        effect_scale
      )

      expected_pairwise(
        rho_query = rho,
        rho_target = rho,
        lookup = reference_lookup,
        study_design = study_design,
        scenario = "perfectly_shared",
        effect_scale = effect_scale,
        lod = lod_reference
      )
    }
  )
)

##################################
# HETEROGENEITY SENSITIVITY      #
##################################

base_rho = synthetic_abs_r

base_positive = probability_matrix(
  rho = base_rho,
  lookup = reference_lookup,
  component = "p_positive"
)

base_negative = probability_matrix(
  rho = base_rho,
  lookup = reference_lookup,
  component = "p_negative"
)

base_detect =
  base_positive +
  base_negative

base_same =
  crossprod(
    base_positive,
    base_positive
  ) +
  crossprod(
    base_negative,
    base_negative
  )

base_any =
  crossprod(
    base_detect,
    base_detect
  )

base_query =
  colSums(
    base_detect
  )

heterogeneity_modes = list(
  attenuated = 0.50 * base_rho,
  target_null = rep(
    0,
    length(base_rho)
  ),
  opposite = -base_rho
)

heterogeneity_fractions = c(
  0.05,
  0.10,
  0.25
)

heterogeneity_results = list()
heterogeneity_index = 1L

for (mode in names(heterogeneity_modes)) {

  message(
    "Heterogeneity mode: ",
    mode
  )

  alt_rho = heterogeneity_modes[[mode]]

  alt_positive = probability_matrix(
    rho = alt_rho,
    lookup = reference_lookup,
    component = "p_positive"
  )

  alt_negative = probability_matrix(
    rho = alt_rho,
    lookup = reference_lookup,
    component = "p_negative"
  )

  alt_detect =
    alt_positive +
    alt_negative

  alt_same =
    crossprod(
      base_positive,
      alt_positive
    ) +
    crossprod(
      base_negative,
      alt_negative
    )

  alt_any =
    crossprod(
      base_detect,
      alt_detect
    )

  for (fraction in heterogeneity_fractions) {

    mixed_same =
      (1 - fraction) *
      base_same +
      fraction *
      alt_same

    mixed_any =
      (1 - fraction) *
      base_any +
      fraction *
      alt_any

    heterogeneity_results[[heterogeneity_index]] =
      pair_table_from_components(
        numerator_same = mixed_same,
        numerator_any = mixed_any,
        denominator_query = base_query,
        study_design = study_design,
        scenario = "heterogeneity_sensitivity",
        effect_scale = 1,
        lod = lod_reference,
        heterogeneity_mode = mode,
        heterogeneity_fraction = fraction
      )

    heterogeneity_index =
      heterogeneity_index + 1L
  }
}

# A full target-null scenario gives the same-locus false-positive floor.
# It is not treated as the genome-wide interval-overlap chance expectation.
heterogeneity_results[[heterogeneity_index]] =
  pair_table_from_components(
    numerator_same = alt_same * 0,
    numerator_any = alt_any * 0,
    denominator_query = base_query,
    study_design = study_design,
    scenario = "placeholder",
    effect_scale = 1,
    lod = lod_reference,
    heterogeneity_mode = "placeholder",
    heterogeneity_fraction = 0
  )

# Replace the placeholder above with the actual all-target-null expectation.
null_rho = rep(
  0,
  length(base_rho)
)

null_positive = probability_matrix(
  rho = null_rho,
  lookup = reference_lookup,
  component = "p_positive"
)

null_negative = probability_matrix(
  rho = null_rho,
  lookup = reference_lookup,
  component = "p_negative"
)

null_detect =
  null_positive +
  null_negative

null_same =
  crossprod(
    base_positive,
    null_positive
  ) +
  crossprod(
    base_negative,
    null_negative
  )

null_any =
  crossprod(
    base_detect,
    null_detect
  )

heterogeneity_results[[heterogeneity_index]] =
  pair_table_from_components(
    numerator_same = null_same,
    numerator_any = null_any,
    denominator_query = base_query,
    study_design = study_design,
    scenario = "same_locus_null_floor",
    effect_scale = 1,
    lod = lod_reference,
    heterogeneity_mode = "target_null",
    heterogeneity_fraction = 1
  )

heterogeneity_rediscovery = do.call(
  rbind,
  heterogeneity_results
)

##################################
# THRESHOLD SENSITIVITY: PAIRS   #
##################################

threshold_pairwise = do.call(
  rbind,
  lapply(
    seq_along(
      threshold_sensitivity_lod
    ),
    function(k) {

      this_lod =
        threshold_sensitivity_lod[[k]]

      message(
        "Pairwise threshold sensitivity: ",
        threshold_sensitivity_label[[k]],
        " (LOD ",
        signif(
          this_lod,
          4
        ),
        ")"
      )

      lookup = make_probability_lookup(
        study_design = study_design,
        lod = this_lod
      )

      out = expected_pairwise(
        rho_query = base_rho,
        rho_target = base_rho,
        lookup = lookup,
        study_design = study_design,
        scenario = "threshold_sensitivity",
        effect_scale = 1,
        lod = this_lod
      )

      out$threshold =
        threshold_sensitivity_label[[k]]

      out
    }
  )
)

########################
# COMPACT SUMMARIES    #
########################

shared_summary = shared_rediscovery %>%
  dplyr::group_by(
    effect_scale,
    comparison
  ) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    mean_rediscovery = mean(
      rediscovery_same_direction
    ),
    median_rediscovery = stats::median(
      rediscovery_same_direction
    ),
    minimum_rediscovery = min(
      rediscovery_same_direction
    ),
    maximum_rediscovery = max(
      rediscovery_same_direction
    ),
    .groups = "drop"
  )

heterogeneity_summary = heterogeneity_rediscovery %>%
  dplyr::filter(
    scenario == "heterogeneity_sensitivity"
  ) %>%
  dplyr::group_by(
    heterogeneity_mode,
    heterogeneity_fraction,
    comparison
  ) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    mean_rediscovery = mean(
      rediscovery_same_direction
    ),
    median_rediscovery = stats::median(
      rediscovery_same_direction
    ),
    .groups = "drop"
  )

########################
# WRITE TABLES         #
########################

write.table(
  study_design,
  file = paste0(
    power_asys,
    "study_design.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  threshold_summary,
  file = paste0(
    power_asys,
    "empirical_threshold_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    lod = empirical_thresholds
  ),
  file = paste0(
    power_asys,
    "empirical_thresholds_all.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  effect_fit,
  file = paste0(
    power_asys,
    "effect_distribution_fit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  fit_quantiles,
  file = paste0(
    power_asys,
    "effect_distribution_quantiles.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  power_curves,
  file = paste0(
    power_asys,
    "study_power_curves.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  threshold_sensitivity,
  file = paste0(
    power_asys,
    "power_threshold_sensitivity.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  validation,
  file = paste0(
    power_asys,
    "analytic_power_validation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  shared_rediscovery,
  file = paste0(
    power_asys,
    "shared_effect_rediscovery.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  shared_summary,
  file = paste0(
    power_asys,
    "shared_effect_rediscovery_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  heterogeneity_rediscovery,
  file = paste0(
    power_asys,
    "heterogeneity_sensitivity.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  heterogeneity_summary,
  file = paste0(
    power_asys,
    "heterogeneity_sensitivity_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  threshold_pairwise,
  file = paste0(
    power_asys,
    "pairwise_threshold_sensitivity.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

review_power = list(
  settings = list(
    seed = seed,
    n_synthetic = n_synthetic,
    n_validation = n_validation,
    reference_study = reference_study,
    reference_lod = lod_reference,
    effect_scales = effect_scales,
    heterogeneity_fractions = heterogeneity_fractions
  ),
  study_design = study_design,
  empirical_thresholds = empirical_thresholds,
  threshold_summary = threshold_summary,
  effect_fit = effect_fit,
  effect_fit_quantiles = fit_quantiles,
  power_curves = power_curves,
  validation = validation,
  shared_rediscovery = shared_rediscovery,
  shared_summary = shared_summary,
  heterogeneity_rediscovery = heterogeneity_rediscovery,
  heterogeneity_summary = heterogeneity_summary,
  threshold_pairwise = threshold_pairwise
)

saveRDS(
  review_power,
  paste0(
    power_imed,
    "review_power.rds"
  )
)

########################
# DIAGNOSTIC FIGURES   #
########################

effect_plot_data = rbind(
  data.frame(
    abs_r = reference_abs_r,
    source = "Observed Albert 2018 QTLs",
    stringsAsFactors = FALSE
  ),
  data.frame(
    abs_r = synthetic_abs_r[
      seq_len(
        min(
          length(synthetic_abs_r),
          100000L
        )
      )
    ],
    source = "Fitted synthetic distribution",
    stringsAsFactors = FALSE
  )
)

effect_plot = ggplot2::ggplot(
  effect_plot_data,
  ggplot2::aes(
    x = abs_r,
    linetype = source
  )
) +
  ggplot2::geom_density(
    linewidth = 0.8
  ) +
  ggplot2::coord_cartesian(
    xlim = c(
      0,
      stats::quantile(
        reference_abs_r,
        0.995
      )
    )
  ) +
  ggplot2::labs(
    x = "Absolute standardized QTL effect (|r|)",
    y = "Density",
    linetype = NULL
  ) +
  ggplot2::theme_bw()

power_plot = ggplot2::ggplot(
  power_curves,
  ggplot2::aes(
    x = rho,
    y = power,
    group = study,
    linetype = analyte
  )
) +
  ggplot2::geom_line(
    linewidth = 0.7,
    alpha = 0.75
  ) +
  ggplot2::facet_wrap(
    ~ study
  ) +
  ggplot2::labs(
    x = "True standardized QTL effect (|r|)",
    y = "Probability of QTL detection",
    linetype = "Gene product"
  ) +
  ggplot2::theme_bw()

shared_heatmap_data = shared_rediscovery %>%
  dplyr::filter(
    effect_scale == 1
  )

shared_heatmap_data$query = factor(
  shared_heatmap_data$query,
  levels = rev(
    study_design$study
  )
)

shared_heatmap_data$target = factor(
  shared_heatmap_data$target,
  levels = study_design$study
)

shared_heatmap = ggplot2::ggplot(
  shared_heatmap_data,
  ggplot2::aes(
    x = target,
    y = query,
    fill = rediscovery_same_direction
  )
) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%.2f",
        rediscovery_same_direction
      )
    ),
    size = 2.5
  ) +
  ggplot2::scale_fill_viridis_c(
    limits = c(
      0,
      1
    ),
    name = "Expected\nrediscovery"
  ) +
  ggplot2::labs(
    x = "Target dataset",
    y = "Query dataset"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

heterogeneity_plot = ggplot2::ggplot(
  heterogeneity_summary,
  ggplot2::aes(
    x = heterogeneity_fraction,
    y = mean_rediscovery,
    linetype = comparison,
    group = comparison
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(
    ~ heterogeneity_mode
  ) +
  ggplot2::labs(
    x = "Fraction of true QTLs altered in target dataset",
    y = "Mean expected query-to-target rediscovery",
    linetype = "Comparison"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  filename = paste0(
    power_asys,
    "effect_distribution_fit.pdf"
  ),
  plot = effect_plot,
  width = 6,
  height = 4
)

ggplot2::ggsave(
  filename = paste0(
    power_asys,
    "effect_distribution_fit.png"
  ),
  plot = effect_plot,
  width = 6,
  height = 4,
  dpi = 300
)

ggplot2::ggsave(
  filename = paste0(
    power_asys,
    "study_power_curves.pdf"
  ),
  plot = power_plot,
  width = 10,
  height = 8
)

ggplot2::ggsave(
  filename = paste0(
    power_asys,
    "study_power_curves.png"
  ),
  plot = power_plot,
  width = 10,
  height = 8,
  dpi = 300
)

ggplot2::ggsave(
  filename = paste0(
    power_asys,
    "shared_effect_rediscovery_heatmap.pdf"
  ),
  plot = shared_heatmap,
  width = 10,
  height = 8
)

ggplot2::ggsave(
  filename = paste0(
    power_asys,
    "shared_effect_rediscovery_heatmap.png"
  ),
  plot = shared_heatmap,
  width = 10,
  height = 8,
  dpi = 300
)

ggplot2::ggsave(
  filename = paste0(
    power_asys,
    "heterogeneity_sensitivity.pdf"
  ),
  plot = heterogeneity_plot,
  width = 10,
  height = 5
)

ggplot2::ggsave(
  filename = paste0(
    power_asys,
    "heterogeneity_sensitivity.png"
  ),
  plot = heterogeneity_plot,
  width = 10,
  height = 5,
  dpi = 300
)

########################
# FINAL SUMMARY        #
########################

message("")
message("Reviewer power simulation complete.")
message(
  "Individual datasets analyzed: ",
  nrow(study_design)
)
message(
  "Reference effect distribution: ",
  reference_study,
  " (",
  format(
    length(reference_abs_r),
    big.mark = ","
  ),
  " mapped QTLs)"
)
message(
  "Reference LOD: ",
  signif(
    lod_reference,
    4
  )
)
message(
  "Outputs: ",
  power_asys
)
