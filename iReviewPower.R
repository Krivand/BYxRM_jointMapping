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

review_imed = paste0(imed, "review/")
review_asys = paste0(asys, "review/")
power_imed = paste0(review_imed, "power/")
power_asys = paste0(review_asys, "power/")

ensure_directory(review_imed)
ensure_directory(review_asys)
ensure_directory(power_imed)
ensure_directory(power_asys)

seed = suppressWarnings(
  as.integer(Sys.getenv("REVIEW_POWER_SEED", "20260818"))
)
if (is.na(seed)) seed = 20260818L

n_synthetic = suppressWarnings(
  as.integer(Sys.getenv("REVIEW_POWER_NQTL", "250000"))
)
if (is.na(n_synthetic) || n_synthetic < 10000L) {
  stop("REVIEW_POWER_NQTL must be at least 10000.")
}

n_validation = suppressWarnings(
  as.integer(Sys.getenv("REVIEW_POWER_NVALID", "5000"))
)
if (is.na(n_validation) || n_validation < 1000L) {
  stop("REVIEW_POWER_NVALID must be at least 1000.")
}

message(
  "Synthetic QTL draws: ",
  format(n_synthetic, big.mark = ",")
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
      stop(
        "Expected object '", object_name,
        "' in ", path
      )
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

correlation_call_prob = function(rho, n, lod) {

  if (length(n) != 1L || n <= 2) {
    stop("n must be a scalar greater than 2.")
  }
  if (length(lod) != 1L || !is.finite(lod)) {
    stop("lod must be one finite scalar.")
  }

  rho = pmax(-0.999999, pmin(0.999999, rho))
  rcrit = lod_to_rcrit(lod = lod, n = n)
  df = n - 2
  tcrit = rcrit * sqrt(df / (1 - rcrit ^ 2))

  # This noncentrality parameter closely matches direct simulation of
  # standardized phenotype-genotype correlations under the mapper's LOD rule.
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
    p_positive = p_positive,
    p_negative = p_negative,
    power = p_positive + p_negative,
    stringsAsFactors = FALSE
  )
}

parse_fdr_thresholds = function(path) {

  if (!file.exists(path)) return(numeric(0))

  lines = readLines(path, warn = FALSE)
  hits = grep(
    "FDR Thresh 5%",
    lines,
    value = TRUE,
    fixed = TRUE
  )

  if (length(hits) == 0L) return(numeric(0))

  values = suppressWarnings(
    as.numeric(
      sub(
        ".*FDR Thresh 5%[[:space:]]+([-+0-9.eE]+).*$",
        "\\1",
        hits
      )
    )
  )

  values[is.finite(values)]
}

make_probability_lookup = function(
  study_design,
  lod,
  rho_grid = seq(-0.95, 0.95, by = 0.0005)
) {

  out = lapply(
    seq_len(nrow(study_design)),
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

probability_matrix = function(rho, lookup, component) {

  if (!(component %in% c(
    "p_positive",
    "p_negative",
    "power"
  ))) {
    stop("Unknown probability component: ", component)
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
    numeric(length(rho))
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

  query_index = match(pairs$query, study_names)
  target_index = match(pairs$target, study_names)

  pairs$expected_query_calls = denominator_query[query_index]
  pairs$expected_both_detected = numerator_any[
    cbind(query_index, target_index)
  ]
  pairs$expected_both_same_direction = numerator_same[
    cbind(query_index, target_index)
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

  metadata_index = match(pairs$query, study_design$study)
  pairs$query_analyte = study_design$analyte[metadata_index]
  pairs$query_medium = study_design$medium[metadata_index]
  pairs$query_n = study_design$n[metadata_index]

  metadata_index = match(pairs$target, study_design$study)
  pairs$target_analyte = study_design$analyte[metadata_index]
  pairs$target_medium = study_design$medium[metadata_index]
  pairs$target_n = study_design$n[metadata_index]

  pairs$comparison = paste0(
    pairs$query_analyte,
    "_to_",
    pairs$target_analyte
  )
  pairs$same_medium =
    pairs$query_medium == pairs$target_medium

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

  query_detect = query_positive + query_negative
  target_detect = target_positive + target_negative

  numerator_same =
    crossprod(query_positive, target_positive) +
    crossprod(query_negative, target_negative)
  numerator_any = crossprod(query_detect, target_detect)
  denominator_query = colSums(query_detect)

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
  lod_hat = -n * log(1 - r_hat ^ 2) / (2 * log(10))

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
    simulated_power = mean(lod_hat > lod),
    analytic_same_direction = analytic$p_positive,
    simulated_same_direction = mean(
      lod_hat > lod & r_hat > 0
    ),
    stringsAsFactors = FALSE
  )
}

fit_observed_overlap_glmm = function(dat) {

  model_data = dat[
    is.finite(dat$observed_rediscovery) &
    dat$n_query_eligible > 0L,
    ,
    drop = FALSE
  ]

  model_data$misses =
    model_data$n_query_eligible -
    model_data$n_matches

  model_data$comparison =
    factor(
      model_data$comparison,
      levels = c(
        "mRNA_to_mRNA",
        "mRNA_to_protein",
        "protein_to_mRNA",
        "protein_to_protein"
      )
    )

  model_data$query =
    factor(model_data$query)

  model_data$target =
    factor(model_data$target)

  fit = lme4::glmer(
    cbind(
      n_matches,
      misses
    ) ~
      0 +
      comparison +
      (1 | query) +
      (1 | target),
    family = stats::binomial(),
    data = model_data,
    control = lme4::glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 200000
      )
    )
  )

  coef_tab =
    summary(fit)$coefficients

  estimate =
    coef_tab[, "Estimate"]

  se =
    coef_tab[, "Std. Error"]

  comparison =
    sub(
      "^comparison",
      "",
      rownames(coef_tab)
    )

  data.frame(
    comparison = comparison,
    estimate_logit = estimate,
    se_logit = se,
    adjusted_agreement = stats::plogis(
      estimate
    ),
    adjusted_lower_95 = stats::plogis(
      estimate -
      1.96 * se
    ),
    adjusted_upper_95 = stats::plogis(
      estimate +
      1.96 * se
    ),
    singular_fit = lme4::isSingular(
      fit,
      tol = 1e-4
    ),
    n_pairs =
      as.integer(
        table(
          model_data$comparison
        )[
          comparison
        ]
      ),
    stringsAsFactors = FALSE
  )
}

safe_cor = function(x, y, method = "pearson") {
  keep = is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) return(NA_real_)
  stats::cor(x[keep], y[keep], method = method)
}

summarize_observed_expected = function(dat, grouping) {

  dat %>%
    dplyr::filter(
      is.finite(observed_rediscovery),
      is.finite(expected_rediscovery)
    ) %>%
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(grouping)
      )
    ) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      mean_observed = mean(observed_rediscovery),
      mean_expected = mean(expected_rediscovery),
      pearson = safe_cor(
        observed_rediscovery,
        expected_rediscovery,
        method = "pearson"
      ),
      spearman = safe_cor(
        observed_rediscovery,
        expected_rediscovery,
        method = "spearman"
      ),
      mae = mean(abs(
        observed_rediscovery -
        expected_rediscovery
      )),
      rmse = sqrt(mean((
        observed_rediscovery -
        expected_rediscovery
      ) ^ 2)),
      .groups = "drop"
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

# AlbertAlt_2018 is the same Albert & Bloom phenotype dataset placed in
# Grossbach marker coordinates for downstream joint analyses. It is not an
# independent experiment and is excluded from individual-study power work.
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
study_design$n_qtl = as.integer(
  qtl_counts[study_design$study]
)
study_design$n_qtl[is.na(study_design$n_qtl)] = 0L
study_design$qtl_per_gene =
  study_design$n_qtl /
  study_design$n_genes

############################
# EMPIRICAL LOD CALIBRATION #
############################

threshold_log = paste0(logs, "iRemapAll.out")
empirical_thresholds = parse_fdr_thresholds(
  threshold_log
)

if (length(empirical_thresholds) == 0L) {

  warning(
    paste0(
      "No saved empirical FDR thresholds were recovered; ",
      "using LOD 4.5 as the reference."
    )
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
    probs = c(0, 0.10, 0.25, 0.50, 0.75, 0.90, 1),
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
  as.numeric(Sys.getenv("REVIEW_POWER_LOD", ""))
)
if (is.finite(lod_override)) {
  lod_reference = lod_override
  message("Using REVIEW_POWER_LOD override: ", lod_reference)
} else {
  message(
    "Reference LOD threshold: ",
    signif(lod_reference, 4)
  )
}

study_design$reference_lod = lod_reference
study_design$reference_rcrit = lod_to_rcrit(
  lod = lod_reference,
  n = study_design$n
)

################################
# EMPIRICAL EFFECT CALIBRATION #
################################

reference_study = "Albert_2018"
reference_qtls = indivFullQTLs[
  indivFullQTLs$study == reference_study,
  ,
  drop = FALSE
]
reference_abs_r = abs(reference_qtls$r)
reference_abs_r = reference_abs_r[
  is.finite(reference_abs_r) &
  reference_abs_r > 0 &
  reference_abs_r < 0.999999
]

if (length(reference_abs_r) < 100L) {
  stop(
    paste0(
      "Too few Albert_2018 QTL effects available for ",
      "empirical effect-size calibration."
    )
  )
}

# Resampling, rather than fitting a parametric distribution, preserves the
# shape and tails of the effect sizes seen in the best-powered individual
# study. Downward effect scaling below is retained as a sensitivity analysis
# because discovered effects are themselves subject to ascertainment.
effect_fit = data.frame(
  reference_study = reference_study,
  n_reference_qtl = length(reference_abs_r),
  fit_family = "empirical_resampling_absolute_r",
  sampling = "with_replacement",
  stringsAsFactors = FALSE
)

set.seed(seed)
synthetic_abs_r = sample(
  reference_abs_r,
  size = n_synthetic,
  replace = TRUE
)
synthetic_abs_r = pmin(synthetic_abs_r, 0.95)

quantile_probs = c(
  0.01,
  0.05,
  0.10,
  0.25,
  0.50,
  0.75,
  0.90,
  0.95,
  0.99
)
fit_quantiles = data.frame(
  quantile = quantile_probs,
  observed_abs_r = as.numeric(
    stats::quantile(
      reference_abs_r,
      probs = quantile_probs,
      names = FALSE
    )
  ),
  resampled_abs_r = as.numeric(
    stats::quantile(
      synthetic_abs_r,
      probs = quantile_probs,
      names = FALSE
    )
  ),
  stringsAsFactors = FALSE
)

############################
# OBSERVED PAIRWISE OVERLAP #
############################


########################################
# REPRODUCE FRANK FIGURE 3 OVERLAPS    #
########################################

normalize_yeast_chr = function(x) {

  x = toupper(
    trimws(
      as.character(x)
    )
  )

  x = sub(
    "^CHR",
    "",
    x
  )

  roman = c(
    "I", "II", "III", "IV",
    "V", "VI", "VII", "VIII",
    "IX", "X", "XI", "XII",
    "XIII", "XIV", "XV", "XVI"
  )

  idx = match(
    x,
    roman
  )

  use = !is.na(idx)

  x[use] =
    as.character(
      idx[use]
    )

  x
}

make_frank_qtl_granges = function(
  qtls,
  pad = 10000L
) {

  required = c(
    "gene_ID",
    "CI_chr",
    "CI_min_base",
    "CI_max_base",
    "cor_coef"
  )

  missing = setdiff(
    required,
    colnames(qtls)
  )

  if (length(missing) > 0L) {
    stop(
      "Frank QTL table missing: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }

  left = suppressWarnings(
    as.numeric(
      qtls$CI_min_base
    )
  )

  right = suppressWarnings(
    as.numeric(
      qtls$CI_max_base
    )
  )

  bad =
    !is.finite(left) |
    !is.finite(right)

  if (any(bad)) {
    stop(
      "Non-finite Frank QTL confidence intervals encountered."
    )
  }

  seq_key = paste(
    normalize_yeast_chr(
      qtls$CI_chr
    ),
    qtls$gene_ID,
    sep = "_"
  )

  strand = ifelse(
    qtls$cor_coef > 0,
    "+",
    ifelse(
      qtls$cor_coef < 0,
      "-",
      "*"
    )
  )

  gr = GenomicRanges::GRanges(
    seqnames = seq_key,
    ranges = IRanges::IRanges(
      start = pmax(
        1,
        pmin(
          left,
          right
        ) -
        as.integer(pad)
      ),
      end =
        pmax(
          left,
          right
        ) +
        as.integer(pad)
    ),
    strand = strand
  )

  S4Vectors::mcols(gr)$gene =
    qtls$gene_ID

  S4Vectors::mcols(gr)$r =
    qtls$cor_coef

  gr
}

make_technical_marker_ranges = function(chr_features) {

  required = c(
    "Feature_Name",
    "Chromosome",
    "Start_Coordinate",
    "Stop_Coordinate"
  )

  missing = setdiff(
    required,
    colnames(chr_features)
  )

  if (length(missing) > 0L) {
    stop(
      paste0(
        "chrFeatures missing columns required for ",
        "Frank technical-marker exclusions: ",
        paste(
          missing,
          collapse = ", "
        )
      )
    )
  }

  tech_orf = c(
    LEU2 = "YCL018W",
    MAT = "YCR039C",
    URA3 = "YEL021W",
    AMN1 = "YBR158W",
    KRE33 = "YNL132W"
  )

  idx = match(
    unname(tech_orf),
    chr_features$Feature_Name
  )

  if (any(is.na(idx))) {
    stop(
      "Technical-marker ORFs missing from chrFeatures: ",
      paste(
        names(tech_orf)[is.na(idx)],
        collapse = ", "
      )
    )
  }

  left = suppressWarnings(
    as.numeric(
      chr_features$Start_Coordinate[idx]
    )
  )

  right = suppressWarnings(
    as.numeric(
      chr_features$Stop_Coordinate[idx]
    )
  )

  gr = GenomicRanges::GRanges(
    seqnames = normalize_yeast_chr(
      chr_features$Chromosome[idx]
    ),
    ranges = IRanges::IRanges(
      start = pmin(
        left,
        right
      ),
      end = pmax(
        left,
        right
      )
    )
  )

  names(gr) =
    names(tech_orf)

  gr
}

frank_technical_markers = function() {

  standard = c(
    "LEU2",
    "URA3",
    "AMN1",
    "MAT"
  )

  list(
    Brem_2002 = standard,
    Brem_2005 = standard,
    SmithEth_2008 = standard,
    SmithGlc_2008 = standard,
    Khan_2009 = standard,
    Foss_2011 = standard,
    Picotti_2013 = standard,
    GrossbachRNA_2022 = standard,
    GrossbachPRT_2022 = standard,
    Albert_2018 = c(
      "MAT",
      "KRE33"
    )
  )
}

build_observed_pairwise_frank = function(
  qtl_by_study,
  studies,
  study_design,
  technical_marker_ranges,
  pad = 10000L
) {

  tech_list =
    frank_technical_markers()

  pair_grid = expand.grid(
    query = study_design$study,
    target = study_design$study,
    stringsAsFactors = FALSE
  )

  pair_grid = pair_grid[
    pair_grid$query !=
    pair_grid$target,
    ,
    drop = FALSE
  ]

  out = lapply(
    seq_len(
      nrow(pair_grid)
    ),
    function(i) {

      query =
        pair_grid$query[[i]]

      target =
        pair_grid$target[[i]]

      query_genes =
        colnames(
          studies[[query]]$pheno
        )

      target_genes =
        colnames(
          studies[[target]]$pheno
        )

      query_qtls =
        qtl_by_study[[query]]

      target_qtls =
        qtl_by_study[[target]]

      # Frank datasetOverlap(): query QTLs must affect genes
      # that were measured in the target dataset.
      query_qtls =
        query_qtls[
          query_qtls$gene_ID %in%
          target_genes,
          ,
          drop = FALSE
        ]

      n_before_technical =
        nrow(query_qtls)

      # Frank datasetOverlap(): remove technical-marker loci
      # present in the query cross but absent from the target.
      exclude_markers = setdiff(
        tech_list[[query]],
        tech_list[[target]]
      )

      n_technical_excluded = 0L

      if (
        nrow(query_qtls) > 0L &&
        length(exclude_markers) > 0L
      ) {

        left = suppressWarnings(
          as.numeric(
            query_qtls$CI_min_base
          )
        )

        right = suppressWarnings(
          as.numeric(
            query_qtls$CI_max_base
          )
        )

        query_chr_gr =
          GenomicRanges::GRanges(
            seqnames =
              normalize_yeast_chr(
                query_qtls$CI_chr
              ),
            ranges = IRanges::IRanges(
              start = pmax(
                1,
                pmin(
                  left,
                  right
                ) -
                as.integer(pad)
              ),
              end =
                pmax(
                  left,
                  right
                ) +
                as.integer(pad)
            )
          )

        technical_overlap =
          suppressWarnings(
            GenomicRanges::countOverlaps(
              query_chr_gr,
              technical_marker_ranges[
                exclude_markers
              ],
              ignore.strand = TRUE
            )
          ) > 0L

        n_technical_excluded =
          sum(
            technical_overlap
          )

        query_qtls =
          query_qtls[
            !technical_overlap,
            ,
            drop = FALSE
          ]
      }

      n_query =
        nrow(query_qtls)

      if (
        n_query == 0L ||
        nrow(target_qtls) == 0L
      ) {

        n_match = 0L

      } else {

        query_gr =
          make_frank_qtl_granges(
            query_qtls,
            pad = pad
          )

        target_gr =
          make_frank_qtl_granges(
            target_qtls,
            pad = pad
          )

        # Exact Figure 3 overlap statistic:
        # chromosome + gene must match and strand preserves
        # the QTL effect direction.
        overlap_count =
          suppressWarnings(
            GenomicRanges::countOverlaps(
              query_gr,
              target_gr,
              ignore.strand = FALSE
            )
          )

        n_match =
          sum(
            overlap_count > 0L
          )
      }

      data.frame(
        query = query,
        target = target,
        n_common_genes =
          length(
            intersect(
              query_genes,
              target_genes
            )
          ),
        n_query_before_technical =
          n_before_technical,
        n_technical_excluded =
          n_technical_excluded,
        n_query_eligible =
          n_query,
        n_target_eligible =
          sum(
            target_qtls$gene_ID %in%
            query_genes
          ),
        n_matches =
          n_match,
        observed_rediscovery =
          if (n_query > 0L) {
            n_match / n_query
          } else {
            NA_real_
          },
        stringsAsFactors = FALSE
      )
    }
  )

  observed =
    do.call(
      rbind,
      out
    )

  idx = match(
    observed$query,
    study_design$study
  )

  observed$query_analyte =
    study_design$analyte[idx]

  observed$query_medium =
    study_design$medium[idx]

  observed$query_n =
    study_design$n[idx]

  idx = match(
    observed$target,
    study_design$study
  )

  observed$target_analyte =
    study_design$analyte[idx]

  observed$target_medium =
    study_design$medium[idx]

  observed$target_n =
    study_design$n[idx]

  observed$comparison = paste0(
    observed$query_analyte,
    "_to_",
    observed$target_analyte
  )

  observed$same_medium =
    observed$query_medium ==
    observed$target_medium

  observed$interval_padding_bp =
    as.integer(pad)

  observed$overlap_definition =
    "Frank_Figure3_overlapsWithStrand"

  observed
}

frank_qtls_list = load_single_rda(
  paste0(
    imed,
    "Frank/qtls_list_260204.rda"
  ),
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

missing_frank = setdiff(
  unname(
    frank_name_map
  ),
  names(
    frank_qtls_list
  )
)

if (length(missing_frank) > 0L) {
  stop(
    "Missing individual datasets from Frank QTL object: ",
    paste(
      missing_frank,
      collapse = ", "
    )
  )
}

frank_individual_qtls = lapply(
  frank_name_map,
  function(name) {
    frank_qtls_list[[name]]
  }
)

names(
  frank_individual_qtls
) =
  names(
    frank_name_map
  )

chrFeatures = load_single_rda(
  paste0(
    raw,
    "chrFeatures.rda"
  ),
  object_name = "chrFeatures"
)

technical_marker_ranges =
  make_technical_marker_ranges(
    chrFeatures
  )

frank_qtl_count_check = data.frame(
  study =
    study_design$study,
  indivFullQTLs =
    study_design$n_qtl,
  frank_qtls_list =
    vapply(
      study_design$study,
      function(study) {
        nrow(
          frank_individual_qtls[[study]]
        )
      },
      integer(1)
    ),
  stringsAsFactors = FALSE
)

frank_qtl_count_check$match =
  frank_qtl_count_check$indivFullQTLs ==
  frank_qtl_count_check$frank_qtls_list

if (!all(
  frank_qtl_count_check$match
)) {
  stop(
    "Frank/current QTL-count mismatch for: ",
    paste(
      frank_qtl_count_check$study[
        !frank_qtl_count_check$match
      ],
      collapse = ", "
    )
  )
}

write.table(
  frank_qtl_count_check,
  file = paste0(
    power_asys,
    "frank_qtl_count_check.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "Frank QTL object matches current remap QTL counts for all 10 datasets."
)

message(
  "Recomputing observed individual-study rediscovery using Frank Figure 3 logic..."
)

observed_pairwise =
  build_observed_pairwise_frank(
    qtl_by_study =
      frank_individual_qtls,
    studies =
      studies,
    study_design =
      study_design,
    technical_marker_ranges =
      technical_marker_ranges,
    pad = 10000L
  )

message(
  "Fitting manuscript-style observed-overlap GLMM..."
)

observed_glmm_summary =
  fit_observed_overlap_glmm(
    observed_pairwise
  )

#########################
# GENERAL POWER CURVES  #
#########################

rho_grid = seq(0, 0.70, by = 0.005)
power_curves = do.call(
  rbind,
  lapply(
    seq_len(nrow(study_design)),
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
      probs = c(0.10, 0.50, 0.90),
      names = FALSE
    )
  )
  threshold_sensitivity_label = c(
    "empirical_q10",
    "empirical_median",
    "empirical_q90"
  )
} else {
  threshold_sensitivity_lod = c(3.5, 4.5, 6)
  threshold_sensitivity_label = c(
    "lod_3.5",
    "lod_4.5",
    "lod_6"
  )
}

threshold_sensitivity = do.call(
  rbind,
  lapply(
    seq_along(threshold_sensitivity_lod),
    function(k) {
      do.call(
        rbind,
        lapply(
          seq_len(nrow(study_design)),
          function(i) {
            rho_values = seq(0.05, 0.70, by = 0.05)
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
    min(study_design$n),
    study_design$n[
      which.min(
        abs(
          study_design$n -
          stats::median(study_design$n)
        )
      )
    ],
    max(study_design$n)
  )
)
validation_rho = c(0.15, 0.30, 0.50)
validation_grid = expand.grid(
  n = validation_n,
  rho = validation_rho,
  stringsAsFactors = FALSE
)

validation = do.call(
  rbind,
  lapply(
    seq_len(nrow(validation_grid)),
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

effect_scales = round(seq(0.20, 1.00, by = 0.05), 2)
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
        synthetic_abs_r * effect_scale
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
# OBSERVED VS EXPECTED           #
##################################

observed_for_join = observed_pairwise[
  ,
  c(
    "query",
    "target",
    "n_common_genes",
    "n_query_eligible",
    "n_target_eligible",
    "n_matches",
    "observed_rediscovery"
  ),
  drop = FALSE
]

observed_vs_expected = merge(
  shared_rediscovery,
  observed_for_join,
  by = c("query", "target"),
  all.x = TRUE,
  sort = FALSE
)
observed_vs_expected$expected_rediscovery =
  observed_vs_expected$rediscovery_same_direction
observed_vs_expected$difference =
  observed_vs_expected$observed_rediscovery -
  observed_vs_expected$expected_rediscovery
observed_vs_expected$absolute_error = abs(
  observed_vs_expected$difference
)
observed_vs_expected$squared_error =
  observed_vs_expected$difference ^ 2

observed_vs_expected_summary =
  summarize_observed_expected(
    observed_vs_expected,
    grouping = "effect_scale"
  )

observed_vs_expected_by_comparison =
  summarize_observed_expected(
    observed_vs_expected,
    grouping = c(
      "effect_scale",
      "comparison"
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
base_detect = base_positive + base_negative
base_same =
  crossprod(base_positive, base_positive) +
  crossprod(base_negative, base_negative)
base_any = crossprod(base_detect, base_detect)
base_query = colSums(base_detect)

heterogeneity_modes = list(
  attenuated = 0.50 * base_rho,
  target_null = rep(0, length(base_rho)),
  opposite = -base_rho
)
heterogeneity_fractions = c(0.05, 0.10, 0.25)
heterogeneity_results = list()
heterogeneity_index = 1L

for (mode in names(heterogeneity_modes)) {

  message("Heterogeneity mode: ", mode)
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
  alt_detect = alt_positive + alt_negative

  alt_same =
    crossprod(base_positive, alt_positive) +
    crossprod(base_negative, alt_negative)
  alt_any = crossprod(base_detect, alt_detect)

  for (fraction in heterogeneity_fractions) {
    mixed_same =
      (1 - fraction) * base_same +
      fraction * alt_same
    mixed_any =
      (1 - fraction) * base_any +
      fraction * alt_any

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

    heterogeneity_index = heterogeneity_index + 1L
  }
}

null_rho = rep(0, length(base_rho))
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
null_detect = null_positive + null_negative
null_same =
  crossprod(base_positive, null_positive) +
  crossprod(base_negative, null_negative)
null_any = crossprod(base_detect, null_detect)

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
    seq_along(threshold_sensitivity_lod),
    function(k) {

      this_lod = threshold_sensitivity_lod[[k]]
      message(
        "Pairwise threshold sensitivity: ",
        threshold_sensitivity_label[[k]],
        " (LOD ", signif(this_lod, 4), ")"
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
      out$threshold = threshold_sensitivity_label[[k]]
      out
    }
  )
)

observed_threshold_sensitivity = merge(
  threshold_pairwise,
  observed_for_join,
  by = c("query", "target"),
  all.x = TRUE,
  sort = FALSE
)
observed_threshold_sensitivity$expected_rediscovery =
  observed_threshold_sensitivity$rediscovery_same_direction
observed_threshold_sensitivity$difference =
  observed_threshold_sensitivity$observed_rediscovery -
  observed_threshold_sensitivity$expected_rediscovery

observed_threshold_sensitivity_summary =
  summarize_observed_expected(
    observed_threshold_sensitivity,
    grouping = "threshold"
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
    mean_rediscovery = mean(rediscovery_same_direction),
    median_rediscovery = stats::median(
      rediscovery_same_direction
    ),
    minimum_rediscovery = min(rediscovery_same_direction),
    maximum_rediscovery = max(rediscovery_same_direction),
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
    mean_rediscovery = mean(rediscovery_same_direction),
    median_rediscovery = stats::median(
      rediscovery_same_direction
    ),
    .groups = "drop"
  )

observed_summary = observed_pairwise %>%
  dplyr::group_by(comparison) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    total_query_qtl = sum(n_query_eligible),
    total_matches = sum(n_matches),
    mean_rediscovery = mean(
      observed_rediscovery,
      na.rm = TRUE
    ),
    median_rediscovery = stats::median(
      observed_rediscovery,
      na.rm = TRUE
    ),
    pooled_rediscovery =
      total_matches /
      total_query_qtl,
    .groups = "drop"
  )

simulation_calibration = data.frame(
  item = c(
    "analysis_universe",
    "covariate_handling",
    "effect_reference",
    "effect_sampling",
    "effect_scale_sensitivity",
    "threshold_reference",
    "threshold_sensitivity",
    "observed_pairwise_definition",
    "observed_pairwise_scope"
  ),
  value = c(
    "10 individual-segregant remapping datasets; AlbertAlt_2018 excluded as duplicate coordinate representation",
    "effects are defined after phenotype adjustment and standardization; covariates are not separately simulated",
    paste0(
      reference_study,
      " individually remapped QTL effects"
    ),
    "empirical resampling of absolute observed r values with replacement; no observed locus identity is reused as simulated truth",
    paste(effect_scales, collapse = ","),
    paste0(
      "median historical 5% FDR LOD calibration = ",
      signif(lod_reference, 6)
    ),
    paste(
      paste0(
        threshold_sensitivity_label,
        "=",
        signif(threshold_sensitivity_lod, 6)
      ),
      collapse = ";"
    ),
    "Frank Figure 3 overlapsWithStrand logic: CI_chr plus gene_ID defines the locus, both confidence intervals are padded 10 kb, and effect sign defines strand",
    "individual-only Figure 3 calibration uses Frank qtls_list_260204, target-gene measurability, and pair-specific LEU2/MAT/URA3/AMN1/KRE33 exclusions"
  ),
  stringsAsFactors = FALSE
)

########################
# WRITE TABLES         #
########################

write_table = function(x, filename) {
  write.table(
    x,
    file = paste0(power_asys, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

write_table(study_design, "study_design.tsv")
write_table(
  threshold_summary,
  "empirical_threshold_summary.tsv"
)
write_table(
  data.frame(lod = empirical_thresholds),
  "empirical_thresholds_all.tsv"
)
write_table(
  effect_fit,
  "effect_distribution_fit.tsv"
)
write_table(
  fit_quantiles,
  "effect_distribution_quantiles.tsv"
)
write_table(
  simulation_calibration,
  "simulation_calibration.tsv"
)
write_table(
  observed_pairwise,
  "observed_individual_rediscovery.tsv"
)
write_table(
  observed_summary,
  "observed_individual_rediscovery_summary.tsv"
)

write_table(
  observed_glmm_summary,
  "observed_individual_rediscovery_glmm.tsv"
)
write_table(
  power_curves,
  "study_power_curves.tsv"
)
write_table(
  threshold_sensitivity,
  "power_threshold_sensitivity.tsv"
)
write_table(
  validation,
  "analytic_power_validation.tsv"
)
write_table(
  shared_rediscovery,
  "shared_effect_rediscovery.tsv"
)
write_table(
  shared_summary,
  "shared_effect_rediscovery_summary.tsv"
)
write_table(
  observed_vs_expected,
  "observed_vs_expected.tsv"
)
write_table(
  observed_vs_expected_summary,
  "observed_vs_expected_summary.tsv"
)
write_table(
  observed_vs_expected_by_comparison,
  "observed_vs_expected_by_comparison.tsv"
)
write_table(
  heterogeneity_rediscovery,
  "heterogeneity_sensitivity.tsv"
)
write_table(
  heterogeneity_summary,
  "heterogeneity_sensitivity_summary.tsv"
)
write_table(
  threshold_pairwise,
  "pairwise_threshold_sensitivity.tsv"
)
write_table(
  observed_threshold_sensitivity,
  "observed_threshold_sensitivity.tsv"
)
write_table(
  observed_threshold_sensitivity_summary,
  "observed_threshold_sensitivity_summary.tsv"
)

review_power = list(
  settings = list(
    seed = seed,
    n_synthetic = n_synthetic,
    n_validation = n_validation,
    reference_study = reference_study,
    reference_lod = lod_reference,
    effect_sampling = "empirical_resampling_absolute_r",
    effect_scales = effect_scales,
    heterogeneity_fractions = heterogeneity_fractions
  ),
  study_design = study_design,
  empirical_thresholds = empirical_thresholds,
  threshold_summary = threshold_summary,
  effect_fit = effect_fit,
  effect_fit_quantiles = fit_quantiles,
  simulation_calibration = simulation_calibration,
  observed_pairwise = observed_pairwise,
  observed_summary = observed_summary,
  observed_glmm_summary = observed_glmm_summary,
  power_curves = power_curves,
  validation = validation,
  shared_rediscovery = shared_rediscovery,
  shared_summary = shared_summary,
  observed_vs_expected = observed_vs_expected,
  observed_vs_expected_summary = observed_vs_expected_summary,
  observed_vs_expected_by_comparison =
    observed_vs_expected_by_comparison,
  heterogeneity_rediscovery = heterogeneity_rediscovery,
  heterogeneity_summary = heterogeneity_summary,
  threshold_pairwise = threshold_pairwise,
  observed_threshold_sensitivity =
    observed_threshold_sensitivity,
  observed_threshold_sensitivity_summary =
    observed_threshold_sensitivity_summary
)

saveRDS(
  review_power,
  paste0(power_imed, "review_power.rds")
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
      seq_len(min(length(synthetic_abs_r), 100000L))
    ],
    source = "Empirically resampled effects",
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
  ggplot2::geom_density(linewidth = 0.8) +
  ggplot2::coord_cartesian(
    xlim = c(
      0,
      stats::quantile(reference_abs_r, 0.995)
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
  ggplot2::facet_wrap(~ study) +
  ggplot2::labs(
    x = "True standardized QTL effect (|r|)",
    y = "Probability of QTL detection",
    linetype = "Gene product"
  ) +
  ggplot2::theme_bw()

shared_heatmap_data = shared_rediscovery %>%
  dplyr::filter(effect_scale == 1)
shared_heatmap_data$query = factor(
  shared_heatmap_data$query,
  levels = rev(study_design$study)
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
    limits = c(0, 1),
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

observed_heatmap_data = observed_pairwise
observed_heatmap_data$query = factor(
  observed_heatmap_data$query,
  levels = rev(study_design$study)
)
observed_heatmap_data$target = factor(
  observed_heatmap_data$target,
  levels = study_design$study
)

observed_heatmap = ggplot2::ggplot(
  observed_heatmap_data,
  ggplot2::aes(
    x = target,
    y = query,
    fill = observed_rediscovery
  )
) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%d/%d",
        n_matches,
        n_query_eligible
      )
    ),
    size = 2.1
  ) +
  ggplot2::scale_fill_viridis_c(
    limits = c(0, 1),
    name = "Observed\nrediscovery"
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

observed_expected_plot = ggplot2::ggplot(
  observed_vs_expected,
  ggplot2::aes(
    x = expected_rediscovery,
    y = observed_rediscovery,
    shape = comparison
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  ggplot2::geom_point(
    alpha = 0.75,
    size = 2
  ) +
  ggplot2::facet_wrap(
    ~ effect_scale,
    labeller = ggplot2::label_both
  ) +
  ggplot2::coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  ggplot2::labs(
    x = "Expected rediscovery under shared effects",
    y = "Observed individual-study rediscovery",
    shape = "Comparison"
  ) +
  ggplot2::theme_bw()

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
  ggplot2::facet_wrap(~ heterogeneity_mode) +
  ggplot2::labs(
    x = "Fraction of true QTLs altered in target dataset",
    y = "Mean expected query-to-target rediscovery",
    linetype = "Comparison"
  ) +
  ggplot2::theme_bw()

save_plot = function(plot, stem, width, height) {
  ggplot2::ggsave(
    filename = paste0(power_asys, stem, ".pdf"),
    plot = plot,
    width = width,
    height = height
  )
  ggplot2::ggsave(
    filename = paste0(power_asys, stem, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

save_plot(
  effect_plot,
  "effect_distribution_fit",
  6,
  4
)
save_plot(
  power_plot,
  "study_power_curves",
  10,
  8
)
save_plot(
  shared_heatmap,
  "shared_effect_rediscovery_heatmap",
  10,
  8
)
save_plot(
  observed_heatmap,
  "observed_individual_rediscovery_heatmap",
  10,
  8
)
save_plot(
  observed_expected_plot,
  "observed_vs_expected",
  10,
  8
)
save_plot(
  heterogeneity_plot,
  "heterogeneity_sensitivity",
  10,
  5
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
  "Reference effect distribution: empirical resampling of ",
  reference_study,
  " (",
  format(length(reference_abs_r), big.mark = ","),
  " mapped QTLs)"
)
message(
  "Observed individual-study pairs: ",
  nrow(observed_pairwise)
)
message(
  "Reference LOD: ",
  signif(lod_reference, 4)
)
message(
  "Outputs: ",
  power_asys
)


########################################
########################################
# FINAL PAIRWISE REVIEW POLISH         #
########################################
########################################

canonicalize_review_seg = function(x) {

  x = trimws(
    as.character(x)
  )

  x = gsub(
    "_",
    ".",
    x,
    fixed = TRUE
  )

  gsub(
    "\\s+",
    "",
    x
  )
}


########################################
# SHARED-SEGREGANT PANEL AUDIT         #
########################################

build_panel_overlap = function(
  study_names,
  studies
) {

  pair_grid = expand.grid(
    query = study_names,
    target = study_names,
    stringsAsFactors = FALSE
  )

  pair_grid = pair_grid[
    pair_grid$query !=
    pair_grid$target,
    ,
    drop = FALSE
  ]

  do.call(
    rbind,
    lapply(
      seq_len(
        nrow(
          pair_grid
        )
      ),
      function(i) {

        query =
          pair_grid$query[[i]]

        target =
          pair_grid$target[[i]]

        query_ids =
          unique(
            canonicalize_review_seg(
              rownames(
                studies[[query]]$pheno
              )
            )
          )

        target_ids =
          unique(
            canonicalize_review_seg(
              rownames(
                studies[[target]]$pheno
              )
            )
          )

        shared =
          intersect(
            query_ids,
            target_ids
          )

        n_query =
          length(
            query_ids
          )

        n_target =
          length(
            target_ids
          )

        n_shared =
          length(
            shared
          )

        data.frame(
          query = query,
          target = target,
          n_query =
            n_query,
          n_target =
            n_target,
          n_shared =
            n_shared,
          fraction_query_shared =
            n_shared /
            n_query,
          fraction_target_shared =
            n_shared /
            n_target,
          fraction_smaller_panel_shared =
            n_shared /
            min(
              n_query,
              n_target
            ),
          jaccard =
            n_shared /
            length(
              union(
                query_ids,
                target_ids
              )
            ),
          major_shared_panel =
            n_shared >= 10L &
            n_shared /
            min(
              n_query,
              n_target
            ) >= 0.50,
          stringsAsFactors = FALSE
        )
      }
    )
  )
}


build_reuse_summary = function(
  study_names,
  studies,
  metadata
) {

  dat =
    do.call(
      rbind,
      lapply(
        study_names,
        function(study) {

          idx =
            match(
              study,
              metadata$study
            )

          data.frame(
            segregant =
              canonicalize_review_seg(
                rownames(
                  studies[[study]]$pheno
                )
              ),
            study = study,
            analyte =
              metadata$analyte[[idx]],
            stringsAsFactors = FALSE
          )
        }
      )
    )

  reuse =
    dat %>%
    dplyr::group_by(
      segregant
    ) %>%
    dplyr::summarise(
      n_studies =
        dplyr::n_distinct(
          study
        ),
      n_analytes =
        dplyr::n_distinct(
          analyte
        ),
      studies =
        paste(
          sort(
            unique(
              study
            )
          ),
          collapse = ";"
        ),
      .groups = "drop"
    )

  summary =
    data.frame(
      total_study_rows =
        nrow(dat),
      unique_segregants =
        dplyr::n_distinct(
          dat$segregant
        ),
      reused_segregants =
        sum(
          reuse$n_studies > 1L
        ),
      fraction_unique_reused =
        mean(
          reuse$n_studies > 1L
        ),
      maximum_studies_per_segregant =
        max(
          reuse$n_studies
        ),
      stringsAsFactors = FALSE
    )

  list(
    reuse = reuse,
    summary = summary
  )
}


individual_panel_overlap =
  build_panel_overlap(
    study_names =
      study_design$study,
    studies =
      studies
  )

joint_study_names =
  intersect(
    scaMetadata$study,
    names(
      studies
    )
  )

joint_panel_overlap =
  build_panel_overlap(
    study_names =
      joint_study_names,
    studies =
      studies
  )

joint_reuse =
  build_reuse_summary(
    study_names =
      joint_study_names,
    studies =
      studies,
    metadata =
      scaMetadata
  )


########################################
# FIGURE 3 SHARED-PANEL SENSITIVITY    #
########################################

observed_panel =
  dplyr::left_join(
    observed_pairwise,
    individual_panel_overlap,
    by = c(
      "query",
      "target"
    )
  )

panel_scenarios =
  list(
    all_pairs =
      rep(
        TRUE,
        nrow(
          observed_panel
        )
      ),
    exclude_major_shared_panels =
      !observed_panel$major_shared_panel
  )

panel_sensitivity_raw =
  do.call(
    rbind,
    lapply(
      names(
        panel_scenarios
      ),
      function(scenario) {

        dat =
          observed_panel[
            panel_scenarios[[scenario]],
            ,
            drop = FALSE
          ]

        dat %>%
          dplyr::group_by(
            comparison
          ) %>%
          dplyr::summarise(
            scenario =
              scenario,
            n_pairs =
              dplyr::n(),
            total_query_qtl =
              sum(
                n_query_eligible
              ),
            total_matches =
              sum(
                n_matches
              ),
            mean_rediscovery =
              mean(
                observed_rediscovery
              ),
            pooled_rediscovery =
              total_matches /
              total_query_qtl,
            .groups = "drop"
          )
      }
    )
  )

panel_sensitivity_glmm =
  do.call(
    rbind,
    lapply(
      names(
        panel_scenarios
      ),
      function(scenario) {

        dat =
          observed_panel[
            panel_scenarios[[scenario]],
            ,
            drop = FALSE
          ]

        fit =
          tryCatch(
            fit_observed_overlap_glmm(
              dat
            ),
            error = function(e) {
              NULL
            }
          )

        if (is.null(fit)) {

          return(
            data.frame(
              comparison =
                NA_character_,
              estimate_logit =
                NA_real_,
              se_logit =
                NA_real_,
              adjusted_agreement =
                NA_real_,
              adjusted_lower_95 =
                NA_real_,
              adjusted_upper_95 =
                NA_real_,
              singular_fit =
                NA,
              n_pairs =
                nrow(dat),
              scenario =
                scenario,
              stringsAsFactors = FALSE
            )
          )
        }

        fit$scenario =
          scenario

        fit
      }
    )
  )


########################################
# ANALYTIC FIGURE 3 CHANCE NULL        #
########################################

make_chr_lengths = function(
  chr_features
) {

  chr =
    normalize_yeast_chr(
      chr_features$Chromosome
    )

  left =
    suppressWarnings(
      as.numeric(
        chr_features$Start_Coordinate
      )
    )

  right =
    suppressWarnings(
      as.numeric(
        chr_features$Stop_Coordinate
      )
    )

  endpoint =
    pmax(
      left,
      right,
      na.rm = TRUE
    )

  keep =
    chr %in%
    as.character(1:16) &
    is.finite(endpoint)

  lengths =
    tapply(
      endpoint[keep],
      chr[keep],
      max,
      na.rm = TRUE
    )

  lengths =
    lengths[
      as.character(1:16)
    ]

  if (
    length(lengths) != 16L ||
    any(!is.finite(lengths))
  ) {
    stop(
      "Could not determine all chromosome lengths."
    )
  }

  lengths
}


prepare_frank_query_for_null = function(
  query,
  target,
  pad = 10000L
) {

  target_genes =
    colnames(
      studies[[target]]$pheno
    )

  qtls =
    frank_individual_qtls[[
      query
    ]]

  qtls =
    qtls[
      qtls$gene_ID %in%
      target_genes,
      ,
      drop = FALSE
    ]

  tech =
    frank_technical_markers()

  exclude =
    setdiff(
      tech[[query]],
      tech[[target]]
    )

  if (
    nrow(qtls) > 0L &&
    length(exclude) > 0L
  ) {

    left =
      pmin(
        as.numeric(
          qtls$CI_min_base
        ),
        as.numeric(
          qtls$CI_max_base
        )
      )

    right =
      pmax(
        as.numeric(
          qtls$CI_min_base
        ),
        as.numeric(
          qtls$CI_max_base
        )
      )

    gr =
      GenomicRanges::GRanges(
        seqnames =
          normalize_yeast_chr(
            qtls$CI_chr
          ),
        ranges =
          IRanges::IRanges(
            start =
              pmax(
                1,
                left -
                as.integer(pad)
              ),
            end =
              right +
              as.integer(pad)
          )
      )

    hit =
      suppressWarnings(
        GenomicRanges::countOverlaps(
          gr,
          technical_marker_ranges[
            exclude
          ],
          ignore.strand = TRUE
        )
      ) > 0L

    qtls =
      qtls[
        !hit,
        ,
        drop = FALSE
      ]
  }

  qtls
}


chance_overlap_pair = function(
  query,
  target,
  chromosome_lengths,
  pad = 10000L
) {

  query_qtls =
    prepare_frank_query_for_null(
      query = query,
      target = target,
      pad = pad
    )

  target_qtls =
    frank_individual_qtls[[
      target
    ]]

  if (
    nrow(query_qtls) == 0L
  ) {
    return(
      data.frame(
        query = query,
        target = target,
        n_query_eligible = 0L,
        expected_chance_matches = 0,
        chance_rediscovery = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  q_chr =
    normalize_yeast_chr(
      query_qtls$CI_chr
    )

  t_chr =
    normalize_yeast_chr(
      target_qtls$CI_chr
    )

  q_key =
    paste(
      q_chr,
      query_qtls$gene_ID,
      sign(
        query_qtls$cor_coef
      ),
      sep = "::"
    )

  t_key =
    paste(
      t_chr,
      target_qtls$gene_ID,
      sign(
        target_qtls$cor_coef
      ),
      sep = "::"
    )

  target_groups =
    split(
      seq_len(
        nrow(
          target_qtls
        )
      ),
      t_key
    )

  q_left =
    pmin(
      as.numeric(
        query_qtls$CI_min_base
      ),
      as.numeric(
        query_qtls$CI_max_base
      )
    )

  q_right =
    pmax(
      as.numeric(
        query_qtls$CI_min_base
      ),
      as.numeric(
        query_qtls$CI_max_base
      )
    )

  t_left =
    pmin(
      as.numeric(
        target_qtls$CI_min_base
      ),
      as.numeric(
        target_qtls$CI_max_base
      )
    )

  t_right =
    pmax(
      as.numeric(
        target_qtls$CI_min_base
      ),
      as.numeric(
        target_qtls$CI_max_base
      )
    )

  target_width =
    t_right -
    t_left

  probability =
    vapply(
      seq_len(
        nrow(
          query_qtls
        )
      ),
      function(i) {

        idx =
          target_groups[[
            q_key[[i]]
          ]]

        if (
          is.null(idx) ||
          length(idx) == 0L
        ) {
          return(0)
        }

        chr_length =
          unname(
            chromosome_lengths[
              q_chr[[i]]
            ]
          )

        widths =
          target_width[idx]

        n_starts =
          floor(
            chr_length -
            widths +
            1
          )

        valid =
          is.finite(n_starts) &
          n_starts >= 1

        if (!any(valid)) {
          return(0)
        }

        widths =
          widths[valid]

        n_starts =
          n_starts[valid]

        # Both query and target are padded 10 kb.
        lower =
          q_left[[i]] -
          widths -
          2 * as.integer(pad)

        upper =
          q_right[[i]] +
          2 * as.integer(pad)

        allowed_low =
          pmax(
            1,
            ceiling(lower)
          )

        allowed_high =
          pmin(
            n_starts,
            floor(upper)
          )

        n_allowed =
          pmax(
            0,
            allowed_high -
            allowed_low +
            1
          )

        p_each =
          pmin(
            1,
            n_allowed /
            n_starts
          )

        1 -
          prod(
            1 -
            p_each
          )
      },
      numeric(1)
    )

  data.frame(
    query = query,
    target = target,
    n_query_eligible =
      nrow(
        query_qtls
      ),
    expected_chance_matches =
      sum(
        probability
      ),
    chance_rediscovery =
      mean(
        probability
      ),
    stringsAsFactors = FALSE
  )
}


chromosome_lengths =
  make_chr_lengths(
    chrFeatures
  )

message(
  "Calculating Figure 3 chance-overlap null..."
)

chance_overlap_null =
  do.call(
    rbind,
    lapply(
      seq_len(
        nrow(
          observed_pairwise
        )
      ),
      function(i) {

        obs =
          observed_pairwise[
            i,
            ,
            drop = FALSE
          ]

        null =
          chance_overlap_pair(
            query =
              obs$query[[1]],
            target =
              obs$target[[1]],
            chromosome_lengths =
              chromosome_lengths,
            pad = 10000L
          )

        if (
          null$n_query_eligible !=
          obs$n_query_eligible
        ) {
          stop(
            "Chance-null denominator mismatch for ",
            obs$query[[1]],
            " -> ",
            obs$target[[1]]
          )
        }

        null$n_matches =
          obs$n_matches

        null$observed_rediscovery =
          obs$observed_rediscovery

        null$comparison =
          obs$comparison

        null
      }
    )
  )

chance_null_definition =
  data.frame(
    item = c(
      "null_model",
      "fixed_query_information",
      "randomized_quantity",
      "preserved_target_information",
      "overlap_definition"
    ),
    value = c(
      "independent target-QTL position randomization within chromosome",
      "Figure 3-eligible query QTLs after target-gene and technical-marker exclusions",
      "target-QTL genomic position",
      "gene, chromosome, effect direction, confidence-interval width, and QTL count",
      "same gene, chromosome, effect direction, and overlapping confidence intervals after 10-kb padding of both intervals"
    ),
    stringsAsFactors = FALSE
  )

chance_overlap_summary =
  chance_overlap_null %>%
  dplyr::group_by(
    comparison
  ) %>%
  dplyr::summarise(
    n_pairs =
      dplyr::n(),
    total_query_qtl =
      sum(
        n_query_eligible
      ),
    observed_matches =
      sum(
        n_matches
      ),
    expected_chance_matches =
      sum(
        expected_chance_matches
      ),
    mean_observed =
      mean(
        observed_rediscovery
      ),
    mean_chance =
      mean(
        chance_rediscovery
      ),
    pooled_observed =
      observed_matches /
      total_query_qtl,
    pooled_chance =
      expected_chance_matches /
      total_query_qtl,
    .groups = "drop"
  )


########################################
# SAVE POLISH OUTPUTS                  #
########################################

write_table(
  individual_panel_overlap,
  "segregant_overlap_individual_studies.tsv"
)

write_table(
  joint_panel_overlap,
  "segregant_overlap_joint_sc_studies.tsv"
)

write_table(
  joint_reuse$reuse,
  "segregant_reuse_joint_sc.tsv"
)

write_table(
  joint_reuse$summary,
  "segregant_reuse_joint_sc_summary.tsv"
)

write_table(
  panel_sensitivity_raw,
  "pairwise_shared_panel_sensitivity_raw.tsv"
)

write_table(
  panel_sensitivity_glmm,
  "pairwise_shared_panel_sensitivity_glmm.tsv"
)

write_table(
  chance_overlap_null,
  "chance_overlap_null.tsv"
)

write_table(
  chance_overlap_summary,
  "chance_overlap_null_summary.tsv"
)

write_table(
  chance_null_definition,
  "chance_overlap_null_definition.tsv"
)

review_power$individual_panel_overlap =
  individual_panel_overlap

review_power$joint_panel_overlap =
  joint_panel_overlap

review_power$joint_reuse =
  joint_reuse

review_power$panel_sensitivity_raw =
  panel_sensitivity_raw

review_power$panel_sensitivity_glmm =
  panel_sensitivity_glmm

review_power$chance_overlap_null =
  chance_overlap_null

review_power$chance_overlap_summary =
  chance_overlap_summary


review_power$chance_null_definition =
  chance_null_definition

saveRDS(
  review_power,
  paste0(
    power_imed,
    "review_power.rds"
  )
)

message(
  "Final pairwise reviewer polish complete."
)
