########################################
########################################
## GENETICS REVIEW POWER ANALYSIS     ##
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
  if (any(!is.finite(rho))) {
    stop("rho contains non-finite values.")
  }
  if (any(abs(rho) >= 1)) {
    stop("rho must satisfy abs(rho) < 1.")
  }

  rcrit = lod_to_rcrit(lod = lod, n = n)
  df = n - 2
  tcrit = rcrit * sqrt(df / (1 - rcrit ^ 2))

  # Analytic approximation to the repeated-sampling distribution of the
  # correlation-test statistic under a nonzero population correlation.
  # This approximation is checked below by explicit phenotype simulation.
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

probability_matrix = function(
  rho,
  study_design,
  lod,
  component
) {

  if (!(component %in% c(
    "p_positive",
    "p_negative",
    "power"
  ))) {
    stop("Unknown probability component: ", component)
  }

  out = vapply(
    seq_len(nrow(study_design)),
    function(i) {
      probs = correlation_call_prob(
        rho = rho,
        n = study_design$n[[i]],
        lod = lod
      )
      probs[[component]]
    },
    numeric(length(rho))
  )

  colnames(out) = study_design$study
  out
}

pair_table_from_shared_effects = function(
  numerator_same,
  denominator_query,
  study_design,
  lod,
  threshold_label = "reference"
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

  pairs$query_detection_weight_sum =
    denominator_query[query_index]
  pairs$joint_same_direction_weight_sum =
    numerator_same[cbind(query_index, target_index)]
  pairs$expected_rediscovery =
    pairs$joint_same_direction_weight_sum /
    pairs$query_detection_weight_sum

  pairs$lod = lod
  pairs$threshold = threshold_label

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

expected_pairwise_shared = function(
  rho,
  study_design,
  lod,
  threshold_label = "reference"
) {

  # Rows are empirical Albert_2018 QTL effect magnitudes. The same row j is
  # propagated through every study; there is no r_j by r_k cross-product of
  # different effect sizes. Matrix cross-products below only sum matched row-j
  # joint detection probabilities across the empirical effect distribution.
  p_positive = probability_matrix(
    rho = rho,
    study_design = study_design,
    lod = lod,
    component = "p_positive"
  )
  p_negative = probability_matrix(
    rho = rho,
    study_design = study_design,
    lod = lod,
    component = "p_negative"
  )
  p_detect = p_positive + p_negative

  numerator_same =
    crossprod(p_positive, p_positive) +
    crossprod(p_negative, p_negative)
  denominator_query = colSums(p_detect)

  pair_table_from_shared_effects(
    numerator_same = numerator_same,
    denominator_query = denominator_query,
    study_design = study_design,
    lod = lod,
    threshold_label = threshold_label
  )
}

weighted_quantile = function(x, w, probs) {

  keep = is.finite(x) & is.finite(w) & w > 0
  x = x[keep]
  w = w[keep]

  if (length(x) == 0L) {
    return(rep(NA_real_, length(probs)))
  }

  ord = order(x)
  x = x[ord]
  w = w[ord]
  cw = cumsum(w) / sum(w)

  vapply(
    probs,
    function(p) x[which(cw >= p)[1]],
    numeric(1)
  )
}

query_ascertainment_summary = function(
  rho,
  study_design,
  lod
) {

  p_detect = probability_matrix(
    rho = rho,
    study_design = study_design,
    lod = lod,
    component = "power"
  )

  do.call(
    rbind,
    lapply(
      seq_len(nrow(study_design)),
      function(i) {
        w = p_detect[, i]
        q = weighted_quantile(
          x = rho,
          w = w,
          probs = c(0.10, 0.50, 0.90)
        )
        data.frame(
          study = study_design$study[[i]],
          analyte = study_design$analyte[[i]],
          n = study_design$n[[i]],
          mean_reference_rho = mean(rho),
          mean_rho_given_detection = sum(rho * w) / sum(w),
          q10_rho_given_detection = q[[1]],
          median_rho_given_detection = q[[2]],
          q90_rho_given_detection = q[[3]],
          detection_weight_sum = sum(w),
          stringsAsFactors = FALSE
        )
      }
    )
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
  threshold_sensitivity_lod = c(3.5, 4.5, 6)
  threshold_sensitivity_label = c(
    "fallback_low",
    "fallback_reference",
    "fallback_high"
  )
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
}

lod_override = suppressWarnings(
  as.numeric(Sys.getenv("REVIEW_POWER_LOD", ""))
)
if (is.finite(lod_override)) {
  lod_reference = lod_override
  threshold_sensitivity_lod[[2]] = lod_override
  threshold_sensitivity_label[[2]] = "reference_override"
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
  reference_abs_r < 1
]

if (length(reference_abs_r) < 100L) {
  stop(
    paste0(
      "Too few Albert_2018 QTL effects available for ",
      "empirical effect-size calibration."
    )
  )
}

# Each Albert_2018 QTL contributes one equally weighted empirical support
# point. No bootstrap resampling, parametric fit, global effect scaling, or
# upper-tail truncation is used in the primary rediscovery calculation.
effect_reference = data.frame(
  reference_study = reference_study,
  n_reference_qtl = length(reference_abs_r),
  reference_distribution =
    "direct_empirical_absolute_r_among_detected_qtls",
  weighting = "one_equal_weight_per_mapped_qtl",
  stringsAsFactors = FALSE
)

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
effect_quantiles = data.frame(
  quantile = quantile_probs,
  observed_abs_r = as.numeric(
    stats::quantile(
      reference_abs_r,
      probs = quantile_probs,
      names = FALSE
    )
  ),
  stringsAsFactors = FALSE
)

############################
# OBSERVED FIGURE 3 OVERLAP #
############################

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

  idx = match(x, roman)
  use = !is.na(idx)
  x[use] = as.character(idx[use])
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

  missing = setdiff(required, colnames(qtls))
  if (length(missing) > 0L) {
    stop(
      "Frank QTL table missing: ",
      paste(missing, collapse = ", ")
    )
  }

  left = suppressWarnings(as.numeric(qtls$CI_min_base))
  right = suppressWarnings(as.numeric(qtls$CI_max_base))
  bad = !is.finite(left) | !is.finite(right)
  if (any(bad)) {
    stop("Non-finite Frank QTL confidence intervals encountered.")
  }

  seq_key = paste(
    normalize_yeast_chr(qtls$CI_chr),
    qtls$gene_ID,
    sep = "_"
  )

  strand = ifelse(
    qtls$cor_coef > 0,
    "+",
    ifelse(qtls$cor_coef < 0, "-", "*")
  )

  gr = GenomicRanges::GRanges(
    seqnames = seq_key,
    ranges = IRanges::IRanges(
      start = pmax(
        1,
        pmin(left, right) - as.integer(pad)
      ),
      end = pmax(left, right) + as.integer(pad)
    ),
    strand = strand
  )

  S4Vectors::mcols(gr)$gene = qtls$gene_ID
  S4Vectors::mcols(gr)$r = qtls$cor_coef
  gr
}

make_technical_marker_ranges = function(chr_features) {

  required = c(
    "Feature_Name",
    "Chromosome",
    "Start_Coordinate",
    "Stop_Coordinate"
  )

  missing = setdiff(required, colnames(chr_features))
  if (length(missing) > 0L) {
    stop(
      paste0(
        "chrFeatures missing columns required for Frank technical-marker ",
        "exclusions: ",
        paste(missing, collapse = ", ")
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

  idx = match(unname(tech_orf), chr_features$Feature_Name)
  if (any(is.na(idx))) {
    stop(
      "Technical-marker ORFs missing from chrFeatures: ",
      paste(names(tech_orf)[is.na(idx)], collapse = ", ")
    )
  }

  left = suppressWarnings(
    as.numeric(chr_features$Start_Coordinate[idx])
  )
  right = suppressWarnings(
    as.numeric(chr_features$Stop_Coordinate[idx])
  )

  gr = GenomicRanges::GRanges(
    seqnames = normalize_yeast_chr(
      chr_features$Chromosome[idx]
    ),
    ranges = IRanges::IRanges(
      start = pmin(left, right),
      end = pmax(left, right)
    )
  )

  names(gr) = names(tech_orf)
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
    Albert_2018 = c("MAT", "KRE33")
  )
}

build_observed_pairwise_frank = function(
  qtl_by_study,
  studies,
  study_design,
  technical_marker_ranges,
  pad = 10000L
) {

  tech_list = frank_technical_markers()

  pair_grid = expand.grid(
    query = study_design$study,
    target = study_design$study,
    stringsAsFactors = FALSE
  )
  pair_grid = pair_grid[
    pair_grid$query != pair_grid$target,
    ,
    drop = FALSE
  ]

  out = lapply(
    seq_len(nrow(pair_grid)),
    function(i) {

      query = pair_grid$query[[i]]
      target = pair_grid$target[[i]]

      query_genes = colnames(studies[[query]]$pheno)
      target_genes = colnames(studies[[target]]$pheno)
      query_qtls = qtl_by_study[[query]]
      target_qtls = qtl_by_study[[target]]

      # Original Figure 3 denominator: query QTLs must affect genes measured
      # in the target dataset.
      query_qtls = query_qtls[
        query_qtls$gene_ID %in% target_genes,
        ,
        drop = FALSE
      ]
      n_before_technical = nrow(query_qtls)

      # Remove technical loci present in the query cross but absent from the
      # target cross, matching the original Figure 3 comparison.
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
          as.numeric(query_qtls$CI_min_base)
        )
        right = suppressWarnings(
          as.numeric(query_qtls$CI_max_base)
        )

        query_chr_gr = GenomicRanges::GRanges(
          seqnames = normalize_yeast_chr(query_qtls$CI_chr),
          ranges = IRanges::IRanges(
            start = pmax(
              1,
              pmin(left, right) - as.integer(pad)
            ),
            end = pmax(left, right) + as.integer(pad)
          )
        )

        technical_overlap = suppressWarnings(
          GenomicRanges::countOverlaps(
            query_chr_gr,
            technical_marker_ranges[exclude_markers],
            ignore.strand = TRUE
          )
        ) > 0L

        n_technical_excluded = sum(technical_overlap)
        query_qtls = query_qtls[
          !technical_overlap,
          ,
          drop = FALSE
        ]
      }

      n_query = nrow(query_qtls)

      if (
        n_query == 0L ||
        nrow(target_qtls) == 0L
      ) {
        n_match = 0L
      } else {
        query_gr = make_frank_qtl_granges(
          query_qtls,
          pad = pad
        )
        target_gr = make_frank_qtl_granges(
          target_qtls,
          pad = pad
        )

        overlap_count = suppressWarnings(
          GenomicRanges::countOverlaps(
            query_gr,
            target_gr,
            ignore.strand = FALSE
          )
        )
        n_match = sum(overlap_count > 0L)
      }

      data.frame(
        query = query,
        target = target,
        n_common_genes = length(
          intersect(query_genes, target_genes)
        ),
        n_query_before_technical = n_before_technical,
        n_technical_excluded = n_technical_excluded,
        n_query_eligible = n_query,
        n_target_eligible = sum(
          target_qtls$gene_ID %in% query_genes
        ),
        n_matches = n_match,
        observed_rediscovery = if (n_query > 0L) {
          n_match / n_query
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  )

  observed = do.call(rbind, out)

  idx = match(observed$query, study_design$study)
  observed$query_analyte = study_design$analyte[idx]
  observed$query_medium = study_design$medium[idx]
  observed$query_n = study_design$n[idx]

  idx = match(observed$target, study_design$study)
  observed$target_analyte = study_design$analyte[idx]
  observed$target_medium = study_design$medium[idx]
  observed$target_n = study_design$n[idx]

  observed$comparison = paste0(
    observed$query_analyte,
    "_to_",
    observed$target_analyte
  )
  observed$same_medium =
    observed$query_medium == observed$target_medium
  observed$interval_padding_bp = as.integer(pad)
  observed$overlap_definition =
    "Frank_Figure3_overlapsWithStrand"

  observed
}

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

missing_frank = setdiff(
  unname(frank_name_map),
  names(frank_qtls_list)
)
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

chrFeatures = load_single_rda(
  paste0(raw, "chrFeatures.rda"),
  object_name = "chrFeatures"
)
technical_marker_ranges = make_technical_marker_ranges(
  chrFeatures
)

frank_qtl_count_check = data.frame(
  study = study_design$study,
  indivFullQTLs = study_design$n_qtl,
  frank_qtls_list = vapply(
    study_design$study,
    function(study) {
      nrow(frank_individual_qtls[[study]])
    },
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
      frank_qtl_count_check$study[
        !frank_qtl_count_check$match
      ],
      collapse = ", "
    )
  )
}

message(
  "Frank QTL object matches current remap QTL counts for all 10 datasets."
)
message(
  "Recomputing observed individual-study rediscovery using Frank Figure 3 logic..."
)

observed_pairwise = build_observed_pairwise_frank(
  qtl_by_study = frank_individual_qtls,
  studies = studies,
  study_design = study_design,
  technical_marker_ranges = technical_marker_ranges,
  pad = 10000L
)

########################
# GENERAL POWER CURVES #
########################

# This grid is only for showing P(detection | rho, n, L). It is not used as
# an effect-size distribution and is not averaged to obtain Figure 3 expected
# rediscovery.
rho_grid = seq(0, 0.95, by = 0.005)
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

##############################
# DIRECT PHENOTYPE VALIDATION #
##############################

validation_n = unique(
  c(
    min(study_design$n),
    study_design$n[
      which.min(
        abs(
          study_design$n - stats::median(study_design$n)
        )
      )
    ],
    max(study_design$n)
  )
)
validation_rho = unique(
  as.numeric(
    stats::quantile(
      reference_abs_r,
      probs = c(0.10, 0.50, 0.90),
      names = FALSE
    )
  )
)
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
        seed = validation_seed + i
      )
    }
  )
)
validation$power_difference =
  validation$simulated_power - validation$analytic_power
validation$direction_difference =
  validation$simulated_same_direction -
  validation$analytic_same_direction

##################################
# PERFECTLY SHARED REDISCOVERY   #
##################################

message(
  "Calculating perfectly shared rediscovery directly over ",
  format(length(reference_abs_r), big.mark = ","),
  " Albert_2018 QTL effects..."
)

shared_rediscovery = expected_pairwise_shared(
  rho = reference_abs_r,
  study_design = study_design,
  lod = lod_reference,
  threshold_label = "reference_median"
)

query_ascertainment = query_ascertainment_summary(
  rho = reference_abs_r,
  study_design = study_design,
  lod = lod_reference
)

##################################
# THRESHOLD SENSITIVITY          #
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
      expected_pairwise_shared(
        rho = reference_abs_r,
        study_design = study_design,
        lod = this_lod,
        threshold_label = threshold_sensitivity_label[[k]]
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
observed_vs_expected$difference =
  observed_vs_expected$observed_rediscovery -
  observed_vs_expected$expected_rediscovery

observed_threshold_sensitivity = merge(
  threshold_pairwise,
  observed_for_join,
  by = c("query", "target"),
  all.x = TRUE,
  sort = FALSE
)
observed_threshold_sensitivity$difference =
  observed_threshold_sensitivity$observed_rediscovery -
  observed_threshold_sensitivity$expected_rediscovery

observed_vs_expected_summary = observed_vs_expected %>%
  dplyr::filter(
    is.finite(observed_rediscovery),
    is.finite(expected_rediscovery)
  ) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    mean_observed = mean(observed_rediscovery),
    mean_expected = mean(expected_rediscovery),
    mean_difference = mean(difference)
  )

observed_vs_expected_by_comparison = observed_vs_expected %>%
  dplyr::filter(
    is.finite(observed_rediscovery),
    is.finite(expected_rediscovery)
  ) %>%
  dplyr::group_by(comparison) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    mean_observed = mean(observed_rediscovery),
    mean_expected = mean(expected_rediscovery),
    mean_difference = mean(difference),
    .groups = "drop"
  )

observed_threshold_sensitivity_summary =
  observed_threshold_sensitivity %>%
  dplyr::filter(
    is.finite(observed_rediscovery),
    is.finite(expected_rediscovery)
  ) %>%
  dplyr::group_by(threshold, lod) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    mean_observed = mean(observed_rediscovery),
    mean_expected = mean(expected_rediscovery),
    mean_difference = mean(difference),
    .groups = "drop"
  )

shared_summary = shared_rediscovery %>%
  dplyr::group_by(comparison) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    mean_rediscovery = mean(expected_rediscovery),
    median_rediscovery = stats::median(expected_rediscovery),
    minimum_rediscovery = min(expected_rediscovery),
    maximum_rediscovery = max(expected_rediscovery),
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
    pooled_rediscovery = total_matches / total_query_qtl,
    .groups = "drop"
  )

########################################
# ANALYTIC FIGURE 3 CHANCE NULL        #
########################################

make_chr_lengths = function(chr_features) {

  chr = normalize_yeast_chr(chr_features$Chromosome)
  left = suppressWarnings(
    as.numeric(chr_features$Start_Coordinate)
  )
  right = suppressWarnings(
    as.numeric(chr_features$Stop_Coordinate)
  )
  endpoint = pmax(left, right, na.rm = TRUE)

  keep =
    chr %in% as.character(1:16) &
    is.finite(endpoint)

  lengths = tapply(
    endpoint[keep],
    chr[keep],
    max,
    na.rm = TRUE
  )
  lengths = lengths[as.character(1:16)]

  if (
    length(lengths) != 16L ||
    any(!is.finite(lengths))
  ) {
    stop("Could not determine all chromosome lengths.")
  }

  lengths
}

prepare_frank_query_for_null = function(
  query,
  target,
  pad = 10000L
) {

  target_genes = colnames(studies[[target]]$pheno)
  qtls = frank_individual_qtls[[query]]
  qtls = qtls[
    qtls$gene_ID %in% target_genes,
    ,
    drop = FALSE
  ]

  tech = frank_technical_markers()
  exclude = setdiff(
    tech[[query]],
    tech[[target]]
  )

  if (
    nrow(qtls) > 0L &&
    length(exclude) > 0L
  ) {

    left = pmin(
      as.numeric(qtls$CI_min_base),
      as.numeric(qtls$CI_max_base)
    )
    right = pmax(
      as.numeric(qtls$CI_min_base),
      as.numeric(qtls$CI_max_base)
    )

    gr = GenomicRanges::GRanges(
      seqnames = normalize_yeast_chr(qtls$CI_chr),
      ranges = IRanges::IRanges(
        start = pmax(
          1,
          left - as.integer(pad)
        ),
        end = right + as.integer(pad)
      )
    )

    hit = suppressWarnings(
      GenomicRanges::countOverlaps(
        gr,
        technical_marker_ranges[exclude],
        ignore.strand = TRUE
      )
    ) > 0L

    qtls = qtls[
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

  query_qtls = prepare_frank_query_for_null(
    query = query,
    target = target,
    pad = pad
  )
  target_qtls = frank_individual_qtls[[target]]

  if (nrow(query_qtls) == 0L) {
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

  q_chr = normalize_yeast_chr(query_qtls$CI_chr)
  t_chr = normalize_yeast_chr(target_qtls$CI_chr)

  q_key = paste(
    q_chr,
    query_qtls$gene_ID,
    sign(query_qtls$cor_coef),
    sep = "::"
  )
  t_key = paste(
    t_chr,
    target_qtls$gene_ID,
    sign(target_qtls$cor_coef),
    sep = "::"
  )

  target_groups = split(
    seq_len(nrow(target_qtls)),
    t_key
  )

  q_left = pmin(
    as.numeric(query_qtls$CI_min_base),
    as.numeric(query_qtls$CI_max_base)
  )
  q_right = pmax(
    as.numeric(query_qtls$CI_min_base),
    as.numeric(query_qtls$CI_max_base)
  )
  t_left = pmin(
    as.numeric(target_qtls$CI_min_base),
    as.numeric(target_qtls$CI_max_base)
  )
  t_right = pmax(
    as.numeric(target_qtls$CI_min_base),
    as.numeric(target_qtls$CI_max_base)
  )
  target_width = t_right - t_left

  probability = vapply(
    seq_len(nrow(query_qtls)),
    function(i) {

      idx = target_groups[[q_key[[i]]]]
      if (is.null(idx) || length(idx) == 0L) {
        return(0)
      }

      chr_length = unname(
        chromosome_lengths[q_chr[[i]]]
      )
      widths = target_width[idx]
      n_starts = floor(chr_length - widths + 1)
      valid = is.finite(n_starts) & n_starts >= 1

      if (!any(valid)) {
        return(0)
      }

      widths = widths[valid]
      n_starts = n_starts[valid]

      # Both query and target confidence intervals receive 10-kb padding.
      lower =
        q_left[[i]] - widths - 2 * as.integer(pad)
      upper =
        q_right[[i]] + 2 * as.integer(pad)

      allowed_low = pmax(1, ceiling(lower))
      allowed_high = pmin(n_starts, floor(upper))
      n_allowed = pmax(
        0,
        allowed_high - allowed_low + 1
      )
      p_each = pmin(1, n_allowed / n_starts)

      1 - prod(1 - p_each)
    },
    numeric(1)
  )

  data.frame(
    query = query,
    target = target,
    n_query_eligible = nrow(query_qtls),
    expected_chance_matches = sum(probability),
    chance_rediscovery = mean(probability),
    stringsAsFactors = FALSE
  )
}

chromosome_lengths = make_chr_lengths(chrFeatures)

message("Calculating Figure 3 chance-overlap null...")

chance_overlap_null = do.call(
  rbind,
  lapply(
    seq_len(nrow(observed_pairwise)),
    function(i) {

      obs = observed_pairwise[
        i,
        ,
        drop = FALSE
      ]

      null = chance_overlap_pair(
        query = obs$query[[1]],
        target = obs$target[[1]],
        chromosome_lengths = chromosome_lengths,
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

      null$n_matches = obs$n_matches
      null$observed_rediscovery = obs$observed_rediscovery
      null$comparison = obs$comparison
      null
    }
  )
)

chance_overlap_summary = chance_overlap_null %>%
  dplyr::group_by(comparison) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    total_query_qtl = sum(n_query_eligible),
    observed_matches = sum(n_matches),
    expected_chance_matches = sum(expected_chance_matches),
    mean_observed = mean(observed_rediscovery),
    mean_chance = mean(chance_rediscovery),
    pooled_observed = observed_matches / total_query_qtl,
    pooled_chance = expected_chance_matches / total_query_qtl,
    .groups = "drop"
  )

chance_null_definition = data.frame(
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

########################
# CALIBRATION SUMMARY  #
########################

simulation_calibration = data.frame(
  item = c(
    "analysis_universe",
    "covariate_handling",
    "effect_reference",
    "effect_weighting",
    "effect_size_sweep",
    "power_curve_grid",
    "threshold_reference",
    "threshold_sensitivity",
    "shared_effect_definition",
    "shared_effect_aggregation",
    "localization_scope",
    "heterogeneity_scope",
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
    "each finite 0 < abs(r) < 1 Albert_2018 mapped QTL enters once with equal empirical weight; no resampling or parametric fit",
    "none; the former 0.20-1.00 global scaling sweep was removed",
    "rho 0-0.95 by 0.005 for visualization of P(detection | rho,n,L) only; grid values are not averaged as an effect distribution",
    paste0(
      "pooled median historical 5% FDR LOD calibration = ",
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
    "the same population correlation rho_j is supplied to query and target for every empirical effect j",
    "same-row joint detection probabilities are summed over empirical effects and divided by the summed query detection probabilities; no r_j by r_k effect-size crossing",
    "expected rediscovery models threshold crossing and sign, not same-locus peak-localization error; observed Figure 3 additionally requires overlapping confidence intervals",
    "no arbitrary attenuation/null/opposite-effect mixture is fit in the primary power analysis; biological/environmental heterogeneity remains a stated limitation",
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
  effect_reference,
  "effect_distribution_reference.tsv"
)
write_table(
  effect_quantiles,
  "effect_distribution_quantiles.tsv"
)
write_table(
  data.frame(abs_r = reference_abs_r),
  "effect_distribution_values.tsv"
)
write_table(
  simulation_calibration,
  "simulation_calibration.tsv"
)
write_table(
  frank_qtl_count_check,
  "frank_qtl_count_check.tsv"
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
  power_curves,
  "study_power_curves.tsv"
)
write_table(
  validation,
  "analytic_power_validation.tsv"
)
write_table(
  query_ascertainment,
  "query_ascertainment_effect_summary.tsv"
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

review_power = list(
  settings = list(
    validation_seed = validation_seed,
    n_validation = n_validation,
    reference_study = reference_study,
    reference_lod = lod_reference,
    effect_reference = "direct_empirical_absolute_r",
    effect_size_sweep = "none",
    threshold_sensitivity = stats::setNames(
      threshold_sensitivity_lod,
      threshold_sensitivity_label
    )
  ),
  study_design = study_design,
  empirical_thresholds = empirical_thresholds,
  threshold_summary = threshold_summary,
  effect_reference = effect_reference,
  effect_quantiles = effect_quantiles,
  simulation_calibration = simulation_calibration,
  observed_pairwise = observed_pairwise,
  observed_summary = observed_summary,
  power_curves = power_curves,
  validation = validation,
  query_ascertainment = query_ascertainment,
  shared_rediscovery = shared_rediscovery,
  shared_summary = shared_summary,
  observed_vs_expected = observed_vs_expected,
  observed_vs_expected_summary = observed_vs_expected_summary,
  observed_vs_expected_by_comparison =
    observed_vs_expected_by_comparison,
  threshold_pairwise = threshold_pairwise,
  observed_threshold_sensitivity =
    observed_threshold_sensitivity,
  observed_threshold_sensitivity_summary =
    observed_threshold_sensitivity_summary,
  chance_overlap_null = chance_overlap_null,
  chance_overlap_summary = chance_overlap_summary,
  chance_null_definition = chance_null_definition
)

saveRDS(
  review_power,
  paste0(power_imed, "review_power.rds")
)

########################
# DIAGNOSTIC FIGURES   #
########################

effect_plot = ggplot2::ggplot(
  data.frame(abs_r = reference_abs_r),
  ggplot2::aes(x = abs_r)
) +
  ggplot2::geom_density(linewidth = 0.8) +
  ggplot2::coord_cartesian(
    xlim = c(
      0,
      stats::quantile(reference_abs_r, 0.995)
    )
  ) +
  ggplot2::labs(
    x = "Absolute Albert 2018 QTL effect (|r|)",
    y = "Density"
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
    x = "True standardized QTL effect (|rho|)",
    y = "Probability of QTL detection",
    linetype = "Gene product"
  ) +
  ggplot2::theme_bw()

shared_heatmap_data = shared_rediscovery
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
    fill = expected_rediscovery
  )
) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf("%.2f", expected_rediscovery)
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
  ggplot2::coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  ggplot2::labs(
    x = "Expected rediscovery under identical shared effects",
    y = "Observed Figure 3 rediscovery",
    shape = "Comparison"
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
  "effect_distribution_reference",
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
  6,
  6
)

########################
# FINAL SUMMARY        #
########################

message("")
message("Reviewer power analysis complete.")
message(
  "Individual datasets analyzed: ",
  nrow(study_design)
)
message(
  "Reference effect distribution: direct empirical ",
  reference_study,
  " QTL effects (",
  format(length(reference_abs_r), big.mark = ","),
  " mapped QTLs)"
)
message(
  "Effect-size scaling sweep: none"
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
