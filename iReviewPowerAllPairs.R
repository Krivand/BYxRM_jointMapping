############################################
############################################
## GENETICS REVIEW POWER ANALYSIS         ##
## ALL ORDERED DATASET PAIRS              ##
############################################
############################################

#################
# PRELIMINARIES #
#################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
power_asys = paste0(review_asys, "power_effect_curve/")
allpair_asys = paste0(power_asys, "all_pairs/")

ensure_directory(review_asys)
ensure_directory(power_asys)

if (dir.exists(allpair_asys)) {
  unlink(allpair_asys, recursive = TRUE, force = TRUE)
}
ensure_directory(allpair_asys)

lod_report = c(2, 3.03, 4.5, 7)
interval_padding_bp = 10000L
power_cores = suppressWarnings(
  as.integer(Sys.getenv("REVIEW_POWER_CORES", "1"))
)
if (is.na(power_cores) || power_cores < 1L) power_cores = 1L

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

write_tsv = function(x, filename) {
  write.table(
    x,
    file = paste0(allpair_asys, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
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

normalize_yeast_chr = function(x) {

  x = toupper(trimws(as.character(x)))
  x = sub("^CHR", "", x)

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

make_frank_qtl_granges = function(qtls, pad = interval_padding_bp) {

  required = c(
    "gene_ID",
    "CI_chr",
    "CI_min_base",
    "CI_max_base",
    "cor_coef"
  )
  missing = setdiff(required, colnames(qtls))
  if (length(missing) > 0L) {
    stop("Frank QTL table missing: ", paste(missing, collapse = ", "))
  }

  left = suppressWarnings(as.numeric(qtls$CI_min_base))
  right = suppressWarnings(as.numeric(qtls$CI_max_base))
  if (any(!is.finite(left) | !is.finite(right))) {
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

  GenomicRanges::GRanges(
    seqnames = seq_key,
    ranges = IRanges::IRanges(
      start = pmax(1, pmin(left, right) - as.integer(pad)),
      end = pmax(left, right) + as.integer(pad)
    ),
    strand = strand
  )
}

make_technical_marker_ranges = function(chr_features) {

  technical_orf = c(
    LEU2 = "YCL018W",
    MAT = "YCR039C",
    URA3 = "YEL021W",
    AMN1 = "YBR158W",
    KRE33 = "YNL132W"
  )

  idx = match(unname(technical_orf), chr_features$Feature_Name)
  if (any(is.na(idx))) {
    stop(
      "Technical-marker ORFs missing from chrFeatures: ",
      paste(names(technical_orf)[is.na(idx)], collapse = ", ")
    )
  }

  left = as.numeric(chr_features$Start_Coordinate[idx])
  right = as.numeric(chr_features$Stop_Coordinate[idx])

  gr = GenomicRanges::GRanges(
    seqnames = normalize_yeast_chr(chr_features$Chromosome[idx]),
    ranges = IRanges::IRanges(
      start = pmin(left, right),
      end = pmax(left, right)
    )
  )
  names(gr) = names(technical_orf)
  gr
}

frank_technical_markers = function() {

  standard = c("LEU2", "URA3", "AMN1", "MAT")

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

prepare_figure3_query_qtls = function(
  query,
  target,
  qtl_by_study,
  studies,
  technical_marker_ranges,
  pad = interval_padding_bp
) {

  target_genes = colnames(studies[[target]]$pheno)
  query_qtls = qtl_by_study[[query]]
  query_qtls = query_qtls[
    query_qtls$gene_ID %in% target_genes,
    ,
    drop = FALSE
  ]

  tech_list = frank_technical_markers()
  exclude_markers = setdiff(
    tech_list[[query]],
    tech_list[[target]]
  )

  if (nrow(query_qtls) > 0L && length(exclude_markers) > 0L) {

    left = as.numeric(query_qtls$CI_min_base)
    right = as.numeric(query_qtls$CI_max_base)

    query_chr_gr = GenomicRanges::GRanges(
      seqnames = normalize_yeast_chr(query_qtls$CI_chr),
      ranges = IRanges::IRanges(
        start = pmax(1, pmin(left, right) - as.integer(pad)),
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

    query_qtls = query_qtls[
      !technical_overlap,
      ,
      drop = FALSE
    ]
  }

  query_qtls
}

same_direction_detection = function(abs_r, n, lod) {

  if (any(!is.finite(abs_r)) || any(abs_r <= 0 | abs_r >= 1)) {
    stop("Effect magnitudes must satisfy 0 < abs_r < 1.")
  }
  if (length(n) != 1L || !is.finite(n) || n <= 2) {
    stop("n must be one finite scalar greater than 2.")
  }

  rcrit = sqrt(pmax(0, 1 - 10 ^ (-2 * lod / n)))
  df = n - 2
  tcrit = rcrit * sqrt(df / (1 - rcrit ^ 2))
  ncp = abs_r * sqrt((n - 1) / (1 - abs_r ^ 2))

  out = stats::pt(
    tcrit,
    df = df,
    ncp = ncp,
    lower.tail = FALSE
  )
  out[!is.finite(out)] = 0
  pmax(0, pmin(1, out))
}

########################
# LOAD CANONICAL INPUTS #
########################

studies = load_single_rda(
  paste0(imed, "studies.rda"),
  object_name = "studies"
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

analysis_studies = setdiff(
  classicMetadata$study,
  "AlbertAlt_2018"
)

if (!setequal(analysis_studies, names(frank_name_map))) {
  stop("The 10 remapped individual studies do not match frank_name_map.")
}

missing_studies = setdiff(analysis_studies, names(studies))
if (length(missing_studies) > 0L) {
  stop(
    "Missing processed studies: ",
    paste(missing_studies, collapse = ", ")
  )
}

missing_qtls = setdiff(unname(frank_name_map), names(frank_qtls_list))
if (length(missing_qtls) > 0L) {
  stop(
    "Missing QTL tables: ",
    paste(missing_qtls, collapse = ", ")
  )
}

qtl_by_study = lapply(
  frank_name_map,
  function(name) frank_qtls_list[[name]]
)
names(qtl_by_study) = names(frank_name_map)

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

chrFeatures = load_single_rda(
  paste0(raw, "chrFeatures.rda"),
  object_name = "chrFeatures"
)
technical_marker_ranges = make_technical_marker_ranges(chrFeatures)

pair_grid = expand.grid(
  query = analysis_studies,
  target = analysis_studies,
  stringsAsFactors = FALSE
)
pair_grid = pair_grid[
  pair_grid$query != pair_grid$target,
  ,
  drop = FALSE
]
rownames(pair_grid) = NULL

if (nrow(pair_grid) != 90L) {
  stop("Expected 90 ordered off-diagonal study pairs.")
}

############################################
# POWER PROJECTION AND OBSERVED REDISCOVERY #
############################################

message("Projecting every query QTL set into every other study...")

pair_results = parallel_lapply(
  seq_len(nrow(pair_grid)),
  function(i) {

    query = pair_grid$query[[i]]
    target = pair_grid$target[[i]]

    query_i = match(query, study_design$study)
    target_i = match(target, study_design$study)
    target_n = study_design$n[[target_i]]
    target_genes = colnames(studies[[target]]$pheno)

    # Power projection: use the query study's mapped QTL effect estimates for
    # genes measured in the target and ask how often effects of those same
    # magnitudes would cross a common LOD threshold at the target sample size.
    # No target QTL calls, positions, or overlap information enter this step.
    power_qtls = qtl_by_study[[query]]
    power_qtls = power_qtls[
      power_qtls$gene_ID %in% target_genes,
      ,
      drop = FALSE
    ]

    if (nrow(power_qtls) == 0L) {
      stop("No query QTLs remain for power projection: ", query, " -> ", target)
    }

    invalid_effect =
      !is.finite(power_qtls$cor_coef) |
      power_qtls$cor_coef == 0 |
      abs(power_qtls$cor_coef) >= 1
    if (any(invalid_effect)) {
      stop(
        "Invalid query effect estimates for ", query, " -> ", target,
        ": ", sum(invalid_effect)
      )
    }

    abs_r = abs(power_qtls$cor_coef)
    predicted = vapply(
      lod_report,
      function(lod) {
        mean(
          same_direction_detection(
            abs_r = abs_r,
            n = target_n,
            lod = lod
          )
        )
      },
      numeric(1)
    )
    names(predicted) = c("lod2", "lod3_03", "lod4_5", "lod7")

    # Observed Figure 3 rediscovery is reconstructed separately using the
    # manuscript matching rule. It is used only to make the diagnostic scatter;
    # it does not enter the power prediction above.
    figure3_query = prepare_figure3_query_qtls(
      query = query,
      target = target,
      qtl_by_study = qtl_by_study,
      studies = studies,
      technical_marker_ranges = technical_marker_ranges,
      pad = interval_padding_bp
    )
    target_qtls = qtl_by_study[[target]]

    if (nrow(figure3_query) == 0L || nrow(target_qtls) == 0L) {
      n_match = 0L
    } else {
      query_gr = make_frank_qtl_granges(figure3_query)
      target_gr = make_frank_qtl_granges(target_qtls)
      overlap_count = suppressWarnings(
        GenomicRanges::countOverlaps(
          query_gr,
          target_gr,
          ignore.strand = FALSE
        )
      )
      n_match = sum(overlap_count > 0L)
    }

    observed = if (nrow(figure3_query) > 0L) {
      n_match / nrow(figure3_query)
    } else {
      NA_real_
    }

    from_name = unname(frank_name_map[[query]])
    to_name = unname(frank_name_map[[target]])
    from_analyte = study_design$analyte[[query_i]]
    to_analyte = study_design$analyte[[target_i]]
    from_medium = study_design$medium[[query_i]]
    to_medium = study_design$medium[[target_i]]

    data.frame(
      query = query,
      target = target,
      fromData = from_name,
      toData = to_name,
      fromData_analyte = from_analyte,
      toData_analyte = to_analyte,
      fromData_medium = from_medium,
      toData_medium = to_medium,
      analyteCombo = paste0(from_analyte, "_", to_analyte),
      mediumSame = from_medium == to_medium,
      fromData_n = study_design$n[[query_i]],
      toData_n = target_n,
      nSourceQTLsUsedForPower = nrow(power_qtls),
      meanAbsSourceEffect = mean(abs_r),
      medianAbsSourceEffect = stats::median(abs_r),
      expectedRediscoveryLOD2 = predicted[["lod2"]],
      expectedRediscoveryLOD3_03 = predicted[["lod3_03"]],
      expectedRediscoveryLOD4_5 = predicted[["lod4_5"]],
      expectedRediscoveryLOD7 = predicted[["lod7"]],
      numberOfQTLs = nrow(figure3_query),
      overlapsWithStrand = n_match,
      fractionOverlapsWithStrand = observed,
      observedMinusExpectedLOD4_5 = observed - predicted[["lod4_5"]],
      stringsAsFactors = FALSE
    )
  }
)

pairwise_comparison = do.call(rbind, pair_results)
rownames(pairwise_comparison) = NULL

if (
  nrow(pairwise_comparison) != 90L ||
  anyDuplicated(pairwise_comparison[, c("fromData", "toData")]) ||
  any(pairwise_comparison$fromData == pairwise_comparison$toData)
) {
  stop("All-pair power table does not contain 90 unique directed pairs.")
}

prediction_columns = c(
  "expectedRediscoveryLOD2",
  "expectedRediscoveryLOD3_03",
  "expectedRediscoveryLOD4_5",
  "expectedRediscoveryLOD7"
)
if (
  any(!is.finite(as.matrix(pairwise_comparison[, prediction_columns]))) ||
  any(as.matrix(pairwise_comparison[, prediction_columns]) < 0) ||
  any(as.matrix(pairwise_comparison[, prediction_columns]) > 1) ||
  any(!is.finite(pairwise_comparison$fractionOverlapsWithStrand))
) {
  stop("All-pair prediction or observed fractions are invalid.")
}

#################################
# MERGE-READY TABLE FOR FRANK   #
#################################

power_join_columns = c(
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
  prediction_columns
)

pairwise_power_for_frank = pairwise_comparison[
  ,
  power_join_columns,
  drop = FALSE
]

#################################
# SUMMARIZE OBSERVED VS EXPECTED #
#################################

fit_summary = do.call(
  rbind,
  lapply(
    prediction_columns,
    function(column) {

      observed = pairwise_comparison$fractionOverlapsWithStrand
      predicted = pairwise_comparison[[column]]
      difference = observed - predicted

      lod = switch(
        column,
        expectedRediscoveryLOD2 = 2,
        expectedRediscoveryLOD3_03 = 3.03,
        expectedRediscoveryLOD4_5 = 4.5,
        expectedRediscoveryLOD7 = 7
      )

      data.frame(
        lod = lod,
        n_pairs = length(observed),
        mean_observed = mean(observed),
        mean_expected = mean(predicted),
        mean_observed_minus_expected = mean(difference),
        mean_abs_difference = mean(abs(difference)),
        rmse = sqrt(mean(difference ^ 2)),
        pearson_r = stats::cor(observed, predicted, method = "pearson"),
        spearman_rho = stats::cor(observed, predicted, method = "spearman"),
        stringsAsFactors = FALSE
      )
    }
  )
)
rownames(fit_summary) = NULL

#################################
# ALBERT SUBSET COMPATIBILITY   #
#################################

legacy_albert = read.delim(
  paste0(power_asys, "albert_reference_crossstudy.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

albert_rows = pairwise_comparison[
  pairwise_comparison$query == "Albert_2018",
  ,
  drop = FALSE
]
legacy_order = match(albert_rows$target, legacy_albert$target)
if (any(is.na(legacy_order)) || nrow(albert_rows) != 9L) {
  stop("Could not align the nine Albert projections to legacy power output.")
}
legacy_albert = legacy_albert[legacy_order, , drop = FALSE]

albert_error = max(abs(c(
  albert_rows$expectedRediscoveryLOD2 -
    legacy_albert$expected_detection_lod2,
  albert_rows$expectedRediscoveryLOD3_03 -
    legacy_albert$expected_detection_lod3_03,
  albert_rows$expectedRediscoveryLOD4_5 -
    legacy_albert$expected_same_direction_lod4_5,
  albert_rows$expectedRediscoveryLOD7 -
    legacy_albert$expected_detection_lod7
)))

# The legacy LOD 2/3.03/7 columns are two-sided detection probabilities,
# whereas the all-pair table is deliberately same-direction throughout.
# Only the LOD 4.5 same-direction column should be numerically identical.
albert_same_direction_error = max(abs(
  albert_rows$expectedRediscoveryLOD4_5 -
    legacy_albert$expected_same_direction_lod4_5
))
if (albert_same_direction_error > 1e-12) {
  stop(
    "All-pair Albert LOD 4.5 projections do not reproduce the existing ",
    "Albert analysis: error = ", albert_same_direction_error
  )
}

######################
# SUPPORTING SCATTER #
######################

plot_data = pairwise_comparison
plot_data$observed_pct = 100 * plot_data$fractionOverlapsWithStrand
plot_data$expected_pct = 100 * plot_data$expectedRediscoveryLOD4_5

lod4 = fit_summary[fit_summary$lod == 4.5, , drop = FALSE]

p = ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = expected_pct,
    y = observed_pct,
    color = analyteCombo,
    size = numberOfQTLs
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  ggplot2::geom_point(alpha = 0.75) +
  ggplot2::scale_size_area(max_size = 5) +
  ggplot2::labs(
    title = "Observed versus power-predicted QTL rediscovery",
    subtitle = sprintf(
      "90 directed individual-study pairs; LOD 4.5; Pearson r = %.2f, Spearman rho = %.2f",
      lod4$pearson_r,
      lod4$spearman_rho
    ),
    x = "Power-predicted same-direction rediscovery (%)",
    y = "Observed Figure 3 rediscovery (%)",
    color = "Comparison",
    size = "Query QTLs"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  filename = paste0(allpair_asys, "all_pairs_observed_vs_power_lod4_5.pdf"),
  plot = p,
  width = 8,
  height = 7
)
ggplot2::ggsave(
  filename = paste0(allpair_asys, "all_pairs_observed_vs_power_lod4_5.png"),
  plot = p,
  width = 8,
  height = 7,
  dpi = 300
)

#################
# WRITE OUTPUTS #
#################

write_tsv(
  pairwise_power_for_frank,
  "pairwise_power_for_frank.tsv"
)
write_tsv(
  pairwise_comparison,
  "pairwise_observed_vs_power.tsv"
)
write_tsv(
  fit_summary,
  "pairwise_power_fit_summary.tsv"
)

message("")
message("All-pair power projection complete.")
message("Directed off-diagonal pairs: ", nrow(pairwise_power_for_frank))
message("Merge keys: fromData + toData")
message("Output directory: ", allpair_asys)
message("")
message("LOD 4.5 all-pair summary:")
print(lod4, row.names = FALSE)
