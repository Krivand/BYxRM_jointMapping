############################################
############################################
## GENETICS REVIEW POWER CALIBRATION      ##
## ALBERT QUERY VS FIGURE 3 REDISCOVERY   ##
############################################
############################################

#################
# PRELIMINARIES #
#################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
power_asys = paste0(review_asys, "power_effect_curve/")
compare_asys = paste0(power_asys, "figure3_calibration/")

ensure_directory(review_asys)
ensure_directory(power_asys)

# This directory contains only outputs from this script. Recreate it on every
# run so that no superseded comparison tables or plots survive.
if (dir.exists(compare_asys)) {
  unlink(compare_asys, recursive = TRUE, force = TRUE)
}
ensure_directory(compare_asys)

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

read_tsv = function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_tsv = function(x, filename) {
  write.table(
    x,
    file = paste0(compare_asys, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

parallel_lapply = function(X, FUN) {
  if (.Platform$OS.type == "unix" && power_cores > 1L && length(X) > 1L) {
    return(
      parallel::mclapply(
        X,
        FUN,
        mc.cores = min(power_cores, length(X)),
        mc.preschedule = TRUE
      )
    )
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

prepare_query_qtls = function(
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

  if (any(!is.finite(abs_r)) || any(abs_r < 0 | abs_r >= 1)) {
    stop("Effect magnitudes must satisfy 0 <= abs_r < 1.")
  }

  rcrit = sqrt(pmax(0, 1 - 10 ^ (-2 * lod / n)))
  df = n - 2
  tcrit = rcrit * sqrt(df / (1 - rcrit ^ 2))
  ncp = abs_r * sqrt((n - 1) / (1 - abs_r ^ 2))

  stats::pt(
    tcrit,
    df = df,
    ncp = ncp,
    lower.tail = FALSE
  )
}

########################
# LOAD CANONICAL INPUTS #
########################

crossstudy = read_tsv(
  paste0(power_asys, "albert_reference_crossstudy.tsv")
)

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

qtl_by_study = lapply(
  frank_name_map,
  function(name) frank_qtls_list[[name]]
)
names(qtl_by_study) = names(frank_name_map)

chrFeatures = load_single_rda(
  paste0(raw, "chrFeatures.rda"),
  object_name = "chrFeatures"
)
technical_marker_ranges = make_technical_marker_ranges(chrFeatures)

query = "Albert_2018"
targets = crossstudy$target

if (
  nrow(crossstudy) != 9L ||
  anyDuplicated(targets) ||
  query %in% targets
) {
  stop("Expected the compact nine-target Albert cross-study table.")
}

########################################
# RECONSTRUCT FIGURE 3 AND PREDICTIONS #
########################################

comparison_rows = parallel_lapply(
  targets,
  function(target) {

    cross_row = crossstudy[crossstudy$target == target, , drop = FALSE]
    if (nrow(cross_row) != 1L) {
      stop("Expected one cross-study row for target ", target)
    }

    # Figure 3 observed rediscovery uses the manuscript's original matching
    # rule: target-measurable query QTLs, pair-specific technical exclusions,
    # same gene/chromosome/direction, and 10-kb-padded confidence intervals.
    query_qtls = prepare_query_qtls(
      query = query,
      target = target,
      qtl_by_study = qtl_by_study,
      studies = studies,
      technical_marker_ranges = technical_marker_ranges,
      pad = interval_padding_bp
    )
    target_qtls = qtl_by_study[[target]]

    if (nrow(query_qtls) == 0L || nrow(target_qtls) == 0L) {
      n_match = 0L
    } else {
      query_gr = make_frank_qtl_granges(query_qtls)
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

    observed = if (nrow(query_qtls) > 0L) {
      n_match / nrow(query_qtls)
    } else {
      NA_real_
    }

    # Power predictions remain external to Figure 3 matching. They use all
    # Albert QTL effects for genes measured in the target and the target sample
    # size, but no target QTL calls or positional information.
    target_genes = colnames(studies[[target]]$pheno)
    reference = qtl_by_study[[query]]
    reference = reference[
      reference$gene_ID %in% target_genes,
      ,
      drop = FALSE
    ]
    reference_abs_r = abs(reference$cor_coef)

    predicted = vapply(
      lod_report,
      function(lod) {
        mean(
          same_direction_detection(
            abs_r = reference_abs_r,
            n = cross_row$target_n[[1]],
            lod = lod
          )
        )
      },
      numeric(1)
    )
    names(predicted) = c("lod2", "lod3_03", "lod4_5", "lod7")

    current_error = abs(
      predicted[["lod4_5"]] -
      cross_row$expected_same_direction_lod4_5[[1]]
    )
    if (current_error > 1e-12) {
      stop(
        "LOD 4.5 prediction does not reproduce the validated cross-study table for ",
        target, ": error = ", current_error
      )
    }

    data.frame(
      target = target,
      target_analyte = cross_row$target_analyte[[1]],
      target_medium = cross_row$target_medium[[1]],
      target_n = cross_row$target_n[[1]],
      n_figure3_query_eligible = nrow(query_qtls),
      n_figure3_matches = n_match,
      observed_figure3_rediscovery = observed,
      predicted_same_direction_lod4_5 = predicted[["lod4_5"]],
      difference_observed_minus_predicted_lod4_5 =
        observed - predicted[["lod4_5"]],
      difference_lod4_5_percentage_points =
        100 * (observed - predicted[["lod4_5"]]),
      predicted_same_direction_lod2 = predicted[["lod2"]],
      predicted_same_direction_lod3_03 = predicted[["lod3_03"]],
      predicted_same_direction_lod7 = predicted[["lod7"]],
      stringsAsFactors = FALSE
    )
  }
)

comparison = do.call(rbind, comparison_rows)
rownames(comparison) = NULL
comparison = comparison[match(targets, comparison$target), , drop = FALSE]

if (
  nrow(comparison) != 9L ||
  any(!is.finite(comparison$observed_figure3_rediscovery))
) {
  stop("Figure 3 calibration did not produce nine finite target comparisons.")
}

#################################
# FIT SUMMARY ACROSS THRESHOLDS #
#################################

prediction_columns = c(
  `2` = "predicted_same_direction_lod2",
  `3.03` = "predicted_same_direction_lod3_03",
  `4.5` = "predicted_same_direction_lod4_5",
  `7` = "predicted_same_direction_lod7"
)

fit_summary = do.call(
  rbind,
  lapply(
    names(prediction_columns),
    function(lod_name) {

      predicted = comparison[[prediction_columns[[lod_name]]]]
      observed = comparison$observed_figure3_rediscovery
      difference = observed - predicted

      data.frame(
        lod = as.numeric(lod_name),
        n_targets = length(observed),
        mean_observed = mean(observed),
        mean_predicted = mean(predicted),
        mean_observed_minus_predicted = mean(difference),
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

########################
# SUPPORTING SCATTER    #
########################

lod4 = fit_summary[fit_summary$lod == 4.5, , drop = FALSE]
plot_data = comparison
plot_data$observed_pct = 100 * plot_data$observed_figure3_rediscovery
plot_data$predicted_pct = 100 * plot_data$predicted_same_direction_lod4_5

p = ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = predicted_pct,
    y = observed_pct,
    label = target
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(
    check_overlap = TRUE,
    nudge_y = 0.25,
    size = 3
  ) +
  ggplot2::labs(
    title = "Albert 2018 query: observed versus power-predicted rediscovery",
    subtitle = sprintf(
      "LOD 4.5; Pearson r = %.2f, Spearman rho = %.2f",
      lod4$pearson_r,
      lod4$spearman_rho
    ),
    x = "Power-predicted same-direction detection (%)",
    y = "Observed Figure 3 rediscovery (%)"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  filename = paste0(
    compare_asys,
    "albert_figure3_observed_vs_power_lod4_5.pdf"
  ),
  plot = p,
  width = 7,
  height = 6
)

ggplot2::ggsave(
  filename = paste0(
    compare_asys,
    "albert_figure3_observed_vs_power_lod4_5.png"
  ),
  plot = p,
  width = 7,
  height = 6,
  dpi = 300
)

#################
# WRITE OUTPUTS #
#################

write_tsv(
  comparison,
  "albert_figure3_observed_vs_power.tsv"
)
write_tsv(
  fit_summary,
  "albert_figure3_power_fit_summary.tsv"
)

message("")
message("Albert/Figure 3 calibration complete.")
message("Output directory: ", compare_asys)
message("")
message("LOD 4.5 summary:")
print(lod4, row.names = FALSE)
message("")
message("Target-level comparison:")
print(
  comparison[, c(
    "target",
    "target_n",
    "observed_figure3_rediscovery",
    "predicted_same_direction_lod4_5",
    "difference_lod4_5_percentage_points"
  )],
  row.names = FALSE
)
