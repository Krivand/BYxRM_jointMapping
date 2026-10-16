############################################
############################################
## GENETICS REVIEW POWER ANALYSIS         ##
## ATTENUATION OBSERVED-VS-POWER PLOTS    ##
############################################
############################################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
power_asys = paste0(review_asys, "power_effect_curve/")
allpair_asys = paste0(power_asys, "all_pairs/")

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

comparison = read_tsv(
  paste0(allpair_asys, "pairwise_observed_vs_power.tsv")
)
attenuation = read_tsv(
  paste0(allpair_asys, "pairwise_power_attenuation.tsv")
)
summary = read_tsv(
  paste0(allpair_asys, "power_attenuation_summary.tsv")
)

required_comparison = c(
  "fromData",
  "toData",
  "analyteCombo",
  "numberOfQTLs",
  "fractionOverlapsWithStrand"
)
required_attenuation = c(
  "fromData",
  "toData",
  "expectedRediscoveryLOD4_5_scale0_8",
  "expectedRediscoveryLOD4_5_scale0_6"
)

if (
  nrow(comparison) != 90L ||
  nrow(attenuation) != 90L ||
  !all(required_comparison %in% names(comparison)) ||
  !all(required_attenuation %in% names(attenuation)) ||
  anyDuplicated(comparison[, c("fromData", "toData")]) ||
  anyDuplicated(attenuation[, c("fromData", "toData")])
) {
  stop("Attenuation plotting inputs failed structural checks.")
}

comparison_key = paste(comparison$fromData, comparison$toData, sep = "\r")
attenuation_key = paste(attenuation$fromData, attenuation$toData, sep = "\r")
idx = match(comparison_key, attenuation_key)
if (any(is.na(idx))) stop("Could not align attenuation rows to comparison rows.")
attenuation = attenuation[idx, , drop = FALSE]

if (
  nrow(summary) != 3L ||
  !all(c("effect_scale", "mean_observed", "mean_expected", "pearson_r", "spearman_rho") %in% names(summary))
) {
  stop("Attenuation summary failed structural checks.")
}

make_attenuation_plot = function(effect_scale, prediction_column, scale_label) {

  summary_row = summary[
    abs(summary$effect_scale - effect_scale) < 1e-12,
    ,
    drop = FALSE
  ]
  if (nrow(summary_row) != 1L) {
    stop("Could not identify summary row for effect scale ", effect_scale)
  }

  plot_data = comparison
  plot_data$observed_pct = 100 * plot_data$fractionOverlapsWithStrand
  plot_data$expected_pct = 100 * attenuation[[prediction_column]]

  if (
    any(!is.finite(plot_data$observed_pct)) ||
    any(!is.finite(plot_data$expected_pct)) ||
    any(plot_data$observed_pct < 0 | plot_data$observed_pct > 100) ||
    any(plot_data$expected_pct < 0 | plot_data$expected_pct > 100)
  ) {
    stop("Plotting values are not finite percentages for scale ", effect_scale)
  }

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
    ggplot2::scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Observed versus power-predicted QTL rediscovery",
      subtitle = sprintf(
        paste0(
          "%s source-effect sensitivity; LOD 4.5; 90 directed pairs\n",
          "Mean predicted %.1f%% vs observed %.1f%%; Pearson r = %.2f, Spearman rho = %.2f"
        ),
        scale_label,
        100 * summary_row$mean_expected,
        100 * summary_row$mean_observed,
        summary_row$pearson_r,
        summary_row$spearman_rho
      ),
      x = "Power-predicted same-direction rediscovery (%)",
      y = "Observed Figure 3 rediscovery (%)",
      color = "Comparison",
      size = "Query QTLs"
    ) +
    ggplot2::theme_bw()

  stem = paste0(
    allpair_asys,
    "all_pairs_observed_vs_power_lod4_5_scale",
    gsub("\\.", "_", format(effect_scale, nsmall = 1))
  )

  ggplot2::ggsave(
    filename = paste0(stem, ".pdf"),
    plot = p,
    width = 8,
    height = 7
  )
  ggplot2::ggsave(
    filename = paste0(stem, ".png"),
    plot = p,
    width = 8,
    height = 7,
    dpi = 300
  )

  c(paste0(stem, ".pdf"), paste0(stem, ".png"))
}

plot_files = c(
  make_attenuation_plot(
    effect_scale = 0.8,
    prediction_column = "expectedRediscoveryLOD4_5_scale0_8",
    scale_label = "80%"
  ),
  make_attenuation_plot(
    effect_scale = 0.6,
    prediction_column = "expectedRediscoveryLOD4_5_scale0_6",
    scale_label = "60%"
  )
)

missing_plots = plot_files[
  !file.exists(plot_files) |
    file.info(plot_files)$size <= 0
]
if (length(missing_plots) > 0L) {
  stop("Attenuation plot generation failed: ", paste(missing_plots, collapse = ", "))
}

cat("\nATTENUATION OBSERVED-VS-POWER PLOTS READY\n")
cat("LOD: 4.5\n")
cat("Effect scales plotted: 0.8, 0.6\n")
cat("Fixed axes: 0-100% on both dimensions\n")
cat("Files:\n")
cat(paste0("  ", plot_files, "\n"), sep = "")
