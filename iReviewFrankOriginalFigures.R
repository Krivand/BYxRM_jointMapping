# Preserve the original two manuscript plot blocks; replace predictions only.
args = commandArgs(TRUE)
out = args[[1]]
old = "/PHShome/kv1026/p05/asys/review/power_effect_curve/all_pairs/"
read_tsv = function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
comparison = read_tsv(paste0(old, "pairwise_observed_vs_power.tsv"))
attenuation = read_tsv(paste0(old, "pairwise_power_attenuation.tsv"))
updated = read_tsv("/PHShome_actual/k/kv1026/tools/eris-jobcards/actions/card_1e080f28d4b1c883cb4f66157701e1a7/power/pairwise_actual_threshold_power.tsv")
key = function(d) paste(d$fromData, d$toData, sep = "\r")
stopifnot(nrow(comparison) == 90L, !anyDuplicated(key(comparison)))
ai = match(key(comparison), key(attenuation))
stopifnot(!anyNA(ai))
attenuation = attenuation[ai, , drop = FALSE]
for (effect_scale in c(1, 0.8)) {
  u = updated[updated$denominator == "historical_gene_overlap" & abs(updated$effect_scale - effect_scale) < 1e-12, , drop = FALSE]
  stopifnot(nrow(u) == 90L, !anyDuplicated(key(u)))
  idx = match(key(comparison), key(u))
  stopifnot(!anyNA(idx))
  u = u[idx, , drop = FALSE]
  fixed = if (effect_scale == 1) comparison$expectedRediscoveryLOD4_5 else attenuation$expectedRediscoveryLOD4_5_scale0_8
  stopifnot(max(abs(fixed - u$expected_fixed_lod4_5)) < 1e-10,
            max(abs(comparison$fractionOverlapsWithStrand - u$observed)) < 1e-10,
            all(comparison$numberOfQTLs == u$n_observed_qtls),
            all(is.finite(u$expected_actual_first_pass)))
  plot_data = comparison
  plot_data$observed_pct = 100 * comparison$fractionOverlapsWithStrand
  plot_data$expected_pct = 100 * u$expected_actual_first_pass
  summary_row = data.frame(mean_expected = mean(u$expected_actual_first_pass),
    mean_observed = mean(u$observed),
    pearson_r = cor(u$observed, u$expected_actual_first_pass),
    spearman_rho = cor(u$observed, u$expected_actual_first_pass, method = "spearman"))
  if (effect_scale == 1) {
    lod4 = summary_row
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
      "90 directed individual-study pairs; actual LODs; Pearson r = %.2f, Spearman rho = %.2f",
      lod4$pearson_r,
      lod4$spearman_rho
    ),
    x = "Power-predicted same-direction rediscovery (%)",
    y = "Observed Figure 3 rediscovery (%)",
    color = "Comparison",
    size = "Query QTLs"
  ) +
  ggplot2::theme_bw()


    filename = "all_pairs_actual_thresholds.png"
  } else {
    scale_label = "80%"
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
          "%s source-effect sensitivity; actual LODs; 90 directed pairs\n",
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


    filename = "all_pairs_actual_thresholds_scale0_8.png"
  }
  ggplot2::ggsave(file.path(out, filename), plot = p, width = 8, height = 7, dpi = 300)
  print(summary_row)
}
