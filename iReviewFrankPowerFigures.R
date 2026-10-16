# Usage: Rscript --vanilla iReviewFrankPowerFigures.R POWER_OUTPUT
args = commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out = normalizePath(args[[1]], mustWork = TRUE)
d = read.delim(file.path(out, "pairwise_actual_threshold_power.tsv"))
d = d[d$denominator == "figure3_matched", ]
stopifnot(nrow(d) == 180L)
protein = c("Khan_2009", "Foss_2011", "Picotti_2013", "GrossbachPRT_2022")
d$comparison = paste(ifelse(d$query %in% protein, "Protein", "mRNA"),
                     ifelse(d$target %in% protein, "Protein", "mRNA"), sep = " to ")
colors = c("mRNA to mRNA" = "#0072B2", "mRNA to Protein" = "#E69F00",
           "Protein to mRNA" = "#009E73", "Protein to Protein" = "#CC79A7")
draw = function(albert = FALSE) {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 4, 1), oma = c(3, 0, 1, 0))
  for (scale in c(1, 0.8)) {
    z = d[d$effect_scale == scale, ]
    if (albert) z = z[z$query == "Albert_2018", ]
    stopifnot(nrow(z) == if (albert) 9L else 90L)
    x = 100 * z$expected_actual_first_pass
    y = 100 * z$observed
    lim = if (albert) c(0, max(10, ceiling(max(c(x, y)) * 1.25))) else c(0, 100)
    plot(x, y, xlim = lim, ylim = lim, pch = 19,
         col = colors[z$comparison], cex = 1.1,
         xlab = "Predicted same-direction detection (%)",
         ylab = "Observed QTL rediscovery (%)",
         main = if (scale == 1) "Observed source effects" else "Source effects reduced by 20%")
    abline(0, 1, lty = 2, col = "grey55")
    mtext(sprintf("n = %d; Pearson r = %.2f; Spearman rho = %.2f",
                  nrow(z), cor(x, y), cor(x, y, method = "spearman")),
          side = 3, line = 0.3, cex = 0.75)
    if (albert) text(x, y, labels = sub("_", " ", z$target), pos = 3, cex = 0.6)
  }
  par(fig = c(0, 1, 0, 0.07), new = TRUE, mar = rep(0, 4))
  plot.new()
  legend("center", legend = names(colors), col = colors, pch = 19,
         horiz = TRUE, bty = "n", cex = 0.85)
}
for (albert in c(FALSE, TRUE)) {
  stem = if (albert) "albert_nine_actual_thresholds" else "all90_actual_thresholds"
  png(file.path(out, paste0(stem, ".png")), width = 2400, height = 1200, res = 200)
  draw(albert)
  dev.off()
  pdf(file.path(out, paste0(stem, ".pdf")), width = 12, height = 6)
  draw(albert)
  dev.off()
}
writeLines(c(
  "Predicted versus observed QTL rediscovery using target-study, chromosome-specific first-pass remapping LOD thresholds.",
  "Each point is a directed study comparison; colors indicate source and target molecular layers.",
  "Left: observed source effects. Right: effects multiplied by 0.8 (sensitivity analysis, not a fitted correction).",
  "The dashed line indicates equality. Both axes use Figure 3 matched source-QTL denominators.",
  "Prediction is a fixed-effect detection benchmark; it does not replay conditional mapping or model localization."
), file.path(out, "figure_captions.txt"))
