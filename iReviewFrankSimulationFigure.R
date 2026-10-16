# Direct phenotype simulation using the original validator functions.
args = commandArgs(TRUE)
out = args[[1]]
independent_rcrit = function(lod, n) {
  sqrt(pmax(0, 1 - 10 ^ (-2 * lod / n)))
}

independent_probability = function(rho, n, lod) {

  rcrit = independent_rcrit(lod = lod, n = n)
  df = n - 2
  tcrit = rcrit * sqrt(df / (1 - rcrit ^ 2))
  ncp = rho * sqrt((n - 1) / (1 - rho ^ 2))

  positive = stats::pt(
    tcrit,
    df = df,
    ncp = ncp,
    lower.tail = FALSE
  )
  negative = stats::pt(
    -tcrit,
    df = df,
    ncp = ncp,
    lower.tail = TRUE
  )

  c(
    same = positive,
    wrong = negative,
    any = positive + negative
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
  lod_hat = -n * log10(1 - r_hat ^ 2) / 2

  analytic = independent_probability(
    rho = rho,
    n = n,
    lod = lod
  )

  data.frame(
    n = n,
    rho = rho,
    lod = lod,
    n_rep = n_rep,
    analytic_any = analytic[["any"]],
    simulated_any = mean(lod_hat > lod),
    analytic_same = analytic[["same"]],
    simulated_same = mean(lod_hat > lod & r_hat > 0),
    stringsAsFactors = FALSE
  )
}


thresholds = read.delim("/PHShome_actual/k/kv1026/tools/eris-jobcards/actions/card_1e080f28d4b1c883cb4f66157701e1a7/thresholds/threshold_full_first_pass_by_study_chr.tsv")
stopifnot(nrow(thresholds) == 160L, all(is.finite(thresholds$lod)))
lods = as.numeric(quantile(thresholds$lod, c(0, 0.5, 1), names = FALSE))
stopifnot(length(unique(lods)) == 3L)
# Retain the original 3 sample sizes, 7 effect sizes and 1,000 replicates.
grid = do.call(rbind, lapply(c(33, 99, 1012), function(n) {
  d = expand.grid(rho = seq(0.10, 0.40, by = 0.05), lod = lods)
  d$n = n
  d
}))
simulation = do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  validate_direct_simulation(grid$n[i], grid$rho[i], grid$lod[i],
                             n_rep = 1000L, seed = 20260818L + 1000L + i)
}))
stopifnot(nrow(simulation) == 63L)
simulation$error = simulation$simulated_any - simulation$analytic_any
simulation$mc_se = sqrt(simulation$analytic_any * (1 - simulation$analytic_any) / simulation$n_rep)
simulation$within_tolerance = abs(simulation$error) <= 0.005 + 4 * simulation$mc_se
write.table(simulation, file.path(out, "simulation_actual_thresholds.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
simulation$sample_size = factor(simulation$n, levels = c(33, 99, 1012),
                                labels = c("33", "99", "1,012"))
simulation$threshold = factor(simulation$lod, levels = lods,
  labels = paste0(c("Low: ", "Median: ", "High: "), formatC(lods, digits = 3, format = "f")))
p = ggplot2::ggplot(simulation, ggplot2::aes(x = analytic_any, y = simulated_any,
                                          color = sample_size, shape = threshold)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  ggplot2::geom_point(size = 3, alpha = 0.85) +
  ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
  ggplot2::scale_shape_manual(values = c(16, 17, 15)) +
  ggplot2::scale_x_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.2)) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.2)) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(x = "Detection probability from power calculation",
    y = "Fraction detected in 1,000 phenotype simulations",
    color = "Sample size", shape = "LOD threshold") +
  ggplot2::theme_classic(base_size = 14)
ggplot2::ggsave(file.path(out, "simulation_actual_thresholds.png"),
               plot = p, width = 8, height = 6, dpi = 300)
caption = sprintf(paste0(
  "Simulation check of the power calculation described in the Second sentence. ",
  "Each point represents one combination of sample size (33, 99 or 1,012), ",
  "effect size (rho = 0.10 to 0.40 in steps of 0.05), and LOD threshold ",
  "(%.3f, %.3f or %.3f: the minimum, median and maximum of 160 recovered ",
  "first-pass study-by-chromosome thresholds). The x-axis gives the detection ",
  "probability calculated from the noncentral t distribution. For the y-axis, ",
  "we generated 1,000 phenotype datasets per setting by adding Gaussian noise ",
  "to a fixed, approximately balanced two-allele genotype effect, then calculated ",
  "the fraction whose genotype-phenotype LOD exceeded the same cutoff. ",
  "The dashed line marks equality. This checks single-locus detection power; ",
  "it does not repeat genome-wide mapping or re-estimate FDR thresholds. ",
  "Maximum absolute difference: %.2f percentage points."),
  lods[1], lods[2], lods[3], 100 * max(abs(simulation$error)))
writeLines(caption, file.path(out, "simulation_caption.txt"))
writeLines(c(sprintf("All 63 settings within Monte Carlo tolerance: %s",
  all(simulation$within_tolerance)), caption), file.path(out, "simulation_summary.txt"))
cat(caption, "\n")
