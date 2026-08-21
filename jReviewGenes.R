########################################
########################################
## GENETICS REVIEW GENE SELECTION     ##
## REPRESENTATIVENESS OF 385 GENES    ##
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

ab2018_imed = paste0(
  review_imed,
  "ab2018/"
)

genes_asys = paste0(
  review_asys,
  "genes/"
)

ensure_directory(review_imed)
ensure_directory(review_asys)
ensure_directory(ab2018_imed)
ensure_directory(genes_asys)

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package openxlsx is required.")
}

review_ab2018 = readRDS(
  paste0(
    ab2018_imed,
    "gene_characteristics_clean.rds"
  )
)

universe =
  review_ab2018$expression_universe

supported =
  review_ab2018$supported

background =
  review_ab2018$expressed_background

########################
# FEATURE DEFINITIONS   #
########################

continuous_features = c(
  "expression",
  "dNdS",
  "dNdS_adjust",
  "PPI",
  "GxGLenient",
  "GxGIntermediate",
  "GxGStrict"
)

binary_features = c(
  "essential",
  "isTF",
  "hasHumanHomolog",
  "hasParalog"
)

feature_order = c(
  "expression",
  "essential",
  "hasHumanHomolog",
  "hasParalog",
  "isTF",
  "PPI",
  "GxGStrict",
  "GxGIntermediate",
  "GxGLenient",
  "dNdS_adjust",
  "dNdS"
)

feature_labels = c(
  expression = "Expression",
  essential = "Essential gene",
  hasHumanHomolog = "Human homolog",
  hasParalog = "Paralog",
  isTF = "Transcription factor",
  PPI = "Protein-protein interactions",
  GxGStrict = "Genetic interactions, strict",
  GxGIntermediate = "Genetic interactions, intermediate",
  GxGLenient = "Genetic interactions, lenient",
  dNdS_adjust = "Adjusted dN/dS",
  dNdS = "dN/dS"
)

# Keep the main-text table compact by retaining one representative
# genetic-interaction measure and one dN/dS measure. All features above
# are still tested individually and retained in the full output.
main_features = c(
  "expression",
  "essential",
  "hasHumanHomolog",
  "hasParalog",
  "isTF",
  "PPI",
  "GxGStrict",
  "dNdS_adjust"
)

########################
# INDIVIDUAL TESTS      #
########################

test_continuous = function(feature) {

  x = supported[[feature]]
  y = background[[feature]]

  x = x[
    !is.na(x)
  ]

  y = y[
    !is.na(y)
  ]

  test = suppressWarnings(
    wilcox.test(
      x,
      y,
      exact = FALSE
    )
  )

  u = as.numeric(
    test$statistic
  )

  rank_biserial =
    2 *
    (
      u /
      (
        length(x) *
        length(y)
      )
    ) -
    1

  data.frame(
    feature = feature,
    feature_type = "continuous",
    test = "Wilcoxon rank-sum",
    effect_metric = "Rank-biserial r",
    n_supported = length(x),
    n_background = length(y),
    supported_n_true = NA_integer_,
    background_n_true = NA_integer_,
    supported_prop = NA_real_,
    background_prop = NA_real_,
    supported_median = median(x),
    supported_q25 = unname(
      quantile(
        x,
        0.25
      )
    ),
    supported_q75 = unname(
      quantile(
        x,
        0.75
      )
    ),
    background_median = median(y),
    background_q25 = unname(
      quantile(
        y,
        0.25
      )
    ),
    background_q75 = unname(
      quantile(
        y,
        0.75
      )
    ),
    effect_estimate = rank_biserial,
    effect_low = NA_real_,
    effect_high = NA_real_,
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

test_binary = function(feature) {

  x = supported[[feature]]
  y = background[[feature]]

  x = x[
    !is.na(x)
  ]

  y = y[
    !is.na(y)
  ]

  table_2x2 = matrix(
    c(
      sum(x),
      sum(!x),
      sum(y),
      sum(!y)
    ),
    nrow = 2,
    byrow = TRUE
  )

  test = fisher.test(
    table_2x2
  )

  data.frame(
    feature = feature,
    feature_type = "binary",
    test = "Fisher exact",
    effect_metric = "Odds ratio",
    n_supported = length(x),
    n_background = length(y),
    supported_n_true = sum(x),
    background_n_true = sum(y),
    supported_prop = mean(x),
    background_prop = mean(y),
    supported_median = NA_real_,
    supported_q25 = NA_real_,
    supported_q75 = NA_real_,
    background_median = NA_real_,
    background_q25 = NA_real_,
    background_q75 = NA_real_,
    effect_estimate = unname(
      test$estimate
    ),
    effect_low = unname(
      test$conf.int[1]
    ),
    effect_high = unname(
      test$conf.int[2]
    ),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

individual_tests = rbind(
  do.call(
    rbind,
    lapply(
      continuous_features,
      test_continuous
    )
  ),
  do.call(
    rbind,
    lapply(
      binary_features,
      test_binary
    )
  )
)

individual_tests$q_value = p.adjust(
  individual_tests$p_value,
  method = "BH"
)

individual_tests$feature_label =
  unname(
    feature_labels[
      individual_tests$feature
    ]
  )

individual_tests$order =
  match(
    individual_tests$feature,
    feature_order
  )

individual_tests =
  individual_tests[
    order(
      individual_tests$order
    ),
    ,
    drop = FALSE
  ]

individual_tests$order = NULL

########################
# REVIEW-READY TABLE    #
########################

format_p = function(x) {
  ifelse(
    x < 0.001,
    formatC(
      x,
      format = "e",
      digits = 2
    ),
    sprintf(
      "%.3f",
      x
    )
  )
}

main_table =
  individual_tests[
    individual_tests$feature %in%
    main_features,
    ,
    drop = FALSE
  ]

main_table$order =
  match(
    main_table$feature,
    main_features
  )

main_table =
  main_table[
    order(
      main_table$order
    ),
    ,
    drop = FALSE
  ]

main_table$supported_genes = ifelse(
  main_table$feature_type == "binary",
  sprintf(
    "%d/%d (%.1f%%)",
    main_table$supported_n_true,
    main_table$n_supported,
    100 * main_table$supported_prop
  ),
  sprintf(
    "%.3f [%.3f, %.3f] (n=%d)",
    main_table$supported_median,
    main_table$supported_q25,
    main_table$supported_q75,
    main_table$n_supported
  )
)

main_table$other_genes = ifelse(
  main_table$feature_type == "binary",
  sprintf(
    "%d/%d (%.1f%%)",
    main_table$background_n_true,
    main_table$n_background,
    100 * main_table$background_prop
  ),
  sprintf(
    "%.3f [%.3f, %.3f] (n=%d)",
    main_table$background_median,
    main_table$background_q25,
    main_table$background_q75,
    main_table$n_background
  )
)

main_table$effect = ifelse(
  main_table$feature_type == "binary",
  sprintf(
    "OR %.2f [%.2f, %.2f]",
    main_table$effect_estimate,
    main_table$effect_low,
    main_table$effect_high
  ),
  sprintf(
    "r = %.3f",
    main_table$effect_estimate
  )
)

main_table$p =
  format_p(
    main_table$p_value
  )

main_table$BH_FDR =
  format_p(
    main_table$q_value
  )

main_table =
  main_table[
    ,
    c(
      "feature_label",
      "supported_genes",
      "other_genes",
      "test",
      "effect_metric",
      "effect",
      "p",
      "BH_FDR"
    ),
    drop = FALSE
  ]

names(main_table) = c(
  "Feature",
  "Supported genes",
  "Other genes",
  "Test",
  "Effect metric",
  "Effect",
  "P",
  "BH FDR"
)

########################
# SAVE OUTPUTS          #
########################

write.table(
  main_table,
  file = paste0(
    genes_asys,
    "reviewer2_gene_characteristics_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  individual_tests,
  file = paste0(
    genes_asys,
    "gene_characteristics_individual_tests.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Keep the legacy unadjusted filename synchronized for compatibility with
# any downstream notes generated before this simplification.
write.table(
  individual_tests,
  file = paste0(
    genes_asys,
    "gene_characteristics_unadjusted.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

results = list(
  main_table = main_table,
  individual_tests = individual_tests,
  feature_coverage =
    review_ab2018$feature_coverage,
  source_summary =
    review_ab2018$source_summary
)

saveRDS(
  results,
  paste0(
    ab2018_imed,
    "reviewer2_gene_characteristics_results.rds"
  )
)

########################
# REVIEW WORKBOOK       #
########################

readme = data.frame(
  item = c(
    "Primary universe",
    "Supported gene set",
    "Continuous features",
    "Binary features",
    "Multiple-testing correction",
    "Main-table selection"
  ),
  description = c(
    "5,720 genes with AB2018 expression measurements",
    "385 genes measured in at least two mRNA and two protein datasets",
    "Wilcoxon rank-sum; effect size is rank-biserial correlation",
    "Fisher exact test; effect size is odds ratio with 95% confidence interval",
    "Benjamini-Hochberg across all 11 individually tested gene characteristics",
    "Main table removes redundant GxG thresholds and raw dN/dS; all 11 tests remain in the full output"
  ),
  stringsAsFactors = FALSE
)

workbook_tables = list(
  README = readme,
  main_table = main_table,
  all_individual_tests = individual_tests,
  feature_coverage =
    review_ab2018$feature_coverage
)

openxlsx::write.xlsx(
  workbook_tables,
  file = paste0(
    genes_asys,
    "reviewer2_gene_characteristics.xlsx"
  ),
  asTable = TRUE,
  overwrite = TRUE
)

message("")
message("Reviewer 2 gene-characteristic analysis complete.")
message("Only individual feature tests were run.")
message(
  "Output directory: ",
  genes_asys
)
message("")
message("Main review-ready table:")
print(
  main_table,
  row.names = FALSE
)
message("")
message("All individual tests:")
print(
  individual_tests[
    ,
    c(
      "feature_label",
      "feature_type",
      "test",
      "effect_metric",
      "effect_estimate",
      "p_value",
      "q_value"
    )
  ],
  row.names = FALSE
)
