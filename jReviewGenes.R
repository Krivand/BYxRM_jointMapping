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

count_features = c(
  "PPI",
  "GxGLenient",
  "GxGIntermediate",
  "GxGStrict"
)

########################
# UNADJUSTED TESTS      #
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
    effect_metric = "rank_biserial",
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
    effect_metric = "odds_ratio",
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

unadjusted = rbind(
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

unadjusted$q_value = p.adjust(
  unadjusted$p_value,
  method = "BH"
)

unadjusted$feature_label =
  unname(
    feature_labels[
      unadjusted$feature
    ]
  )

unadjusted$order =
  match(
    unadjusted$feature,
    feature_order
  )

unadjusted =
  unadjusted[
    order(
      unadjusted$order
    ),
    ,
    drop = FALSE
  ]

unadjusted$order = NULL

########################
# EXPRESSION ADJUSTMENT #
########################

adjusted_models = list()
adjusted_rows = list()

adjustment_features = setdiff(
  feature_order,
  "expression"
)

for (feature in adjustment_features) {

  dat = universe[
    ,
    c(
      "in_supported_analyte",
      "expression",
      feature
    ),
    drop = FALSE
  ]

  names(dat)[3] =
    "raw_feature"

  dat = dat[
    complete.cases(dat),
    ,
    drop = FALSE
  ]

  if (feature %in% binary_features) {

    dat$feature_value =
      as.integer(
        dat$raw_feature
      )

    feature_scale =
      "TRUE vs FALSE"

  } else {

    transformed =
      dat$raw_feature

    if (feature %in% count_features) {
      transformed =
        log1p(
          transformed
        )
    }

    dat$feature_value =
      as.numeric(
        scale(
          transformed
        )
      )

    if (feature %in% count_features) {
      feature_scale =
        "per SD of log1p(feature)"
    } else {
      feature_scale =
        "per SD of feature"
    }
  }

  fit = glm(
    in_supported_analyte ~
      splines::ns(
        expression,
        df = 4
      ) +
      feature_value,
    data = dat,
    family = binomial()
  )

  coefficient_table =
    summary(fit)$coefficients

  beta =
    coefficient_table[
      "feature_value",
      "Estimate"
    ]

  se =
    coefficient_table[
      "feature_value",
      "Std. Error"
    ]

  p_value =
    coefficient_table[
      "feature_value",
      "Pr(>|z|)"
    ]

  adjusted_models[[feature]] =
    fit

  adjusted_rows[[feature]] =
    data.frame(
      feature = feature,
      feature_label = unname(
        feature_labels[
          feature
        ]
      ),
      feature_type = ifelse(
        feature %in% binary_features,
        "binary",
        "continuous"
      ),
      feature_scale = feature_scale,
      n_model = nrow(dat),
      n_supported = sum(
        dat$in_supported_analyte
      ),
      n_background = sum(
        !dat$in_supported_analyte
      ),
      log_odds_estimate = beta,
      standard_error = se,
      odds_ratio = exp(beta),
      odds_ratio_low = exp(
        beta -
        1.96 *
        se
      ),
      odds_ratio_high = exp(
        beta +
        1.96 *
        se
      ),
      p_value = p_value,
      stringsAsFactors = FALSE
    )
}

expression_adjusted = do.call(
  rbind,
  adjusted_rows
)

rownames(expression_adjusted) =
  NULL

expression_adjusted$q_value =
  p.adjust(
    expression_adjusted$p_value,
    method = "BH"
  )

expression_adjusted$order =
  match(
    expression_adjusted$feature,
    feature_order
  )

expression_adjusted =
  expression_adjusted[
    order(
      expression_adjusted$order
    ),
    ,
    drop = FALSE
  ]

expression_adjusted$order =
  NULL

########################
# CORE MULTIVARIABLE    #
########################

core_features = c(
  "essential",
  "hasHumanHomolog",
  "hasParalog",
  "isTF"
)

core_dat =
  universe[
    ,
    c(
      "in_supported_analyte",
      "expression",
      core_features
    ),
    drop = FALSE
  ]

for (feature in core_features) {

  core_dat[[feature]] =
    as.integer(
      core_dat[[feature]]
    )
}

core_fit = glm(
  in_supported_analyte ~
    splines::ns(
      expression,
      df = 4
    ) +
    essential +
    hasHumanHomolog +
    hasParalog +
    isTF,
  data = core_dat,
  family = binomial()
)

core_coefficients =
  summary(core_fit)$coefficients

core_rows = lapply(
  core_features,
  function(feature) {

    beta =
      core_coefficients[
        feature,
        "Estimate"
      ]

    se =
      core_coefficients[
        feature,
        "Std. Error"
      ]

    data.frame(
      feature = feature,
      feature_label = unname(
        feature_labels[
          feature
        ]
      ),
      n_model = nrow(core_dat),
      odds_ratio = exp(beta),
      odds_ratio_low = exp(
        beta -
        1.96 *
        se
      ),
      odds_ratio_high = exp(
        beta +
        1.96 *
        se
      ),
      p_value = core_coefficients[
        feature,
        "Pr(>|z|)"
      ],
      stringsAsFactors = FALSE
    )
  }
)

core_model_table =
  do.call(
    rbind,
    core_rows
  )

core_model_table$q_value =
  p.adjust(
    core_model_table$p_value,
    method = "BH"
  )

########################
# EXPRESSION SUMMARY    #
########################

expression_row =
  unadjusted[
    unadjusted$feature ==
    "expression",
    ,
    drop = FALSE
  ]

expression_summary = data.frame(
  n_expression_universe = nrow(universe),
  n_supported = nrow(supported),
  n_background = nrow(background),
  universe_expression_median =
    median(
      universe$expression
    ),
  supported_expression_median =
    median(
      supported$expression
    ),
  background_expression_median =
    median(
      background$expression
    ),
  supported_median_percentile =
    median(
      supported$expression_percentile
    ),
  supported_fraction_above_50th =
    mean(
      supported$expression_percentile >
      0.50
    ),
  supported_fraction_above_75th =
    mean(
      supported$expression_percentile >
      0.75
    ),
  supported_fraction_above_90th =
    mean(
      supported$expression_percentile >
      0.90
    ),
  expression_rank_biserial =
    expression_row$effect_estimate,
  expression_wilcoxon_p =
    expression_row$p_value,
  stringsAsFactors = FALSE
)

########################
# EXPRESSION DECILES    #
########################

expression_deciles = do.call(
  rbind,
  lapply(
    1:10,
    function(decile) {

      this_decile =
        universe[
          universe$expression_decile ==
          decile,
          ,
          drop = FALSE
        ]

      data.frame(
        expression_decile = decile,
        n_universe = nrow(
          this_decile
        ),
        n_supported = sum(
          this_decile$in_supported_analyte
        ),
        n_background = sum(
          !this_decile$in_supported_analyte
        ),
        fraction_of_decile_supported =
          mean(
            this_decile$in_supported_analyte
          ),
        fraction_of_supported_genes =
          sum(
            this_decile$in_supported_analyte
          ) /
          nrow(supported),
        stringsAsFactors = FALSE
      )
    }
  )
)

##############################
# EXPRESSION-STANDARDIZED    #
# BINARY CHARACTERISTICS     #
##############################

binary_standardized_rows = list()

for (feature in binary_features) {

  supported_complete =
    supported[
      !is.na(
        supported[[feature]]
      ),
      ,
      drop = FALSE
    ]

  decile_rows = lapply(
    1:10,
    function(decile) {

      s =
        supported_complete[
          supported_complete$expression_decile ==
          decile,
          ,
          drop = FALSE
        ]

      b =
        background[
          background$expression_decile ==
          decile &
          !is.na(
            background[[feature]]
          ),
          ,
          drop = FALSE
        ]

      data.frame(
        expression_decile = decile,
        n_supported = nrow(s),
        n_background = nrow(b),
        background_rate = ifelse(
          nrow(b) > 0L,
          mean(
            b[[feature]]
          ),
          NA_real_
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  decile_table =
    do.call(
      rbind,
      decile_rows
    )

  if (
    any(
      decile_table$n_supported > 0L &
      is.na(
        decile_table$background_rate
      )
    )
  ) {
    standardized_rate =
      NA_real_
  } else {
    standardized_rate =
      sum(
        decile_table$n_supported *
        decile_table$background_rate,
        na.rm = TRUE
      ) /
      sum(
        decile_table$n_supported
      )
  }

  observed_rate =
    mean(
      supported_complete[[feature]]
    )

  binary_standardized_rows[[feature]] =
    data.frame(
      feature = feature,
      feature_label = unname(
        feature_labels[
          feature
        ]
      ),
      n_supported = nrow(
        supported_complete
      ),
      observed_supported_rate =
        observed_rate,
      expected_rate_given_expression =
        standardized_rate,
      observed_minus_expected =
        observed_rate -
        standardized_rate,
      observed_over_expected =
        observed_rate /
        standardized_rate,
      stringsAsFactors = FALSE
    )
}

binary_expression_standardized =
  do.call(
    rbind,
    binary_standardized_rows
  )

rownames(
  binary_expression_standardized
) = NULL

########################
# REVIEWER SUMMARY      #
########################

reviewer_summary =
  unadjusted

reviewer_summary$supported_summary =
  ifelse(
    reviewer_summary$feature_type ==
    "binary",
    sprintf(
      "%d/%d (%.1f%%)",
      reviewer_summary$supported_n_true,
      reviewer_summary$n_supported,
      100 *
      reviewer_summary$supported_prop
    ),
    sprintf(
      "%.3f [%.3f, %.3f]",
      reviewer_summary$supported_median,
      reviewer_summary$supported_q25,
      reviewer_summary$supported_q75
    )
  )

reviewer_summary$background_summary =
  ifelse(
    reviewer_summary$feature_type ==
    "binary",
    sprintf(
      "%d/%d (%.1f%%)",
      reviewer_summary$background_n_true,
      reviewer_summary$n_background,
      100 *
      reviewer_summary$background_prop
    ),
    sprintf(
      "%.3f [%.3f, %.3f]",
      reviewer_summary$background_median,
      reviewer_summary$background_q25,
      reviewer_summary$background_q75
    )
  )

adjusted_index = match(
  reviewer_summary$feature,
  expression_adjusted$feature
)

reviewer_summary$expression_adjusted_n_model =
  expression_adjusted$n_model[
    adjusted_index
  ]

reviewer_summary$expression_adjusted_scale =
  expression_adjusted$feature_scale[
    adjusted_index
  ]

reviewer_summary$expression_adjusted_or =
  expression_adjusted$odds_ratio[
    adjusted_index
  ]

reviewer_summary$expression_adjusted_or_low =
  expression_adjusted$odds_ratio_low[
    adjusted_index
  ]

reviewer_summary$expression_adjusted_or_high =
  expression_adjusted$odds_ratio_high[
    adjusted_index
  ]

reviewer_summary$expression_adjusted_p =
  expression_adjusted$p_value[
    adjusted_index
  ]

reviewer_summary$expression_adjusted_q =
  expression_adjusted$q_value[
    adjusted_index
  ]

reviewer_summary =
  reviewer_summary[
    ,
    c(
      "feature",
      "feature_label",
      "feature_type",
      "n_supported",
      "n_background",
      "supported_summary",
      "background_summary",
      "effect_metric",
      "effect_estimate",
      "effect_low",
      "effect_high",
      "p_value",
      "q_value",
      "expression_adjusted_n_model",
      "expression_adjusted_scale",
      "expression_adjusted_or",
      "expression_adjusted_or_low",
      "expression_adjusted_or_high",
      "expression_adjusted_p",
      "expression_adjusted_q"
    )
  ]

########################
# SAVE OUTPUTS          #
########################

write.table(
  reviewer_summary,
  file = paste0(
    genes_asys,
    "reviewer2_gene_characteristics_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  unadjusted,
  file = paste0(
    genes_asys,
    "gene_characteristics_unadjusted.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  expression_adjusted,
  file = paste0(
    genes_asys,
    "gene_characteristics_expression_adjusted.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  core_model_table,
  file = paste0(
    genes_asys,
    "gene_characteristics_core_model.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  binary_expression_standardized,
  file = paste0(
    genes_asys,
    "gene_characteristics_binary_expression_standardized.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  expression_summary,
  file = paste0(
    genes_asys,
    "gene_characteristics_expression_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  expression_deciles,
  file = paste0(
    genes_asys,
    "gene_characteristics_expression_deciles.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

results = list(
  reviewer_summary = reviewer_summary,
  unadjusted = unadjusted,
  expression_adjusted = expression_adjusted,
  core_model_table = core_model_table,
  binary_expression_standardized =
    binary_expression_standardized,
  expression_summary = expression_summary,
  expression_deciles = expression_deciles,
  feature_coverage =
    review_ab2018$feature_coverage,
  adjusted_models = adjusted_models,
  core_model = core_fit
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
    "Unadjusted continuous tests",
    "Unadjusted binary tests",
    "Expression adjustment",
    "Count-feature adjustment",
    "Multiple-testing correction"
  ),
  description = c(
    "5,720 genes with AB2018 expression measurements",
    "385 genes measured in at least two mRNA and two protein datasets",
    "Wilcoxon rank-sum; effect size is rank-biserial correlation",
    "Fisher exact test; effect size is odds ratio",
    "Logistic regression of supported membership with expression modeled using a natural spline with 4 df",
    "PPI and GxG counts are log1p-transformed and standardized before adjusted logistic regression",
    "Benjamini-Hochberg within each reported feature-test family"
  ),
  stringsAsFactors = FALSE
)

workbook_tables = list(
  README = readme,
  source_summary = review_ab2018$source_summary,
  reviewer_summary = reviewer_summary,
  expression_summary = expression_summary,
  expression_deciles = expression_deciles,
  unadjusted_tests = unadjusted,
  expression_adjusted = expression_adjusted,
  binary_standardized =
    binary_expression_standardized,
  core_model = core_model_table,
  feature_coverage =
    review_ab2018$feature_coverage,
  supported_genes = supported,
  gene_level = universe
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
message(
  "Output directory: ",
  genes_asys
)
message("")
print(
  reviewer_summary
)
