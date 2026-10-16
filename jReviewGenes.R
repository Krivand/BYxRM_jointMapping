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

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package readxl is required.")
}

genes_asys = paste0(asys, "review/genes/")
ensure_directory(genes_asys)

########################
# LOAD AB2018 FEATURES #
########################

source_file = "raw/Albert_2018_gene_characteristics.xlsx"
source_sheet = "SI_Data_06_geneFeatures"

features = readxl::read_excel(
  source_file,
  sheet = source_sheet,
  na = c("", "NA")
) %>%
  as.data.frame(stringsAsFactors = FALSE)

required_columns = c(
  "gene",
  "expression",
  "essential",
  "hasHumanHomolog",
  "PPI"
)

missing_columns = setdiff(
  required_columns,
  names(features)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing expected AB2018 columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

features = features[, required_columns, drop = FALSE]

if (anyDuplicated(features$gene)) {
  stop("AB2018 feature table contains duplicated genes.")
}

supported_genes = readLines(
  paste0(imed, "supported_analyte.txt")
)

if (
  length(supported_genes) != 385L ||
  anyDuplicated(supported_genes)
) {
  stop("Expected 385 unique supported analyte genes.")
}

expression_universe = features[
  !is.na(features$expression),
  ,
  drop = FALSE
]

if (nrow(expression_universe) != 5720L) {
  stop(
    "Expected 5,720 genes with AB2018 expression values, found ",
    nrow(expression_universe),
    "."
  )
}

if (!all(supported_genes %in% expression_universe$gene)) {
  stop("At least one supported gene is absent from the expression universe.")
}

supported = expression_universe[
  expression_universe$gene %in% supported_genes,
  ,
  drop = FALSE
]

background = expression_universe[
  !(expression_universe$gene %in% supported_genes),
  ,
  drop = FALSE
]

if (
  nrow(supported) != 385L ||
  nrow(background) != 5335L
) {
  stop("Unexpected supported/background gene counts.")
}

####################
# HELPER FUNCTIONS #
####################

format_p = function(x) {
  ifelse(
    x < 0.01,
    formatC(x, format = "e", digits = 2),
    sprintf("%.2f", x)
  )
}

continuous_row = function(
  feature,
  label,
  fold_scale = c("linear", "log2")
) {

  fold_scale = match.arg(fold_scale)

  x = as.numeric(supported[[feature]])
  y = as.numeric(background[[feature]])

  x = x[!is.na(x)]
  y = y[!is.na(y)]

  test = suppressWarnings(
    wilcox.test(x, y, exact = FALSE)
  )

  median_x = median(x)
  median_y = median(y)

  if (fold_scale == "log2") {
    median_difference = median_x - median_y
    fold_difference = 2^abs(median_difference)
    direction = ifelse(
      median_difference >= 0,
      "higher",
      "lower"
    )
  } else {
    if (median_x <= 0 || median_y <= 0) {
      stop("Linear-scale fold differences require positive medians.")
    }
    fold_difference = max(median_x, median_y) / min(median_x, median_y)
    direction = ifelse(
      median_x >= median_y,
      "higher",
      "lower"
    )
  }

  effect_size = sprintf(
    "%.2f-fold %s",
    fold_difference,
    direction
  )

  data.frame(
    Characteristic = label,
    `Supported genes` = sprintf("%.2f", median_x),
    `Other genes` = sprintf("%.2f", median_y),
    `Effect size` = effect_size,
    P = format_p(test$p.value),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

binary_row = function(feature, label) {

  x = supported[[feature]]
  y = background[[feature]]

  x = as.logical(x[!is.na(x)])
  y = as.logical(y[!is.na(y)])

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

  test = fisher.test(table_2x2)

  data.frame(
    Characteristic = label,
    `Supported genes` = sprintf(
      "%d/%d",
      sum(x),
      length(x)
    ),
    `Other genes` = sprintf(
      "%d/%d",
      sum(y),
      length(y)
    ),
    `Effect size` = sprintf(
      "OR %.2f",
      unname(test$estimate)
    ),
    P = format_p(test$p.value),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

########################
# REVIEW-READY TABLE    #
########################

review_table = rbind(
  continuous_row(
    "expression",
    "Mean log2 TPM",
    fold_scale = "log2"
  ),
  binary_row(
    "essential",
    "Essential gene"
  ),
  binary_row(
    "hasHumanHomolog",
    "Human homolog"
  ),
  continuous_row(
    "PPI",
    "PPI partners",
    fold_scale = "linear"
  )
)

################
# SAVE OUTPUT  #
################

output_file = paste0(
  genes_asys,
  "supported_gene_characteristics_summary.tsv"
)

write.table(
  review_table,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("")
message("Supported-gene representativeness analysis complete.")
message("Output: ", output_file)
message("")
print(review_table, row.names = FALSE)
