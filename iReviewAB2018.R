########################################
########################################
## GENETICS REVIEW AB2018 FEATURES    ##
## CANONICAL GENE CHARACTERISTICS     ##
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

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package readxl is required.")
}

####################
# SOURCE AB2018    #
####################

source_file =
  "raw/Albert_2018_gene_characteristics.xlsx"

source_sheet =
  "SI_Data_06_geneFeatures"

message(
  "Reading Albert et al. 2018 gene characteristics..."
)

ab2018 = readxl::read_excel(
  source_file,
  sheet = source_sheet,
  na = c(
    "",
    "NA"
  )
)

ab2018 = as.data.frame(
  ab2018,
  stringsAsFactors = FALSE
)

expected_columns = c(
  "gene",
  "expression",
  "essential",
  "dNdS",
  "dNdS_adjust",
  "PPI",
  "GxGLenient",
  "GxGIntermediate",
  "GxGStrict",
  "isTF",
  "hasHumanHomolog",
  "hasParalog"
)

missing_columns = setdiff(
  expected_columns,
  names(ab2018)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing expected AB2018 columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

if (nrow(ab2018) != 7126L) {
  stop(
    "Expected 7,126 AB2018 rows, found ",
    nrow(ab2018),
    "."
  )
}

if (anyDuplicated(ab2018$gene)) {
  stop(
    "AB2018 feature table contains duplicated genes."
  )
}

########################
# SUPPORTED GENE SET    #
########################

supported_analyte = readLines(
  paste0(
    imed,
    "supported_analyte.txt"
  )
)

if (length(supported_analyte) != 385L) {
  stop(
    "Expected 385 supported analyte genes, found ",
    length(supported_analyte),
    "."
  )
}

if (anyDuplicated(supported_analyte)) {
  stop(
    "supported_analyte contains duplicated genes."
  )
}

missing_supported = setdiff(
  supported_analyte,
  ab2018$gene
)

if (length(missing_supported) > 0L) {
  stop(
    "Supported genes absent from AB2018: ",
    paste(
      missing_supported,
      collapse = ", "
    )
  )
}

########################
# DERIVED MEMBERSHIP    #
########################

ab2018$in_supported_analyte =
  ab2018$gene %in%
  supported_analyte

ab2018$expression_measured =
  !is.na(
    ab2018$expression
  )

expression_universe =
  ab2018[
    ab2018$expression_measured,
    ,
    drop = FALSE
  ]

if (nrow(expression_universe) != 5720L) {
  stop(
    "Expected 5,720 expression-measured genes, found ",
    nrow(expression_universe),
    "."
  )
}

supported =
  expression_universe[
    expression_universe$in_supported_analyte,
    ,
    drop = FALSE
  ]

expressed_background =
  expression_universe[
    !expression_universe$in_supported_analyte,
    ,
    drop = FALSE
  ]

if (nrow(supported) != 385L) {
  stop(
    "Expected all 385 supported genes in expression universe."
  )
}

########################
# EXPRESSION POSITION   #
########################

expression_ecdf = ecdf(
  expression_universe$expression
)

ab2018$expression_percentile =
  NA_real_

ab2018$expression_percentile[
  ab2018$expression_measured
] =
  expression_ecdf(
    ab2018$expression[
      ab2018$expression_measured
    ]
  )

expression_rank = rank(
  expression_universe$expression,
  ties.method = "average"
)

expression_decile = ceiling(
  10 *
  expression_rank /
  nrow(expression_universe)
)

expression_decile =
  pmax(
    1L,
    pmin(
      10L,
      expression_decile
    )
  )

ab2018$expression_decile =
  NA_integer_

ab2018$expression_decile[
  match(
    expression_universe$gene,
    ab2018$gene
  )
] =
  expression_decile

expression_universe =
  ab2018[
    ab2018$expression_measured,
    ,
    drop = FALSE
  ]

supported =
  expression_universe[
    expression_universe$in_supported_analyte,
    ,
    drop = FALSE
  ]

expressed_background =
  expression_universe[
    !expression_universe$in_supported_analyte,
    ,
    drop = FALSE
  ]

########################
# FEATURE COVERAGE      #
########################

feature_columns = setdiff(
  expected_columns,
  "gene"
)

coverage_for_group = function(
  dat,
  group_name
) {

  data.frame(
    group = group_name,
    feature = feature_columns,
    n_group = nrow(dat),
    n_nonmissing = vapply(
      feature_columns,
      function(feature) {
        sum(
          !is.na(
            dat[[feature]]
          )
        )
      },
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
}

feature_coverage = rbind(
  coverage_for_group(
    ab2018,
    "all_annotated"
  ),
  coverage_for_group(
    expression_universe,
    "expression_universe"
  ),
  coverage_for_group(
    supported,
    "supported"
  ),
  coverage_for_group(
    expressed_background,
    "expressed_background"
  )
)

feature_coverage$fraction_nonmissing =
  feature_coverage$n_nonmissing /
  feature_coverage$n_group

########################
# SOURCE SUMMARY        #
########################

source_summary = data.frame(
  item = c(
    "source_file",
    "source_sheet",
    "n_all_annotated",
    "n_expression_universe",
    "n_supported",
    "n_expressed_background",
    "readxl_version"
  ),
  value = c(
    source_file,
    source_sheet,
    nrow(ab2018),
    nrow(expression_universe),
    nrow(supported),
    nrow(expressed_background),
    as.character(
      packageVersion(
        "readxl"
      )
    )
  ),
  stringsAsFactors = FALSE
)

########################
# SAVE CANONICAL OBJECT #
########################

review_ab2018 = list(
  source_summary = source_summary,
  all_gene_features = ab2018,
  expression_universe = expression_universe,
  supported = supported,
  expressed_background = expressed_background,
  feature_coverage = feature_coverage,
  supported_genes = supported_analyte
)

saveRDS(
  review_ab2018,
  paste0(
    ab2018_imed,
    "gene_characteristics_clean.rds"
  )
)

writeLines(
  expression_universe$gene,
  paste0(
    ab2018_imed,
    "expression_universe_genes.txt"
  )
)

writeLines(
  expressed_background$gene,
  paste0(
    ab2018_imed,
    "expressed_background_genes.txt"
  )
)

writeLines(
  supported$gene,
  paste0(
    ab2018_imed,
    "supported_genes.txt"
  )
)

write.table(
  feature_coverage,
  file = paste0(
    genes_asys,
    "AB2018_feature_coverage.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  source_summary,
  file = paste0(
    genes_asys,
    "AB2018_source_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  expression_universe,
  file = paste0(
    genes_asys,
    "AB2018_expression_universe.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("")
message("AB2018 preparation complete.")
message(
  "All annotated genes: ",
  nrow(ab2018)
)
message(
  "Expression universe: ",
  nrow(expression_universe)
)
message(
  "Supported genes: ",
  nrow(supported)
)
message(
  "Expressed comparison genes: ",
  nrow(expressed_background)
)
