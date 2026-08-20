########################################
########################################
## GENETICS REVIEW MISSINGNESS AUDIT  ##
## PHENOTYPE FILTERING                ##
########################################
########################################

#################
# PRELIMINARIES #
#################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(
  asys,
  "review/missingness/"
)

ensure_directory(review_asys)

summary_rows = list()
removed_row_rows = list()
removed_col_rows = list()
threshold_check_rows = list()

###################
# AUDIT FUNCTIONS #
###################

append_removed = function(
  destination,
  dataset,
  margin,
  unit,
  values
) {

  if (length(values) == 0L) {
    return(destination)
  }

  destination[[length(destination) + 1L]] = data.frame(
    dataset = dataset,
    margin = margin,
    unit = unit,
    id = as.character(values),
    stringsAsFactors = FALSE
  )

  destination
}

audit_rmIterNA = function(
  x,
  dataset,
  row_unit,
  col_unit
) {

  if (is.null(dim(x)) || length(dim(x)) != 2L) {
    stop("Expected a two-dimensional object for ", dataset, ".")
  }

  before_rows = rownames(x)
  before_cols = colnames(x)
  before_dim = dim(x)

  filtered = rmIterNA(
    x,
    t = threshNA
  )

  after_rows = rownames(filtered)
  after_cols = colnames(filtered)
  after_dim = dim(filtered)

  removed_rows =
    if (!is.null(before_rows) && !is.null(after_rows)) {
      before_rows[!(before_rows %in% after_rows)]
    } else {
      character(0)
    }

  removed_cols =
    if (!is.null(before_cols) && !is.null(after_cols)) {
      before_cols[!(before_cols %in% after_cols)]
    } else {
      character(0)
    }

  summary_rows[[length(summary_rows) + 1L]] <<- data.frame(
    dataset = dataset,
    filter_applied = TRUE,
    row_unit = row_unit,
    col_unit = col_unit,
    n_rows_before = before_dim[1],
    n_rows_after = after_dim[1],
    n_rows_removed = before_dim[1] - after_dim[1],
    n_cols_before = before_dim[2],
    n_cols_after = after_dim[2],
    n_cols_removed = before_dim[2] - after_dim[2],
    stringsAsFactors = FALSE
  )

  removed_row_rows <<- append_removed(
    removed_row_rows,
    dataset = dataset,
    margin = "row",
    unit = row_unit,
    values = removed_rows
  )

  removed_col_rows <<- append_removed(
    removed_col_rows,
    dataset = dataset,
    margin = "column",
    unit = col_unit,
    values = removed_cols
  )

  filtered
}

audit_threshold_only = function(
  x,
  dataset,
  row_unit,
  col_unit
) {

  x = as.matrix(x)

  row_frac = rowMeans(is.na(x))
  col_frac = colMeans(is.na(x))

  row_hits = which(row_frac >= threshNA)
  col_hits = which(col_frac >= threshNA)

  threshold_check_rows[[length(threshold_check_rows) + 1L]] <<- data.frame(
    dataset = dataset,
    filter_applied = FALSE,
    row_unit = row_unit,
    col_unit = col_unit,
    n_rows = nrow(x),
    n_rows_at_or_above_25pct = length(row_hits),
    n_cols = ncol(x),
    n_cols_at_or_above_25pct = length(col_hits),
    stringsAsFactors = FALSE
  )

  invisible(NULL)
}

prepare_brem_gds = function(
  filename
) {

  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop(
      "GEOquery is required to audit Brem SOFT files. ",
      "Install/load it in the project environment and rerun."
    )
  }

  datGDS = GEOquery::getGEO(
    filename = paste0(
      raw,
      filename
    )
  ) %>%
    GEOquery::Table()

  datGDS[, 3:ncol(datGDS)] = apply(
    as.matrix(
      datGDS[, -c(1:2)]
    ),
    2,
    function(x) {
      x - mean(
        x[which(datGDS$IDENTIFIER == "blank")],
        na.rm = TRUE
      )
    }
  )

  datGDS$IDENTIFIER = dictReplace(
    datGDS$IDENTIFIER,
    chrFeatures[c(
      "Standard_Gene_Name",
      "Feature_Name"
    )]
  )

  datGDS = datGDS[
    datGDS$IDENTIFIER %in% chrFeatures$Feature_Name,
    ,
    drop = FALSE
  ]

  rownames(datGDS) = make.unique(
    as.character(
      datGDS$IDENTIFIER
    )
  )

  datGDS
}

########################
# LOAD FEATURE METADATA #
########################

load(
  paste0(
    raw,
    "chrFeatures.rda"
  )
)

########################
# BREM 2002 AND 2005   #
########################

brem_files = list(
  Brem_2002_GDS91 = "GDS91.soft",
  Brem_2002_GDS92 = "GDS92.soft",
  Brem_2005_GDS1115 = "GDS1115.soft",
  Brem_2005_GDS1116 = "GDS1116.soft"
)

for (dataset in names(brem_files)) {
  dat = prepare_brem_gds(
    brem_files[[dataset]]
  )

  audit_rmIterNA(
    dat,
    dataset = dataset,
    row_unit = "ORF/probe row",
    col_unit = "array"
  )

  rm(dat)
  gc()
}

########################
# SMITH 2008           #
########################

Smith_2008_pheno = readxl::read_excel(
  paste0(
    raw,
    "pbio.0060083.sd002.xls"
  )
) %>%
  as.data.frame()

rownames(Smith_2008_pheno) = sapply(
  Smith_2008_pheno$NAME,
  function(x) strsplit(
    x,
    "\\|\\|"
  )
) %>%
  (function(x) sapply(
    x,
    function(y) y[4]
  )) %>%
  (function(x) sapply(
    x,
    function(y) gsub(
      " ",
      "",
      y
    )
  ))

Smith_2008_pheno = Smith_2008_pheno[, -c(1:2)]
nuCol = rownames(Smith_2008_pheno)
nuRow = colnames(Smith_2008_pheno)

Smith_2008_pheno = apply(
  Smith_2008_pheno,
  c(1, 2),
  as.numeric
) %>%
  as.matrix() %>%
  t()

colnames(Smith_2008_pheno) = nuCol
rownames(Smith_2008_pheno) = nuRow

SmithMeta = sapply(
  nuRow,
  function(x) strsplit(
    x,
    " "
  )
) %>%
  (function(x) do.call(
    rbind,
    x
  )) %>%
  `colnames<-`(c(
    "id",
    "segregant",
    "condition"
  )) %>%
  as.data.frame()

SmithEth_2008_pheno = Smith_2008_pheno[
  which(SmithMeta$condition == "ethanol"),
  ,
  drop = FALSE
]
rownames(SmithEth_2008_pheno) = SmithMeta$segregant[
  which(SmithMeta$condition == "ethanol")
]

SmithGlc_2008_pheno = Smith_2008_pheno[
  which(SmithMeta$condition == "glucose"),
  ,
  drop = FALSE
]
rownames(SmithGlc_2008_pheno) = SmithMeta$segregant[
  which(SmithMeta$condition == "glucose")
]

audit_rmIterNA(
  SmithEth_2008_pheno,
  dataset = "SmithEth_2008",
  row_unit = "segregant",
  col_unit = "gene"
)

audit_rmIterNA(
  SmithGlc_2008_pheno,
  dataset = "SmithGlc_2008",
  row_unit = "segregant",
  col_unit = "gene"
)

rm(
  Smith_2008_pheno,
  SmithMeta,
  SmithEth_2008_pheno,
  SmithGlc_2008_pheno,
  nuCol,
  nuRow
)

gc()

########################
# KHAN 2009            #
########################

Khan_2009_pheno = read.delim(
  paste0(
    raw,
    "Khan_2009_protein.tsv"
  ),
  sep = "\t",
  row.names = 1,
  check.names = FALSE
) %>%
  as.matrix()

audit_rmIterNA(
  Khan_2009_pheno,
  dataset = "Khan_2009",
  row_unit = "gene",
  col_unit = "segregant"
)

rm(Khan_2009_pheno)
gc()

########################
# FOSS 2011            #
########################

Foss_2011_pheno = readxl::read_excel(
  paste0(
    raw,
    "foss_11_protein_measures.xls"
  )
) %>%
  as.data.frame()

rownames(Foss_2011_pheno) =
  Foss_2011_pheno$`systematic name`

Foss_2011_pheno = Foss_2011_pheno[, -c(1:3)] %>%
  mutate_all(as.numeric) %>%
  as.matrix()

colnames(Foss_2011_pheno) = colnames(Foss_2011_pheno) %>%
  (function(x) sapply(
    x,
    function(y) gsub(
      "_",
      ".",
      y
    )
  ))

audit_rmIterNA(
  Foss_2011_pheno,
  dataset = "Foss_2011",
  row_unit = "gene",
  col_unit = "segregant"
)

rm(Foss_2011_pheno)
gc()

########################
# PICOTTI 2013         #
########################

Picotti_2013_pheno = read.csv(
  paste0(
    raw,
    "Picotti_2013_SRM_dataset.csv"
  ),
  check.names = FALSE
)

picotti_numeric = Picotti_2013_pheno[, -1, drop = FALSE] %>%
  as.matrix()

storage.mode(picotti_numeric) = "numeric"

audit_threshold_only(
  picotti_numeric,
  dataset = "Picotti_2013",
  row_unit = "protein measurement",
  col_unit = "strain measurement"
)

rm(
  Picotti_2013_pheno,
  picotti_numeric
)
gc()

########################
# ALBERT 2018          #
########################

load(
  paste0(
    raw,
    "log2_t.tpm.matrix.RData"
  )
)

audit_threshold_only(
  t.tpm.matrix,
  dataset = "Albert_2018",
  row_unit = "segregant",
  col_unit = "gene"
)

rm(t.tpm.matrix)
gc()

########################
# GROSSBACH 2022 RNA   #
########################

grossbach_rna = readxl::read_excel(
  paste0(
    raw,
    "msb202110712-sup-0003-datasetev1.xlsx"
  ),
  sheet = "transcriptome"
) %>%
  as.data.frame()

grossbach_rna_genes = grossbach_rna$gene
grossbach_rna = t(
  grossbach_rna[, -1]
)
colnames(grossbach_rna) = grossbach_rna_genes

audit_threshold_only(
  grossbach_rna,
  dataset = "GrossbachRNA_2022",
  row_unit = "culture measurement",
  col_unit = "gene"
)

rm(
  grossbach_rna,
  grossbach_rna_genes
)
gc()

########################
# GROSSBACH 2022 PRT   #
########################

grossbach_prt = readxl::read_excel(
  paste0(
    raw,
    "msb202110712-sup-0003-datasetev1.xlsx"
  ),
  sheet = "proteome",
  na = "NA"
) %>%
  as.data.frame()

grossbach_prt_genes = grossbach_prt$gene
grossbach_prt = t(
  grossbach_prt[, -1]
)
colnames(grossbach_prt) = grossbach_prt_genes

grossbach_prt_filtered = audit_rmIterNA(
  grossbach_prt,
  dataset = "GrossbachPRT_2022",
  row_unit = "culture measurement",
  col_unit = "gene"
)

Grossbach_2022_meta = readxl::read_excel(
  paste0(
    raw,
    "msb202110712-sup-0005-datasetev3.xlsx"
  ),
  sheet = "metadata"
) %>%
  as.data.frame()

Grossbach_2022_meta$strain = Grossbach_2022_meta$strain %>%
  unifyChain(c(
    "BY",
    "RM"
  ))

removed_grossbach_prt_cultures = setdiff(
  rownames(grossbach_prt),
  rownames(grossbach_prt_filtered)
)

removed_grossbach_prt_strains = Grossbach_2022_meta$strain[
  match(
    removed_grossbach_prt_cultures,
    Grossbach_2022_meta$culture
  )
]

removed_grossbach_prt_strains = unique(
  removed_grossbach_prt_strains[
    !is.na(removed_grossbach_prt_strains) &
      !(removed_grossbach_prt_strains %in% c("BY", "RM"))
  ]
)

grossbach_prt_strain_summary = data.frame(
  dataset = "GrossbachPRT_2022",
  n_culture_measurements_removed_by_missingness =
    length(removed_grossbach_prt_cultures),
  n_unique_nonparent_strains_touched =
    length(removed_grossbach_prt_strains),
  stringsAsFactors = FALSE
)

rm(
  grossbach_prt,
  grossbach_prt_genes,
  grossbach_prt_filtered,
  Grossbach_2022_meta
)
gc()

############################
# FINAL MAPPING DIMENSIONS #
############################

mapping_studies = c(
  "Brem_2002",
  "Brem_2005",
  "SmithEth_2008",
  "SmithGlc_2008",
  "Khan_2009",
  "Foss_2011",
  "Picotti_2013",
  "Albert_2018",
  "GrossbachRNA_2022",
  "GrossbachPRT_2022"
)

mapping_dimensions = lapply(
  mapping_studies,
  function(study_name) {
    env = new.env(
      parent = emptyenv()
    )

    load(
      paste0(
        imed,
        study_name,
        ".rda"
      ),
      envir = env
    )

    study = env[[study_name]]

    data.frame(
      dataset = study_name,
      final_n_segregants = nrow(study$pheno),
      final_n_genes = ncol(study$pheno),
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

################
# WRITE OUTPUT #
################

missingness_summary = bind_rows(
  summary_rows
)

removed_rows =
  if (length(removed_row_rows) > 0L) {
    bind_rows(removed_row_rows)
  } else {
    data.frame(
      dataset = character(0),
      margin = character(0),
      unit = character(0),
      id = character(0)
    )
  }

removed_columns =
  if (length(removed_col_rows) > 0L) {
    bind_rows(removed_col_rows)
  } else {
    data.frame(
      dataset = character(0),
      margin = character(0),
      unit = character(0),
      id = character(0)
    )
  }

threshold_checks = bind_rows(
  threshold_check_rows
)

write.table(
  missingness_summary,
  file = paste0(
    review_asys,
    "phenotype_missingness_filter_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  removed_rows,
  file = paste0(
    review_asys,
    "phenotype_missingness_removed_rows.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  removed_columns,
  file = paste0(
    review_asys,
    "phenotype_missingness_removed_columns.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  threshold_checks,
  file = paste0(
    review_asys,
    "phenotype_missingness_threshold_checks.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  grossbach_prt_strain_summary,
  file = paste0(
    review_asys,
    "GrossbachPRT_missingness_strain_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  mapping_dimensions,
  file = paste0(
    review_asys,
    "final_mapping_dimensions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("")
message("Phenotype missingness audit complete.")
message("")
print(missingness_summary)
message("")
message("Datasets without an explicit rmIterNA call:")
print(threshold_checks)
message("")
message("Grossbach protein strain-level consequence:")
print(grossbach_prt_strain_summary)
message("")
message("Final mapping dimensions:")
print(mapping_dimensions)
message("")
message("Outputs written to: ", review_asys)
