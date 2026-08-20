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

review_asys = paste0(asys, "review/missingness/")
ensure_directory(review_asys)

summary_rows = list()
removed_rows = list()

###################
# AUDIT FUNCTIONS #
###################

record_removed = function(dataset, margin, unit, ids) {
  if (length(ids) == 0L) {
    return(invisible(NULL))
  }

  removed_rows[[length(removed_rows) + 1L]] <<- data.frame(
    dataset = dataset,
    margin = margin,
    unit = unit,
    id = as.character(ids),
    stringsAsFactors = FALSE
  )

  invisible(NULL)
}

audit_filter = function(
  x,
  dataset,
  row_unit,
  col_unit,
  apply_filter = TRUE
) {
  if (is.null(dim(x)) || length(dim(x)) != 2L) {
    stop("Expected a two-dimensional object for ", dataset, ".")
  }

  before_rows = rownames(x)
  before_cols = colnames(x)
  before_dim = dim(x)

  if (isTRUE(apply_filter)) {
    filtered = rmIterNA(x, t = threshNA)
  } else {
    filtered = x
  }

  after_rows = rownames(filtered)
  after_cols = colnames(filtered)
  after_dim = dim(filtered)

  removed_row_ids = if (!is.null(before_rows) && !is.null(after_rows)) {
    before_rows[!(before_rows %in% after_rows)]
  } else {
    character(0)
  }

  removed_col_ids = if (!is.null(before_cols) && !is.null(after_cols)) {
    before_cols[!(before_cols %in% after_cols)]
  } else {
    character(0)
  }

  row_missing_fraction = if (is.data.frame(x)) {
    numeric_cols = vapply(x, is.numeric, logical(1))
    if (sum(numeric_cols) > 0L) {
      rowMeans(is.na(x[, numeric_cols, drop = FALSE]))
    } else {
      rep(NA_real_, nrow(x))
    }
  } else {
    rowMeans(is.na(x))
  }

  col_missing_fraction = if (is.data.frame(x)) {
    numeric_cols = vapply(x, is.numeric, logical(1))
    out = rep(NA_real_, ncol(x))
    out[numeric_cols] = colMeans(is.na(x[, numeric_cols, drop = FALSE]))
    out
  } else {
    colMeans(is.na(x))
  }

  summary_rows[[length(summary_rows) + 1L]] <<- data.frame(
    dataset = dataset,
    filter_applied = isTRUE(apply_filter),
    row_unit = row_unit,
    col_unit = col_unit,
    n_rows_before = before_dim[1],
    n_rows_after = after_dim[1],
    n_rows_removed = before_dim[1] - after_dim[1],
    n_cols_before = before_dim[2],
    n_cols_after = after_dim[2],
    n_cols_removed = before_dim[2] - after_dim[2],
    n_rows_at_or_above_25pct_before = sum(
      row_missing_fraction >= threshNA,
      na.rm = TRUE
    ),
    n_cols_at_or_above_25pct_before = sum(
      col_missing_fraction >= threshNA,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )

  record_removed(
    dataset,
    "row",
    row_unit,
    removed_row_ids
  )

  record_removed(
    dataset,
    "column",
    col_unit,
    removed_col_ids
  )

  list(
    filtered = filtered,
    removed_rows = removed_row_ids,
    removed_cols = removed_col_ids
  )
}

read_gds_soft_table = function(path) {
  lines = readLines(path, warn = FALSE)

  begin = grep(
    "^!dataset_table_begin",
    lines,
    ignore.case = TRUE
  )
  end = grep(
    "^!dataset_table_end",
    lines,
    ignore.case = TRUE
  )

  if (length(begin) != 1L || length(end) != 1L || end <= begin) {
    stop("Could not identify a single dataset table in ", path, ".")
  }

  table_lines = lines[(begin + 1L):(end - 1L)]

  dat = read.delim(
    text = paste(table_lines, collapse = "\n"),
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "null", "NULL")
  )

  if (!("IDENTIFIER" %in% colnames(dat))) {
    stop("IDENTIFIER column not found in ", path, ".")
  }

  if (ncol(dat) < 3L) {
    stop("Unexpected SOFT table dimensions in ", path, ".")
  }

  for (j in 3:ncol(dat)) {
    dat[[j]] = suppressWarnings(as.numeric(dat[[j]]))
  }

  dat
}

prepare_brem_gds = function(filename) {
  datGDS = read_gds_soft_table(
    paste0(raw, filename)
  )

  blank_rows = which(datGDS$IDENTIFIER == "blank")

  if (length(blank_rows) > 0L) {
    datGDS[, 3:ncol(datGDS)] = apply(
      as.matrix(datGDS[, -c(1:2), drop = FALSE]),
      2,
      function(x) {
        x - mean(x[blank_rows], na.rm = TRUE)
      }
    )
  }

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
    as.character(datGDS$IDENTIFIER)
  )

  datGDS
}

########################
# LOAD FEATURE METADATA #
########################

load(paste0(raw, "chrFeatures.rda"))

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
  message("Auditing ", dataset, "...")

  dat = prepare_brem_gds(
    brem_files[[dataset]]
  )

  audit_filter(
    dat,
    dataset = dataset,
    row_unit = "ORF/probe row",
    col_unit = "array",
    apply_filter = TRUE
  )

  rm(dat)
  gc()
}

########################
# SMITH 2008           #
########################

message("Auditing Smith 2008...")

Smith_2008_pheno = readxl::read_excel(
  paste0(raw, "pbio.0060083.sd002.xls")
) %>%
  as.data.frame()

rownames(Smith_2008_pheno) = sapply(
  Smith_2008_pheno$NAME,
  function(x) strsplit(x, "\\|\\|")
) %>%
  (function(x) sapply(x, function(y) y[4])) %>%
  (function(x) sapply(x, function(y) gsub(" ", "", y)))

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
  function(x) strsplit(x, " ")
) %>%
  (function(x) do.call(rbind, x)) %>%
  `colnames<-`(c("id", "segregant", "condition")) %>%
  as.data.frame()

SmithEth_2008_pheno = Smith_2008_pheno[
  SmithMeta$condition == "ethanol",
  ,
  drop = FALSE
]
rownames(SmithEth_2008_pheno) = SmithMeta$segregant[
  SmithMeta$condition == "ethanol"
]

SmithGlc_2008_pheno = Smith_2008_pheno[
  SmithMeta$condition == "glucose",
  ,
  drop = FALSE
]
rownames(SmithGlc_2008_pheno) = SmithMeta$segregant[
  SmithMeta$condition == "glucose"
]

audit_filter(
  SmithEth_2008_pheno,
  "SmithEth_2008",
  "segregant",
  "gene"
)

audit_filter(
  SmithGlc_2008_pheno,
  "SmithGlc_2008",
  "segregant",
  "gene"
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

message("Auditing Khan 2009...")

Khan_2009_pheno = read.delim(
  paste0(raw, "Khan_2009_protein.tsv"),
  sep = "\t",
  row.names = 1,
  check.names = FALSE
) %>%
  as.matrix()

audit_filter(
  Khan_2009_pheno,
  "Khan_2009",
  "gene",
  "segregant"
)

rm(Khan_2009_pheno)
gc()

########################
# FOSS 2011            #
########################

message("Auditing Foss 2011...")

Foss_2011_pheno = readxl::read_excel(
  paste0(raw, "foss_11_protein_measures.xls")
) %>%
  as.data.frame()

rownames(Foss_2011_pheno) =
  Foss_2011_pheno$`systematic name`

Foss_2011_pheno = Foss_2011_pheno[, -c(1:3)] %>%
  mutate_all(as.numeric) %>%
  as.matrix()

colnames(Foss_2011_pheno) = gsub(
  "_",
  ".",
  colnames(Foss_2011_pheno)
)

audit_filter(
  Foss_2011_pheno,
  "Foss_2011",
  "gene",
  "segregant"
)

rm(Foss_2011_pheno)
gc()

########################
# PICOTTI 2013         #
########################

message("Checking Picotti 2013...")

Picotti_2013_pheno = read.csv(
  paste0(raw, "Picotti_2013_SRM_dataset.csv"),
  check.names = FALSE
)

picotti_numeric = as.matrix(
  Picotti_2013_pheno[, -1, drop = FALSE]
)
storage.mode(picotti_numeric) = "numeric"

audit_filter(
  picotti_numeric,
  "Picotti_2013",
  "protein measurement",
  "strain measurement",
  apply_filter = FALSE
)

rm(Picotti_2013_pheno, picotti_numeric)
gc()

########################
# ALBERT 2018          #
########################

message("Checking Albert 2018...")

load(paste0(raw, "log2_t.tpm.matrix.RData"))

audit_filter(
  t.tpm.matrix,
  "Albert_2018",
  "segregant",
  "gene",
  apply_filter = FALSE
)

rm(t.tpm.matrix)
gc()

########################
# GROSSBACH 2022 RNA   #
########################

message("Checking Grossbach 2022 RNA...")

grossbach_rna = readxl::read_excel(
  paste0(raw, "msb202110712-sup-0003-datasetev1.xlsx"),
  sheet = "transcriptome"
) %>%
  as.data.frame()

grossbach_rna_genes = grossbach_rna$gene
grossbach_rna = t(grossbach_rna[, -1])
colnames(grossbach_rna) = grossbach_rna_genes

audit_filter(
  grossbach_rna,
  "GrossbachRNA_2022",
  "culture measurement",
  "gene",
  apply_filter = FALSE
)

rm(grossbach_rna, grossbach_rna_genes)
gc()

########################
# GROSSBACH 2022 PRT   #
########################

message("Auditing Grossbach 2022 protein...")

grossbach_prt = readxl::read_excel(
  paste0(raw, "msb202110712-sup-0003-datasetev1.xlsx"),
  sheet = "proteome",
  na = "NA"
) %>%
  as.data.frame()

grossbach_prt_genes = grossbach_prt$gene
grossbach_prt = t(grossbach_prt[, -1])
colnames(grossbach_prt) = grossbach_prt_genes

grossbach_prt_audit = audit_filter(
  grossbach_prt,
  "GrossbachPRT_2022",
  "culture measurement",
  "gene",
  apply_filter = TRUE
)

Grossbach_2022_meta = readxl::read_excel(
  paste0(raw, "msb202110712-sup-0005-datasetev3.xlsx"),
  sheet = "metadata"
) %>%
  as.data.frame()

Grossbach_2022_meta$strain = unifyChain(
  Grossbach_2022_meta$strain,
  c("BY", "RM")
)

before_strains = Grossbach_2022_meta$strain[
  match(rownames(grossbach_prt), Grossbach_2022_meta$culture)
]
after_strains = Grossbach_2022_meta$strain[
  match(
    rownames(grossbach_prt_audit$filtered),
    Grossbach_2022_meta$culture
  )
]

before_strains = unique(
  before_strains[
    !is.na(before_strains) &
      !(before_strains %in% c("BY", "RM"))
  ]
)
after_strains = unique(
  after_strains[
    !is.na(after_strains) &
      !(after_strains %in% c("BY", "RM"))
  ]
)

grossbach_prt_strain_loss = data.frame(
  dataset = "GrossbachPRT_2022",
  n_nonparent_strains_before_missingness = length(before_strains),
  n_nonparent_strains_after_missingness = length(after_strains),
  n_nonparent_strains_lost_to_missingness = length(
    setdiff(before_strains, after_strains)
  ),
  stringsAsFactors = FALSE
)

rm(
  grossbach_prt,
  grossbach_prt_genes,
  grossbach_prt_audit,
  Grossbach_2022_meta,
  before_strains,
  after_strains
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
    env = new.env(parent = emptyenv())
    load(
      paste0(imed, study_name, ".rda"),
      envir = env
    )
    study = env[[study_name]]

    data.frame(
      dataset = study_name,
      final_n_segregants = nrow(study$pheno),
      final_n_features = ncol(study$pheno),
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

################
# WRITE OUTPUT #
################

missingness_summary = bind_rows(summary_rows)

removed_ids = if (length(removed_rows) > 0L) {
  bind_rows(removed_rows)
} else {
  data.frame(
    dataset = character(0),
    margin = character(0),
    unit = character(0),
    id = character(0),
    stringsAsFactors = FALSE
  )
}

write.table(
  missingness_summary,
  file = paste0(
    review_asys,
    "missingness_filter_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  removed_ids,
  file = paste0(
    review_asys,
    "missingness_removed_ids.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  grossbach_prt_strain_loss,
  file = paste0(
    review_asys,
    "GrossbachPRT_strain_loss.tsv"
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
message("Missingness audit complete.")
message("")
message("Direct effects of the 25% missingness filter:")
print(missingness_summary)
message("")
message("Grossbach protein strain-level loss:")
print(grossbach_prt_strain_loss)
message("")
message("Final mapping dimensions after all preprocessing/harmonization:")
print(mapping_dimensions)
message("")
message("Outputs written to: ", review_asys)
