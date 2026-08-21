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

load(paste0(raw, "chrFeatures.rda"))

summary_rows = list()

####################
# HELPER FUNCTIONS #
####################

add_summary = function(
  dataset,
  feature_unit,
  n_features_before,
  n_features_after,
  sample_unit,
  n_samples_before,
  n_samples_after
) {

  summary_rows[[length(summary_rows) + 1L]] <<- data.frame(
    dataset = dataset,
    feature_unit = feature_unit,
    n_features_before = n_features_before,
    n_features_after = n_features_after,
    n_features_removed = n_features_before - n_features_after,
    sample_unit = sample_unit,
    n_samples_before = n_samples_before,
    n_samples_after = n_samples_after,
    n_samples_removed = n_samples_before - n_samples_after,
    stringsAsFactors = FALSE
  )
}

audit_matrix = function(
  x,
  dataset,
  feature_margin,
  feature_unit,
  sample_unit,
  apply_filter = TRUE
) {

  if (is.null(dim(x)) || length(dim(x)) != 2L) {
    stop("Expected a two-dimensional object for ", dataset, ".")
  }

  before = dim(x)

  filtered = if (isTRUE(apply_filter)) {
    rmIterNA(x, t = threshNA)
  } else {
    x
  }

  after = dim(filtered)

  sample_margin = if (feature_margin == 1L) 2L else 1L

  add_summary(
    dataset = dataset,
    feature_unit = feature_unit,
    n_features_before = before[feature_margin],
    n_features_after = after[feature_margin],
    sample_unit = sample_unit,
    n_samples_before = before[sample_margin],
    n_samples_after = after[sample_margin]
  )

  invisible(NULL)
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

  if (
    length(begin) != 1L ||
    length(end) != 1L ||
    end <= begin
  ) {
    stop("Could not identify a single dataset table in ", path, ".")
  }

  dat = read.delim(
    text = paste(
      lines[(begin + 1L):(end - 1L)],
      collapse = "\n"
    ),
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "null", "NULL")
  )

  for (j in 3:ncol(dat)) {
    dat[[j]] = suppressWarnings(as.numeric(dat[[j]]))
  }

  dat
}

prepare_brem_gds = function(filename) {

  dat = read_gds_soft_table(
    paste0(raw, filename)
  )

  blank_rows = which(dat$IDENTIFIER == "blank")

  if (length(blank_rows) > 0L) {
    dat[, 3:ncol(dat)] = apply(
      as.matrix(dat[, -c(1:2), drop = FALSE]),
      2,
      function(x) {
        x - mean(x[blank_rows], na.rm = TRUE)
      }
    )
  }

  dat$IDENTIFIER = dictReplace(
    dat$IDENTIFIER,
    chrFeatures[c(
      "Standard_Gene_Name",
      "Feature_Name"
    )]
  )

  dat[
    dat$IDENTIFIER %in% chrFeatures$Feature_Name,
    ,
    drop = FALSE
  ]
}

brem_channel_sets = function(
  filename,
  pattern,
  old = "-",
  new = "."
) {

  dat = prepare_brem_gds(filename)

  strain_key = getStrainInfo(
    paste0(raw, filename)
  )

  strain_key$strain = strain_key$strain %>%
    (function(x) sapply(
      x,
      function(y) sub(pattern, "", y)
    )) %>%
    unifyChain(c("BY", "RM"))

  strain_key$strain = gsub(
    old,
    new,
    strain_key$strain
  )

  before_genes = unique(as.character(dat$IDENTIFIER))
  before_arrays = colnames(dat)[3:ncol(dat)]
  before_strains = unique(
    strain_key$strain[
      match(before_arrays, strain_key$gsmID)
    ]
  )
  before_strains = before_strains[!is.na(before_strains)]

  filtered = rmIterNA(
    dat,
    t = threshNA
  )

  after_genes = unique(as.character(filtered$IDENTIFIER))
  after_arrays = colnames(filtered)[3:ncol(filtered)]
  after_strains = unique(
    strain_key$strain[
      match(after_arrays, strain_key$gsmID)
    ]
  )
  after_strains = after_strains[!is.na(after_strains)]

  list(
    before_genes = before_genes,
    after_genes = after_genes,
    before_strains = before_strains,
    after_strains = after_strains
  )
}

audit_brem = function(
  study,
  file1,
  pattern1,
  file2,
  pattern2
) {

  ch1 = brem_channel_sets(
    file1,
    pattern1
  )

  ch2 = brem_channel_sets(
    file2,
    pattern2
  )

  genes_before = intersect(
    ch1$before_genes,
    ch2$before_genes
  )

  genes_after = intersect(
    ch1$after_genes,
    ch2$after_genes
  )

  strains_before = intersect(
    ch1$before_strains,
    ch2$before_strains
  )

  strains_after = intersect(
    ch1$after_strains,
    ch2$after_strains
  )

  add_summary(
    dataset = study,
    feature_unit = "gene",
    n_features_before = length(genes_before),
    n_features_after = length(genes_after),
    sample_unit = "segregant",
    n_samples_before = length(strains_before),
    n_samples_after = length(strains_after)
  )
}

########################
# BREM 2002 AND 2005   #
########################

audit_brem(
  study = "Brem_2002",
  file1 = "GDS91.soft",
  pattern1 = "BY-",
  file2 = "GDS92.soft",
  pattern2 = "-BY"
)

audit_brem(
  study = "Brem_2005",
  file1 = "GDS1115.soft",
  pattern1 = "Cy5",
  file2 = "GDS1116.soft",
  pattern2 = "Cy3"
)

########################
# SMITH 2008           #
########################

smith = readxl::read_excel(
  paste0(raw, "pbio.0060083.sd002.xls")
) %>%
  as.data.frame()

rownames(smith) = sapply(
  smith$NAME,
  function(x) strsplit(x, "\\|\\|")
) %>%
  (function(x) sapply(x, function(y) y[4])) %>%
  (function(x) sapply(x, function(y) gsub(" ", "", y)))

smith = smith[, -c(1:2)]
genes = rownames(smith)
sample_names = colnames(smith)

smith = apply(
  smith,
  c(1, 2),
  as.numeric
) %>%
  as.matrix() %>%
  t()

colnames(smith) = genes
rownames(smith) = sample_names

smith_meta = sapply(
  sample_names,
  function(x) strsplit(x, " ")
) %>%
  (function(x) do.call(rbind, x)) %>%
  `colnames<-`(c("id", "segregant", "condition")) %>%
  as.data.frame()

smith_eth = smith[
  smith_meta$condition == "ethanol",
  ,
  drop = FALSE
]

smith_glc = smith[
  smith_meta$condition == "glucose",
  ,
  drop = FALSE
]

rownames(smith_eth) = smith_meta$segregant[
  smith_meta$condition == "ethanol"
]

rownames(smith_glc) = smith_meta$segregant[
  smith_meta$condition == "glucose"
]

audit_matrix(
  smith_eth,
  "SmithEth_2008",
  feature_margin = 2L,
  feature_unit = "gene",
  sample_unit = "segregant"
)

audit_matrix(
  smith_glc,
  "SmithGlc_2008",
  feature_margin = 2L,
  feature_unit = "gene",
  sample_unit = "segregant"
)

rm(
  smith,
  smith_meta,
  smith_eth,
  smith_glc,
  genes,
  sample_names
)
gc()

########################
# KHAN 2009            #
########################

khan = read.delim(
  paste0(raw, "Khan_2009_protein.tsv"),
  sep = "\t",
  row.names = 1,
  check.names = FALSE
) %>%
  as.matrix()

audit_matrix(
  khan,
  "Khan_2009",
  feature_margin = 1L,
  feature_unit = "gene",
  sample_unit = "segregant"
)

rm(khan)
gc()

########################
# FOSS 2011            #
########################

foss = readxl::read_excel(
  paste0(raw, "foss_11_protein_measures.xls")
) %>%
  as.data.frame()

rownames(foss) = foss$`systematic name`

foss = foss[, -c(1:3)] %>%
  mutate_all(as.numeric) %>%
  as.matrix()

colnames(foss) = gsub(
  "_",
  ".",
  colnames(foss)
)

audit_matrix(
  foss,
  "Foss_2011",
  feature_margin = 1L,
  feature_unit = "gene",
  sample_unit = "segregant"
)

rm(foss)
gc()

########################
# DATASETS WITH NO LOSS #
########################

picotti = read.csv(
  paste0(raw, "Picotti_2013_SRM_dataset.csv"),
  check.names = FALSE
)

picotti = as.matrix(
  picotti[, -1, drop = FALSE]
)
storage.mode(picotti) = "numeric"

audit_matrix(
  picotti,
  "Picotti_2013",
  feature_margin = 1L,
  feature_unit = "protein measurement",
  sample_unit = "strain measurement",
  apply_filter = FALSE
)

rm(picotti)
gc()

load(paste0(raw, "log2_t.tpm.matrix.RData"))

audit_matrix(
  t.tpm.matrix,
  "Albert_2018",
  feature_margin = 2L,
  feature_unit = "gene",
  sample_unit = "segregant",
  apply_filter = FALSE
)

rm(t.tpm.matrix)
gc()

grossbach_rna = readxl::read_excel(
  paste0(raw, "msb202110712-sup-0003-datasetev1.xlsx"),
  sheet = "transcriptome"
) %>%
  as.data.frame()

grossbach_rna_genes = grossbach_rna$gene
grossbach_rna = t(grossbach_rna[, -1])
colnames(grossbach_rna) = grossbach_rna_genes

audit_matrix(
  grossbach_rna,
  "GrossbachRNA_2022",
  feature_margin = 2L,
  feature_unit = "gene",
  sample_unit = "culture measurement",
  apply_filter = FALSE
)

rm(grossbach_rna, grossbach_rna_genes)
gc()

grossbach_prt = readxl::read_excel(
  paste0(raw, "msb202110712-sup-0003-datasetev1.xlsx"),
  sheet = "proteome",
  na = "NA"
) %>%
  as.data.frame()

grossbach_prt_genes = grossbach_prt$gene
grossbach_prt = t(grossbach_prt[, -1])
colnames(grossbach_prt) = grossbach_prt_genes

audit_matrix(
  grossbach_prt,
  "GrossbachPRT_2022",
  feature_margin = 2L,
  feature_unit = "gene",
  sample_unit = "culture measurement"
)

rm(grossbach_prt, grossbach_prt_genes)
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
# SAVE OUTPUTS #
################

missingness_summary = bind_rows(summary_rows)

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
print(missingness_summary)
message("")
message("Final mapping dimensions:")
print(mapping_dimensions)
message("")
message("Outputs written to: ", review_asys)
