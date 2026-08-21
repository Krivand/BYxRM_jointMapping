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

####################
# HELPER FUNCTIONS #
####################

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

  after_arrays = colnames(filtered)[3:ncol(filtered)]
  after_strains = unique(
    strain_key$strain[
      match(after_arrays, strain_key$gsmID)
    ]
  )
  after_strains = after_strains[!is.na(after_strains)]

  list(
    before_genes = unique(as.character(dat$IDENTIFIER)),
    after_genes = unique(as.character(filtered$IDENTIFIER)),
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

  segregants_before = intersect(
    ch1$before_strains,
    ch2$before_strains
  )

  segregants_after = intersect(
    ch1$after_strains,
    ch2$after_strains
  )

  data.frame(
    Dataset = study,
    `Genes before filter` = length(genes_before),
    `Genes removed` = length(
      setdiff(genes_before, genes_after)
    ),
    `Segregants before filter` = length(segregants_before),
    `Segregants removed` = length(
      setdiff(segregants_before, segregants_after)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

audit_standard = function(
  phenotype,
  study
) {

  genes_before = rownames(phenotype)
  segregants_before = colnames(phenotype)

  filtered = rmIterNA(
    phenotype,
    t = threshNA
  )

  data.frame(
    Dataset = study,
    `Genes before filter` = nrow(phenotype),
    `Genes removed` = nrow(phenotype) - nrow(filtered),
    `Segregants before filter` = ncol(phenotype),
    `Segregants removed` = ncol(phenotype) - ncol(filtered),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

########################
# BREM 2002 AND 2005   #
########################

brem_2002 = audit_brem(
  study = "Brem 2002",
  file1 = "GDS91.soft",
  pattern1 = "BY-",
  file2 = "GDS92.soft",
  pattern2 = "-BY"
)

brem_2005 = audit_brem(
  study = "Brem 2005",
  file1 = "GDS1115.soft",
  pattern1 = "Cy5",
  file2 = "GDS1116.soft",
  pattern2 = "Cy3"
)

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

khan_2009 = audit_standard(
  khan,
  "Khan 2009"
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

foss_2011 = audit_standard(
  foss,
  "Foss 2011"
)

rm(foss)
gc()

########################
# REVIEW-READY SUMMARY #
########################

missingness_summary = rbind(
  brem_2002,
  brem_2005,
  khan_2009,
  foss_2011
)

expected_genes_before = c(
  6167L,
  6167L,
  539L,
  354L
)

expected_genes_removed = c(
  490L,
  263L,
  159L,
  116L
)

expected_segregants_removed = c(
  0L,
  3L,
  2L,
  0L
)

if (
  !identical(
    as.integer(missingness_summary$`Genes before filter`),
    expected_genes_before
  ) ||
  !identical(
    as.integer(missingness_summary$`Genes removed`),
    expected_genes_removed
  ) ||
  !identical(
    as.integer(missingness_summary$`Segregants removed`),
    expected_segregants_removed
  )
) {
  stop("Missingness counts differ from the reviewed results.")
}

################
# SAVE OUTPUT  #
################

output_file = paste0(
  review_asys,
  "missingness_filter_summary.tsv"
)

write.table(
  missingness_summary,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("")
message("Missingness audit complete.")
message("Output: ", output_file)
message("")
print(missingness_summary, row.names = FALSE)
