########################################
########################################
## GENETICS REVIEW BREM MISSINGNESS   ##
## FINAL DYE-MERGED EFFECTS           ##
########################################
########################################

source("/PHShome/kv1026/p05/remap_functions.R")

load(paste0(raw, "chrFeatures.rda"))

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

  datGDS
}

channel_sets = function(
  filename,
  pat,
  old = "-",
  new = "."
) {
  dat = prepare_brem_gds(filename)

  strainKey = getStrainInfo(
    paste0(raw, filename)
  )

  strainKey$strain = strainKey$strain %>%
    (function(x) sapply(
      x,
      function(y) sub(pat, "", y)
    )) %>%
    unifyChain(c("BY", "RM"))

  if (is.character(old) && is.character(new)) {
    strainKey$strain = gsub(
      old,
      new,
      strainKey$strain
    )
  }

  before_genes = unique(
    as.character(dat$IDENTIFIER)
  )

  before_arrays = colnames(dat)[3:ncol(dat)]
  before_strains = unique(
    strainKey$strain[
      match(
        before_arrays,
        strainKey$gsmID
      )
    ]
  )
  before_strains = before_strains[!is.na(before_strains)]

  filtered = rmIterNA(
    dat,
    t = threshNA
  )

  after_genes = unique(
    as.character(filtered$IDENTIFIER)
  )

  after_arrays = colnames(filtered)[3:ncol(filtered)]
  after_strains = unique(
    strainKey$strain[
      match(
        after_arrays,
        strainKey$gsmID
      )
    ]
  )
  after_strains = after_strains[!is.na(after_strains)]

  list(
    before_genes = before_genes,
    after_genes = after_genes,
    before_strains = before_strains,
    after_strains = after_strains,
    removed_arrays = setdiff(
      before_arrays,
      after_arrays
    )
  )
}

audit_brem = function(
  study,
  file1,
  pat1,
  file2,
  pat2
) {
  ch1 = channel_sets(
    file1,
    pat1
  )
  ch2 = channel_sets(
    file2,
    pat2
  )

  merged_genes_before = intersect(
    ch1$before_genes,
    ch2$before_genes
  )
  merged_genes_after = intersect(
    ch1$after_genes,
    ch2$after_genes
  )

  merged_strains_before = intersect(
    ch1$before_strains,
    ch2$before_strains
  )
  merged_strains_after = intersect(
    ch1$after_strains,
    ch2$after_strains
  )

  cat("\n", study, "\n", sep = "")
  cat(
    "Raw arrays removed by channel: ",
    length(ch1$removed_arrays),
    " and ",
    length(ch2$removed_arrays),
    "\n",
    sep = ""
  )
  cat(
    "Dye-merged strains before missingness filter: ",
    length(merged_strains_before),
    "\n",
    sep = ""
  )
  cat(
    "Dye-merged strains after missingness filter: ",
    length(merged_strains_after),
    "\n",
    sep = ""
  )
  cat(
    "Dye-merged strains removed by missingness: ",
    length(setdiff(
      merged_strains_before,
      merged_strains_after
    )),
    "\n",
    sep = ""
  )
  cat(
    "Dye-merged genes before missingness filter: ",
    length(merged_genes_before),
    "\n",
    sep = ""
  )
  cat(
    "Dye-merged genes after missingness filter: ",
    length(merged_genes_after),
    "\n",
    sep = ""
  )
  cat(
    "Dye-merged genes removed by missingness: ",
    length(setdiff(
      merged_genes_before,
      merged_genes_after
    )),
    "\n",
    sep = ""
  )

  removed_strains = setdiff(
    merged_strains_before,
    merged_strains_after
  )

  if (length(removed_strains) > 0L) {
    cat(
      "Removed strain IDs: ",
      paste(
        removed_strains,
        collapse = ", "
      ),
      "\n",
      sep = ""
    )
  }
}

audit_brem(
  study = "Brem 2002",
  file1 = "GDS91.soft",
  pat1 = "BY-",
  file2 = "GDS92.soft",
  pat2 = "-BY"
)

audit_brem(
  study = "Brem 2005",
  file1 = "GDS1115.soft",
  pat1 = "Cy5",
  file2 = "GDS1116.soft",
  pat2 = "Cy3"
)
