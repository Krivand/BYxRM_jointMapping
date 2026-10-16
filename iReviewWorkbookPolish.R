########################################
########################################
## GENETICS REVIEW WORKBOOK POLISH    ##
## PLAIN FORMATTING AND NORMAL LABELS ##
########################################
########################################

source("/PHShome/kv1026/p05/remap_functions.R")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package openxlsx is required to polish reviewer workbooks.")
}

review_asys = paste0(asys, "review/")
final_asys = paste0(review_asys, "final/")
polish_asys = paste0(review_asys, ".polished_workbooks/")

final_files = c(
  "replicate_filtering.xlsx",
  "supported_gene_characteristics.xlsx",
  "interaction_GO.xlsx",
  "power_analysis.xlsx",
  "segregant_dependence.xlsx"
)

sheet_plan = list(
  replicate_filtering.xlsx = c(
    Summary = "Summary",
    Removed_items = "Removed data"
  ),
  supported_gene_characteristics.xlsx = c(
    Summary = "Summary",
    Gene_level = "Gene data"
  ),
  interaction_GO.xlsx = c(
    Summary = "Summary",
    All_strict_64 = "All 64 interactions",
    Discordant_37 = "Discordant 37",
    Protein_specific_14 = "Protein specific 14",
    mRNA_specific_9 = "mRNA specific 9",
    Same_direction_4 = "Same direction 4",
    Significant_all = "Significant terms",
    Ribosome_contrast = "Ribosome comparison",
    Gene_sets = "Gene sets"
  ),
  power_analysis.xlsx = c(
    Study_summary = "Study summary",
    Validation_summary = "Sensitivity summary",
    Power_grid = "Power by study",
    Albert_crossstudy = "Albert to other studies",
    Figure3_targets = "Figure 3 by study",
    Figure3_fit = "Figure 3 summary",
    Albert_benchmark = "Albert benchmark",
    Analytic_engine_grid = "Power method comparison",
    Direct_simulation = "Simulation results"
  ),
  segregant_dependence.xlsx = c(
    Summary = "Summary",
    Locus_level = "Locus results",
    Panel_overlap = "Panel overlap"
  )
)

purpose = list(
  replicate_filtering.xlsx = c(
    "Summary" = "Replicate and missingness counts",
    "Removed data" = "Genes and segregants removed by the iterative 25% missingness filter"
  ),
  supported_gene_characteristics.xlsx = c(
    "Summary" = "Summary of the 385-gene representativeness comparisons",
    "Gene data" = "Albert et al. 2018 gene characteristics for the 5,720-gene comparison set"
  ),
  interaction_GO.xlsx = c(
    "Summary" = "GO terms highlighted in the response and the ribosomal-subunit class comparison",
    "All 64 interactions" = "All tested GO terms for the 64 strict interaction-QTL genes",
    "Discordant 37" = "All tested GO terms for the 37 discordant interaction-QTL genes",
    "Protein specific 14" = "All tested GO terms for the 14 protein-specific interaction-QTL genes",
    "mRNA specific 9" = "All tested GO terms for the 9 mRNA-specific interaction-QTL genes",
    "Same direction 4" = "All tested GO terms for the 4 same-direction interaction-QTL genes",
    "Significant terms" = "All GO terms passing FDR < 0.05 across the gene-set and ontology analyses",
    "Ribosome comparison" = "Discordant versus other interaction-QTL genes for six ribosome-related GO terms",
    "Gene sets" = "Genes included in the full strict-interaction and subclass GO analyses"
  ),
  power_analysis.xlsx = c(
    "Study summary" = "Sample size, detectable PVE, observed rediscovery, and predicted rediscovery",
    "Sensitivity summary" = "Attenuation, published-benchmark, and simulation summary measures",
    "Power by study" = "Full 10-study by four-LOD detectable-effect results",
    "Albert to other studies" = "Albert-reference predictions for each smaller study and sensitivity scales",
    "Figure 3 by study" = "Observed Figure 3 rediscovery and power predictions for each target study",
    "Figure 3 summary" = "Across-study summaries at each stipulated LOD",
    "Albert benchmark" = "Reproduction of the Albert et al. 2018 public power calculation",
    "Power method comparison" = "Published-code and current analytic power across 0.5-10% PVE",
    "Simulation results" = "Direct phenotype simulations across 27 representative settings"
  ),
  segregant_dependence.xlsx = c(
    "Summary" = "Shared-segregant sensitivity results",
    "Locus results" = "All 928 targeted original, segregant-clustered, and crossed-random-effect results",
    "Panel overlap" = "Shared-segregant counts for all 21 pairs of joint-mapping datasets"
  )
)

clean_visible_terms = function(dat) {
  names(dat) = gsub("contrast", "comparison", names(dat), ignore.case = TRUE)
  names(dat) = gsub("engine", "method", names(dat), ignore.case = TRUE)

  for (j in seq_along(dat)) {
    if (is.factor(dat[[j]])) dat[[j]] = as.character(dat[[j]])
    if (is.character(dat[[j]])) {
      dat[[j]] = gsub("contrast", "comparison", dat[[j]], ignore.case = TRUE)
      dat[[j]] = gsub("engine", "method", dat[[j]], ignore.case = TRUE)
    }
  }

  dat
}

sheet_widths = function(dat) {
  vapply(
    seq_along(dat),
    function(j) {
      txt = as.character(dat[[j]])
      txt[is.na(txt)] = "NA"
      observed = if (length(txt) > 0L) max(nchar(txt), na.rm = TRUE) else 0L
      max(10, min(45, max(nchar(names(dat)[[j]]), observed) + 2))
    },
    numeric(1)
  )
}

write_plain_workbook = function(sheets, filename) {
  wb = openxlsx::createWorkbook()
  header_style = openxlsx::createStyle(
    textDecoration = "bold",
    border = "bottom",
    borderStyle = "thin",
    halign = "left",
    valign = "top"
  )
  body_style = openxlsx::createStyle(
    wrapText = TRUE,
    valign = "top"
  )

  for (sheet_name in names(sheets)) {
    dat = clean_visible_terms(sheets[[sheet_name]])
    if (ncol(dat) < 1L) stop("Sheet has no columns: ", filename, " / ", sheet_name)

    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = dat,
      startRow = 1,
      startCol = 1,
      headerStyle = header_style,
      withFilter = nrow(dat) > 0L
    )
    openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
    openxlsx::setColWidths(
      wb,
      sheet = sheet_name,
      cols = seq_len(ncol(dat)),
      widths = sheet_widths(dat)
    )
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = body_style,
      rows = seq_len(max(1L, nrow(dat) + 1L)),
      cols = seq_len(ncol(dat)),
      gridExpand = TRUE,
      stack = TRUE
    )
  }

  openxlsx::saveWorkbook(
    wb,
    file = paste0(polish_asys, filename),
    overwrite = TRUE
  )
}

present = sort(list.files(final_asys, pattern = "\\.xlsx$"))
if (!identical(present, sort(final_files))) {
  stop(
    "Expected the five reviewer workbooks before polishing; found: ",
    paste(present, collapse = ", ")
  )
}

if (dir.exists(polish_asys)) unlink(polish_asys, recursive = TRUE, force = TRUE)
ensure_directory(polish_asys)

polished_data = list()

for (filename in final_files) {
  input_path = paste0(final_asys, filename)
  input_sheets = openxlsx::getSheetNames(input_path)
  plan = sheet_plan[[filename]]

  missing = setdiff(names(plan), input_sheets)
  if (length(missing) > 0L) {
    stop(
      "Missing expected sheet(s) in ", filename, ": ",
      paste(missing, collapse = ", ")
    )
  }

  sheets = lapply(
    names(plan),
    function(old_name) {
      openxlsx::read.xlsx(
        input_path,
        sheet = old_name,
        check.names = FALSE,
        na.strings = c("NA")
      )
    }
  )
  names(sheets) = unname(plan)

  polished_data[[filename]] = sheets
  write_plain_workbook(sheets, filename)
}

expected_sheet_names = lapply(sheet_plan, unname)
for (filename in final_files) {
  observed = openxlsx::getSheetNames(paste0(polish_asys, filename))
  expected = expected_sheet_names[[filename]]
  if (!identical(observed, expected)) {
    stop("Polished workbook sheet mismatch: ", filename)
  }
  if (any(grepl("validation|engine|contrast", observed, ignore.case = TRUE))) {
    stop("Unwanted technical wording remains in sheet names: ", filename)
  }
}

if (dir.exists(final_asys)) unlink(final_asys, recursive = TRUE, force = TRUE)
if (!file.rename(polish_asys, final_asys)) {
  stop("Could not promote polished reviewer workbooks to final/.")
}

manifest_rows = list()
manifest_i = 1L

for (filename in final_files) {
  for (sheet_name in names(polished_data[[filename]])) {
    dat = clean_visible_terms(polished_data[[filename]][[sheet_name]])
    manifest_rows[[manifest_i]] = data.frame(
      workbook = filename,
      sheet = sheet_name,
      rows = nrow(dat),
      columns = ncol(dat),
      purpose = unname(purpose[[filename]][[sheet_name]]),
      stringsAsFactors = FALSE
    )
    manifest_i = manifest_i + 1L
  }
}

manifest = do.call(rbind, manifest_rows)
rownames(manifest) = NULL

write.table(
  manifest,
  file = paste0(review_asys, "review_output_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

if (nrow(manifest) != 25L) {
  stop("Expected 25 final workbook sheets after removing the QA-only sheet.")
}

cat("\nFINAL REVIEWER WORKBOOKS POLISHED\n")
cat("Workbooks: 5\n")
cat("Sheets: 25\n")
cat("Formatting: plain white sheets with bold headers and filters\n")
cat("QA-only validation sheet removed\n")
cat("Technical labels 'contrast' and 'engine' removed from the visible workbooks\n")
cat("Directory: ", final_asys, "\n", sep = "")
