########################################
########################################
## GENETICS REVIEW PRE-POWER SUPPORT  ##
## FINAL WORKBOOK CONSISTENCY PASS     ##
########################################
########################################

source("/PHShome/kv1026/p05/remap_functions.R")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package openxlsx is required for reviewer workbook support.")
}

review_asys = paste0(asys, "review/")
final_asys = paste0(review_asys, "final/")
manifest_path = paste0(review_asys, "review_output_manifest.tsv")
go_workbook_path = paste0(final_asys, "interaction_GO.xlsx")

if (!exists("go_env", inherits = TRUE)) {
  stop("Pre-power finalizer must run after iReviewWorkbookBundle.R.")
}
if (!exists("hotspot_summary", envir = go_env, inherits = FALSE)) {
  stop("GO analysis did not generate hotspot_summary.")
}
if (!file.exists(go_workbook_path)) {
  stop("Missing finalized interaction GO workbook: ", go_workbook_path)
}
if (!file.exists(manifest_path)) {
  stop("Missing reviewer output manifest: ", manifest_path)
}

hotspot_summary = get(
  "hotspot_summary",
  envir = go_env,
  inherits = FALSE
)

########################################
# HELPERS                              #
########################################

sheet_widths = function(dat) {
  if (ncol(dat) == 0L) return(numeric(0))

  vapply(
    seq_along(dat),
    function(j) {
      txt = as.character(dat[[j]])
      txt[is.na(txt)] = "NA"
      observed = if (length(txt) > 0L) {
        max(nchar(txt), na.rm = TRUE)
      } else {
        0L
      }
      max(
        10,
        min(
          45,
          max(nchar(names(dat)[[j]]), observed) + 2
        )
      )
    },
    numeric(1)
  )
}

clean_name = function(x) {
  x = gsub("contrast", "comparison", x, ignore.case = TRUE)
  x = gsub("engine", "method", x, ignore.case = TRUE)
  x = gsub("[ .]+", "_", x)
  x = gsub("_+", "_", x)
  x = gsub("^_|_$", "", x)
  x
}

clean_table = function(dat) {
  dat = as.data.frame(
    dat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(dat) = clean_name(names(dat))
  dat
}

########################################
# ADD OBVIOUS MKT1/BUL2 SUPPORT SHEET #
########################################

hotspot_sheet = "MKT1 BUL2 summary"
hotspot_summary = clean_table(hotspot_summary)

wb = openxlsx::loadWorkbook(go_workbook_path)
existing_sheets = openxlsx::sheets(wb)

if (hotspot_sheet %in% existing_sheets) {
  openxlsx::removeWorksheet(wb, hotspot_sheet)
}

openxlsx::addWorksheet(wb, hotspot_sheet)

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

openxlsx::writeData(
  wb,
  sheet = hotspot_sheet,
  x = hotspot_summary,
  startRow = 1,
  startCol = 1,
  headerStyle = header_style,
  keepNA = TRUE,
  na.string = "NA"
)

if (nrow(hotspot_summary) > 0L) {
  openxlsx::addFilter(
    wb,
    sheet = hotspot_sheet,
    row = 1,
    cols = seq_len(ncol(hotspot_summary))
  )
}

openxlsx::freezePane(
  wb,
  sheet = hotspot_sheet,
  firstRow = TRUE
)
openxlsx::setColWidths(
  wb,
  sheet = hotspot_sheet,
  cols = seq_len(ncol(hotspot_summary)),
  widths = sheet_widths(hotspot_summary)
)
openxlsx::addStyle(
  wb,
  sheet = hotspot_sheet,
  style = body_style,
  rows = seq_len(max(1L, nrow(hotspot_summary) + 1L)),
  cols = seq_len(ncol(hotspot_summary)),
  gridExpand = TRUE,
  stack = TRUE
)

workbook_tmp = paste0(go_workbook_path, ".tmp.xlsx")
if (file.exists(workbook_tmp)) unlink(workbook_tmp, force = TRUE)
openxlsx::saveWorkbook(wb, workbook_tmp, overwrite = TRUE)

hotspot_check = openxlsx::read.xlsx(
  workbook_tmp,
  sheet = hotspot_sheet,
  check.names = FALSE
)

expected_hotspot = data.frame(
  hotspot = c("MKT1", "BUL2"),
  strict_interaction_QTLs = c(36L, 11L),
  discordant_QTLs = c(26L, 1L),
  mRNA_specific_QTLs = c(8L, 0L),
  protein_specific_QTLs = c(2L, 9L),
  same_direction_QTLs = c(0L, 1L),
  discordant_ribosome_associated = c(21L, 0L),
  discordant_mRNA_down_protein_up = c(26L, 1L),
  discordant_mRNA_up_protein_down = c(0L, 0L),
  stringsAsFactors = FALSE
)

hotspot_order = match(expected_hotspot$hotspot, hotspot_check$hotspot)
if (any(is.na(hotspot_order))) {
  stop("Final hotspot support sheet is missing MKT1 or BUL2.")
}
hotspot_check = hotspot_check[
  hotspot_order,
  names(expected_hotspot),
  drop = FALSE
]
for (column in setdiff(names(expected_hotspot), "hotspot")) {
  if (!identical(
    as.integer(hotspot_check[[column]]),
    expected_hotspot[[column]]
  )) {
    stop("Unexpected finalized hotspot count in column: ", column)
  }
}

if (!file.rename(workbook_tmp, go_workbook_path)) {
  stop("Could not promote interaction GO workbook with hotspot summary.")
}

########################################
# UPDATE FINAL MANIFEST                #
########################################

manifest = read.delim(
  manifest_path,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

manifest = manifest[
  !(
    manifest$workbook == "interaction_GO.xlsx" &
    manifest$sheet == hotspot_sheet
  ),
  ,
  drop = FALSE
]

manifest$purpose[
  manifest$workbook == "segregant_dependence.xlsx" &
  manifest$sheet == "Summary"
] = "Crossed-segregant random-intercept sensitivity results"
manifest$purpose[
  manifest$workbook == "segregant_dependence.xlsx" &
  manifest$sheet == "Locus results"
] = "All 928 targeted original and crossed-segregant-random-effect results"

hotspot_manifest = data.frame(
  workbook = "interaction_GO.xlsx",
  sheet = hotspot_sheet,
  rows = nrow(hotspot_summary),
  columns = ncol(hotspot_summary),
  purpose = paste(
    "MKT1 and BUL2 strict-interaction counts, interaction classes,",
    "ribosome membership, effect orientation, and class-enrichment tests"
  ),
  stringsAsFactors = FALSE
)

manifest = rbind(manifest, hotspot_manifest)
rownames(manifest) = NULL

if (nrow(manifest) != 26L) {
  stop("Expected exactly 26 final supporting-data sheets.")
}

write.table(
  manifest,
  file = manifest_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

########################################
# VERIFY CONTINUOUS EFFECT DISPLAY     #
########################################

gene_summary = openxlsx::read.xlsx(
  paste0(final_asys, "supported_gene_characteristics.xlsx"),
  sheet = "Summary",
  check.names = FALSE
)

expression_row = gene_summary[
  gene_summary$Characteristic == "Mean log2 TPM",
  ,
  drop = FALSE
]
ppi_row = gene_summary[
  gene_summary$Characteristic == "PPI partners",
  ,
  drop = FALSE
]

if (
  nrow(expression_row) != 1L ||
  nrow(ppi_row) != 1L ||
  expression_row$Effect_size[[1]] != "6.52-fold higher" ||
  ppi_row$Effect_size[[1]] != "2.21-fold higher"
) {
  stop("Continuous representativeness effects are not both shown as fold differences.")
}

########################################
# VERIFY ROBUST ANALYSIS IS GONE       #
########################################

dependence_summary = openxlsx::read.xlsx(
  paste0(final_asys, "segregant_dependence.xlsx"),
  sheet = "Summary",
  check.names = FALSE
)
dependence_locus = openxlsx::read.xlsx(
  paste0(final_asys, "segregant_dependence.xlsx"),
  sheet = "Locus results",
  check.names = FALSE
)

if (
  any(grepl("robust", names(dependence_summary), ignore.case = TRUE)) ||
  any(grepl("robust", names(dependence_locus), ignore.case = TRUE))
) {
  stop("Cluster-robust output remains in final segregant-dependence workbook.")
}

required_dependence_columns = c(
  "segregant_random_fit_success",
  "segregant_random_full_singular",
  "segregant_random_interaction_effect",
  "segregant_random_interaction_p"
)
if (!all(required_dependence_columns %in% names(dependence_locus))) {
  stop("Crossed-segregant diagnostic columns are missing from the final locus table.")
}
if (sum(dependence_locus$segregant_random_full_singular %in% TRUE) != 928L) {
  stop("Expected all 928 crossed-segregant full models to be singular.")
}

########################################
# FINAL STATE                         #
########################################

cat("\nFINAL PRE-POWER SUPPORT CHECKS PASSED\n")
cat("Workbooks: 5\n")
cat("Supporting-data sheets: 26\n")
cat("Continuous effects: fold differences for expression and PPI\n")
cat("Hotspot support rows: 2 (MKT1, BUL2)\n")
cat("Dependence sensitivity: crossed-segregant random intercept only\n")
cat("Directory: ", final_asys, "\n", sep = "")
