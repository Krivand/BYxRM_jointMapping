########################################
########################################
## GENETICS REVIEW ALL-PAIR POWER     ##
## FINAL OUTPUT SUPPORT               ##
########################################
########################################

source("/PHShome/kv1026/p05/remap_functions.R")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package openxlsx is required for reviewer workbook support.")
}

review_asys = paste0(asys, "review/")
final_asys = paste0(review_asys, "final/")
manifest_path = paste0(review_asys, "review_output_manifest.tsv")
power_workbook_path = paste0(final_asys, "power_analysis.xlsx")
final_tsv_path = paste0(final_asys, "pairwise_power_for_frank.tsv")

if (!exists("allpairs_env", inherits = TRUE)) {
  stop("All-pair finalizer must run after iReviewPowerAllPairs.R.")
}
if (!exists(
  "pairwise_power_for_frank",
  envir = allpairs_env,
  inherits = FALSE
)) {
  stop("All-pair analysis did not generate pairwise_power_for_frank.")
}
if (!file.exists(power_workbook_path)) {
  stop("Missing finalized power workbook: ", power_workbook_path)
}
if (!file.exists(manifest_path)) {
  stop("Missing reviewer output manifest: ", manifest_path)
}

pairwise_power = get(
  "pairwise_power_for_frank",
  envir = allpairs_env,
  inherits = FALSE
)
pairwise_power = as.data.frame(
  pairwise_power,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns = c(
  "fromData",
  "toData",
  "fromData_analyte",
  "toData_analyte",
  "fromData_medium",
  "toData_medium",
  "analyteCombo",
  "mediumSame",
  "fromData_n",
  "toData_n",
  "nSourceQTLsUsedForPower",
  "meanAbsSourceEffect",
  "medianAbsSourceEffect",
  "expectedRediscoveryLOD2",
  "expectedRediscoveryLOD3_03",
  "expectedRediscoveryLOD4_5",
  "expectedRediscoveryLOD7"
)

if (
  nrow(pairwise_power) != 90L ||
  !all(required_columns %in% names(pairwise_power)) ||
  anyDuplicated(pairwise_power[, c("fromData", "toData")]) ||
  any(pairwise_power$fromData == pairwise_power$toData)
) {
  stop("Merge-ready all-pair power table failed final checks.")
}

########################
# WORKBOOK SUPPORT SHEET #
########################

sheet_name = "All pair power"
wb = openxlsx::loadWorkbook(power_workbook_path)
existing_sheets = openxlsx::sheets(wb)
if (sheet_name %in% existing_sheets) {
  openxlsx::removeWorksheet(wb, sheet_name)
}
openxlsx::addWorksheet(wb, sheet_name)

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
  sheet = sheet_name,
  x = pairwise_power,
  startRow = 1,
  startCol = 1,
  headerStyle = header_style,
  keepNA = TRUE,
  na.string = "NA"
)
openxlsx::addFilter(
  wb,
  sheet = sheet_name,
  row = 1,
  cols = seq_len(ncol(pairwise_power))
)
openxlsx::freezePane(
  wb,
  sheet = sheet_name,
  firstRow = TRUE
)

widths = vapply(
  seq_along(pairwise_power),
  function(j) {
    txt = as.character(pairwise_power[[j]])
    txt[is.na(txt)] = "NA"
    observed = max(nchar(txt), na.rm = TRUE)
    max(10, min(40, max(nchar(names(pairwise_power)[[j]]), observed) + 2))
  },
  numeric(1)
)
openxlsx::setColWidths(
  wb,
  sheet = sheet_name,
  cols = seq_len(ncol(pairwise_power)),
  widths = widths
)
openxlsx::addStyle(
  wb,
  sheet = sheet_name,
  style = body_style,
  rows = seq_len(nrow(pairwise_power) + 1L),
  cols = seq_len(ncol(pairwise_power)),
  gridExpand = TRUE,
  stack = TRUE
)

workbook_tmp = paste0(power_workbook_path, ".tmp.xlsx")
if (file.exists(workbook_tmp)) unlink(workbook_tmp, force = TRUE)
openxlsx::saveWorkbook(wb, workbook_tmp, overwrite = TRUE)

sheet_check = openxlsx::read.xlsx(
  workbook_tmp,
  sheet = sheet_name,
  check.names = FALSE
)
if (
  nrow(sheet_check) != 90L ||
  !all(required_columns %in% names(sheet_check)) ||
  anyDuplicated(sheet_check[, c("fromData", "toData")])
) {
  stop("Final all-pair power workbook sheet failed verification.")
}

if (!file.rename(workbook_tmp, power_workbook_path)) {
  stop("Could not promote power workbook with all-pair support sheet.")
}

########################
# MERGE-READY TSV       #
########################

write.table(
  pairwise_power,
  file = final_tsv_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

roundtrip = read.delim(
  final_tsv_path,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (
  nrow(roundtrip) != 90L ||
  !identical(roundtrip$fromData, pairwise_power$fromData) ||
  !identical(roundtrip$toData, pairwise_power$toData)
) {
  stop("Final merge-ready TSV failed round-trip verification.")
}

########################
# UPDATE MANIFEST       #
########################

manifest = read.delim(
  manifest_path,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
manifest = manifest[
  !(
    manifest$workbook == "power_analysis.xlsx" &
    manifest$sheet == sheet_name
  ),
  ,
  drop = FALSE
]

pairwise_manifest = data.frame(
  workbook = "power_analysis.xlsx",
  sheet = sheet_name,
  rows = nrow(pairwise_power),
  columns = ncol(pairwise_power),
  purpose = paste(
    "Merge-ready power predictions for all 90 directed pairs among the 10",
    "individual-segregant datasets; join on fromData and toData"
  ),
  stringsAsFactors = FALSE
)
manifest = rbind(manifest, pairwise_manifest)
rownames(manifest) = NULL

if (nrow(manifest) != 27L) {
  stop("Expected exactly 27 final supporting-data sheets.")
}

write.table(
  manifest,
  file = manifest_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("\nFINAL ALL-PAIR POWER SUPPORT READY\n")
cat("Directed pairs: 90\n")
cat("Workbook sheet: All pair power\n")
cat("Merge keys: fromData + toData\n")
cat("TSV: ", final_tsv_path, "\n", sep = "")
cat("Supporting-data sheets: 27\n")
