############################################
############################################
## GENETICS REVIEW ALL-PAIR ATTENUATION   ##
## FINAL OUTPUT SUPPORT                   ##
############################################
############################################

source("/PHShome/kv1026/p05/remap_functions.R")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package openxlsx is required for reviewer workbook support.")
}

review_asys = paste0(asys, "review/")
power_asys = paste0(review_asys, "power_effect_curve/")
allpair_asys = paste0(power_asys, "all_pairs/")
final_asys = paste0(review_asys, "final/")
manifest_path = paste0(review_asys, "review_output_manifest.tsv")
power_workbook_path = paste0(final_asys, "power_analysis.xlsx")
final_tsv_path = paste0(final_asys, "pairwise_power_predictions.tsv")
attenuation_path = paste0(allpair_asys, "pairwise_power_attenuation.tsv")

required_files = c(
  power_workbook_path,
  final_tsv_path,
  manifest_path,
  attenuation_path
)
missing_files = required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

read_tsv = function(path) {
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

pairwise_power = read_tsv(final_tsv_path)
attenuation = read_tsv(attenuation_path)

if (
  nrow(pairwise_power) != 90L ||
  nrow(attenuation) != 90L ||
  anyDuplicated(pairwise_power[, c("fromData", "toData")]) ||
  anyDuplicated(attenuation[, c("fromData", "toData")])
) {
  stop("Expected 90 unique directed rows in both power tables.")
}

new_columns = c(
  "expectedRediscoveryLOD4_5_scale0_8",
  "expectedRediscoveryLOD4_5_scale0_6"
)
if (!all(new_columns %in% names(attenuation))) {
  stop("Attenuation table is missing required scaled prediction columns.")
}

key = paste(pairwise_power$fromData, pairwise_power$toData, sep = "\r")
attenuation_key = paste(attenuation$fromData, attenuation$toData, sep = "\r")
idx = match(key, attenuation_key)
if (any(is.na(idx))) stop("Could not align attenuation rows to final power rows.")
attenuation = attenuation[idx, , drop = FALSE]

pairwise_power = pairwise_power[
  ,
  setdiff(names(pairwise_power), new_columns),
  drop = FALSE
]

insert_after = match("expectedRediscoveryLOD4_5", names(pairwise_power))
if (is.na(insert_after)) {
  stop("Final power table lacks expectedRediscoveryLOD4_5.")
}

left_names = names(pairwise_power)[seq_len(insert_after)]
right_names = names(pairwise_power)[
  seq.int(insert_after + 1L, ncol(pairwise_power))
]
if (insert_after == ncol(pairwise_power)) right_names = character(0)

pairwise_power = data.frame(
  pairwise_power[, left_names, drop = FALSE],
  expectedRediscoveryLOD4_5_scale0_8 = attenuation$expectedRediscoveryLOD4_5_scale0_8,
  expectedRediscoveryLOD4_5_scale0_6 = attenuation$expectedRediscoveryLOD4_5_scale0_6,
  pairwise_power[, right_names, drop = FALSE],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  nrow(pairwise_power) != 90L ||
  anyDuplicated(pairwise_power[, c("fromData", "toData")]) ||
  !all(new_columns %in% names(pairwise_power)) ||
  any(!is.finite(as.matrix(pairwise_power[, new_columns]))) ||
  any(as.matrix(pairwise_power[, new_columns]) < 0) ||
  any(as.matrix(pairwise_power[, new_columns]) > 1) ||
  any(
    pairwise_power$expectedRediscoveryLOD4_5 <
      pairwise_power$expectedRediscoveryLOD4_5_scale0_8
  ) ||
  any(
    pairwise_power$expectedRediscoveryLOD4_5_scale0_8 <
      pairwise_power$expectedRediscoveryLOD4_5_scale0_6
  )
) {
  stop("Final attenuation columns failed validation.")
}

########################
# WORKBOOK SUPPORT SHEET #
########################

sheet_name = "All pair power"
wb = openxlsx::loadWorkbook(power_workbook_path)
existing_sheets = openxlsx::sheets(wb)
if (!(sheet_name %in% existing_sheets)) {
  stop("Power workbook lacks the All pair power sheet.")
}
openxlsx::removeWorksheet(wb, sheet_name)
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
openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)

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

workbook_tmp = paste0(power_workbook_path, ".attenuation.tmp.xlsx")
if (file.exists(workbook_tmp)) unlink(workbook_tmp, force = TRUE)
openxlsx::saveWorkbook(wb, workbook_tmp, overwrite = TRUE)

sheet_check = openxlsx::read.xlsx(
  workbook_tmp,
  sheet = sheet_name,
  check.names = FALSE
)
if (
  nrow(sheet_check) != 90L ||
  !all(new_columns %in% names(sheet_check)) ||
  anyDuplicated(sheet_check[, c("fromData", "toData")])
) {
  stop("Updated All pair power workbook sheet failed verification.")
}

if (!file.rename(workbook_tmp, power_workbook_path)) {
  stop("Could not promote attenuation-updated power workbook.")
}

########################
# FINAL TSV             #
########################

write.table(
  pairwise_power,
  file = final_tsv_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

roundtrip = read_tsv(final_tsv_path)
if (
  nrow(roundtrip) != 90L ||
  !identical(roundtrip$fromData, pairwise_power$fromData) ||
  !identical(roundtrip$toData, pairwise_power$toData) ||
  !all(new_columns %in% names(roundtrip))
) {
  stop("Updated final power TSV failed round-trip verification.")
}

########################
# MANIFEST              #
########################

manifest = read_tsv(manifest_path)
row_i = which(
  manifest$workbook == "power_analysis.xlsx" &
    manifest$sheet == sheet_name
)
if (length(row_i) != 1L) {
  stop("Could not identify exactly one All pair power manifest row.")
}
manifest$rows[[row_i]] = nrow(pairwise_power)
manifest$columns[[row_i]] = ncol(pairwise_power)
manifest$purpose[[row_i]] = paste(
  "Merge-ready power predictions for all 90 directed pairs among the 10",
  "individual-segregant datasets, including 0.8 and 0.6 effect-scale",
  "sensitivities; join on fromData and toData"
)
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

cat("\nFINAL ALL-PAIR ATTENUATION SUPPORT READY\n")
cat("Directed pairs: 90\n")
cat("Added columns: ", paste(new_columns, collapse = ", "), "\n", sep = "")
cat("Workbook sheet: All pair power\n")
cat("TSV: ", final_tsv_path, "\n", sep = "")
cat("Supporting-data sheets: 27\n")
