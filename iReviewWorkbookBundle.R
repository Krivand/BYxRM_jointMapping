###########################################
###########################################
## GENETICS REVIEW WORKBOOK OUTPUT BUNDLE ##
## FIVE WORKBOOKS WITH SUPPORTING TABS     ##
###########################################
###########################################

source("/PHShome/kv1026/p05/remap_functions.R")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package openxlsx is required to build reviewer workbooks.")
}

repo_wd = "/PHShome/kv1026/p05/"
review_asys = paste0(asys, "review/")
review_imed = paste0(imed, "review/")
final_asys = paste0(review_asys, "final/")
build_asys = paste0(review_asys, ".final_workbook_build/")

ensure_directory(review_asys)

final_files = c(
  "replicate_filtering.xlsx",
  "supported_gene_characteristics.xlsx",
  "interaction_GO.xlsx",
  "power_analysis.xlsx",
  "segregant_dependence.xlsx"
)

########################################
# HELPERS                              #
########################################

source_review_script = function(filename) {
  path = paste0(repo_wd, filename)
  if (!file.exists(path)) stop("Missing reviewer script: ", path)
  env = new.env(parent = globalenv())
  sys.source(path, envir = env)
  env
}

as_plain_df = function(x) {
  x = as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  for (j in seq_along(x)) {
    if (is.factor(x[[j]])) x[[j]] = as.character(x[[j]])
  }
  x
}

sheet_widths = function(dat) {
  if (ncol(dat) == 0L) return(numeric(0))
  vapply(
    seq_along(dat),
    function(j) {
      x = dat[[j]]
      txt = as.character(x)
      txt[is.na(txt)] = "NA"
      observed = if (length(txt) > 0L) max(nchar(txt), na.rm = TRUE) else 0L
      max(10, min(45, max(nchar(names(dat)[[j]]), observed) + 2))
    },
    numeric(1)
  )
}

write_workbook = function(sheets, filename) {
  if (is.null(names(sheets)) || any(names(sheets) == "")) {
    stop("Every workbook sheet must have a name.")
  }
  if (any(nchar(names(sheets)) > 31L)) {
    stop("Excel sheet names must not exceed 31 characters: ", paste(names(sheets), collapse = ", "))
  }
  if (anyDuplicated(names(sheets))) {
    stop("Workbook contains duplicated sheet names: ", paste(names(sheets), collapse = ", "))
  }

  wb = openxlsx::createWorkbook()

  for (i in seq_along(sheets)) {
    sheet_name = names(sheets)[[i]]
    dat = as_plain_df(sheets[[i]])

    if (ncol(dat) == 0L) {
      stop("Workbook sheet has zero columns: ", sheet_name)
    }

    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeDataTable(
      wb,
      sheet = sheet_name,
      x = dat,
      startRow = 1,
      startCol = 1,
      tableName = paste0("ReviewerTable", sprintf("%02d", i)),
      tableStyle = "TableStyleMedium2",
      withFilter = TRUE
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
      style = openxlsx::createStyle(wrapText = TRUE, valign = "top"),
      rows = seq_len(max(1L, nrow(dat) + 1L)),
      cols = seq_len(ncol(dat)),
      gridExpand = TRUE,
      stack = TRUE
    )
  }

  path = paste0(build_asys, filename)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

make_clean_power_table = function(
  study,
  n,
  pve50,
  pve80,
  observed,
  predicted
) {
  out = data.frame(
    study = as.character(study),
    n = as.numeric(n),
    pve_50pct_power_lod3_03_percent = as.numeric(pve50),
    pve_80pct_power_lod3_03_percent = as.numeric(pve80),
    figure3_observed_albert_query_percent = as.numeric(observed),
    power_predicted_same_direction_lod4_5_percent = as.numeric(predicted),
    observed_minus_predicted_lod4_5_pp = as.numeric(observed) - as.numeric(predicted),
    stringsAsFactors = FALSE
  )

  target =
    out$study != "Albert_2018" &
    is.finite(out$figure3_observed_albert_query_percent) &
    is.finite(out$power_predicted_same_direction_lod4_5_percent)

  if (sum(target) != 9L) {
    stop("Expected nine finite Albert-to-target rediscovery comparisons.")
  }

  mean_row = data.frame(
    study = "MEAN_9_ALBERT_TARGETS",
    n = NA_real_,
    pve_50pct_power_lod3_03_percent = NA_real_,
    pve_80pct_power_lod3_03_percent = NA_real_,
    figure3_observed_albert_query_percent = mean(
      out$figure3_observed_albert_query_percent[target]
    ),
    power_predicted_same_direction_lod4_5_percent = mean(
      out$power_predicted_same_direction_lod4_5_percent[target]
    ),
    observed_minus_predicted_lod4_5_pp = mean(
      out$observed_minus_predicted_lod4_5_pp[target]
    ),
    stringsAsFactors = FALSE
  )

  rbind(out, mean_row)
}

make_power_support_table = function(
  observed,
  predicted,
  scale08,
  scale06,
  benchmark,
  simulation
) {
  keep =
    is.finite(observed) &
    is.finite(predicted) &
    is.finite(scale08) &
    is.finite(scale06)

  if (sum(keep) != 9L || nrow(benchmark) != 1L) {
    stop("Unexpected power-validation inputs.")
  }

  simulation_error = max(
    c(simulation$any_error, simulation$same_error),
    na.rm = TRUE
  )

  data.frame(
    metric = c(
      "mean_observed_figure3_rediscovery",
      "mean_predicted_same_direction_lod4_5",
      "mean_observed_minus_predicted_lod4_5",
      "mean_predicted_lod4_5_effect_scale_0_8",
      "mean_predicted_lod4_5_effect_scale_0_6",
      "observed_vs_predicted_pearson_r",
      "observed_vs_predicted_spearman_rho",
      "albert_public_power_at_2_5pct_pve",
      "our_power_at_2_5pct_pve",
      "max_analytic_engine_difference_0_5_to_10pct_pve",
      "max_direct_simulation_difference"
    ),
    value = c(
      100 * mean(observed[keep]),
      100 * mean(predicted[keep]),
      100 * mean(observed[keep] - predicted[keep]),
      100 * mean(scale08[keep]),
      100 * mean(scale06[keep]),
      stats::cor(observed[keep], predicted[keep], method = "pearson"),
      stats::cor(observed[keep], predicted[keep], method = "spearman"),
      100 * benchmark$public_code_power_t_test[[1]],
      100 * benchmark$current_engine_any_detection[[1]],
      100 * benchmark$engine_grid_max_abs_difference[[1]],
      100 * simulation_error
    ),
    unit = c(
      "percent",
      "percent",
      "percentage_points",
      "percent",
      "percent",
      "correlation",
      "correlation",
      "percent",
      "percent",
      "percentage_points",
      "percentage_points"
    ),
    stringsAsFactors = FALSE
  )
}

removed_rows = function(dataset, entity_type, before, after) {
  removed = sort(setdiff(as.character(before), as.character(after)))
  if (length(removed) == 0L) {
    return(data.frame(
      Dataset = character(0),
      Entity = character(0),
      ID = character(0),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    Dataset = dataset,
    Entity = entity_type,
    ID = removed,
    stringsAsFactors = FALSE
  )
}

########################################
# START FROM A CLEAN NONFINAL STAGING  #
########################################

staging_dirs = c(
  paste0(review_asys, "missingness/"),
  paste0(review_asys, "genes/"),
  paste0(review_asys, "go/"),
  paste0(review_asys, "power_effect_curve/"),
  paste0(review_asys, "interaction_dependence/")
)

for (path in staging_dirs) {
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
}
if (dir.exists(build_asys)) unlink(build_asys, recursive = TRUE, force = TRUE)
ensure_directory(build_asys)

########################################
# REGENERATE AND VALIDATE ANALYSES     #
########################################

message("Running replicate/missingness review analysis...")
missing_env = source_review_script("iReviewMissingness.R")

message("Running supported-gene review analysis...")
genes_env = source_review_script("jReviewGenes.R")

message("Running GO review analysis...")
go_env = source_review_script("jReviewGO.R")

message("Running individual-study power review analysis...")
power_env = source_review_script("iReviewPowerEffectCurve.R")
power_validation_env = source_review_script("iReviewPowerEffectCurveValidate.R")
figure3_env = source_review_script("iReviewPowerFigure3Compare.R")

message("Running shared-segregant review sensitivity...")
dependence_env = source_review_script("iReviewInteractionDependence.R")
invisible(source_review_script("iReviewInteractionDependenceValidate.R"))

########################################
# 1. REPLICATE FILTERING WORKBOOK      #
########################################

b02_1 = missing_env$brem_channel_sets("GDS91.soft", "BY-")
b02_2 = missing_env$brem_channel_sets("GDS92.soft", "-BY")
b05_1 = missing_env$brem_channel_sets("GDS1115.soft", "Cy5")
b05_2 = missing_env$brem_channel_sets("GDS1116.soft", "Cy3")

b02_genes_before = intersect(b02_1$before_genes, b02_2$before_genes)
b02_genes_after = intersect(b02_1$after_genes, b02_2$after_genes)
b02_segs_before = intersect(b02_1$before_strains, b02_2$before_strains)
b02_segs_after = intersect(b02_1$after_strains, b02_2$after_strains)

b05_genes_before = intersect(b05_1$before_genes, b05_2$before_genes)
b05_genes_after = intersect(b05_1$after_genes, b05_2$after_genes)
b05_segs_before = intersect(b05_1$before_strains, b05_2$before_strains)
b05_segs_after = intersect(b05_1$after_strains, b05_2$after_strains)

khan_raw = read.delim(
  paste0(raw, "Khan_2009_protein.tsv"),
  sep = "\t",
  row.names = 1,
  check.names = FALSE
) %>%
  as.matrix()
khan_filtered = rmIterNA(khan_raw, t = threshNA)

foss_raw = readxl::read_excel(
  paste0(raw, "foss_11_protein_measures.xls")
) %>%
  as.data.frame()
rownames(foss_raw) = foss_raw$`systematic name`
foss_raw = foss_raw[, -c(1:3)] %>%
  mutate_all(as.numeric) %>%
  as.matrix()
colnames(foss_raw) = gsub("_", ".", colnames(foss_raw))
foss_filtered = rmIterNA(foss_raw, t = threshNA)

missingness_removed = do.call(
  rbind,
  list(
    removed_rows("Brem 2002", "gene", b02_genes_before, b02_genes_after),
    removed_rows("Brem 2002", "segregant", b02_segs_before, b02_segs_after),
    removed_rows("Brem 2005", "gene", b05_genes_before, b05_genes_after),
    removed_rows("Brem 2005", "segregant", b05_segs_before, b05_segs_after),
    removed_rows("Khan 2009", "gene", rownames(khan_raw), rownames(khan_filtered)),
    removed_rows("Khan 2009", "segregant", colnames(khan_raw), colnames(khan_filtered)),
    removed_rows("Foss 2011", "gene", rownames(foss_raw), rownames(foss_filtered)),
    removed_rows("Foss 2011", "segregant", colnames(foss_raw), colnames(foss_filtered))
  )
)
rownames(missingness_removed) = NULL

replicate_sheets = list(
  Summary = missing_env$missingness_summary,
  Removed_items = missingness_removed
)
write_workbook(replicate_sheets, "replicate_filtering.xlsx")

rm(
  b02_1,
  b02_2,
  b05_1,
  b05_2,
  khan_raw,
  khan_filtered,
  foss_raw,
  foss_filtered
)
gc()

########################################
# 2. SUPPORTED-GENE WORKBOOK           #
########################################

gene_level = genes_env$expression_universe
gene_level$supported_for_joint_mapping =
  gene_level$gene %in% genes_env$supported_genes

gene_level = gene_level[
  ,
  c(
    "gene",
    "supported_for_joint_mapping",
    "expression",
    "essential",
    "hasHumanHomolog",
    "PPI"
  ),
  drop = FALSE
]

names(gene_level) = c(
  "gene",
  "supported_for_joint_mapping",
  "mean_log2_TPM",
  "essential",
  "has_human_homolog",
  "PPI_partners"
)

gene_sheets = list(
  Summary = genes_env$review_table,
  Gene_level = gene_level
)
write_workbook(gene_sheets, "supported_gene_characteristics.xlsx")

########################################
# 3. GO WORKBOOK                       #
########################################

go_summary = go_env$key_go
go_contrast = go_env$ribosome_contrast

contrast_row = go_contrast[
  go_contrast$GO_term == "ribosomal subunit",
  ,
  drop = FALSE
]
key_row = which(go_summary$`GO term` == "ribosomal subunit")

if (nrow(contrast_row) != 1L || length(key_row) != 1L) {
  stop("Could not identify the ribosomal-subunit GO contrast.")
}

go_summary$`Discordant vs other interaction QTLs` = NA_character_
go_summary$`Class contrast FDR` = NA_character_
go_summary$`Discordant vs other interaction QTLs`[[key_row]] = paste(
  contrast_row$discordant[[1]],
  "vs",
  contrast_row$other_interactions[[1]]
)
go_summary$`Class contrast FDR`[[key_row]] = formatC(
  contrast_row$fisher_q[[1]],
  format = "e",
  digits = 2
)

go_gene_sets = do.call(
  rbind,
  lapply(
    names(go_env$gene_sets),
    function(set_name) {
      data.frame(
        gene_set = set_name,
        gene = go_env$gene_sets[[set_name]],
        stringsAsFactors = FALSE
      )
    }
  )
)
rownames(go_gene_sets) = NULL

go_all = go_env$all_go_results

go_sheets = list(
  Summary = go_summary,
  All_strict_64 = go_all[go_all$gene_set == "all_interaction", , drop = FALSE],
  Discordant_37 = go_all[go_all$gene_set == "discordant", , drop = FALSE],
  Protein_specific_14 = go_all[go_all$gene_set == "protein_specific", , drop = FALSE],
  mRNA_specific_9 = go_all[go_all$gene_set == "mRNA_specific", , drop = FALSE],
  Same_direction_4 = go_all[go_all$gene_set == "same_direction", , drop = FALSE],
  Significant_all = go_env$significant_go,
  Ribosome_contrast = go_contrast,
  Gene_sets = go_gene_sets
)
write_workbook(go_sheets, "interaction_GO.xlsx")

########################################
# 4. POWER WORKBOOK                    #
########################################

critical = power_env$critical_r_power_table
crossstudy = power_env$albert_reference_crossstudy_table
benchmark = power_env$albert_2018_power_benchmark
figure3 = figure3_env$comparison
figure3_fit = figure3_env$fit_summary
simulation = power_validation_env$simulation
validation_checks = power_validation_env$validation_summary

critical_303 = critical[abs(critical$lod - 3.03) < 1e-10, , drop = FALSE]
if (nrow(critical_303) != 10L || nrow(figure3) != 9L) {
  stop("Unexpected study counts in power workbook inputs.")
}

observed = rep(NA_real_, nrow(critical_303))
predicted = rep(NA_real_, nrow(critical_303))
figure_index = match(figure3$target, critical_303$study)
if (any(is.na(figure_index))) {
  stop("Could not align Figure 3 targets to power-study rows.")
}
observed[figure_index] = figure3$observed_figure3_rediscovery
predicted[figure_index] = figure3$predicted_same_direction_lod4_5

power_summary = make_clean_power_table(
  study = critical_303$study,
  n = critical_303$n,
  pve50 = 100 * critical_303$pve_50pct_power,
  pve80 = 100 * critical_303$pve_80pct_power,
  observed = 100 * observed,
  predicted = 100 * predicted
)

cross_index = match(figure3$target, crossstudy$target)
if (any(is.na(cross_index))) {
  stop("Could not align Figure 3 targets to Albert-reference power rows.")
}

power_validation_summary = make_power_support_table(
  observed = figure3$observed_figure3_rediscovery,
  predicted = figure3$predicted_same_direction_lod4_5,
  scale08 = crossstudy$expected_detection_lod4_5_scale0_8[cross_index],
  scale06 = crossstudy$expected_detection_lod4_5_scale0_6[cross_index],
  benchmark = benchmark,
  simulation = simulation
)

analytic_engine_grid = data.frame(
  pve = power_validation_env$engine_grid_pve,
  public_code_power = power_validation_env$engine_public,
  current_engine_power = power_validation_env$engine_current,
  current_minus_public = power_validation_env$engine_diff,
  absolute_difference = abs(power_validation_env$engine_diff),
  stringsAsFactors = FALSE
)

power_sheets = list(
  Study_summary = power_summary,
  Validation_summary = power_validation_summary,
  Power_grid = critical,
  Albert_crossstudy = crossstudy,
  Figure3_targets = figure3,
  Figure3_fit = figure3_fit,
  Albert_benchmark = benchmark,
  Analytic_engine_grid = analytic_engine_grid,
  Direct_simulation = simulation,
  Validation_checks = validation_checks
)
write_workbook(power_sheets, "power_analysis.xlsx")

########################################
# 5. SHARED-SEGREGANT WORKBOOK         #
########################################

dependence_sheets = list(
  Summary = dependence_env$segregant_dependence_summary,
  Locus_level = dependence_env$segregant_dependence_locus_results,
  Panel_overlap = dependence_env$segregant_panel_overlap
)
write_workbook(dependence_sheets, "segregant_dependence.xlsx")

########################################
# VALIDATE WORKBOOK BUILD              #
########################################

present = sort(list.files(build_asys, pattern = "\\.xlsx$"))
if (!identical(present, sort(final_files))) {
  stop(
    "Temporary reviewer workbook bundle is incomplete: ",
    paste(present, collapse = ", ")
  )
}

expected_sheets = list(
  replicate_filtering.xlsx = names(replicate_sheets),
  supported_gene_characteristics.xlsx = names(gene_sheets),
  interaction_GO.xlsx = names(go_sheets),
  power_analysis.xlsx = names(power_sheets),
  segregant_dependence.xlsx = names(dependence_sheets)
)

for (filename in final_files) {
  path = paste0(build_asys, filename)
  observed_sheets = openxlsx::getSheetNames(path)
  if (!identical(observed_sheets, expected_sheets[[filename]])) {
    stop(
      "Unexpected workbook sheets in ", filename,
      ": ", paste(observed_sheets, collapse = ", ")
    )
  }

  first_sheet = observed_sheets[[1]]
  check = openxlsx::read.xlsx(
    path,
    sheet = first_sheet,
    check.names = FALSE
  )
  if (nrow(check) < 1L || ncol(check) < 1L) {
    stop("First workbook sheet is empty: ", filename, " / ", first_sheet)
  }
}

########################################
# PROMOTE TRANSACTIONALLY              #
########################################

if (dir.exists(final_asys)) {
  unlink(final_asys, recursive = TRUE, force = TRUE)
}

if (!file.rename(build_asys, final_asys)) {
  stop("Could not promote validated reviewer workbook bundle to final/.")
}

########################################
# MANIFEST                             #
########################################

purpose = list(
  replicate_filtering.xlsx = c(
    Summary = "Reviewer-facing replicate and missingness counts",
    Removed_items = "Exact genes and segregants removed by the iterative 25% missingness filter"
  ),
  supported_gene_characteristics.xlsx = c(
    Summary = "Reviewer-facing representativeness comparisons",
    Gene_level = "Gene-level Albert et al. 2018 characteristics for the full 5,720-gene comparison universe"
  ),
  interaction_GO.xlsx = c(
    Summary = "Reviewer-facing GO terms and ribosomal-subunit class contrast",
    All_strict_64 = "All tested GO terms for the 64 strict interaction-QTL genes",
    Discordant_37 = "All tested GO terms for the 37 discordant strict interaction-QTL genes",
    Protein_specific_14 = "All tested GO terms for the 14 protein-specific strict interaction-QTL genes",
    mRNA_specific_9 = "All tested GO terms for the 9 mRNA-specific strict interaction-QTL genes",
    Same_direction_4 = "All tested GO terms for the 4 same-direction strict interaction-QTL genes",
    Significant_all = "All GO terms passing FDR < 0.05 across gene-set by ontology analyses",
    Ribosome_contrast = "Direct discordant-versus-other-interaction contrasts for six ribosome-focused GO terms",
    Gene_sets = "Genes entering the total strict-interaction and subclass GO analyses"
  ),
  power_analysis.xlsx = c(
    Study_summary = "Reviewer-facing sample-size, detectable-PVE, observed rediscovery, and predicted rediscovery summary",
    Validation_summary = "Reviewer-facing power validation and attenuation sensitivity metrics",
    Power_grid = "Full 10-study by four-LOD detectable-effect grid",
    Albert_crossstudy = "Albert-reference cross-study detection predictions and sensitivity scales",
    Figure3_targets = "Per-target observed Figure 3 rediscovery and power predictions",
    Figure3_fit = "Across-target calibration summaries at each stipulated LOD",
    Albert_benchmark = "Reproduction of the Albert et al. 2018 public power benchmark",
    Analytic_engine_grid = "Public-code versus current analytic power across 0.5-10% PVE",
    Direct_simulation = "Direct phenotype-simulation validation across 27 settings",
    Validation_checks = "Programmatic validation checks for the power analysis"
  ),
  segregant_dependence.xlsx = c(
    Summary = "Reviewer-facing shared-segregant sensitivity summaries",
    Locus_level = "All 928 targeted locus-level original, cluster-robust, and crossed-random-effect results",
    Panel_overlap = "Pairwise shared-segregant counts for the seven joint-mapping datasets"
  )
)

manifest_rows = list()
manifest_i = 1L
for (filename in final_files) {
  sheets = expected_sheets[[filename]]
  for (sheet_name in sheets) {
    dat = switch(
      filename,
      "replicate_filtering.xlsx" = replicate_sheets[[sheet_name]],
      "supported_gene_characteristics.xlsx" = gene_sheets[[sheet_name]],
      "interaction_GO.xlsx" = go_sheets[[sheet_name]],
      "power_analysis.xlsx" = power_sheets[[sheet_name]],
      "segregant_dependence.xlsx" = dependence_sheets[[sheet_name]]
    )
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

########################################
# REMOVE NONFINAL REVIEW ARTIFACTS     #
########################################

review_entries = list.files(
  review_asys,
  full.names = TRUE,
  recursive = FALSE,
  all.files = FALSE
)

keep_names = c(
  basename(sub("/$", "", final_asys)),
  "review_output_manifest.tsv"
)
stale = review_entries[!(basename(review_entries) %in% keep_names)]
if (length(stale) > 0L) {
  unlink(stale, recursive = TRUE, force = TRUE)
}

if (dir.exists(review_imed)) {
  unlink(review_imed, recursive = TRUE, force = TRUE)
}

final_present = sort(list.files(final_asys, pattern = "\\.xlsx$"))
if (!identical(final_present, sort(final_files))) {
  stop("Final reviewer workbook bundle changed unexpectedly during cleanup.")
}

cat("\nFINAL REVIEWER WORKBOOK BUNDLE\n")
cat("Substantive workbooks: 5\n")
cat("Supporting sheets are retained inside each workbook.\n")
cat("Directory: ", final_asys, "\n", sep = "")
cat("Manifest: ", paste0(review_asys, "review_output_manifest.tsv"), "\n", sep = "")
cat("All nonfinal review outputs and intermediates removed.\n")
