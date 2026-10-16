########################################
########################################
## GENETICS REVIEW WORKBOOK FINALIZER ##
## CLEAN SUPPORTING DATA PRESENTATION ##
########################################
########################################

source("/PHShome/kv1026/p05/remap_functions.R")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package openxlsx is required to finalize reviewer workbooks.")
}

review_asys = paste0(asys, "review/")
final_asys = paste0(review_asys, "final/")
build_asys = paste0(review_asys, ".final_presentation_build/")

final_files = c(
  "replicate_filtering.xlsx",
  "supported_gene_characteristics.xlsx",
  "interaction_GO.xlsx",
  "power_analysis.xlsx",
  "segregant_dependence.xlsx"
)

required_objects = c(
  "replicate_sheets",
  "gene_sheets",
  "go_env",
  "critical",
  "crossstudy",
  "benchmark",
  "figure3",
  "figure3_fit",
  "power_summary",
  "power_validation_env",
  "dependence_sheets",
  "make_power_support_table"
)

missing_objects = required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects) > 0L) {
  stop(
    "Workbook finalizer must run immediately after iReviewWorkbookBundle.R. ",
    "Missing objects: ", paste(missing_objects, collapse = ", ")
  )
}

########################################
# HELPERS                              #
########################################

clean_name = function(x) {
  x = gsub("contrast", "comparison", x, ignore.case = TRUE)
  x = gsub("engine", "method", x, ignore.case = TRUE)
  x = gsub("[ .]+", "_", x)
  x = gsub("_+", "_", x)
  x = gsub("^_|_$", "", x)
  x
}

clean_visible_terms = function(dat) {
  dat = as.data.frame(
    dat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(dat) = clean_name(names(dat))

  for (j in seq_along(dat)) {
    if (is.factor(dat[[j]])) dat[[j]] = as.character(dat[[j]])
    if (is.character(dat[[j]])) {
      dat[[j]] = gsub(
        "contrast",
        "comparison",
        dat[[j]],
        ignore.case = TRUE
      )
      dat[[j]] = gsub(
        "engine",
        "method",
        dat[[j]],
        ignore.case = TRUE
      )
    }
  }

  dat
}

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

write_plain_workbook = function(sheets, filename) {
  if (is.null(names(sheets)) || any(names(sheets) == "")) {
    stop("Every workbook sheet must have a name.")
  }
  if (any(nchar(names(sheets)) > 31L)) {
    stop("Excel sheet names must not exceed 31 characters.")
  }
  if (anyDuplicated(names(sheets))) {
    stop("Workbook contains duplicated sheet names: ", filename)
  }

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

    if (ncol(dat) < 1L) {
      stop("Sheet has no columns: ", filename, " / ", sheet_name)
    }

    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = dat,
      startRow = 1,
      startCol = 1,
      headerStyle = header_style,
      withFilter = nrow(dat) > 0L,
      keepNA = TRUE,
      na.string = "NA"
    )
    openxlsx::freezePane(
      wb,
      sheet = sheet_name,
      firstRow = TRUE
    )
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
    file = paste0(build_asys, filename),
    overwrite = TRUE
  )
}

########################################
# REGULAR DIRECT-SIMULATION GRID       #
########################################

study_sizes = unique(
  critical[
    ,
    c("study", "n"),
    drop = FALSE
  ]
)

if (
  nrow(study_sizes) != 10L ||
  anyDuplicated(study_sizes$study)
) {
  stop("Expected one sample size for each of the 10 studies.")
}

study_sizes = study_sizes[
  order(study_sizes$n, study_sizes$study),
  ,
  drop = FALSE
]

smallest_n = min(study_sizes$n)
median_n = stats::median(study_sizes$n)
largest_n = max(study_sizes$n)

size_table = data.frame(
  n = c(smallest_n, median_n, largest_n),
  data_size = c("smallest", "median", "largest"),
  stringsAsFactors = FALSE
)

if (
  nrow(size_table) != 3L ||
  anyDuplicated(size_table$n) ||
  !identical(size_table$n, c(33, 99, 1012))
) {
  stop(
    "Unexpected minimum, median, or maximum study sizes: ",
    paste(size_table$n, collapse = ", ")
  )
}

simulation_rho = seq(0.10, 0.40, by = 0.05)
simulation_lod = c(2, 4.5, 7)
simulation_n_rep = 1000L

simulation_grid = do.call(
  rbind,
  lapply(
    seq_len(nrow(size_table)),
    function(i) {
      grid = expand.grid(
        rho = simulation_rho,
        lod = simulation_lod,
        stringsAsFactors = FALSE
      )
      grid$n = size_table$n[[i]]
      grid$data_size = size_table$data_size[[i]]
      grid[
        ,
        c("n", "data_size", "rho", "lod"),
        drop = FALSE
      ]
    }
  )
)
rownames(simulation_grid) = NULL

simulation_worker = function(i) {
  row = simulation_grid[i, , drop = FALSE]
  out = power_validation_env$validate_direct_simulation(
    n = row$n[[1]],
    rho = row$rho[[1]],
    lod = row$lod[[1]],
    n_rep = simulation_n_rep,
    seed = power_validation_env$validation_seed + 1000L + i
  )
  out$data_size = row$data_size[[1]]
  out[
    ,
    c(
      "n",
      "data_size",
      "rho",
      "lod",
      "n_rep",
      "analytic_any",
      "simulated_any",
      "analytic_same",
      "simulated_same"
    ),
    drop = FALSE
  ]
}

if (
  .Platform$OS.type == "unix" &&
  power_validation_env$power_cores > 1L
) {
  simulation_list = parallel::mclapply(
    seq_len(nrow(simulation_grid)),
    simulation_worker,
    mc.cores = min(
      power_validation_env$power_cores,
      nrow(simulation_grid)
    ),
    mc.preschedule = TRUE
  )
} else {
  simulation_list = lapply(
    seq_len(nrow(simulation_grid)),
    simulation_worker
  )
}

simulation = do.call(rbind, simulation_list)
rownames(simulation) = NULL

simulation$any_error = abs(
  simulation$analytic_any - simulation$simulated_any
)
simulation$same_error = abs(
  simulation$analytic_same - simulation$simulated_same
)

any_se = sqrt(
  pmax(
    0,
    simulation$analytic_any * (1 - simulation$analytic_any) /
      simulation$n_rep
  )
)
same_se = sqrt(
  pmax(
    0,
    simulation$analytic_same * (1 - simulation$analytic_same) /
      simulation$n_rep
  )
)

simulation$any_tolerance = 0.005 + 4 * any_se
simulation$same_tolerance = 0.005 + 4 * same_se

if (
  nrow(simulation) != 63L ||
  !identical(sort(unique(simulation$rho)), simulation_rho) ||
  !identical(sort(unique(simulation$lod)), simulation_lod) ||
  !identical(unique(simulation$n_rep), simulation_n_rep) ||
  !all(simulation$any_error <= simulation$any_tolerance) ||
  !all(simulation$same_error <= simulation$same_tolerance)
) {
  stop("Regular direct-simulation grid failed its checks.")
}

########################################
# GO SUMMARY: COMPLETE CLASS COMPARISON#
########################################

go_summary = go_env$key_go
go_comparison = go_env$ribosome_contrast

comparison_index = match(
  tolower(go_summary$`GO term`),
  tolower(go_comparison$GO_term)
)

if (any(is.na(comparison_index))) {
  stop("Could not match every highlighted GO term to its class comparison.")
}

go_summary$`Discordant vs other interaction QTLs` = paste(
  go_comparison$discordant[comparison_index],
  "vs",
  go_comparison$other_interactions[comparison_index]
)
go_summary$`Class comparison FDR` = formatC(
  go_comparison$fisher_q[comparison_index],
  format = "e",
  digits = 2
)

if (
  any(is.na(go_summary$`Discordant vs other interaction QTLs`)) ||
  any(is.na(go_summary$`Class comparison FDR`))
) {
  stop("GO summary class-comparison columns are incomplete.")
}

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

go_sheets_final = list(
  "Summary" = go_summary,
  "All 64 interactions" = go_all[
    go_all$gene_set == "all_interaction",
    ,
    drop = FALSE
  ],
  "Discordant 37" = go_all[
    go_all$gene_set == "discordant",
    ,
    drop = FALSE
  ],
  "Protein specific 14" = go_all[
    go_all$gene_set == "protein_specific",
    ,
    drop = FALSE
  ],
  "mRNA specific 9" = go_all[
    go_all$gene_set == "mRNA_specific",
    ,
    drop = FALSE
  ],
  "Same direction 4" = go_all[
    go_all$gene_set == "same_direction",
    ,
    drop = FALSE
  ],
  "Significant terms" = go_env$significant_go,
  "Ribosome comparison" = go_comparison,
  "Gene sets" = go_gene_sets
)

########################################
# POWER WORKBOOK                       #
########################################

cross_index = match(figure3$target, crossstudy$target)
if (any(is.na(cross_index))) {
  stop("Could not align Figure 3 targets with power results.")
}

power_sensitivity = make_power_support_table(
  observed = figure3$observed_figure3_rediscovery,
  predicted = figure3$predicted_same_direction_lod4_5,
  scale08 = crossstudy$expected_detection_lod4_5_scale0_8[cross_index],
  scale06 = crossstudy$expected_detection_lod4_5_scale0_6[cross_index],
  benchmark = benchmark,
  simulation = simulation
)

benchmark_row = benchmark[1, , drop = FALSE]
benchmark_long = data.frame(
  measure = clean_name(names(benchmark_row)),
  value = vapply(
    benchmark_row,
    function(x) {
      value = x[[1]]
      if (is.na(value)) return(NA_character_)
      as.character(value)
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)
benchmark_long = clean_visible_terms(benchmark_long)

power_method_comparison = data.frame(
  pve = power_validation_env$engine_grid_pve,
  published_code_power = power_validation_env$engine_public,
  current_method_power = power_validation_env$engine_current,
  current_minus_published = power_validation_env$engine_diff,
  absolute_difference = abs(power_validation_env$engine_diff),
  stringsAsFactors = FALSE
)

power_sheets_final = list(
  "Study summary" = power_summary,
  "Sensitivity summary" = power_sensitivity,
  "Power by study" = critical,
  "Albert to other studies" = crossstudy,
  "Figure 3 by study" = figure3,
  "Figure 3 summary" = figure3_fit,
  "Albert benchmark" = benchmark_long,
  "Power method comparison" = power_method_comparison,
  "Simulation results" = simulation
)

########################################
# OTHER WORKBOOKS                      #
########################################

replicate_sheets_final = list(
  "Summary" = replicate_sheets$Summary,
  "Removed data" = replicate_sheets$Removed_items
)

gene_sheets_final = list(
  "Summary" = gene_sheets$Summary,
  "Gene data" = gene_sheets$Gene_level
)

dependence_sheets_final = list(
  "Summary" = dependence_sheets$Summary,
  "Locus results" = dependence_sheets$Locus_level,
  "Panel overlap" = dependence_sheets$Panel_overlap
)

########################################
# BUILD PLAIN WORKBOOKS                #
########################################

if (dir.exists(build_asys)) {
  unlink(build_asys, recursive = TRUE, force = TRUE)
}
ensure_directory(build_asys)

write_plain_workbook(
  replicate_sheets_final,
  "replicate_filtering.xlsx"
)
write_plain_workbook(
  gene_sheets_final,
  "supported_gene_characteristics.xlsx"
)
write_plain_workbook(
  go_sheets_final,
  "interaction_GO.xlsx"
)
write_plain_workbook(
  power_sheets_final,
  "power_analysis.xlsx"
)
write_plain_workbook(
  dependence_sheets_final,
  "segregant_dependence.xlsx"
)

########################################
# VALIDATE PRESENTATION                #
########################################

expected_sheets = list(
  replicate_filtering.xlsx = names(replicate_sheets_final),
  supported_gene_characteristics.xlsx = names(gene_sheets_final),
  interaction_GO.xlsx = names(go_sheets_final),
  power_analysis.xlsx = names(power_sheets_final),
  segregant_dependence.xlsx = names(dependence_sheets_final)
)

present = sort(list.files(build_asys, pattern = "\\.xlsx$"))
if (!identical(present, sort(final_files))) {
  stop("Final workbook build is incomplete.")
}

for (filename in final_files) {
  path = paste0(build_asys, filename)
  observed = openxlsx::getSheetNames(path)

  if (!identical(observed, expected_sheets[[filename]])) {
    stop("Unexpected sheets in final workbook: ", filename)
  }

  if (any(grepl("validation|engine|contrast", observed, ignore.case = TRUE))) {
    stop("Unwanted technical wording remains in sheet names: ", filename)
  }
}

simulation_check = openxlsx::read.xlsx(
  paste0(build_asys, "power_analysis.xlsx"),
  sheet = "Simulation results",
  check.names = FALSE
)
if (
  nrow(simulation_check) != 63L ||
  !identical(
    names(simulation_check)[1:5],
    c("n", "data_size", "rho", "lod", "n_rep")
  )
) {
  stop("Simulation-results sheet has unexpected structure.")
}

benchmark_check = openxlsx::read.xlsx(
  paste0(build_asys, "power_analysis.xlsx"),
  sheet = "Albert benchmark",
  check.names = FALSE
)
if (
  ncol(benchmark_check) != 2L ||
  !identical(names(benchmark_check), c("measure", "value"))
) {
  stop("Albert benchmark must be a two-column table.")
}

go_summary_check = openxlsx::read.xlsx(
  paste0(build_asys, "interaction_GO.xlsx"),
  sheet = "Summary",
  check.names = FALSE
)
comparison_col = "Discordant_vs_other_interaction_QTLs"
comparison_fdr_col = "Class_comparison_FDR"
if (
  !(comparison_col %in% names(go_summary_check)) ||
  !(comparison_fdr_col %in% names(go_summary_check)) ||
  any(is.na(go_summary_check[[comparison_col]])) ||
  any(is.na(go_summary_check[[comparison_fdr_col]]))
) {
  stop("GO summary does not contain complete class-comparison results.")
}

########################################
# PROMOTE                              #
########################################

if (dir.exists(final_asys)) {
  unlink(final_asys, recursive = TRUE, force = TRUE)
}
if (!file.rename(build_asys, final_asys)) {
  stop("Could not promote final reviewer workbooks.")
}

########################################
# MANIFEST                             #
########################################

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
    "Summary" = "GO terms highlighted in the response with class comparisons",
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
    "Albert benchmark" = "Albert et al. 2018 public power calculation and matched current calculation",
    "Power method comparison" = "Published-code and current analytic power across 0.5-10% PVE",
    "Simulation results" = "Direct phenotype simulations over regular sample-size, effect-size, and LOD grids"
  ),
  segregant_dependence.xlsx = c(
    "Summary" = "Shared-segregant sensitivity results",
    "Locus results" = "All 928 targeted original, segregant-clustered, and crossed-random-effect results",
    "Panel overlap" = "Shared-segregant counts for all 21 pairs of joint-mapping datasets"
  )
)

workbook_data = list(
  replicate_filtering.xlsx = replicate_sheets_final,
  supported_gene_characteristics.xlsx = gene_sheets_final,
  interaction_GO.xlsx = go_sheets_final,
  power_analysis.xlsx = power_sheets_final,
  segregant_dependence.xlsx = dependence_sheets_final
)

manifest_rows = list()
manifest_i = 1L
for (filename in final_files) {
  for (sheet_name in names(workbook_data[[filename]])) {
    dat = clean_visible_terms(
      workbook_data[[filename]][[sheet_name]]
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

if (nrow(manifest) != 25L) {
  stop("Expected exactly 25 final supporting-data sheets.")
}

write.table(
  manifest,
  file = paste0(review_asys, "review_output_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("\nFINAL REVIEWER WORKBOOKS READY\n")
cat("Workbooks: 5\n")
cat("Supporting-data sheets: 25\n")
cat("Simulation settings: 63 (3 sample sizes x 7 rho values x 3 LOD values)\n")
cat("Simulation replicates per setting: 1000\n")
cat("Simulation rho values: ", paste(simulation_rho, collapse = ", "), "\n", sep = "")
cat("Formatting: plain white sheets with bold headers and filters\n")
cat("Directory: ", final_asys, "\n", sep = "")