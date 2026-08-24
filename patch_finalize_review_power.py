#!/usr/bin/env python3
from pathlib import Path

main_path = Path("iReviewPowerEffectCurve.R")
val_path = Path("iReviewPowerEffectCurveValidate.R")

main = main_path.read_text()
val = val_path.read_text()

# -----------------------------------------------------------------------------
# Primary analysis: compact the cross-study output, add full Albert-engine
# agreement diagnostics, and expose optional fork-based parallelism.
# -----------------------------------------------------------------------------

old = 'albert_scales = c(1.00, 0.80, 0.60)\n'
new = '''albert_scales = c(1.00, 0.80, 0.60)\nengine_grid_pve = seq(0.005, 0.10, by = 0.001)\n\npower_cores = suppressWarnings(\n  as.integer(Sys.getenv("REVIEW_POWER_CORES", "1"))\n)\nif (is.na(power_cores) || power_cores < 1L) power_cores = 1L\n'''
if old not in main:
    raise SystemExit("STOP: primary albert_scales anchor not found")
main = main.replace(old, new, 1)

old = '''solve_r_for_power = function(n, lod, target_power) {\n'''
insert = '''parallel_lapply = function(X, FUN) {\n  if (\n    .Platform$OS.type == "unix" &&\n    power_cores > 1L &&\n    length(X) > 1L\n  ) {\n    return(parallel::mclapply(\n      X,\n      FUN,\n      mc.cores = min(power_cores, length(X)),\n      mc.preschedule = TRUE\n    ))\n  }\n  lapply(X, FUN)\n}\n\nsolve_r_for_power = function(n, lod, target_power) {\n'''
if old not in main:
    raise SystemExit("STOP: primary solve_r_for_power anchor not found")
main = main.replace(old, insert, 1)

old = '''message(\n  "Albert-reference effect scales: ",\n  paste(albert_scales, collapse = ", ")\n)\n'''
new = '''message(\n  "Albert-reference effect scales: ",\n  paste(albert_scales, collapse = ", ")\n)\nmessage("Parallel worker limit: ", power_cores)\n'''
if old not in main:
    raise SystemExit("STOP: primary message anchor not found")
main = main.replace(old, new, 1)

old = '''albert_reference_crossstudy = do.call(\n  rbind,\n  lapply(\n    albert_targets,\n'''
new = '''albert_reference_crossstudy = do.call(\n  rbind,\n  parallel_lapply(\n    albert_targets,\n'''
if old not in main:
    raise SystemExit("STOP: primary cross-study lapply anchor not found")
main = main.replace(old, new, 1)

anchor = 'rownames(albert_reference_crossstudy) = NULL\n\n################################################\n# TABLE 3: ALBERT ET AL. 2018 POWER BENCHMARK #\n################################################\n'
compact = '''rownames(albert_reference_crossstudy) = NULL\n\n# Keep the full target x scale x LOD grid in memory for validation, but emit\n# only the compact table needed for interpretation and reviewer-facing use.\nextract_crossstudy_power = function(effect_scale, lod, column) {\n  rows = albert_reference_crossstudy[\n    albert_reference_crossstudy$effect_scale == effect_scale &\n      albert_reference_crossstudy$lod == lod,\n    ,\n    drop = FALSE\n  ]\n  rows = rows[match(albert_targets, rows$target), , drop = FALSE]\n  rows[[column]]\n}\n\nbase_crossstudy = albert_reference_crossstudy[\n  albert_reference_crossstudy$effect_scale == 1 &\n    albert_reference_crossstudy$lod == 3.03,\n  c(\n    "target",\n    "target_analyte",\n    "target_medium",\n    "target_n",\n    "target_n_genes",\n    "n_albert_qtl_on_measured_genes",\n    "reference_mean_abs_r",\n    "reference_median_abs_r",\n    "reference_mean_pve",\n    "reference_median_pve"\n  ),\n  drop = FALSE\n]\nbase_crossstudy = base_crossstudy[\n  match(albert_targets, base_crossstudy$target),\n  ,\n  drop = FALSE\n]\n\nalbert_reference_crossstudy_table = base_crossstudy\nalbert_reference_crossstudy_table$expected_detection_lod2 =\n  extract_crossstudy_power(1, 2, "expected_any_detection")\nalbert_reference_crossstudy_table$expected_detection_lod3_03 =\n  extract_crossstudy_power(1, 3.03, "expected_any_detection")\nalbert_reference_crossstudy_table$expected_detection_lod4_5 =\n  extract_crossstudy_power(1, 4.5, "expected_any_detection")\nalbert_reference_crossstudy_table$expected_detection_lod7 =\n  extract_crossstudy_power(1, 7, "expected_any_detection")\nalbert_reference_crossstudy_table$expected_detection_lod4_5_scale0_8 =\n  extract_crossstudy_power(0.8, 4.5, "expected_any_detection")\nalbert_reference_crossstudy_table$expected_detection_lod4_5_scale0_6 =\n  extract_crossstudy_power(0.6, 4.5, "expected_any_detection")\nalbert_reference_crossstudy_table$expected_same_direction_lod4_5 =\n  extract_crossstudy_power(1, 4.5, "expected_same_direction_detection")\nrownames(albert_reference_crossstudy_table) = NULL\n\n################################################\n# TABLE 3: ALBERT ET AL. 2018 POWER BENCHMARK #\n################################################\n'''
if anchor not in main:
    raise SystemExit("STOP: primary compact-table anchor not found")
main = main.replace(anchor, compact, 1)

anchor = '''our_ab2018_power = correlation_detection_prob(\n  rho = ab2018_abs_r,\n  n = ab2018_total_n,\n  lod = ab2018_lod\n)\n\nalbert_2018_power_benchmark = data.frame(\n'''
engine = '''our_ab2018_power = correlation_detection_prob(\n  rho = ab2018_abs_r,\n  n = ab2018_total_n,\n  lod = ab2018_lod\n)\n\n# Compare the two analytic engines across a range of effect sizes, not only\n# at the published 2.5%-PVE benchmark. This is deterministic; stochastic\n# phenotype simulation remains an independent validator-only check.\nab2018_public_grid_power = vapply(\n  engine_grid_pve,\n  function(pve) {\n    delta = 2 * sqrt(pve / (1 - pve))\n    stats::power.t.test(\n      n = ab2018_n_per_group,\n      delta = delta,\n      sd = 1,\n      sig.level = ab2018_sig_level,\n      type = "two.sample",\n      alternative = "two.sided"\n    )$power\n  },\n  numeric(1)\n)\nab2018_current_grid_power = vapply(\n  sqrt(engine_grid_pve),\n  function(rho) {\n    correlation_detection_prob(\n      rho = rho,\n      n = ab2018_total_n,\n      lod = ab2018_lod\n    )$p_any_detection\n  },\n  numeric(1)\n)\nab2018_grid_difference =\n  ab2018_current_grid_power - ab2018_public_grid_power\n\nalbert_2018_power_benchmark = data.frame(\n'''
if anchor not in main:
    raise SystemExit("STOP: primary benchmark-grid anchor not found")
main = main.replace(anchor, engine, 1)

old = '''  current_minus_public_code =\n    our_ab2018_power$p_any_detection - ab2018_power_t_test,\n  agrees_with_published_gt90 =\n'''
new = '''  current_minus_public_code =\n    our_ab2018_power$p_any_detection - ab2018_power_t_test,\n  engine_grid_pve_min = min(engine_grid_pve),\n  engine_grid_pve_max = max(engine_grid_pve),\n  engine_grid_pve_step = unique(diff(engine_grid_pve))[[1]],\n  engine_grid_n = length(engine_grid_pve),\n  engine_grid_max_abs_difference = max(abs(ab2018_grid_difference)),\n  engine_grid_mean_abs_difference = mean(abs(ab2018_grid_difference)),\n  engine_grid_rmse = sqrt(mean(ab2018_grid_difference ^ 2)),\n  agrees_with_published_gt90 =\n'''
if old not in main:
    raise SystemExit("STOP: primary benchmark columns anchor not found")
main = main.replace(old, new, 1)

old = '''write_table(\n  albert_reference_crossstudy,\n  "albert_reference_crossstudy.tsv"\n)\n'''
new = '''write_table(\n  albert_reference_crossstudy_table,\n  "albert_reference_crossstudy.tsv"\n)\n'''
if old not in main:
    raise SystemExit("STOP: primary cross-study write anchor not found")
main = main.replace(old, new, 1)

# -----------------------------------------------------------------------------
# Validator: understand the compact cross-study table, validate the full public
# vs current engine curve, and parallelize the independent phenotype simulations.
# -----------------------------------------------------------------------------

old = '''if (is.na(n_validation) || n_validation < 1000L) {\n  stop("REVIEW_POWER_NVALID must be at least 1000.")\n}\n'''
new = '''if (is.na(n_validation) || n_validation < 1000L) {\n  stop("REVIEW_POWER_NVALID must be at least 1000.")\n}\n\npower_cores = suppressWarnings(\n  as.integer(Sys.getenv("REVIEW_POWER_CORES", "1"))\n)\nif (is.na(power_cores) || power_cores < 1L) power_cores = 1L\n'''
if old not in val:
    raise SystemExit("STOP: validator n_validation anchor not found")
val = val.replace(old, new, 1)

old = '''is_nonincreasing = function(x, tolerance = 1e-12) {\n  if (length(x) <= 1L) return(TRUE)\n  all(diff(x) <= tolerance)\n}\n'''
new = '''is_nonincreasing = function(x, tolerance = 1e-12) {\n  if (length(x) <= 1L) return(TRUE)\n  all(diff(x) <= tolerance)\n}\n\nparallel_lapply = function(X, FUN) {\n  if (\n    .Platform$OS.type == "unix" &&\n    power_cores > 1L &&\n    length(X) > 1L\n  ) {\n    return(parallel::mclapply(\n      X,\n      FUN,\n      mc.cores = min(power_cores, length(X)),\n      mc.preschedule = TRUE\n    ))\n  }\n  lapply(X, FUN)\n}\n'''
if old not in val:
    raise SystemExit("STOP: validator parallel helper anchor not found")
val = val.replace(old, new, 1)

# Replace the existing cross-study validator block wholesale.
start = val.index('########################################\n# ALBERT CROSS-STUDY TABLE             #\n########################################\n')
end = val.index('########################################\n# DIRECT PHENOTYPE SIMULATION          #\n########################################\n')
new_cross = '''########################################\n# ALBERT CROSS-STUDY TABLE             #\n########################################\n\nexpected_targets = setdiff(expected_studies, "Albert_2018")\nexpected_cross_columns = c(\n  "expected_detection_lod2",\n  "expected_detection_lod3_03",\n  "expected_detection_lod4_5",\n  "expected_detection_lod7",\n  "expected_detection_lod4_5_scale0_8",\n  "expected_detection_lod4_5_scale0_6",\n  "expected_same_direction_lod4_5"\n)\n\nchecks[[check_i]] = add_check(\n  "crossstudy_table_is_compact_one_row_per_target",\n  nrow(crossstudy) == length(expected_targets) &&\n    setequal(crossstudy$target, expected_targets) &&\n    all(expected_cross_columns %in% colnames(crossstudy)),\n  paste(\n    "rows =", nrow(crossstudy),\n    "; targets =", length(unique(crossstudy$target))\n  )\n)\ncheck_i = check_i + 1L\n\nstudies = load_single_rda(\n  paste0(imed, "studies.rda"),\n  object_name = "studies"\n)\nfrank_qtls_list = load_single_rda(\n  paste0(imed, "Frank/qtls_list_260204.rda"),\n  object_name = "qtls_list"\n)\n\nalbert = frank_qtls_list[["Albert_2018"]]\nif (is.null(albert)) {\n  stop("Albert_2018 missing from Frank QTL object.")\n}\nalbert_abs_r = abs(albert$cor_coef)\n\ncross_recomputed = parallel_lapply(\n  expected_targets,\n  function(target) {\n    target_genes = colnames(studies[[target]]$pheno)\n    keep = albert$gene_ID %in% target_genes\n    ref_r = albert_abs_r[keep]\n    row = crossstudy[crossstudy$target == target, , drop = FALSE]\n    if (nrow(row) != 1L) stop("Expected one compact row for ", target)\n\n    configs = data.frame(\n      effect_scale = c(1, 1, 1, 1, 0.8, 0.6, 1),\n      lod = c(2, 3.03, 4.5, 7, 4.5, 4.5, 4.5),\n      direction = c("any", "any", "any", "any", "any", "any", "same"),\n      output_col = expected_cross_columns,\n      stringsAsFactors = FALSE\n    )\n\n    expected = vapply(\n      seq_len(nrow(configs)),\n      function(i) {\n        scaled_r = pmin(0.999, ref_r * configs$effect_scale[[i]])\n        probs = vapply(\n          scaled_r,\n          function(rho) {\n            independent_probability(\n              rho = rho,\n              n = row$target_n[[1]],\n              lod = configs$lod[[i]]\n            )[[configs$direction[[i]]]]\n          },\n          numeric(1)\n        )\n        mean(probs)\n      },\n      numeric(1)\n    )\n    names(expected) = configs$output_col\n\n    c(\n      reference_error = max(\n        abs(row$n_albert_qtl_on_measured_genes[[1]] - length(ref_r)),\n        abs(row$reference_mean_abs_r[[1]] - mean(ref_r)),\n        abs(row$reference_median_abs_r[[1]] - stats::median(ref_r)),\n        abs(row$reference_mean_pve[[1]] - mean(ref_r ^ 2)),\n        abs(row$reference_median_pve[[1]] - stats::median(ref_r ^ 2))\n      ),\n      power_error = max(abs(as.numeric(row[1, expected_cross_columns]) - expected))\n    )\n  }\n)\ncross_recomputed = do.call(rbind, cross_recomputed)\n\nchecks[[check_i]] = add_check(\n  "crossstudy_reference_and_power_recompute_exactly",\n  max(cross_recomputed[, "reference_error"]) < 1e-12 &&\n    max(cross_recomputed[, "power_error"]) < 1e-12,\n  paste(\n    "maximum reference-summary error =",\n    signif(max(cross_recomputed[, "reference_error"]), 7),\n    "; maximum power error =",\n    signif(max(cross_recomputed[, "power_error"]), 7)\n  )\n)\ncheck_i = check_i + 1L\n\nlod_monotone = apply(\n  crossstudy[, c(\n    "expected_detection_lod2",\n    "expected_detection_lod3_03",\n    "expected_detection_lod4_5",\n    "expected_detection_lod7"\n  )],\n  1,\n  is_nonincreasing\n)\nscale_monotone = apply(\n  crossstudy[, c(\n    "expected_detection_lod4_5_scale0_6",\n    "expected_detection_lod4_5_scale0_8",\n    "expected_detection_lod4_5"\n  )],\n  1,\n  is_nondecreasing\n)\n\nchecks[[check_i]] = add_check(\n  "crossstudy_sensitivities_are_monotone",\n  all(lod_monotone) && all(scale_monotone),\n  paste(\n    "LOD-monotone targets =", sum(lod_monotone), "/", length(lod_monotone),\n    "; scale-monotone targets =", sum(scale_monotone), "/", length(scale_monotone)\n  )\n)\ncheck_i = check_i + 1L\n\n'''
val = val[:start] + new_cross + val[end:]

# Extend the Albert benchmark validator to the full effect-size curve summary.
anchor = '''ab_method_difference = abs(ab_current[["any"]] - ab_public_power)\nchecks[[check_i]] = add_check(\n  "albert_2018_and_current_power_engines_agree",\n  ab_method_difference < 0.002,\n  paste(\n    "absolute power difference at PVE 2.5% =",\n    signif(ab_method_difference, 7)\n  )\n)\ncheck_i = check_i + 1L\n\n'''
replacement = '''ab_method_difference = abs(ab_current[["any"]] - ab_public_power)\nchecks[[check_i]] = add_check(\n  "albert_2018_and_current_power_engines_agree",\n  ab_method_difference < 0.002,\n  paste(\n    "absolute power difference at PVE 2.5% =",\n    signif(ab_method_difference, 7)\n  )\n)\ncheck_i = check_i + 1L\n\nengine_grid_pve = seq(0.005, 0.10, by = 0.001)\nengine_public = vapply(\n  engine_grid_pve,\n  function(pve) {\n    delta = 2 * sqrt(pve / (1 - pve))\n    stats::power.t.test(\n      n = ab_n_per_group,\n      delta = delta,\n      sd = 1,\n      sig.level = ab_alpha,\n      type = "two.sample",\n      alternative = "two.sided"\n    )$power\n  },\n  numeric(1)\n)\nengine_current = vapply(\n  sqrt(engine_grid_pve),\n  function(rho) {\n    independent_probability(\n      rho = rho,\n      n = ab_total_n,\n      lod = ab_lod\n    )[["any"]]\n  },\n  numeric(1)\n)\nengine_diff = engine_current - engine_public\nengine_summary_error = max(\n  abs(benchmark$engine_grid_pve_min - min(engine_grid_pve)),\n  abs(benchmark$engine_grid_pve_max - max(engine_grid_pve)),\n  abs(benchmark$engine_grid_pve_step - unique(diff(engine_grid_pve))[[1]]),\n  abs(benchmark$engine_grid_n - length(engine_grid_pve)),\n  abs(benchmark$engine_grid_max_abs_difference - max(abs(engine_diff))),\n  abs(benchmark$engine_grid_mean_abs_difference - mean(abs(engine_diff))),\n  abs(benchmark$engine_grid_rmse - sqrt(mean(engine_diff ^ 2)))\n)\nchecks[[check_i]] = add_check(\n  "albert_2018_engines_agree_across_effect_grid",\n  engine_summary_error < 1e-12 && max(abs(engine_diff)) < 0.01,\n  paste(\n    "PVE grid = 0.005 to 0.10; maximum absolute difference =",\n    signif(max(abs(engine_diff)), 7),\n    "; table-summary error =", signif(engine_summary_error, 7)\n  )\n)\ncheck_i = check_i + 1L\n\n'''
if anchor not in val:
    raise SystemExit("STOP: validator engine-grid anchor not found")
val = val.replace(anchor, replacement, 1)

old = '''simulation = do.call(\n  rbind,\n  lapply(\n    seq_len(nrow(validation_grid)),\n'''
new = '''simulation = do.call(\n  rbind,\n  parallel_lapply(\n    seq_len(nrow(validation_grid)),\n'''
if old not in val:
    raise SystemExit("STOP: validator simulation lapply anchor not found")
val = val.replace(old, new, 1)

# Apply the already-tested self-referential style-check fix if the local copy
# still has the older form. If it is already fixed, leave it unchanged.
old = '''checks[[check_i]] = add_check(\n  "review_power_scripts_follow_equals_style",\n  !grepl("<-", analysis_text, fixed = TRUE) &&\n    !grepl("<-", validator_text, fixed = TRUE),\n  "No '<-' assignments detected in primary or validator script."\n)\n'''
new = '''assignment_token = paste0("<", "-")\n\nchecks[[check_i]] = add_check(\n  "review_power_scripts_follow_equals_style",\n  !grepl(assignment_token, analysis_text, fixed = TRUE) &&\n    !grepl(assignment_token, validator_text, fixed = TRUE),\n  "No left-arrow assignments detected in primary or validator script."\n)\n'''
if old in val:
    val = val.replace(old, new, 1)

main_path.write_text(main)
val_path.write_text(val)

print("Final reviewer-power patch applied.")
print("  - compact 9-row Albert cross-study table")
print("  - full Albert-2018 analytic engine agreement summary")
print("  - REVIEW_POWER_CORES fork parallelism")
print("  - updated independent validator")
