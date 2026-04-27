# prelims
source("/PHShome/kv1026/p05/remap_functions.R")

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
# idir = stable input location for study objects and mash inputs
# odir = stable output location for mash-generated arm files
# ------------------------------------------------------------------------------

idir = imed
odir = paste0(imed, "mashr/")
ncores = 48
# ncores = max(1L, parallel::detectCores())

# ------------------------------------------------------------------
# IMPORTANT DELINEATION
# ------------------------------------------------------------------
# metadata_support_template:
#   only used to define the analyte-supported gene universe
#   (the 385 genes supported in >=2 analyte datasets)
#
# metadata_qtl:
#   built fresh from the full study set actually used in a mash arm
#   and enriched with per-study segregant IDs after any arm-specific
#   harmonization / downsampling
#
# In other words:
#   support genes come from scaMetadata
#   mash columns come from full studyMetadata / studies.rda
# ------------------------------------------------------------------
metadata_support_template = scaMetadata

# mashr generation settings
min_present = 4L
max_covariance_rows = 20000L
max_mixture_rows = 200000L
alpha = 1
overwrite = FALSE

# covariance-row selection settings
signal_fraction = 0.40
signal_quantile = 0.80

# LD-pruning settings for mash training rows
ld_abs_cor_threshold = 0.8

# choose any combination of arms to run
prefix = c("fs", "ds", "na")

# analysis-arm settings
analysis_arms = list(
  ds = list(
    label = "downsampled",
    prefix = "ds",
    downsample = TRUE,
    target_study = "AlbertAlt_2018",
    target_n = 100L,
    seed = 343L,
    exclude_studies = character(0)
  ),
  fs = list(
    label = "fullscan",
    prefix = "fs",
    downsample = FALSE,
    seed = 454L,
    exclude_studies = character(0)
  ),
  na = list(
    label = "no_albert",
    prefix = "na",
    downsample = FALSE,
    seed = 565L,
    exclude_studies = "AlbertAlt_2018"
  )
)

valid_prefixes = names(analysis_arms)
if (length(prefix) < 1L) {
  stop("prefix must contain at least one of: ", paste(valid_prefixes, collapse = ", "))
}
if (anyDuplicated(prefix)) {
  stop("prefix contains duplicates: ", paste(unique(prefix[duplicated(prefix)]), collapse = ", "))
}
if (!all(prefix %in% valid_prefixes)) {
  stop(
    "Invalid prefix value(s): ",
    paste(setdiff(prefix, valid_prefixes), collapse = ", "),
    ". Choose any combination of: ",
    paste(valid_prefixes, collapse = ", ")
  )
}

analysis_arms = analysis_arms[prefix]

#######################################
# HELPERS FOR MASHR OBJECT GENERATION #
#######################################

load_or_compute_rda = function(path, object_name, expr, overwrite = FALSE, message_prefix = NULL) {
  if (file.exists(path) && !isTRUE(overwrite)) {
    if (!is.null(message_prefix)) {
      message(message_prefix, "Loading existing ", object_name, " from: ", path)
    } else {
      message("Loading existing ", object_name, " from: ", path)
    }
    return(load_rda_object(path, object_name = object_name))
  }

  if (!is.null(message_prefix)) {
    message(message_prefix, "Generating ", object_name, " ...")
  } else {
    message("Generating ", object_name, " ...")
  }

  object = force(expr)
  assign(object_name, object)
  save(list = object_name, file = path)
  object
}

load_or_compute_rds = function(path, expr, overwrite = FALSE, label = NULL) {
  if (file.exists(path) && !isTRUE(overwrite)) {
    if (is.null(label)) {
      message("Loading existing RDS from: ", path)
    } else {
      message("Loading existing ", label, " from: ", path)
    }
    return(readRDS(path))
  }

  if (is.null(label)) {
    message("Generating object for: ", path)
  } else {
    message("Generating ", label, " ...")
  }

  object = force(expr)
  saveRDS(object, path)
  object
}

write_lines_if_needed = function(x, path, overwrite = FALSE) {
  if (!file.exists(path) || isTRUE(overwrite)) {
    writeLines(x, path)
  }
  invisible(path)
}

infile = function(filename) file.path(idir, filename)
outfile = function(filename) file.path(odir, filename)
prefixed_path = function(prefix, filename) outfile(paste0(prefix, "_", filename))

is_albert_study = function(x) x %in% c("Albert_2018", "AlbertAlt_2018")

canonicalize_seg_ids = function(x) {
  x = as.character(x)
  x = trimws(x)
  x = gsub("_", ".", x)
  x = gsub("\\s+", "", x)
  x
}

get_primary_seg_ids = function(study) {
  if (!is.null(study$pheno) && !is.null(rownames(study$pheno))) return(canonicalize_seg_ids(rownames(study$pheno)))
  if (!is.null(study$pclean) && !is.null(rownames(study$pclean))) return(canonicalize_seg_ids(rownames(study$pclean)))
  if (!is.null(study$gdata) && !is.null(rownames(study$gdata))) return(canonicalize_seg_ids(rownames(study$gdata)))
  stop("Could not determine primary segregant IDs for study.")
}

rekey_one_matrix = function(mat, seg_ids, mat_name, study_name) {
  if (is.null(mat)) return(NULL)

  if (nrow(mat) != length(seg_ids)) {
    rn = rownames(mat)
    if (!is.null(rn)) {
      rn_c = canonicalize_seg_ids(rn)
      idx = match(seg_ids, rn_c)
      if (all(!is.na(idx)) && length(unique(idx)) == length(idx)) {
        mat = mat[idx, , drop = FALSE]
        rownames(mat) = seg_ids
        return(mat)
      }
    }
    stop("Matrix ", mat_name, " in study ", study_name, " could not be re-keyed to the target segregant IDs.")
  }

  rn = rownames(mat)
  if (!is.null(rn)) {
    rn_c = canonicalize_seg_ids(rn)
    idx = match(seg_ids, rn_c)
    if (all(!is.na(idx)) && length(unique(idx)) == length(idx)) {
      mat = mat[idx, , drop = FALSE]
      rownames(mat) = seg_ids
      return(mat)
    }
  }

  rownames(mat) = seg_ids
  mat
}

harmonize_study_rowkeys = function(study, study_name = NA_character_) {
  seg_ids = get_primary_seg_ids(study)

  if (length(seg_ids) == 0L) stop("Study ", study_name, " has zero segregants after key extraction.")
  if (anyDuplicated(seg_ids)) {
    dup = unique(seg_ids[duplicated(seg_ids)])
    stop("Study ", study_name, " has duplicated canonical segregant IDs: ", paste(head(dup, 10L), collapse = ", "))
  }

  if (!is.null(study$gdata)) study$gdata = rekey_one_matrix(study$gdata, seg_ids, "gdata", study_name)
  if (!is.null(study$gscaled)) study$gscaled = rekey_one_matrix(study$gscaled, seg_ids, "gscaled", study_name)
  if (!is.null(study$pclean)) study$pclean = rekey_one_matrix(study$pclean, seg_ids, "pclean", study_name)
  if (!is.null(study$pheno)) study$pheno = rekey_one_matrix(study$pheno, seg_ids, "pheno", study_name)

  attr(study, "seg_ids") = seg_ids
  study
}

subset_study_rows = function(study, keep_rows) {
  seg_ids = attr(study, "seg_ids")
  if (is.null(seg_ids)) seg_ids = get_primary_seg_ids(study)

  keep_rows = canonicalize_seg_ids(keep_rows)
  keep_rows = intersect(keep_rows, seg_ids)

  if (length(keep_rows) == 0L) stop("Downsampling produced zero retained rows after enforcing segregant keys.")

  if (!is.null(study$gdata)) study$gdata = study$gdata[keep_rows, , drop = FALSE]
  if (!is.null(study$gscaled)) study$gscaled = study$gscaled[keep_rows, , drop = FALSE]
  if (!is.null(study$pclean)) study$pclean = study$pclean[keep_rows, , drop = FALSE]
  if (!is.null(study$pheno)) study$pheno = study$pheno[keep_rows, , drop = FALSE]

  attr(study, "seg_ids") = keep_rows
  study
}

downsample_study_object = function(study, target_n, seed = 1L, study_name = NA_character_) {
  seg_ids = attr(study, "seg_ids")
  if (is.null(seg_ids)) {
    seg_ids = get_primary_seg_ids(study)
    attr(study, "seg_ids") = seg_ids
  }

  n0 = length(seg_ids)
  message("Downsampling keyed row universe for ", study_name, ": ", n0, " segregants.")

  if (n0 == 0L) stop("Study ", study_name, " has zero keyed segregants after harmonization.")

  if (n0 <= target_n) {
    study = subset_study_rows(study, seg_ids)
    if (length(attr(study, "seg_ids")) != n0) stop("Failed to preserve keyed segregants for study ", study_name, ".")
    return(study)
  }

  set.seed(seed)
  keep_rows = sort(sample(seg_ids, size = target_n, replace = FALSE))
  study = subset_study_rows(study, keep_rows)

  if (length(attr(study, "seg_ids")) != target_n) {
    stop("Downsampling failed for study ", study_name, ": expected ", target_n, " keyed segregants, got ", length(attr(study, "seg_ids")), ".")
  }

  study
}

build_metadata_qtl_enriched = function(studies) {
  study_names = names(studies)
  md = studyMetadata[studyMetadata$study %in% study_names, , drop = FALSE]
  md = md[match(study_names, md$study), , drop = FALSE]
  if (any(is.na(md$study))) stop("Some mash study names were not found in studyMetadata.")
  md$segs = lapply(study_names, function(x) attr(studies[[x]], "seg_ids"))
  md$n_segs = vapply(md$segs, length, integer(1))
  md
}

prepare_qtl_studies = function(arm_cfg) {
  message("Loading full QTL studies and scaling genotypes...")

  studies = load_rda_object(infile("studies.rda"), "studies")
  if ("Albert_2018" %in% names(studies)) studies = studies[names(studies) != "Albert_2018"]

  exclude_studies = arm_cfg$exclude_studies
  if (is.null(exclude_studies)) exclude_studies = character(0)

  if (length(exclude_studies) > 0L) {
    message("Excluding studies for ", arm_cfg$label, ": ", paste(exclude_studies, collapse = ", "))
    studies = studies[setdiff(names(studies), exclude_studies)]
  }

  if (length(studies) == 0L) stop("No studies remain after exclusion.")

  studies = lapply(names(studies), function(x) harmonize_study_rowkeys(studies[[x]], study_name = x)) |> `names<-`(names(studies))

  if (isTRUE(arm_cfg$downsample)) {
    target_study = arm_cfg$target_study
    if (!(target_study %in% names(studies))) stop("Target study for downsampling not found in studies: ", target_study)

    message("Downsampling ", target_study, " to ", arm_cfg$target_n, " segregants for ", arm_cfg$label, ".")
    studies[[target_study]] = downsample_study_object(
      studies[[target_study]],
      target_n = arm_cfg$target_n,
      seed = arm_cfg$seed,
      study_name = target_study
    )
  }

  studies = lapply(studies, function(y) {
    y$gscaled = scale(y$gdata)
    attr(y, "seg_ids") = rownames(y$gdata)
    y
  })

  studies
}

prepare_support_studies = function(metadata_support, arm_cfg) {
  message("Loading support-only study objects...")

  exclude_studies = arm_cfg$exclude_studies
  if (is.null(exclude_studies)) exclude_studies = character(0)
  keep_studies = setdiff(metadata_support$study, exclude_studies)

  out = vector("list", length(keep_studies))
  names(out) = keep_studies

  pb = txtProgressBar(min = 0, max = length(keep_studies), style = 3)
  for (i in seq_along(keep_studies)) {
    x = keep_studies[i]
    env = new.env(parent = emptyenv())
    load(infile(paste0(x, ".rda")), envir = env)
    out[[x]] = harmonize_study_rowkeys(env[[x]], study_name = x)
    rm(env)

    if (isTRUE(arm_cfg$downsample) && identical(x, arm_cfg$target_study)) {
      out[[x]] = downsample_study_object(
        out[[x]],
        target_n = arm_cfg$target_n,
        seed = arm_cfg$seed,
        study_name = x
      )
    }

    setTxtProgressBar(pb, i)
    if (i %% 5L == 0L) gc()
  }
  close(pb)

  out
}

build_support_metadata = function(metadata_support, support_studies) {
  keep_studies = names(support_studies)

  md = metadata_support
  md = md[md$study %in% keep_studies, , drop = FALSE]
  md = md[match(keep_studies, md$study), , drop = FALSE]

  md$genes = vector("list", length(md$study))
  md$segs = vector("list", length(md$study))
  md$n_segs = integer(length(md$study))

  pb = txtProgressBar(min = 0, max = length(md$study), style = 3)
  for (i in seq_along(md$study)) {
    x = md$study[i]
    md$genes[[i]] = colnames(support_studies[[x]]$pheno)
    md$segs[[i]] = attr(support_studies[[x]], "seg_ids")
    md$n_segs[i] = length(md$segs[[i]])
    setTxtProgressBar(pb, i)
  }
  close(pb)

  md
}

subsample_training_rows = function(Bhat, Shat, max_rows, signal_fraction = 0.40, signal_quantile = 0.80) {
  null_fraction = 1 - signal_fraction
  n = nrow(Bhat)

  if (n <= max_rows) {
    return(list(
      Bhat = Bhat,
      Shat = Shat,
      idx = seq_len(n),
      row_score = apply(abs(Bhat / Shat), 1, max, na.rm = TRUE),
      threshold = NA_real_
    ))
  }

  Z = Bhat / Shat
  row_score = apply(abs(Z), 1, max, na.rm = TRUE)
  row_score[!is.finite(row_score)] = 0

  thr = as.numeric(stats::quantile(row_score, probs = signal_quantile, na.rm = TRUE, names = FALSE))
  signal_idx = which(row_score >= thr)
  null_idx = which(row_score < thr)

  n_null = round(max_rows * null_fraction)
  n_signal = max_rows - n_null

  if (length(null_idx) < n_null) {
    deficit = n_null - length(null_idx)
    n_null = length(null_idx)
    n_signal = min(length(signal_idx), n_signal + deficit)
  }

  if (length(signal_idx) < n_signal) {
    deficit = n_signal - length(signal_idx)
    n_signal = length(signal_idx)
    n_null = min(length(null_idx), n_null + deficit)
  }

  idx_null = if (n_null > 0L) sample(null_idx, n_null, replace = FALSE) else integer(0)
  idx_signal = if (n_signal > 0L) sample(signal_idx, n_signal, replace = FALSE) else integer(0)
  idx = sort(c(idx_null, idx_signal))

  list(
    Bhat = Bhat[idx, , drop = FALSE],
    Shat = Shat[idx, , drop = FALSE],
    idx = idx,
    row_score = row_score,
    threshold = thr
  )
}

choose_ld_reference_study = function(studies) {
  preferred = intersect(c("GrossbachPRT_2022", "GrossbachRNA_2022"), names(studies))

  if (length(preferred) > 0L) {
    n_by_study = vapply(preferred, function(x) nrow(studies[[x]]$gdata), integer(1))
    return(preferred[which.max(n_by_study)])
  }

  fallback = setdiff(names(studies), c("Albert_2018", "AlbertAlt_2018"))
  if (length(fallback) == 0L) stop("Could not find a non-Albert study to use as the LD reference.")

  n_by_study = vapply(fallback, function(x) nrow(studies[[x]]$gdata), integer(1))
  fallback[which.max(n_by_study)]
}

global_ld_prune_markers = function(markers, ld_submat, abs_cor_threshold = 0.8) {
  markers = as.character(markers)
  if (length(markers) == 0L) return(character(0))

  keep = character(0)
  blocked = setNames(rep(FALSE, length(markers)), markers)

  for (marker_i in markers) {
    if (blocked[[marker_i]]) next
    keep = c(keep, marker_i)
    linked = abs(ld_submat[marker_i, markers]) >= abs_cor_threshold
    linked[!is.finite(linked)] = FALSE
    blocked[markers[linked]] = TRUE
  }

  keep
}

###########################
# MASHR OBJECT GENERATION #
###########################

run_mash_generation = function(arm_name, arm_cfg, ld_reference) {

  prefix = arm_cfg$prefix
  odir_arm = file.path(odir, paste0(prefix, "_training"))
  arm_file = function(filename) file.path(odir_arm, paste0(prefix, "_", filename))
  arm_overwrite = overwrite || isTRUE(arm_cfg$downsample)

  message("============================================================")
  message("Running mashr arm: ", arm_name, " [", arm_cfg$label, "]")
  if (isTRUE(arm_cfg$downsample)) {
    message("Downsampled arm detected: forcing recomputation of study-dependent cached objects.")
  }
  message("============================================================")

  prepare_output_dirs(odir_arm)

  studies = prepare_qtl_studies(arm_cfg)
  markers = colnames(studies[[1]]$gdata)
  study_names = names(studies)
  ld_submat = align_ld_matrix(ld_reference, markers)
  metadata_qtl = build_metadata_qtl_enriched(studies)

  if (isTRUE(arm_cfg$downsample) && arm_cfg$target_study %in% metadata_qtl$study) {
    n_target = metadata_qtl$n_segs[match(arm_cfg$target_study, metadata_qtl$study)]
    if (!identical(as.integer(n_target), as.integer(arm_cfg$target_n))) {
      stop("Downsampled arm ", prefix, " expected ", arm_cfg$target_n, " segregants for ", arm_cfg$target_study, " but metadata_qtl reports ", n_target, ".")
    }
  }

  write_tsv_if_needed(metadata_qtl, arm_file("metadata_qtl.tsv"), overwrite = TRUE)
  saveRDS(metadata_qtl, arm_file("metadata_qtl.rds"))

  support_studies = load_or_compute_rds(
    arm_file("studies_analyte.rds"),
    expr = prepare_support_studies(metadata_support_template, arm_cfg),
    overwrite = arm_overwrite,
    label = paste0(prefix, "_support_studies [", arm_name, "]")
  )

  metadata_support = load_or_compute_rds(
    arm_file("metadata_analyte_enriched.rds"),
    expr = {
      message("Collecting support-only gene and segregant metadata...")
      build_support_metadata(metadata_support_template, support_studies)
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_metadata_support_enriched [", arm_name, "]")
  )

  supported_analyte = load_or_compute_rds(
    arm_file("supported_analyte.rds"),
    expr = {
      message("Computing supported analytes from support-only metadata...")
      suppTwo(metadata_support)
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_supported_analyte [", arm_name, "]")
  )

  write_lines_if_needed(supported_analyte, arm_file("supported_analyte.txt"), overwrite = TRUE)
  message("Analyte-supported genes [", arm_name, "]: ", length(supported_analyte))

  genes_any_qtl = unique(unlist(lapply(studies, function(x) colnames(x$pclean))))
  run_genes = intersect(supported_analyte, genes_any_qtl)

  dropped_support_genes = setdiff(supported_analyte, genes_any_qtl)
  if (length(dropped_support_genes) > 0L) {
    warning(
      "Some support-defined genes were absent from all full QTL studies in this arm: ",
      paste(head(dropped_support_genes, 10L), collapse = ", "),
      if (length(dropped_support_genes) > 10L) " ..."
    )
  }

  message("Genes to run after intersecting support genes with full QTL studies [", arm_name, "]: ", length(run_genes))
  saveRDS(run_genes, arm_file("run_genes.rds"))
  write_lines_if_needed(run_genes, arm_file("run_genes.txt"), overwrite = TRUE)

  gene_fits = load_or_compute_rda(
    arm_file("gene_fits_supported_analyte.rda"),
    object_name = "gene_fits",
    expr = {
      message("Fitting all supported genes across full QTL studies...")
      gf = parallel::mclapply(
        X = run_genes,
        FUN = function(g) {
          out = lapply(studies, fit_gene_in_study, gene = g)
          names(out) = study_names
          out
        },
        mc.cores = ncores
      )
      names(gf) = run_genes
      gf
    },
    overwrite = arm_overwrite,
    message_prefix = paste0("[", prefix, "] ")
  )

  gene_present = load_or_compute_rds(
    arm_file("gene_present_supported_analyte.rds"),
    expr = {
      message("Summarizing per-gene dataset presence across full QTL studies...")
      gp = integer(length(run_genes))
      names(gp) = run_genes
      pb = txtProgressBar(min = 0, max = length(run_genes), style = 3)
      for (i in seq_along(run_genes)) {
        g = run_genes[i]
        fits = gene_fits[[g]]
        gp[i] = sum(vapply(fits, function(x) any(is.finite(x$beta) & is.finite(x$se)), logical(1)))
        setTxtProgressBar(pb, i)
      }
      close(pb)
      gp
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_gene_present [", arm_name, "]")
  )

  albert_survival = load_or_compute_rds(
    arm_file("albert_survival_diagnostic.rds"),
    expr = {
      target = "AlbertAlt_2018"

      if (!(target %in% study_names)) {
        out = list(
          study = target,
          present_in_arm = FALSE,
          gene_summary = data.frame(),
          overall = data.frame(
            study = target,
            present_in_arm = FALSE,
            n_run_genes = length(run_genes),
            n_genes_with_any_finite = NA_integer_,
            frac_genes_with_any_finite = NA_real_,
            mean_fraction_marker_rows_finite = NA_real_,
            median_fraction_marker_rows_finite = NA_real_,
            stringsAsFactors = FALSE
          )
        )
      } else {
        gene_summary = lapply(run_genes, function(g) {
          fit = gene_fits[[g]][[target]]
          beta_ok = is.finite(fit$beta)
          se_ok = is.finite(fit$se)
          row_ok = beta_ok & se_ok
          data.frame(
            gene = g,
            n_marker_rows = length(fit$beta),
            n_marker_rows_finite = sum(row_ok),
            fraction_marker_rows_finite = mean(row_ok),
            any_finite = any(row_ok),
            n = fit$n,
            stringsAsFactors = FALSE
          )
        }) |> do.call(what = rbind)

        out = list(
          study = target,
          present_in_arm = TRUE,
          gene_summary = gene_summary,
          overall = data.frame(
            study = target,
            present_in_arm = TRUE,
            n_run_genes = nrow(gene_summary),
            n_genes_with_any_finite = sum(gene_summary$any_finite),
            frac_genes_with_any_finite = mean(gene_summary$any_finite),
            mean_fraction_marker_rows_finite = mean(gene_summary$fraction_marker_rows_finite),
            median_fraction_marker_rows_finite = stats::median(gene_summary$fraction_marker_rows_finite),
            stringsAsFactors = FALSE
          )
        )
      }

      out
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_albert_survival_diagnostic [", arm_name, "]")
  )

  if (isTRUE(albert_survival$present_in_arm)) {
    write_tsv_if_needed(albert_survival$gene_summary, arm_file("albert_survival_by_gene.tsv"), overwrite = TRUE)
    write_tsv_if_needed(albert_survival$overall, arm_file("albert_survival_summary.tsv"), overwrite = TRUE)

    msg = albert_survival$overall
    message(
      "[", prefix, "] AlbertAlt survival: ",
      msg$n_genes_with_any_finite, "/", msg$n_run_genes,
      " genes with any finite data (", round(100 * msg$frac_genes_with_any_finite, 1), "%); ",
      "mean finite marker-row fraction = ", round(msg$mean_fraction_marker_rows_finite, 3),
      ", median = ", round(msg$median_fraction_marker_rows_finite, 3)
    )

    if (is.na(msg$n_genes_with_any_finite) || msg$n_genes_with_any_finite == 0L) {
      stop("[", prefix, "] AlbertAlt_2018 is included in this arm but has zero genes with any finite beta+SE after gene fitting. Aborting mash training.")
    }
  } else {
    message("[", prefix, "] AlbertAlt_2018 is not present in this arm.")
  }

  global_ld_pruned_markers = load_or_compute_rds(
    arm_file("global_ld_pruned_markers.rds"),
    expr = {
      message("Computing one global LD-pruned marker panel for this arm...")
      out = global_ld_prune_markers(markers = markers, ld_submat = ld_submat, abs_cor_threshold = ld_abs_cor_threshold)
      message("Global LD-pruned marker panel: ", length(out), " of ", length(markers), " markers retained.")
      out
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_global_ld_pruned_markers [", arm_name, "]")
  )

  write_lines_if_needed(global_ld_pruned_markers, arm_file("global_ld_pruned_markers.txt"), overwrite = TRUE)

  global_marker_inventory = data.frame(
    arm = arm_name,
    prefix = prefix,
    n_markers_before_global_ld_prune = length(markers),
    n_markers_after_global_ld_prune = length(global_ld_pruned_markers),
    frac_markers_retained_after_global_ld_prune = length(global_ld_pruned_markers) / length(markers),
    ld_abs_cor_threshold = ld_abs_cor_threshold,
    stringsAsFactors = FALSE
  )

  write_tsv_if_needed(global_marker_inventory, arm_file("global_marker_inventory.tsv"), overwrite = TRUE)
  write_lines_if_needed(
    c(
      paste0("arm\t", global_marker_inventory$arm),
      paste0("prefix\t", global_marker_inventory$prefix),
      paste0("n_markers_before_global_ld_prune\t", global_marker_inventory$n_markers_before_global_ld_prune),
      paste0("n_markers_after_global_ld_prune\t", global_marker_inventory$n_markers_after_global_ld_prune),
      paste0("frac_markers_retained_after_global_ld_prune\t", round(global_marker_inventory$frac_markers_retained_after_global_ld_prune, 6)),
      paste0("ld_abs_cor_threshold\t", global_marker_inventory$ld_abs_cor_threshold)
    ),
    arm_file("global_marker_inventory.txt"),
    overwrite = TRUE
  )

  training_rows = load_or_compute_rds(
    arm_file("training_rows_supported_analyte_global_ldpruned.rds"),
    expr = {
      message("Building pooled training rows using one global LD-pruned marker panel...")
      tr = vector("list", length(run_genes))
      ps = vector("list", length(run_genes))
      marker_catalog_before = vector("list", length(run_genes))
      marker_catalog_after = vector("list", length(run_genes))

      pb = txtProgressBar(min = 0, max = length(run_genes), style = 3)
      for (i in seq_along(run_genes)) {
        gene = run_genes[i]
        mats = build_gene_mats(gene, gene_fits, study_names, markers)
        keep_present = rowSums(is.finite(mats$Bhat) & is.finite(mats$Shat)) >= min_present

        if (any(keep_present)) {
          B0 = mats$Bhat[keep_present, , drop = FALSE]
          S0 = mats$Shat[keep_present, , drop = FALSE]
          gene_markers = rownames(B0)
          keep_global = gene_markers %in% global_ld_pruned_markers

          tr[[i]] = list(Bhat = B0[keep_global, , drop = FALSE], Shat = S0[keep_global, , drop = FALSE])
          ps[[i]] = data.frame(gene = gene, n_rows_before = nrow(B0), n_rows_after = sum(keep_global), stringsAsFactors = FALSE)
          marker_catalog_before[[i]] = gene_markers
          marker_catalog_after[[i]] = gene_markers[keep_global]

          rm(B0, S0, gene_markers, keep_global)
        }

        setTxtProgressBar(pb, i)
        if (i %% 50L == 0L) gc()
      }
      close(pb)

      tr = tr[!vapply(tr, is.null, logical(1))]
      ps = ps[!vapply(ps, is.null, logical(1))]
      marker_catalog_before = marker_catalog_before[!vapply(marker_catalog_before, is.null, logical(1))]
      marker_catalog_after = marker_catalog_after[!vapply(marker_catalog_after, is.null, logical(1))]

      if (length(tr) == 0L) stop("No training rows were collected after global LD pruning.")

      list(
        blocks = tr,
        prune_stats = do.call(rbind, ps),
        marker_catalog_before = marker_catalog_before,
        marker_catalog_after = marker_catalog_after
      )
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_training_rows_global_ldpruned [", arm_name, "]")
  )

  prune_stats = training_rows$prune_stats
  write_tsv_if_needed(prune_stats, arm_file("training_row_prune_stats.tsv"), overwrite = TRUE)

  training_row_marker_ids_before_ld_prune = unlist(training_rows$marker_catalog_before, use.names = FALSE)
  training_row_marker_ids_after_ld_prune = unlist(training_rows$marker_catalog_after, use.names = FALSE)

  training_row_inventory = data.frame(
    arm = arm_name,
    prefix = prefix,
    n_training_blocks = length(training_rows$blocks),
    n_rows_before_global_ld_prune = sum(prune_stats$n_rows_before),
    n_rows_after_global_ld_prune = sum(prune_stats$n_rows_after),
    n_unique_markers_before_global_ld_prune = length(unique(training_row_marker_ids_before_ld_prune)),
    n_unique_markers_after_global_ld_prune = length(unique(training_row_marker_ids_after_ld_prune)),
    mean_rows_before_global_ld_prune_per_gene = mean(prune_stats$n_rows_before),
    mean_rows_after_global_ld_prune_per_gene = mean(prune_stats$n_rows_after),
    median_rows_before_global_ld_prune_per_gene = median(prune_stats$n_rows_before),
    median_rows_after_global_ld_prune_per_gene = median(prune_stats$n_rows_after),
    stringsAsFactors = FALSE
  )

  write_tsv_if_needed(training_row_inventory, arm_file("training_row_inventory.tsv"), overwrite = TRUE)
  write_lines_if_needed(
    c(
      paste0("arm\t", training_row_inventory$arm),
      paste0("prefix\t", training_row_inventory$prefix),
      paste0("n_training_blocks\t", training_row_inventory$n_training_blocks),
      paste0("n_rows_before_global_ld_prune\t", training_row_inventory$n_rows_before_global_ld_prune),
      paste0("n_rows_after_global_ld_prune\t", training_row_inventory$n_rows_after_global_ld_prune),
      paste0("n_unique_markers_before_global_ld_prune\t", training_row_inventory$n_unique_markers_before_global_ld_prune),
      paste0("n_unique_markers_after_global_ld_prune\t", training_row_inventory$n_unique_markers_after_global_ld_prune),
      paste0("mean_rows_before_global_ld_prune_per_gene\t", round(training_row_inventory$mean_rows_before_global_ld_prune_per_gene, 3)),
      paste0("mean_rows_after_global_ld_prune_per_gene\t", round(training_row_inventory$mean_rows_after_global_ld_prune_per_gene, 3)),
      paste0("median_rows_before_global_ld_prune_per_gene\t", round(training_row_inventory$median_rows_before_global_ld_prune_per_gene, 3)),
      paste0("median_rows_after_global_ld_prune_per_gene\t", round(training_row_inventory$median_rows_after_global_ld_prune_per_gene, 3))
    ),
    arm_file("training_row_inventory.txt"),
    overwrite = TRUE
  )

  message("Collected globally LD-pruned training blocks from ", length(training_rows$blocks), " genes [", arm_name, "].")
  message("Mean rows per gene before global prune: ", round(mean(prune_stats$n_rows_before), 2), " | after global prune: ", round(mean(prune_stats$n_rows_after), 2))

  training_mats_file = arm_file("mash_training_mats_supported_analyte_global_ldpruned.rda")
  if (file.exists(training_mats_file) && !isTRUE(arm_overwrite)) {
    message("Loading existing training matrices from: ", training_mats_file)
    env = new.env(parent = emptyenv())
    loaded_names = load(training_mats_file, envir = env)

    if ("training_mats" %in% loaded_names) {
      training_mats = env$training_mats
    } else if (all(c("cov_Bhat", "cov_Shat", "cov_ss", "mix_Bhat", "mix_Shat", "mix_idx", "n_rows_total") %in% loaded_names)) {
      training_mats = list(
        cov_Bhat = env$cov_Bhat,
        cov_Shat = env$cov_Shat,
        cov_ss = env$cov_ss,
        mix_Bhat = env$mix_Bhat,
        mix_Shat = env$mix_Shat,
        mix_idx = env$mix_idx,
        n_rows_total = env$n_rows_total
      )
    } else {
      stop(
        "Existing training mats file does not contain either 'training_mats' or all of ",
        "'cov_Bhat', 'cov_Shat', 'cov_ss', 'mix_Bhat', 'mix_Shat', 'mix_idx', and 'n_rows_total': ",
        training_mats_file
      )
    }

    rm(env, loaded_names)
    gc()
  } else {
    message("Row-binding pooled globally LD-pruned training matrices...")
    all_Bhat = do.call(rbind, lapply(training_rows$blocks, `[[`, "Bhat"))
    all_Shat = do.call(rbind, lapply(training_rows$blocks, `[[`, "Shat"))
    rownames(all_Bhat) = NULL
    rownames(all_Shat) = NULL

    n_rows_total = nrow(all_Bhat)
    message("Training matrix dimensions before any sampling [", arm_name, "]: ", n_rows_total, " rows x ", ncol(all_Bhat), " studies")

    message("Sampling covariance-learning rows to ", max_covariance_rows, " with enrichment on max |Z|...")
    set.seed(arm_cfg$seed)
    cov_ss = subsample_training_rows(
      all_Bhat,
      all_Shat,
      max_rows = max_covariance_rows,
      signal_fraction = signal_fraction,
      signal_quantile = signal_quantile
    )
    cov_Bhat = cov_ss$Bhat
    cov_Shat = cov_ss$Shat

    message("Covariance signal threshold (", signal_quantile * 100, "th pct of max |Z|): ", round(cov_ss$threshold, 3))
    message("Covariance rows retained: ", nrow(cov_Bhat), " [target null-ish fraction=", 1 - signal_fraction, ", signal-ish fraction=", signal_fraction, "]")

    set.seed(arm_cfg$seed + 1000L)
    mix_idx = if (n_rows_total > max_mixture_rows) sort(sample.int(n_rows_total, max_mixture_rows, replace = FALSE)) else seq_len(n_rows_total)
    mix_Bhat = all_Bhat[mix_idx, , drop = FALSE]
    mix_Shat = all_Shat[mix_idx, , drop = FALSE]

    message("Mixture-training rows retained: ", nrow(mix_Bhat), " of ", n_rows_total, " (uniform random sample after global LD pruning)")

    training_mats = list(
      cov_Bhat = cov_Bhat,
      cov_Shat = cov_Shat,
      cov_ss = cov_ss,
      mix_Bhat = mix_Bhat,
      mix_Shat = mix_Shat,
      mix_idx = mix_idx,
      n_rows_total = n_rows_total
    )

    save(training_mats, file = training_mats_file)
    rm(all_Bhat, all_Shat, cov_Bhat, cov_Shat, mix_Bhat, mix_Shat, mix_idx, n_rows_total)
    gc()
  }

  rm(training_rows)
  gc()

  cov_Bhat = training_mats$cov_Bhat
  cov_Shat = training_mats$cov_Shat
  mix_Bhat = training_mats$mix_Bhat
  mix_Shat = training_mats$mix_Shat

  n_rows_available_for_sampling = length(training_row_marker_ids_after_ld_prune)
  n_rows_available_for_covariance = training_mats$n_rows_total
  n_rows_available_for_mixture = training_mats$n_rows_total
  n_rows_used_for_covariance = nrow(cov_Bhat)
  n_rows_used_for_mixture = nrow(mix_Bhat)

  cov_sampled_marker_ids = if (!is.null(training_mats$cov_ss$idx)) training_row_marker_ids_after_ld_prune[training_mats$cov_ss$idx] else training_row_marker_ids_after_ld_prune
  mix_sampled_marker_ids = if (!is.null(training_mats$mix_idx)) training_row_marker_ids_after_ld_prune[training_mats$mix_idx] else training_row_marker_ids_after_ld_prune

  training_sampling_inventory = data.frame(
    arm = arm_name,
    prefix = prefix,
    max_covariance_rows = max_covariance_rows,
    max_mixture_rows = max_mixture_rows,
    n_rows_available_for_sampling = n_rows_available_for_sampling,
    n_rows_available_for_covariance = n_rows_available_for_covariance,
    n_rows_available_for_mixture = n_rows_available_for_mixture,
    n_rows_used_for_covariance = n_rows_used_for_covariance,
    n_rows_used_for_mixture = n_rows_used_for_mixture,
    n_rows_dropped_by_covariance_sampling = n_rows_available_for_covariance - n_rows_used_for_covariance,
    n_rows_dropped_by_mixture_sampling = n_rows_available_for_mixture - n_rows_used_for_mixture,
    frac_rows_retained_for_covariance = n_rows_used_for_covariance / n_rows_available_for_covariance,
    frac_rows_retained_for_mixture = n_rows_used_for_mixture / n_rows_available_for_mixture,
    n_unique_markers_available_for_sampling = length(unique(training_row_marker_ids_after_ld_prune)),
    n_unique_markers_used_for_covariance = length(unique(cov_sampled_marker_ids)),
    n_unique_markers_used_for_mixture = length(unique(mix_sampled_marker_ids)),
    n_studies_in_training_matrix = ncol(cov_Bhat),
    signal_fraction_target = signal_fraction,
    signal_quantile = signal_quantile,
    covariance_signal_threshold = if (!is.null(training_mats$cov_ss$threshold)) training_mats$cov_ss$threshold else NA_real_,
    stringsAsFactors = FALSE
  )

  write_tsv_if_needed(training_sampling_inventory, arm_file("training_sampling_inventory.tsv"), overwrite = TRUE)
  write_lines_if_needed(
    c(
      paste0("arm\t", training_sampling_inventory$arm),
      paste0("prefix\t", training_sampling_inventory$prefix),
      paste0("max_covariance_rows\t", training_sampling_inventory$max_covariance_rows),
      paste0("max_mixture_rows\t", training_sampling_inventory$max_mixture_rows),
      paste0("n_rows_available_for_sampling\t", training_sampling_inventory$n_rows_available_for_sampling),
      paste0("n_rows_available_for_covariance\t", training_sampling_inventory$n_rows_available_for_covariance),
      paste0("n_rows_available_for_mixture\t", training_sampling_inventory$n_rows_available_for_mixture),
      paste0("n_rows_used_for_covariance\t", training_sampling_inventory$n_rows_used_for_covariance),
      paste0("n_rows_used_for_mixture\t", training_sampling_inventory$n_rows_used_for_mixture),
      paste0("n_rows_dropped_by_covariance_sampling\t", training_sampling_inventory$n_rows_dropped_by_covariance_sampling),
      paste0("n_rows_dropped_by_mixture_sampling\t", training_sampling_inventory$n_rows_dropped_by_mixture_sampling),
      paste0("frac_rows_retained_for_covariance\t", round(training_sampling_inventory$frac_rows_retained_for_covariance, 6)),
      paste0("frac_rows_retained_for_mixture\t", round(training_sampling_inventory$frac_rows_retained_for_mixture, 6)),
      paste0("n_unique_markers_available_for_sampling\t", training_sampling_inventory$n_unique_markers_available_for_sampling),
      paste0("n_unique_markers_used_for_covariance\t", training_sampling_inventory$n_unique_markers_used_for_covariance),
      paste0("n_unique_markers_used_for_mixture\t", training_sampling_inventory$n_unique_markers_used_for_mixture),
      paste0("n_studies_in_training_matrix\t", training_sampling_inventory$n_studies_in_training_matrix),
      paste0("signal_fraction_target\t", training_sampling_inventory$signal_fraction_target),
      paste0("signal_quantile\t", training_sampling_inventory$signal_quantile),
      paste0("covariance_signal_threshold\t", round(training_sampling_inventory$covariance_signal_threshold, 6))
    ),
    arm_file("training_sampling_inventory.txt"),
    overwrite = TRUE
  )

  if ("AlbertAlt_2018" %in% study_names) {
    albert_cov_ok = is.finite(cov_Bhat[, "AlbertAlt_2018"]) & is.finite(cov_Shat[, "AlbertAlt_2018"])
    albert_mix_ok = is.finite(mix_Bhat[, "AlbertAlt_2018"]) & is.finite(mix_Shat[, "AlbertAlt_2018"])

    albert_training_summary = data.frame(
      study = "AlbertAlt_2018",
      n_covariance_rows = nrow(cov_Bhat),
      n_covariance_rows_finite = sum(albert_cov_ok),
      fraction_covariance_rows_finite = mean(albert_cov_ok),
      n_mixture_rows = nrow(mix_Bhat),
      n_mixture_rows_finite = sum(albert_mix_ok),
      fraction_mixture_rows_finite = mean(albert_mix_ok),
      stringsAsFactors = FALSE
    )

    write_tsv_if_needed(albert_training_summary, arm_file("albert_training_row_summary.tsv"), overwrite = TRUE)

    message(
      "[", prefix, "] AlbertAlt covariance-row survival: ",
      albert_training_summary$n_covariance_rows_finite, "/",
      albert_training_summary$n_covariance_rows,
      " rows finite (", round(100 * albert_training_summary$fraction_covariance_rows_finite, 1), "%); ",
      "mixture-row survival: ",
      albert_training_summary$n_mixture_rows_finite, "/",
      albert_training_summary$n_mixture_rows,
      " rows finite (", round(100 * albert_training_summary$fraction_mixture_rows_finite, 1), "%)"
    )

    if (albert_training_summary$n_covariance_rows_finite == 0L) {
      stop("[", prefix, "] AlbertAlt_2018 is included in this arm but has zero finite rows in the covariance-learning matrices. Aborting mash training.")
    }
    if (albert_training_summary$n_mixture_rows_finite == 0L) {
      stop("[", prefix, "] AlbertAlt_2018 is included in this arm but has zero finite rows in the mixture-training matrices. Aborting mash training.")
    }
  }

  Vhat = load_or_compute_rds(
    arm_file("null_correlation.rds"),
    expr = {
      message("Using identity null correlation Vhat.")
      diag(length(study_names))
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_Vhat [", arm_name, "]")
  )

  write_tsv_if_needed(as.data.frame(Vhat), arm_file("null_correlation.tsv"), overwrite = TRUE)

  cov_data0 = load_or_compute_rds(
    arm_file("cov_data0.rds"),
    expr = {
      message("Constructing mash covariance-learning data...")
      mash_set_data(cov_Bhat, cov_Shat, alpha = alpha)
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_cov_data0 [", arm_name, "]")
  )

  cov_data = load_or_compute_rds(
    arm_file("cov_data.rds"),
    expr = mash_update_data(cov_data0, V = Vhat),
    overwrite = arm_overwrite,
    label = paste0(prefix, "_cov_data [", arm_name, "]")
  )

  mix_data0 = load_or_compute_rds(
    arm_file("mix_data0.rds"),
    expr = {
      message("Constructing mash mixture-training data...")
      mash_set_data(mix_Bhat, mix_Shat, alpha = alpha)
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_mix_data0 [", arm_name, "]")
  )

  mix_data = load_or_compute_rds(
    arm_file("mix_data.rds"),
    expr = mash_update_data(mix_data0, V = Vhat),
    overwrite = arm_overwrite,
    label = paste0(prefix, "_mix_data [", arm_name, "]")
  )

  npc = load_or_compute_rds(
    arm_file("npc.rds"),
    expr = min(5L, length(study_names)),
    overwrite = arm_overwrite,
    label = paste0(prefix, "_npc [", arm_name, "]")
  )

  message("Using npc = ", npc, " [", arm_name, "]")

  U_c = load_or_compute_rds(
    arm_file("cov_canonical.rds"),
    expr = {
      message("Computing canonical covariance candidates...")
      t0 = Sys.time()
      out = cov_canonical(cov_data)
      message("Done canonical in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
      out
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_U_c [", arm_name, "]")
  )

  U_pca = load_or_compute_rds(
    arm_file("cov_pca.rds"),
    expr = {
      message("Computing PCA covariance candidates...")
      t0 = Sys.time()
      out = cov_pca(cov_data, npc = npc)
      message("Done PCA in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
      out
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_U_pca [", arm_name, "]")
  )

  U_ed = load_or_compute_rds(
    arm_file("cov_ed.rds"),
    expr = {
      message("Computing ED covariance candidates...")
      t0 = Sys.time()
      out = cov_ed(cov_data, U_pca)
      message("Done ED in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
      out
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_U_ed [", arm_name, "]")
  )

  Ulist = load_or_compute_rds(
    arm_file("covariance_candidates.rds"),
    expr = {
      out = c(U_c, U_pca, U_ed)
      message("Total covariance candidates: ", length(out))
      out
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_Ulist [", arm_name, "]")
  )

  message("Total covariance candidates: ", length(Ulist), " [", arm_name, "]")

  global_mash = load_or_compute_rds(
    arm_file("global_mash_model.rds"),
    expr = {
      message("Fitting global mash model on representative pooled rows...")
      t0 = Sys.time()
      out = mash(mix_data, Ulist, outputlevel = 1)
      message("Done global mash in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
      out
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_global_mash [", arm_name, "]")
  )

  mix_df = load_or_compute_rds(
    arm_file("mixture_weights.rds"),
    expr = {
      data.frame(
        covariance = names(get_estimated_pi(global_mash)),
        pi = as.numeric(get_estimated_pi(global_mash)),
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    },
    overwrite = arm_overwrite,
    label = paste0(prefix, "_mixture_weights [", arm_name, "]")
  )

  write_tsv_if_needed(mix_df, arm_file("mixture_weights.tsv"), overwrite = TRUE)

  summary_info = data.frame(
    arm = arm_name,
    prefix = prefix,
    label = arm_cfg$label,
    n_qtl_studies = length(study_names),
    n_support_studies = nrow(metadata_support),
    n_supported_analyte_genes = length(supported_analyte),
    n_run_genes = length(run_genes),
    ld_reference_study = ld_reference$study,
    ld_abs_cor_threshold = ld_abs_cor_threshold,
    n_markers_before_global_ld_prune = length(markers),
    n_markers_after_global_ld_prune = length(global_ld_pruned_markers),
    n_rows_before_global_ld_prune = sum(prune_stats$n_rows_before),
    n_rows_after_global_ld_prune = sum(prune_stats$n_rows_after),
    n_rows_used_for_covariance = training_sampling_inventory$n_rows_used_for_covariance,
    n_rows_used_for_mixture = training_sampling_inventory$n_rows_used_for_mixture,
    n_unique_markers_used_for_covariance = training_sampling_inventory$n_unique_markers_used_for_covariance,
    n_unique_markers_used_for_mixture = training_sampling_inventory$n_unique_markers_used_for_mixture,
    stringsAsFactors = FALSE
  )

  write_tsv_if_needed(summary_info, arm_file("analysis_summary.tsv"), overwrite = TRUE)

  out = list(
    arm_name = arm_name,
    prefix = prefix,
    analysis_dir = odir_arm,
    study_names = study_names,
    n_run_genes = length(run_genes),
    ld_reference_study = ld_reference$study,
    mixture_weights = mix_df,
    summary_info = summary_info,
    files = list(
      metadata_qtl = arm_file("metadata_qtl.rds"),
      metadata_support = arm_file("metadata_analyte_enriched.rds"),
      supported_analyte = arm_file("supported_analyte.rds"),
      run_genes = arm_file("run_genes.rds"),
      gene_fits = arm_file("gene_fits_supported_analyte.rda"),
      gene_present = arm_file("gene_present_supported_analyte.rds"),
      albert_survival = arm_file("albert_survival_diagnostic.rds"),
      albert_survival_gene_tsv = arm_file("albert_survival_by_gene.tsv"),
      albert_survival_summary_tsv = arm_file("albert_survival_summary.tsv"),
      albert_training_row_summary_tsv = arm_file("albert_training_row_summary.tsv"),
      global_ld_pruned_markers_rds = arm_file("global_ld_pruned_markers.rds"),
      global_ld_pruned_markers_txt = arm_file("global_ld_pruned_markers.txt"),
      global_marker_inventory_tsv = arm_file("global_marker_inventory.tsv"),
      global_marker_inventory_txt = arm_file("global_marker_inventory.txt"),
      training_rows = arm_file("training_rows_supported_analyte_global_ldpruned.rds"),
      training_row_inventory_tsv = arm_file("training_row_inventory.tsv"),
      training_row_inventory_txt = arm_file("training_row_inventory.txt"),
      training_sampling_inventory_tsv = arm_file("training_sampling_inventory.tsv"),
      training_sampling_inventory_txt = arm_file("training_sampling_inventory.txt"),
      training_mats = arm_file("mash_training_mats_supported_analyte_global_ldpruned.rda"),
      cov_data = arm_file("cov_data.rds"),
      mix_data = arm_file("mix_data.rds"),
      covariance_candidates = arm_file("covariance_candidates.rds"),
      global_mash = arm_file("global_mash_model.rds"),
      mixture_weights = arm_file("mixture_weights.rds")
    )
  )

  rm(
    studies,
    support_studies,
    metadata_qtl,
    metadata_support,
    supported_analyte,
    genes_any_qtl,
    run_genes,
    gene_fits,
    gene_present,
    albert_survival,
    ld_submat,
    global_ld_pruned_markers,
    global_marker_inventory,
    training_mats,
    cov_Bhat,
    cov_Shat,
    mix_Bhat,
    mix_Shat,
    Vhat,
    cov_data0,
    cov_data,
    mix_data0,
    mix_data,
    U_c,
    U_pca,
    U_ed,
    Ulist,
    global_mash,
    prune_stats,
    summary_info,
    mix_df,
    training_row_inventory,
    training_sampling_inventory,
    training_row_marker_ids_before_ld_prune,
    training_row_marker_ids_after_ld_prune,
    cov_sampled_marker_ids,
    mix_sampled_marker_ids
  )
  gc()

  out
}

################
# RUN ALL ARMS #
################

prepare_output_dirs(odir)

message("Running mash generation for arm(s): ", paste(names(analysis_arms), collapse = ", "))

shared_ld_reference = build_shared_ld_reference(overwrite = overwrite)

analysis_results = list()
for (arm_name in names(analysis_arms)) {
  analysis_results[[arm_name]] = run_mash_generation(
    arm_name = arm_name,
    arm_cfg = analysis_arms[[arm_name]],
    ld_reference = shared_ld_reference
  )
  gc()
}

rm(shared_ld_reference)
gc()

Message("Done.")