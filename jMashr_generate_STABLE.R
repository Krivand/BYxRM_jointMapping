# prelims
source("/PHShome/kv1026/p05/remap_functions.R")

odir = paste0(imed, "mashr/")
ncores = max(1L, parallel::detectCores())

# ------------------------------------------------------------------
# IMPORTANT DELINEATION
# ------------------------------------------------------------------
# metadata_support_template:
#   only used to define the analyte-supported gene universe
#   (the 385 genes supported in >=2 analyte datasets)
#
# metadata_qtl:
#   always built from the full study set actually used in a mash arm
#   and therefore may include AlbertAlt_2018
#
# In other words:
#   support genes come from scaMetadata
#   mash columns come from full studyMetadata / studies.rda
# ------------------------------------------------------------------
metadata_support_template = scaMetadata

# mashr generation settings
min_present = 1L
max_training_rows = 10000L
alpha = 1
overwrite = FALSE

# training-row selection settings
null_fraction = 0.70
signal_fraction = 0.30
signal_quantile = 0.80

# LD-pruning settings for mash training rows
# We prune rows before pooling so that one LD block does not flood the training set
# with many near-duplicate marker rows.
ld_abs_cor_threshold = 0.8
max_ld_representatives_per_gene = Inf

# analysis-arm settings
analysis_arms = list(
  fs = list(
    label = "fullscan",
    prefix = "fs",
    downsample = FALSE,
    seed = 1L,
    exclude_studies = character(0)
  ),
  ds = list(
    label = "downsampled",
    prefix = "ds",
    downsample = TRUE,
    target_study = "AlbertAlt_2018",
    target_n = 100L,
    seed = 1L,
    exclude_studies = character(0)
  ),
  na = list(
    label = "no_albert",
    prefix = "na",
    downsample = FALSE,
    seed = 1L,
    exclude_studies = "AlbertAlt_2018"
  )
)

#######################################
# HELPERS FOR MASHR OBJECT GENERATION #
#######################################

load_rda_object = function(path, object_name = NULL) {
  env = new.env(parent = emptyenv())
  nm = load(path, envir = env)

  if (!is.null(object_name)) {
    if (!(object_name %in% nm)) {
      stop("Object '", object_name, "' not found in ", path)
    }
    return(env[[object_name]])
  }

  if (length(nm) != 1L) {
    stop(
      "load_rda_object() expected exactly one object in ",
      path,
      " when object_name is NULL, but found: ",
      paste(nm, collapse = ", ")
    )
  }

  env[[nm]]
}

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

write_tsv_if_needed = function(dt, path, overwrite = FALSE, sep = "\t") {
  if (!file.exists(path) || isTRUE(overwrite)) {
    data.table::fwrite(dt, path, sep = sep)
  }
  invisible(path)
}

write_lines_if_needed = function(x, path, overwrite = FALSE) {
  if (!file.exists(path) || isTRUE(overwrite)) {
    writeLines(x, path)
  }
  invisible(path)
}

prepare_output_dirs = function(analysis_dir) {
  dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
}

subset_study_rows = function(study, keep_rows) {
  keep_rows = intersect(keep_rows, rownames(study$gdata))

  if (length(keep_rows) == 0L) {
    stop("Downsampling produced zero retained rows.")
  }

  if (!is.null(study$gdata)) {
    study$gdata = study$gdata[keep_rows, , drop = FALSE]
  }
  if (!is.null(study$gscaled)) {
    study$gscaled = study$gscaled[keep_rows, , drop = FALSE]
  }
  if (!is.null(study$pclean)) {
    keep_pclean = intersect(keep_rows, rownames(study$pclean))
    study$pclean = study$pclean[keep_pclean, , drop = FALSE]
  }
  if (!is.null(study$pheno)) {
    keep_pheno = intersect(keep_rows, rownames(study$pheno))
    study$pheno = study$pheno[keep_pheno, , drop = FALSE]
  }

  study
}

downsample_study_object = function(study, target_n, seed = 1L) {
  n0 = nrow(study$gdata)
  if (is.null(n0) || n0 <= target_n) {
    return(study)
  }

  set.seed(seed)
  keep_rows = sort(sample(rownames(study$gdata), size = target_n, replace = FALSE))
  subset_study_rows(study, keep_rows)
}

# This helper handles the Albert naming inconsistency in one place.
is_albert_study = function(x) {
  x %in% c("Albert_2018", "AlbertAlt_2018")
}

# ------------------------------------------------------------------
# FULL QTL STUDIES FOR THE ARM
# ------------------------------------------------------------------
# These are the studies that define the columns of mash.
# This can include AlbertAlt_2018.
# ------------------------------------------------------------------
prepare_qtl_studies = function(arm_cfg) {
  message("Loading full QTL studies and scaling genotypes...")

  studies_file = paste0(imed, "studies.rda")
  studies = load_rda_object(studies_file, "studies")

  # Drop legacy Albert_2018 if present.
  if ("Albert_2018" %in% names(studies)) {
    studies = studies[names(studies) != "Albert_2018"]
  }

  exclude_studies = arm_cfg$exclude_studies
  if (is.null(exclude_studies)) {
    exclude_studies = character(0)
  }

  if (length(exclude_studies) > 0L) {
    message("Excluding studies for ", arm_cfg$label, ": ", paste(exclude_studies, collapse = ", "))
    studies = studies[setdiff(names(studies), exclude_studies)]
  }

  if (length(studies) == 0L) {
    stop("No studies remain after exclusion.")
  }

  if (isTRUE(arm_cfg$downsample)) {
    target_study = arm_cfg$target_study
    if (!(target_study %in% names(studies))) {
      stop("Target study for downsampling not found in studies: ", target_study)
    }

    message(
      "Downsampling ", target_study, " to ", arm_cfg$target_n,
      " segregants for ", arm_cfg$label, "."
    )

    studies[[target_study]] = downsample_study_object(
      studies[[target_study]],
      target_n = arm_cfg$target_n,
      seed = arm_cfg$seed
    )
  }

  studies = lapply(studies, function(y) {
    y$gscaled = scale(y$gdata)
    y
  })

  studies
}

# ------------------------------------------------------------------
# SUPPORT STUDIES ONLY
# ------------------------------------------------------------------
# These are used only to define the analyte-supported gene set.
# They intentionally come from scaMetadata and therefore do NOT need
# to include AlbertAlt_2018.
# ------------------------------------------------------------------
prepare_support_studies = function(metadata_support, arm_cfg) {
  message("Loading support-only study objects...")

  exclude_studies = arm_cfg$exclude_studies
  if (is.null(exclude_studies)) {
    exclude_studies = character(0)
  }

  keep_studies = setdiff(metadata_support$study, exclude_studies)

  out = vector("list", length(keep_studies))
  names(out) = keep_studies

  pb = txtProgressBar(min = 0, max = length(keep_studies), style = 3)
  for (i in seq_along(keep_studies)) {
    x = keep_studies[i]
    env = new.env(parent = emptyenv())
    load(paste0(imed, x, ".rda"), envir = env)
    out[[x]] = env[[x]]
    rm(env)

    # Downsampling only matters if a support study is ever targeted.
    if (isTRUE(arm_cfg$downsample) && identical(x, arm_cfg$target_study)) {
      out[[x]] = downsample_study_object(
        out[[x]],
        target_n = arm_cfg$target_n,
        seed = arm_cfg$seed
      )
    }

    setTxtProgressBar(pb, i)
    if (i %% 5L == 0L) gc()
  }
  close(pb)

  out
}

# Build enriched support metadata for suppTwo().
# This metadata is ONLY for support-gene discovery.
build_support_metadata = function(metadata_support, support_studies) {
  keep_studies = names(support_studies)

  md = metadata_support
  md = md[md$study %in% keep_studies, , drop = FALSE]
  md = md[match(keep_studies, md$study), , drop = FALSE]

  md$genes = vector("list", length(md$study))
  md$segs = vector("list", length(md$study))

  pb = txtProgressBar(min = 0, max = length(md$study), style = 3)
  for (i in seq_along(md$study)) {
    x = md$study[i]
    md$genes[[i]] = colnames(support_studies[[x]]$pheno)
    md$segs[[i]] = rownames(support_studies[[x]]$pheno)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  md
}

subsample_training_rows = function(
  Bhat,
  Shat,
  max_rows,
  null_fraction = 0.70,
  signal_fraction = 0.30,
  signal_quantile = 0.80
) {
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

  thr = as.numeric(stats::quantile(
    row_score,
    probs = signal_quantile,
    na.rm = TRUE,
    names = FALSE
  ))

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

fit_gene_in_study = function(study, gene) {
  if (!(gene %in% colnames(study$pclean))) {
    m = ncol(study$gscaled)
    return(list(
      beta = rep(NA_real_, m),
      se = rep(NA_real_, m),
      n = NA_integer_
    ))
  }

  y = as.numeric(study$pclean[, gene])
  G = study$gscaled

  ok = is.finite(y)
  y = y[ok]
  G = G[ok, , drop = FALSE]

  n = length(y)
  m = ncol(G)

  if (n < 3L) {
    return(list(
      beta = rep(NA_real_, m),
      se = rep(NA_real_, m),
      n = n
    ))
  }

  numer = as.numeric(crossprod(y, G))
  denom = colSums(G^2)
  beta = numer / denom

  ssy = sum(y^2)
  rss = ssy - (numer^2 / denom)
  se = sqrt((rss / (n - 2)) / denom)

  beta[!is.finite(beta)] = NA_real_
  se[!is.finite(se)] = NA_real_

  list(beta = beta, se = se, n = n)
}

build_gene_mats = function(gene, gene_fits, study_names, markers) {
  fits = gene_fits[[gene]]

  Bhat = do.call(cbind, lapply(fits, `[[`, "beta"))
  Shat = do.call(cbind, lapply(fits, `[[`, "se"))

  colnames(Bhat) = study_names
  colnames(Shat) = study_names
  rownames(Bhat) = markers
  rownames(Shat) = markers

  ns = vapply(fits, `[[`, integer(1), "n")

  list(Bhat = Bhat, Shat = Shat, n = ns)
}

prefixed_path = function(prefix, filename) {
  file.path(odir, paste0(prefix, "_", filename))
}

# Pick one reference study for LD pruning.
# Prefer whichever Grossbach study has more segregants.
# If neither is present, fall back to the non-Albert study with the most segregants.
choose_ld_reference_study = function(studies) {
  preferred = intersect(c("GrossbachPRT_2022", "GrossbachRNA_2022"), names(studies))

  if (length(preferred) > 0L) {
    n_by_study = vapply(
      preferred,
      function(x) nrow(studies[[x]]$gdata),
      integer(1)
    )
    return(preferred[which.max(n_by_study)])
  }

  fallback = setdiff(names(studies), c("Albert_2018", "AlbertAlt_2018"))
  if (length(fallback) == 0L) {
    stop("Could not find a non-Albert study to use as the LD reference.")
  }

  n_by_study = vapply(
    fallback,
    function(x) nrow(studies[[x]]$gdata),
    integer(1)
  )
  fallback[which.max(n_by_study)]
}

# Build the shared LD matrix once before all arms.
# We do not do this per arm because the marker map is the same and both fs and ds
# should use the exact same pruning rule.
build_shared_ld_reference = function(overwrite = FALSE) {
  path = file.path(odir, "shared_ld_reference.rds")

  load_or_compute_rds(
    path,
    expr = {
      message("Loading studies once to compute shared LD reference...")
      studies_file = paste0(imed, "studies.rda")
      studies_all = load_rda_object(studies_file, "studies")

      if ("Albert_2018" %in% names(studies_all)) {
        studies_all = studies_all[names(studies_all) != "Albert_2018"]
      }

      ref_study = choose_ld_reference_study(studies_all)
      message("Using ", ref_study, " as the LD reference study.")

      G = studies_all[[ref_study]]$gdata
      if (is.null(G)) {
        stop("Reference study ", ref_study, " does not contain gdata.")
      }

      ld_r = stats::cor(G, use = "pairwise.complete.obs")

      out = list(
        study = ref_study,
        markers = colnames(G),
        ld_r = ld_r
      )

      rm(studies_all, G, ld_r)
      gc()

      out
    },
    overwrite = overwrite,
    label = "shared LD reference"
  )
}

# Reorder and subset the shared LD matrix to a requested marker order.
# This makes sure gene-level pruning lines up with the marker rows in Bhat/Shat.
align_ld_matrix = function(ld_reference, markers) {
  if (!all(markers %in% ld_reference$markers)) {
    missing_markers = setdiff(markers, ld_reference$markers)
    stop(
      "Some markers are missing from the shared LD matrix: ",
      paste(head(missing_markers, 10L), collapse = ", ")
    )
  }

  ld_reference$ld_r[markers, markers, drop = FALSE]
}

# Greedy LD pruning for one gene's marker rows.
# The purpose is simple: keep one representative row from each tightly linked block,
# using the strongest row in that block as the representative.
ld_prune_gene_rows = function(
  Bhat,
  Shat,
  ld_submat,
  abs_cor_threshold = 0.95,
  max_keep = Inf
) {
  stopifnot(identical(rownames(Bhat), rownames(Shat)))
  stopifnot(identical(rownames(Bhat), rownames(ld_submat)))
  stopifnot(identical(colnames(ld_submat), rownames(ld_submat)))

  Z = Bhat / Shat
  row_score = apply(abs(Z), 1, max, na.rm = TRUE)
  row_score[!is.finite(row_score)] = -Inf

  ord = order(row_score, decreasing = TRUE)
  marker_names = rownames(Bhat)
  keep = integer(0)
  blocked = rep(FALSE, length(marker_names))
  names(blocked) = marker_names

  for (idx in ord) {
    marker_i = marker_names[idx]

    if (isTRUE(blocked[marker_i])) {
      next
    }

    keep = c(keep, idx)

    if (length(keep) >= max_keep) {
      break
    }

    linked = abs(ld_submat[marker_i, marker_names]) >= abs_cor_threshold
    linked[!is.finite(linked)] = FALSE
    blocked[marker_names[linked]] = TRUE
  }

  sort(unique(keep))
}

###########################
# MASHR OBJECT GENERATION #
###########################

run_mash_generation = function(arm_name, arm_cfg, ld_reference) {
  prefix = arm_cfg$prefix

  message("============================================================")
  message("Running mashr arm: ", arm_name, " [", arm_cfg$label, "]")
  message("============================================================")

  prepare_output_dirs(odir)

  # ------------------------------------------------------------
  # 1. FULL QTL STUDIES = actual mash columns
  # ------------------------------------------------------------
  studies = prepare_qtl_studies(arm_cfg)
  markers = colnames(studies[[1]]$gdata)
  study_names = names(studies)

  ld_submat = align_ld_matrix(ld_reference, markers)

  metadata_qtl = studyMetadata[studyMetadata$study %in% study_names, , drop = FALSE]
  metadata_qtl = metadata_qtl[match(study_names, metadata_qtl$study), , drop = FALSE]

  if (any(is.na(metadata_qtl$study))) {
    stop("Some mash study names were not found in studyMetadata.")
  }

  write_tsv_if_needed(
    metadata_qtl,
    prefixed_path(prefix, "metadata_qtl.tsv"),
    overwrite = overwrite
  )
  saveRDS(metadata_qtl, prefixed_path(prefix, "metadata_qtl.rds"))

  # ------------------------------------------------------------
  # 2. SUPPORT STUDIES = only for supported_analyte
  # ------------------------------------------------------------
  support_studies = load_or_compute_rds(
    prefixed_path(prefix, "studies_analyte.rds"),
    expr = prepare_support_studies(metadata_support_template, arm_cfg),
    overwrite = overwrite,
    label = paste0(prefix, "_support_studies [", arm_name, "]")
  )

  metadata_support = load_or_compute_rds(
    prefixed_path(prefix, "metadata_analyte_enriched.rds"),
    expr = {
      message("Collecting support-only gene and segregant metadata...")
      build_support_metadata(metadata_support_template, support_studies)
    },
    overwrite = overwrite,
    label = paste0(prefix, "_metadata_support_enriched [", arm_name, "]")
  )

  supported_analyte = load_or_compute_rds(
    prefixed_path(prefix, "supported_analyte.rds"),
    expr = {
      message("Computing supported analytes from support-only metadata...")
      suppTwo(metadata_support)
    },
    overwrite = overwrite,
    label = paste0(prefix, "_supported_analyte [", arm_name, "]")
  )

  write_lines_if_needed(
    supported_analyte,
    prefixed_path(prefix, "supported_analyte.txt"),
    overwrite = overwrite
  )

  message("Analyte-supported genes [", arm_name, "]: ", length(supported_analyte))

  # ------------------------------------------------------------
  # 3. RUN GENE SET
  # ------------------------------------------------------------
  # Define from support first.
  # Then only sanity-check against the full QTL studies used in mash.
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

  message(
    "Genes to run after intersecting support genes with full QTL studies [",
    arm_name, "]: ", length(run_genes)
  )

  saveRDS(run_genes, prefixed_path(prefix, "run_genes.rds"))
  write_lines_if_needed(
    run_genes,
    prefixed_path(prefix, "run_genes.txt"),
    overwrite = overwrite
  )

  # ------------------------------------------------------------
  # 4. FIT EFFECTS ACROSS THE FULL QTL STUDY SET
  # ------------------------------------------------------------
  gene_fits = load_or_compute_rda(
    prefixed_path(prefix, "gene_fits_supported_analyte.rda"),
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
    overwrite = overwrite,
    message_prefix = paste0("[", prefix, "] ")
  )

  gene_present = load_or_compute_rds(
    prefixed_path(prefix, "gene_present_supported_analyte.rds"),
    expr = {
      message("Summarizing per-gene dataset presence across full QTL studies...")
      gp = integer(length(run_genes))
      names(gp) = run_genes

      pb = txtProgressBar(min = 0, max = length(run_genes), style = 3)
      for (i in seq_along(run_genes)) {
        g = run_genes[i]
        fits = gene_fits[[g]]
        gp[i] = sum(vapply(
          fits,
          function(x) any(is.finite(x$beta) & is.finite(x$se)),
          logical(1)
        ))
        setTxtProgressBar(pb, i)
      }
      close(pb)

      gp
    },
    overwrite = overwrite,
    label = paste0(prefix, "_gene_present [", arm_name, "]")
  )

  # Build training rows gene by gene, but now prune markers within each gene
  # using the shared LD matrix before they are allowed into the pooled training set.
  training_rows = load_or_compute_rds(
    prefixed_path(prefix, "training_rows_supported_analyte_ldpruned.rds"),
    expr = {
      message("Building pooled training rows for covariance learning with LD pruning...")
      tr = vector("list", length(run_genes))
      ps = vector("list", length(run_genes))

      pb = txtProgressBar(min = 0, max = length(run_genes), style = 3)
      for (i in seq_along(run_genes)) {
        gene = run_genes[i]
        mats = build_gene_mats(gene, gene_fits, study_names, markers)

        keep_present = rowSums(is.finite(mats$Bhat) & is.finite(mats$Shat)) >= min_present

        if (any(keep_present)) {
          B0 = mats$Bhat[keep_present, , drop = FALSE]
          S0 = mats$Shat[keep_present, , drop = FALSE]
          gene_markers = rownames(B0)

          gene_ld = ld_submat[gene_markers, gene_markers, drop = FALSE]

          keep_ld = ld_prune_gene_rows(
            Bhat = B0,
            Shat = S0,
            ld_submat = gene_ld,
            abs_cor_threshold = ld_abs_cor_threshold,
            max_keep = max_ld_representatives_per_gene
          )

          tr[[i]] = list(
            Bhat = B0[keep_ld, , drop = FALSE],
            Shat = S0[keep_ld, , drop = FALSE]
          )

          ps[[i]] = data.frame(
            gene = gene,
            n_rows_before = nrow(B0),
            n_rows_after = length(keep_ld),
            stringsAsFactors = FALSE
          )

          rm(B0, S0, gene_markers, gene_ld, keep_ld)
        }

        setTxtProgressBar(pb, i)
        if (i %% 50L == 0L) gc()
      }
      close(pb)

      tr = tr[!vapply(tr, is.null, logical(1))]
      ps = ps[!vapply(ps, is.null, logical(1))]

      if (length(tr) == 0L) {
        stop("No training rows were collected after LD pruning.")
      }

      list(
        blocks = tr,
        prune_stats = do.call(rbind, ps)
      )
    },
    overwrite = overwrite,
    label = paste0(prefix, "_training_rows_ldpruned [", arm_name, "]")
  )

  prune_stats = training_rows$prune_stats
  write_tsv_if_needed(
    prune_stats,
    prefixed_path(prefix, "training_row_prune_stats.tsv"),
    overwrite = overwrite
  )

  message("Collected LD-pruned training blocks from ", length(training_rows$blocks), " genes [", arm_name, "].")
  message(
    "Mean rows per gene before pruning: ",
    round(mean(prune_stats$n_rows_before), 2),
    " | after pruning: ",
    round(mean(prune_stats$n_rows_after), 2)
  )

  training_mats_file = prefixed_path(prefix, "mash_training_mats_supported_analyte_ldpruned.rda")
  if (file.exists(training_mats_file) && !isTRUE(overwrite)) {
    message("Loading existing training matrices from: ", training_mats_file)

    env = new.env(parent = emptyenv())
    loaded_names = load(training_mats_file, envir = env)

    if ("training_mats" %in% loaded_names) {
      training_mats = env$training_mats
    } else if (all(c("train_Bhat", "train_Shat", "ss") %in% loaded_names)) {
      training_mats = list(
        train_Bhat = env$train_Bhat,
        train_Shat = env$train_Shat,
        ss = env$ss
      )
    } else {
      stop(
        "Existing training mats file does not contain either 'training_mats' or all of ",
        "'train_Bhat', 'train_Shat', and 'ss': ",
        training_mats_file
      )
    }

    rm(env, loaded_names)
    gc()
  } else {
    message("Row-binding pooled LD-pruned training matrices...")
    train_Bhat = do.call(rbind, lapply(training_rows$blocks, `[[`, "Bhat"))
    train_Shat = do.call(rbind, lapply(training_rows$blocks, `[[`, "Shat"))

    rownames(train_Bhat) = NULL
    rownames(train_Shat) = NULL

    message(
      "Training matrix dimensions before subsampling [", arm_name, "]: ",
      nrow(train_Bhat), " rows x ", ncol(train_Bhat), " studies"
    )

    if (nrow(train_Bhat) > max_training_rows) {
      message("Subsampling training rows to ", max_training_rows, " with stratification on max |Z|...")
      set.seed(arm_cfg$seed)

      ss = subsample_training_rows(
        train_Bhat,
        train_Shat,
        max_rows = max_training_rows,
        null_fraction = null_fraction,
        signal_fraction = signal_fraction,
        signal_quantile = signal_quantile
      )

      train_Bhat = ss$Bhat
      train_Shat = ss$Shat

      message("Signal threshold (", signal_quantile * 100, "th pct of max |Z|): ", round(ss$threshold, 3))
      message(
        "Sampled rows retained: ", nrow(train_Bhat),
        " [target null-ish fraction=", null_fraction,
        ", signal-ish fraction=", signal_fraction, "]"
      )
    } else {
      ss = list(
        idx = seq_len(nrow(train_Bhat)),
        row_score = apply(abs(train_Bhat / train_Shat), 1, max, na.rm = TRUE),
        threshold = NA_real_
      )
    }

    message(
      "Training matrix dimensions used for mash covariance learning [", arm_name, "]: ",
      nrow(train_Bhat), " rows x ", ncol(train_Bhat), " studies"
    )

    training_mats = list(
      train_Bhat = train_Bhat,
      train_Shat = train_Shat,
      ss = ss
    )

    save(training_mats, file = training_mats_file)
  }

  rm(training_rows)
  gc()

  train_Bhat = training_mats$train_Bhat
  train_Shat = training_mats$train_Shat

  Vhat = load_or_compute_rds(
    prefixed_path(prefix, "null_correlation.rds"),
    expr = {
      message("Using identity null correlation Vhat.")
      diag(length(study_names))
    },
    overwrite = overwrite,
    label = paste0(prefix, "_Vhat [", arm_name, "]")
  )

  write_tsv_if_needed(
    as.data.frame(Vhat),
    prefixed_path(prefix, "null_correlation.tsv"),
    overwrite = overwrite
  )

  train_data0 = load_or_compute_rds(
    prefixed_path(prefix, "train_data0.rds"),
    expr = {
      message("Constructing mash training data...")
      mash_set_data(train_Bhat, train_Shat, alpha = alpha)
    },
    overwrite = overwrite,
    label = paste0(prefix, "_train_data0 [", arm_name, "]")
  )

  train_data = load_or_compute_rds(
    prefixed_path(prefix, "train_data.rds"),
    expr = {
      mash_update_data(train_data0, V = Vhat)
    },
    overwrite = overwrite,
    label = paste0(prefix, "_train_data [", arm_name, "]")
  )

  npc = load_or_compute_rds(
    prefixed_path(prefix, "npc.rds"),
    expr = {
      min(5L, length(study_names))
    },
    overwrite = overwrite,
    label = paste0(prefix, "_npc [", arm_name, "]")
  )

  message("Using npc = ", npc, " [", arm_name, "]")

  U_c = load_or_compute_rds(
    prefixed_path(prefix, "cov_canonical.rds"),
    expr = {
      message("Computing canonical covariance candidates...")
      t0 = Sys.time()
      out = cov_canonical(train_data)
      message("Done canonical in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
      out
    },
    overwrite = overwrite,
    label = paste0(prefix, "_U_c [", arm_name, "]")
  )

  U_pca = load_or_compute_rds(
    prefixed_path(prefix, "cov_pca.rds"),
    expr = {
      message("Computing PCA covariance candidates...")
      t0 = Sys.time()
      out = cov_pca(train_data, npc = npc)
      message("Done PCA in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
      out
    },
    overwrite = overwrite,
    label = paste0(prefix, "_U_pca [", arm_name, "]")
  )

  U_ed = load_or_compute_rds(
    prefixed_path(prefix, "cov_ed.rds"),
    expr = {
      message("Computing ED covariance candidates...")
      t0 = Sys.time()
      out = cov_ed(train_data, U_pca)
      message("Done ED in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
      out
    },
    overwrite = overwrite,
    label = paste0(prefix, "_U_ed [", arm_name, "]")
  )

  Ulist = load_or_compute_rds(
    prefixed_path(prefix, "covariance_candidates.rds"),
    expr = {
      out = c(U_c, U_pca, U_ed)
      message("Total covariance candidates: ", length(out))
      out
    },
    overwrite = overwrite,
    label = paste0(prefix, "_Ulist [", arm_name, "]")
  )

  message("Total covariance candidates: ", length(Ulist), " [", arm_name, "]")

  global_mash = load_or_compute_rds(
    prefixed_path(prefix, "global_mash_model.rds"),
    expr = {
      message("Fitting global mash model on pooled rows...")
      t0 = Sys.time()
      out = mash(train_data, Ulist, outputlevel = 1)
      message("Done global mash in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
      out
    },
    overwrite = overwrite,
    label = paste0(prefix, "_global_mash [", arm_name, "]")
  )

  mix_df = load_or_compute_rds(
    prefixed_path(prefix, "mixture_weights.rds"),
    expr = {
      data.frame(
        covariance = names(get_estimated_pi(global_mash)),
        pi = as.numeric(get_estimated_pi(global_mash)),
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    },
    overwrite = overwrite,
    label = paste0(prefix, "_mixture_weights [", arm_name, "]")
  )

  write_tsv_if_needed(
    mix_df,
    prefixed_path(prefix, "mixture_weights.tsv"),
    overwrite = overwrite
  )

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
    mean_rows_before_ld_prune = mean(prune_stats$n_rows_before),
    mean_rows_after_ld_prune = mean(prune_stats$n_rows_after),
    stringsAsFactors = FALSE
  )

  write_tsv_if_needed(
    summary_info,
    prefixed_path(prefix, "analysis_summary.tsv"),
    overwrite = TRUE
  )

  out = list(
    arm_name = arm_name,
    prefix = prefix,
    analysis_dir = odir,
    study_names = study_names,
    n_run_genes = length(run_genes),
    ld_reference_study = ld_reference$study,
    mixture_weights = mix_df,
    summary_info = summary_info,
    files = list(
      metadata_qtl = prefixed_path(prefix, "metadata_qtl.rds"),
      metadata_support = prefixed_path(prefix, "metadata_analyte_enriched.rds"),
      supported_analyte = prefixed_path(prefix, "supported_analyte.rds"),
      run_genes = prefixed_path(prefix, "run_genes.rds"),
      gene_fits = prefixed_path(prefix, "gene_fits_supported_analyte.rda"),
      gene_present = prefixed_path(prefix, "gene_present_supported_analyte.rds"),
      training_rows = prefixed_path(prefix, "training_rows_supported_analyte_ldpruned.rds"),
      training_mats = prefixed_path(prefix, "mash_training_mats_supported_analyte_ldpruned.rda"),
      train_data = prefixed_path(prefix, "train_data.rds"),
      covariance_candidates = prefixed_path(prefix, "covariance_candidates.rds"),
      global_mash = prefixed_path(prefix, "global_mash_model.rds"),
      mixture_weights = prefixed_path(prefix, "mixture_weights.rds")
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
    ld_submat,
    training_mats,
    train_Bhat,
    train_Shat,
    Vhat,
    train_data0,
    train_data,
    U_c,
    U_pca,
    U_ed,
    Ulist,
    global_mash,
    prune_stats,
    summary_info,
    mix_df
  )
  gc()

  out
}

################
# RUN ALL ARMS #
################

prepare_output_dirs(odir)

# Build the shared LD matrix once.
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

message("Done.")