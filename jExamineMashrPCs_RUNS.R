#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# PURPOSE
# ------------------------------------------------------------------------------
# This script interprets the covariance candidates learned by mash.
#
# It now does two related things:
#
#   1. PCA-focused interpretation:
#      - Which PCA covariances line up with analyte / medium / study contrasts?
#      - Which PCA covariances actually carry nonzero pi?
#
#   2. Full covariance interpretation:
#      - Where does the mash mixture weight actually sit across ALL covariance
#        candidates (canonical, PCA, ED)?
#      - For the top-weighted covariance candidates, what metadata-derived
#        contrasts do their leading loading directions resemble most?
#
# Important caution:
#   - These are descriptive structure analyses, not formal tests of interaction.
#   - Because the design is not fully crossed across analyte and medium, any
#     analyte-versus-medium interpretation must be framed carefully.
#   - A covariance can align with a contrast without carrying much mixture mass.
#     Alignment tells us geometric similarity; pi tells us how often mash uses
#     that pattern.
# ------------------------------------------------------------------------------

source("/PHShome/kv1026/p05/remap_functions.R")
library(data.table)
library(ggplot2)

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------

overwrite = FALSE
n_perm = 10000L
perm_seed = 1L
active_pi_threshold = 1e-3
top_n_covariances_for_plots = 12L

prefix = "na"   # choose one: "na", "ds", or "fs"

valid_prefixes = c("na", "ds", "fs")
if (!(prefix %in% valid_prefixes)) {
  stop(
    "Invalid prefix: ", prefix,
    ". Choose one of: ",
    paste(valid_prefixes, collapse = ", ")
  )
}

message("Running covariance interpretation for arm: ", prefix)

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
# idir = stable input location for mash objects produced earlier
# odir = this script's output location
# ------------------------------------------------------------------------------

idir = paste0(imed, "mashr/", prefix, "_training/")
odir = paste0(imed, "mashr/", prefix, "_pc_contrast_analysis/")

dir.create(odir, recursive = TRUE, showWarnings = FALSE)

infile = function(filename) {
  file.path(idir, paste0(prefix, "_", filename))
}

outfile = function(filename) {
  file.path(odir, filename)
}

analysis_dir = odir

write_tsv_if_needed = function(dt, path, overwrite = FALSE, sep = "\t") {
  if (!file.exists(path) || isTRUE(overwrite)) {
    data.table::fwrite(dt, path, sep = sep)
  }
  invisible(path)
}

sanitize_for_filename = function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

# ------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------

load_rda_object = function(path, object_name = NULL) {
  env = new.env(parent = emptyenv())
  nm = load(path, envir = env)

  if (!is.null(object_name)) {
    return(env[[object_name]])
  }

  env[[nm]]
}

resolve_metadata_col = function(dt, candidates, required = TRUE) {
  nms = colnames(dt)
  hit = candidates[candidates %in% nms]

  if (length(hit) >= 1L) {
    return(hit[1])
  }

  if (isTRUE(required)) {
    stop(
      "Could not find required metadata column. Tried: ",
      paste(candidates, collapse = ", ")
    )
  }

  NULL
}

center_cosine = function(x, y) {
  ok = is.finite(x) & is.finite(y)
  x = x[ok]
  y = y[ok]

  if (length(x) < 3L) {
    return(NA_real_)
  }

  x = x - mean(x)
  y = y - mean(y)

  sx = sqrt(sum(x^2))
  sy = sqrt(sum(y^2))

  if (sx == 0 || sy == 0) {
    return(NA_real_)
  }

  sum(x * y) / (sx * sy)
}

permute_alignment_p = function(x, y, n_perm = 10000L, seed = 1L) {
  ok = is.finite(x) & is.finite(y)
  x = x[ok]
  y = y[ok]

  if (length(x) < 3L) {
    return(NA_real_)
  }

  x = x - mean(x)
  y = y - mean(y)

  sx = sqrt(sum(x^2))
  sy = sqrt(sum(y^2))

  if (sx == 0 || sy == 0) {
    return(NA_real_)
  }

  obs = sum(x * y) / (sx * sy)
  obs_abs = abs(obs)

  set.seed(seed)
  ge = 0L

  for (b in seq_len(n_perm)) {
    yp = sample(y, replace = FALSE)
    stat = sum(x * yp) / (sx * sqrt(sum(yp^2)))
    if (abs(stat) >= obs_abs) {
      ge = ge + 1L
    }
  }

  (ge + 1) / (n_perm + 1)
}

extract_covariance_loading = function(S) {
  ee = eigen(S, symmetric = TRUE)
  data.frame(
    loading = as.numeric(ee$vectors[, 1]),
    leading_eigenvalue = as.numeric(ee$values[1]),
    stringsAsFactors = FALSE
  )
}

make_subset_contrast = function(
  md,
  mask,
  group,
  positive_level,
  negative_level,
  name,
  family,
  details = NA_character_
) {
  out = rep(0, nrow(md))

  use = which(mask & !is.na(group) & group %in% c(positive_level, negative_level))
  if (length(use) < 3L) {
    return(NULL)
  }

  g = group[use]
  if (length(unique(g)) < 2L) {
    return(NULL)
  }

  out[use[g == positive_level]] = 1
  out[use[g == negative_level]] = -1

  data.frame(
    contrast = name,
    family = family,
    details = details,
    vector = I(list(out)),
    support_n = sum(out != 0),
    pos_n = sum(out == 1),
    neg_n = sum(out == -1),
    stringsAsFactors = FALSE
  )
}

make_one_vs_rest_contrast = function(
  md,
  mask,
  group,
  focal_level,
  name,
  family,
  details = NA_character_
) {
  out = rep(0, nrow(md))

  use = which(mask & !is.na(group))
  if (length(use) < 3L) {
    return(NULL)
  }

  g = group[use]
  if (!(focal_level %in% g)) {
    return(NULL)
  }

  if (length(unique(g)) < 2L) {
    return(NULL)
  }

  out[use[g == focal_level]] = 1
  out[use[g != focal_level]] = -1

  if (length(unique(out[use])) < 2L) {
    return(NULL)
  }

  data.frame(
    contrast = name,
    family = family,
    details = details,
    vector = I(list(out)),
    support_n = sum(out != 0),
    pos_n = sum(out == 1),
    neg_n = sum(out == -1),
    stringsAsFactors = FALSE
  )
}

weighted_alignment_score = function(abs_alignment, pi) {
  ifelse(is.finite(abs_alignment) & is.finite(pi), abs_alignment * pi, NA_real_)
}

escape_regex = function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

classify_covariance_source = function(cov_name, canonical_names, pca_names, ed_names) {
  if (cov_name %in% canonical_names) {
    return("canonical")
  }
  if (cov_name %in% pca_names) {
    return("pca")
  }
  if (cov_name %in% ed_names) {
    return("ed")
  }
  "unmatched"
}

collapse_mixture_to_covariances = function(mix_df, covariance_names, canonical_names, pca_names, ed_names) {
  out = data.table(
    covariance = covariance_names,
    source_family = vapply(
      covariance_names,
      classify_covariance_source,
      character(1),
      canonical_names = canonical_names,
      pca_names = pca_names,
      ed_names = ed_names
    ),
    pi = 0
  )

  mix_dt = as.data.table(mix_df)
  if (!all(c("covariance", "pi") %in% colnames(mix_dt))) {
    stop("mixture weights file must contain columns named 'covariance' and 'pi'")
  }

  unmatched_rows = vector("list", nrow(mix_dt))
  unmatched_idx = 1L

  cov_names_sorted = covariance_names[order(nchar(covariance_names), decreasing = TRUE)]

  for (i in seq_len(nrow(mix_dt))) {
    raw_name = as.character(mix_dt$covariance[i])
    raw_pi = as.numeric(mix_dt$pi[i])

    match_name = NULL

    if (raw_name %in% covariance_names) {
      match_name = raw_name
    } else {
      for (cov_name in cov_names_sorted) {
        patt = paste0("^", escape_regex(cov_name), "(_|\\.|:|$)")
        if (grepl(patt, raw_name)) {
          match_name = cov_name
          break
        }
      }
    }

    if (is.null(match_name)) {
      unmatched_rows[[unmatched_idx]] = data.table(
        raw_covariance = raw_name,
        pi = raw_pi
      )
      unmatched_idx = unmatched_idx + 1L
    } else {
      out[covariance == match_name, pi := pi + raw_pi]
    }
  }

  unmatched_rows = unmatched_rows[!vapply(unmatched_rows, is.null, logical(1))]
  unmatched_dt = if (length(unmatched_rows) > 0L) {
    rbindlist(unmatched_rows)
  } else {
    data.table(raw_covariance = character(), pi = numeric())
  }

  list(
    collapsed = out[order(-pi, covariance)],
    unmatched = unmatched_dt[order(-pi, raw_covariance)]
  )
}

build_covariance_long = function(Ulist, metadata_qtl, study_col) {
  cov_names = names(Ulist)
  if (is.null(cov_names)) {
    cov_names = paste0("cov_", seq_along(Ulist))
    names(Ulist) = cov_names
  }

  cov_long = vector("list", length(Ulist))
  cov_eigen = data.table(
    covariance = cov_names,
    leading_eigenvalue = NA_real_
  )

  for (i in seq_along(Ulist)) {
    tmp = extract_covariance_loading(Ulist[[i]])
    cov_long[[i]] = data.table(
      covariance = cov_names[i],
      study = metadata_qtl[[study_col]],
      loading = tmp$loading
    )
    cov_eigen$leading_eigenvalue[i] = tmp$leading_eigenvalue[1]
  }

  list(
    long = rbindlist(cov_long),
    eigen = cov_eigen
  )
}

run_alignment_against_contrasts = function(
  loading_wide,
  cov_names,
  contrast_def,
  n_perm = 10000L,
  perm_seed = 1L,
  cov_col = "covariance"
) {
  alignment_res = vector("list", length(cov_names) * nrow(contrast_def))
  idx = 1L

  for (cov_name in cov_names) {
    v_full = loading_wide[[cov_name]]

    for (j in seq_len(nrow(contrast_def))) {
      c_full = as.numeric(contrast_def$vector[[j]])
      mask = c_full != 0 & is.finite(v_full)

      if (sum(mask) < 3L) {
        next
      }

      v = v_full[mask]
      c = c_full[mask]

      stat = center_cosine(v, c)
      p_perm = permute_alignment_p(
        x = v,
        y = c,
        n_perm = n_perm,
        seed = perm_seed
      )

      alignment_res[[idx]] = data.table(
        covariance = cov_name,
        contrast = contrast_def$contrast[j],
        family = contrast_def$family[j],
        details = contrast_def$details[j],
        support_n = sum(mask),
        alignment = stat,
        abs_alignment = abs(stat),
        p_perm = p_perm
      )
      idx = idx + 1L
    }
  }

  alignment_res = alignment_res[!vapply(alignment_res, is.null, logical(1))]
  alignment_dt = rbindlist(alignment_res)
  alignment_dt[, fdr_global := p.adjust(p_perm, method = "BH")]
  alignment_dt[, fdr_within_family := p.adjust(p_perm, method = "BH"), by = family]
  setorder(alignment_dt, -abs_alignment, p_perm)

  alignment_dt
}

orient_loadings_to_best_contrast = function(loading_wide, cov_names, contrast_def, alignment_dt) {
  best_by_cov = alignment_dt[order(-abs_alignment, p_perm), .SD[1], by = covariance]

  orientation = data.table(
    covariance = cov_names,
    orient_to_contrast = NA_character_,
    sign_flip = 1L
  )

  oriented = copy(loading_wide)

  for (i in seq_len(nrow(best_by_cov))) {
    cov_name = best_by_cov$covariance[i]
    contrast_name = best_by_cov$contrast[i]

    v = oriented[[cov_name]]
    c = as.numeric(contrast_def$vector[[match(contrast_name, contrast_def$contrast)]])
    mask = c != 0 & is.finite(v)

    s = center_cosine(v[mask], c[mask])
    flip = ifelse(is.na(s) || s >= 0, 1L, -1L)

    orientation[covariance == cov_name, orient_to_contrast := contrast_name]
    orientation[covariance == cov_name, sign_flip := flip]
    oriented[[cov_name]] = flip * oriented[[cov_name]]
  }

  list(
    orientation = orientation,
    oriented_wide = oriented
  )
}

make_covariance_matrix_long = function(Ulist, metadata_qtl, study_col, covariance_names) {
  out = vector("list", length(covariance_names))

  for (i in seq_along(covariance_names)) {
    cov_name = covariance_names[i]
    S = Ulist[[cov_name]]
    dt = as.data.table(as.table(S))
    setnames(dt, c("i", "j", "value"))

    dt[, i := as.integer(i)]
    dt[, j := as.integer(j)]
    dt[, study_i := metadata_qtl[[study_col]][i]]
    dt[, study_j := metadata_qtl[[study_col]][j]]
    dt[, covariance := cov_name]

    out[[i]] = dt[, .(covariance, study_i, study_j, value)]
  }

  rbindlist(out)
}

# ------------------------------------------------------------------------------
# load arm-specific mash outputs
# ------------------------------------------------------------------------------

metadata_qtl = readRDS(infile("metadata_qtl.rds"))
U_c = readRDS(infile("cov_canonical.rds"))
U_pca = readRDS(infile("cov_pca.rds"))
U_ed = readRDS(infile("cov_ed.rds"))
Ulist = readRDS(infile("covariance_candidates.rds"))
mix_df = readRDS(infile("mixture_weights.rds"))

global_mash_path = infile("global_mash_model.rds")
if (!file.exists(global_mash_path)) {
  stop("Could not find global mash model: ", global_mash_path)
}
global_mash = readRDS(global_mash_path)

if (!is.data.frame(metadata_qtl)) {
  stop("metadata_qtl is not a data.frame")
}

if (!is.list(U_pca) || length(U_pca) == 0L) {
  stop("U_pca is empty or not a list")
}
if (!is.list(Ulist) || length(Ulist) == 0L) {
  stop("Ulist is empty or not a list")
}

study_col = resolve_metadata_col(metadata_qtl, c("study"))
analyte_col = resolve_metadata_col(metadata_qtl, c("analyte"), required = FALSE)
medium_col = resolve_metadata_col(metadata_qtl, c("medium"), required = FALSE)

K = nrow(Ulist[[1]])
if (nrow(metadata_qtl) != K) {
  stop(
    "metadata_qtl rows (", nrow(metadata_qtl),
    ") do not match covariance dimension (", K, ")"
  )
}

if (is.null(names(U_c))) {
  names(U_c) = paste0("canonical_", seq_along(U_c))
}
if (is.null(names(U_pca))) {
  names(U_pca) = c(paste0("PCA_", seq_len(length(U_pca) - 1L)), "tPCA")
}
if (is.null(names(U_ed))) {
  names(U_ed) = paste0("ED_", seq_along(U_ed))
}
if (is.null(names(Ulist))) {
  names(Ulist) = c(names(U_c), names(U_pca), names(U_ed))
}

canonical_names = names(U_c)
pca_names = names(U_pca)
ed_names = names(U_ed)
covariance_names = names(Ulist)

message("Loaded metadata_qtl with ", nrow(metadata_qtl), " studies")
message("Loaded ", length(canonical_names), " canonical covariance candidates")
message("Loaded ", length(pca_names), " PCA covariance candidates")
message("Loaded ", length(ed_names), " ED covariance candidates")
message("Loaded ", length(covariance_names), " total covariance candidates")

# ------------------------------------------------------------------------------
# analysis 0: extract one loading vector per PCA covariance
# ------------------------------------------------------------------------------
# explicit PCA on training Z matrix
# ------------------------------------------------------------------------------
# We compute PCA ourselves from the arm-specific mash training matrix:
#   Z = Bhat / Shat
#
# The pooled training matrix is sparse across studies, so we:
#   1. center each study column using observed values only
#   2. replace remaining missing entries with 0 after centering
#      (equivalent to mean imputation)
#   3. run prcomp on the centered, mean-imputed matrix
# ------------------------------------------------------------------------------

training_mats = load_rda_object(
  infile("mash_training_mats_supported_analyte_ldpruned.rda"),
  object_name = "training_mats"
)

train_Bhat = training_mats$train_Bhat
train_Shat = training_mats$train_Shat

Z = train_Bhat / train_Shat
Z[!is.finite(Z)] = NA_real_

if (ncol(Z) != nrow(metadata_qtl)) {
  stop(
    "Number of study columns in Z (", ncol(Z),
    ") does not match metadata_qtl rows (", nrow(metadata_qtl), ")."
  )
}

# useful diagnostics on sparsity / overlap
obs_mask = is.finite(Z)

study_nonmissing = data.table(
  study = metadata_qtl[[study_col]],
  n_nonmissing = colSums(obs_mask),
  fraction_nonmissing = colMeans(obs_mask)
)

pairwise_overlap = obs_mask * 1L
pairwise_overlap = t(pairwise_overlap) %*% pairwise_overlap
pairwise_overlap = as.data.table(pairwise_overlap, keep.rownames = "study")
setnames(pairwise_overlap, c("study", metadata_qtl[[study_col]]))

write_tsv_if_needed(
  study_nonmissing,
  outfile("global_pca_nonmissing_by_study.tsv"),
  overwrite = overwrite
)

write_tsv_if_needed(
  pairwise_overlap,
  outfile("global_pca_pairwise_overlap.tsv"),
  overwrite = overwrite
)

# center each column using observed values only
col_means = colMeans(Z, na.rm = TRUE)
for (j in seq_len(ncol(Z))) {
  ok = is.finite(Z[, j])
  Z[ok, j] = Z[ok, j] - col_means[j]
}

# mean imputation after centering
Z[!is.finite(Z)] = 0

if (nrow(Z) < 2L || ncol(Z) < 2L) {
  stop("Not enough data to run PCA.")
}

pca = prcomp(Z, center = FALSE, scale. = FALSE)

var_explained = (pca$sdev^2) / sum(pca$sdev^2)

pca_variance = data.table(
  pc = paste0("PC_", seq_along(var_explained)),
  variance_explained = var_explained,
  cumulative_variance = cumsum(var_explained)
)

pca_loadings = data.table(
  study = metadata_qtl[[study_col]]
)

for (k in seq_len(ncol(pca$rotation))) {
  pca_loadings[[paste0("PC_", k)]] = pca$rotation[, k]
}

write_tsv_if_needed(
  pca_variance,
  outfile("global_pca_variance_explained.tsv"),
  overwrite = overwrite
)

write_tsv_if_needed(
  pca_loadings,
  outfile("global_pca_loadings_by_study.tsv"),
  overwrite = overwrite
)

message("\nGlobal PCA variance explained:")
print(pca_variance[1:min(10, nrow(pca_variance))])

# ------------------------------------------------------------------------------
# analysis 1: build metadata-driven contrasts
# ------------------------------------------------------------------------------

contrast_rows = list()

if (!is.null(analyte_col)) {
  analyte_vals = as.character(metadata_qtl[[analyte_col]])
  analyte_levels = sort(unique(analyte_vals[!is.na(analyte_vals)]))

  if (length(analyte_levels) >= 2L) {
    cc = utils::combn(analyte_levels, 2, simplify = FALSE)
    for (pair in cc) {
      nm = paste0("analyte__", pair[1], "_vs_", pair[2])
      contrast_rows[[length(contrast_rows) + 1L]] = make_subset_contrast(
        md = metadata_qtl,
        mask = rep(TRUE, nrow(metadata_qtl)),
        group = analyte_vals,
        positive_level = pair[1],
        negative_level = pair[2],
        name = nm,
        family = "analyte_pairwise",
        details = "pairwise analyte contrast"
      )
    }
  }
}

if (!is.null(medium_col)) {
  medium_vals = as.character(metadata_qtl[[medium_col]])
  medium_levels = sort(unique(medium_vals[!is.na(medium_vals)]))

  if (length(medium_levels) >= 2L) {
    cc = utils::combn(medium_levels, 2, simplify = FALSE)
    for (pair in cc) {
      nm = paste0("medium__", pair[1], "_vs_", pair[2])
      contrast_rows[[length(contrast_rows) + 1L]] = make_subset_contrast(
        md = metadata_qtl,
        mask = rep(TRUE, nrow(metadata_qtl)),
        group = medium_vals,
        positive_level = pair[1],
        negative_level = pair[2],
        name = nm,
        family = "medium_pairwise",
        details = "pairwise medium contrast across all studies"
      )
    }
  }
}

if (!is.null(analyte_col) && !is.null(medium_col)) {
  analyte_vals = as.character(metadata_qtl[[analyte_col]])
  medium_vals = as.character(metadata_qtl[[medium_col]])
  analyte_levels = sort(unique(analyte_vals[!is.na(analyte_vals)]))
  medium_levels = sort(unique(medium_vals[!is.na(medium_vals)]))

  for (a in analyte_levels) {
    mask_a = analyte_vals == a
    med_a = sort(unique(medium_vals[mask_a & !is.na(medium_vals)]))

    if (length(med_a) >= 2L) {
      cc = utils::combn(med_a, 2, simplify = FALSE)
      for (pair in cc) {
        nm = paste0("within_", a, "__", pair[1], "_vs_", pair[2])
        contrast_rows[[length(contrast_rows) + 1L]] = make_subset_contrast(
          md = metadata_qtl,
          mask = mask_a,
          group = medium_vals,
          positive_level = pair[1],
          negative_level = pair[2],
          name = nm,
          family = "medium_within_analyte",
          details = paste0("medium contrast within analyte=", a)
        )
      }
    }
  }

  for (m in medium_levels) {
    mask_m = medium_vals == m
    analyte_m = sort(unique(analyte_vals[mask_m & !is.na(analyte_vals)]))

    if (length(analyte_m) >= 2L) {
      cc = utils::combn(analyte_m, 2, simplify = FALSE)
      for (pair in cc) {
        nm = paste0("within_", m, "__", pair[1], "_vs_", pair[2])
        contrast_rows[[length(contrast_rows) + 1L]] = make_subset_contrast(
          md = metadata_qtl,
          mask = mask_m,
          group = analyte_vals,
          positive_level = pair[1],
          negative_level = pair[2],
          name = nm,
          family = "analyte_within_medium",
          details = paste0("analyte contrast within medium=", m)
        )
      }
    }
  }
}

study_vals = as.character(metadata_qtl[[study_col]])
study_levels = sort(unique(study_vals[!is.na(study_vals)]))

for (s in study_levels) {
  nm = paste0("study__", s, "_vs_rest")
  contrast_rows[[length(contrast_rows) + 1L]] = make_one_vs_rest_contrast(
    md = metadata_qtl,
    mask = rep(TRUE, nrow(metadata_qtl)),
    group = study_vals,
    focal_level = s,
    name = nm,
    family = "study_one_vs_rest",
    details = "study-specific one-vs-rest contrast"
  )
}

contrast_rows = contrast_rows[!vapply(contrast_rows, is.null, logical(1))]
if (length(contrast_rows) == 0L) {
  stop("No usable contrasts could be generated from metadata_qtl")
}

contrast_def = rbindlist(contrast_rows, fill = TRUE)
contrast_mat = do.call(cbind, lapply(contrast_def$vector, as.numeric))
colnames(contrast_mat) = contrast_def$contrast
rownames(contrast_mat) = metadata_qtl[[study_col]]

write_tsv_if_needed(
  contrast_def[, .(contrast, family, details, support_n, pos_n, neg_n)],
  outfile("contrast_definitions.tsv"),
  overwrite = overwrite
)

contrast_vectors_dt = data.table(
  study = metadata_qtl[[study_col]],
  contrast_mat,
  check.names = FALSE
)

write_tsv_if_needed(
  contrast_vectors_dt,
  outfile("contrast_vectors_by_study.tsv"),
  overwrite = overwrite
)

# ------------------------------------------------------------------------------
# analysis 2: full mixture-weight analysis across all covariances
# ------------------------------------------------------------------------------

pi_map = collapse_mixture_to_covariances(
  mix_df = mix_df,
  covariance_names = covariance_names,
  canonical_names = canonical_names,
  pca_names = pca_names,
  ed_names = ed_names
)

covariance_pi_dt = pi_map$collapsed
unmatched_mixture_dt = pi_map$unmatched

covariance_pi_dt[, is_active := is.finite(pi) & pi > active_pi_threshold]

write_tsv_if_needed(
  covariance_pi_dt,
  outfile("covariance_pi_collapsed.tsv"),
  overwrite = overwrite
)

write_tsv_if_needed(
  unmatched_mixture_dt,
  outfile("covariance_pi_unmatched_raw.tsv"),
  overwrite = overwrite
)

pi_by_source_family = covariance_pi_dt[
  , .(
    n_covariance = .N,
    n_active_covariance = sum(is_active, na.rm = TRUE),
    total_pi = sum(pi, na.rm = TRUE),
    active_pi = sum(pi[is_active], na.rm = TRUE)
  ),
  by = source_family
][order(-active_pi, -total_pi, source_family)]

write_tsv_if_needed(
  pi_by_source_family,
  outfile("pi_by_source_family.tsv"),
  overwrite = overwrite
)

top_covariances_by_pi = covariance_pi_dt[order(-pi, covariance)]
write_tsv_if_needed(
  top_covariances_by_pi,
  outfile("top_covariances_by_pi.tsv"),
  overwrite = overwrite
)

# ------------------------------------------------------------------------------
# analysis 3: extract leading loading vector for every covariance candidate
# ------------------------------------------------------------------------------

cov_loading_obj = build_covariance_long(Ulist, metadata_qtl, study_col)
cov_loadings_long = cov_loading_obj$long
cov_eigen_dt = cov_loading_obj$eigen

cov_loadings_wide = dcast(cov_loadings_long, study ~ covariance, value.var = "loading")
cov_loadings_wide = merge(
  data.table(metadata_qtl),
  cov_loadings_wide,
  by.x = study_col,
  by.y = "study",
  all.x = TRUE,
  sort = FALSE
)

write_tsv_if_needed(
  cov_eigen_dt,
  outfile("covariance_leading_eigenvalues.tsv"),
  overwrite = overwrite
)

write_tsv_if_needed(
  cov_loadings_wide,
  outfile("covariance_loadings_by_study.tsv"),
  overwrite = overwrite
)

# ------------------------------------------------------------------------------
# analysis 4: align ALL covariance candidates to metadata-derived contrasts
# ------------------------------------------------------------------------------

cov_alignment_dt = run_alignment_against_contrasts(
  loading_wide = cov_loadings_wide,
  cov_names = covariance_names,
  contrast_def = contrast_def,
  n_perm = n_perm,
  perm_seed = perm_seed
)

cov_alignment_dt = merge(
  cov_alignment_dt,
  covariance_pi_dt,
  by = "covariance",
  all.x = TRUE,
  sort = FALSE
)

cov_alignment_dt[, weighted_alignment := weighted_alignment_score(abs_alignment, pi)]

write_tsv_if_needed(
  cov_alignment_dt,
  outfile("covariance_contrast_alignment.tsv"),
  overwrite = overwrite
)

top_contrast_by_covariance = cov_alignment_dt[order(-abs_alignment, p_perm), .SD[1], by = covariance]
top_contrast_by_covariance[, is_active := is.finite(pi) & pi > active_pi_threshold]
top_contrast_by_covariance[, weighted_alignment := weighted_alignment_score(abs_alignment, pi)]
setorder(top_contrast_by_covariance, -pi, -weighted_alignment, -abs_alignment, covariance)

write_tsv_if_needed(
  top_contrast_by_covariance,
  outfile("top_contrast_per_covariance.tsv"),
  overwrite = overwrite
)

top_covariances_interpreted = top_contrast_by_covariance[
  order(-pi, -weighted_alignment, -abs_alignment),
  .(
    covariance,
    source_family,
    pi,
    is_active,
    best_contrast = contrast,
    best_family = family,
    details,
    support_n,
    alignment,
    abs_alignment,
    weighted_alignment,
    p_perm,
    fdr_global,
    fdr_within_family
  )
]

write_tsv_if_needed(
  top_covariances_interpreted,
  outfile("top_covariances_by_pi_interpreted.tsv"),
  overwrite = overwrite
)

pi_by_top_contrast_family = top_contrast_by_covariance[
  , .(
    n_covariance = .N,
    n_active_covariance = sum(is_active, na.rm = TRUE),
    total_pi = sum(pi, na.rm = TRUE),
    active_pi = sum(pi[is_active], na.rm = TRUE),
    total_weighted_alignment = sum(weighted_alignment, na.rm = TRUE)
  ),
  by = .(best_family = family)
][order(-active_pi, -total_pi, -total_weighted_alignment)]

write_tsv_if_needed(
  pi_by_top_contrast_family,
  outfile("pi_by_top_contrast_family.tsv"),
  overwrite = overwrite
)

# ------------------------------------------------------------------------------
# analysis 5: collapse pi by best-matching contrast family
# ------------------------------------------------------------------------------

alignment_threshold = 0.5

top_contrast_by_covariance[
  , best_family_with_other := ifelse(
    abs_alignment >= alignment_threshold,
    family,
    "other"
  )
]

pi_by_best_family = top_contrast_by_covariance[
  , .(
    n_covariance = .N,
    n_active_covariance = sum(is_active, na.rm = TRUE),
    total_pi = sum(pi, na.rm = TRUE),
    active_pi = sum(pi[is_active], na.rm = TRUE)
  ),
  by = .(best_family = best_family_with_other)
][order(-total_pi)]

write_tsv_if_needed(
  pi_by_best_family,
  outfile("pi_by_best_family.tsv"),
  overwrite = overwrite
)

message("\nPi by best-matching contrast family:")
print(pi_by_best_family)

medium_like_pi = pi_by_best_family[
  best_family %in% c("medium_pairwise", "medium_within_analyte"),
  sum(total_pi, na.rm = TRUE)
]

analyte_like_pi = pi_by_best_family[
  best_family %in% c("analyte_pairwise", "analyte_within_medium"),
  sum(total_pi, na.rm = TRUE)
]

other_like_pi = pi_by_best_family[
  best_family == "other",
  sum(total_pi, na.rm = TRUE)
]

message("\nTotal medium-like pi: ", sprintf("%.4f", medium_like_pi))
message("Total analyte-like pi: ", sprintf("%.4f", analyte_like_pi))
message("Total other pi: ", sprintf("%.4f", other_like_pi))

# ------------------------------------------------------------------------------
# analysis 6: orient ALL covariance loadings for display
# ------------------------------------------------------------------------------

orient_all = orient_loadings_to_best_contrast(
  loading_wide = cov_loadings_wide,
  cov_names = covariance_names,
  contrast_def = contrast_def,
  alignment_dt = cov_alignment_dt
)

cov_orientation = orient_all$orientation
cov_loadings_wide_oriented = orient_all$oriented_wide

write_tsv_if_needed(
  cov_orientation,
  outfile("covariance_orientation.tsv"),
  overwrite = overwrite
)

write_tsv_if_needed(
  cov_loadings_wide_oriented,
  outfile("covariance_loadings_by_study_oriented.tsv"),
  overwrite = overwrite
)

cov_alignment_oriented_dt = run_alignment_against_contrasts(
  loading_wide = cov_loadings_wide_oriented,
  cov_names = covariance_names,
  contrast_def = contrast_def,
  n_perm = 0L,
  perm_seed = perm_seed
)

cov_alignment_oriented_dt = cov_alignment_oriented_dt[
  , .(covariance, contrast, family, details, support_n, alignment, abs_alignment)
]
cov_alignment_oriented_dt = merge(
  cov_alignment_oriented_dt,
  covariance_pi_dt,
  by = "covariance",
  all.x = TRUE,
  sort = FALSE
)
cov_alignment_oriented_dt[, weighted_alignment := weighted_alignment_score(abs(alignment), pi)]

write_tsv_if_needed(
  cov_alignment_oriented_dt,
  outfile("covariance_contrast_alignment_oriented.tsv"),
  overwrite = overwrite
)

# ------------------------------------------------------------------------------
# analysis 7: PCA-only convenience tables for continuity with earlier workflow
# ------------------------------------------------------------------------------

pc_names = pca_names

pc_pi_dt = covariance_pi_dt[covariance %in% pc_names]
setnames(pc_pi_dt, "covariance", "pc")

pc_alignment_dt = cov_alignment_dt[covariance %in% pc_names]
setnames(pc_alignment_dt, "covariance", "pc")

pc_alignment_oriented_dt = cov_alignment_oriented_dt[covariance %in% pc_names]
setnames(pc_alignment_oriented_dt, "covariance", "pc")

pc_orientation = cov_orientation[covariance %in% pc_names]
setnames(pc_orientation, "covariance", "pc")

pc_loadings_wide = cov_loadings_wide[, c(colnames(metadata_qtl), pc_names), with = FALSE]
pc_loadings_wide_oriented = cov_loadings_wide_oriented[, c(colnames(metadata_qtl), pc_names), with = FALSE]

write_tsv_if_needed(
  pc_pi_dt,
  outfile("pc_pi_collapsed.tsv"),
  overwrite = overwrite
)

write_tsv_if_needed(
  pc_alignment_dt,
  outfile("pc_contrast_alignment_with_weights.tsv"),
  overwrite = overwrite
)

write_tsv_if_needed(
  pc_alignment_oriented_dt,
  outfile("pc_contrast_alignment_oriented_with_weights.tsv"),
  overwrite = overwrite
)

write_tsv_if_needed(
  pc_orientation,
  outfile("pc_orientation.tsv"),
  overwrite = overwrite
)

write_tsv_if_needed(
  pc_loadings_wide,
  outfile("pc_loadings_by_study.tsv"),
  overwrite = overwrite
)

write_tsv_if_needed(
  pc_loadings_wide_oriented,
  outfile("pc_loadings_by_study_oriented.tsv"),
  overwrite = overwrite
)

top_by_pc = pc_alignment_dt[order(-abs_alignment, p_perm), .SD[1], by = pc]
top_by_pc[, is_active := is.finite(pi) & pi > active_pi_threshold]
top_by_pc[, weighted_alignment := weighted_alignment_score(abs_alignment, pi)]
setorder(top_by_pc, -pi, -weighted_alignment, -abs_alignment, pc)

write_tsv_if_needed(
  top_by_pc,
  outfile("top_contrast_per_pc.tsv"),
  overwrite = overwrite
)

top_by_contrast = pc_alignment_dt[order(-abs_alignment, p_perm), .SD[1], by = contrast]
write_tsv_if_needed(
  top_by_contrast,
  outfile("top_pc_per_contrast.tsv"),
  overwrite = overwrite
)

pi_by_top_contrast = top_by_pc[
  , .(
    n_pc = .N,
    n_active_pc = sum(is_active, na.rm = TRUE),
    total_pi = sum(pi, na.rm = TRUE),
    active_pi = sum(pi[is_active], na.rm = TRUE),
    total_weighted_alignment = sum(weighted_alignment, na.rm = TRUE)
  ),
  by = .(contrast, family, details)
][order(-active_pi, -total_pi, -total_weighted_alignment)]

write_tsv_if_needed(
  pi_by_top_contrast,
  outfile("pi_by_top_contrast.tsv"),
  overwrite = overwrite
)

pi_by_family = top_by_pc[
  , .(
    n_pc = .N,
    n_active_pc = sum(is_active, na.rm = TRUE),
    total_pi = sum(pi, na.rm = TRUE),
    active_pi = sum(pi[is_active], na.rm = TRUE),
    median_abs_alignment = median(abs_alignment, na.rm = TRUE),
    total_weighted_alignment = sum(weighted_alignment, na.rm = TRUE)
  ),
  by = .(family)
][order(-active_pi, -total_pi, -total_weighted_alignment)]

write_tsv_if_needed(
  pi_by_family,
  outfile("pi_by_family.tsv"),
  overwrite = overwrite
)

# ------------------------------------------------------------------------------
# analysis 8: targeted summaries
# ------------------------------------------------------------------------------

target_contrasts = c(
  "analyte__mRNA_vs_protein",
  "within_SC__mRNA_vs_protein",
  "within_mRNA__SC_vs_YNB2G",
  "within_mRNA__SC_vs_YNB1EtOH"
)

target_summary = cov_alignment_dt[
  contrast %in% target_contrasts,
  .(
    best_abs_alignment = max(abs_alignment, na.rm = TRUE),
    best_weighted_alignment = max(weighted_alignment, na.rm = TRUE),
    total_pi_from_matching_top_covariance = sum(
      top_contrast_by_covariance[contrast == .BY$contrast]$pi,
      na.rm = TRUE
    )
  ),
  by = contrast
]

write_tsv_if_needed(
  target_summary,
  outfile("target_contrast_summary.tsv"),
  overwrite = overwrite
)

if ("analyte__mRNA_vs_protein" %in% cov_alignment_dt$contrast) {
  mrna_protein_rank = cov_alignment_dt[contrast == "analyte__mRNA_vs_protein"]
  setorder(mrna_protein_rank, -pi, -abs_alignment, -weighted_alignment, p_perm)

  write_tsv_if_needed(
    mrna_protein_rank,
    outfile("mRNA_vs_protein_covariance_rank.tsv"),
    overwrite = overwrite
  )
}

if ("within_SC__mRNA_vs_protein" %in% cov_alignment_dt$contrast) {
  sc_analyte_rank = cov_alignment_dt[contrast == "within_SC__mRNA_vs_protein"]
  setorder(sc_analyte_rank, -pi, -abs_alignment, -weighted_alignment, p_perm)

  write_tsv_if_needed(
    sc_analyte_rank,
    outfile("within_SC_mRNA_vs_protein_covariance_rank.tsv"),
    overwrite = overwrite
  )
}

# ------------------------------------------------------------------------------
# figures
# ------------------------------------------------------------------------------

# 1. PCA alignment heatmap
heat_dt = copy(pc_alignment_oriented_dt)
heat_dt[, pc := factor(pc, levels = pc_names)]
contrast_levels = unique(pc_alignment_dt[order(family, contrast)]$contrast)
heat_dt[, contrast := factor(contrast, levels = contrast_levels)]

p_heat = ggplot(heat_dt, aes(x = contrast, y = pc, fill = alignment)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", alignment)), size = 3) +
  scale_fill_gradient2(limits = c(-1, 1)) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = paste0(prefix, ": PCA loading alignment with metadata-derived contrasts"),
    subtitle = "Centered cosine alignment on relevant study subsets",
    x = "Contrast",
    y = "PCA covariance"
  )

ggsave(
  filename = outfile("pc_contrast_alignment_heatmap.pdf"),
  plot = p_heat,
  width = max(10, 0.45 * length(unique(heat_dt$contrast))),
  height = max(4, 0.50 * length(unique(heat_dt$pc)))
)

# 2. PCA pi bar plot
p_pc_pi = ggplot(pc_pi_dt[order(pi)], aes(x = reorder(pc, pi), y = pi)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = paste0(prefix, ": collapsed pi by PCA covariance"),
    subtitle = "Pi collapsed across all effect-size scales for each PCA covariance",
    x = "PCA covariance",
    y = "Collapsed pi"
  )

ggsave(
  filename = outfile("pc_pi_barplot.pdf"),
  plot = p_pc_pi,
  width = 7,
  height = 4.5
)

# 3. Full source-family pi plot
p_source = ggplot(pi_by_source_family, aes(x = reorder(source_family, active_pi), y = active_pi)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = paste0(prefix, ": active mixture weight by covariance source family"),
    subtitle = "Collapsed pi across canonical, PCA, and ED covariance candidates",
    x = "Source family",
    y = "Active pi"
  )

ggsave(
  filename = outfile("pi_by_source_family.pdf"),
  plot = p_source,
  width = 7,
  height = 4.5
)

# 4. Top weighted covariances
top_cov_plot_dt = top_covariances_interpreted[seq_len(min(top_n_covariances_for_plots, .N))]
top_cov_plot_dt[, covariance_label := paste0(covariance, "\n", best_contrast)]

p_top_cov = ggplot(
  top_cov_plot_dt[order(pi)],
  aes(x = reorder(covariance_label, pi), y = pi, fill = source_family)
) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = paste0(prefix, ": top covariance candidates by collapsed pi"),
    subtitle = "Labels show best-matching metadata contrast for each covariance",
    x = "Covariance candidate",
    y = "Collapsed pi",
    fill = "Source"
  )

ggsave(
  filename = outfile("top_covariances_by_pi.pdf"),
  plot = p_top_cov,
  width = 9,
  height = max(5, 0.45 * nrow(top_cov_plot_dt))
)

# 5. Heatmap of covariance matrices for top weighted covariances
top_cov_names = top_covariances_interpreted$covariance[seq_len(min(top_n_covariances_for_plots, nrow(top_covariances_interpreted)))]
top_cov_matrix_dt = make_covariance_matrix_long(Ulist, metadata_qtl, study_col, top_cov_names)
top_cov_matrix_dt[, covariance := factor(covariance, levels = top_cov_names)]

p_cov_mats = ggplot(top_cov_matrix_dt, aes(x = study_j, y = study_i, fill = value)) +
  geom_tile() +
  facet_wrap(~ covariance, scales = "free", ncol = 3) +
  scale_fill_gradient2() +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = paste0(prefix, ": covariance matrices for top pi-weighted candidates"),
    x = "Study",
    y = "Study",
    fill = "Covariance"
  )

ggsave(
  filename = outfile("top_covariance_matrices_heatmap.pdf"),
  plot = p_cov_mats,
  width = 13,
  height = max(6, 2.8 * ceiling(length(top_cov_names) / 3))
)

# 6. PCA loading plots
plot_dt = as.data.table(pc_loadings_wide_oriented)

for (pc_name in pc_names) {
  sub = top_by_pc[pc == pc_name]
  if (nrow(sub) == 0L) {
    next
  }

  subtitle_txt = paste0(
    "best contrast = ", sub$contrast[1],
    "; |alignment| = ", sprintf("%.3f", sub$abs_alignment[1]),
    "; pi = ", sprintf("%.4f", sub$pi[1]),
    "; weighted = ", sprintf("%.4f", sub$weighted_alignment[1])
  )

  this_dt = copy(plot_dt)
  this_dt[, plot_y := get(pc_name)]
  this_dt[, plot_x := factor(get(study_col), levels = get(study_col)[order(plot_y)])]

  p = ggplot(this_dt, aes(x = plot_x, y = plot_y)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(size = 3) +
    coord_flip() +
    theme_bw(base_size = 11) +
    labs(
      title = paste0(prefix, ": ", pc_name, " loadings across studies"),
      subtitle = subtitle_txt,
      x = "Study",
      y = paste0(pc_name, " loading")
    )

  if (!is.null(analyte_col) && !is.null(medium_col)) {
    p = p + geom_point(aes(color = get(analyte_col), shape = get(medium_col)), size = 3)
  } else if (!is.null(analyte_col)) {
    p = p + geom_point(aes(color = get(analyte_col)), size = 3)
  } else if (!is.null(medium_col)) {
    p = p + geom_point(aes(color = get(medium_col)), size = 3)
  }

  ggsave(
    filename = outfile(paste0(sanitize_for_filename(pc_name), "_loadings.pdf")),
    plot = p,
    width = 8,
    height = 5
  )
}

plot_top_contrast_membership = function(
  contrast_name,
  filename_stub,
  point_size = 3
) {
  if (!(contrast_name %in% contrast_def$contrast)) {
    return(invisible(NULL))
  }

  sub = cov_alignment_dt[contrast == contrast_name][order(-pi, -abs_alignment, -weighted_alignment, p_perm)]
  if (nrow(sub) == 0L) {
    return(invisible(NULL))
  }

  top_cov = sub$covariance[1]
  cvec = as.numeric(contrast_def$vector[[match(contrast_name, contrast_def$contrast)]])
  dt = copy(cov_loadings_wide_oriented)
  dt[, contrast_member := cvec]
  dt = dt[contrast_member != 0]

  dt[, plot_y := get(top_cov)]
  dt[, plot_x := factor(get(study_col), levels = get(study_col)[order(plot_y)])]

  p = ggplot(
    dt,
    aes(
      x = plot_x,
      y = plot_y,
      color = factor(contrast_member)
    )
  ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(size = point_size) +
    coord_flip() +
    theme_bw(base_size = 11) +
    labs(
      title = paste0(prefix, ": top weighted covariance for ", contrast_name),
      subtitle = paste0("top covariance = ", top_cov),
      x = "Study",
      y = paste0(top_cov, " loading"),
      color = "Contrast group"
    )

  ggsave(
    filename = outfile(paste0(filename_stub, ".pdf")),
    plot = p,
    width = 8,
    height = 5
  )

  invisible(NULL)
}

if ("analyte__mRNA_vs_protein" %in% contrast_def$contrast) {
  plot_top_contrast_membership(
    contrast_name = "analyte__mRNA_vs_protein",
    filename_stub = "top_covariance_mRNA_vs_protein"
  )
}

if ("within_SC__mRNA_vs_protein" %in% contrast_def$contrast) {
  plot_top_contrast_membership(
    contrast_name = "within_SC__mRNA_vs_protein",
    filename_stub = "top_covariance_within_SC_mRNA_vs_protein"
  )
}

# ------------------------------------------------------------------------------
# plain-English summary table
# ------------------------------------------------------------------------------

summary_lines = c(
  paste0("Arm analyzed: ", prefix),
  paste0("Number of total covariance candidates: ", length(covariance_names)),
  paste0("Number of active covariance candidates (pi > ", active_pi_threshold, "): ", sum(covariance_pi_dt$is_active, na.rm = TRUE)),
  paste0(
    "Top covariance by pi: ",
    if (nrow(top_covariances_interpreted) > 0) top_covariances_interpreted$covariance[1] else NA_character_
  ),
  paste0(
    "Top covariance pi: ",
    if (nrow(top_covariances_interpreted) > 0) sprintf("%.4f", top_covariances_interpreted$pi[1]) else NA_character_
  ),
  paste0(
    "Top covariance best contrast: ",
    if (nrow(top_covariances_interpreted) > 0) top_covariances_interpreted$best_contrast[1] else NA_character_
  ),
  paste0(
    "Top source family by active pi: ",
    if (nrow(pi_by_source_family) > 0) pi_by_source_family$source_family[1] else NA_character_
  ),
  paste0(
    "Active pi in top source family: ",
    if (nrow(pi_by_source_family) > 0) sprintf("%.4f", pi_by_source_family$active_pi[1]) else NA_character_
  ),
  paste0(
    "Top PCA covariance by pi: ",
    if (nrow(pc_pi_dt[order(-pi)]) > 0) pc_pi_dt[order(-pi)]$pc[1] else NA_character_
  ),
  paste0(
    "Top PCA pi: ",
    if (nrow(pc_pi_dt[order(-pi)]) > 0) sprintf("%.4f", pc_pi_dt[order(-pi)]$pi[1]) else NA_character_
  )
)

if ("analyte__mRNA_vs_protein" %in% target_summary$contrast) {
  tmp = target_summary[contrast == "analyte__mRNA_vs_protein"]
  summary_lines = c(
    summary_lines,
    paste0(
      "Global mRNA_vs_protein: best |alignment| = ",
      sprintf("%.3f", tmp$best_abs_alignment),
      "; best weighted alignment = ",
      sprintf("%.4f", tmp$best_weighted_alignment)
    )
  )
}

if ("within_SC__mRNA_vs_protein" %in% target_summary$contrast) {
  tmp = target_summary[contrast == "within_SC__mRNA_vs_protein"]
  summary_lines = c(
    summary_lines,
    paste0(
      "Within-SC mRNA_vs_protein: best |alignment| = ",
      sprintf("%.3f", tmp$best_abs_alignment),
      "; best weighted alignment = ",
      sprintf("%.4f", tmp$best_weighted_alignment)
    )
  )
}

if ("within_mRNA__SC_vs_YNB2G" %in% target_summary$contrast) {
  tmp = target_summary[contrast == "within_mRNA__SC_vs_YNB2G"]
  summary_lines = c(
    summary_lines,
    paste0(
      "Within-mRNA SC_vs_YNB2G: best |alignment| = ",
      sprintf("%.3f", tmp$best_abs_alignment),
      "; best weighted alignment = ",
      sprintf("%.4f", tmp$best_weighted_alignment)
    )
  )
}

writeLines(summary_lines, con = outfile("plain_english_summary.txt"))

# ------------------------------------------------------------------------------
# console summary
# ------------------------------------------------------------------------------

message("\nTop covariances by pi:")
print(top_covariances_interpreted)

message("\nMixture weight by source family:")
print(pi_by_source_family)

message("\nTop contrasts per PCA covariance:")
print(top_by_pc)

if ("analyte__mRNA_vs_protein" %in% cov_alignment_dt$contrast) {
  message("\nTop covariance candidates for mRNA vs protein:")
  print(cov_alignment_dt[contrast == "analyte__mRNA_vs_protein"][order(-pi, -abs_alignment, -weighted_alignment, p_perm)])
}

if ("within_SC__mRNA_vs_protein" %in% cov_alignment_dt$contrast) {
  message("\nTop covariance candidates for mRNA vs protein within SC:")
  print(cov_alignment_dt[contrast == "within_SC__mRNA_vs_protein"][order(-pi, -abs_alignment, -weighted_alignment, p_perm)])
}

message("\nDone. Outputs written to: ", analysis_dir)