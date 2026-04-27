source("/PHShome/kv1026/p05/remap_functions.R")

load(paste0(imed, "qtls_aug_analyte_ixn.rda"))
load(paste0(imed, "qtls_aug_analyte_geno.rda"))

qtls_aug_merge =
  bind_rows(
    qtls_aug_analyte_ixn %>%
      mutate(
        bin = recode(
          bin,
          "line cross" = "analyte_discordant",
          "mRNA specific" = "mRNA_specific",
          "protein specific" = "protein_specific",
          "nd in singles" = "analyte_undetected",
          "same direction" = "analyte_concordant"
        )
      ) %>%
      select(gene, marker, flavor, bin) %>%
      mutate(partition = "analyte"),

    qtls_aug_analyte_geno %>%
      mutate(
        bin = recode(
          bin,
          "line cross" = "analyte_discordant",
          "mRNA specific" = "mRNA_specific",
          "protein specific" = "protein_specific",
          "nd in singles" = "analyte_undetected",
          "same direction" = "analyte_concordant"
        )
      ) %>%
      select(gene, marker, flavor, bin) %>%
      mutate(partition = "analyte")
  ) %>%
  mutate(
    chr = match(
      sub(":.*", "", marker),
      c(
        "chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII",
        "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"
      )
    ),
    pos = as.numeric(sub(".*:", "", marker))
  ) %>%
  select(gene, marker, chr, pos, partition, flavor, bin) %>%
  arrange(chr, pos)

saveRDS(qtls_aug_merge, paste0(imed, "mashr/qtls_aug_merge.rds"))

qtl_markers = unique(qtls_aug_merge$marker)

# choose any combination of: "fs", "ds", "na"
prefix = c("fs", "ds", "na")

valid_prefixes = c("fs", "ds", "na")
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

alpha = 0.05
min_present = 5L
ld_abs_cor_threshold = 0.8
calibration_iterations = 1000L
calibration_seed = 343L
n_cores = 48
# n_cores = max(1L, parallel::detectCores())

load_rda_object = function(path, object_name = NULL) {
  env = new.env(parent = emptyenv())
  nm = load(path, envir = env)

  if (!is.null(object_name)) {
    return(env[[object_name]])
  }

  env[[nm]]
}

build_pair_table = function(tbl, flag_name) {
  tbl %>%
    dplyr::select(gene, marker) %>%
    dplyr::distinct() %>%
    dplyr::mutate(!!flag_name := TRUE)
}

build_rows_for_pairs = function(pair_df, gene_fits, study_names, min_present = 5L) {
  out_rows = vector("list", nrow(pair_df))
  keep = rep(FALSE, nrow(pair_df))

  for (i in seq_len(nrow(pair_df))) {
    gene = pair_df$gene[i]
    marker = pair_df$marker[i]

    if (!(gene %in% names(gene_fits))) next
    fits = gene_fits[[gene]]
    if (is.null(fits) || !(marker %in% names(fits[[1]]$beta))) next

    beta_row = rep(NA_real_, length(study_names))
    se_row = rep(NA_real_, length(study_names))

    for (j in seq_along(study_names)) {
      st = study_names[j]
      beta_row[j] = fits[[st]]$beta[marker]
      se_row[j] = fits[[st]]$se[marker]
    }

    out_rows[[i]] = list(
      gene = gene,
      marker = marker,
      Bhat = beta_row,
      Shat = se_row
    )
    keep[i] = TRUE
  }

  pair_df = pair_df[keep, , drop = FALSE]
  out_rows = out_rows[keep]

  if (length(out_rows) == 0L) {
    return(list(
      pairs = pair_df[0, , drop = FALSE],
      Bhat = matrix(
        numeric(0),
        nrow = 0,
        ncol = length(study_names),
        dimnames = list(NULL, study_names)
      ),
      Shat = matrix(
        numeric(0),
        nrow = 0,
        ncol = length(study_names),
        dimnames = list(NULL, study_names)
      )
    ))
  }

  Bhat = do.call(rbind, lapply(out_rows, `[[`, "Bhat"))
  Shat = do.call(rbind, lapply(out_rows, `[[`, "Shat"))

  colnames(Bhat) = study_names
  colnames(Shat) = study_names
  rownames(Bhat) = paste(pair_df$gene, pair_df$marker, sep = "::")
  rownames(Shat) = rownames(Bhat)

  finite_counts = rowSums(is.finite(Bhat) & is.finite(Shat))
  keep2 = finite_counts >= min_present

  list(
    pairs = pair_df[keep2, , drop = FALSE],
    Bhat = Bhat[keep2, , drop = FALSE],
    Shat = Shat[keep2, , drop = FALSE]
  )
}

compute_mash_outputs = function(Bhat, Shat, global_mash, Vhat) {
  data0 = mash_set_data(Bhat, Shat, alpha = 1)
  data = mash_update_data(data0, V = Vhat)
  post = mash_compute_posterior_matrices(global_mash, data)

  rep_mash = global_mash
  rep_mash$result = post
  rep_mash$data = data

  list(
    pm = get_pm(rep_mash),
    psd = get_psd(rep_mash),
    lfsr = get_lfsr(rep_mash)
  )
}

classify_two_context_signal = function(
  pos_a,
  neg_a,
  pos_b,
  neg_b,
  a_specific_label,
  b_specific_label,
  shared_prefix
) {
  has_a = (pos_a + neg_a) > 0
  has_b = (pos_b + neg_b) > 0

  consistent_a = xor(pos_a > 0, neg_a > 0)
  consistent_b = xor(pos_b > 0, neg_b > 0)

  if (!has_a && !has_b) {
    return(paste0(shared_prefix, "_undetected"))
  }

  if (has_a && !has_b) {
    if (consistent_a) {
      return(a_specific_label)
    } else {
      return(paste0(shared_prefix, "_mixed"))
    }
  }

  if (!has_a && has_b) {
    if (consistent_b) {
      return(b_specific_label)
    } else {
      return(paste0(shared_prefix, "_mixed"))
    }
  }

  if (!consistent_a || !consistent_b) {
    return(paste0(shared_prefix, "_mixed"))
  }

  same_pos = pos_a > 0 && pos_b > 0
  same_neg = neg_a > 0 && neg_b > 0
  cross_opp = (pos_a > 0 && neg_b > 0) || (neg_a > 0 && pos_b > 0)

  if (same_pos || same_neg) {
    return(paste0(shared_prefix, "_concordant"))
  }

  if (cross_opp) {
    return(paste0(shared_prefix, "_discordant"))
  }

  paste0(shared_prefix, "_mixed")
}

summarise_pair_signal = function(posteriors, alpha = 0.05) {
  posteriors %>%
    mutate(
      is_strong = is.finite(lfsr) & lfsr < alpha,
      is_pos = is_strong & is.finite(pm) & pm > 0,
      is_neg = is_strong & is.finite(pm) & pm < 0
    ) %>%
    group_by(gene, marker) %>%
    summarise(
      n_pos_strong = sum(is_pos, na.rm = TRUE),
      n_neg_strong = sum(is_neg, na.rm = TRUE),
      has_signal = (n_pos_strong + n_neg_strong) > 0,
      has_consistent_signal = xor(n_pos_strong > 0, n_neg_strong > 0),
      .groups = "drop"
    )
}

build_final_summary_table = function(posteriors, alpha = 0.05) {
  posteriors = posteriors %>%
    mutate(
      is_strong = is.finite(lfsr) & lfsr < alpha,
      is_pos = is_strong & is.finite(pm) & pm > 0,
      is_neg = is_strong & is.finite(pm) & pm < 0
    )

  final_tbl =
    posteriors %>%
    group_by(gene, marker) %>%
    summarise(
      pos_analyte_protein_strong = sum(analyte == "protein" & is_pos, na.rm = TRUE),
      pos_analyte_mRNA_strong = sum(analyte == "mRNA" & is_pos, na.rm = TRUE),
      neg_analyte_protein_strong = sum(analyte == "protein" & is_neg, na.rm = TRUE),
      neg_analyte_mRNA_strong = sum(analyte == "mRNA" & is_neg, na.rm = TRUE),
      analyte_mRNA_has_signal =
        (pos_analyte_mRNA_strong + neg_analyte_mRNA_strong) > 0,
      analyte_protein_has_signal =
        (pos_analyte_protein_strong + neg_analyte_protein_strong) > 0,
      analyte_mRNA_consistent_direction =
        xor(pos_analyte_mRNA_strong > 0, neg_analyte_mRNA_strong > 0),
      analyte_protein_consistent_direction =
        xor(pos_analyte_protein_strong > 0, neg_analyte_protein_strong > 0),
      .groups = "drop"
    ) %>%
    mutate(
      analyte_bin = purrr::pmap_chr(
        list(
          pos_analyte_mRNA_strong,
          neg_analyte_mRNA_strong,
          pos_analyte_protein_strong,
          neg_analyte_protein_strong
        ),
        function(pos_mrna, neg_mrna, pos_protein, neg_protein) {
          classify_two_context_signal(
            pos_a = pos_mrna,
            neg_a = neg_mrna,
            pos_b = pos_protein,
            neg_b = neg_protein,
            a_specific_label = "mRNA_specific",
            b_specific_label = "protein_specific",
            shared_prefix = "analyte"
          )
        }
      )
    ) %>%
    arrange(gene, marker)

  final_tbl
}

draw_heatmap = function(mat, main, cex_axis = 0.8, cex_text = 0.65) {
  stopifnot(is.matrix(mat) || is.data.frame(mat))
  mat = as.matrix(mat)

  pal = rev(colorRampPalette(
    rev(RColorBrewer::brewer.pal(n = 9, name = "YlOrRd"))
  )(100))

  nr = nrow(mat)
  nc = ncol(mat)
  gap_col = "grey90"

  image(
    x = seq_len(nc),
    y = seq_len(nr),
    z = t(mat),
    col = pal,
    zlim = c(0, max(mat, na.rm = TRUE)),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main
  )

  axis(
    side = 1,
    at = seq_len(nc),
    labels = colnames(mat),
    las = 2,
    cex.axis = cex_axis
  )

  axis(
    side = 2,
    at = seq_len(nr),
    labels = rownames(mat),
    las = 2,
    cex.axis = cex_axis
  )

  abline(v = seq(0.5, nc + 0.5, by = 1), col = gap_col)
  abline(h = seq(0.5, nr + 0.5, by = 1), col = gap_col)

  text(
    x = rep(seq_len(nc), each = nr),
    y = rep(seq_len(nr), times = nc),
    labels = sprintf("%.2f", as.vector(mat)),
    col = gap_col,
    cex = cex_text
  )

  box()
}

save_heatmap_dual = function(mat, main, base_path, cex_axis = 0.8, cex_text = 0.65,
                             width = 7, height = 6) {

  pdf(paste0(base_path, ".pdf"), width = width, height = height)
  par(mar = c(10, 10, 4, 2))
  draw_heatmap(mat, main = main, cex_axis = cex_axis, cex_text = cex_text)
  dev.off()

  png(paste0(base_path, ".png"), width = width, height = height, units = "in", res = 300)
  par(mar = c(10, 10, 4, 2))
  draw_heatmap(mat, main = main, cex_axis = cex_axis, cex_text = cex_text)
  dev.off()
}

write_offdiag_summary_csv = function(mat, out_csv, matrix_name) {
  stopifnot(is.matrix(mat) || is.data.frame(mat))
  mat = as.matrix(mat)

  if (nrow(mat) != ncol(mat)) {
    stop("Expected a square matrix for off-diagonal summary: ", matrix_name)
  }

  offdiag_vals = mat[row(mat) != col(mat)]
  offdiag_vals = offdiag_vals[is.finite(offdiag_vals)]

  summary_tbl =
    data.frame(
      matrix = matrix_name,
      n_offdiag = length(offdiag_vals),
      min = if (length(offdiag_vals) > 0L) min(offdiag_vals) else NA_real_,
      max = if (length(offdiag_vals) > 0L) max(offdiag_vals) else NA_real_,
      median = if (length(offdiag_vals) > 0L) median(offdiag_vals) else NA_real_,
      mean = if (length(offdiag_vals) > 0L) mean(offdiag_vals) else NA_real_,
      stringsAsFactors = FALSE
    )

  write.csv(summary_tbl, file = out_csv, row.names = FALSE)

  invisible(summary_tbl)
}

shared_ld_reference = readRDS(paste0(imed, "mashr/shared_ld_reference.rds"))

candidate_nonsig = setdiff(
  shared_ld_reference$markers,
  intersect(qtl_markers, shared_ld_reference$markers)
)

linked_to_qtl =
  apply(
    abs(
      shared_ld_reference$ld_r[
        candidate_nonsig,
        intersect(qtl_markers, shared_ld_reference$markers),
        drop = FALSE
      ]
    ) >= ld_abs_cor_threshold,
    1,
    any
  )

nonsig_candidates = candidate_nonsig[!linked_to_qtl]

message("QTL markers in LD reference: ", length(intersect(qtl_markers, shared_ld_reference$markers)))
message(
  "Non-QTL candidates not linked to any QTL marker at |r| >= ",
  ld_abs_cor_threshold,
  ": ",
  length(nonsig_candidates)
)
message("Calibration iterations will use ", n_cores, " cores.")

run_one_prefix = function(prefix) {

  message("============================================================")
  message("Running mash strong-pair analysis for arm: ", prefix)
  message("============================================================")

  idir = paste0(imed, "mashr/", prefix, "_training/")
  odir = paste0(imed, "mashr/", prefix, "_sharing_analyses/")

  if (!dir.exists(odir)) dir.create(odir, recursive = TRUE)

  infile = function(filename) file.path(idir, paste0(prefix, "_", filename))
  outfile = function(filename) file.path(odir, paste0(prefix, "_", filename))

  pairs_analyte_geno = build_pair_table(qtls_aug_analyte_geno, "from_analyte_geno")
  pairs_analyte_ixn = build_pair_table(qtls_aug_analyte_ixn, "from_analyte_ixn")

  selected_pairs =
    list(
      pairs_analyte_geno,
      pairs_analyte_ixn
    ) %>%
    purrr::reduce(full_join, by = c("gene", "marker")) %>%
    dplyr::mutate(
      across(starts_with("from_"), ~ifelse(is.na(.x), FALSE, .x))
    ) %>%
    dplyr::mutate(
      n_source_flags = rowSums(dplyr::select(., starts_with("from_")) * 1L),
      source_class = dplyr::case_when(
        n_source_flags > 1 ~ "multi_source",
        from_analyte_geno ~ "analyte_geno_only",
        from_analyte_ixn ~ "analyte_ixn_only",
        TRUE ~ "none"
      )
    ) %>%
    dplyr::arrange(gene, marker)

  global_mash = readRDS(infile("global_mash_model.rds"))
  gene_fits = load_rda_object(infile("gene_fits_supported_analyte.rda"), "gene_fits")
  metadata_qtl = readRDS(infile("metadata_qtl.rds"))
  Vhat = readRDS(infile("null_correlation.rds"))
  study_names = metadata_qtl$study

  selected_mats = build_rows_for_pairs(
    pair_df = selected_pairs,
    gene_fits = gene_fits,
    study_names = study_names,
    min_present = min_present
  )

  selected_pairs = selected_mats$pairs
  Bhat = selected_mats$Bhat
  Shat = selected_mats$Shat

  mash_selected = compute_mash_outputs(Bhat, Shat, global_mash, Vhat)
  pm = mash_selected$pm
  psd = mash_selected$psd
  lfsr = mash_selected$lfsr

  metadata_qtl_flat = metadata_qtl %>% dplyr::select(where(function(x) !is.list(x)))

  posteriors =
    data.frame(
      gene = rep(selected_pairs$gene, times = ncol(pm)),
      marker = rep(selected_pairs$marker, times = ncol(pm)),
      study = rep(colnames(pm), each = nrow(pm)),
      pm = as.vector(pm),
      psd = as.vector(psd),
      lfsr = as.vector(lfsr),
      stringsAsFactors = FALSE
    ) %>%
    left_join(metadata_qtl_flat, by = "study") %>%
    left_join(selected_pairs, by = c("gene", "marker"))

  pair_signal_summary = summarise_pair_signal(posteriors, alpha = alpha)

  strong_pairs =
    posteriors %>%
    filter(is.finite(lfsr) & lfsr < alpha) %>%
    left_join(
      pair_signal_summary %>%
        select(gene, marker, has_consistent_signal),
      by = c("gene", "marker")
    ) %>%
    transmute(
      gene = gene,
      marker = marker,
      study = study,
      pm = pm,
      psd = psd,
      lfsr = lfsr,
      analyte = analyte,
      medium = medium,
      sample_size = n_segs,
      pair_has_consistent_signal = has_consistent_signal
    ) %>%
    arrange(gene, marker, lfsr, study)

  write.table(
    strong_pairs,
    file = outfile("strong_pairs.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # --------------------------------------------------------------------------
  # LD-AWARE CALIBRATION, NOW PARALLELIZED ACROSS ITERATIONS
  # --------------------------------------------------------------------------

  calibration_gene_pool = sort(unique(selected_pairs$gene))
  calibration_gene_pool = calibration_gene_pool[calibration_gene_pool %in% names(gene_fits)]

  calibration_qtl_markers = sort(intersect(unique(selected_pairs$marker), shared_ld_reference$markers))
  calibration_null_markers = sort(intersect(nonsig_candidates, shared_ld_reference$markers))

  calibration_null_gene_availability = setNames(
    vector("list", length(calibration_null_markers)),
    calibration_null_markers
  )

  for (marker_i in calibration_null_markers) {
    calibration_null_gene_availability[[marker_i]] =
      calibration_gene_pool[
        vapply(
          calibration_gene_pool,
          function(gene_i) marker_i %in% names(gene_fits[[gene_i]][[1]]$beta),
          logical(1)
        )
      ]
  }

  calibration_null_markers =
    calibration_null_markers[
      vapply(calibration_null_gene_availability, length, integer(1)) > 0L
    ]

  if (length(calibration_qtl_markers) > 0L && length(calibration_null_markers) > 0L) {

    calibration_one_iter = function(iter) {
      set.seed(calibration_seed + iter)

      qtl_order = sample(calibration_qtl_markers)
      qtl_keep = character(0)
      qtl_blocked = setNames(rep(FALSE, length(qtl_order)), qtl_order)

      for (marker_i in qtl_order) {
        if (qtl_blocked[[marker_i]]) next

        qtl_keep = c(qtl_keep, marker_i)

        linked =
          abs(shared_ld_reference$ld_r[marker_i, qtl_order]) >= ld_abs_cor_threshold
        linked[!is.finite(linked)] = FALSE

        qtl_blocked[qtl_order[linked]] = TRUE
      }

      null_order = sample(calibration_null_markers)
      null_keep = character(0)
      null_blocked = setNames(rep(FALSE, length(null_order)), null_order)

      for (marker_i in null_order) {
        if (null_blocked[[marker_i]]) next

        null_keep = c(null_keep, marker_i)

        linked =
          abs(shared_ld_reference$ld_r[marker_i, null_order]) >= ld_abs_cor_threshold
        linked[!is.finite(linked)] = FALSE

        null_blocked[null_order[linked]] = TRUE
      }

      qtl_pairs_iter =
        selected_pairs %>%
        filter(marker %in% qtl_keep) %>%
        select(gene, marker) %>%
        distinct()

      if (nrow(qtl_pairs_iter) == 0L || length(null_keep) == 0L) {
        return(NULL)
      }

      null_pairs_iter =
        data.frame(
          marker = null_keep,
          gene = vapply(
            null_keep,
            function(marker_i) sample(calibration_null_gene_availability[[marker_i]], size = 1),
            character(1)
          ),
          stringsAsFactors = FALSE
        ) %>%
        select(gene, marker)

      null_mats = build_rows_for_pairs(
        pair_df = null_pairs_iter,
        gene_fits = gene_fits,
        study_names = study_names,
        min_present = min_present
      )

      null_pairs_iter = null_mats$pairs
      null_Bhat = null_mats$Bhat
      null_Shat = null_mats$Shat

      if (nrow(null_Bhat) == 0L) {
        return(NULL)
      }

      null_post = compute_mash_outputs(null_Bhat, null_Shat, global_mash, Vhat)

      null_posteriors =
        data.frame(
          gene = rep(null_pairs_iter$gene, times = ncol(null_post$pm)),
          marker = rep(null_pairs_iter$marker, times = ncol(null_post$pm)),
          study = rep(colnames(null_post$pm), each = nrow(null_post$pm)),
          pm = as.vector(null_post$pm),
          lfsr = as.vector(null_post$lfsr),
          stringsAsFactors = FALSE
        )

      qtl_posteriors_iter =
        posteriors %>%
        semi_join(qtl_pairs_iter, by = c("gene", "marker")) %>%
        select(gene, marker, study, pm, lfsr)

      qtl_pair_signal =
        summarise_pair_signal(qtl_posteriors_iter, alpha = alpha)

      null_pair_signal =
        summarise_pair_signal(null_posteriors, alpha = alpha)

      qtl_marker_signal =
        qtl_pair_signal %>%
        group_by(marker) %>%
        summarise(
          any_signal = any(has_consistent_signal, na.rm = TRUE),
          .groups = "drop"
        )

      null_marker_signal =
        null_pair_signal %>%
        group_by(marker) %>%
        summarise(
          any_signal = any(has_consistent_signal, na.rm = TRUE),
          .groups = "drop"
        )

      bind_rows(
        data.frame(
          set = "qtl",
          unlinked_markers = length(unique(qtl_keep)),
          unlinked_pairs = nrow(qtl_pairs_iter),
          frac_signal_markers = mean(qtl_marker_signal$any_signal, na.rm = TRUE),
          frac_signal_pairs = mean(qtl_pair_signal$has_consistent_signal, na.rm = TRUE),
          stringsAsFactors = FALSE
        ),
        data.frame(
          set = "nonqtl",
          unlinked_markers = length(unique(null_keep)),
          unlinked_pairs = nrow(null_pairs_iter),
          frac_signal_markers = mean(null_marker_signal$any_signal, na.rm = TRUE),
          frac_signal_pairs = mean(null_pair_signal$has_consistent_signal, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      )
    }

    calibration_rows =
      parallel::mclapply(
        X = seq_len(calibration_iterations),
        FUN = calibration_one_iter,
        mc.cores = n_cores
      )

    calibration_rows = calibration_rows[!vapply(calibration_rows, is.null, logical(1))]

    if (length(calibration_rows) > 0L) {
      ld_calibration =
        bind_rows(calibration_rows) %>%
        group_by(set) %>%
        summarise(
          median_unlinked_markers = median(unlinked_markers, na.rm = TRUE),
          median_unlinked_pairs = median(unlinked_pairs, na.rm = TRUE),
          frac_signal_markers = median(frac_signal_markers, na.rm = TRUE),
          frac_signal_pairs = median(frac_signal_pairs, na.rm = TRUE),
          .groups = "drop"
        )

      write.table(
        ld_calibration,
        file = outfile("ld_calibration.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
  }

  # --------------------------------------------------------------------------
  # FINAL SUMMARY TABLE
  # --------------------------------------------------------------------------

  final_summary = build_final_summary_table(posteriors = posteriors, alpha = alpha)

  write.table(
    final_summary,
    file = outfile("final_strong_effect_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # --------------------------------------------------------------------------
  # STRONG EFFECT CROSS-COMPARISON HEATMAPS
  # --------------------------------------------------------------------------

  metadata_qtl$group = paste(metadata_qtl$analyte, metadata_qtl$medium, sep = "_")
  study_ids = metadata_qtl$study
  group_ids = unique(metadata_qtl$group)

  study_signal = lfsr < alpha
  study_signal[!is.finite(lfsr)] = FALSE

  study_pos = study_signal & is.finite(pm) & (pm > 0)
  study_neg = study_signal & is.finite(pm) & (pm < 0)

  group_pos = matrix(
    FALSE,
    nrow = nrow(lfsr),
    ncol = length(group_ids),
    dimnames = list(rownames(lfsr), group_ids)
  )

  group_neg = matrix(
    FALSE,
    nrow = nrow(lfsr),
    ncol = length(group_ids),
    dimnames = list(rownames(lfsr), group_ids)
  )

  for (g in seq_along(group_ids)) {
    idx = metadata_qtl$group == group_ids[g]

    if (sum(idx) == 1) {
      group_pos[, g] = study_pos[, idx]
      group_neg[, g] = study_neg[, idx]
    } else {
      group_pos[, g] = apply(study_pos[, idx, drop = FALSE], 1, any, na.rm = TRUE)
      group_neg[, g] = apply(study_neg[, idx, drop = FALSE], 1, any, na.rm = TRUE)
    }
  }

  group_signal = xor(group_pos, group_neg)

  study_overlap = matrix(
    NA_real_,
    nrow = length(study_ids),
    ncol = length(study_ids),
    dimnames = list(study_ids, study_ids)
  )

  for (i in seq_along(study_ids)) {
    for (j in seq_along(study_ids)) {
      study_overlap[i, j] = mean(study_signal[, i] & study_signal[, j], na.rm = TRUE)
    }
  }

  group_overlap = matrix(
    NA_real_,
    nrow = length(group_ids),
    ncol = length(group_ids),
    dimnames = list(group_ids, group_ids)
  )

  for (i in seq_along(group_ids)) {
    for (j in seq_along(group_ids)) {
      group_overlap[i, j] = mean(group_signal[, i] & group_signal[, j], na.rm = TRUE)
    }
  }

  study_sharing = matrix(
    NA_real_,
    nrow = length(study_ids),
    ncol = length(study_ids),
    dimnames = list(study_ids, study_ids)
  )

  for (i in seq_along(study_ids)) {
    Ai = study_signal[, i]
    denom = sum(Ai, na.rm = TRUE)

    for (j in seq_along(study_ids)) {
      Bj = study_signal[, j]
      study_sharing[i, j] =
        if (denom == 0) NA_real_ else sum(Ai & Bj, na.rm = TRUE) / denom
    }
  }

  group_sharing = matrix(
    NA_real_,
    nrow = length(group_ids),
    ncol = length(group_ids),
    dimnames = list(group_ids, group_ids)
  )

  for (i in seq_along(group_ids)) {
    Ai = group_signal[, i]
    denom = sum(Ai, na.rm = TRUE)

    for (j in seq_along(group_ids)) {
      Bj = group_signal[, j]
      group_sharing[i, j] =
        if (denom == 0) NA_real_ else sum(Ai & Bj, na.rm = TRUE) / denom
    }
  }

  group_ids_noEtOH = setdiff(group_ids, "mRNA_YNB1EtOH")
  group_signal_noEtOH = group_signal[, group_ids_noEtOH, drop = FALSE]

  group_sharing_noEtOH = matrix(
    NA_real_,
    nrow = length(group_ids_noEtOH),
    ncol = length(group_ids_noEtOH),
    dimnames = list(group_ids_noEtOH, group_ids_noEtOH)
  )

  for (i in seq_along(group_ids_noEtOH)) {
    Ai = group_signal_noEtOH[, i]
    denom = sum(Ai, na.rm = TRUE)

    for (j in seq_along(group_ids_noEtOH)) {
      Bj = group_signal_noEtOH[, j]
      group_sharing_noEtOH[i, j] =
        if (denom == 0) NA_real_ else sum(Ai & Bj, na.rm = TRUE) / denom
    }
  }

  write_offdiag_summary_csv(
    study_overlap,
    outfile("study_overlap_offdiag_summary.csv"),
    "study_overlap"
  )

  write_offdiag_summary_csv(
    study_sharing,
    outfile("study_sharing_offdiag_summary.csv"),
    "study_sharing"
  )

  write_offdiag_summary_csv(
    group_overlap,
    outfile("group_overlap_offdiag_summary.csv"),
    "group_overlap"
  )

  write_offdiag_summary_csv(
    group_sharing,
    outfile("group_sharing_offdiag_summary.csv"),
    "group_sharing"
  )

  write_offdiag_summary_csv(
    group_sharing_noEtOH,
    outfile("group_sharing_noEtOH_offdiag_summary.csv"),
    "group_sharing_noEtOH"
  )

  write.table(
    data.frame(study = rownames(study_overlap), study_overlap, check.names = FALSE),
    file = outfile("study_overlap_heatmap_matrix.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    data.frame(group = rownames(group_overlap), group_overlap, check.names = FALSE),
    file = outfile("group_overlap_heatmap_matrix.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    data.frame(study = rownames(study_sharing), study_sharing, check.names = FALSE),
    file = outfile("study_sharing_heatmap_matrix.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    data.frame(group = rownames(group_sharing), group_sharing, check.names = FALSE),
    file = outfile("group_sharing_heatmap_matrix.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    data.frame(group = rownames(group_sharing_noEtOH), group_sharing_noEtOH, check.names = FALSE),
    file = outfile("group_sharing_noEtOH_heatmap_matrix.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  save_heatmap_dual(
    study_overlap,
    main = "Study-level overlap",
    base_path = file.path(odir, paste0(prefix, "_study_overlap_heatmap")),
    cex_axis = 0.75,
    cex_text = 0.55,
    width = 8,
    height = 7
  )

  save_heatmap_dual(
    group_overlap,
    main = "Group-level overlap",
    base_path = file.path(odir, paste0(prefix, "_group_overlap_heatmap")),
    cex_axis = 0.9,
    cex_text = 0.8,
    width = 7,
    height = 6
  )

  save_heatmap_dual(
    study_sharing,
    main = "Study-level sharing of strong effects",
    base_path = file.path(odir, paste0(prefix, "_study_sharing_heatmap")),
    cex_axis = 0.75,
    cex_text = 0.55,
    width = 8,
    height = 7
  )

  save_heatmap_dual(
    group_sharing,
    main = "Group-level sharing of strong effects",
    base_path = file.path(odir, paste0(prefix, "_group_sharing_heatmap")),
    cex_axis = 0.9,
    cex_text = 0.8,
    width = 7,
    height = 6
  )

  save_heatmap_dual(
    group_sharing_noEtOH,
    main = "Group-level sharing of strong effects",
    base_path = file.path(odir, paste0(prefix, "_group_sharing_noEtOH_heatmap")),
    cex_axis = 0.9,
    cex_text = 0.8,
    width = 6,
    height = 5.5
  )

  invisible(
    list(
      posteriors = posteriors,
      strong_pairs = strong_pairs,
      final_summary = final_summary
    )
  )
}

results = lapply(prefix, run_one_prefix)
names(results) = prefix