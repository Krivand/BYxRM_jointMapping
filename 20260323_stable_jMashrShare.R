source("/PHShome/kv1026/p05/remap_functions.R")

load(paste0(imed, "qtls_aug_analyte_ixn.rda"))
load(paste0(imed, "qtls_aug_medium_ixn.rda"))
load(paste0(imed, "qtls_aug_analyte_geno.rda"))
load(paste0(imed, "qtls_aug_medium_geno.rda"))

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
eps = 1e-12

load_rda_object = function(path, object_name = NULL) {
  env = new.env(parent = emptyenv())
  nm = load(path, envir = env)

  if (!is.null(object_name)) {
    return(env[[object_name]])
  }

  env[[nm]]
}

run_one_prefix = function(prefix) {

  message("============================================================")
  message("Running mash strong-pair analysis for arm: ", prefix)
  message("============================================================")

  # ----------------------------------------------------------------------------
  # paths
  # ----------------------------------------------------------------------------
  # idir = stable input location for mash objects produced earlier
  # odir = this script's output location
  # ----------------------------------------------------------------------------

  idir = paste0(imed, "mashr/", prefix, "_training/")
  odir = paste0(imed, "mashr/", prefix, "_sharing_analyses/")

  if (!dir.exists(odir)) dir.create(odir, recursive = TRUE)

  infile = function(filename) {
    file.path(idir, paste0(prefix, "_", filename))
  }

  outfile = function(filename) {
    file.path(odir, paste0(prefix, "_", filename))
  }

  # ----------------------------------------------------------------------------
  # strong marker sets
  # ----------------------------------------------------------------------------

  analyte_geno_markers = qtls_aug_analyte_geno$marker %>% unique()
  medium_ixn_markers = qtls_aug_medium_ixn$marker %>% unique()
  strong_markers = c(analyte_geno_markers, medium_ixn_markers) %>% unique()

  length(analyte_geno_markers)
  length(medium_ixn_markers)
  length(strong_markers)

  # ----------------------------------------------------------------------------
  # selected strong gene-marker pairs
  # keep provenance, but collapse to unique gene-marker rows
  # ----------------------------------------------------------------------------

  pairs_analyte_geno = qtls_aug_analyte_geno %>%
    dplyr::select(gene, marker) %>%
    dplyr::distinct() %>%
    dplyr::mutate(from_analyte_geno = TRUE)

  pairs_medium_ixn = qtls_aug_medium_ixn %>%
    dplyr::select(gene, marker) %>%
    dplyr::distinct() %>%
    dplyr::mutate(from_medium_ixn = TRUE)

  selected_pairs = full_join(
    pairs_analyte_geno,
    pairs_medium_ixn,
    by = c("gene", "marker")
  ) %>%
    dplyr::mutate(
      from_analyte_geno = ifelse(is.na(from_analyte_geno), FALSE, from_analyte_geno),
      from_medium_ixn = ifelse(is.na(from_medium_ixn), FALSE, from_medium_ixn)
    ) %>%
    dplyr::arrange(gene, marker)

  write.table(
    selected_pairs,
    file = outfile("selected_strong_pairs.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # load existing mash fit + inputs
  # ----------------------------------------------------------------------------

  global_mash = readRDS(infile("global_mash_model.rds"))
  gene_fits = load_rda_object(
    infile("gene_fits_supported_analyte.rda"),
    "gene_fits"
  )
  metadata_qtl = readRDS(infile("metadata_qtl.rds"))
  Vhat = readRDS(infile("null_correlation.rds"))

  study_names = metadata_qtl$study

  # ----------------------------------------------------------------------------
  # build matrices for just the selected strong pairs
  # ----------------------------------------------------------------------------

  selected_rows = vector("list", nrow(selected_pairs))
  keep = rep(FALSE, nrow(selected_pairs))

  for (i in seq_len(nrow(selected_pairs))) {
    gene = selected_pairs$gene[i]
    marker = selected_pairs$marker[i]

    if (!(gene %in% names(gene_fits))) next

    fits = gene_fits[[gene]]
    if (is.null(fits)) next
    if (!(marker %in% names(fits[[1]]$beta))) next

    beta_row = rep(NA_real_, length(study_names))
    se_row = rep(NA_real_, length(study_names))

    for (j in seq_along(study_names)) {
      st = study_names[j]
      beta_row[j] = fits[[st]]$beta[marker]
      se_row[j] = fits[[st]]$se[marker]
    }

    selected_rows[[i]] = list(
      gene = gene,
      marker = marker,
      Bhat = beta_row,
      Shat = se_row
    )
    keep[i] = TRUE
  }

  selected_pairs = selected_pairs[keep, , drop = FALSE]
  selected_rows = selected_rows[keep]

  Bhat = do.call(rbind, lapply(selected_rows, `[[`, "Bhat"))
  Shat = do.call(rbind, lapply(selected_rows, `[[`, "Shat"))

  colnames(Bhat) = study_names
  colnames(Shat) = study_names
  rownames(Bhat) = paste(selected_pairs$gene, selected_pairs$marker, sep = "::")
  rownames(Shat) = rownames(Bhat)

  finite_counts = rowSums(is.finite(Bhat) & is.finite(Shat))
  keep2 = finite_counts >= 2

  selected_pairs = selected_pairs[keep2, , drop = FALSE]
  Bhat = Bhat[keep2, , drop = FALSE]
  Shat = Shat[keep2, , drop = FALSE]

  # ----------------------------------------------------------------------------
  # posterior summaries from existing mash fit
  # ----------------------------------------------------------------------------

  data0 = mash_set_data(Bhat, Shat, alpha = 1)
  data = mash_update_data(data0, V = Vhat)

  post = mash_compute_posterior_matrices(global_mash, data)

  rep_mash = global_mash
  rep_mash$result = post
  rep_mash$data = data

  pm = get_pm(rep_mash)
  psd = get_psd(rep_mash)
  lfsr = get_lfsr(rep_mash)

  # ----------------------------------------------------------------------------
  # long posterior table
  # ----------------------------------------------------------------------------

  metadata_qtl_flat = metadata_qtl %>%
    dplyr::select(where(function(x) !is.list(x)))

  posteriors = data.frame(
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

  vapply(posteriors, is.list, logical(1))
  vapply(metadata_qtl, is.list, logical(1))

  write.table(
    posteriors,
    file = outfile("selected_strong_posteriors.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # analyte angle: compare mRNA_SC vs protein_SC
  # ----------------------------------------------------------------------------

  mRNA_SC = metadata_qtl$study[
    metadata_qtl$analyte == "mRNA" &
    metadata_qtl$medium == "SC"
  ]

  protein_SC = metadata_qtl$study[
    metadata_qtl$analyte == "protein" &
    metadata_qtl$medium == "SC"
  ]

  summary_rows = vector("list", nrow(selected_pairs))

  for (i in seq_len(nrow(selected_pairs))) {
  l_mrna = lfsr[i, mRNA_SC, drop = TRUE]
  p_mrna = pm[i, mRNA_SC, drop = TRUE]

  l_prot = lfsr[i, protein_SC, drop = TRUE]
  p_prot = pm[i, protein_SC, drop = TRUE]

  mrna_signal = any(l_mrna < alpha, na.rm = TRUE)
  prot_signal = any(l_prot < alpha, na.rm = TRUE)

  mrna_best_idx = which.min(ifelse(is.finite(l_mrna), l_mrna, Inf))
  prot_best_idx = which.min(ifelse(is.finite(l_prot), l_prot, Inf))

  mrna_best_study = if (length(mRNA_SC) > 0) mRNA_SC[mrna_best_idx] else NA
  prot_best_study = if (length(protein_SC) > 0) protein_SC[prot_best_idx] else NA

  mrna_best_pm = if (length(p_mrna) > 0) p_mrna[mrna_best_idx] else NA_real_
  prot_best_pm = if (length(p_prot) > 0) p_prot[prot_best_idx] else NA_real_

  mrna_best_lfsr = if (length(l_mrna) > 0) l_mrna[mrna_best_idx] else NA_real_
  prot_best_lfsr = if (length(l_prot) > 0) l_prot[prot_best_idx] else NA_real_

  state = "neither"
  sharing = "not_shared"
  direction = NA_character_

  if (mrna_signal & !prot_signal) state = "mRNA_only"
  if (!mrna_signal & prot_signal) state = "protein_only"

  if (mrna_signal & prot_signal) {
    sharing = "shared"

    if (sign(mrna_best_pm) == sign(prot_best_pm)) {
      state = "shared_same_sign"
      direction = "same_sign"
    } else {
      state = "shared_opposite_sign"
      direction = "opposite_sign"
    }
  }

  log2_abs_ratio = NA_real_
  if (mrna_signal & prot_signal) {
    log2_abs_ratio = log2(
      pmax(abs(mrna_best_pm), eps) /
      pmax(abs(prot_best_pm), eps)
    )
  }

  summary_rows[[i]] = data.frame(
    gene = selected_pairs$gene[i],
    marker = selected_pairs$marker[i],
    from_analyte_geno = selected_pairs$from_analyte_geno[i],
    from_medium_ixn = selected_pairs$from_medium_ixn[i],
    mRNA_SC_signal = mrna_signal,
    protein_SC_signal = prot_signal,
    mRNA_SC_best_study = mrna_best_study,
    protein_SC_best_study = prot_best_study,
    mRNA_SC_best_pm = mrna_best_pm,
    protein_SC_best_pm = prot_best_pm,
    mRNA_SC_best_lfsr = mrna_best_lfsr,
    protein_SC_best_lfsr = prot_best_lfsr,
    analyte_sharing_SC = sharing,
    analyte_direction_SC = direction,
    analyte_state_SC = state,
    log2_abs_ratio_mRNA_over_protein = log2_abs_ratio,
    stringsAsFactors = FALSE
  )
}

  analyte_summary = do.call(rbind, summary_rows)

  write.table(
    analyte_summary,
    file = outfile("selected_strong_analyte_SC_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # compact counts for the story
  # ----------------------------------------------------------------------------

  headline = analyte_summary %>%
    count(analyte_state_SC, name = "n") %>%
    arrange(desc(n))

  write.table(
    headline,
    file = outfile("selected_strong_analyte_SC_headline.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  sharing_counts = analyte_summary %>%
    count(analyte_sharing_SC, name = "n") %>%
    mutate(frac = n / sum(n)) %>%
    arrange(desc(n))

  write.table(
    sharing_counts,
    file = outfile("selected_strong_analyte_SC_sharing.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  sharing_direction_counts = analyte_summary %>%
    filter(analyte_sharing_SC == "shared") %>%
    count(analyte_direction_SC, name = "n") %>%
    mutate(frac_within_shared = n / sum(n)) %>%
    arrange(desc(n))

  write.table(
    sharing_direction_counts,
    file = outfile("selected_strong_analyte_SC_sharing_direction.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  sharing_by_state = analyte_summary %>%
    count(analyte_sharing_SC, analyte_state_SC, name = "n") %>%
    arrange(analyte_sharing_SC, desc(n))

  write.table(
    sharing_by_state,
    file = outfile("selected_strong_analyte_SC_sharing_by_state.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  shared_hits = analyte_summary %>%
    filter(analyte_sharing_SC == "shared") %>%
    arrange(protein_SC_best_lfsr, mRNA_SC_best_lfsr)

  write.table(
    shared_hits,
    file = outfile("selected_strong_analyte_SC_shared.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  shared_same_sign_hits = analyte_summary %>%
    filter(analyte_state_SC == "shared_same_sign") %>%
    arrange(protein_SC_best_lfsr, mRNA_SC_best_lfsr)

  write.table(
    shared_same_sign_hits,
    file = outfile("selected_strong_analyte_SC_shared_same_sign.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  shared_opposite_sign_hits = analyte_summary %>%
    filter(analyte_state_SC == "shared_opposite_sign") %>%
    arrange(protein_SC_best_lfsr, mRNA_SC_best_lfsr)

  write.table(
    shared_opposite_sign_hits,
    file = outfile("selected_strong_analyte_SC_shared_opposite_sign.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # quick peeking + plots (non-interactive)
  # ----------------------------------------------------------------------------

  asys_dir = file.path(odir, "asys")
  if (!dir.exists(asys_dir)) dir.create(asys_dir, recursive = TRUE)

  # ----------------------------------------------------------------------------
  # console summaries
  # ----------------------------------------------------------------------------

  cat("\n====================\n")
  cat("ANALYTE STATE COUNTS\n")
  cat("====================\n")
  print(
    analyte_summary %>%
      count(analyte_state_SC, name = "n") %>%
      mutate(frac = n / sum(n)) %>%
      arrange(desc(n))
  )

  cat("\n====================\n")
  cat("SHARING COUNTS\n")
  cat("====================\n")
  print(
    analyte_summary %>%
      count(analyte_sharing_SC, name = "n") %>%
      mutate(frac = n / sum(n))
  )

  cat("\n====================\n")
  cat("DIRECTION (WITHIN SHARED)\n")
  cat("====================\n")
  print(
    analyte_summary %>%
      filter(analyte_sharing_SC == "shared") %>%
      count(analyte_direction_SC, name = "n") %>%
      mutate(frac = n / sum(n))
  )

  # ----------------------------------------------------------------------------
  # save compact TSV summaries
  # ----------------------------------------------------------------------------

  write.table(
    analyte_summary %>%
      count(analyte_state_SC, name = "n") %>%
      mutate(frac = n / sum(n)),
    file = file.path(asys_dir, paste0(prefix, "_barplot_state.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    analyte_summary %>%
      count(analyte_sharing_SC, name = "n") %>%
      mutate(frac = n / sum(n)),
    file = file.path(asys_dir, paste0(prefix, "_barplot_sharing.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    analyte_summary %>%
      filter(analyte_sharing_SC == "shared") %>%
      count(analyte_direction_SC, name = "n") %>%
      mutate(frac = n / sum(n)),
    file = file.path(asys_dir, paste0(prefix, "_barplot_direction.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # base R barplots (saved to file)
  # ----------------------------------------------------------------------------

  state_counts = analyte_summary %>%
    count(analyte_state_SC, name = "n") %>%
    arrange(desc(n))

  pdf(file.path(asys_dir, paste0(prefix, "_barplot_state.pdf")), width = 6, height = 4)
  barplot(
    height = state_counts$n,
    names.arg = state_counts$analyte_state_SC,
    las = 2,
    main = "Analyte states (mRNA_SC vs protein_SC)",
    ylab = "Count"
  )
  dev.off()

  sharing_counts = analyte_summary %>%
    count(analyte_sharing_SC, name = "n")

  pdf(file.path(asys_dir, paste0(prefix, "_barplot_sharing.pdf")), width = 5, height = 4)
  barplot(
    height = sharing_counts$n,
    names.arg = sharing_counts$analyte_sharing_SC,
    main = "Sharing (mRNA_SC vs protein_SC)",
    ylab = "Count"
  )
  dev.off()

  direction_counts = analyte_summary %>%
    filter(analyte_sharing_SC == "shared") %>%
    count(analyte_direction_SC, name = "n")

  pdf(file.path(asys_dir, paste0(prefix, "_barplot_direction.pdf")), width = 5, height = 4)
  barplot(
    height = direction_counts$n,
    names.arg = direction_counts$analyte_direction_SC,
    main = "Direction among shared loci",
    ylab = "Count"
  )
  dev.off()

  cat("\nPlots written to:", asys_dir, "\n")

  # ----------------------------------------------------------------------------
  # writing-ready summaries
  # ----------------------------------------------------------------------------

  analyte_summary_onepergene = analyte_summary %>%
    mutate(
      best_rank = pmin(mRNA_SC_best_lfsr, protein_SC_best_lfsr, na.rm = TRUE),
      best_rank = ifelse(is.finite(best_rank), best_rank, Inf)
    ) %>%
    arrange(gene, best_rank, desc(abs(mRNA_SC_best_pm) + abs(protein_SC_best_pm))) %>%
    group_by(gene) %>%
    slice(1) %>%
    ungroup()

  write.table(
    analyte_summary_onepergene,
    file = outfile("selected_strong_analyte_SC_summary_onepergene.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # effect-size summaries
  # ----------------------------------------------------------------------------

  effect_size_summary_pair = analyte_summary %>%
    mutate(
      abs_mRNA_SC_pm = abs(mRNA_SC_best_pm),
      abs_protein_SC_pm = abs(protein_SC_best_pm)
    ) %>%
    group_by(analyte_state_SC) %>%
    summarise(
      n = n(),
      median_abs_mRNA_SC_pm = median(abs_mRNA_SC_pm, na.rm = TRUE),
      median_abs_protein_SC_pm = median(abs_protein_SC_pm, na.rm = TRUE),
      median_abs_ratio_mRNA_over_protein = median(
        log2_abs_ratio_mRNA_over_protein,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    arrange(desc(n))

  effect_size_summary_gene = analyte_summary_onepergene %>%
    mutate(
      abs_mRNA_SC_pm = abs(mRNA_SC_best_pm),
      abs_protein_SC_pm = abs(protein_SC_best_pm)
    ) %>%
    group_by(analyte_state_SC) %>%
    summarise(
      n = n(),
      median_abs_mRNA_SC_pm = median(abs_mRNA_SC_pm, na.rm = TRUE),
      median_abs_protein_SC_pm = median(abs_protein_SC_pm, na.rm = TRUE),
      median_abs_ratio_mRNA_over_protein = median(
        log2_abs_ratio_mRNA_over_protein,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    arrange(desc(n))

  write.table(
    effect_size_summary_pair,
    file = file.path(asys_dir, paste0(prefix, "_writing_effect_size_summary_pair.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    effect_size_summary_gene,
    file = file.path(asys_dir, paste0(prefix, "_writing_effect_size_summary_onepergene.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # simple within-analyte consistency checks across studies
  # ----------------------------------------------------------------------------

  consistency_rows = vector("list", nrow(selected_pairs))

  for (i in seq_len(nrow(selected_pairs))) {
    pm_mrna = pm[i, mRNA_SC, drop = TRUE]
    lfsr_mrna = lfsr[i, mRNA_SC, drop = TRUE]

    pm_prot = pm[i, protein_SC, drop = TRUE]
    lfsr_prot = lfsr[i, protein_SC, drop = TRUE]

    mrna_supported = is.finite(lfsr_mrna) & (lfsr_mrna < alpha)
    prot_supported = is.finite(lfsr_prot) & (lfsr_prot < alpha)

    n_mrna_supported = sum(mrna_supported, na.rm = TRUE)
    n_prot_supported = sum(prot_supported, na.rm = TRUE)

    mrna_same_sign_supported = NA
    prot_same_sign_supported = NA

    if (n_mrna_supported > 0) {
      mrna_signs = sign(pm_mrna[mrna_supported])
      mrna_same_sign_supported = length(unique(mrna_signs)) == 1
    }

    if (n_prot_supported > 0) {
      prot_signs = sign(pm_prot[prot_supported])
      prot_same_sign_supported = length(unique(prot_signs)) == 1
    }

    consistency_rows[[i]] = data.frame(
      gene = selected_pairs$gene[i],
      marker = selected_pairs$marker[i],
      n_mRNA_SC_supported = n_mrna_supported,
      n_protein_SC_supported = n_prot_supported,
      mRNA_SC_same_sign_supported = mrna_same_sign_supported,
      protein_SC_same_sign_supported = prot_same_sign_supported,
      both_analytes_ge2_supported = (n_mrna_supported >= 2) & (n_prot_supported >= 2),
      stringsAsFactors = FALSE
    )
  }

  consistency_summary = do.call(rbind, consistency_rows)

  analyte_summary = analyte_summary %>%
    left_join(consistency_summary, by = c("gene", "marker"))

  analyte_summary_onepergene = analyte_summary_onepergene %>%
    left_join(consistency_summary, by = c("gene", "marker"))

  write.table(
    analyte_summary,
    file = outfile("selected_strong_analyte_SC_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    analyte_summary_onepergene,
    file = outfile("selected_strong_analyte_SC_summary_onepergene.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  consistency_pair = analyte_summary %>%
    summarise(
      n_total = n(),
      n_shared = sum(analyte_sharing_SC == "shared"),
      frac_shared = mean(analyte_sharing_SC == "shared"),
      n_shared_ge2_each = sum(
        analyte_sharing_SC == "shared" & both_analytes_ge2_supported,
        na.rm = TRUE
      ),
      frac_shared_ge2_each = mean(
        both_analytes_ge2_supported[analyte_sharing_SC == "shared"],
        na.rm = TRUE
      )
    )

  consistency_gene = analyte_summary_onepergene %>%
    summarise(
      n_total = n(),
      n_shared = sum(analyte_sharing_SC == "shared"),
      frac_shared = mean(analyte_sharing_SC == "shared"),
      n_shared_ge2_each = sum(
        analyte_sharing_SC == "shared" & both_analytes_ge2_supported,
        na.rm = TRUE
      ),
      frac_shared_ge2_each = mean(
        both_analytes_ge2_supported[analyte_sharing_SC == "shared"],
        na.rm = TRUE
      )
    )

  write.table(
    consistency_pair,
    file = file.path(asys_dir, paste0(prefix, "_writing_consistency_pair.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    consistency_gene,
    file = file.path(asys_dir, paste0(prefix, "_writing_consistency_onepergene.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # focused tables for writing
  # ----------------------------------------------------------------------------

  opposite_sign_hits_onepergene = analyte_summary_onepergene %>%
    filter(analyte_state_SC == "shared_opposite_sign") %>%
    arrange(protein_SC_best_lfsr, mRNA_SC_best_lfsr)

  shared_same_sign_hits_onepergene = analyte_summary_onepergene %>%
    filter(analyte_state_SC == "shared_same_sign") %>%
    arrange(protein_SC_best_lfsr, mRNA_SC_best_lfsr)

  write.table(
    opposite_sign_hits_onepergene,
    file = file.path(asys_dir, paste0(prefix, "_writing_opposite_sign_hits_onepergene.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    shared_same_sign_hits_onepergene,
    file = file.path(asys_dir, paste0(prefix, "_writing_shared_same_sign_hits_onepergene.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # console output for section writing
  # ----------------------------------------------------------------------------

  cat("\n====================\n")
  cat("PAIR-LEVEL TOTALS\n")
  cat("====================\n")
  cat("Selected gene-marker pairs:", nrow(analyte_summary), "\n")
  cat("Unique genes in selected set:", dplyr::n_distinct(analyte_summary$gene), "\n")

  cat("\n====================\n")
  cat("PAIR-LEVEL STATE COUNTS\n")
  cat("====================\n")
  print(
    analyte_summary %>%
      count(analyte_state_SC, name = "n") %>%
      mutate(frac = n / sum(n)) %>%
      arrange(desc(n))
  )

  cat("\n====================\n")
  cat("PAIR-LEVEL SHARING\n")
  cat("====================\n")
  print(
    analyte_summary %>%
      count(analyte_sharing_SC, name = "n") %>%
      mutate(frac = n / sum(n)) %>%
      arrange(desc(n))
  )

  cat("\n====================\n")
  cat("PAIR-LEVEL DIRECTION WITHIN SHARED\n")
  cat("====================\n")
  print(
    analyte_summary %>%
      filter(analyte_sharing_SC == "shared") %>%
      count(analyte_direction_SC, name = "n") %>%
      mutate(frac = n / sum(n)) %>%
      arrange(desc(n))
  )

  cat("\n====================\n")
  cat("ONE-PER-GENE TOTALS\n")
  cat("====================\n")
  cat("Representative genes:", nrow(analyte_summary_onepergene), "\n")

  cat("\n====================\n")
  cat("ONE-PER-GENE STATE COUNTS\n")
  cat("====================\n")
  print(
    analyte_summary_onepergene %>%
      count(analyte_state_SC, name = "n") %>%
      mutate(frac = n / sum(n)) %>%
      arrange(desc(n))
  )

  cat("\n====================\n")
  cat("ONE-PER-GENE SHARING\n")
  cat("====================\n")
  print(
    analyte_summary_onepergene %>%
      count(analyte_sharing_SC, name = "n") %>%
      mutate(frac = n / sum(n)) %>%
      arrange(desc(n))
  )

  cat("\n====================\n")
  cat("ONE-PER-GENE DIRECTION WITHIN SHARED\n")
  cat("====================\n")
  print(
    analyte_summary_onepergene %>%
      filter(analyte_sharing_SC == "shared") %>%
      count(analyte_direction_SC, name = "n") %>%
      mutate(frac = n / sum(n)) %>%
      arrange(desc(n))
  )

  cat("\n====================\n")
  cat("PAIR-LEVEL EFFECT SIZE SUMMARY\n")
  cat("====================\n")
  print(effect_size_summary_pair)

  cat("\n====================\n")
  cat("ONE-PER-GENE EFFECT SIZE SUMMARY\n")
  cat("====================\n")
  print(effect_size_summary_gene)

  cat("\n====================\n")
  cat("PAIR-LEVEL CONSISTENCY SUMMARY\n")
  cat("====================\n")
  print(consistency_pair)

  cat("\n====================\n")
  cat("ONE-PER-GENE CONSISTENCY SUMMARY\n")
  cat("====================\n")
  print(consistency_gene)

  cat("\n====================\n")
  cat("TOP OPPOSITE-SIGN HITS (PAIR-LEVEL)\n")
  cat("====================\n")
  print(
    analyte_summary %>%
      filter(analyte_state_SC == "shared_opposite_sign") %>%
      select(
        gene, marker,
        mRNA_SC_best_study, protein_SC_best_study,
        mRNA_SC_best_pm, protein_SC_best_pm,
        mRNA_SC_best_lfsr, protein_SC_best_lfsr
      ) %>%
      arrange(protein_SC_best_lfsr, mRNA_SC_best_lfsr)
  )

  cat("\n====================\n")
  cat("TOP OPPOSITE-SIGN HITS (ONE-PER-GENE)\n")
  cat("====================\n")
  print(
    opposite_sign_hits_onepergene %>%
      select(
        gene, marker,
        mRNA_SC_best_study, protein_SC_best_study,
        mRNA_SC_best_pm, protein_SC_best_pm,
        mRNA_SC_best_lfsr, protein_SC_best_lfsr
      )
  )

  cat("\nWriting tables saved in:\n", asys_dir, "\n")

  # ----------------------------------------------------------------------------
  # figure: effect overlap and sharing heatmaps
  # ----------------------------------------------------------------------------

  asys_dir = file.path(odir, "asys")
  if (!dir.exists(asys_dir)) dir.create(asys_dir, recursive = TRUE)

  metadata_qtl$group = paste(metadata_qtl$analyte, metadata_qtl$medium, sep = "_")

  study_ids = metadata_qtl$study
  group_ids = unique(metadata_qtl$group)

  study_signal = lfsr < alpha
  study_signal[!is.finite(lfsr)] = FALSE

  group_signal = matrix(
    FALSE,
    nrow = nrow(lfsr),
    ncol = length(group_ids),
    dimnames = list(rownames(lfsr), group_ids)
  )

  for (g in seq_along(group_ids)) {
    idx = metadata_qtl$group == group_ids[g]
    if (sum(idx) == 1) {
      group_signal[, g] = study_signal[, idx]
    } else {
      group_signal[, g] = apply(study_signal[, idx, drop = FALSE], 1, any, na.rm = TRUE)
    }
  }

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

      if (denom == 0) {
        study_sharing[i, j] = NA_real_
      } else {
        study_sharing[i, j] = sum(Ai & Bj, na.rm = TRUE) / denom
      }
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

      if (denom == 0) {
        group_sharing[i, j] = NA_real_
      } else {
        group_sharing[i, j] = sum(Ai & Bj, na.rm = TRUE) / denom
      }
    }
  }

  write.table(
    data.frame(study = rownames(study_overlap), study_overlap, check.names = FALSE),
    file = file.path(asys_dir, paste0(prefix, "_study_overlap_heatmap_matrix.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    data.frame(group = rownames(group_overlap), group_overlap, check.names = FALSE),
    file = file.path(asys_dir, paste0(prefix, "_group_overlap_heatmap_matrix.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    data.frame(study = rownames(study_sharing), study_sharing, check.names = FALSE),
    file = file.path(asys_dir, paste0(prefix, "_study_sharing_heatmap_matrix.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    data.frame(group = rownames(group_sharing), group_sharing, check.names = FALSE),
    file = file.path(asys_dir, paste0(prefix, "_group_sharing_heatmap_matrix.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("\n====================\n")
  cat("STUDY-LEVEL OVERLAP MATRIX\n")
  cat("====================\n")
  print(round(study_overlap, 3))

  cat("\n====================\n")
  cat("GROUP-LEVEL OVERLAP MATRIX\n")
  cat("====================\n")
  print(round(group_overlap, 3))

  cat("\n====================\n")
  cat("STUDY-LEVEL SHARING MATRIX\n")
  cat("====================\n")
  print(round(study_sharing, 3))

  cat("\n====================\n")
  cat("GROUP-LEVEL SHARING MATRIX\n")
  cat("====================\n")
  print(round(group_sharing, 3))

  draw_heatmap = function(mat, main, cex_axis = 0.8, cex_text = 0.65) {
  stopifnot(is.matrix(mat) || is.data.frame(mat))
  mat = as.matrix(mat)

  pal = viridis::viridis(100)

  nr = nrow(mat)
  nc = ncol(mat)

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

  abline(v = seq(0.5, nc + 0.5, by = 1), col = "grey90")
  abline(h = seq(0.5, nr + 0.5, by = 1), col = "grey90")

  text(
    x = rep(seq_len(nc), each = nr),
    y = rep(seq_len(nr), times = nc),
    labels = sprintf("%.2f", as.vector(mat)),
    cex = cex_text
  )

  box()
  }

  pdf(
    file.path(asys_dir, paste0(prefix, "_overlap_heatmaps.pdf")),
    width = 14,
    height = 7
  )

  layout(matrix(c(1, 2), nrow = 1))

  par(mar = c(10, 10, 4, 2))
  draw_heatmap(
    study_overlap,
    main = "Study-level overlap\nfraction of selected loci with signal in both studies",
    cex_axis = 0.75,
    cex_text = 0.55
  )

  par(mar = c(10, 10, 4, 2))
  draw_heatmap(
    group_overlap,
    main = "Group-level overlap\nfraction of selected loci with signal in both groups",
    cex_axis = 0.9,
    cex_text = 0.8
  )

  dev.off()

  cat("\nOverlap heatmap figure written to:\n")
  cat(file.path(asys_dir, paste0(prefix, "_overlap_heatmaps.pdf")), "\n")

  pdf(
    file.path(asys_dir, paste0(prefix, "_sharing_heatmaps.pdf")),
    width = 14,
    height = 7
  )

  layout(matrix(c(1, 2), nrow = 1))

  par(mar = c(10, 10, 4, 2))
  draw_heatmap(
    study_sharing,
    main = "Study-level sharing\nfraction of row-study loci also significant in column-study",
    cex_axis = 0.75,
    cex_text = 0.55
  )

  par(mar = c(10, 10, 4, 2))
  draw_heatmap(
    group_sharing,
    main = "Group-level sharing\nfraction of row-group loci also significant in column-group",
    cex_axis = 0.9,
    cex_text = 0.8
  )

  dev.off()

  cat("\nSharing heatmap figure written to:\n")
  cat(file.path(asys_dir, paste0(prefix, "_sharing_heatmaps.pdf")), "\n")

  message("Done with arm: ", prefix)

  invisible(list(
    prefix = prefix,
    idir = idir,
    odir = odir,
    asys_dir = asys_dir
  ))
}

analysis_results = list()
for (px in prefix) {
  analysis_results[[px]] = run_one_prefix(px)
}

message("Done.")