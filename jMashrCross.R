source("/PHShome/kv1026/p05/remap_functions.R")

valid_prefixes = c("fs", "ds", "na")
prefixes = c("ds")

if (length(prefixes) < 1L) {
  stop("prefixes must contain at least one of: ", paste(valid_prefixes, collapse = ", "))
}
if (anyDuplicated(prefixes)) {
  stop("prefixes contains duplicates: ", paste(unique(prefixes[duplicated(prefixes)]), collapse = ", "))
}
if (!all(prefixes %in% valid_prefixes)) {
  stop(
    "Invalid prefix value(s): ",
    paste(setdiff(prefixes, valid_prefixes), collapse = ", "),
    ". Choose any combination of: ",
    paste(valid_prefixes, collapse = ", ")
  )
}

alpha = 0.05
min_present = 1L
interval_sd_mult = 2
figure_width = 8.5
figure_height = 7.5
hard_xlim = c(-1.5, 1.5)

selected_pairs = tibble::tribble(
  ~gene,      ~marker,          ~example_class,
  "YML126C",  "chrXII:648649",  "Concordant",
  "YDR454C",  "chrXIV:465187",  "Discordant",
  "YML028W",  "chrXIV:465005",  "mRNA-specific",
  "YKL081W",  "chrXIII:48814",  "Protein-specific"
)

load_rda_object = function(path, object_name = NULL) {
  env = new.env(parent = emptyenv())
  nm = load(path, envir = env)

  if (!is.null(object_name)) {
    if (!(object_name %in% nm)) {
      stop("Object '", object_name, "' not found in ", path)
    }
    return(env[[object_name]])
  }

  if (length(nm) == 1L) {
    return(env[[nm]])
  }

  mget(nm, envir = env)
}

build_rows_for_pairs = function(pair_df, gene_fits, study_names, min_present = 1L) {
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
      Bhat = matrix(numeric(0), nrow = 0, ncol = length(study_names), dimnames = list(NULL, study_names)),
      Shat = matrix(numeric(0), nrow = 0, ncol = length(study_names), dimnames = list(NULL, study_names))
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

build_selected_qtl_mashr_estimates = function(selected_pairs, prefix, imed, min_present = 1L) {
  train_dir = file.path(imed, "mashr", paste0(prefix, "_training"))

  global_mash = readRDS(file.path(train_dir, paste0(prefix, "_global_mash_model.rds")))
  metadata_qtl = readRDS(file.path(train_dir, paste0(prefix, "_metadata_qtl.rds")))
  Vhat = readRDS(file.path(train_dir, paste0(prefix, "_null_correlation.rds")))
  gene_fits = load_rda_object(
    file.path(train_dir, paste0(prefix, "_gene_fits_supported_analyte.rda")),
    "gene_fits"
  )

  study_names = metadata_qtl$study

  mats = build_rows_for_pairs(
    pair_df = selected_pairs,
    gene_fits = gene_fits,
    study_names = study_names,
    min_present = min_present
  )

  if (nrow(mats$Bhat) == 0L) {
    stop("No mashr rows could be built for ", prefix)
  }

  missing_pairs = dplyr::anti_join(selected_pairs, mats$pairs, by = c("gene", "marker"))
  if (nrow(missing_pairs) > 0L) {
    warning(
      "Skipping ",
      nrow(missing_pairs),
      " selected gene-marker pair(s) not found for ",
      prefix,
      ": ",
      paste(paste0(missing_pairs$gene, " / ", missing_pairs$marker), collapse = ", ")
    )
  }

  mash_out = compute_mash_outputs(
    Bhat = mats$Bhat,
    Shat = mats$Shat,
    global_mash = global_mash,
    Vhat = Vhat
  )

  pm_long = data.frame(
    gene = mats$pairs$gene,
    marker = mats$pairs$marker,
    mash_out$pm,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  psd_long = data.frame(
    gene = mats$pairs$gene,
    marker = mats$pairs$marker,
    mash_out$psd,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  lfsr_long = data.frame(
    gene = mats$pairs$gene,
    marker = mats$pairs$marker,
    mash_out$lfsr,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(pm_long)[-(1:2)] = paste0(colnames(pm_long)[-(1:2)], "_pm")
  colnames(psd_long)[-(1:2)] = paste0(colnames(psd_long)[-(1:2)], "_psd")
  colnames(lfsr_long)[-(1:2)] = paste0(colnames(lfsr_long)[-(1:2)], "_lfsr")

  selected_pairs %>%
    dplyr::select(gene, marker, example_class) %>%
    dplyr::inner_join(mats$pairs, by = c("gene", "marker")) %>%
    dplyr::left_join(pm_long, by = c("gene", "marker")) %>%
    dplyr::left_join(psd_long, by = c("gene", "marker")) %>%
    dplyr::left_join(lfsr_long, by = c("gene", "marker"))
}

get_study_order = function(metadata_qtl, exclude_studies = character(0)) {
  metadata_qtl %>%
    dplyr::filter(!(study %in% exclude_studies)) %>%
    dplyr::mutate(
      original_order = dplyr::row_number(),
      analyte = factor(analyte, levels = c("mRNA", "protein"))
    ) %>%
    dplyr::arrange(analyte, original_order) %>%
    dplyr::pull(study)
}

build_plot_df = function(
  one_row,
  metadata_qtl,
  study_order,
  alpha = 0.05,
  interval_sd_mult = 2,
  exclude_studies = character(0),
  hard_xlim = c(-1.5, 1.5)
) {
  study_names = metadata_qtl$study
  keep_studies = setdiff(study_names, exclude_studies)

  pm_cols = paste0(keep_studies, "_pm")
  psd_cols = paste0(keep_studies, "_psd")
  lfsr_cols = paste0(keep_studies, "_lfsr")

  missing_pm = setdiff(pm_cols, colnames(one_row))
  missing_psd = setdiff(psd_cols, colnames(one_row))
  if (length(missing_pm) > 0L || length(missing_psd) > 0L) {
    stop("Missing pm/psd columns needed for plotting: ", paste(c(missing_pm, missing_psd), collapse = ", "))
  }

  metadata_qtl %>%
    dplyr::filter(study %in% keep_studies) %>%
    dplyr::select(study, analyte) %>%
    dplyr::mutate(
      pm = as.numeric(one_row[1, pm_cols, drop = TRUE]),
      psd = as.numeric(one_row[1, psd_cols, drop = TRUE]),
      lfsr = as.numeric(one_row[1, lfsr_cols, drop = TRUE]),
      xmin_raw = pm - interval_sd_mult * psd,
      xmax_raw = pm + interval_sd_mult * psd,
      xmin = pmax(xmin_raw, hard_xlim[1]),
      xmax = pmin(xmax_raw, hard_xlim[2]),
      sig_label = dplyr::case_when(
        is.finite(lfsr) & lfsr < alpha ~ "lfsr < 0.05",
        TRUE ~ "lfsr ≥ 0.05"
      ),
      sig_label = factor(sig_label, levels = c("lfsr ≥ 0.05", "lfsr < 0.05")),
      study_label = factor(study, levels = rev(study_order))
    )
}

build_fourpanel_plot = function(
  selected_estimates,
  metadata_qtl,
  selected_pairs,
  exclude_studies = character(0),
  alpha = 0.05,
  interval_sd_mult = 2,
  hard_xlim = c(-1.5, 1.5)
) {
  study_order = get_study_order(
    metadata_qtl = metadata_qtl,
    exclude_studies = exclude_studies
  )

  plot_rows = vector("list", nrow(selected_pairs))

  for (i in seq_len(nrow(selected_pairs))) {
    one_row = selected_estimates %>%
      dplyr::filter(
        gene == selected_pairs$gene[i],
        marker == selected_pairs$marker[i]
      )

    if (nrow(one_row) != 1L) {
      stop(
        "Expected exactly one row for ",
        selected_pairs$gene[i],
        " / ",
        selected_pairs$marker[i],
        ", found ",
        nrow(one_row),
        "."
      )
    }

    plot_df = build_plot_df(
      one_row = one_row,
      metadata_qtl = metadata_qtl,
      study_order = study_order,
      alpha = alpha,
      interval_sd_mult = interval_sd_mult,
      exclude_studies = exclude_studies,
      hard_xlim = hard_xlim
    )

    plot_df =
      plot_df %>%
      dplyr::mutate(
        pair_label = factor(
          paste0(
            selected_pairs$example_class[i],
            ": ",
            one_row$gene[1],
            "-",
            one_row$marker[1]
          ),
          levels = paste0(
            selected_pairs$example_class,
            ": ",
            selected_pairs$gene,
            "-",
            selected_pairs$marker
          )
        )
      )

    plot_rows[[i]] = plot_df
  }

  combined_df = dplyr::bind_rows(plot_rows)

  ggplot(
    combined_df,
    aes(
      x = pm,
      y = study_label
    )
  ) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "black") +
    geom_segment(
      aes(x = xmin, xend = xmax, yend = study_label, color = analyte, alpha = sig_label),
      linewidth = 0.45,
      lineend = "butt",
      show.legend = FALSE
    ) +
    geom_point(
      aes(color = analyte, alpha = sig_label),
      shape = 16,
      size = 2.2,
      show.legend = TRUE
    ) +
    facet_wrap(
      ~ pair_label,
      ncol = 2,
      scales = "fixed"
    ) +
    scale_color_manual(
      name = "Analyte",
      values = colAnalyte,
      drop = FALSE,
      guide = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2.5))
    ) +
    scale_alpha_manual(
      name = "Significance",
      values = c("lfsr ≥ 0.05" = 0.40, "lfsr < 0.05" = 1),
      drop = FALSE,
      guide = guide_legend(order = 2, override.aes = list(color = "grey35", size = 2.5))
    ) +
    scale_x_continuous(
      limits = hard_xlim,
      breaks = seq(hard_xlim[1], hard_xlim[2], by = 0.5)
    ) +
    labs(
      x = "mash posterior mean; intervals show +/- 2 posterior SD",
      y = NULL
    ) +
    coord_cartesian(xlim = hard_xlim, clip = "off") +
    theme_classic(base_size = 10.5) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      strip.text = element_text(color = "black", face = "plain", size = 10),
      strip.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = grid::unit(0.10, "cm")
    )
}

for (prefix in prefixes) {
  train_dir = file.path(imed, "mashr", paste0(prefix, "_training"))
  out_dir = file.path(imed, "mashr", paste0(prefix, "_cross_analyses"))

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  metadata_qtl = readRDS(file.path(train_dir, paste0(prefix, "_metadata_qtl.rds")))

  selected_estimates = build_selected_qtl_mashr_estimates(
    selected_pairs = selected_pairs,
    prefix = prefix,
    imed = imed,
    min_present = min_present
  )

  combined_plot = build_fourpanel_plot(
    selected_estimates = selected_estimates,
    metadata_qtl = metadata_qtl,
    selected_pairs = selected_pairs,
    exclude_studies = character(0),
    alpha = alpha,
    interval_sd_mult = interval_sd_mult,
    hard_xlim = hard_xlim
  )

  out_stem = file.path(
    out_dir,
    paste0(prefix, "_selected_analyte_qtl_fourpanel")
  )

  ggplot2::ggsave(
    filename = paste0(out_stem, ".pdf"),
    plot = combined_plot,
    width = figure_width,
    height = figure_height,
    units = "in",
    device = grDevices::cairo_pdf
  )

  write.table(
    selected_estimates,
    file = paste0(out_stem, ".tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}