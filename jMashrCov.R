# jMashrCov.R
# Build pi-weighted mash covariance summaries and Urbut-style figures
# for any combination of mashr arms (fs, ds, na)

source("/PHShome/kv1026/p05/remap_functions.R")

# choose any combination of arms to run
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

write_tsv_if_needed = function(dt, path, overwrite = TRUE, sep = "\t") {
  if (!file.exists(path) || isTRUE(overwrite)) {
    data.table::fwrite(dt, path, sep = sep)
  }
  invisible(path)
}

symmetrize_matrix = function(mat) {
  mat = as.matrix(mat)
  0.5 * (mat + t(mat))
}

safe_cor_from_cov = function(U) {
  U = symmetrize_matrix(U)

  d = diag(U)
  d[!is.finite(d)] = NA_real_

  if (any(d < 0, na.rm = TRUE)) {
    warning("Negative diagonal entries detected in covariance matrix. Setting tiny floor before cov2cor.")
    d[d < 0] = 0
    diag(U) = d
  }

  eps = 1e-12
  zero_diag = is.na(diag(U)) | diag(U) <= eps
  if (any(zero_diag)) {
    diag(U)[zero_diag] = eps
  }

  cor_mat = stats::cov2cor(U)
  cor_mat = symmetrize_matrix(cor_mat)
  cor_mat
}

scale_vec_for_plot = function(x) {
  mx = max(abs(x), na.rm = TRUE)
  if (!is.finite(mx) || mx == 0) {
    return(x)
  }
  x / mx
}

orient_pc_for_plot = function(x) {
  x = scale_vec_for_plot(x)

  if (all(!is.finite(x))) {
    return(x)
  }

  tallest_idx = which.max(abs(x))
  if (length(tallest_idx) == 1L && is.finite(x[tallest_idx]) && x[tallest_idx] < 0) {
    x = -x
  }

  x
}

base_plot_theme = function() {
  theme_bw(base_size = 11) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(hjust = 0, face = "plain", color = "black"),
      legend.title = element_text(color = "black"),
      legend.text = element_text(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(color = "black"),
      plot.background = element_blank()
    )
}

build_heatmap_df = function(mat, clustering = NULL, upper_only = TRUE) {
  mat = as.matrix(mat)
  stopifnot(nrow(mat) == ncol(mat))

  ord =
    if (!is.null(clustering)) {
      clustering$order
    } else {
      seq_len(nrow(mat))
    }

  mat = mat[ord, ord, drop = FALSE]
  labs = colnames(mat)

  df = expand.grid(
    row = seq_len(nrow(mat)),
    col = seq_len(ncol(mat)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      study_row = labs[row],
      study_col = labs[col],
      value = mat[cbind(row, col)]
    )

  if (isTRUE(upper_only)) {
    df = df %>% filter(col >= row)
  }

  df %>%
    mutate(
      study_row = factor(study_row, levels = rev(labs)),
      study_col = factor(study_col, levels = labs),
      label = sprintf("%.2f", value),
      text_col = "white"
    )
}

make_heatmap_plot = function(mat, main, clustering = NULL) {
  df = build_heatmap_df(mat = mat, clustering = clustering, upper_only = TRUE)
  pal = rev(colorRampPalette(rev(RColorBrewer::brewer.pal(9, "YlOrRd")))(100))

  ggplot(df, aes(x = study_col, y = study_row, fill = value)) +
    geom_tile(color = "grey90", linewidth = 0.35) +
    geom_text(aes(label = label, color = text_col), size = 4) +
    scale_color_identity() +
    scale_fill_gradientn(
      colours = pal,
      limits = c(0, 1),
      name = NULL,
      guide = guide_colorbar(
        barheight = unit(6, "cm"),
        barwidth = unit(0.6, "cm"),
        frame.colour = NA,
        ticks.colour = "black"
      )
    ) +
    coord_fixed() +
    labs(title = main, x = NULL, y = NULL) +
    base_plot_theme() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      legend.position = "right"
    )
}

make_pc_barplot = function(eig, pc_index, study_order = NULL, show_x_labels = FALSE) {
  vals = eig$vectors[, pc_index]

  if (!is.null(study_order)) {
    vals = vals[study_order]
    labs = rownames(eig$vectors)[study_order]
  } else {
    labs = rownames(eig$vectors)
  }

  vals = orient_pc_for_plot(vals)
  pve = eig$values / sum(eig$values)

  df = data.frame(
    study = factor(labs, levels = labs),
    loading = vals,
    stringsAsFactors = FALSE
  )

  p =
    ggplot(df, aes(x = study, y = loading)) +
    geom_col() +
    labs(
      title = paste0("PC", pc_index, " (PVE = ", sprintf("%.1f", 100 * pve[pc_index]), "%)"),
      x = "Study",
      y = "Scaled loading"
    ) +
    base_plot_theme() +
    theme(
      axis.text.y = element_text(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      axis.line.x = element_line(color = "black", linewidth = 0.4),
      axis.line.y = element_line(color = "black", linewidth = 0.4),
      panel.border = element_blank()
    )

  if (isTRUE(show_x_labels)) {
    p = p +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black")
      )
  } else {
    p = p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
      )
  }

  p
}

save_all_pc_barplots = function(eig, outdir, study_order = NULL, width = 10, height_pdf = 5, height_png = 7) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  n_pcs = ncol(eig$vectors)

  for (pc_index in seq_len(n_pcs)) {
    p_pdf = make_pc_barplot(
      eig = eig,
      pc_index = pc_index,
      study_order = study_order,
      show_x_labels = TRUE
    )

    p_png =
      p_pdf +
      theme(
        plot.margin = margin(t = 16, r = 16, b = 80, l = 40),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          color = "black",
          size = 10
        )
      ) +
      coord_cartesian(clip = "off")

    pdf_file = file.path(outdir, paste0(sprintf("PC%02d", pc_index), "_barplot.pdf"))
    png_file = file.path(outdir, paste0(sprintf("PC%02d", pc_index), "_barplot.png"))

    ggplot2::ggsave(
      filename = pdf_file,
      plot = p_pdf,
      width = width,
      height = height_pdf,
      units = "in",
      device = cairo_pdf
    )

    ggplot2::ggsave(
      filename = png_file,
      plot = p_png,
      width = width,
      height = height_png,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }
}

plot_figure = function(cor_mat, eig, main_title, pdf_file, png_file, clustering = NULL) {
  ord =
    if (!is.null(clustering)) {
      clustering$order
    } else {
      seq_len(nrow(cor_mat))
    }

  p_heat = make_heatmap_plot(
    mat = cor_mat,
    main = main_title,
    clustering = clustering
  )

  p1 = make_pc_barplot(
    eig = eig,
    pc_index = 1,
    study_order = ord,
    show_x_labels = FALSE
  )

  p2 = make_pc_barplot(
    eig = eig,
    pc_index = 2,
    study_order = ord,
    show_x_labels = FALSE
  )

  p3 = make_pc_barplot(
    eig = eig,
    pc_index = 3,
    study_order = ord,
    show_x_labels = TRUE
  )

  right_col = p1 / p2 / p3
  full_plot = p_heat | right_col
  full_plot = full_plot + patchwork::plot_layout(widths = c(1.8, 1))

  ggplot2::ggsave(
    filename = pdf_file,
    plot = full_plot,
    width = 16,
    height = 10.5,
    units = "in",
    device = cairo_pdf
  )

  ggplot2::ggsave(
    filename = png_file,
    plot = full_plot,
    width = 16,
    height = 10.5,
    units = "in",
    dpi = 300
  )
}

run_one_prefix = function(prefix) {
  idir = file.path(imed, "mashr", paste0(prefix, "_training"))
  odir = file.path(imed, "mashr", paste0(prefix, "_cov"))

  if (!dir.exists(odir)) {
    dir.create(odir, recursive = TRUE)
  }

  infile = function(filename) file.path(idir, paste0(prefix, "_", filename))
  outfile = function(filename) file.path(odir, paste0(prefix, "_", filename))

  message("============================================================")
  message("Running covariance summary for arm: ", prefix)
  message("Input dir: ", idir)
  message("Output dir: ", odir)
  message("============================================================")

  global_mash_path = infile("global_mash_model.rds")
  metadata_qtl_path = infile("metadata_qtl.rds")

  if (!file.exists(global_mash_path)) {
    stop("Missing mash model for prefix ", prefix, ": ", global_mash_path)
  }
  if (!file.exists(metadata_qtl_path)) {
    stop("Missing metadata_qtl for prefix ", prefix, ": ", metadata_qtl_path)
  }

  message("Loading ", prefix, " mash objects...")
  global_mash = readRDS(global_mash_path)
  metadata_qtl = readRDS(metadata_qtl_path)

  study_names = metadata_qtl$study
  K = length(study_names)

  pi_vec = get_estimated_pi(global_mash)
  Ulist = global_mash$fitted_g$Ulist

  if (is.null(Ulist)) {
    stop("No fitted Ulist found in global_mash for prefix ", prefix, ".")
  }

  if (is.null(names(Ulist))) {
    stop("Ulist does not have names for prefix ", prefix, ", so pi weights cannot be mapped back to covariance patterns.")
  }

  message("Constructing pi-weighted covariance matrix for ", prefix, "...")

  base_names = names(Ulist)

  pi_by_base = sapply(base_names, function(u_name) {
    sum(pi_vec[startsWith(names(pi_vec), u_name)])
  })

  if (all(pi_by_base == 0)) {
    stop("Failed to map estimated mixture weights back to base covariance names for prefix ", prefix, ".")
  }

  Sigma = matrix(0, nrow = K, ncol = K)

  for (u_name in base_names) {
    Ui = as.matrix(Ulist[[u_name]])

    if (nrow(Ui) != K || ncol(Ui) != K) {
      stop(
        "Component ", u_name, " has dimensions ", nrow(Ui), " x ", ncol(Ui),
        " but expected ", K, " x ", K, " for prefix ", prefix, "."
      )
    }

    Sigma = Sigma + as.numeric(pi_by_base[[u_name]]) * Ui
  }

  Sigma = symmetrize_matrix(Sigma)
  rownames(Sigma) = study_names
  colnames(Sigma) = study_names

  cor_mat = safe_cor_from_cov(Sigma)

  eig = eigen(Sigma, symmetric = TRUE)
  rownames(eig$vectors) = study_names

  eig_values = data.frame(
    component = seq_along(eig$values),
    eigenvalue = eig$values,
    proportion_variance_explained = eig$values / sum(eig$values),
    cumulative_pve = cumsum(eig$values / sum(eig$values))
  )

  summary_df = data.frame(
    prefix = prefix,
    n_studies = K,
    n_mixture_components = length(pi_vec),
    n_base_covariance_patterns = length(Ulist),
    stringsAsFactors = FALSE
  )

  mixture_weights = data.frame(
    component = names(pi_vec),
    pi = as.numeric(pi_vec),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(pi)) %>%
    mutate(rank = seq_len(n()))

  base_covariance_weights = data.frame(
    covariance_pattern = names(pi_by_base),
    pi_collapsed = as.numeric(pi_by_base),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(pi_collapsed)) %>%
    mutate(rank = seq_len(n()))

  eig_df = data.frame(
    study = rownames(eig$vectors),
    eig$vectors,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (!("analyte" %in% colnames(metadata_qtl))) {
    stop("metadata_qtl does not contain an 'analyte' column for prefix ", prefix, ".")
  }

  analyte_contrast = ifelse(metadata_qtl$analyte == "protein", 1, 0)

  pc_analyte_cor =
    lapply(seq_len(ncol(eig$vectors)), function(i) {
      x = eig$vectors[, i]
      keep = is.finite(x) & is.finite(analyte_contrast)

      if (sum(keep) < 3L || length(unique(analyte_contrast[keep])) < 2L) {
        return(data.frame(
          component = i,
          pc_name = paste0("PC", i),
          correlation_with_protein_vs_mRNA = NA_real_,
          abs_correlation_with_protein_vs_mRNA = NA_real_,
          p_value = NA_real_,
          proportion_variance_explained = eig_values$proportion_variance_explained[i],
          cumulative_pve = eig_values$cumulative_pve[i],
          stringsAsFactors = FALSE
        ))
      }

      ct = suppressWarnings(cor.test(x[keep], analyte_contrast[keep], method = "pearson"))

      data.frame(
        component = i,
        pc_name = paste0("PC", i),
        correlation_with_protein_vs_mRNA = unname(ct$estimate),
        abs_correlation_with_protein_vs_mRNA = abs(unname(ct$estimate)),
        p_value = ct$p.value,
        proportion_variance_explained = eig_values$proportion_variance_explained[i],
        cumulative_pve = eig_values$cumulative_pve[i],
        stringsAsFactors = FALSE
      )
    }) %>%
    bind_rows() %>%
    mutate(
      p_value_bonf = p.adjust(p_value, method = "bonferroni")
    ) %>%
    arrange(desc(abs_correlation_with_protein_vs_mRNA))

  hc = hclust(as.dist(1 - cor_mat), method = "average")

  pc_plot_dir_unclustered = file.path(odir, "pc_barplots_unclustered")
  pc_plot_dir_clustered = file.path(odir, "pc_barplots_clustered")

  save_all_pc_barplots(
    eig = eig,
    outdir = pc_plot_dir_unclustered,
    study_order = seq_len(nrow(cor_mat))
  )

  save_all_pc_barplots(
    eig = eig,
    outdir = pc_plot_dir_clustered,
    study_order = hc$order
  )

  clustered_study_order =
    data.frame(
      study = study_names,
      analyte = metadata_qtl$analyte,
      stringsAsFactors = FALSE
    ) %>%
    mutate(cluster_rank = match(study, study_names[hc$order])) %>%
    arrange(cluster_rank)

  write_tsv_if_needed(summary_df, outfile("pi_weighted_summary.tsv"))
  write_tsv_if_needed(eig_df, outfile("pi_weighted_eigenvectors.tsv"))
  write_tsv_if_needed(eig_values, outfile("pi_weighted_eigenvalues.tsv"))
  write_tsv_if_needed(as.data.frame(cor_mat), outfile("pi_weighted_correlation_matrix.tsv"))
  write_tsv_if_needed(as.data.frame(Sigma), outfile("pi_weighted_covariance_matrix.tsv"))
  write_tsv_if_needed(mixture_weights, outfile("mixture_weights_ranked.tsv"))
  write_tsv_if_needed(base_covariance_weights, outfile("base_covariance_weights.tsv"))
  write_tsv_if_needed(pc_analyte_cor, outfile("pc_vs_protein_mrna_contrast.tsv"))
  write_tsv_if_needed(clustered_study_order, outfile("clustered_study_order.tsv"))

  saveRDS(Sigma, outfile("pi_weighted_covariance_matrix.rds"))
  saveRDS(cor_mat, outfile("pi_weighted_correlation_matrix.rds"))
  saveRDS(hc, outfile("pi_weighted_correlation_hclust.rds"))

  figure_title = "\u03c0-weighted mash correlation"

  plot_figure(
    cor_mat = cor_mat,
    eig = eig,
    main_title = figure_title,
    pdf_file = outfile("figure3_style_pi_unclustered.pdf"),
    png_file = outfile("figure3_style_pi_unclustered.png"),
    clustering = NULL
  )

  plot_figure(
    cor_mat = cor_mat,
    eig = eig,
    main_title = figure_title,
    pdf_file = outfile("figure3_style_pi_clustered.pdf"),
    png_file = outfile("figure3_style_pi_clustered.png"),
    clustering = hc
  )

  message("Done with ", prefix, ".")
  message("Constructed pi-weighted covariance from ", length(pi_vec), " mixture components.")
  message("Saved unclustered and clustered figures plus PC/analyte correlation summaries to: ", odir)
}

invisible(lapply(prefix, run_one_prefix))
message("Script Finished!")