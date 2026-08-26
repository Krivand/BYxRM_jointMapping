############################################
############################################
## GENETICS REVIEW POWER ANALYSIS         ##
## ALL-PAIR EFFECT ATTENUATION SENSITIVITY ##
############################################
############################################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
power_asys = paste0(review_asys, "power_effect_curve/")
allpair_asys = paste0(power_asys, "all_pairs/")

ensure_directory(review_asys)
ensure_directory(power_asys)
ensure_directory(allpair_asys)

lod = 4.5
effect_scales = c(1.0, 0.8, 0.6)
power_cores = suppressWarnings(
  as.integer(Sys.getenv("REVIEW_POWER_CORES", "1"))
)
if (is.na(power_cores) || power_cores < 1L) power_cores = 1L

load_single_rda = function(path, object_name = NULL) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  env = new.env(parent = emptyenv())
  loaded = load(path, envir = env)
  if (!is.null(object_name)) {
    if (!(object_name %in% loaded)) {
      stop("Expected object '", object_name, "' in ", path)
    }
    return(env[[object_name]])
  }
  if (length(loaded) != 1L) {
    stop("Expected one object in ", path, "; found: ", paste(loaded, collapse = ", "))
  }
  env[[loaded[[1]]]]
}

read_tsv = function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_tsv = function(x, filename) {
  write.table(
    x,
    file = paste0(allpair_asys, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

parallel_lapply = function(X, FUN) {
  if (
    .Platform$OS.type == "unix" &&
    power_cores > 1L &&
    length(X) > 1L
  ) {
    return(parallel::mclapply(
      X,
      FUN,
      mc.cores = min(power_cores, length(X)),
      mc.preschedule = TRUE
    ))
  }
  lapply(X, FUN)
}

same_direction_detection = function(abs_r, n, lod) {
  if (any(!is.finite(abs_r)) || any(abs_r <= 0 | abs_r >= 1)) {
    stop("Effect magnitudes must satisfy 0 < abs_r < 1.")
  }
  if (length(n) != 1L || !is.finite(n) || n <= 2) {
    stop("n must be one finite scalar greater than 2.")
  }
  rcrit = sqrt(pmax(0, 1 - 10 ^ (-2 * lod / n)))
  df = n - 2
  tcrit = rcrit * sqrt(df / (1 - rcrit ^ 2))
  ncp = abs_r * sqrt((n - 1) / (1 - abs_r ^ 2))
  out = stats::pt(tcrit, df = df, ncp = ncp, lower.tail = FALSE)
  out[!is.finite(out)] = 0
  pmax(0, pmin(1, out))
}

studies = load_single_rda(
  paste0(imed, "studies.rda"),
  object_name = "studies"
)
qtl_list = load_single_rda(
  paste0(imed, "Frank/qtls_list_260204.rda"),
  object_name = "qtls_list"
)

name_map = c(
  Brem_2002 = "Brem_2002",
  Brem_2005 = "Brem_2005",
  SmithEth_2008 = "Smith_2008_Ethanol",
  SmithGlc_2008 = "Smith_2008_Glucose",
  Khan_2009 = "Khan_2009",
  Foss_2011 = "Foss_2011",
  Picotti_2013 = "Picotti_2013",
  Albert_2018 = "Albert_2018",
  GrossbachRNA_2022 = "Grossbach_2022_mRNA",
  GrossbachPRT_2022 = "Grossbach_2022_protein"
)
reverse_name_map = setNames(names(name_map), unname(name_map))

comparison = read_tsv(
  paste0(allpair_asys, "pairwise_observed_vs_power.tsv")
)

if (
  nrow(comparison) != 90L ||
  anyDuplicated(comparison[, c("fromData", "toData")]) ||
  any(comparison$fromData == comparison$toData)
) {
  stop("Expected 90 unique directed all-pair comparison rows.")
}

missing_names = setdiff(
  unique(c(comparison$fromData, comparison$toData)),
  names(reverse_name_map)
)
if (length(missing_names) > 0L) {
  stop("Unrecognized dataset names: ", paste(missing_names, collapse = ", "))
}

message(
  "Evaluating all 90 directed pairs at effect scales ",
  paste(effect_scales, collapse = ", "),
  " and LOD ", lod, "..."
)

pair_results = parallel_lapply(
  seq_len(nrow(comparison)),
  function(i) {
    from_name = comparison$fromData[[i]]
    to_name = comparison$toData[[i]]
    source = unname(reverse_name_map[[from_name]])
    target = unname(reverse_name_map[[to_name]])

    if (!(source %in% names(studies)) || !(target %in% names(studies))) {
      stop("Missing processed study for ", from_name, " -> ", to_name)
    }
    if (!(from_name %in% names(qtl_list))) {
      stop("Missing QTL table for source dataset ", from_name)
    }

    target_genes = colnames(studies[[target]]$pheno)
    source_qtls = qtl_list[[from_name]]
    source_qtls = source_qtls[
      source_qtls$gene_ID %in% target_genes,
      ,
      drop = FALSE
    ]

    if (nrow(source_qtls) == 0L) {
      stop("No source QTLs remain for ", from_name, " -> ", to_name)
    }

    abs_r = abs(source_qtls$cor_coef)
    invalid = !is.finite(abs_r) | abs_r <= 0 | abs_r >= 1
    if (any(invalid)) {
      stop("Invalid source effects for ", from_name, " -> ", to_name)
    }

    target_n = nrow(studies[[target]]$pheno)
    expected = vapply(
      effect_scales,
      function(effect_scale) {
        mean(
          same_direction_detection(
            abs_r = abs_r * effect_scale,
            n = target_n,
            lod = lod
          )
        )
      },
      numeric(1)
    )

    data.frame(
      fromData = from_name,
      toData = to_name,
      expectedRediscoveryLOD4_5 = expected[[1]],
      expectedRediscoveryLOD4_5_scale0_8 = expected[[2]],
      expectedRediscoveryLOD4_5_scale0_6 = expected[[3]],
      stringsAsFactors = FALSE
    )
  }
)

attenuation = do.call(rbind, pair_results)
rownames(attenuation) = NULL

key = paste(comparison$fromData, comparison$toData, sep = "\r")
attenuation_key = paste(attenuation$fromData, attenuation$toData, sep = "\r")
idx = match(key, attenuation_key)
if (any(is.na(idx))) stop("Could not align attenuation rows to all-pair rows.")
attenuation = attenuation[idx, , drop = FALSE]

base_error = max(abs(
  attenuation$expectedRediscoveryLOD4_5 -
    comparison$expectedRediscoveryLOD4_5
))
if (base_error > 1e-12) {
  stop("Scale 1.0 does not reproduce existing all-pair power: ", base_error)
}

monotone =
  attenuation$expectedRediscoveryLOD4_5 >=
    attenuation$expectedRediscoveryLOD4_5_scale0_8 &
  attenuation$expectedRediscoveryLOD4_5_scale0_8 >=
    attenuation$expectedRediscoveryLOD4_5_scale0_6
if (!all(monotone)) {
  stop("Effect attenuation increased predicted power for at least one pair.")
}

observed = comparison$fractionOverlapsWithStrand
prediction_columns = c(
  "expectedRediscoveryLOD4_5",
  "expectedRediscoveryLOD4_5_scale0_8",
  "expectedRediscoveryLOD4_5_scale0_6"
)

attenuation_summary = do.call(
  rbind,
  lapply(
    seq_along(effect_scales),
    function(i) {
      predicted = attenuation[[prediction_columns[[i]]]]
      difference = observed - predicted
      data.frame(
        effect_scale = effect_scales[[i]],
        n_pairs = length(observed),
        mean_observed = mean(observed),
        mean_expected = mean(predicted),
        mean_observed_minus_expected = mean(difference),
        mean_abs_difference = mean(abs(difference)),
        rmse = sqrt(mean(difference ^ 2)),
        pearson_r = stats::cor(observed, predicted, method = "pearson"),
        spearman_rho = stats::cor(observed, predicted, method = "spearman"),
        stringsAsFactors = FALSE
      )
    }
  )
)
rownames(attenuation_summary) = NULL

source_summary = do.call(
  rbind,
  lapply(
    seq_along(effect_scales),
    function(i) {
      predicted = attenuation[[prediction_columns[[i]]]]
      source_rows = split(seq_len(nrow(comparison)), comparison$fromData)
      do.call(
        rbind,
        lapply(
          names(source_rows),
          function(source_name) {
            rows = source_rows[[source_name]]
            data.frame(
              fromData = source_name,
              effect_scale = effect_scales[[i]],
              n_pairs = length(rows),
              mean_observed = mean(observed[rows]),
              mean_expected = mean(predicted[rows]),
              mean_observed_minus_expected = mean(observed[rows] - predicted[rows]),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)
rownames(source_summary) = NULL

attenuation_output = data.frame(
  fromData = comparison$fromData,
  toData = comparison$toData,
  fractionOverlapsWithStrand = comparison$fractionOverlapsWithStrand,
  attenuation[, prediction_columns, drop = FALSE],
  stringsAsFactors = FALSE
)

write_tsv(attenuation_output, "pairwise_power_attenuation.tsv")
write_tsv(attenuation_summary, "power_attenuation_summary.tsv")
write_tsv(source_summary, "power_attenuation_by_source.tsv")

cat("\nALL-PAIR EFFECT ATTENUATION SENSITIVITY PASSED\n")
cat("Directed pairs: 90\n")
cat("LOD: 4.5\n")
cat("Effect scales: 1.0, 0.8, 0.6\n")
cat("Scale 1.0 maximum reproduction error: ", signif(base_error, 7), "\n", sep = "")
cat("\nOverall attenuation summary:\n")
print(attenuation_summary, row.names = FALSE)
cat("\nBy-source summary:\n")
print(source_summary, row.names = FALSE)
