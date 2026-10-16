########################################
## FRANK REVIEW ACTUAL LOD POWER       ##
########################################

# Standalone: never source remap_functions.R (it installs packages).
# Usage: Rscript --vanilla iReviewFrankPower.R P05_ROOT NEW_OUTPUT_DIRECTORY [THRESHOLD_AUDIT_DIRECTORY]
args = commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(2L, 3L))) stop("Supply root, new output and optional audit directory.")
root = normalizePath(args[[1]], mustWork = TRUE)
out = args[[2]]
if (file.exists(out)) stop("Output already exists; refusing to overwrite.")
for (pkg in c("GenomicRanges", "IRanges")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg)
}
read_tsv = function(path) {
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}
load_object = function(path, name) {
  env = new.env(parent = emptyenv())
  found = load(path, envir = env)
  if (!(name %in% found)) stop("Missing object ", name, " in ", path)
  env[[name]]
}
write_tsv = function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}
interval_padding_bp = 10000L

normalize_yeast_chr = function(x) {

  x = toupper(trimws(as.character(x)))
  x = sub("^CHR", "", x)

  roman = c(
    "I", "II", "III", "IV",
    "V", "VI", "VII", "VIII",
    "IX", "X", "XI", "XII",
    "XIII", "XIV", "XV", "XVI"
  )

  idx = match(x, roman)
  use = !is.na(idx)
  x[use] = as.character(idx[use])
  x
}

make_frank_qtl_granges = function(qtls, pad = interval_padding_bp) {

  required = c(
    "gene_ID",
    "CI_chr",
    "CI_min_base",
    "CI_max_base",
    "cor_coef"
  )
  missing = setdiff(required, colnames(qtls))
  if (length(missing) > 0L) {
    stop("Frank QTL table missing: ", paste(missing, collapse = ", "))
  }

  left = suppressWarnings(as.numeric(qtls$CI_min_base))
  right = suppressWarnings(as.numeric(qtls$CI_max_base))
  if (any(!is.finite(left) | !is.finite(right))) {
    stop("Non-finite Frank QTL confidence intervals encountered.")
  }

  seq_key = paste(
    normalize_yeast_chr(qtls$CI_chr),
    qtls$gene_ID,
    sep = "_"
  )
  strand = ifelse(
    qtls$cor_coef > 0,
    "+",
    ifelse(qtls$cor_coef < 0, "-", "*")
  )

  GenomicRanges::GRanges(
    seqnames = seq_key,
    ranges = IRanges::IRanges(
      start = pmax(1, pmin(left, right) - as.integer(pad)),
      end = pmax(left, right) + as.integer(pad)
    ),
    strand = strand
  )
}

make_technical_marker_ranges = function(chr_features) {

  technical_orf = c(
    LEU2 = "YCL018W",
    MAT = "YCR039C",
    URA3 = "YEL021W",
    AMN1 = "YBR158W",
    KRE33 = "YNL132W"
  )

  idx = match(unname(technical_orf), chr_features$Feature_Name)
  if (any(is.na(idx))) {
    stop(
      "Technical-marker ORFs missing from chrFeatures: ",
      paste(names(technical_orf)[is.na(idx)], collapse = ", ")
    )
  }

  left = as.numeric(chr_features$Start_Coordinate[idx])
  right = as.numeric(chr_features$Stop_Coordinate[idx])

  gr = GenomicRanges::GRanges(
    seqnames = normalize_yeast_chr(chr_features$Chromosome[idx]),
    ranges = IRanges::IRanges(
      start = pmin(left, right),
      end = pmax(left, right)
    )
  )
  names(gr) = names(technical_orf)
  gr
}

frank_technical_markers = function() {

  standard = c("LEU2", "URA3", "AMN1", "MAT")

  list(
    Brem_2002 = standard,
    Brem_2005 = standard,
    SmithEth_2008 = standard,
    SmithGlc_2008 = standard,
    Khan_2009 = standard,
    Foss_2011 = standard,
    Picotti_2013 = standard,
    GrossbachRNA_2022 = standard,
    GrossbachPRT_2022 = standard,
    Albert_2018 = c("MAT", "KRE33")
  )
}

prepare_figure3_query_qtls = function(
  query,
  target,
  qtl_by_study,
  studies,
  technical_marker_ranges,
  pad = interval_padding_bp
) {

  target_genes = colnames(studies[[target]]$pheno)
  query_qtls = qtl_by_study[[query]]
  query_qtls = query_qtls[
    query_qtls$gene_ID %in% target_genes,
    ,
    drop = FALSE
  ]

  tech_list = frank_technical_markers()
  exclude_markers = setdiff(
    tech_list[[query]],
    tech_list[[target]]
  )

  if (nrow(query_qtls) > 0L && length(exclude_markers) > 0L) {

    left = as.numeric(query_qtls$CI_min_base)
    right = as.numeric(query_qtls$CI_max_base)

    query_chr_gr = GenomicRanges::GRanges(
      seqnames = normalize_yeast_chr(query_qtls$CI_chr),
      ranges = IRanges::IRanges(
        start = pmax(1, pmin(left, right) - as.integer(pad)),
        end = pmax(left, right) + as.integer(pad)
      )
    )

    technical_overlap = suppressWarnings(
      GenomicRanges::countOverlaps(
        query_chr_gr,
        technical_marker_ranges[exclude_markers],
        ignore.strand = TRUE
      )
    ) > 0L

    query_qtls = query_qtls[
      !technical_overlap,
      ,
      drop = FALSE
    ]
  }

  query_qtls
}

same_direction_detection = function(abs_r, n, lod) {

  if (any(!is.finite(abs_r)) || any(abs_r <= 0 | abs_r >= 1)) {
    stop("Effect magnitudes must satisfy 0 < abs_r < 1.")
  }
  if (length(n) != 1L || !is.finite(n) || n <= 2) {
    stop("n must be one finite scalar greater than 2.")
  }

  if (!(length(lod) %in% c(1L, length(abs_r))) ||
      any(!is.finite(lod)) || any(lod <= 0)) stop("Invalid LOD thresholds.")
  rcrit = sqrt(pmax(0, 1 - 10 ^ (-2 * lod / n)))
  df = n - 2
  tcrit = rcrit * sqrt(df / (1 - rcrit ^ 2))
  ncp = abs_r * sqrt((n - 1) / (1 - abs_r ^ 2))

  out = stats::pt(
    tcrit,
    df = df,
    ncp = ncp,
    lower.tail = FALSE
  )
  if (any(!is.finite(out))) stop("Non-finite detection probability.")
  pmax(0, pmin(1, out))
}


metadata = read.csv(file.path(root, "raw/studyMetadata.csv"),
                     stringsAsFactors = FALSE)
metadata = metadata[metadata$mapping_approach == "QTL" &
                    metadata$study != "AlbertAlt_2018", ]
if (nrow(metadata) != 10L || anyDuplicated(metadata$study)) {
  stop("Expected 10 unique individual datasets.")
}
studies = load_object(file.path(root, "imed/studies.rda"), "studies")
frank = load_object(file.path(root, "imed/Frank/qtls_list_260204.rda"), "qtls_list")
name_map = c(
  Brem_2002 = "Brem_2002", Brem_2005 = "Brem_2005",
  SmithEth_2008 = "Smith_2008_Ethanol", SmithGlc_2008 = "Smith_2008_Glucose",
  Khan_2009 = "Khan_2009", Foss_2011 = "Foss_2011",
  Picotti_2013 = "Picotti_2013", Albert_2018 = "Albert_2018",
  GrossbachRNA_2022 = "Grossbach_2022_mRNA",
  GrossbachPRT_2022 = "Grossbach_2022_protein"
)
stopifnot(setequal(metadata$study, names(name_map)),
          all(metadata$study %in% names(studies)),
          all(name_map %in% names(frank)))
qtl_by_study = lapply(name_map, function(x) frank[[x]])
technical = make_technical_marker_ranges(
  load_object(file.path(root, "raw/chrFeatures.rda"), "chrFeatures")
)

# Reuse the existing audit, but verify its exact log records and completed calls.
audit_dir = if (length(args) == 3L) args[[3]] else file.path(root, "asys/review/power")
thresholds = read_tsv(file.path(audit_dir, "threshold_full_first_pass_by_study_chr.tsv"))
calls = read_tsv(file.path(audit_dir, "threshold_selected_full_calls.tsv"))
required = c("study", "chromosome_number", "lod", "line_number", "call", "available")
stopifnot(all(required %in% names(thresholds)))
key = paste(thresholds$study, thresholds$chromosome_number, sep = "::")
expected = expand.grid(study = metadata$study, chromosome_number = 1:16)
stopifnot(nrow(thresholds) == 160L, !anyDuplicated(key),
          setequal(key, paste(expected$study, expected$chromosome_number, sep = "::")),
          all(thresholds$available), all(is.finite(thresholds$lod)),
          all(thresholds$lod > 0),
          nrow(calls) == 10L, !anyDuplicated(calls$study),
          setequal(calls$study, metadata$study),
          all(calls$mode == "f"), all(is.finite(calls$finish_line)))
log_lines = readLines(file.path(root, "logs/iRemapAll.out"), warn = FALSE)
stopifnot(all(thresholds$line_number >= 1),
          all(thresholds$line_number <= length(log_lines)))
audit_lines = log_lines[thresholds$line_number]
stopifnot(all(grepl("FDR Thresh 5%", audit_lines, fixed = TRUE)))
logged_lod = suppressWarnings(as.numeric(sub(
  ".*FDR Thresh 5%[[:space:]]+([-+0-9.eE]+).*$", "\\1", audit_lines
)))
stopifnot(all(is.finite(logged_lod)),
          max(abs(logged_lod - thresholds$lod)) < 1e-10)
ci = match(thresholds$study, calls$study)
stopifnot(all(thresholds$call == calls$call[ci]),
          all(thresholds$line_number > calls$start_line[ci]),
          all(thresholds$line_number < calls$finish_line[ci]))

# Fixed-effect power engine checks. These verify identities, not real-data output.
stopifnot(abs(same_direction_detection(sqrt(0.025), 1012, 3.03) -
              0.91038) < 0.001)
stopifnot(same_direction_detection(0.2, 100, 7) <
          same_direction_detection(0.2, 100, 2))

pairs = expand.grid(query = metadata$study, target = metadata$study,
                    stringsAsFactors = FALSE)
pairs = pairs[pairs$query != pairs$target, ]
results = list()
loci = list()
for (i in seq_len(nrow(pairs))) {
  query = pairs$query[[i]]
  target = pairs$target[[i]]
  n = nrow(studies[[target]]$pheno)
  historical = qtl_by_study[[query]]
  historical = historical[historical$gene_ID %in%
                            colnames(studies[[target]]$pheno), ]
  matched = prepare_figure3_query_qtls(
    query, target, qtl_by_study, studies, technical
  )
  observed = if (nrow(matched)) {
    mean(suppressWarnings(GenomicRanges::countOverlaps(
      make_frank_qtl_granges(matched),
      make_frank_qtl_granges(qtl_by_study[[target]]),
      ignore.strand = FALSE
    )) > 0L)
  } else NA_real_
  for (denominator in c("historical_gene_overlap", "figure3_matched")) {
    q = if (denominator == "figure3_matched") matched else historical
    if (!nrow(q)) stop("Empty source QTL set: ", query, " -> ", target)
    chr = normalize_yeast_chr(q$CI_chr)
    ix = match(paste(target, chr, sep = "::"), key)
    if (anyNA(ix)) stop("Missing target chromosome threshold.")
    lod = thresholds$lod[ix]
    effect = abs(q$cor_coef)
    for (scale in c(1, 0.8)) {
      p_actual = same_direction_detection(effect * scale, n, lod)
      p_fixed = same_direction_detection(effect * scale, n, 4.5)
      k = length(results) + 1L
      results[[k]] = data.frame(
        query = query, target = target,
        fromData = name_map[[query]], toData = name_map[[target]],
        denominator = denominator, effect_scale = scale,
        target_n = n, n_power_qtls = nrow(q), n_observed_qtls = nrow(matched),
        observed = observed, expected_actual_first_pass = mean(p_actual),
        expected_fixed_lod4_5 = mean(p_fixed),
        min_target_lod = min(lod), max_target_lod = max(lod)
      )
      loci[[k]] = data.frame(query = query, target = target,
        denominator = denominator, effect_scale = scale,
        source_row = seq_len(nrow(q)), gene = q$gene_ID,
        chromosome = chr, abs_source_r = effect, target_n = n,
        target_lod = lod, threshold_log_line = thresholds$line_number[ix],
        expected_actual_first_pass = p_actual, expected_fixed_lod4_5 = p_fixed)
    }
  }
}
result = do.call(rbind, results)
stopifnot(nrow(result) == 360L, all(is.finite(result$observed)))
summaries = list()
for (denom in unique(result$denominator)) {
  for (scale in c(1, 0.8)) {
    for (subset in c("all_90_pairs", "albert_9_targets")) {
      d = result[result$denominator == denom & result$effect_scale == scale, ]
      if (subset == "albert_9_targets") d = d[d$query == "Albert_2018", ]
      for (method in c("expected_actual_first_pass", "expected_fixed_lod4_5")) {
        pred = d[[method]]
        delta = pred - d$observed
        summaries[[length(summaries) + 1L]] = data.frame(
          denominator = denom, effect_scale = scale, subset = subset,
          method = method, n_pairs = nrow(d),
          mean_observed = mean(d$observed), mean_expected = mean(pred),
          mean_error = mean(delta), mean_absolute_error = mean(abs(delta)),
          rmse = sqrt(mean(delta^2)), pearson_r = cor(pred, d$observed),
          spearman_rho = cor(pred, d$observed, method = "spearman")
        )
      }
    }
  }
}
# Outputs are created only after input checks and calculations succeed.
if (!dir.create(out, recursive = TRUE)) stop("Could not create new output directory.")
write_tsv(result, "pairwise_actual_threshold_power.tsv")
write_tsv(do.call(rbind, summaries), "actual_threshold_fit_summary.tsv")
write_tsv(do.call(rbind, loci), "actual_threshold_locus_predictions.tsv")
write_tsv(thresholds, "thresholds_used.tsv")
writeLines(c(
  "Uses recovered FULL-map first-pass thresholds per target study and chromosome.",
  "This is an analytic fixed-effect detection benchmark, not a replay of iterative mapping.",
  "Later conditional passes, localization, environment and measurement changes are not simulated.",
  "Verify that the selected mapping log belongs to the QTL tables used in Figure 3 before interpretation.",
  "Both historical and Figure 3 matched denominators are shown so threshold and denominator changes can be separated.",
  "Effect scale 0.8 is a prespecified sensitivity, not a fitted winner's-curse correction."
), file.path(out, "interpretation.txt"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
print(do.call(rbind, summaries), row.names = FALSE)

