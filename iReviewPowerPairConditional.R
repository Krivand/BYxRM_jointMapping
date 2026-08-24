############################################
############################################
## GENETICS REVIEW POWER ANALYSIS         ##
## PAIR-CONDITIONAL REDISCOVERY           ##
############################################
############################################

#################
# PRELIMINARIES #
#################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
review_imed = paste0(imed, "review/")
power_asys = paste0(review_asys, "power_pair_conditional/")
power_imed = paste0(review_imed, "power_pair_conditional/")

ensure_directory(review_asys)
ensure_directory(review_imed)
ensure_directory(power_asys)
ensure_directory(power_imed)

validation_seed = suppressWarnings(
  as.integer(Sys.getenv("REVIEW_POWER_SEED", "20260818"))
)
if (is.na(validation_seed)) validation_seed = 20260818L

n_validation = suppressWarnings(
  as.integer(Sys.getenv("REVIEW_POWER_NVALID", "5000"))
)
if (is.na(n_validation) || n_validation < 1000L) {
  stop("REVIEW_POWER_NVALID must be at least 1000.")
}

lod_sweep = seq(2, 7, by = 0.25)
lod_display = 4.5
interval_padding_bp = 10000L

message(
  "Pair-conditional reviewer power analysis: LOD ",
  min(lod_sweep), " to ", max(lod_sweep),
  " by ", unique(diff(lod_sweep))
)
message(
  "Direct-validation replicates per setting: ",
  format(n_validation, big.mark = ",")
)

####################
# HELPER FUNCTIONS #
####################

load_single_rda = function(path, object_name = NULL) {

  if (!file.exists(path)) {
    stop("Missing required file: ", path)
  }

  env = new.env(parent = emptyenv())
  loaded = load(path, envir = env)

  if (!is.null(object_name)) {
    if (!(object_name %in% loaded)) {
      stop("Expected object '", object_name, "' in ", path)
    }
    return(env[[object_name]])
  }

  if (length(loaded) != 1L) {
    stop(
      "Expected one object in ", path,
      "; found: ", paste(loaded, collapse = ", ")
    )
  }

  env[[loaded[[1]]]]
}

write_table = function(x, filename) {
  write.table(
    x,
    file = paste0(power_asys, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

write_table_gz = function(x, filename) {
  con = gzfile(paste0(power_asys, filename), open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(
    x,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

lod_to_rcrit = function(lod, n) {

  if (any(n <= 2)) {
    stop("Sample sizes must exceed 2.")
  }

  sqrt(
    pmax(
      0,
      1 - 10 ^ (-2 * lod / n)
    )
  )
}

correlation_call_prob = function(rho, n, lod) {

  if (length(n) != 1L || n <= 2) {
    stop("n must be a scalar greater than 2.")
  }
  if (length(lod) != 1L || !is.finite(lod)) {
    stop("lod must be one finite scalar.")
  }
  if (any(!is.finite(rho))) {
    stop("rho contains non-finite values.")
  }
  if (any(abs(rho) >= 1)) {
    stop("rho must satisfy abs(rho) < 1.")
  }

  rcrit = lod_to_rcrit(lod = lod, n = n)
  df = n - 2
  tcrit = rcrit * sqrt(df / (1 - rcrit ^ 2))
  ncp = rho * sqrt((n - 1) / (1 - rho ^ 2))

  p_positive = stats::pt(
    tcrit,
    df = df,
    ncp = ncp,
    lower.tail = FALSE
  )
  p_negative = stats::pt(
    -tcrit,
    df = df,
    ncp = ncp,
    lower.tail = TRUE
  )

  p_positive[!is.finite(p_positive)] = 0
  p_negative[!is.finite(p_negative)] = 0
  p_positive = pmax(0, pmin(1, p_positive))
  p_negative = pmax(0, pmin(1, p_negative))

  data.frame(
    rho = rho,
    p_positive = p_positive,
    p_negative = p_negative,
    power = p_positive + p_negative,
    stringsAsFactors = FALSE
  )
}

same_direction_detection_prob = function(rho, n, lod) {

  probs = correlation_call_prob(
    rho = rho,
    n = n,
    lod = lod
  )

  same = ifelse(
    rho > 0,
    probs$p_positive,
    ifelse(rho < 0, probs$p_negative, 0)
  )
  wrong = ifelse(
    rho > 0,
    probs$p_negative,
    ifelse(rho < 0, probs$p_positive, 0)
  )

  data.frame(
    rho = rho,
    p_same_direction = same,
    p_wrong_direction = wrong,
    p_any_detection = probs$power,
    stringsAsFactors = FALSE
  )
}

validate_direct_simulation = function(
  n,
  rho,
  lod,
  n_rep,
  seed
) {

  set.seed(seed)

  x = rep(c(-1, 1), length.out = n)
  x = as.numeric(scale(x))
  residual_sd = sqrt(1 - rho ^ 2)

  noise = matrix(
    stats::rnorm(n * n_rep),
    nrow = n,
    ncol = n_rep
  )
  y = rho * x + residual_sd * noise
  y = sweep(y, 2, colMeans(y), "-")

  y_sd = sqrt(colSums(y ^ 2) / (n - 1))
  r_hat = drop(crossprod(x, y)) / ((n - 1) * y_sd)
  r_hat = pmax(-0.999999, pmin(0.999999, r_hat))
  lod_hat = -n * log(1 - r_hat ^ 2) / (2 * log(10))

  analytic = same_direction_detection_prob(
    rho = rho,
    n = n,
    lod = lod
  )

  data.frame(
    n = n,
    rho = rho,
    lod = lod,
    n_rep = n_rep,
    analytic_power = analytic$p_any_detection,
    simulated_power = mean(lod_hat > lod),
    analytic_same_direction = analytic$p_same_direction,
    simulated_same_direction = mean(
      lod_hat > lod & sign(r_hat) == sign(rho)
    ),
    stringsAsFactors = FALSE
  )
}

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

prepare_pair_query_qtls = function(
  query,
  target,
  qtl_by_study,
  studies,
  technical_marker_ranges,
  pad = interval_padding_bp
) {

  query_genes = colnames(studies[[query]]$pheno)
  target_genes = colnames(studies[[target]]$pheno)

  query_qtls = qtl_by_study[[query]]
  query_qtls$query_qtl_id = seq_len(nrow(query_qtls))
  query_qtls = query_qtls[
    query_qtls$gene_ID %in% target_genes,
    ,
    drop = FALSE
  ]
  n_before_technical = nrow(query_qtls)

  tech_list = frank_technical_markers()
  exclude_markers = setdiff(
    tech_list[[query]],
    tech_list[[target]]
  )
  n_technical_excluded = 0L

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

    n_technical_excluded = sum(technical_overlap)
    query_qtls = query_qtls[
      !technical_overlap,
      ,
      drop = FALSE
    ]
  }

  list(
    qtls = query_qtls,
    n_common_genes = length(intersect(query_genes, target_genes)),
    n_before_technical = n_before_technical,
    n_technical_excluded = n_technical_excluded,
    query_genes = query_genes,
    target_genes = target_genes
  )
}

########################
# LOAD INDIVIDUAL DATA #
########################

studies = load_single_rda(
  paste0(imed, "studies.rda"),
  object_name = "studies"
)
indivFullQTLs = load_single_rda(
  paste0(imed, "indivFullQTLs.rda"),
  object_name = "indivFullQTLs"
)

analysis_studies = setdiff(
  classicMetadata$study,
  "AlbertAlt_2018"
)

missing_studies = setdiff(analysis_studies, names(studies))
if (length(missing_studies) > 0L) {
  stop(
    "Missing processed individual studies: ",
    paste(missing_studies, collapse = ", ")
  )
}

study_design = classicMetadata[
  match(analysis_studies, classicMetadata$study),
  c("study", "analyte", "medium"),
  drop = FALSE
]
study_design$n = vapply(
  study_design$study,
  function(study) nrow(studies[[study]]$pheno),
  integer(1)
)
study_design$n_genes = vapply(
  study_design$study,
  function(study) ncol(studies[[study]]$pheno),
  integer(1)
)
qtl_counts = table(indivFullQTLs$study)
study_design$n_qtl = as.integer(qtl_counts[study_design$study])
study_design$n_qtl[is.na(study_design$n_qtl)] = 0L

frank_qtls_list = load_single_rda(
  paste0(imed, "Frank/qtls_list_260204.rda"),
  object_name = "qtls_list"
)
frank_name_map = c(
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

missing_frank = setdiff(unname(frank_name_map), names(frank_qtls_list))
if (length(missing_frank) > 0L) {
  stop(
    "Missing individual datasets from Frank QTL object: ",
    paste(missing_frank, collapse = ", ")
  )
}

frank_individual_qtls = lapply(
  frank_name_map,
  function(name) frank_qtls_list[[name]]
)
names(frank_individual_qtls) = names(frank_name_map)

frank_qtl_count_check = data.frame(
  study = study_design$study,
  indivFullQTLs = study_design$n_qtl,
  frank_qtls_list = vapply(
    study_design$study,
    function(study) nrow(frank_individual_qtls[[study]]),
    integer(1)
  ),
  stringsAsFactors = FALSE
)
frank_qtl_count_check$match =
  frank_qtl_count_check$indivFullQTLs == frank_qtl_count_check$frank_qtls_list
if (!all(frank_qtl_count_check$match)) {
  stop(
    "Frank/current QTL-count mismatch for: ",
    paste(
      frank_qtl_count_check$study[!frank_qtl_count_check$match],
      collapse = ", "
    )
  )
}

chrFeatures = load_single_rda(
  paste0(raw, "chrFeatures.rda"),
  object_name = "chrFeatures"
)
technical_marker_ranges = make_technical_marker_ranges(chrFeatures)

#############################################
# FIX THE ACTUAL FIGURE 3 QUERY DENOMINATOR #
#############################################

pair_grid = expand.grid(
  query = study_design$study,
  target = study_design$study,
  stringsAsFactors = FALSE
)
pair_grid = pair_grid[pair_grid$query != pair_grid$target, , drop = FALSE]

pair_results = lapply(
  seq_len(nrow(pair_grid)),
  function(i) {

    query = pair_grid$query[[i]]
    target = pair_grid$target[[i]]

    prepared = prepare_pair_query_qtls(
      query = query,
      target = target,
      qtl_by_study = frank_individual_qtls,
      studies = studies,
      technical_marker_ranges = technical_marker_ranges,
      pad = interval_padding_bp
    )

    query_qtls = prepared$qtls
    target_qtls = frank_individual_qtls[[target]]
    n_query = nrow(query_qtls)

    if (n_query == 0L || nrow(target_qtls) == 0L) {
      n_match = 0L
    } else {
      query_gr = make_frank_qtl_granges(query_qtls)
      target_gr = make_frank_qtl_granges(target_qtls)
      overlap_count = suppressWarnings(
        GenomicRanges::countOverlaps(
          query_gr,
          target_gr,
          ignore.strand = FALSE
        )
      )
      n_match = sum(overlap_count > 0L)
    }

    eligible = query_qtls
    if (nrow(eligible) > 0L) {
      eligible$query = query
      eligible$target = target
      eligible$query_n = study_design$n[
        match(query, study_design$study)
      ]
      eligible$target_n = study_design$n[
        match(target, study_design$study)
      ]
    }

    summary = data.frame(
      query = query,
      target = target,
      n_common_genes = prepared$n_common_genes,
      n_query_before_technical = prepared$n_before_technical,
      n_technical_excluded = prepared$n_technical_excluded,
      n_query_eligible = n_query,
      n_matches = n_match,
      observed_rediscovery = if (n_query > 0L) n_match / n_query else NA_real_,
      query_n = study_design$n[match(query, study_design$study)],
      target_n = study_design$n[match(target, study_design$study)],
      stringsAsFactors = FALSE
    )

    list(summary = summary, eligible = eligible)
  }
)

observed_pairwise = do.call(
  rbind,
  lapply(pair_results, function(x) x$summary)
)

eligible_nonempty = lapply(pair_results, function(x) x$eligible)
eligible_nonempty = eligible_nonempty[
  vapply(eligible_nonempty, nrow, integer(1)) > 0L
]
pair_eligible_query_qtls = do.call(rbind, eligible_nonempty)
rownames(pair_eligible_query_qtls) = NULL

invalid_r =
  !is.finite(pair_eligible_query_qtls$cor_coef) |
  pair_eligible_query_qtls$cor_coef == 0 |
  abs(pair_eligible_query_qtls$cor_coef) >= 1
if (any(invalid_r)) {
  stop(
    "Eligible query QTLs contain ", sum(invalid_r),
    " invalid correlation coefficients."
  )
}

message(
  "Fixed Figure 3 denominator reconstructed for ",
  nrow(observed_pairwise), " ordered pairs."
)

############################################
# PAIR-CONDITIONAL PERFECT-SHARING POWER  #
############################################

pair_power_at_lod = function(query, target, lod) {

  dat = pair_eligible_query_qtls[
    pair_eligible_query_qtls$query == query &
    pair_eligible_query_qtls$target == target,
    ,
    drop = FALSE
  ]

  if (nrow(dat) == 0L) {
    return(data.frame(
      query = query,
      target = target,
      lod = lod,
      n_query_eligible = 0L,
      mean_abs_query_r = NA_real_,
      median_abs_query_r = NA_real_,
      expected_rediscovery = NA_real_,
      expected_any_detection = NA_real_,
      expected_wrong_direction = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  target_n = study_design$n[match(target, study_design$study)]
  probs = same_direction_detection_prob(
    rho = dat$cor_coef,
    n = target_n,
    lod = lod
  )

  data.frame(
    query = query,
    target = target,
    lod = lod,
    n_query_eligible = nrow(dat),
    mean_abs_query_r = mean(abs(dat$cor_coef)),
    median_abs_query_r = stats::median(abs(dat$cor_coef)),
    expected_rediscovery = mean(probs$p_same_direction),
    expected_any_detection = mean(probs$p_any_detection),
    expected_wrong_direction = mean(probs$p_wrong_direction),
    stringsAsFactors = FALSE
  )
}

power_sweep = do.call(
  rbind,
  lapply(
    lod_sweep,
    function(lod) {
      message("Pair-conditional power: LOD ", lod)
      do.call(
        rbind,
        lapply(
          seq_len(nrow(pair_grid)),
          function(i) {
            pair_power_at_lod(
              query = pair_grid$query[[i]],
              target = pair_grid$target[[i]],
              lod = lod
            )
          }
        )
      )
    }
  )
)

power_sweep = merge(
  power_sweep,
  observed_pairwise[, c(
    "query",
    "target",
    "n_matches",
    "observed_rediscovery",
    "query_n",
    "target_n"
  )],
  by = c("query", "target"),
  all.x = TRUE,
  sort = FALSE
)
power_sweep$difference =
  power_sweep$observed_rediscovery - power_sweep$expected_rediscovery

power_summary = power_sweep %>%
  dplyr::group_by(lod) %>%
  dplyr::summarise(
    n_pairs = sum(is.finite(expected_rediscovery)),
    mean_observed = mean(observed_rediscovery, na.rm = TRUE),
    mean_expected = mean(expected_rediscovery, na.rm = TRUE),
    mean_difference = mean(difference, na.rm = TRUE),
    .groups = "drop"
  )

power_envelope = power_sweep %>%
  dplyr::group_by(query, target) %>%
  dplyr::summarise(
    n_query_eligible = dplyr::first(n_query_eligible),
    observed_rediscovery = dplyr::first(observed_rediscovery),
    expected_minimum = min(expected_rediscovery, na.rm = TRUE),
    expected_maximum = max(expected_rediscovery, na.rm = TRUE),
    lod_at_expected_minimum = lod[which.min(expected_rediscovery)][1],
    lod_at_expected_maximum = lod[which.max(expected_rediscovery)][1],
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    relative_to_expected_range = dplyr::case_when(
      observed_rediscovery < expected_minimum ~ "below",
      observed_rediscovery > expected_maximum ~ "above",
      TRUE ~ "within"
    )
  )

pair_effect_summary = power_sweep %>%
  dplyr::filter(lod == min(lod_sweep)) %>%
  dplyr::select(
    query,
    target,
    n_query_eligible,
    mean_abs_query_r,
    median_abs_query_r,
    query_n,
    target_n
  )

################################
# DIRECT ANALYTIC VALIDATION   #
################################

validation_abs_r = unlist(
  lapply(
    frank_individual_qtls,
    function(x) abs(x$cor_coef)
  ),
  use.names = FALSE
)
validation_abs_r = validation_abs_r[
  is.finite(validation_abs_r) &
  validation_abs_r > 0 &
  validation_abs_r < 1
]

validation_n = unique(c(
  min(study_design$n),
  study_design$n[
    which.min(abs(study_design$n - stats::median(study_design$n)))
  ],
  max(study_design$n)
))
validation_rho = unique(as.numeric(stats::quantile(
  validation_abs_r,
  probs = c(0.10, 0.50, 0.90),
  names = FALSE
)))
validation_grid = expand.grid(
  n = validation_n,
  rho = validation_rho,
  lod = c(min(lod_sweep), lod_display, max(lod_sweep)),
  stringsAsFactors = FALSE
)

validation = do.call(
  rbind,
  lapply(
    seq_len(nrow(validation_grid)),
    function(i) {
      validate_direct_simulation(
        n = validation_grid$n[[i]],
        rho = validation_grid$rho[[i]],
        lod = validation_grid$lod[[i]],
        n_rep = n_validation,
        seed = validation_seed + i
      )
    }
  )
)
validation$power_difference =
  validation$simulated_power - validation$analytic_power
validation$direction_difference =
  validation$simulated_same_direction - validation$analytic_same_direction

#########################
# MODEL DEFINITION TABLE #
#########################

model_definition = data.frame(
  item = c(
    "query_set",
    "effect_input",
    "sharing_null",
    "target_power",
    "query_reascertainment",
    "albert_special_role",
    "lod_strategy",
    "localization",
    "winner_curse"
  ),
  definition = c(
    "actual Figure 3 eligible query QTLs for each ordered pair",
    "observed signed query-study QTL correlation for each eligible QTL",
    "target population correlation equals the observed query correlation",
    "analytic same-direction detection probability at actual target n",
    "none; the observed query denominator is fixed",
    "none; Albert_2018 is treated like every other query study",
    "common target LOD sweep 2-7 by 0.25; no historical threshold averaging",
    "not simulated; detected shared QTLs are assumed positionally matched",
    "observed query effects are used without shrinkage, making power optimistic"
  ),
  stringsAsFactors = FALSE
)

#################
# WRITE OUTPUTS #
#################

write_table(study_design, "study_design.tsv")
write_table(frank_qtl_count_check, "frank_qtl_count_check.tsv")
write_table(observed_pairwise, "observed_individual_rediscovery.tsv")
write_table_gz(pair_eligible_query_qtls, "pair_eligible_query_qtls.tsv.gz")
write_table(pair_effect_summary, "pair_effect_summary.tsv")
write_table(power_sweep, "pair_conditional_power_sweep.tsv")
write_table(power_summary, "pair_conditional_power_summary.tsv")
write_table(power_envelope, "pair_conditional_power_envelope.tsv")
write_table(validation, "analytic_power_validation.tsv")
write_table(model_definition, "model_definition.tsv")

saveRDS(
  list(
    study_design = study_design,
    observed_pairwise = observed_pairwise,
    pair_eligible_query_qtls = pair_eligible_query_qtls,
    power_sweep = power_sweep,
    power_summary = power_summary,
    power_envelope = power_envelope,
    validation = validation,
    model_definition = model_definition
  ),
  file = paste0(power_imed, "pair_conditional_power.rds")
)

################
# SIMPLE PLOTS #
################

mean_observed = mean(observed_pairwise$observed_rediscovery, na.rm = TRUE)

p_summary = ggplot2::ggplot(
  power_summary,
  ggplot2::aes(x = lod, y = mean_expected)
) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::geom_hline(
    yintercept = mean_observed,
    linetype = 2
  ) +
  ggplot2::labs(
    x = "Target LOD threshold",
    y = "Mean rediscovery fraction",
    title = "Pair-conditional expected rediscovery",
    subtitle = "Solid: perfect-sharing expectation; dashed: observed mean"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  paste0(power_asys, "pair_conditional_mean_expected_vs_lod.png"),
  p_summary,
  width = 7,
  height = 5
)
ggplot2::ggsave(
  paste0(power_asys, "pair_conditional_mean_expected_vs_lod.pdf"),
  p_summary,
  width = 7,
  height = 5
)

p_display_data = power_sweep[
  abs(power_sweep$lod - lod_display) < 1e-12,
  ,
  drop = FALSE
]
p_display = ggplot2::ggplot(
  p_display_data,
  ggplot2::aes(
    x = expected_rediscovery,
    y = observed_rediscovery
  )
) +
  ggplot2::geom_point() +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    x = "Expected rediscovery under perfect sharing",
    y = "Observed Figure 3 rediscovery",
    title = paste0("Observed vs expected at display LOD ", lod_display)
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  paste0(power_asys, "pair_conditional_observed_vs_expected_lod4.5.png"),
  p_display,
  width = 6,
  height = 6
)
ggplot2::ggsave(
  paste0(power_asys, "pair_conditional_observed_vs_expected_lod4.5.pdf"),
  p_display,
  width = 6,
  height = 6
)

#################
# FINAL SUMMARY #
#################

message("")
message("Pair-conditional reviewer power analysis complete.")
message("Individual datasets: ", nrow(study_design))
message("Ordered Figure 3 pairs: ", nrow(observed_pairwise))
message("Albert-specific effect reference: none")
message("Query re-ascertainment in power model: none")
message("Query denominator: actual Figure 3 eligible QTL set")
message("LOD sweep: 2 to 7 by 0.25")
message("Outputs: ", power_asys)
