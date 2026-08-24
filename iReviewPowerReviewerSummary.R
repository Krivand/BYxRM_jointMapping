############################################
############################################
## GENETICS REVIEW POWER SUMMARY          ##
## REVIEWER-FACING CALIBRATION TABLES     ##
############################################
############################################

#################
# PRELIMINARIES #
#################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
power_asys = paste0(review_asys, "power_pair_conditional/")
ensure_directory(review_asys)
ensure_directory(power_asys)

interval_padding_bp = 10000L
lod_low = 2
lod_display = 4.5
lod_high = 7

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

write_table = function(x, filename) {
  write.table(
    x,
    file = paste0(power_asys, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
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

make_chr_lengths = function(chr_features) {

  chr = normalize_yeast_chr(chr_features$Chromosome)
  left = suppressWarnings(
    as.numeric(chr_features$Start_Coordinate)
  )
  right = suppressWarnings(
    as.numeric(chr_features$Stop_Coordinate)
  )
  endpoint = pmax(left, right, na.rm = TRUE)

  keep =
    chr %in% as.character(1:16) &
    is.finite(endpoint)

  lengths = tapply(
    endpoint[keep],
    chr[keep],
    max,
    na.rm = TRUE
  )

  lengths[as.character(1:16)]
}

chance_overlap_pair = function(
  query_qtls,
  target_qtls,
  chromosome_lengths,
  pad = interval_padding_bp
) {

  required = c(
    "gene_ID",
    "CI_chr",
    "CI_min_base",
    "CI_max_base",
    "cor_coef"
  )

  missing_query = setdiff(required, colnames(query_qtls))
  missing_target = setdiff(required, colnames(target_qtls))
  if (length(missing_query) > 0L || length(missing_target) > 0L) {
    stop(
      "Chance-null QTL tables are missing required columns. Query: ",
      paste(missing_query, collapse = ", "),
      "; target: ", paste(missing_target, collapse = ", ")
    )
  }

  if (nrow(query_qtls) == 0L) {
    return(data.frame(
      n_query_eligible = 0L,
      expected_chance_matches = 0,
      chance_rediscovery = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  q_chr = normalize_yeast_chr(query_qtls$CI_chr)
  t_chr = normalize_yeast_chr(target_qtls$CI_chr)

  q_key = paste(
    q_chr,
    query_qtls$gene_ID,
    sign(query_qtls$cor_coef),
    sep = "::"
  )
  t_key = paste(
    t_chr,
    target_qtls$gene_ID,
    sign(target_qtls$cor_coef),
    sep = "::"
  )
  target_groups = split(
    seq_len(nrow(target_qtls)),
    t_key
  )

  q_left = pmin(
    as.numeric(query_qtls$CI_min_base),
    as.numeric(query_qtls$CI_max_base)
  )
  q_right = pmax(
    as.numeric(query_qtls$CI_min_base),
    as.numeric(query_qtls$CI_max_base)
  )
  t_left = pmin(
    as.numeric(target_qtls$CI_min_base),
    as.numeric(target_qtls$CI_max_base)
  )
  t_right = pmax(
    as.numeric(target_qtls$CI_min_base),
    as.numeric(target_qtls$CI_max_base)
  )

  if (any(!is.finite(q_left) | !is.finite(q_right))) {
    stop("Non-finite query confidence intervals in chance null.")
  }
  if (any(!is.finite(t_left) | !is.finite(t_right))) {
    stop("Non-finite target confidence intervals in chance null.")
  }

  target_width = t_right - t_left

  probability = vapply(
    seq_len(nrow(query_qtls)),
    function(i) {

      idx = target_groups[[q_key[[i]]]]
      if (is.null(idx) || length(idx) == 0L) {
        return(0)
      }

      chr_length = unname(
        chromosome_lengths[q_chr[[i]]]
      )
      widths = target_width[idx]
      n_starts = floor(chr_length - widths + 1)
      valid = is.finite(n_starts) & n_starts >= 1

      if (!any(valid)) {
        return(0)
      }

      widths = widths[valid]
      n_starts = n_starts[valid]

      # The Figure 3 matching rule pads both query and target confidence
      # intervals by 10 kb. The target interval start is randomized uniformly
      # within its chromosome while its gene, chromosome, sign, width and QTL
      # count are held fixed. For multiple target QTLs in the same group, the
      # probability is that at least one independently randomized interval
      # overlaps the fixed query interval.
      lower =
        q_left[[i]] - widths - 2 * as.integer(pad)
      upper =
        q_right[[i]] + 2 * as.integer(pad)

      allowed_low = pmax(1, ceiling(lower))
      allowed_high = pmin(n_starts, floor(upper))
      n_allowed = pmax(
        0,
        allowed_high - allowed_low + 1
      )
      p_each = pmin(1, n_allowed / n_starts)

      1 - prod(1 - p_each)
    },
    numeric(1)
  )

  data.frame(
    n_query_eligible = nrow(query_qtls),
    expected_chance_matches = sum(probability),
    chance_rediscovery = mean(probability),
    stringsAsFactors = FALSE
  )
}

summarize_group = function(
  group,
  pair_rows,
  chance_null,
  sweep
) {

  pair_rows = pair_rows[, c("query", "target"), drop = FALSE]
  pair_rows$key = paste(pair_rows$query, pair_rows$target, sep = "|")

  chance = chance_null[
    paste(chance_null$query, chance_null$target, sep = "|") %in%
      pair_rows$key,
    ,
    drop = FALSE
  ]
  expected = sweep[
    paste(sweep$query, sweep$target, sep = "|") %in%
      pair_rows$key,
    ,
    drop = FALSE
  ]

  at_lod = function(lod) {
    dat = expected[abs(expected$lod - lod) < 1e-12, , drop = FALSE]
    mean(dat$expected_rediscovery, na.rm = TRUE)
  }

  data.frame(
    group = group,
    n_pairs = nrow(chance),
    mean_chance = mean(chance$chance_rediscovery, na.rm = TRUE),
    mean_observed = mean(chance$observed_rediscovery, na.rm = TRUE),
    expected_lod2 = at_lod(lod_low),
    expected_lod4_5 = at_lod(lod_display),
    expected_lod7 = at_lod(lod_high),
    residual_lod4_5 =
      mean(chance$observed_rediscovery, na.rm = TRUE) -
      at_lod(lod_display),
    stringsAsFactors = FALSE
  )
}

################
# READ OUTPUTS #
################

study_design = read_tsv(
  paste0(power_asys, "study_design.tsv")
)
observed = read_tsv(
  paste0(power_asys, "observed_individual_rediscovery.tsv")
)
eligible = read_tsv(
  gzfile(paste0(power_asys, "pair_eligible_query_qtls.tsv.gz"))
)
sweep = read_tsv(
  paste0(power_asys, "pair_conditional_power_sweep.tsv")
)
analytic_validation = read_tsv(
  paste0(power_asys, "analytic_power_validation.tsv")
)

required_lod = c(lod_low, lod_display, lod_high)
missing_lod = setdiff(required_lod, unique(sweep$lod))
if (length(missing_lod) > 0L) {
  stop(
    "Power sweep does not contain required reviewer LOD values: ",
    paste(missing_lod, collapse = ", ")
  )
}

########################################
# RECONSTRUCT ANALYTIC POSITIONAL NULL #
########################################

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
frank_individual_qtls = lapply(
  frank_name_map,
  function(name) frank_qtls_list[[name]]
)
names(frank_individual_qtls) = names(frank_name_map)

chrFeatures = load_single_rda(
  paste0(raw, "chrFeatures.rda"),
  object_name = "chrFeatures"
)
chromosome_lengths = make_chr_lengths(chrFeatures)

message("Calculating analytic Figure 3 positional chance null...")

chance_overlap_null = do.call(
  rbind,
  lapply(
    seq_len(nrow(observed)),
    function(i) {

      query = observed$query[[i]]
      target = observed$target[[i]]
      query_qtls = eligible[
        eligible$query == query &
        eligible$target == target,
        ,
        drop = FALSE
      ]
      target_qtls = frank_individual_qtls[[target]]

      null = chance_overlap_pair(
        query_qtls = query_qtls,
        target_qtls = target_qtls,
        chromosome_lengths = chromosome_lengths,
        pad = interval_padding_bp
      )

      data.frame(
        query = query,
        target = target,
        n_query_eligible = null$n_query_eligible,
        n_matches = observed$n_matches[[i]],
        observed_rediscovery = observed$observed_rediscovery[[i]],
        expected_chance_matches = null$expected_chance_matches,
        chance_rediscovery = null$chance_rediscovery,
        stringsAsFactors = FALSE
      )
    }
  )
)

qmeta = study_design[, c("study", "analyte", "medium"), drop = FALSE]
names(qmeta) = c("query", "query_analyte", "query_medium")
tmeta = study_design[, c("study", "analyte", "medium"), drop = FALSE]
names(tmeta) = c("target", "target_analyte", "target_medium")

chance_overlap_null = merge(
  chance_overlap_null,
  qmeta,
  by = "query",
  all.x = TRUE,
  sort = FALSE
)
chance_overlap_null = merge(
  chance_overlap_null,
  tmeta,
  by = "target",
  all.x = TRUE,
  sort = FALSE
)
chance_overlap_null$analyte_pair = paste0(
  chance_overlap_null$query_analyte,
  "_to_",
  chance_overlap_null$target_analyte
)

########################################
# REVIEWER-FACING PRIMARY SUMMARY      #
########################################

all_pairs = chance_overlap_null[, c("query", "target"), drop = FALSE]
albert_pairs = all_pairs[all_pairs$query == "Albert_2018", , drop = FALSE]
other_pairs = all_pairs[all_pairs$query != "Albert_2018", , drop = FALSE]

reviewer_power_summary = do.call(
  rbind,
  list(
    summarize_group(
      "all_ordered_pairs",
      all_pairs,
      chance_overlap_null,
      sweep
    ),
    summarize_group(
      "Albert_2018_query",
      albert_pairs,
      chance_overlap_null,
      sweep
    ),
    summarize_group(
      "other_studies_query",
      other_pairs,
      chance_overlap_null,
      sweep
    ),
    summarize_group(
      "mRNA_to_mRNA",
      chance_overlap_null[
        chance_overlap_null$analyte_pair == "mRNA_to_mRNA",
        c("query", "target"),
        drop = FALSE
      ],
      chance_overlap_null,
      sweep
    ),
    summarize_group(
      "mRNA_to_protein",
      chance_overlap_null[
        chance_overlap_null$analyte_pair == "mRNA_to_protein",
        c("query", "target"),
        drop = FALSE
      ],
      chance_overlap_null,
      sweep
    ),
    summarize_group(
      "protein_to_mRNA",
      chance_overlap_null[
        chance_overlap_null$analyte_pair == "protein_to_mRNA",
        c("query", "target"),
        drop = FALSE
      ],
      chance_overlap_null,
      sweep
    ),
    summarize_group(
      "protein_to_protein",
      chance_overlap_null[
        chance_overlap_null$analyte_pair == "protein_to_protein",
        c("query", "target"),
        drop = FALSE
      ],
      chance_overlap_null,
      sweep
    )
  )
)
rownames(reviewer_power_summary) = NULL

reviewer_power_table = data.frame(
  summary = c(
    "All ordered comparisons",
    "Albert & Bloom 2018 as query",
    "Other studies as query",
    "mRNA to mRNA",
    "mRNA to protein",
    "protein to mRNA",
    "protein to protein"
  ),
  n_pairs = reviewer_power_summary$n_pairs,
  chance = sprintf(
    "%.2f%%",
    100 * reviewer_power_summary$mean_chance
  ),
  observed = sprintf(
    "%.1f%%",
    100 * reviewer_power_summary$mean_observed
  ),
  expected_lod2 = sprintf(
    "%.1f%%",
    100 * reviewer_power_summary$expected_lod2
  ),
  expected_lod4_5 = sprintf(
    "%.1f%%",
    100 * reviewer_power_summary$expected_lod4_5
  ),
  expected_lod7 = sprintf(
    "%.1f%%",
    100 * reviewer_power_summary$expected_lod7
  ),
  stringsAsFactors = FALSE
)

########################################
# RECIPROCAL DIRECTIONAL ASYMMETRY     #
########################################

reverse = data.frame(
  query = sweep$target,
  target = sweep$query,
  lod = sweep$lod,
  observed_reverse = sweep$observed_rediscovery,
  expected_reverse = sweep$expected_rediscovery,
  stringsAsFactors = FALSE
)

asymmetry = merge(
  sweep[, c(
    "query",
    "target",
    "lod",
    "observed_rediscovery",
    "expected_rediscovery"
  )],
  reverse,
  by = c("query", "target", "lod"),
  all = FALSE,
  sort = FALSE
)
asymmetry = asymmetry[
  asymmetry$query < asymmetry$target,
  ,
  drop = FALSE
]
asymmetry$observed_asymmetry =
  asymmetry$observed_rediscovery -
  asymmetry$observed_reverse
asymmetry$expected_asymmetry =
  asymmetry$expected_rediscovery -
  asymmetry$expected_reverse

asymmetry_summary = do.call(
  rbind,
  lapply(
    sort(unique(asymmetry$lod)),
    function(lod) {
      dat = asymmetry[abs(asymmetry$lod - lod) < 1e-12, , drop = FALSE]
      data.frame(
        lod = lod,
        n_unordered_pairs = nrow(dat),
        spearman = stats::cor(
          dat$observed_asymmetry,
          dat$expected_asymmetry,
          method = "spearman"
        ),
        pearson = stats::cor(
          dat$observed_asymmetry,
          dat$expected_asymmetry,
          method = "pearson"
        ),
        direction_agreement = mean(
          sign(dat$observed_asymmetry) ==
          sign(dat$expected_asymmetry)
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

asymmetry_key_result = data.frame(
  metric = c(
    "spearman_reciprocal_asymmetry",
    "pearson_reciprocal_asymmetry",
    "reciprocal_direction_agreement"
  ),
  minimum = c(
    min(asymmetry_summary$spearman),
    min(asymmetry_summary$pearson),
    min(asymmetry_summary$direction_agreement)
  ),
  maximum = c(
    max(asymmetry_summary$spearman),
    max(asymmetry_summary$pearson),
    max(asymmetry_summary$direction_agreement)
  ),
  lod_min = lod_low,
  lod_max = lod_high,
  stringsAsFactors = FALSE
)

########################################
# MINIMAL SUPPORTING VALIDATION TABLE  #
########################################

pair_validation_path = paste0(
  power_asys,
  "pair_conditional_validation_summary.tsv"
)
pair_validation = if (file.exists(pair_validation_path)) {
  read_tsv(pair_validation_path)
} else {
  NULL
}

max_validation_error = max(
  abs(
    analytic_validation$simulated_power -
    analytic_validation$analytic_power
  ),
  abs(
    analytic_validation$simulated_same_direction -
    analytic_validation$analytic_same_direction
  ),
  na.rm = TRUE
)

reviewer_supporting_validation = data.frame(
  metric = c(
    "pair_conditional_validator",
    "max_analytic_vs_direct_simulation_error",
    "reciprocal_asymmetry_spearman_range",
    "reciprocal_direction_agreement_range"
  ),
  value = c(
    if (is.null(pair_validation)) {
      "not available"
    } else if (all(pair_validation$passed)) {
      "all checks passed"
    } else {
      "one or more checks failed"
    },
    sprintf("%.4f", max_validation_error),
    sprintf(
      "%.3f-%.3f",
      min(asymmetry_summary$spearman),
      max(asymmetry_summary$spearman)
    ),
    sprintf(
      "%.1f%%-%.1f%%",
      100 * min(asymmetry_summary$direction_agreement),
      100 * max(asymmetry_summary$direction_agreement)
    )
  ),
  stringsAsFactors = FALSE
)

#################
# WRITE OUTPUTS #
#################

write_table(chance_overlap_null, "chance_overlap_null.tsv")
write_table(reviewer_power_summary, "reviewer_power_summary.tsv")
write_table(reviewer_power_table, "reviewer_power_table.tsv")
write_table(asymmetry, "reviewer_reciprocal_asymmetry.tsv")
write_table(asymmetry_summary, "reviewer_reciprocal_asymmetry_summary.tsv")
write_table(asymmetry_key_result, "reviewer_reciprocal_asymmetry_key.tsv")
write_table(
  reviewer_supporting_validation,
  "reviewer_supporting_validation.tsv"
)

message("")
message("Reviewer-facing power summary complete.")
message(
  "Primary table: ",
  paste0(power_asys, "reviewer_power_table.tsv")
)
message(
  "Raw summary: ",
  paste0(power_asys, "reviewer_power_summary.tsv")
)
message(
  "Asymmetry summary: ",
  paste0(power_asys, "reviewer_reciprocal_asymmetry_summary.tsv")
)
message(
  "Chance null: ",
  paste0(power_asys, "chance_overlap_null.tsv")
)
