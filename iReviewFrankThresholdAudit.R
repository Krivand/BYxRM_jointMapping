########################################
########################################
## GENETICS REVIEW POWER THRESHOLDS   ##
## FULL-MAP PERMUTATION AUDIT         ##
########################################
########################################

#################
# PRELIMINARIES #
#################

# Usage: Rscript --vanilla iReviewFrankThresholdAudit.R P05_ROOT NEW_OUTPUT
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply p05 root and new audit output directory.")
root = normalizePath(args[[1]], mustWork = TRUE)
power_asys = paste0(args[[2]], "/")
if (file.exists(args[[2]])) stop("Output exists; refusing overwrite.")
threshold_log = file.path(root, "logs/iRemapAll.out")
if (!file.exists(threshold_log)) stop("Missing remapping log: ", threshold_log)
metadata = read.csv(file.path(root, "raw/studyMetadata.csv"),
                    stringsAsFactors = FALSE)
all_studies = metadata$study[metadata$mapping_approach == "QTL"]
analysis_studies = setdiff(all_studies, "AlbertAlt_2018")
if (length(analysis_studies) != 10L) stop("Expected ten studies.")
# Load only the seed function definition, never the package-installing preamble.
definitions = parse(file.path(root, "remap_functions.R"))
seed_defs = Filter(function(x) is.call(x) &&
  as.character(x[[1]]) %in% c("=", "<-") &&
  is.symbol(x[[2]]) && as.character(x[[2]]) == "seed_from" &&
  is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function")), definitions)
if (length(seed_defs) != 1L) stop("Expected one seed_from definition.")
seed_env = new.env(parent = globalenv())
eval(seed_defs[[1]], envir = seed_env)
study_seed = vapply(all_studies, seed_env$seed_from, integer(1))
names(study_seed) = all_studies
if (!dir.create(args[[2]], recursive = TRUE)) stop("Cannot create audit directory.")

####################
# HELPER FUNCTIONS #
####################

extract_first_match = function(pattern, x) {

  hit = regexec(pattern, x, perl = TRUE)
  matched = regmatches(x, hit)[[1]]

  if (length(matched) < 2L) {
    return(NA_character_)
  }

  matched[[2]]
}

summarize_thresholds = function(x) {

  x = x[is.finite(x)]

  if (length(x) == 0L) {
    return(
      data.frame(
        n_thresholds = 0L,
        minimum = NA_real_,
        q10 = NA_real_,
        q25 = NA_real_,
        median = NA_real_,
        mean = NA_real_,
        q75 = NA_real_,
        q90 = NA_real_,
        maximum = NA_real_,
        sd = NA_real_,
        iqr = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  q = stats::quantile(
    x,
    probs = c(0, 0.10, 0.25, 0.50, 0.75, 0.90, 1),
    names = FALSE,
    na.rm = TRUE
  )

  data.frame(
    n_thresholds = length(x),
    minimum = q[[1]],
    q10 = q[[2]],
    q25 = q[[3]],
    median = q[[4]],
    mean = mean(x),
    q75 = q[[5]],
    q90 = q[[6]],
    maximum = q[[7]],
    sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
    iqr = stats::IQR(x),
    stringsAsFactors = FALSE
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

normalize_chr = function(x) {

  x = trimws(as.character(x))
  x = sub("^chr", "", x, ignore.case = TRUE)
  x = toupper(x)

  roman = c(
    "I", "II", "III", "IV",
    "V", "VI", "VII", "VIII",
    "IX", "X", "XI", "XII",
    "XIII", "XIV", "XV", "XVI"
  )

  idx = match(x, roman)
  out = rep(NA_integer_, length(x))
  out[!is.na(idx)] = idx[!is.na(idx)]

  numeric_x = suppressWarnings(as.integer(x))
  use_numeric = is.na(out) & is.finite(numeric_x)
  out[use_numeric] = numeric_x[use_numeric]
  out
}

########################
# PARSE REMAPPING LOG  #
########################

lines = readLines(threshold_log, warn = FALSE)

current_call = NA_integer_
current_study = NA_character_
current_mode = NA_character_
current_mode_source = NA_character_
current_chr = NA_character_
current_pass = 0L
pending_study = NA_character_
pending_mode = NA_character_
call_counter = 0L

calls = list()
records = list()
record_index = 1L

start_call = function(
  line_number,
  study = NA_character_,
  mode = NA_character_,
  mode_source = NA_character_,
  seed = NA_integer_
) {

  call_counter <<- call_counter + 1L
  current_call <<- call_counter
  current_study <<- study
  current_mode <<- mode
  current_mode_source <<- mode_source
  current_chr <<- NA_character_
  current_pass <<- 0L

  calls[[current_call]] <<- data.frame(
    call = current_call,
    start_line = line_number,
    finish_line = NA_integer_,
    study = study,
    seed = seed,
    mode = mode,
    mode_source = mode_source,
    has_mm_signature = FALSE,
    stringsAsFactors = FALSE
  )
}

for (line_number in seq_along(lines)) {

  line = lines[[line_number]]

  if (grepl("Mapping FULL", line, ignore.case = TRUE) ||
      grepl("FULL mapping", line, ignore.case = TRUE)) {
    pending_mode = "f"
  }
  if (grepl("Mapping SIMPLE", line, ignore.case = TRUE) ||
      grepl("SIMPLE mapping", line, ignore.case = TRUE)) {
    pending_mode = "s"
  }

  if (grepl(
    "^Mapping QTLs for .+ \\([fs]\\)\\.\\.\\.$",
    line,
    perl = TRUE
  )) {

    pending_study = extract_first_match(
      "^Mapping QTLs for (.+) \\([fs]\\)\\.\\.\\.$",
      line
    )
    pending_mode = extract_first_match(
      "^Mapping QTLs for .+ \\(([fs])\\)\\.\\.\\.$",
      line
    )
  }

  if (grepl("^mapQTL: seed set to [0-9]+$", line)) {

    this_seed = suppressWarnings(
      as.integer(
        sub("^mapQTL: seed set to ([0-9]+)$", "\\1", line)
      )
    )
    seed_hit = names(study_seed)[study_seed == this_seed]
    seed_study = if (length(seed_hit) == 1L) {
      seed_hit[[1]]
    } else {
      NA_character_
    }

    this_study = if (!is.na(pending_study)) {
      if (!is.na(seed_study) && pending_study != seed_study) {
        stop("Conflicting study label and seed at log line ", line_number)
      }
      pending_study
    } else {
      seed_study
    }

    this_mode_source = if (!is.na(pending_mode)) {
      "explicit_log_context"
    } else {
      NA_character_
    }

    start_call(
      line_number = line_number,
      study = this_study,
      mode = pending_mode,
      mode_source = this_mode_source,
      seed = this_seed
    )

    pending_study = NA_character_
    pending_mode = NA_character_
    next
  }

  if (grepl("^Mapping:[[:space:]]+", line)) {

    if (is.na(current_call)) {
      start_call(
        line_number = line_number,
        study = pending_study,
        mode = pending_mode,
        mode_source = if (!is.na(pending_mode)) {
          "explicit_log_context"
        } else {
          NA_character_
        }
      )
      pending_study = NA_character_
    }

    current_chr = sub(
      "^Mapping:[[:space:]]+([^[:space:]]+).*$",
      "\\1",
      line
    )
    current_pass = 0L
  }

  if (!is.na(current_call) &&
      grepl("MM residualization", line, fixed = TRUE)) {
    calls[[current_call]]$has_mm_signature = TRUE
  }

  if (grepl("FDR Thresh 5%", line, fixed = TRUE)) {

    value = suppressWarnings(
      as.numeric(
        sub(
          ".*FDR Thresh 5%[[:space:]]+([-+0-9.eE]+).*$",
          "\\1",
          line
        )
      )
    )

    current_pass = current_pass + 1L

    records[[record_index]] = data.frame(
      line_number = line_number,
      call = current_call,
      study = current_study,
      seed = if (!is.na(current_call)) calls[[current_call]]$seed else NA_integer_,
      mode = current_mode,
      mode_source = current_mode_source,
      chromosome = current_chr,
      chromosome_number = normalize_chr(current_chr),
      pass = current_pass,
      lod = value,
      stringsAsFactors = FALSE
    )
    record_index = record_index + 1L
  }

  if (grepl("^mapQTL: total elapsed = ", line)) {
    if (!is.na(current_call)) {
      calls[[current_call]]$finish_line = line_number
    }
    current_call = NA_integer_
    current_study = NA_character_
    current_mode = NA_character_
    current_mode_source = NA_character_
    current_chr = NA_character_
    current_pass = 0L
  }
}

if (length(records) == 0L) {
  stop(
    "No 'FDR Thresh 5%' records were found in ",
    threshold_log
  )
}

threshold_records_all = do.call(rbind, records)
threshold_calls = if (length(calls) > 0L) {
  do.call(rbind, calls)
} else {
  data.frame()
}

# Historical seed_from collapses unsigned hashes above the signed integer limit
# to seed 1. Resolve collisions using independently printed phenotype dimensions
# and the final MM residualization column, never the call's position in the log.
logged_genes = setNames(rep(NA_integer_, length(all_studies)), all_studies)
for (study in all_studies) {
  headers = which(lines == paste0("$", study))
  values = integer()
  for (h in headers[headers < min(threshold_calls$start_line)]) {
    if (h < length(lines) && grepl("^\\[1\\][[:space:]]+[0-9]+[[:space:]]+[0-9]+[[:space:]]*$", lines[[h + 1L]])) {
      values = c(values, as.integer(tail(strsplit(trimws(sub("^\\[1\\]", "", lines[[h + 1L]])), "[[:space:]]+")[[1]], 1)))
    }
  }
  if (length(unique(values)) != 1L) stop("Ambiguous logged dimensions: ", study)
  logged_genes[[study]] = unique(values)
}
threshold_calls$identity_source = ifelse(is.na(threshold_calls$study), NA_character_, "unique_seed_or_label")
threshold_calls$n_genes_logged = NA_integer_
for (i in which(threshold_calls$has_mm_signature)) {
  if (is.na(threshold_calls$finish_line[[i]])) stop("Incomplete FULL call.")
  segment = lines[threshold_calls$start_line[[i]]:threshold_calls$finish_line[[i]]]
  mm = segment[grepl("MM residualization batch", segment, fixed = TRUE)]
  endpoints = suppressWarnings(as.integer(sub(".*cols [0-9]+[^0-9]+([0-9]+)\\).*", "\\1", mm)))
  if (!length(endpoints) || anyNA(endpoints)) stop("Cannot parse MM column counts.")
  n_genes = max(endpoints)
  candidates = all_studies[study_seed == threshold_calls$seed[[i]] & logged_genes == n_genes]
  if (length(candidates) != 1L) stop("Seed and gene count do not uniquely identify call ", i)
  previous = threshold_calls$study[[i]]
  if (!is.na(previous) && previous != candidates[[1]]) stop("Conflicting study identity.")
  threshold_calls$study[[i]] = candidates[[1]]
  threshold_calls$identity_source[[i]] = "seed_and_logged_gene_count"
  threshold_calls$n_genes_logged[[i]] = n_genes
}
write_table(data.frame(study = all_studies, seed = study_seed,
                       n_genes = logged_genes), "threshold_identity_evidence.tsv")

# A FULL map uses the LOCO mixed-model residualization. When an older wrapper
# did not print an explicit FULL/SIMPLE label, the mapQTL MM log signature
# identifies FULL calls without relying on run order.
if (nrow(threshold_calls) > 0L) {
  infer_full = is.na(threshold_calls$mode) &
    threshold_calls$has_mm_signature
  threshold_calls$mode[infer_full] = "f"
  threshold_calls$mode_source[infer_full] = "MM_residualization_signature"

  call_match = match(
    threshold_records_all$call,
    threshold_calls$call
  )
  threshold_records_all$study = threshold_calls$study[call_match]
  threshold_records_all$mode = threshold_calls$mode[call_match]
  threshold_records_all$mode_source = threshold_calls$mode_source[call_match]
}

threshold_records_all$assigned =
  !is.na(threshold_records_all$study) &
  !is.na(threshold_records_all$mode) &
  !is.na(threshold_records_all$chromosome_number) &
  is.finite(threshold_records_all$lod)

# Always write raw diagnostics before selecting calls.
context_index = which(grepl("FDR Thresh 5%", lines, fixed = TRUE))
threshold_context = data.frame(
  line_number = context_index,
  previous_3 = vapply(
    context_index,
    function(i) paste(lines[pmax(1, i - 3):(i - 1)], collapse = " | "),
    character(1)
  ),
  threshold_line = lines[context_index],
  next_2 = vapply(
    context_index,
    function(i) paste(lines[(i + 1):pmin(length(lines), i + 2)], collapse = " | "),
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_table(threshold_calls, "threshold_mapqtl_calls.tsv")
write_table(threshold_records_all, "threshold_records_all.tsv")
write_table(threshold_context, "threshold_context.tsv")

###############################
# SELECT FULL-MAP CALLS ONLY  #
###############################

full_calls = threshold_calls[
  threshold_calls$mode == "f" &
  threshold_calls$study %in% analysis_studies,
  ,
  drop = FALSE
]

if (nrow(full_calls) == 0L) {
  stop(
    "No FULL mapQTL calls could be assigned. Raw diagnostics were written to ",
    power_asys,
    "threshold_mapqtl_calls.tsv, threshold_records_all.tsv, and threshold_context.tsv."
  )
}

# If a study appears in more than one historical FULL call, use its latest
# call. This avoids mixing repeated reruns of the same mapping job.
selected_calls = do.call(
  rbind,
  lapply(
    analysis_studies,
    function(study) {
      candidates = full_calls[
        full_calls$study == study,
        ,
        drop = FALSE
      ]
      if (nrow(candidates) == 0L) return(NULL)
      candidates[which.max(candidates$call), , drop = FALSE]
    }
  )
)

if (is.null(selected_calls) || nrow(selected_calls) == 0L) {
  stop("Could not select any FULL mapQTL calls.")
}

threshold_records = threshold_records_all[
  threshold_records_all$call %in% selected_calls$call &
  threshold_records_all$assigned &
  threshold_records_all$mode == "f",
  ,
  drop = FALSE
]

missing_studies = setdiff(
  analysis_studies,
  unique(threshold_records$study)
)

############################
# SUMMARIZE THRESHOLD SHAPE #
############################

by_study_all_pass = do.call(
  rbind,
  lapply(
    analysis_studies,
    function(study) {
      dat = threshold_records[
        threshold_records$study == study,
        ,
        drop = FALSE
      ]
      out = summarize_thresholds(dat$lod)
      out$study = study
      out$n_chromosomes = length(unique(dat$chromosome_number))
      out$n_passes_max = if (nrow(dat) > 0L) max(dat$pass) else 0L
      out
    }
  )
)

first_pass = threshold_records[
  threshold_records$pass == 1L,
  ,
  drop = FALSE
]

by_study_first_pass = do.call(
  rbind,
  lapply(
    analysis_studies,
    function(study) {
      dat = first_pass[first_pass$study == study, , drop = FALSE]
      out = summarize_thresholds(dat$lod)
      out$study = study
      out$n_chromosomes = length(unique(dat$chromosome_number))
      out
    }
  )
)

first_pass_grid = expand.grid(
  study = analysis_studies,
  chromosome_number = 1:16,
  stringsAsFactors = FALSE
)
first_pass_lookup = first_pass[
  ,
  c(
    "study",
    "chromosome_number",
    "chromosome",
    "lod",
    "line_number",
    "call",
    "mode_source"
  ),
  drop = FALSE
]

if (any(duplicated(first_pass_lookup[, c("study", "chromosome_number")]))) {
  stop(
    "More than one selected first-pass threshold was recovered for a study/chromosome. ",
    "Inspect threshold_full_first_pass_records.tsv."
  )
}

first_pass_by_study_chr = merge(
  first_pass_grid,
  first_pass_lookup,
  by = c("study", "chromosome_number"),
  all.x = TRUE,
  sort = FALSE
)
first_pass_by_study_chr = first_pass_by_study_chr[
  order(
    match(first_pass_by_study_chr$study, analysis_studies),
    first_pass_by_study_chr$chromosome_number
  ),
  ,
  drop = FALSE
]
first_pass_by_study_chr$available =
  is.finite(first_pass_by_study_chr$lod)

########################
# WRITE AUDIT OUTPUTS  #
########################

write_table(selected_calls, "threshold_selected_full_calls.tsv")
write_table(threshold_records, "threshold_full_records.tsv")
write_table(first_pass, "threshold_full_first_pass_records.tsv")
write_table(by_study_all_pass, "threshold_full_all_pass_by_study.tsv")
write_table(by_study_first_pass, "threshold_full_first_pass_by_study.tsv")
write_table(first_pass_by_study_chr, "threshold_full_first_pass_by_study_chr.tsv")

########################
# FINAL AUDIT SUMMARY  #
########################

cat("\n")
cat("FULL-map permutation-threshold audit\n")
cat("Log: ", threshold_log, "\n", sep = "")
cat(
  "Total FDR-threshold lines: ",
  nrow(threshold_records_all),
  "\n",
  sep = ""
)
cat(
  "Assigned FULL-map threshold lines: ",
  nrow(threshold_records),
  "\n",
  sep = ""
)
cat(
  "Assigned FULL first-pass thresholds: ",
  nrow(first_pass),
  " / ",
  length(analysis_studies) * 16L,
  " possible study/chromosome cells\n",
  sep = ""
)

cat("\nSelected FULL calls:\n")
print(selected_calls, row.names = FALSE)

cat("\nPer-study FIRST-PASS threshold summary:\n")
print(by_study_first_pass, row.names = FALSE)

if (length(missing_studies) > 0L) {
  stop(
    "No usable FULL-map thresholds recovered for: ",
    paste(missing_studies, collapse = ", "),
    ". Inspect threshold_mapqtl_calls.tsv and threshold_context.tsv."
  )
}

missing_first_pass = first_pass_by_study_chr[
  !first_pass_by_study_chr$available,
  c("study", "chromosome_number"),
  drop = FALSE
]

if (nrow(missing_first_pass) > 0L) {
  stop(
    nrow(missing_first_pass),
    " study/chromosome first-pass thresholds were not recoverable. ",
    "Do not run the primary power model until these cells are reviewed."
  )
}

message(
  "Threshold audit complete. The primary power model should use recovered ",
  "FULL-map first-pass thresholds by study and chromosome rather than a pooled ",
  "median or q10/q90 threshold summary."
)
