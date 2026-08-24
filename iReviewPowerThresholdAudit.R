########################################
########################################
## GENETICS REVIEW POWER THRESHOLDS   ##
## FULL-MAP PERMUTATION AUDIT         ##
########################################
########################################

#################
# PRELIMINARIES #
#################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
power_asys = paste0(review_asys, "power/")

ensure_directory(review_asys)
ensure_directory(power_asys)

threshold_log = paste0(logs, "iRemapAll.out")

if (!file.exists(threshold_log)) {
  stop("Missing remapping log: ", threshold_log)
}

analysis_studies = setdiff(
  classicMetadata$study,
  "AlbertAlt_2018"
)

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

########################
# PARSE REMAPPING LOG  #
########################

lines = readLines(threshold_log, warn = FALSE)

current_study = NA_character_
current_method = NA_character_
current_chr = NA_character_
current_block = NA_integer_
block_counter = 0L

blocks = list()
block_index = 1L
records = list()
record_index = 1L

for (line_number in seq_along(lines)) {

  line = lines[[line_number]]

  if (grepl(
    "^Mapping QTLs for .+ \\([fs]\\)\\.\\.\\.$",
    line,
    perl = TRUE
  )) {

    study = extract_first_match(
      "^Mapping QTLs for (.+) \\([fs]\\)\\.\\.\\.$",
      line
    )
    method = extract_first_match(
      "^Mapping QTLs for .+ \\(([fs])\\)\\.\\.\\.$",
      line
    )

    block_counter = block_counter + 1L
    current_block = block_counter
    current_study = study
    current_method = method
    current_chr = NA_character_

    blocks[[block_index]] = data.frame(
      block = current_block,
      study = current_study,
      method = current_method,
      start_line = line_number,
      finish_line = NA_integer_,
      completed = FALSE,
      stringsAsFactors = FALSE
    )
    block_index = block_index + 1L
    next
  }

  if (grepl("^Mapping:[[:space:]]+", line)) {
    current_chr = sub(
      "^Mapping:[[:space:]]+([^[:space:]]+).*$",
      "\\1",
      line
    )
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

    records[[record_index]] = data.frame(
      line_number = line_number,
      block = current_block,
      study = current_study,
      method = current_method,
      chromosome = current_chr,
      lod = value,
      stringsAsFactors = FALSE
    )
    record_index = record_index + 1L
  }

  if (grepl("^Finished mapping QTLs for ", line)) {

    finished_study = extract_first_match(
      "^Finished mapping QTLs for (.+) in [0-9.]+ s$",
      line
    )

    if (
      is.finite(current_block) &&
      !is.na(finished_study) &&
      identical(finished_study, current_study)
    ) {
      hit = which(vapply(
        blocks,
        function(x) x$block[[1]] == current_block,
        logical(1)
      ))
      if (length(hit) == 1L) {
        blocks[[hit]]$finish_line = line_number
        blocks[[hit]]$completed = TRUE
      }
    }

    current_study = NA_character_
    current_method = NA_character_
    current_chr = NA_character_
    current_block = NA_integer_
  }
}

if (length(records) == 0L) {
  stop(
    "No 'FDR Thresh 5%' records were found in ",
    threshold_log
  )
}

threshold_records_all = do.call(rbind, records)
threshold_blocks = if (length(blocks) > 0L) {
  do.call(rbind, blocks)
} else {
  data.frame(
    block = integer(0),
    study = character(0),
    method = character(0),
    start_line = integer(0),
    finish_line = integer(0),
    completed = logical(0),
    stringsAsFactors = FALSE
  )
}

threshold_records_all$assigned =
  !is.na(threshold_records_all$study) &
  !is.na(threshold_records_all$method) &
  !is.na(threshold_records_all$chromosome) &
  is.finite(threshold_records_all$lod)

###############################
# SELECT FULL-MAP BLOCKS ONLY #
###############################

full_blocks = threshold_blocks[
  threshold_blocks$method == "f" &
  threshold_blocks$study %in% analysis_studies,
  ,
  drop = FALSE
]

if (nrow(full_blocks) == 0L) {
  write_table(
    threshold_records_all,
    "threshold_records_all.tsv"
  )
  write_table(
    threshold_blocks,
    "threshold_mapping_blocks.tsv"
  )
  stop(
    "No structured FULL-map blocks were recovered. ",
    "Audit files were written for inspection."
  )
}

selected_blocks = do.call(
  rbind,
  lapply(
    analysis_studies,
    function(study) {

      candidates = full_blocks[
        full_blocks$study == study,
        ,
        drop = FALSE
      ]

      if (nrow(candidates) == 0L) {
        return(NULL)
      }

      completed = candidates[candidates$completed, , drop = FALSE]
      if (nrow(completed) > 0L) {
        candidates = completed
      }

      candidates[
        which.max(candidates$block),
        ,
        drop = FALSE
      ]
    }
  )
)

if (is.null(selected_blocks) || nrow(selected_blocks) == 0L) {
  stop("Could not select any FULL-map threshold blocks.")
}

threshold_records = threshold_records_all[
  threshold_records_all$block %in% selected_blocks$block &
  threshold_records_all$assigned,
  ,
  drop = FALSE
]

threshold_records = threshold_records[
  threshold_records$study %in% analysis_studies &
  threshold_records$method == "f",
  ,
  drop = FALSE
]

threshold_records$pass = ave(
  threshold_records$line_number,
  threshold_records$study,
  threshold_records$chromosome,
  FUN = seq_along
)
threshold_records$pass = as.integer(threshold_records$pass)

missing_studies = setdiff(
  analysis_studies,
  unique(threshold_records$study)
)

############################
# SUMMARIZE THRESHOLD SHAPE #
############################

by_study = do.call(
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
      out$n_chromosomes = length(unique(dat$chromosome))
      out$n_passes_max = if (nrow(dat) > 0L) {
        max(dat$pass)
      } else {
        0L
      }
      out$selected_block = if (study %in% selected_blocks$study) {
        selected_blocks$block[match(study, selected_blocks$study)]
      } else {
        NA_integer_
      }
      out$block_completed = if (study %in% selected_blocks$study) {
        selected_blocks$completed[match(study, selected_blocks$study)]
      } else {
        NA
      }
      out
    }
  )
)

by_study = by_study[
  match(analysis_studies, by_study$study),
  c(
    "study",
    "n_thresholds",
    "n_chromosomes",
    "n_passes_max",
    "minimum",
    "q10",
    "q25",
    "median",
    "mean",
    "q75",
    "q90",
    "maximum",
    "sd",
    "iqr",
    "selected_block",
    "block_completed"
  ),
  drop = FALSE
]

study_chr_groups = split(
  threshold_records,
  interaction(
    threshold_records$study,
    threshold_records$chromosome,
    drop = TRUE
  )
)

by_study_chr = do.call(
  rbind,
  lapply(
    study_chr_groups,
    function(dat) {
      summary = summarize_thresholds(dat$lod)
      data.frame(
        study = dat$study[[1]],
        chromosome = dat$chromosome[[1]],
        n_passes = nrow(dat),
        minimum = summary$minimum,
        median = summary$median,
        maximum = summary$maximum,
        stringsAsFactors = FALSE
      )
    }
  )
)

block_threshold_counts = merge(
  threshold_blocks,
  aggregate(
    lod ~ block,
    data = threshold_records_all[
      is.finite(threshold_records_all$lod),
      ,
      drop = FALSE
    ],
    FUN = length
  ),
  by = "block",
  all.x = TRUE,
  sort = FALSE
)
names(block_threshold_counts)[
  names(block_threshold_counts) == "lod"
] = "n_thresholds"
block_threshold_counts$n_thresholds[
  is.na(block_threshold_counts$n_thresholds)
] = 0L

########################
# WRITE AUDIT OUTPUTS  #
########################

write_table(
  threshold_records_all,
  "threshold_records_all.tsv"
)
write_table(
  threshold_blocks,
  "threshold_mapping_blocks.tsv"
)
write_table(
  block_threshold_counts,
  "threshold_mapping_block_counts.tsv"
)
write_table(
  selected_blocks,
  "threshold_selected_full_blocks.tsv"
)
write_table(
  threshold_records,
  "threshold_full_records.tsv"
)
write_table(
  by_study,
  "threshold_full_by_study.tsv"
)
write_table(
  by_study_chr,
  "threshold_full_by_study_chromosome.tsv"
)

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
  "Unassigned threshold lines: ",
  sum(!threshold_records_all$assigned),
  "\n",
  sep = ""
)
cat(
  "Selected FULL-map records: ",
  nrow(threshold_records),
  "\n",
  sep = ""
)
cat("\nPer-study FULL-map threshold summary:\n")
print(by_study, row.names = FALSE)

if (length(missing_studies) > 0L) {
  stop(
    "No usable FULL-map thresholds recovered for: ",
    paste(missing_studies, collapse = ", "),
    ". Inspect threshold_mapping_blocks.tsv and threshold_records_all.tsv."
  )
}

if (any(!by_study$block_completed)) {
  warning(
    "At least one selected FULL-map block lacks a recorded completion line; ",
    "inspect threshold_selected_full_blocks.tsv before using its thresholds."
  )
}

message(
  "Threshold audit complete. Primary candidate cutoff is the per-study ",
  "median; sensitivity bounds should be chosen only after inspecting the ",
  "actual per-study threshold counts and ranges."
)
