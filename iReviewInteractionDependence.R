########################################
########################################
# REVIEWER SEGREGANT DEPENDENCE       #
# TARGETED JOINT-QTL SENSITIVITY      #
########################################
########################################

source("/PHShome/kv1026/p05/remap_functions.R")

review_asys = paste0(asys, "review/")
dependence_asys = paste0(review_asys, "interaction_dependence/")
ensure_directory(review_asys)
ensure_directory(dependence_asys)

cores = suppressWarnings(
  as.integer(
    Sys.getenv(
      "REVIEW_DEPENDENCE_CORES",
      Sys.getenv("SLURM_CPUS_PER_TASK", "8")
    )
  )
)
if (is.na(cores) || cores < 1L) cores = 8L

message("Shared-segregant sensitivity using ", cores, " worker(s).")

####################
# HELPER FUNCTIONS #
####################

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
    stop(
      "Expected one object in ", path,
      "; found: ", paste(loaded, collapse = ", ")
    )
  }

  env[[loaded[[1]]]]
}

write_tsv = function(x, filename) {
  write.table(
    x,
    file = paste0(dependence_asys, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

parallel_lapply = function(X, FUN) {
  if (.Platform$OS.type == "unix" && cores > 1L && length(X) > 1L) {
    return(
      parallel::mclapply(
        X,
        FUN,
        mc.cores = min(cores, length(X)),
        mc.preschedule = TRUE
      )
    )
  }
  lapply(X, FUN)
}

canonical_segregant = function(x) {
  x = trimws(as.character(x))
  x = gsub("_", ".", x, fixed = TRUE)
  gsub("\\s+", "", x)
}

interaction_term = function(beta_names) {
  term = grep(
    "(^geno:analyte|^analyte.*:geno$)",
    beta_names,
    value = TRUE
  )
  if (length(term) != 1L) {
    stop(
      "Expected one genotype-by-analyte coefficient; found: ",
      paste(term, collapse = ", ")
    )
  }
  term[[1]]
}

fit_segregant_random_lrt = function(model_frame, segregant) {

  dat = model_frame
  dat$segregant = factor(segregant)
  dat$study = factor(dat$study)
  dat$analyte = factor(dat$analyte)

  f0 = stats::as.formula(
    "pheno ~ (1 | study) + analyte + geno + (1 | segregant)"
  )
  f1 = stats::as.formula(
    "pheno ~ (1 | study) + analyte + geno + geno:analyte + (1 | segregant)"
  )

  ctrl = lme4::lmerControl(
    check.conv.singular = "ignore",
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )

  warning_text = character(0)

  fit_result = tryCatch(
    withCallingHandlers(
      {
        m0 = lme4::lmer(
          f0,
          data = dat,
          control = ctrl,
          REML = FALSE
        )
        m1 = lme4::lmer(
          f1,
          data = dat,
          control = ctrl,
          REML = FALSE
        )
        list(m0 = m0, m1 = m1)
      },
      warning = function(w) {
        warning_text <<- c(warning_text, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit_result, "error")) {
    return(
      data.frame(
        segregant_random_fit_success = FALSE,
        segregant_random_warning = paste(unique(warning_text), collapse = " | "),
        segregant_random_error = conditionMessage(fit_result),
        segregant_random_interaction_effect = NA_real_,
        segregant_random_interaction_p = NA_real_,
        segregant_random_lod = NA_real_,
        segregant_random_null_singular = NA,
        segregant_random_full_singular = NA,
        segregant_random_variance = NA_real_,
        segregant_random_study_variance = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  m0 = fit_result$m0
  m1 = fit_result$m1
  term = interaction_term(names(lme4::fixef(m1)))

  lrt = tryCatch(
    stats::anova(m0, m1),
    error = function(e) NULL
  )

  p = if (is.null(lrt)) NA_real_ else as.numeric(lrt[2, "Pr(>Chisq)"])
  lod = (
    as.numeric(stats::logLik(m1)) -
    as.numeric(stats::logLik(m0))
  ) / log(10)

  vc = as.data.frame(lme4::VarCorr(m1))
  segregant_variance = vc$vcov[
    vc$grp == "segregant" & vc$var1 == "(Intercept)"
  ]
  study_variance = vc$vcov[
    vc$grp == "study" & vc$var1 == "(Intercept)"
  ]

  data.frame(
    segregant_random_fit_success = is.finite(p),
    segregant_random_warning = paste(unique(warning_text), collapse = " | "),
    segregant_random_error = if (is.finite(p)) "" else "LRT p-value was not finite",
    segregant_random_interaction_effect = unname(lme4::fixef(m1)[term]),
    segregant_random_interaction_p = p,
    segregant_random_lod = lod,
    segregant_random_null_singular = lme4::isSingular(m0, tol = 1e-5),
    segregant_random_full_singular = lme4::isSingular(m1, tol = 1e-5),
    segregant_random_variance = if (
      length(segregant_variance) == 1L
    ) segregant_variance else NA_real_,
    segregant_random_study_variance = if (
      length(study_variance) == 1L
    ) study_variance else NA_real_,
    stringsAsFactors = FALSE
  )
}

########################
# LOAD JOINT DATA      #
########################

metadata_analyte = scaMetadata

studies_analyte = lapply(
  metadata_analyte$study,
  function(study) {
    load_single_rda(
      paste0(imed, study, ".rda"),
      object_name = study
    )
  }
)
names(studies_analyte) = metadata_analyte$study

metadata_analyte$genes = lapply(
  metadata_analyte$study,
  function(study) colnames(studies_analyte[[study]]$pheno)
)
metadata_analyte$segs = lapply(
  metadata_analyte$study,
  function(study) rownames(studies_analyte[[study]]$pheno)
)

supported_analyte = suppTwo(metadata_analyte)
unionBG_analyte = load_single_rda(
  paste0(imed, "unionAnalyteBG.rda"),
  object_name = "unionAnalyteBG"
)

jBH_analyte_geno = load_single_rda(
  paste0(imed, "jBH_analyte_geno.rda"),
  object_name = "jBH_analyte_geno"
)
jBH_analyte_ixn = load_single_rda(
  paste0(imed, "jBH_analyte_ixn.rda"),
  object_name = "jBH_analyte_ixn"
)
qtls_aug_analyte_geno = load_single_rda(
  paste0(imed, "qtls_aug_analyte_geno.rda"),
  object_name = "qtls_aug_analyte_geno"
)
qtls_aug_analyte_ixn = load_single_rda(
  paste0(imed, "qtls_aug_analyte_ixn.rda"),
  object_name = "qtls_aug_analyte_ixn"
)

if (nrow(jBH_analyte_geno) != 868L) {
  stop("Expected 868 joint QTLs; found ", nrow(jBH_analyte_geno), ".")
}
if (nrow(jBH_analyte_ixn) != 64L) {
  stop("Expected 64 strict interaction QTLs; found ", nrow(jBH_analyte_ixn), ".")
}

joint_keys = paste(jBH_analyte_geno$gene, jBH_analyte_geno$marker, sep = "::")
strict_keys = paste(jBH_analyte_ixn$gene, jBH_analyte_ixn$marker, sep = "::")

target_pairs = unique(
  rbind(
    jBH_analyte_geno[, c("gene", "marker"), drop = FALSE],
    jBH_analyte_ixn[, c("gene", "marker"), drop = FALSE]
  )
)
target_pairs$gene = as.character(target_pairs$gene)
target_pairs$marker = as.character(target_pairs$marker)

target_keys = paste(target_pairs$gene, target_pairs$marker, sep = "::")
message("Targeted loci: ", nrow(target_pairs))

########################################
# DOCUMENT SEGREGANT REUSE             #
########################################

study_segregant_ids = lapply(
  studies_analyte,
  function(study) canonical_segregant(rownames(study$gdata))
)

if (any(vapply(study_segregant_ids, anyDuplicated, integer(1)) > 0L)) {
  stop("At least one joint-analysis dataset has duplicated canonical segregant IDs.")
}

pair_index = utils::combn(names(study_segregant_ids), 2, simplify = FALSE)
segregant_panel_overlap = do.call(
  rbind,
  lapply(
    pair_index,
    function(pair) {
      x = study_segregant_ids[[pair[[1]]]]
      y = study_segregant_ids[[pair[[2]]]]
      shared = intersect(x, y)
      data.frame(
        study_1 = pair[[1]],
        study_2 = pair[[2]],
        n_study_1 = length(x),
        n_study_2 = length(y),
        n_shared = length(shared),
        fraction_study_1_shared = length(shared) / length(x),
        fraction_study_2_shared = length(shared) / length(y),
        stringsAsFactors = FALSE
      )
    }
  )
)

if (max(segregant_panel_overlap$n_shared) == 0L) {
  stop("No segregant reuse detected among the joint-analysis datasets.")
}

########################################
# REFIT ORIGINAL TARGETED LOCI         #
########################################

message("Refitting targeted loci with the original joint model...")

modelQTL = drillQTL(
  studies = studies_analyte,
  dfQTL = target_pairs,
  unionBg = unionBG_analyte,
  smetadata = metadata_analyte,
  suppGenes = supported_analyte,
  fixedEffect = "analyte",
  correctBG = TRUE,
  correctGRM = FALSE,
  cores = cores,
  chromosomes = NULL
)

if (length(modelQTL) != nrow(target_pairs)) {
  stop(
    "Expected ", nrow(target_pairs),
    " targeted refits; received ", length(modelQTL), "."
  )
}

########################################
# RECONSTRUCT SEGREGANT IDS PER MODEL  #
########################################

reconstruct_rows = function(gene, marker, model_frame) {

  support = metadata_analyte$study[
    vapply(
      metadata_analyte$genes,
      function(x) gene %in% x,
      logical(1)
    )
  ]
  support = support[support %in% names(studies_analyte)]

  study_expected = unlist(
    lapply(
      support,
      function(study) {
        rep(study, nrow(studies_analyte[[study]]$gdata))
      }
    ),
    use.names = FALSE
  )

  segregant_expected = unlist(
    lapply(
      support,
      function(study) {
        canonical_segregant(rownames(studies_analyte[[study]]$gdata))
      }
    ),
    use.names = FALSE
  )

  genotype_raw = unlist(
    lapply(
      support,
      function(study) studies_analyte[[study]]$gdata[, marker]
    ),
    use.names = FALSE
  )
  genotype_expected = as.numeric(scale(genotype_raw))

  if (length(study_expected) != nrow(model_frame)) {
    stop(
      "Row-count mismatch for ", gene, " / ", marker,
      ": reconstructed ", length(study_expected),
      ", model ", nrow(model_frame)
    )
  }

  if (!all(as.character(model_frame$study) == study_expected)) {
    stop("Study-order mismatch for ", gene, " / ", marker, ".")
  }

  genotype_difference = max(
    abs(as.numeric(model_frame$geno) - genotype_expected),
    na.rm = TRUE
  )

  if (!is.finite(genotype_difference) || genotype_difference > 1e-8) {
    stop(
      "Genotype-order mismatch for ", gene, " / ", marker,
      "; maximum difference = ", genotype_difference
    )
  }

  list(
    support = support,
    segregant = segregant_expected,
    genotype_difference = genotype_difference
  )
}

########################################
# LOCUS-LEVEL SENSITIVITY              #
########################################

extract_locus = function(x) {

  gene = as.character(x$gene)
  marker = as.character(x$marker)
  key = paste(gene, marker, sep = "::")

  m = x$interaction
  if (is.null(m) || !inherits(m, "lmerMod")) {
    stop("Missing original interaction model for ", key, ".")
  }

  mf = stats::model.frame(m)
  alignment = reconstruct_rows(gene, marker, mf)
  segregant = alignment$segregant

  segregant_counts = table(segregant)
  n_segregant_clusters = length(segregant_counts)
  n_reused_segregants = sum(segregant_counts > 1L)

  original_beta = lme4::fixef(m)
  original_term = interaction_term(names(original_beta))
  original_se = sqrt(diag(stats::vcov(m)))[original_term]

  segregant_random = fit_segregant_random_lrt(mf, segregant)

  data.frame(
    gene = gene,
    marker = marker,
    key = key,
    in_joint_868 = key %in% joint_keys,
    in_strict_64 = key %in% strict_keys,
    n_observations = nrow(mf),
    n_studies = length(alignment$support),
    n_segregant_clusters = n_segregant_clusters,
    n_reused_segregants = n_reused_segregants,
    reused_segregant_fraction = n_reused_segregants / n_segregant_clusters,
    maximum_observations_per_segregant = max(segregant_counts),
    maximum_genotype_alignment_difference = alignment$genotype_difference,
    original_interaction_effect = unname(original_beta[original_term]),
    original_interaction_se = unname(original_se),
    original_interaction_p = as.numeric(x$interaction_pval),
    original_model_singular = lme4::isSingular(m, tol = 1e-5),
    segregant_random,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

message("Calculating crossed-segregant random-effect sensitivity...")

locus_list = parallel_lapply(modelQTL, extract_locus)
segregant_dependence_locus_results = do.call(rbind, locus_list)
rownames(segregant_dependence_locus_results) = NULL

if (anyDuplicated(segregant_dependence_locus_results$key)) {
  stop("Targeted dependence output contains duplicated gene-marker keys.")
}
if (!setequal(segregant_dependence_locus_results$key, target_keys)) {
  stop("Targeted dependence output does not contain the exact requested loci.")
}

########################################
# VERIFY ORIGINAL REFITS               #
########################################

joint_reference = qtls_aug_analyte_geno[
  ,
  c("gene", "marker", "interaction_effect", "interaction_pval"),
  drop = FALSE
]
names(joint_reference)[3:4] = c(
  "stored_joint_interaction_effect",
  "stored_joint_interaction_p"
)

strict_reference = qtls_aug_analyte_ixn[
  ,
  c("gene", "marker", "interaction_effect", "interaction_pval"),
  drop = FALSE
]
names(strict_reference)[3:4] = c(
  "stored_strict_interaction_effect",
  "stored_strict_interaction_p"
)

segregant_dependence_locus_results = merge(
  segregant_dependence_locus_results,
  joint_reference,
  by = c("gene", "marker"),
  all.x = TRUE,
  sort = FALSE
)
segregant_dependence_locus_results = merge(
  segregant_dependence_locus_results,
  strict_reference,
  by = c("gene", "marker"),
  all.x = TRUE,
  sort = FALSE
)
segregant_dependence_locus_results$key = paste(
  segregant_dependence_locus_results$gene,
  segregant_dependence_locus_results$marker,
  sep = "::"
)
segregant_dependence_locus_results = segregant_dependence_locus_results[
  match(target_keys, segregant_dependence_locus_results$key),
  ,
  drop = FALSE
]

joint_check = segregant_dependence_locus_results[
  segregant_dependence_locus_results$in_joint_868,
  ,
  drop = FALSE
]
strict_check = segregant_dependence_locus_results[
  segregant_dependence_locus_results$in_strict_64,
  ,
  drop = FALSE
]

joint_p_error = max(
  abs(joint_check$original_interaction_p - joint_check$stored_joint_interaction_p),
  na.rm = TRUE
)
strict_p_error = max(
  abs(strict_check$original_interaction_p - strict_check$stored_strict_interaction_p),
  na.rm = TRUE
)

if (joint_p_error > 1e-8 || strict_p_error > 1e-8) {
  stop(
    "Targeted refits do not reproduce stored interaction p-values. ",
    "Joint max error = ", joint_p_error,
    "; strict max error = ", strict_p_error, "."
  )
}

if (sum(joint_check$original_interaction_p < 0.05) != 149L) {
  stop("Expected 149/868 nominal post-hoc interaction signals in the original model.")
}

########################################
# REVIEWER-FACING SUMMARY              #
########################################

summarize_set = function(label, dat) {

  random_ok = dat$segregant_random_fit_success &
    is.finite(dat$segregant_random_interaction_p)

  data.frame(
    set = label,
    n_loci = nrow(dat),
    median_segregant_clusters = stats::median(dat$n_segregant_clusters),
    median_reused_segregant_fraction = stats::median(dat$reused_segregant_fraction),
    original_nominal_p_lt_0_05 = sum(dat$original_interaction_p < 0.05),
    segregant_random_fit_success = sum(random_ok),
    segregant_random_full_singular = sum(
      dat$segregant_random_full_singular %in% TRUE,
      na.rm = TRUE
    ),
    segregant_random_nominal_p_lt_0_05 = sum(
      dat$segregant_random_interaction_p < 0.05,
      na.rm = TRUE
    ),
    original_and_segregant_random_nominal_p_lt_0_05 = sum(
      dat$original_interaction_p < 0.05 &
      dat$segregant_random_interaction_p < 0.05,
      na.rm = TRUE
    ),
    median_segregant_random_variance = stats::median(
      dat$segregant_random_variance[random_ok],
      na.rm = TRUE
    ),
    segregant_random_effect_pearson_r = stats::cor(
      dat$original_interaction_effect[random_ok],
      dat$segregant_random_interaction_effect[random_ok],
      method = "pearson"
    ),
    segregant_random_effect_spearman_rho = stats::cor(
      dat$original_interaction_effect[random_ok],
      dat$segregant_random_interaction_effect[random_ok],
      method = "spearman"
    ),
    segregant_random_sign_concordance = mean(
      sign(dat$original_interaction_effect[random_ok]) ==
      sign(dat$segregant_random_interaction_effect[random_ok])
    ),
    stringsAsFactors = FALSE
  )
}

segregant_dependence_summary = do.call(
  rbind,
  list(
    summarize_set(
      "all_unique_targeted_loci",
      segregant_dependence_locus_results
    ),
    summarize_set(
      "868_joint_QTLs",
      segregant_dependence_locus_results[
        segregant_dependence_locus_results$in_joint_868,
        ,
        drop = FALSE
      ]
    ),
    summarize_set(
      "64_strict_interaction_QTLs",
      segregant_dependence_locus_results[
        segregant_dependence_locus_results$in_strict_64,
        ,
        drop = FALSE
      ]
    )
  )
)
rownames(segregant_dependence_summary) = NULL

#################
# WRITE OUTPUTS #
#################

write_tsv(
  segregant_panel_overlap,
  "segregant_panel_overlap.tsv"
)
write_tsv(
  segregant_dependence_locus_results,
  "segregant_dependence_locus_results.tsv"
)
write_tsv(
  segregant_dependence_summary,
  "segregant_dependence_summary.tsv"
)

message("")
message("Shared-segregant sensitivity complete.")
message("Output directory: ", dependence_asys)
message("Original refit max p error, joint: ", signif(joint_p_error, 4))
message("Original refit max p error, strict: ", signif(strict_p_error, 4))
message("")
print(segregant_dependence_summary, row.names = FALSE)
