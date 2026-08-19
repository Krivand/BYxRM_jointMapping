########################################
########################################
# REVIEWER INTERACTION POWER ANALYSIS  #
########################################
########################################

source(
  "/PHShome/kv1026/p05/remap_functions.R"
)

review_asys =
  paste0(
    asys,
    "review/"
  )

interaction_asys =
  paste0(
    review_asys,
    "interaction_power/"
  )

interaction_imed =
  paste0(
    imed,
    "review/interaction_power/"
  )

ensure_directory(
  interaction_asys
)

ensure_directory(
  interaction_imed
)

cores =
  suppressWarnings(
    as.integer(
      Sys.getenv(
        "SLURM_CPUS_PER_TASK",
        "8"
      )
    )
  )

if (
  is.na(cores) ||
  cores < 1L
) {
  cores = 8L
}


########################################
# LOAD THE EXACT JOINT ANALYSIS INPUTS #
########################################

metadata_analyte =
  scaMetadata

studies_analyte =
  lapply(
    metadata_analyte$study,
    function(study) {

      env =
        new.env(
          parent = emptyenv()
        )

      load(
        paste0(
          imed,
          study,
          ".rda"
        ),
        envir = env
      )

      env[[study]]
    }
  )

names(
  studies_analyte
) =
  metadata_analyte$study

metadata_analyte$genes =
  lapply(
    metadata_analyte$study,
    function(study) {
      colnames(
        studies_analyte[[
          study
        ]]$pheno
      )
    }
  )

metadata_analyte$segs =
  lapply(
    metadata_analyte$study,
    function(study) {
      rownames(
        studies_analyte[[
          study
        ]]$pheno
      )
    }
  )

supported_analyte =
  suppTwo(
    metadata_analyte
  )

load(
  paste0(
    imed,
    "unionAnalyteBG.rda"
  )
)

unionBG_analyte =
  unionAnalyteBG

rm(
  unionAnalyteBG
)

load(
  paste0(
    imed,
    "jBH_analyte_geno.rda"
  )
)

load(
  paste0(
    imed,
    "jBH_analyte_ixn.rda"
  )
)

load(
  paste0(
    imed,
    "qtls_aug_analyte_geno.rda"
  )
)

if (
  nrow(
    jBH_analyte_geno
  ) != 868L
) {
  stop(
    "Expected 868 discovered joint genotype QTLs; found ",
    nrow(
      jBH_analyte_geno
    )
  )
}


########################################
# TARGETED REFIT OF THE 868 JOINT QTLS #
########################################

message(
  "Refitting interaction models at 868 known joint-QTL loci..."
)

modelQTL =
  drillQTL(
    studies =
      studies_analyte,
    dfQTL =
      jBH_analyte_geno,
    unionBg =
      unionBG_analyte,
    smetadata =
      metadata_analyte,
    suppGenes =
      supported_analyte,
    fixedEffect =
      "analyte",
    correctBG =
      TRUE,
    correctGRM =
      FALSE,
    cores =
      cores,
    chromosomes =
      NULL
  )


########################################
# EXTRACT INTERACTION INFORMATION      #
########################################

extract_one = function(x) {

  gene =
    if (!is.null(x$gene)) {
      as.character(
        x$gene
      )
    } else {
      NA_character_
    }

  marker =
    if (!is.null(x$marker)) {
      as.character(
        x$marker
      )
    } else {
      NA_character_
    }

  error =
    if (!is.null(x$error)) {
      as.character(
        x$error
      )
    } else {
      NA_character_
    }

  m =
    x$interaction

  if (
    is.null(m) ||
    !inherits(
      m,
      "lmerMod"
    )
  ) {

    return(
      data.frame(
        gene = gene,
        marker = marker,
        error = error,
        interaction_term =
          NA_character_,
        interaction_effect =
          NA_real_,
        interaction_se =
          NA_real_,
        interaction_wald_z =
          NA_real_,
        interaction_wald_p =
          NA_real_,
        interaction_lrt_p =
          if (!is.null(
            x$interaction_pval
          )) {
            as.numeric(
              x$interaction_pval
            )
          } else {
            NA_real_
          },
        n_observations =
          NA_integer_,
        n_studies =
          NA_integer_,
        n_mRNA =
          NA_integer_,
        n_protein =
          NA_integer_,
        singular =
          NA,
        residual_sigma =
          NA_real_,
        convergence_message =
          NA_character_,
        stringsAsFactors = FALSE
      )
    )
  }

  beta =
    lme4::fixef(
      m
    )

  ixn_terms =
    grep(
      "(^geno:|:geno$)",
      names(beta),
      value = TRUE
    )

  if (
    length(ixn_terms) == 0L
  ) {

    interaction_term =
      NA_character_

    interaction_effect =
      NA_real_

    interaction_se =
      NA_real_

  } else {

    if (
      length(ixn_terms) == 1L
    ) {

      interaction_term =
        ixn_terms[[1]]

    } else {

      interaction_term =
        ixn_terms[
          which.max(
            abs(
              beta[
                ixn_terms
              ]
            )
          )
        ]
    }

    interaction_effect =
      unname(
        beta[
          interaction_term
        ]
      )

    vc =
      tryCatch(
        as.matrix(
          stats::vcov(
            m
          )
        ),
        error = function(e) {
          NULL
        }
      )

    if (
      is.null(vc) ||
      !interaction_term %in%
      rownames(vc)
    ) {

      interaction_se =
        NA_real_

    } else {

      interaction_se =
        sqrt(
          vc[
            interaction_term,
            interaction_term
          ]
        )
    }
  }

  if (
    is.finite(
      interaction_effect
    ) &&
    is.finite(
      interaction_se
    ) &&
    interaction_se > 0
  ) {

    interaction_wald_z =
      interaction_effect /
      interaction_se

    interaction_wald_p =
      2 *
      stats::pnorm(
        -abs(
          interaction_wald_z
        )
      )

  } else {

    interaction_wald_z =
      NA_real_

    interaction_wald_p =
      NA_real_
  }

  mf =
    stats::model.frame(
      m
    )

  analyte =
    if (
      "analyte" %in%
      colnames(mf)
    ) {
      as.character(
        mf$analyte
      )
    } else {
      rep(
        NA_character_,
        nrow(mf)
      )
    }

  study =
    if (
      "study" %in%
      colnames(mf)
    ) {
      as.character(
        mf$study
      )
    } else {
      rep(
        NA_character_,
        nrow(mf)
      )
    }

  conv =
    m@optinfo$conv$lme4$messages

  if (
    is.null(conv)
  ) {
    conv =
      ""
  } else {
    conv =
      paste(
        conv,
        collapse = "; "
      )
  }

  data.frame(
    gene = gene,
    marker = marker,
    error = error,
    interaction_term =
      interaction_term,
    interaction_effect =
      interaction_effect,
    interaction_se =
      interaction_se,
    interaction_wald_z =
      interaction_wald_z,
    interaction_wald_p =
      interaction_wald_p,
    interaction_lrt_p =
      if (!is.null(
        x$interaction_pval
      )) {
        as.numeric(
          x$interaction_pval
        )
      } else {
        NA_real_
      },
    n_observations =
      stats::nobs(
        m
      ),
    n_studies =
      length(
        unique(
          study[
            !is.na(study)
          ]
        )
      ),
    n_mRNA =
      sum(
        analyte == "mRNA",
        na.rm = TRUE
      ),
    n_protein =
      sum(
        analyte == "protein",
        na.rm = TRUE
      ),
    singular =
      lme4::isSingular(
        m,
        tol = 1e-4
      ),
    residual_sigma =
      stats::sigma(
        m
      ),
    convergence_message =
      conv,
    stringsAsFactors = FALSE
  )
}

interaction_design =
  do.call(
    rbind,
    lapply(
      modelQTL,
      extract_one
    )
  )


########################################
# VERIFY AGAINST EXISTING AUGMENTED QTL#
########################################

existing =
  qtls_aug_analyte_geno[
    ,
    c(
      "gene",
      "marker",
      "interaction_effect",
      "interaction_pval",
      "bin"
    ),
    drop = FALSE
  ]

names(existing)[
  names(existing) ==
  "interaction_effect"
] =
  "stored_interaction_effect"

names(existing)[
  names(existing) ==
  "interaction_pval"
] =
  "stored_interaction_pval"

interaction_design =
  merge(
    interaction_design,
    existing,
    by = c(
      "gene",
      "marker"
    ),
    all.x = TRUE,
    sort = FALSE
  )

interaction_design$effect_difference =
  interaction_design$interaction_effect -
  interaction_design$stored_interaction_effect

interaction_design$p_difference =
  interaction_design$interaction_lrt_p -
  interaction_design$stored_interaction_pval


########################################
# CALIBRATE LRT VS WALD APPROXIMATION  #
########################################

valid =
  is.finite(
    interaction_design$interaction_lrt_p
  ) &
  is.finite(
    interaction_design$interaction_wald_p
  ) &
  interaction_design$interaction_lrt_p > 0 &
  interaction_design$interaction_wald_p > 0

lrt_wald_calibration =
  data.frame(
    n =
      sum(valid),
    pearson_mlog10p =
      stats::cor(
        -log10(
          interaction_design$interaction_lrt_p[
            valid
          ]
        ),
        -log10(
          interaction_design$interaction_wald_p[
            valid
          ]
        ),
        method = "pearson"
      ),
    spearman_mlog10p =
      stats::cor(
        -log10(
          interaction_design$interaction_lrt_p[
            valid
          ]
        ),
        -log10(
          interaction_design$interaction_wald_p[
            valid
          ]
        ),
        method = "spearman"
      ),
    median_abs_effect_difference =
      stats::median(
        abs(
          interaction_design$effect_difference
        ),
        na.rm = TRUE
      ),
    median_abs_p_difference =
      stats::median(
        abs(
          interaction_design$p_difference
        ),
        na.rm = TRUE
      ),
    singular_fraction =
      mean(
        interaction_design$singular,
        na.rm = TRUE
      ),
    stringsAsFactors = FALSE
  )


########################################
# POWER AS FUNCTION OF TRUE IXN EFFECT #
########################################

strict_raw_p_boundary =
  max(
    jBH_analyte_ixn$p_value,
    na.rm = TRUE
  )

thresholds =
  data.frame(
    threshold =
      c(
        "nominal_0.05",
        "observed_5pct_FDR_raw_p_boundary"
      ),
    alpha =
      c(
        0.05,
        strict_raw_p_boundary
      ),
    stringsAsFactors = FALSE
  )

delta_grid =
  seq(
    0,
    0.60,
    by = 0.025
  )

valid_se =
  is.finite(
    interaction_design$interaction_se
  ) &
  interaction_design$interaction_se > 0

power_rows =
  lapply(
    seq_len(
      nrow(
        thresholds
      )
    ),
    function(t) {

      alpha =
        thresholds$alpha[[t]]

      critical =
        stats::qchisq(
          1 - alpha,
          df = 1
        )

      do.call(
        rbind,
        lapply(
          delta_grid,
          function(delta) {

            ncp =
              (
                delta /
                interaction_design$interaction_se[
                  valid_se
                ]
              )^2

            power =
              stats::pchisq(
                critical,
                df = 1,
                ncp = ncp,
                lower.tail = FALSE
              )

            data.frame(
              threshold =
                thresholds$threshold[[t]],
              alpha =
                alpha,
              true_abs_interaction =
                delta,
              n_loci =
                length(power),
              mean_power =
                mean(power),
              median_power =
                stats::median(power),
              p10_power =
                as.numeric(
                  stats::quantile(
                    power,
                    0.10
                  )
                ),
              p90_power =
                as.numeric(
                  stats::quantile(
                    power,
                    0.90
                  )
                ),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )

interaction_power =
  do.call(
    rbind,
    power_rows
  )


########################################
# MAGNITUDE x PREVALENCE SENSITIVITY   #
########################################

prevalence_grid =
  c(
    0.05,
    0.10,
    0.25,
    0.50,
    1.00
  )

interaction_prevalence =
  do.call(
    rbind,
    lapply(
      seq_len(
        nrow(
          interaction_power
        )
      ),
      function(i) {

        row =
          interaction_power[i, ]

        data.frame(
          threshold =
            row$threshold,
          alpha =
            row$alpha,
          true_abs_interaction =
            row$true_abs_interaction,
          prevalence =
            prevalence_grid,
          n_loci =
            row$n_loci,
          expected_true_discordant =
            prevalence_grid *
            row$n_loci,
          mean_power =
            row$mean_power,
          expected_true_detected =
            prevalence_grid *
            row$n_loci *
            row$mean_power,
          stringsAsFactors = FALSE
        )
      }
    )
  )


########################################
# OBSERVED EFFECT DISTRIBUTION          #
########################################

observed_interaction_summary =
  interaction_design %>%
  dplyr::summarise(
    n_loci =
      dplyr::n(),
    n_valid_se =
      sum(
        is.finite(
          interaction_se
        ) &
        interaction_se > 0
      ),
    median_abs_interaction =
      stats::median(
        abs(
          interaction_effect
        ),
        na.rm = TRUE
      ),
    q25_abs_interaction =
      as.numeric(
        stats::quantile(
          abs(
            interaction_effect
          ),
          0.25,
          na.rm = TRUE
        )
      ),
    q75_abs_interaction =
      as.numeric(
        stats::quantile(
          abs(
            interaction_effect
          ),
          0.75,
          na.rm = TRUE
        )
      ),
    median_se =
      stats::median(
        interaction_se,
        na.rm = TRUE
      ),
    q25_se =
      as.numeric(
        stats::quantile(
          interaction_se,
          0.25,
          na.rm = TRUE
        )
      ),
    q75_se =
      as.numeric(
        stats::quantile(
          interaction_se,
          0.75,
          na.rm = TRUE
        )
      ),
    nominal_p_lt_0.05 =
      sum(
        interaction_lrt_p < 0.05,
        na.rm = TRUE
      ),
    stringsAsFactors = FALSE
  )


########################################
# WRITE OUTPUTS                        #
########################################

write.table(
  interaction_design,
  file =
    paste0(
      interaction_asys,
      "interaction_design_868.tsv"
    ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  lrt_wald_calibration,
  file =
    paste0(
      interaction_asys,
      "lrt_wald_calibration.tsv"
    ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  thresholds,
  file =
    paste0(
      interaction_asys,
      "interaction_thresholds.tsv"
    ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  interaction_power,
  file =
    paste0(
      interaction_asys,
      "interaction_power_by_magnitude.tsv"
    ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  interaction_prevalence,
  file =
    paste0(
      interaction_asys,
      "interaction_power_by_magnitude_prevalence.tsv"
    ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  observed_interaction_summary,
  file =
    paste0(
      interaction_asys,
      "observed_interaction_summary.tsv"
    ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  list(
    interaction_design =
      interaction_design,
    lrt_wald_calibration =
      lrt_wald_calibration,
    thresholds =
      thresholds,
    interaction_power =
      interaction_power,
    interaction_prevalence =
      interaction_prevalence,
    observed_interaction_summary =
      observed_interaction_summary
  ),
  file =
    paste0(
      interaction_imed,
      "review_interaction_power.rds"
    )
)

rm(
  modelQTL
)

invisible(
  gc()
)

message("")
message(
  "Reviewer interaction-power analysis complete."
)
message(
  "Valid locus-specific SEs: ",
  sum(valid_se),
  " / ",
  nrow(
    interaction_design
  )
)
message(
  "Observed raw-p boundary among original 5% FDR discoveries: ",
  signif(
    strict_raw_p_boundary,
    5
  )
)
message(
  "Outputs: ",
  interaction_asys
)
