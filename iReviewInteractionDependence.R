########################################
########################################
# REVIEWER CLONE-DEPENDENCE ANALYSIS   #
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

dependence_asys =
  paste0(
    review_asys,
    "interaction_dependence/"
  )

dependence_imed =
  paste0(
    imed,
    "review/interaction_dependence/"
  )

ensure_directory(
  dependence_asys
)

ensure_directory(
  dependence_imed
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
# LOAD JOINT-ANALYSIS OBJECTS          #
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
        studies_analyte[[study]]$pheno
      )
    }
  )

metadata_analyte$segs =
  lapply(
    metadata_analyte$study,
    function(study) {
      rownames(
        studies_analyte[[study]]$pheno
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

load(
  paste0(
    imed,
    "qtls_aug_analyte_ixn.rda"
  )
)


########################################
# COMBINE UNIQUE TARGETED LOCI         #
########################################

joint_keys =
  paste(
    jBH_analyte_geno$gene,
    jBH_analyte_geno$marker,
    sep = "::"
  )

strict_keys =
  paste(
    jBH_analyte_ixn$gene,
    jBH_analyte_ixn$marker,
    sep = "::"
  )

target_pairs =
  unique(
    rbind(
      jBH_analyte_geno[
        ,
        c(
          "gene",
          "marker"
        ),
        drop = FALSE
      ],
      jBH_analyte_ixn[
        ,
        c(
          "gene",
          "marker"
        ),
        drop = FALSE
      ]
    )
  )

message(
  "Unique targeted loci: ",
  nrow(
    target_pairs
  )
)


########################################
# TARGETED MODEL REFITS                #
########################################

message(
  "Refitting known QTLs for clone-dependence sensitivity..."
)

modelQTL =
  drillQTL(
    studies =
      studies_analyte,
    dfQTL =
      target_pairs,
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
# CLONE-ID RECONSTRUCTION              #
########################################

canonical_clone = function(x) {

  x =
    trimws(
      as.character(x)
    )

  x =
    gsub(
      "_",
      ".",
      x,
      fixed = TRUE
    )

  gsub(
    "\\s+",
    "",
    x
  )
}


reconstruct_rows = function(
  gene,
  marker,
  model_frame
) {

  support =
    metadata_analyte$study[
      vapply(
        metadata_analyte$genes,
        function(x) {
          gene %in% x
        },
        logical(1)
      )
    ]

  support =
    support[
      support %in%
      names(
        studies_analyte
      )
    ]

  study_expected =
    unlist(
      lapply(
        support,
        function(study) {
          rep(
            study,
            nrow(
              studies_analyte[[study]]$gdata
            )
          )
        }
      ),
      use.names = FALSE
    )

  clone_expected =
    unlist(
      lapply(
        support,
        function(study) {
          canonical_clone(
            rownames(
              studies_analyte[[study]]$gdata
            )
          )
        }
      ),
      use.names = FALSE
    )

  genotype_raw =
    unlist(
      lapply(
        support,
        function(study) {
          studies_analyte[[study]]$gdata[
            ,
            marker
          ]
        }
      ),
      use.names = FALSE
    )

  genotype_expected =
    as.numeric(
      scale(
        genotype_raw
      )
    )

  if (
    length(study_expected) !=
    nrow(model_frame)
  ) {
    stop(
      "Row-count mismatch for ",
      gene,
      " / ",
      marker,
      ": reconstructed ",
      length(study_expected),
      ", model ",
      nrow(model_frame)
    )
  }

  study_match =
    all(
      as.character(
        model_frame$study
      ) ==
      study_expected
    )

  genotype_difference =
    max(
      abs(
        as.numeric(
          model_frame$geno
        ) -
        genotype_expected
      ),
      na.rm = TRUE
    )

  if (!study_match) {
    stop(
      "Study-order mismatch for ",
      gene,
      " / ",
      marker
    )
  }

  if (
    !is.finite(
      genotype_difference
    ) ||
    genotype_difference >
    1e-8
  ) {
    stop(
      "Genotype-order mismatch for ",
      gene,
      " / ",
      marker,
      "; maximum difference = ",
      genotype_difference
    )
  }

  list(
    support = support,
    study =
      study_expected,
    clone =
      clone_expected,
    genotype_difference =
      genotype_difference
  )
}


########################################
# CR1 CLUSTER-ROBUST COVARIANCE        #
########################################

cluster_vcov_cr1 = function(
  fit,
  cluster
) {

  X =
    stats::model.matrix(
      fit
    )

  residual =
    stats::residuals(
      fit
    )

  cluster =
    factor(
      cluster
    )

  N =
    nrow(X)

  K =
    ncol(X)

  G =
    nlevels(
      cluster
    )

  if (G <= K) {
    stop(
      "Too few clusters for CR1 covariance."
    )
  }

  bread =
    solve(
      crossprod(X)
    )

  score =
    X *
    residual

  score_cluster =
    rowsum(
      score,
      group = cluster,
      reorder = FALSE
    )

  meat =
    crossprod(
      score_cluster
    )

  correction =
    (
      G /
      (G - 1)
    ) *
    (
      (N - 1) /
      (N - K)
    )

  vcov =
    correction *
    bread %*%
    meat %*%
    bread

  list(
    vcov = vcov,
    n_clusters = G,
    df = G - 1
  )
}


########################################
# EXTRACT LOCUS-SPECIFIC SENSITIVITY   #
########################################

extract_clustered = function(x) {

  gene =
    as.character(
      x$gene
    )

  marker =
    as.character(
      x$marker
    )

  key =
    paste(
      gene,
      marker,
      sep = "::"
    )

  m =
    x$interaction

  if (
    is.null(m) ||
    !inherits(
      m,
      "lmerMod"
    )
  ) {
    stop(
      "Missing interaction model for ",
      key
    )
  }

  mf =
    stats::model.frame(
      m
    )

  alignment =
    reconstruct_rows(
      gene = gene,
      marker = marker,
      model_frame = mf
    )

  lm_fit =
    stats::lm(
      pheno ~
        analyte +
        geno +
        geno:analyte,
      data = mf
    )

  lm_beta =
    stats::coef(
      lm_fit
    )

  lmer_beta =
    lme4::fixef(
      m
    )

  ixn_term =
    grep(
      "(^geno:analyte|^analyte.*:geno$)",
      names(
        lm_beta
      ),
      value = TRUE
    )

  if (
    length(ixn_term) != 1L
  ) {
    stop(
      "Could not identify one interaction term for ",
      key
    )
  }

  ixn_term =
    ixn_term[[1]]

  lmer_ixn_term =
    grep(
      "(^geno:analyte|^analyte.*:geno$)",
      names(
        lmer_beta
      ),
      value = TRUE
    )

  if (
    length(lmer_ixn_term) != 1L
  ) {
    stop(
      "Could not identify lmer interaction term for ",
      key
    )
  }

  lmer_ixn_term =
    lmer_ixn_term[[1]]

  beta_lm =
    unname(
      lm_beta[
        ixn_term
      ]
    )

  beta_lmer =
    unname(
      lmer_beta[
        lmer_ixn_term
      ]
    )

  ordinary_lm_se =
    summary(
      lm_fit
    )$coefficients[
      ixn_term,
      "Std. Error"
    ]

  lmer_vcov =
    as.matrix(
      stats::vcov(
        m
      )
    )

  lmer_se =
    sqrt(
      lmer_vcov[
        lmer_ixn_term,
        lmer_ixn_term
      ]
    )

  clustered =
    cluster_vcov_cr1(
      fit = lm_fit,
      cluster =
        alignment$clone
    )

  clustered_se =
    sqrt(
      clustered$vcov[
        ixn_term,
        ixn_term
      ]
    )

  clustered_t =
    beta_lm /
    clustered_se

  clustered_p =
    2 *
    stats::pt(
      -abs(
        clustered_t
      ),
      df =
        clustered$df
    )

  random_effect =
    as.data.frame(
      lme4::VarCorr(
        m
      )
    )

  study_variance =
    random_effect$vcov[
      random_effect$grp ==
      "study"
    ]

  if (
    length(study_variance) == 0L
  ) {
    study_variance =
      NA_real_
  } else {
    study_variance =
      study_variance[[1]]
  }

  data.frame(
    gene = gene,
    marker = marker,
    in_joint =
      key %in%
      joint_keys,
    in_strict =
      key %in%
      strict_keys,
    n_observations =
      nrow(mf),
    n_studies =
      length(
        alignment$support
      ),
    n_clone_clusters =
      clustered$n_clusters,
    cluster_df =
      clustered$df,
    maximum_genotype_alignment_difference =
      alignment$genotype_difference,
    study_random_variance =
      study_variance,
    interaction_effect_lmer =
      beta_lmer,
    interaction_effect_lm =
      beta_lm,
    effect_difference =
      beta_lm -
      beta_lmer,
    lmer_se =
      lmer_se,
    ordinary_lm_se =
      ordinary_lm_se,
    clustered_se =
      clustered_se,
    clustered_se_ratio =
      clustered_se /
      lmer_se,
    clustered_t =
      clustered_t,
    clustered_p =
      clustered_p,
    original_targeted_lrt_p =
      as.numeric(
        x$interaction_pval
      ),
    stringsAsFactors = FALSE
  )
}


message(
  "Calculating clone-clustered covariance estimates..."
)

clustered_results =
  do.call(
    rbind,
    lapply(
      modelQTL,
      extract_clustered
    )
  )


########################################
# MERGE STORED POST-HOC RESULTS        #
########################################

joint_stored =
  qtls_aug_analyte_geno[
    ,
    c(
      "gene",
      "marker",
      "interaction_effect",
      "interaction_pval"
    ),
    drop = FALSE
  ]

names(
  joint_stored
)[
  names(
    joint_stored
  ) ==
  "interaction_effect"
] =
  "stored_joint_interaction_effect"

names(
  joint_stored
)[
  names(
    joint_stored
  ) ==
  "interaction_pval"
] =
  "stored_joint_interaction_p"

clustered_results =
  merge(
    clustered_results,
    joint_stored,
    by = c(
      "gene",
      "marker"
    ),
    all.x = TRUE,
    sort = FALSE
  )


########################################
# SUMMARY OF ALIGNMENT AND SE INFLATION#
########################################

joint_results =
  clustered_results[
    clustered_results$in_joint,
    ,
    drop = FALSE
  ]

strict_results =
  clustered_results[
    clustered_results$in_strict,
    ,
    drop = FALSE
  ]

dependence_summary =
  data.frame(
    set =
      c(
        "868_joint_QTLs",
        "64_strict_interaction_QTLs"
      ),
    n_loci =
      c(
        nrow(
          joint_results
        ),
        nrow(
          strict_results
        )
      ),
    minimum_clone_clusters =
      c(
        min(
          joint_results$n_clone_clusters
        ),
        min(
          strict_results$n_clone_clusters
        )
      ),
    median_clone_clusters =
      c(
        stats::median(
          joint_results$n_clone_clusters
        ),
        stats::median(
          strict_results$n_clone_clusters
        )
      ),
    maximum_clone_clusters =
      c(
        max(
          joint_results$n_clone_clusters
        ),
        max(
          strict_results$n_clone_clusters
        )
      ),
    median_lmer_se =
      c(
        stats::median(
          joint_results$lmer_se
        ),
        stats::median(
          strict_results$lmer_se
        )
      ),
    median_clustered_se =
      c(
        stats::median(
          joint_results$clustered_se
        ),
        stats::median(
          strict_results$clustered_se
        )
      ),
    median_se_ratio =
      c(
        stats::median(
          joint_results$clustered_se_ratio
        ),
        stats::median(
          strict_results$clustered_se_ratio
        )
      ),
    q25_se_ratio =
      c(
        as.numeric(
          stats::quantile(
            joint_results$clustered_se_ratio,
            0.25
          )
        ),
        as.numeric(
          stats::quantile(
            strict_results$clustered_se_ratio,
            0.25
          )
        )
      ),
    q75_se_ratio =
      c(
        as.numeric(
          stats::quantile(
            joint_results$clustered_se_ratio,
            0.75
          )
        ),
        as.numeric(
          stats::quantile(
            strict_results$clustered_se_ratio,
            0.75
          )
        )
      ),
    nominal_original_p_lt_0.05 =
      c(
        sum(
          joint_results$original_targeted_lrt_p <
          0.05
        ),
        sum(
          strict_results$original_targeted_lrt_p <
          0.05
        )
      ),
    nominal_clustered_p_lt_0.05 =
      c(
        sum(
          joint_results$clustered_p <
          0.05
        ),
        sum(
          strict_results$clustered_p <
          0.05
        )
      ),
    maximum_abs_effect_difference =
      c(
        max(
          abs(
            joint_results$effect_difference
          )
        ),
        max(
          abs(
            strict_results$effect_difference
          )
        )
      ),
    maximum_genotype_alignment_difference =
      c(
        max(
          joint_results$maximum_genotype_alignment_difference
        ),
        max(
          strict_results$maximum_genotype_alignment_difference
        )
      ),
    maximum_study_random_variance =
      c(
        max(
          joint_results$study_random_variance,
          na.rm = TRUE
        ),
        max(
          strict_results$study_random_variance,
          na.rm = TRUE
        )
      ),
    stringsAsFactors = FALSE
  )


########################################
# POWER WITH CLONE-CLUSTERED SES       #
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

power_rows =
  list()

power_index =
  1L

for (
  threshold_index in
  seq_len(
    nrow(
      thresholds
    )
  )
) {

  alpha =
    thresholds$alpha[
      threshold_index
    ]

  for (
    delta in
    delta_grid
  ) {

    ####################################
    # ORIGINAL SE APPROXIMATION        #
    ####################################

    critical_z =
      stats::qnorm(
        1 -
        alpha /
        2
      )

    ncp_original =
      delta /
      joint_results$lmer_se

    original_power =
      stats::pnorm(
        -critical_z,
        mean =
          ncp_original
      ) +
      stats::pnorm(
        critical_z,
        mean =
          ncp_original,
        lower.tail = FALSE
      )

    ####################################
    # CLONE-CLUSTERED APPROXIMATION    #
    ####################################

    critical_t =
      stats::qt(
        1 -
        alpha /
        2,
      df =
        joint_results$cluster_df
      )

    ncp_cluster =
      delta /
      joint_results$clustered_se

    clustered_power =
      stats::pt(
        -critical_t,
        df =
          joint_results$cluster_df,
        ncp =
          ncp_cluster
      ) +
      stats::pt(
        critical_t,
        df =
          joint_results$cluster_df,
        ncp =
          ncp_cluster,
        lower.tail = FALSE
      )

    power_rows[[
      power_index
    ]] =
      data.frame(
        threshold =
          thresholds$threshold[
            threshold_index
          ],
        alpha =
          alpha,
        true_abs_interaction =
          delta,
        inference =
          c(
            "original",
            "clone_clustered"
          ),
        mean_power =
          c(
            mean(
              original_power
            ),
            mean(
              clustered_power
            )
          ),
        median_power =
          c(
            stats::median(
              original_power
            ),
            stats::median(
              clustered_power
            )
          ),
        p10_power =
          c(
            as.numeric(
              stats::quantile(
                original_power,
                0.10
              )
            ),
            as.numeric(
              stats::quantile(
                clustered_power,
                0.10
              )
            )
          ),
        p90_power =
          c(
            as.numeric(
              stats::quantile(
                original_power,
                0.90
              )
            ),
            as.numeric(
              stats::quantile(
                clustered_power,
                0.90
              )
            )
          ),
        stringsAsFactors = FALSE
      )

    power_index =
      power_index +
      1L
  }
}

clustered_power =
  do.call(
    rbind,
    power_rows
  )


########################################
# PREVALENCE SENSITIVITY               #
########################################

prevalence =
  c(
    0.05,
    0.10,
    0.25,
    0.50,
    1.00
  )

clustered_prevalence =
  do.call(
    rbind,
    lapply(
      seq_len(
        nrow(
          clustered_power
        )
      ),
      function(i) {

        row =
          clustered_power[
            i,
            ,
            drop = FALSE
          ]

        data.frame(
          threshold =
            row$threshold,
          alpha =
            row$alpha,
          true_abs_interaction =
            row$true_abs_interaction,
          inference =
            row$inference,
          prevalence =
            prevalence,
          n_joint_loci =
            nrow(
              joint_results
            ),
          expected_true_discordant =
            prevalence *
            nrow(
              joint_results
            ),
          mean_power =
            row$mean_power,
          expected_detected =
            prevalence *
            nrow(
              joint_results
            ) *
            row$mean_power,
          stringsAsFactors = FALSE
        )
      }
    )
  )


########################################
# STRICT-QTL TARGETED SENSITIVITY      #
########################################

strict_sensitivity =
  data.frame(
    n_strict =
      nrow(
        strict_results
      ),
    clustered_nominal_p_lt_0.05 =
      sum(
        strict_results$clustered_p <
        0.05
      ),
    clustered_p_below_observed_5pct_FDR_raw_p_boundary =
      sum(
        strict_results$clustered_p <=
        strict_raw_p_boundary
      ),
    observed_5pct_FDR_raw_p_boundary =
      strict_raw_p_boundary,
    median_clustered_p =
      stats::median(
        strict_results$clustered_p
      ),
    stringsAsFactors = FALSE
  )


########################################
# SAVE OUTPUTS                         #
########################################

write.table(
  clustered_results,
  file =
    paste0(
      dependence_asys,
      "clone_clustered_locus_results.tsv"
    ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  dependence_summary,
  file =
    paste0(
      dependence_asys,
      "clone_clustered_summary.tsv"
    ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  clustered_power,
  file =
    paste0(
      dependence_asys,
      "clone_clustered_power.tsv"
    ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  clustered_prevalence,
  file =
    paste0(
      dependence_asys,
      "clone_clustered_prevalence.tsv"
    ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  strict_sensitivity,
  file =
    paste0(
      dependence_asys,
      "strict_qtl_clone_sensitivity.tsv"
    ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  list(
    locus_results =
      clustered_results,
    summary =
      dependence_summary,
    power =
      clustered_power,
    prevalence =
      clustered_prevalence,
    strict_sensitivity =
      strict_sensitivity
  ),
  file =
    paste0(
      dependence_imed,
      "review_interaction_dependence.rds"
    )
)

message("")
message(
  "Clone-dependence sensitivity complete."
)
message(
  "Joint loci: ",
  nrow(
    joint_results
  )
)
message(
  "Strict loci: ",
  nrow(
    strict_results
  )
)
message(
  "Median clustered/original SE ratio, joint: ",
  signif(
    stats::median(
      joint_results$clustered_se_ratio
    ),
    4
  )
)
