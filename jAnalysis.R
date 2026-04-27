#################
# preliminaries #
#################

# source functions and constants
source("/PHShome/kv1026/p05/remap_functions.R")[]

# rename for symmetry
metadata_analyte = scaMetadata
metadata_medium  = medMeta

# designate number of cores
cores = suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
if (is.na(cores) || cores < 1L) cores = 8L

#####################################
# load previously generated objects #
#####################################

chunk_in = c(
  "/PHShome/kv1026/p05/imed/jBH_analyte_geno",
  "/PHShome/kv1026/p05/imed/jBH_analyte_analyte",
  "/PHShome/kv1026/p05/imed/jBH_medium_geno",
  "/PHShome/kv1026/p05/imed/jBH_medium_medium"
)

merge_out = c(
  "/PHShome/kv1026/p05/imed/jBH_analyte_geno.rda",
  "/PHShome/kv1026/p05/imed/jBH_analyte_ixn.rda",
  "/PHShome/kv1026/p05/imed/jBH_medium_geno.rda",
  "/PHShome/kv1026/p05/imed/jBH_medium_ixn.rda"
)

for (i in seq_along(chunk_in)) {

  if (!file.exists(merge_out[[i]])) {
    aggregate_chunks(chunk_in[[i]], merge_out[[i]])
  }
  
  load(merge_out[[i]])

}

#########################
# build mapping objects #
#########################

# ANALYTE-SPECIFIC OBJECTS

studies_analyte = lapply(metadata_analyte$study, function(x) {
  load(paste0(imed, x, ".rda"))
  get(x)
}) %>% `names<-`(metadata_analyte$study)

metadata_analyte$genes = lapply(metadata_analyte$study,
  function(x) colnames(studies_analyte[[x]]$pheno))
metadata_analyte$segs = lapply(metadata_analyte$study,
  function(x) rownames(studies_analyte[[x]]$pheno))

supported_analyte = suppTwo(metadata_analyte)
writeLines(supported_analyte, paste0(imed, "supported_analyte.txt"))
message("Analyte-supported genes: ", length(supported_analyte))

analyte_bg_file = paste0(imed, "unionAnalyteBG.rda")
if (!file.exists(analyte_bg_file)) {

  bgAnalyte = lapply(metadata_analyte$study, function(x) {
    load(paste0(imed, x, "_Background.rda"))
    get(paste0(x, "_Background"))
  }) %>% `names<-`(metadata_analyte$study)

  unionBG_analyte = collapseBgQTLs(bgAnalyte)
  save(unionBG_analyte, file = analyte_bg_file)

} else {
  message("Loading existing unionBG_analyte: ", analyte_bg_file)
  load(analyte_bg_file) # expects unionAnalyteBG
  unionBG_analyte = unionAnalyteBG
  rm(unionAnalyteBG)
}

# MEDIUM-SPECIFIC OBJECTS

studies_medium = lapply(metadata_medium$study, function(x) {
  load(paste0(imed, x, ".rda"))
  get(x)
}) %>% `names<-`(metadata_medium$study)

metadata_medium$genes = lapply(metadata_medium$study,
  function(x) colnames(studies_medium[[x]]$pheno))
metadata_medium$segs = lapply(metadata_medium$study,
  function(x) rownames(studies_medium[[x]]$pheno))

# medium-only support (pre-intersection)
proto_medium = suppTwo(metadata_medium, type = "medium")
message("Medium-supported genes (pre-intersection): ", length(proto_medium))

medium_bg_file = paste0(imed, "unionMediumBG.rda")
if (!file.exists(medium_bg_file)) {

  bgMedium = lapply(metadata_medium$study, function(x) {
    load(paste0(imed, x, "_Background.rda"))
    get(paste0(x, "_Background"))
  }) %>% `names<-`(metadata_medium$study)

  unionBG_medium = collapseBgQTLs(bgMedium)
  save(unionBG_medium, file = medium_bg_file)

} else {
  message("Loading existing unionBG_medium: ", medium_bg_file)
  load(medium_bg_file) # expects unionMediumBG
  unionBG_medium = unionMediumBG
  rm(unionMediumBG)
}

# intersect to create final support set for medium
supported_medium = intersect(supported_analyte, proto_medium)
writeLines(supported_medium,  paste0(imed, "supported_medium.txt"))
message("Genes supported by analyte and medium data: ", length(supported_medium))

#########################
# summarize QTLs nicely #
#########################

levels = c("analyte", "medium")
modes  = c("geno", "ixn")

for (lvl in levels) {
  for (md in modes) {

    # NEW output (augmented dfQTL, with modelQTL attached as attr)
    out_name = paste0("qtls_aug_", lvl, "_", md)
    out_file = paste0(imed, out_name, ".rda")

    if (!file.exists(out_file)) {

      # fixedEffect must be a COLUMN NAME in smetadata
      fixedEffect = if (lvl == "medium") "medium" else "analyte"

      # levels used for binning
      fixed_levels = if (lvl == "medium") c("SC", "YNB2G") else c("mRNA", "protein")

      df_in = get(paste0("jBH_", lvl, "_", md))

      qtls_aug = attach_drill_info_and_save(
        dfQTL = df_in,
        studies = get(paste0("studies_", lvl)),
        unionBg = get(paste0("unionBG_", lvl)),
        smetadata = get(paste0("metadata_", lvl)),
        suppGenes = get(paste0("supported_", lvl)),
        fixedEffect = fixedEffect,
        correctBG = TRUE,
        correctGRM = FALSE,
        cores = cores,
        chromosomes = NULL,
        do_bin = TRUE,
        bin_level1 = fixed_levels[[1]],
        bin_level2 = fixed_levels[[2]],
        bin_alpha = 0.05,
        out_name = out_name,
        out_file = out_file
      )

    } else {

      load(out_file)  # loads object named out_name into .GlobalEnv
      qtls_aug = get(out_name, envir = .GlobalEnv)

    }

    message("OK: ", out_name)
  }
}
