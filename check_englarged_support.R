source("/PHShome/kv1026/p05/remap_functions.R")

md =
  studyMetadata[
    studyMetadata$mapping_approach == "QTL" &
      studyMetadata$analyte %in% c("mRNA", "protein"),
    ,
    drop = FALSE
  ]

studies =
  lapply(md$study, function(x) {
    load(paste0(imed, x, ".rda"))
    get(x)
  }) %>%
  `names<-`(md$study)

md$genes =
  lapply(
    md$study,
    function(x) colnames(studies[[x]]$pheno)
  )

mrna_counts =
  table(unlist(md$genes[md$analyte == "mRNA"]))

protein_counts =
  table(unlist(md$genes[md$analyte == "protein"]))

mrna_supported = names(mrna_counts[mrna_counts >= 2])
protein_supported = names(protein_counts[protein_counts >= 2])

doubleSupp_all_mrna = intersect(mrna_supported, protein_supported)

message("Genes in >=2 mRNA datasets: ", length(mrna_supported))
message("Genes in >=2 protein datasets: ", length(protein_supported))
message("Genes in both universes: ", length(doubleSupp_all_mrna))