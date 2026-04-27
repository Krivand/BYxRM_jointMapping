#####################################
#####################################
## JOINT DATASET REMAPPING ANALYTE ##
#####################################
#####################################

##############################
## PULL AND PARSE ARGUMENTS ##
##############################

# source functions and constants
source(paste0("/PHShome/kv1026/p05/", "remap_functions.R"))

# args
args = commandArgs(trailingOnly = TRUE)
chunk_to_run = as.integer(args[[1]])
chunk_size = if (length(args) >= 2) as.integer(args[[2]]) else 50L
fx_arg = if (length(args) >= 3) args[[3]] else ""

# parse fixed effects
fx_arg = trimws(fx_arg)
fixedEffect = if (nzchar(fx_arg)) fx_arg else NULL
mode_tag = if (is.null(fixedEffect)) "geno" else fixedEffect

# designate number of cores
cores = suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
if (is.na(cores) || cores < 1L) cores = 8L

#########################
# BUILD MAPPING OBJECTS #
#########################

studies = lapply(scaMetadata$study, function(x) {
  load(paste0(imed, x, ".rda"))
  get(x)
}) %>% `names<-`(scaMetadata$study)

scaMetadata$genes = lapply(scaMetadata$study,
  function(x) colnames(studies[[x]]$pheno))
scaMetadata$segs  = lapply(scaMetadata$study,
  function(x) rownames(studies[[x]]$pheno))

doubleSupp = suppTwo(scaMetadata)
message("Double-supported genes: ", length(doubleSupp))

# union background QTLs
union_bg_file = paste0(imed, "unionAnalyteBG.rda")

if (!file.exists(union_bg_file)) {

  # get the background QTLs
  bgAnalyte =
    lapply(scaMetadata$study, function(x) {
      load(paste0(imed, x, "_Background.rda"))
      get(paste0(x, "_Background"))
    }) %>% `names<-`(scaMetadata$study)

  # collapse and save into the union
  unionAnalyteBG = collapseBgQTLs(bgAnalyte)
  save(unionAnalyteBG, file = union_bg_file)

} else {
  message("Loading existing unionAnalyteBG: ", union_bg_file)
  load(union_bg_file)
}

################
# Output setup #
################

# setup outputs
base_name = paste0("jBH_analyte_", mode_tag)  # geno or fixef
out_dir   = paste0(imed, base_name)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# chunking
gene_chunks = split(doubleSupp, ceiling(seq_along(doubleSupp) / chunk_size))
n_chunks    = length(gene_chunks)

message("Mode: ", mode_tag)
message("Total genes: ", length(doubleSupp),
        " | chunk size: ", chunk_size,
        " | number of chunks: ", n_chunks)

if (is.na(chunk_to_run) || chunk_to_run < 1L || chunk_to_run > n_chunks)
  stop("chunk_index out of range: 1..", n_chunks)

chunk_genes = gene_chunks[[chunk_to_run]]
chunk_file = paste0(out_dir, "/chunk", chunk_to_run, ".rda")

if (file.exists(chunk_file)) {
  message("Chunk ", chunk_to_run, " already exists at ", chunk_file, " — exiting.")
  quit(save = "no", status = 0)
}

###################
# Run the mapping #
###################

message("Running jointMarkerFDR mode=", mode_tag,
        " chunk ", chunk_to_run, "/", n_chunks,
        " (", length(chunk_genes), " genes) with ", cores, " threads.")

chunkAssocs = jointMarkerFDR(
  studies       = studies,
  targetMarkers = colnames(studies$Brem_2002$gdata),
  unionBg       = unionAnalyteBG,
  smetadata     = scaMetadata,
  suppGenes     = chunk_genes,
  FDR           = 0.05,
  fixedEffect   = fixedEffect,
  fn            = base_name,
  correctBG     = TRUE,
  correctGRM    = FALSE,
  cores         = cores,
  chromosomes   = NULL
)

save(chunkAssocs, file = chunk_file)
message("Saved chunk ", chunk_to_run, " to ", chunk_file)
