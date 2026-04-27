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

# ANALYTE PROCESSING

# determine supported genes
studiesAnalyte = lapply(scaMetadata$study, function(x) {
  load(paste0(imed, x, ".rda"))
  get(x)
}) %>% `names<-`(scaMetadata$study)

scaMetadata$genes = lapply(scaMetadata$study,
  function(x) colnames(studiesAnalyte[[x]]$pheno))
scaMetadata$segs  = lapply(scaMetadata$study,
  function(x) rownames(studiesAnalyte[[x]]$pheno))

suppAnalyte = suppTwo(scaMetadata)
message("Analyte-supported genes: ", length(suppAnalyte))

# unite background QTLs
analyte_bg_file = paste0(imed, "unionAnalyteBG.rda")

if (!file.exists(analyte_bg_file)) {

  # get the background QTLs
  bgAnalyte =
    lapply(scaMetadata$study, function(x) {
      load(paste0(imed, x, "_Background.rda"))
      get(paste0(x, "_Background"))
    }) %>% `names<-`(scaMetadata$study)

  # collapse and save into the union
  unionAnalyteBG = collapseBgQTLs(bgAnalyte)
  save(unionAnalyteBG, file = analyte_bg_file)

} else {
  message("Loading existing unionAnalyteBG: ", analyte_bg_file)
  load(analyte_bg_file)
}

# MEDIUM PROCESSING

# determine supported genes
studiesMedium = lapply(medMeta$study, function(x) {
  load(paste0(imed, x, ".rda"))
  get(x)
}) %>% `names<-`(medMeta$study)

medMeta$genes = lapply(medMeta$study,
  function(x) colnames(studiesMedium[[x]]$pheno))
medMeta$segs  = lapply(medMeta$study,
  function(x) rownames(studiesMedium[[x]]$pheno))

suppMedium = suppTwo(medMeta,type="medium")
message("Medium-supported genes: ", length(suppMedium))

# unite background QTLs
medium_bg_file = paste0(imed, "unionMediumBG.rda")

if (!file.exists(medium_bg_file)) {

  # get the background QTLs
  bgMedium =
    lapply(medMeta$study, function(x) {
      load(paste0(imed, x, "_Background.rda"))
      get(paste0(x, "_Background"))
    }) %>% `names<-`(medMeta$study)

  # collapse and save into the union
  unionMediumBG = collapseBgQTLs(bgMedium)
  save(unionMediumBG, file = medium_bg_file)

} else {
  message("Loading existing unionMediumBG: ", medium_bg_file)
  load(medium_bg_file)
}

##########################
# Intersect support sets #
##########################

doubleSupp = intersect(suppAnalyte, suppMedium)
message("Genes supported by analyte and medium data: ", length(doubleSupp))

##################################
# downsample Albert for analysis #
##################################

study_name = "AlbertAlt_2018"
n_segs_keep = 100L

# ensure we are all in the same order
ok = with(studiesMedium[[study_name]], {
  g = rownames(gdata)

  all(
    sub("-.*", "", rownames(pheno))  == g,
    sub("-.*", "", rownames(pclean)) == g,
    sub("-.*", "", rownames(covars)) == g,
    sub("-.*", "", medMeta$segs[medMeta$study == study_name][[1]]) == g
  )
})

cat(
  "Are Albert segregants all correctly in the same order? ",
  if (ok) "YES" else "NO",
  "\n"
)

if (!ok) stop("Sample order mismatch: refusing to downsample.", call. = FALSE)

# deterministic subset: stable across chunks/runs (no RNG)
g_ids = rownames(studiesMedium[[study_name]]$gdata)

seg_subset_idx = order(
  vapply(g_ids, function(s) sum(as.integer(charToRaw(s))), numeric(1)),
  g_ids
)[seq_len(n_segs_keep)]

# apply in place
studiesMedium[[study_name]]$gdata = 
  studiesMedium[[study_name]]$gdata[seg_subset_idx, , drop = FALSE]
studiesMedium[[study_name]]$pheno = 
  studiesMedium[[study_name]]$pheno[seg_subset_idx, , drop = FALSE]
studiesMedium[[study_name]]$pclean = 
  studiesMedium[[study_name]]$pclean[seg_subset_idx, , drop = FALSE]
studiesMedium[[study_name]]$covars = 
  studiesMedium[[study_name]]$covars[seg_subset_idx, , drop = FALSE]

meta_row = which(medMeta$study == study_name)
if (length(meta_row) != 1) stop("Expected exactly one medMeta row for ",
  study_name, call. = FALSE)

medMeta$segs[[meta_row]] = medMeta$segs[[meta_row]][seg_subset_idx]

message("Downsampled ", study_name, " to ", n_segs_keep, " segregants.")

# update medMeta metadata to match the downsample
medMeta$sample_size[meta_row] = n_segs_keep

study_name = "AlbertAlt_2018"

cat("\nPost-downsample dimensions for", study_name, ":\n")

cat(
  "  gdata rows : ", nrow(studiesMedium[[study_name]]$gdata),  "\n",
  "  pheno rows : ", nrow(studiesMedium[[study_name]]$pheno),  "\n",
  "  pclean rows: ", nrow(studiesMedium[[study_name]]$pclean), "\n",
  "  covars rows: ", nrow(studiesMedium[[study_name]]$covars), "\n",
  sep = ""
)

meta_row = which(medMeta$study == study_name)
cat(
  "  medMeta segs: ", length(medMeta$segs[[meta_row]]), "\n\n",
  sep = ""
)

################
# Output setup #
################

# setup outputs
base_name = paste0("jBH_medium_", mode_tag)  # geno or fixef
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
  studies       = studiesMedium,
  targetMarkers = colnames(studiesMedium[[1]]$gdata),
  unionBg       = unionMediumBG,
  smetadata     = medMeta,
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
