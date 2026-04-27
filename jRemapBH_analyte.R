#############################
#############################
## JOINT DATASET REMAPPING ##
#############################
#############################

#########################################
# RETRIEVE AND FORMAT NECESSARY OBJECTS #
#########################################

# SOURCE FXNS AND CONSTANTS
source(paste0("/PHShome/kv1026/p05/",
              "remap_functions.R"))

# create an object for easier joint mapping
studies = lapply(scaMetadata$study, function(x) {
  load(paste0(imed, x, ".rda"))
  get(x)
}) %>% `names<-`(scaMetadata$study)

# retrieve background QTLs across studies
bgQTLs = lapply(scaMetadata$study, function(x) {
  load(paste0(imed, x, "_Background.rda"))
  get(paste0(x, "_Background"))
}) %>% `names<-`(scaMetadata$study)

# add gene support to metadata
scaMetadata$genes = lapply(scaMetadata$study, function(x)
  colnames(studies[[x]]$pheno))

# add segregant support to metadata
scaMetadata$segs = lapply(scaMetadata$study, function(x)
  rownames(studies[[x]]$pheno))

# within the synthetic complete media (sc) datasets
# find which genes are supported by at least two datasets
doubleSupp = suppTwo(scaMetadata)
message("Within the SC media with genotype A subset there are ",
        length(doubleSupp),
        " genes supported by two mRNA and two protein datasets")

# take the union of all background QTLs discovered across all datasets
union_bg_file = file.path(imed, "unionBg.rda")

if (!file.exists(union_bg_file)) {
  unionBg = collapseBgQTLs(bgQTLs)
  save(unionBg, file = union_bg_file)
} else {
  message("Union of background QTLs already processed and loading now...")
  load(union_bg_file)
}

################################################################
# JOINT FDR-BASED QTL MAPPING ACROSS ALL GENE-MARKER PAIRS     #
################################################################

## Base folder where everything for this analysis will live
base_name = "jointAssocsAnalyteFDR"
out_dir   = file.path(imed, base_name)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

## Chunking setup
chunk_size = 50L

gene_chunks = split(
  doubleSupp,
  ceiling(seq_along(doubleSupp) / chunk_size)
)

message("Total genes: ", length(doubleSupp),
        " | chunk size: ", chunk_size,
        " | number of chunks: ", length(gene_chunks))

## Specify which chunks to run
chunks_to_run = NULL

if (is.null(chunks_to_run)) {
  chunks_to_run = seq_along(gene_chunks)
}

if (any(chunks_to_run < 1L | chunks_to_run > length(gene_chunks))) {
  stop("chunks_to_run out of range: 1..", length(gene_chunks))
}

for (i in chunks_to_run) {

  chunk_genes = gene_chunks[[i]]
  chunk_file  = file.path(out_dir, paste0("chunk", i, ".rda"))

  if (file.exists(chunk_file)) {
    message("Chunk ", i, " already exists at ", chunk_file, " — skipping.")
    next
  }

  message("Running jointMarkerFDR for chunk ", i, " of ", length(gene_chunks),
          " (", length(chunk_genes), " genes).")

  chunkAssocs = jointMarkerFDR(
    studies       = studies,
    targetMarkers = colnames(studies$Brem_2002$gdata),
    unionBg       = unionBg,
    smetadata     = scaMetadata,
    suppGenes     = chunk_genes,
    FDR           = 0.05,
    fixedEffect   = "analyte",
    fn            = base_name,
    correctBG     = TRUE,
    correctGRM    = FALSE,
    chromosomes   = NULL
  )

  save(chunkAssocs, file = chunk_file)
  message("Saved chunk ", i, " to ", chunk_file)
}

##############################
## MERGE ALL CHUNKS AT END  ##
##############################

chunk_pattern = "^chunk[0-9]+\\.rda$"

chunk_files = list.files(
  path       = out_dir,   # FIXED: search correct directory
  pattern    = chunk_pattern,
  full.names = TRUE
)

# numeric sort (safe for chunk10+)
chunk_files = chunk_files[
  order(as.integer(sub("^chunk([0-9]+)\\.rda$", "\\1", basename(chunk_files))))
]

if (length(chunk_files) == length(gene_chunks)) {

  message("All ", length(gene_chunks), " chunk files found. Merging.")

  jointAssocsAnalyteFDR = NULL

  for (f in chunk_files) {
    message("Loading ", f)
    load(f)  # creates object 'chunkAssocs'

    if (is.null(jointAssocsAnalyteFDR)) {
      jointAssocsAnalyteFDR = chunkAssocs
    } else {
      for (nm in names(chunkAssocs)) {
        jointAssocsAnalyteFDR[[nm]] = rbind(
          jointAssocsAnalyteFDR[[nm]],
          chunkAssocs[[nm]]
        )
      }
    }
  }

  final_file = file.path(out_dir, "jointAssocsAnalyteFDR.rda")
  save(jointAssocsAnalyteFDR, file = final_file)
  message("Saved merged jointAssocsAnalyteFDR to: ", final_file)

} else {
  message("Not all chunks are present (expected ",
          length(gene_chunks),
          ", found ",
          length(chunk_files),
          "). Skipping merge.")
}