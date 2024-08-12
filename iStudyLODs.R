# SOURCE FXNS AND CONSTANTS
source("/home/albertf/vandy147/p05/remap_functions.R")

# Get the studies object
studies = lapply(classicMetadata$study, function(x) {
  load(paste0(imed, x, '.rda'))
  get(x)
})
names(studies) = classicMetadata$study

# Detect the number of cores
num_cores = detectCores()

# Create gene_marker_associations_object and save it (LODS and Rs)
gene_marker_associations = mclapply(
  studies,function(x)fasterLOD(n.pheno = nrow(x$pheno),
                               pheno.s = scale(x$pheno),
                               gdata.s = scale(x$gdata),
                               betas = TRUE), mc.cores = num_cores)

save(gene_marker_associations,file=paste0(imed,'gene_marker_associations.rda'))