# SOURCE FXNS AND CONSTANTS
source("/PHShome/kv1026/p05/remap_functions.R")  # imed is defined inside

# Get the studies object
studies = lapply(classicMetadata$study, function(x) {
  load(paste0(imed, x, ".rda"))
  get(x)
})
names(studies) = classicMetadata$study

# Cores: prefer SLURM_CPUS_PER_TASK, fall back to detectCores()
num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                  parallel::detectCores()))
if (is.na(num_cores) || num_cores < 1) num_cores = 1

# Compute gene_marker_associations (LODs and betas) in parallel
gene_marker_associations = parallel::mclapply(
  studies,
  function(x) fasterLOD(
    n.pheno = nrow(x$pheno),
    pheno.s = scale(x$pheno),
    gdata.s = scale(x$gdata),
    betas   = TRUE
  ),
  mc.cores = num_cores
)
save(gene_marker_associations, file = paste0(imed, "gene_marker_associations.rda"))
rm(gene_marker_associations)

# Compute gene_marker_pvals in parallel
gene_marker_pvals = parallel::mclapply(
  studies,
  function(x) fasterP(
    n.pheno = nrow(x$pheno),
    pheno.s = scale(x$pheno),
    gdata.s = scale(x$gdata)
  ),
  mc.cores = num_cores
)
save(gene_marker_pvals, file = paste0(imed, "gene_marker_pvals.rda"))
rm(gene_marker_pvals)

# # SANITY CHECK

# cor_checks = lapply(names(gene_marker_pvals), function(study) {
#   rmat = as.matrix(gene_marker_associations[[study]]$r)
#   pmat = as.matrix(gene_marker_pvals[[study]])

#   # total elements
#   n_total = length(rmat)
#   n_samp  = min(100, n_total)  # sample up to 100

#   idx = sample.int(n_total, n_samp)
#   rvals = abs(rmat[idx])
#   pvals = pmat[idx]

#   cor_val = cor(rvals, -log10(pvals))
#   data.frame(study = study, cor = cor_val)
# })

# cor_df = do.call(rbind, cor_checks)
# print(cor_df)