########################################
## FRANK REVIEW CIS TRANS COMPARISON   ##
########################################

# Usage: Rscript --vanilla jReviewGeneRegulation.R P05_ROOT NEW_OUTPUT
# Uses the existing Albert_2018 QTL table and its cis annotation.
# The tested universe comes from the processed phenotype matrix, not QTL hits.
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply p05 root and new output directory.")
root = normalizePath(args[[1]], mustWork = TRUE)
out = args[[2]]
if (file.exists(out)) stop("Output already exists; refusing to overwrite.")
if (!requireNamespace("readxl", quietly = TRUE)) stop("Package readxl required.")
load_object = function(path, name) {
  e = new.env(parent = emptyenv())
  found = load(path, envir = e)
  if (!(name %in% found)) stop("Missing object: ", name)
  e[[name]]
}
qtls = load_object(file.path(root, "imed/Frank/qtls_list_260204.rda"), "qtls_list")
ab18 = qtls$Albert_2018
stopifnot(is.data.frame(ab18), nrow(ab18) > 0L,
          all(c("gene_ID", "marker_name", "cis") %in% names(ab18)))
cis_value = tolower(as.character(ab18$cis))
if (anyNA(cis_value) || !all(cis_value %in% c("true", "false", "1", "0"))) {
  stop("Unrecognized cis encoding; inspect annotation before proceeding.")
}
q = data.frame(qtl_id = paste(ab18$gene_ID, ab18$marker_name, sep = "::"),
               gene = as.character(ab18$gene_ID),
               class = ifelse(cis_value %in% c("true", "1"), "cis", "trans"))
stopifnot(!anyNA(ab18$gene_ID), !anyNA(ab18$marker_name),
          !anyDuplicated(q$qtl_id), all(nzchar(q$gene)))
study = load_object(file.path(root, "imed/Albert_2018.rda"), "Albert_2018")
tested = colnames(study$pheno)
stopifnot(length(tested) > 0L, !anyDuplicated(tested), all(nzchar(tested)),
          all(q$gene %in% tested))
supported = readLines(file.path(root, "imed/supported_analyte.txt"), warn = FALSE)
stopifnot(length(supported) == 385L, !anyDuplicated(supported))
features = as.data.frame(readxl::read_excel(
  file.path(root, "raw/Albert_2018_gene_characteristics.xlsx"),
  sheet = "SI_Data_06_geneFeatures", na = c("", "NA")
))
stopifnot(all(c("gene", "expression") %in% names(features)),
          !anyDuplicated(features$gene))
expression_genes = features$gene[!is.na(features$expression)]
stopifnot(length(expression_genes) == 5720L,
          all(supported %in% expression_genes))
universe = intersect(expression_genes, tested)
if (!all(supported %in% universe)) {
  stop("Some of the 385 supported genes were not tested in AB18; review denominator.")
}
if (length(setdiff(universe, supported)) == 0L) stop("No comparison genes.")
gene_counts = data.frame(gene = universe, supported = universe %in% supported)
summary_rows = list()
count_rows = list()
for (cls in c("cis", "trans")) {
  counts = table(factor(q$gene[q$class == cls], levels = universe))
  gene_counts[[paste0("n_", cls)]] = as.integer(counts)
  a = as.integer(counts[gene_counts$supported])
  b = as.integer(counts[!gene_counts$supported])
  contingency = matrix(c(sum(a > 0L), sum(a == 0L),
                         sum(b > 0L), sum(b == 0L)), nrow = 2L, byrow = TRUE)
  fit = fisher.test(contingency)
  summary_rows[[cls]] = data.frame(
    class = cls, supported_with_qtl = sum(a > 0L), supported_total = length(a),
    other_with_qtl = sum(b > 0L), other_total = length(b),
    supported_fraction = mean(a > 0L), other_fraction = mean(b > 0L),
    odds_ratio = unname(fit$estimate), ci_lower = fit$conf.int[[1]],
    ci_upper = fit$conf.int[[2]], p = fit$p.value
  )
  count_fit = suppressWarnings(wilcox.test(a, b, exact = FALSE))
  count_rows[[cls]] = data.frame(
    class = cls, supported_mean = mean(a), other_mean = mean(b),
    supported_median = median(a), other_median = median(b),
    p = count_fit$p.value
  )
}
summary_table = do.call(rbind, summary_rows)
summary_table$fdr_two_presence_tests = p.adjust(summary_table$p, method = "BH")
count_table = do.call(rbind, count_rows)
count_table$fdr_two_count_tests = p.adjust(count_table$p, method = "BH")
if (!dir.create(out, recursive = TRUE)) stop("Could not create new output.")
write_tsv = function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE)
}
write_tsv(summary_table, "supported_gene_cis_trans_presence.tsv")
write_tsv(count_table, "supported_gene_cis_trans_counts.tsv")
write_tsv(gene_counts, "supported_gene_regulation_per_gene.tsv")
write_tsv(data.frame(
  expression_universe = length(expression_genes),
  ab18_tested_genes = length(tested), comparison_universe = length(universe),
  supported = length(supported), other = length(setdiff(universe, supported)),
  expression_genes_not_tested = length(setdiff(expression_genes, tested)),
  qtl_rows_outside_comparison_universe = sum(!(q$gene %in% universe))
), "regulation_denominators.tsv")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
print(summary_table, row.names = FALSE)
