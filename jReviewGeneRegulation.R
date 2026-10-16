########################################
## FRANK REVIEW CIS TRANS COMPARISON   ##
########################################

# Usage: Rscript --vanilla jReviewGeneRegulation.R P05_ROOT AB18_QTLS_TSV TESTED_GENES_TXT NEW_OUTPUT
# AB18_QTLS_TSV must contain one row per published significant AB18 eQTL,
# with unique qtl_id, gene (systematic ORF) and class (cis or trans).
# Preserve the published local/distant classification in the reviewed adapter.
# TESTED_GENES_TXT must include every gene tested, including genes with no QTL.
# Do not construct that universe from the significant-QTL list.
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Supply root, AB18 QTL table, tested genes and new output.")
root = normalizePath(args[[1]], mustWork = TRUE)
out = args[[4]]
if (file.exists(out)) stop("Output already exists; refusing to overwrite.")
if (!requireNamespace("readxl", quietly = TRUE)) stop("Package readxl required.")
q = read.delim(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("qtl_id", "gene", "class") %in% names(q)),
          nrow(q) > 0L, !anyNA(q[, c("qtl_id", "gene", "class")]),
          !anyDuplicated(q$qtl_id), all(nzchar(q$qtl_id)),
          all(nzchar(q$gene)), all(q$class %in% c("cis", "trans")))
tested = readLines(args[[3]], warn = FALSE)
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
