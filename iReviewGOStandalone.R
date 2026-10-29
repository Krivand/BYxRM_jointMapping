# Rebuild only GO results in a new output directory. No power/mapping reruns.
args = commandArgs(TRUE)
if (length(args) != 1L) stop("Supply one new GO output directory.")
out = args[[1]]
if (file.exists(out) || dir.exists(out)) stop("Output already exists: ", out)
source_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(source_arg) == 1L)
script_dir = dirname(normalizePath(sub("^--file=", "", source_arg)))
dir.create(out, recursive = TRUE)
out = normalizePath(out)
Sys.setenv(REVIEW_GO_OUTPUT_DIR = file.path(out, "tables"))

go_env = new.env(parent = globalenv())
sys.source(file.path(script_dir, "jReviewGO.R"), envir = go_env)

# Reuse only the writer definitions, never execute the full finalizer.
writer = new.env(parent = globalenv())
writer$build_asys = paste0(out, "/")
needed = c("clean_name", "clean_visible_terms", "sheet_widths", "write_plain_workbook")
found = character(0)
for (expr in parse(file.path(script_dir, "iReviewWorkbookFinalize.R"))) {
  if (is.call(expr) && as.character(expr[[1]]) %in% c("=", "<-") &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% needed) {
    eval(expr, envir = writer)
    found = c(found, as.character(expr[[2]]))
  }
}
stopifnot(setequal(found, needed))

# Boundary checks for empty/full term membership and separated classes.
compare = go_env$compare_go_term
stopifnot(compare(character(0), c("a", "b"), c("c", "d"))$fisher_p == 1,
          !compare(c("a", "b", "c", "d"), c("a", "b"), c("c", "d"))$estimable,
          is.infinite(compare(c("a", "b"), c("a", "b"), c("c", "d"))$odds_ratio))
comparison = go_env$class_comparison
for (ontology in go_env$ontologies) {
  d = comparison[comparison$ontology == ontology, , drop = FALSE]
  primary = go_env$all_go_results
  primary = primary[primary$gene_set == "discordant" & primary$ontology == ontology, ]
  stopifnot(setequal(d$GO.ID, primary$GO.ID),
            all(d$discordant_members + d$discordant_nonmembers == 37L),
            all(d$other_members + d$other_nonmembers == 27L),
            isTRUE(all.equal(d$fisher_q, p.adjust(d$fisher_p, "BH"))))
  # Independent hypergeometric upper-tail calculation of one-sided Fisher P.
  expected_p = phyper(d$discordant_members - 1,
                     d$discordant_members + d$other_members,
                     d$discordant_nonmembers + d$other_nonmembers,
                     37, lower.tail = FALSE)
  stopifnot(max(abs(d$fisher_p - expected_p)) < 1e-12)
}

all_go = go_env$all_go_results
sheets = list(
  "Summary" = go_env$discordant_summary,
  "All 64 interactions" = all_go[all_go$gene_set == "all_interaction", ],
  "Discordant 37" = all_go[all_go$gene_set == "discordant", ],
  "Protein specific 14" = all_go[all_go$gene_set == "protein_specific", ],
  "mRNA specific 9" = all_go[all_go$gene_set == "mRNA_specific", ],
  "Same direction 4" = all_go[all_go$gene_set == "same_direction", ],
  "Significant terms" = go_env$significant_go,
  "Discordant vs other" = comparison,
  "Gene sets" = data.frame(
    gene_set = rep(names(go_env$gene_sets), lengths(go_env$gene_sets)),
    gene = unlist(go_env$gene_sets, use.names = FALSE)
  ),
  "MKT1 BUL2 summary" = go_env$hotspot_summary
)
writer$write_plain_workbook(sheets, "interaction_GO.xlsx")
path = file.path(out, "interaction_GO.xlsx")
stopifnot(identical(openxlsx::getSheetNames(path), names(sheets)))
for (name in names(sheets)) {
  observed = openxlsx::read.xlsx(path, sheet = name, check.names = FALSE,
                                 skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  expected = writer$clean_visible_terms(sheets[[name]])
  stopifnot(identical(names(observed), names(expected)), nrow(observed) == nrow(expected))
  for (j in seq_along(expected)) {
    a = as.character(observed[[j]])
    b = as.character(expected[[j]])
    a[is.na(a) | a == ""] = "NA"
    b[is.na(b) | b == ""] = "NA"
    if (is.numeric(expected[[j]])) {
      stopifnot(isTRUE(all.equal(suppressWarnings(as.numeric(a)),
                                suppressWarnings(as.numeric(b)), tolerance = 1e-12)))
      infinite = is.infinite(expected[[j]])
      stopifnot(identical(a[infinite], b[infinite]))
    } else stopifnot(identical(a, b))
  }
}
write.table(all_go, file.path(out, "tables", "GO_all_terms.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(annotation = go_env$go_annotation, gene_sets = go_env$gene_sets,
             all_go_results = all_go, class_comparison = comparison),
        file.path(out, "GO_results.rds"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
report = do.call(rbind, lapply(go_env$ontologies, function(ontology) {
  d = comparison[comparison$ontology == ontology, ]
  summary = go_env$discordant_summary
  summary = summary[summary$ontology == ontology, ]
  data.frame(ontology = ontology, tested = nrow(d),
             first_pass_discordant_hits = nrow(summary),
             followup_hits = sum(d$fisher_q < 0.05),
             significant_both = sum(summary$significant_both))
}))
write.table(report, file.path(out, "GO_run_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(report, row.names = FALSE)
cat("GO-only workbook and statistical checks passed. Output:", out, "\n")
