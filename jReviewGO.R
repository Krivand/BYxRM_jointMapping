########################################
########################################
## GENETICS REVIEW GO ENRICHMENT      ##
## JOINT ANALYTE INTERACTION GENE SETS ##
########################################
########################################

#################
# PRELIMINARIES #
#################

source("/PHShome/kv1026/p05/remap_functions.R")

review_pkgs = c(
  "AnnotationDbi",
  "org.Sc.sgd.db",
  "topGO"
)

for (pkg in review_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop("Package ", pkg, " is required.")
  }
}

rm(review_pkgs)

go_asys = paste0(sub("/+$", "", Sys.getenv(
  "REVIEW_GO_OUTPUT_DIR", unset = paste0(asys, "review/go/")
)), "/")
ensure_directory(go_asys)

########################
# LOAD CANONICAL INPUT #
########################

supported_analyte = readLines(
  paste0(imed, "supported_analyte.txt")
)

load(
  paste0(imed, "qtls_aug_analyte_ixn.rda")
)

if (
  length(supported_analyte) != 385L ||
  anyDuplicated(supported_analyte)
) {
  stop("Expected 385 unique supported analyte genes.")
}

if (
  nrow(qtls_aug_analyte_ixn) != 64L ||
  anyDuplicated(qtls_aug_analyte_ixn$gene)
) {
  stop("Expected 64 strict interaction QTLs, one per gene.")
}

required_qtl_columns = c(
  "gene",
  "marker",
  "bin",
  "mRNA_effect",
  "protein_effect"
)

missing_qtl_columns = setdiff(
  required_qtl_columns,
  names(qtls_aug_analyte_ixn)
)

if (length(missing_qtl_columns) > 0L) {
  stop(
    "Interaction QTL table is missing columns: ",
    paste(missing_qtl_columns, collapse = ", ")
  )
}

###################
# BUILD GENE SETS #
###################

gene_sets = list(
  all_interaction = sort(
    qtls_aug_analyte_ixn$gene
  ),
  discordant = sort(
    qtls_aug_analyte_ixn$gene[
      qtls_aug_analyte_ixn$bin == "line cross"
    ]
  ),
  protein_specific = sort(
    qtls_aug_analyte_ixn$gene[
      qtls_aug_analyte_ixn$bin == "protein specific"
    ]
  ),
  mRNA_specific = sort(
    qtls_aug_analyte_ixn$gene[
      qtls_aug_analyte_ixn$bin == "mRNA specific"
    ]
  ),
  same_direction = sort(
    qtls_aug_analyte_ixn$gene[
      qtls_aug_analyte_ixn$bin == "same direction"
    ]
  )
)

expected_sizes = c(
  all_interaction = 64L,
  discordant = 37L,
  protein_specific = 14L,
  mRNA_specific = 9L,
  same_direction = 4L
)

observed_sizes = vapply(
  gene_sets,
  length,
  integer(1)
)

if (!identical(observed_sizes, expected_sizes)) {
  stop("Unexpected interaction-class sizes.")
}

########################
# BUILD SGD GO MAPPING #
########################

go_annotation = AnnotationDbi::select(
  org.Sc.sgd.db::org.Sc.sgd.db,
  keys = supported_analyte,
  keytype = "ORF",
  columns = c("GO", "ONTOLOGY")
) %>%
  as.data.frame() %>%
  filter(
    !is.na(ORF),
    !is.na(GO),
    !is.na(ONTOLOGY),
    ORF %in% supported_analyte,
    ONTOLOGY %in% c("BP", "CC", "MF")
  ) %>%
  distinct(ORF, GO, ONTOLOGY)

gene2GO = go_annotation %>%
  dplyr::select(ORF, GO) %>%
  distinct()

gene2GO = split(
  gene2GO$GO,
  gene2GO$ORF
)

if (length(intersect(supported_analyte, names(gene2GO))) != 385L) {
  stop("Not all 385 supported genes have GO annotations.")
}

####################
# HELPER FUNCTIONS #
####################

run_one_go = function(
  interesting_genes,
  set_name,
  ontology
) {

  gene_score = as.integer(
    supported_analyte %in% interesting_genes
  )

  names(gene_score) = supported_analyte

  go_data = new(
    "topGOdata",
    ontology = ontology,
    allGenes = gene_score,
    geneSel = function(x) x == 1L,
    annot = annFUN.gene2GO,
    gene2GO = gene2GO,
    nodeSize = 3L
  )

  go_test = runTest(
    go_data,
    algorithm = "classic",
    statistic = "fisher"
  )

  go_table = GenTable(
    go_data,
    classicFisher = go_test,
    topNodes = length(score(go_test)),
    numChar = 200
  ) %>%
    as.data.frame()

  term_members = genesInTerm(
    go_data,
    whichGO = go_table$GO.ID
  )

  go_table$member_genes = vapply(
    go_table$GO.ID,
    function(this_go) {
      paste(
        sort(
          intersect(
            term_members[[this_go]],
            interesting_genes
          )
        ),
        collapse = ";"
      )
    },
    character(1)
  )

  raw_p = score(go_test)[go_table$GO.ID]
  stopifnot(length(raw_p) == nrow(go_table), all(is.finite(raw_p)))

  go_table = go_table %>%
    mutate(
      gene_set = set_name,
      ontology = ontology,
      n_gene_set = length(interesting_genes),
      n_universe = length(supported_analyte),
      fisher_p = unname(raw_p),
      fisher_q = p.adjust(fisher_p, method = "BH"),
      Annotated = as.numeric(Annotated),
      Significant = as.numeric(Significant),
      Expected = Annotated * length(interesting_genes) / length(supported_analyte),
      fold_enrichment = Significant / Expected
    ) %>%
    dplyr::select(
      gene_set,
      ontology,
      GO.ID,
      Term,
      n_gene_set,
      n_universe,
      Annotated,
      Significant,
      Expected,
      fold_enrichment,
      fisher_p,
      fisher_q,
      member_genes
    ) %>%
    arrange(fisher_p, GO.ID)

  list(
    data = go_data,
    table = go_table
  )
}

#####################
# RUN GO ENRICHMENT #
#####################

ontologies = c("BP", "CC", "MF")
go_objects = list()
go_tables = list()

for (set_name in names(gene_sets)) {
  for (ontology in ontologies) {

    result = run_one_go(
      interesting_genes = gene_sets[[set_name]],
      set_name = set_name,
      ontology = ontology
    )

    result_name = paste(
      set_name,
      ontology,
      sep = "_"
    )

    go_objects[[result_name]] = result$data
    go_tables[[result_name]] = result$table
  }
}

all_go_results = bind_rows(go_tables)

significant_go = all_go_results %>%
  filter(fisher_q < 0.05) %>%
  arrange(gene_set, ontology, fisher_q, fisher_p)

########################################
# ALL-TERM DISCORDANT CLASS COMPARISON  #
########################################

# Eligibility is fixed by the original 385-gene background, not significance
# or the number of discordant genes. Use topGO-propagated term membership.
discordant_genes = gene_sets$discordant
other_interaction_genes = setdiff(gene_sets$all_interaction, discordant_genes)
stopifnot(length(other_interaction_genes) == 27L,
          all(gene_sets$all_interaction %in% supported_analyte))

compare_go_term = function(term_members, discordant_genes, other_genes) {
  a = sum(discordant_genes %in% term_members)
  b = length(discordant_genes) - a
  c = sum(other_genes %in% term_members)
  d = length(other_genes) - c
  informative = (a + c > 0L) && (b + d > 0L)
  # Empty/full columns contain no class contrast. Retain them in the BH family
  # with P = 1 and an explicitly non-estimable OR.
  fit = if (informative) {
    fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
                alternative = "greater")
  } else NULL
  data.frame(
    discordant_members = a, discordant_nonmembers = b,
    other_members = c, other_nonmembers = d,
    discordant = paste0(a, "/", length(discordant_genes)),
    other_interactions = paste0(c, "/", length(other_genes)),
    odds_ratio = if (informative) unname(fit$estimate) else NA_real_,
    fisher_p = if (informative) fit$p.value else 1,
    estimable = informative,
    discordant_genes = paste(sort(intersect(discordant_genes, term_members)), collapse = ";"),
    other_interaction_genes = paste(sort(intersect(other_genes, term_members)), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

comparison_tables = lapply(ontologies, function(ontology) {
  primary = go_tables[[paste("discordant", ontology, sep = "_")]]
  members = genesInTerm(go_objects[[paste("discordant", ontology, sep = "_")]],
                        whichGO = primary$GO.ID)
  counts = bind_rows(lapply(primary$GO.ID, function(id) {
    compare_go_term(members[[id]], discordant_genes, other_interaction_genes)
  }))
  out = cbind(primary[, c("ontology", "GO.ID", "Term", "Annotated")], counts)
  out$fisher_q = p.adjust(out$fisher_p, method = "BH")
  out$n_tests = nrow(out)
  out
})
class_comparison = bind_rows(comparison_tables) %>%
  arrange(ontology, fisher_q, fisher_p, GO.ID)

# Report every first-pass discordant hit, including nonsignificant follow-ups.
# Join on stable GO ID plus ontology, never on a display term name.
discordant_summary = significant_go %>%
  filter(gene_set == "discordant") %>%
  left_join(
    class_comparison %>% dplyr::select(
      ontology, GO.ID, other_interactions, odds_ratio, estimable,
      comparison_p = fisher_p, comparison_q = fisher_q,
      comparison_n_tests = n_tests
    ),
    by = c("ontology", "GO.ID")
  ) %>%
  mutate(significant_both = fisher_q < 0.05 & comparison_q < 0.05) %>%
  arrange(ontology, fisher_q, GO.ID)
stopifnot(nrow(discordant_summary) == sum(significant_go$gene_set == "discordant"),
          !anyNA(discordant_summary$comparison_q),
          !anyDuplicated(class_comparison[, c("ontology", "GO.ID")]))

########################################
# MKT1 AND BUL2 DESCRIPTIVE SUMMARY     #
########################################

marker_parts = strsplit(
  as.character(qtls_aug_analyte_ixn$marker),
  ":",
  fixed = TRUE
)
if (any(lengths(marker_parts) != 2L)) {
  stop("Could not parse all strict-interaction peak markers.")
}

qtl_chr = vapply(
  marker_parts,
  function(x) x[[1]],
  character(1)
)
qtl_position = suppressWarnings(
  as.numeric(
    vapply(
      marker_parts,
      function(x) x[[2]],
      character(1)
    )
  )
)
if (any(!is.finite(qtl_position))) {
  stop("Could not parse strict-interaction peak positions.")
}

hotspot_windows = data.frame(
  hotspot = c("MKT1", "BUL2"),
  chromosome = c("chrXIV", "chrXIII"),
  locus_start = c(431063L, 28519L),
  locus_end = c(532790L, 79575L),
  stringsAsFactors = FALSE
)

hotspot_summary_rows = lapply(
  seq_len(nrow(hotspot_windows)),
  function(i) {
    window = hotspot_windows[i, , drop = FALSE]
    keep =
      qtl_chr == window$chromosome &
      qtl_position >= window$locus_start &
      qtl_position <= window$locus_end

    dat = qtls_aug_analyte_ixn[keep, , drop = FALSE]

    n_discordant = sum(dat$bin == "line cross")
    n_mRNA_specific = sum(dat$bin == "mRNA specific")
    n_protein_specific = sum(dat$bin == "protein specific")
    n_same_direction = sum(dat$bin == "same direction")

    discordant = dat[dat$bin == "line cross", , drop = FALSE]

    data.frame(
      hotspot = window$hotspot,
      chromosome = window$chromosome,
      locus_start = window$locus_start,
      locus_end = window$locus_end,
      strict_interaction_QTLs = nrow(dat),
      percent_of_strict_64 = 100 * nrow(dat) / 64,
      discordant_QTLs = n_discordant,
      percent_of_discordant_37 = 100 * n_discordant / 37,
      mRNA_specific_QTLs = n_mRNA_specific,
      protein_specific_QTLs = n_protein_specific,
      same_direction_QTLs = n_same_direction,
      discordant_mRNA_down_protein_up = sum(
        discordant$mRNA_effect < 0 & discordant$protein_effect > 0
      ),
      discordant_mRNA_up_protein_down = sum(
        discordant$mRNA_effect > 0 & discordant$protein_effect < 0
      ),
      stringsAsFactors = FALSE
    )
  }
)

hotspot_summary = do.call(
  rbind,
  hotspot_summary_rows
)
rownames(hotspot_summary) = NULL

expected_hotspot_counts = data.frame(
  hotspot = c("MKT1", "BUL2"),
  strict_interaction_QTLs = c(36L, 11L),
  discordant_QTLs = c(26L, 1L),
  mRNA_specific_QTLs = c(8L, 0L),
  protein_specific_QTLs = c(2L, 9L),
  same_direction_QTLs = c(0L, 1L),
  discordant_mRNA_down_protein_up = c(26L, 1L),
  discordant_mRNA_up_protein_down = c(0L, 0L),
  stringsAsFactors = FALSE
)

hotspot_order = match(
  expected_hotspot_counts$hotspot,
  hotspot_summary$hotspot
)
if (any(is.na(hotspot_order))) {
  stop("MKT1/BUL2 hotspot summary is incomplete.")
}
hotspot_check = hotspot_summary[
  hotspot_order,
  names(expected_hotspot_counts),
  drop = FALSE
]

for (column in setdiff(names(expected_hotspot_counts), "hotspot")) {
  if (!identical(
    as.integer(hotspot_check[[column]]),
    expected_hotspot_counts[[column]]
  )) {
    stop("Unexpected MKT1/BUL2 count in column: ", column)
  }
}

################
# SAVE OUTPUTS #
################

write.table(
  significant_go,
  file = paste0(
    go_asys,
    "GO_significant_terms.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  discordant_summary,
  file = paste0(
    go_asys,
    "GO_discordant_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  class_comparison,
  file = paste0(
    go_asys,
    "GO_discordant_vs_other_interactions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  hotspot_summary,
  file = paste0(
    go_asys,
    "MKT1_BUL2_hotspot_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("")
message("GO review analysis complete.")
message("Output directory: ", go_asys)
message("")
message("All significant first-pass discordant terms with class comparisons:")
print(discordant_summary, row.names = FALSE)
message("")
message("MKT1/BUL2 support summary:")
print(hotspot_summary, row.names = FALSE)

