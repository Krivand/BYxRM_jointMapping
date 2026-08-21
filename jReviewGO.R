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

go_asys = paste0(asys, "review/go/")
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
  "bin"
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

parse_topgo_p = function(x) {
  x = as.character(x)
  x = gsub("^<\\s*", "", x)
  suppressWarnings(as.numeric(x))
}

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

  go_table = go_table %>%
    mutate(
      gene_set = set_name,
      ontology = ontology,
      n_gene_set = length(interesting_genes),
      n_universe = length(supported_analyte),
      fisher_p = parse_topgo_p(classicFisher),
      fisher_q = p.adjust(fisher_p, method = "BH"),
      Annotated = as.numeric(Annotated),
      Significant = as.numeric(Significant),
      Expected = as.numeric(Expected),
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

############################
# REVIEWER-FOCUSED GO ROWS #
############################

key_terms = c(
  "GO:0006412",
  "GO:0042254",
  "GO:0044391",
  "GO:0003735"
)

key_go = all_go_results %>%
  filter(
    gene_set == "discordant",
    GO.ID %in% key_terms
  ) %>%
  mutate(
    term_order = match(GO.ID, key_terms),
    `GO term` = Term,
    Ontology = ontology,
    `Discordant genes` = paste0(Significant, "/37"),
    Expected = sprintf("%.2f", Expected),
    `Fold enrichment` = sprintf("%.2f", fold_enrichment),
    FDR = formatC(fisher_q, format = "e", digits = 2)
  ) %>%
  arrange(term_order) %>%
  dplyr::select(
    `GO term`,
    Ontology,
    `Discordant genes`,
    Expected,
    `Fold enrichment`,
    FDR
  )

if (
  nrow(key_go) != 4L ||
  !identical(
    key_go$`Discordant genes`,
    c("25/37", "15/37", "21/37", "20/37")
  )
) {
  stop("Unexpected counts for reviewer-focused GO terms.")
}

##############################
# RIBOSOME CLASS CONTRAST    #
##############################

ribosome_focus_terms = data.frame(
  ontology = c(
    "CC",
    "CC",
    "MF",
    "BP",
    "BP",
    "BP"
  ),
  GO.ID = c(
    "GO:0005840",
    "GO:0044391",
    "GO:0003735",
    "GO:0006412",
    "GO:0002181",
    "GO:0042254"
  ),
  Term = c(
    "ribosome",
    "ribosomal subunit",
    "structural constituent of ribosome",
    "translation",
    "cytoplasmic translation",
    "ribosome biogenesis"
  ),
  stringsAsFactors = FALSE
)

discordant_genes = gene_sets$discordant
other_interaction_genes = setdiff(
  gene_sets$all_interaction,
  discordant_genes
)

ribosome_contrast_rows = lapply(
  seq_len(nrow(ribosome_focus_terms)),
  function(i) {

    this_term = ribosome_focus_terms[i, , drop = FALSE]
    result_name = paste(
      "all_interaction",
      this_term$ontology,
      sep = "_"
    )

    term_members = genesInTerm(
      go_objects[[result_name]],
      whichGO = this_term$GO.ID
    )[[this_term$GO.ID]]

    discordant_members = sort(
      intersect(discordant_genes, term_members)
    )

    other_members = sort(
      intersect(other_interaction_genes, term_members)
    )

    comparison_table = matrix(
      c(
        length(discordant_members),
        length(discordant_genes) - length(discordant_members),
        length(other_members),
        length(other_interaction_genes) - length(other_members)
      ),
      nrow = 2,
      byrow = TRUE
    )

    fisher_result = fisher.test(
      comparison_table,
      alternative = "greater"
    )

    data.frame(
      GO.ID = this_term$GO.ID,
      GO_term = this_term$Term,
      discordant = paste0(length(discordant_members), "/37"),
      other_interactions = paste0(length(other_members), "/27"),
      odds_ratio = unname(fisher_result$estimate),
      fisher_p = fisher_result$p.value,
      stringsAsFactors = FALSE
    )
  }
)

ribosome_contrast = do.call(
  rbind,
  ribosome_contrast_rows
)

ribosome_contrast$fisher_q = p.adjust(
  ribosome_contrast$fisher_p,
  method = "BH"
)

ribosomal_subunit_row = ribosome_contrast[
  ribosome_contrast$GO.ID == "GO:0044391",
  ,
  drop = FALSE
]

if (
  ribosomal_subunit_row$discordant != "21/37" ||
  ribosomal_subunit_row$other_interactions != "1/27"
) {
  stop("Unexpected ribosomal-subunit class contrast.")
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
  key_go,
  file = paste0(
    go_asys,
    "GO_key_terms.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ribosome_contrast,
  file = paste0(
    go_asys,
    "GO_discordant_vs_other_interactions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("")
message("GO review analysis complete.")
message("Output directory: ", go_asys)
message("")
message("Reviewer-focused GO table:")
print(key_go, row.names = FALSE)
