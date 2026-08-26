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
# RIBOSOME CLASS COMPARISON  #
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

ribosome_comparison_rows = lapply(
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
      discordant_genes = paste(discordant_members, collapse = ";"),
      other_interaction_genes = paste(other_members, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
)

ribosome_contrast = do.call(
  rbind,
  ribosome_comparison_rows
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
  stop("Unexpected ribosomal-subunit class comparison.")
}

########################################
# MKT1 AND BUL2 SUPPORT SUMMARY        #
########################################

ribosome_member_text = c(
  ribosome_contrast$discordant_genes,
  ribosome_contrast$other_interaction_genes
)
ribosome_member_text = ribosome_member_text[
  !is.na(ribosome_member_text) & nzchar(ribosome_member_text)
]
ribosome_related_genes = sort(
  unique(
    unlist(
      strsplit(
        ribosome_member_text,
        ";",
        fixed = TRUE
      )
    )
  )
)
ribosome_related_genes = ribosome_related_genes[
  nzchar(ribosome_related_genes)
]

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

class_enrichment_p = function(in_hotspot, n_hotspot, class_total) {
  outside_class = class_total - in_hotspot
  outside_total = 64L - n_hotspot
  table_2x2 = matrix(
    c(
      in_hotspot,
      n_hotspot - in_hotspot,
      outside_class,
      outside_total - outside_class
    ),
    nrow = 2,
    byrow = TRUE
  )
  fisher.test(
    table_2x2,
    alternative = "greater"
  )$p.value
}

hotspot_summary_rows = lapply(
  seq_len(nrow(hotspot_windows)),
  function(i) {
    window = hotspot_windows[i, , drop = FALSE]
    keep =
      qtl_chr == window$chromosome &
      qtl_position >= window$locus_start &
      qtl_position <= window$locus_end

    dat = qtls_aug_analyte_ixn[keep, , drop = FALSE]
    dat$ribosome_related = dat$gene %in% ribosome_related_genes

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
      discordant_ribosome_associated = sum(discordant$ribosome_related),
      discordant_nonribosome = sum(!discordant$ribosome_related),
      discordant_mRNA_down_protein_up = sum(
        discordant$mRNA_effect < 0 & discordant$protein_effect > 0
      ),
      discordant_mRNA_up_protein_down = sum(
        discordant$mRNA_effect > 0 & discordant$protein_effect < 0
      ),
      discordant_enrichment_fisher_p = class_enrichment_p(
        n_discordant,
        nrow(dat),
        expected_sizes[["discordant"]]
      ),
      mRNA_specific_enrichment_fisher_p = class_enrichment_p(
        n_mRNA_specific,
        nrow(dat),
        expected_sizes[["mRNA_specific"]]
      ),
      protein_specific_enrichment_fisher_p = class_enrichment_p(
        n_protein_specific,
        nrow(dat),
        expected_sizes[["protein_specific"]]
      ),
      ribosome_definition = paste0(
        "Member of >=1 of the six ribosome/translation GO terms in Ribosome comparison"
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
  discordant_ribosome_associated = c(21L, 0L),
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
message("Reviewer-focused GO table:")
print(key_go, row.names = FALSE)
message("")
message("MKT1/BUL2 support summary:")
print(hotspot_summary, row.names = FALSE)
