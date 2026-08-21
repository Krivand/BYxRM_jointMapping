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
  "mRNA_effect",
  "protein_effect",
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

normalize_chr = function(x) {

  x = toupper(trimws(as.character(x)))
  x = sub("^CHR", "", x)

  roman = c(
    "I", "II", "III", "IV",
    "V", "VI", "VII", "VIII",
    "IX", "X", "XI", "XII",
    "XIII", "XIV", "XV", "XVI"
  )

  out = suppressWarnings(as.integer(x))
  use_roman = is.na(out) & x %in% roman
  out[use_roman] = match(x[use_roman], roman)
  out
}

parse_marker = function(x) {

  x = as.character(x)

  chr = sub(":.*$", "", x)
  position = sub("^[^:]+:", "", x)
  position = sub("_.*$", "", position)

  data.frame(
    chr = normalize_chr(chr),
    position = suppressWarnings(as.numeric(position)),
    stringsAsFactors = FALSE
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
    term_order = match(GO.ID, key_terms)
  ) %>%
  arrange(term_order) %>%
  dplyr::select(-term_order)

if (
  nrow(key_go) != 4L ||
  !identical(
    key_go$Significant,
    c(25, 15, 21, 20)
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

ribosome_term_members = vector(
  "list",
  nrow(ribosome_focus_terms)
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

    ribosome_term_members[[i]] <<- term_members

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
      ontology = this_term$ontology,
      GO.ID = this_term$GO.ID,
      Term = this_term$Term,
      n_discordant = length(discordant_genes),
      n_discordant_member = length(discordant_members),
      n_other_interaction = length(other_interaction_genes),
      n_other_interaction_member = length(other_members),
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
  ribosomal_subunit_row$n_discordant_member != 21L ||
  ribosomal_subunit_row$n_other_interaction_member != 1L
) {
  stop("Unexpected ribosomal-subunit class contrast.")
}

ribosome_related_genes = intersect(
  gene_sets$all_interaction,
  unique(unlist(ribosome_term_members))
)

########################
# MKT1 / BUL2 SUMMARY  #
########################

# These are the 100-kb hotspot bins used for the interaction-QTL
# hotspot summaries in the manuscript.

peak = parse_marker(
  qtls_aug_analyte_ixn$marker
)

if (
  any(!is.finite(peak$chr)) ||
  any(!is.finite(peak$position))
) {
  stop("Could not parse at least one interaction-QTL peak marker.")
}

hotspot = rep(
  NA_character_,
  nrow(qtls_aug_analyte_ixn)
)

hotspot[
  peak$chr == 14 &
  peak$position >= 400000 &
  peak$position < 500000
] = "MKT1"

hotspot[
  peak$chr == 13 &
  peak$position >= 0 &
  peak$position < 100000
] = "BUL2"

summarize_hotspot = function(hotspot_name) {

  x = qtls_aug_analyte_ixn[
    hotspot == hotspot_name,
    ,
    drop = FALSE
  ]

  discordant = x$bin == "line cross"

  data.frame(
    hotspot = hotspot_name,
    n_strict_interaction = nrow(x),
    n_discordant = sum(discordant),
    n_mRNA_specific = sum(x$bin == "mRNA specific"),
    n_protein_specific = sum(x$bin == "protein specific"),
    n_same_direction = sum(x$bin == "same direction"),
    n_discordant_ribosome_related = sum(
      discordant &
      x$gene %in% ribosome_related_genes
    ),
    n_discordant_RM_RNA_down_protein_up = sum(
      discordant &
      x$mRNA_effect < 0 &
      x$protein_effect > 0
    ),
    stringsAsFactors = FALSE
  )
}

hotspot_summary = rbind(
  summarize_hotspot("MKT1"),
  summarize_hotspot("BUL2")
)

mkt1 = hotspot_summary[
  hotspot_summary$hotspot == "MKT1",
  ,
  drop = FALSE
]

bul2 = hotspot_summary[
  hotspot_summary$hotspot == "BUL2",
  ,
  drop = FALSE
]

if (
  mkt1$n_strict_interaction != 36L ||
  mkt1$n_discordant != 26L ||
  mkt1$n_discordant_ribosome_related != 21L ||
  mkt1$n_discordant_RM_RNA_down_protein_up != 26L
) {
  stop("Unexpected MKT1 interaction summary.")
}

if (
  bul2$n_strict_interaction != 11L ||
  bul2$n_protein_specific != 9L
) {
  stop("Unexpected BUL2 interaction summary.")
}

################
# SAVE OUTPUTS #
################

write.table(
  significant_go,
  file = paste0(
    go_asys,
    "analyte_interaction_GO_report.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  key_go,
  file = paste0(
    go_asys,
    "analyte_interaction_GO_key.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ribosome_contrast,
  file = paste0(
    go_asys,
    "analyte_interaction_ribosome_class_contrast.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  hotspot_summary,
  file = paste0(
    go_asys,
    "analyte_interaction_hotspot_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("")
message("GO review analysis complete.")
message("Production tables written to: ", go_asys)
message("")
message("Reviewer-focused GO terms:")
print(key_go)
message("")
message("Hotspot summary:")
print(hotspot_summary)
