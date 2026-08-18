########################################
########################################
## GENETICS REVIEW GO ENRICHMENT      ##
## JOINT ANALYTE INTERACTION GENE SETS ##
########################################
########################################

#################
# PRELIMINARIES #
#################

# source functions and constants
source("/PHShome/kv1026/p05/remap_functions.R")

# packages specific to reviewer analyses
review_pkgs = c(
  "AnnotationDbi",
  "org.Sc.sgd.db",
  "topGO"
)

for (pkg in review_pkgs) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    message("Installing ", pkg, "...")

    suppressMessages(
      BiocManager::install(
        pkg,
        lib = personal_lib_path,
        ask = FALSE,
        update = FALSE,
        quiet = TRUE
      )
    )
  }

  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop("Package ", pkg, " could not be loaded.")
  }
}

rm(review_pkgs)

#################
# OUTPUT SETUP  #
#################

review_imed = paste0(imed, "review/")
review_asys = paste0(asys, "review/")

ensure_directory(review_imed)
ensure_directory(review_asys)

go_imed = paste0(review_imed, "go/")
go_asys = paste0(review_asys, "go/")

ensure_directory(go_imed)
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

#################
# INPUT CHECKS  #
#################

if (length(supported_analyte) != 385L) {
  stop(
    "Expected 385 supported analyte genes, found ",
    length(supported_analyte),
    "."
  )
}

if (anyDuplicated(supported_analyte)) {
  stop("supported_analyte contains duplicated genes.")
}

if (nrow(qtls_aug_analyte_ixn) != 64L) {
  stop(
    "Expected 64 analyte interaction QTLs, found ",
    nrow(qtls_aug_analyte_ixn),
    "."
  )
}

ixn_genes = unique(qtls_aug_analyte_ixn$gene)

if (!all(ixn_genes %in% supported_analyte)) {
  stop(
    "Some analyte interaction genes are outside ",
    "the supported analyte universe."
  )
}

gene_bin_check =
  qtls_aug_analyte_ixn %>%
  distinct(gene, bin) %>%
  count(gene, name = "n_bins")

if (any(gene_bin_check$n_bins != 1L)) {
  stop("At least one interaction gene belongs to multiple bins.")
}

###################
# BUILD GENE SETS #
###################

gene_sets = list(
  all_interaction = sort(
    unique(qtls_aug_analyte_ixn$gene)
  ),
  discordant = sort(
    unique(
      qtls_aug_analyte_ixn$gene[
        qtls_aug_analyte_ixn$bin == "line cross"
      ]
    )
  ),
  protein_specific = sort(
    unique(
      qtls_aug_analyte_ixn$gene[
        qtls_aug_analyte_ixn$bin == "protein specific"
      ]
    )
  ),
  mRNA_specific = sort(
    unique(
      qtls_aug_analyte_ixn$gene[
        qtls_aug_analyte_ixn$bin == "mRNA specific"
      ]
    )
  ),
  same_direction = sort(
    unique(
      qtls_aug_analyte_ixn$gene[
        qtls_aug_analyte_ixn$bin == "same direction"
      ]
    )
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

  stop(
    "Unexpected review gene-set sizes:\n",
    paste(
      names(observed_sizes),
      observed_sizes,
      collapse = "\n"
    )
  )
}

message("Review gene sets:")
print(observed_sizes)

########################
# BUILD SGD GO MAPPING #
########################

message("Collecting SGD GO annotations...")

go_annotation =
  AnnotationDbi::select(
    org.Sc.sgd.db::org.Sc.sgd.db,
    keys = supported_analyte,
    keytype = "ORF",
    columns = c(
      "GO",
      "ONTOLOGY"
    )
  ) %>%
  as.data.frame() %>%
  filter(
    !is.na(ORF),
    !is.na(GO),
    !is.na(ONTOLOGY),
    ORF %in% supported_analyte,
    ONTOLOGY %in% c("BP", "CC", "MF")
  ) %>%
  distinct(
    ORF,
    GO,
    ONTOLOGY
  )

gene2GO =
  go_annotation %>%
  dplyr::select(
    ORF,
    GO
  ) %>%
  distinct()

gene2GO = split(
  gene2GO$GO,
  gene2GO$ORF
)

annotated_genes = intersect(
  supported_analyte,
  names(gene2GO)
)

message(
  "GO annotated genes in 385-gene universe: ",
  length(annotated_genes),
  " / ",
  length(supported_analyte)
)

saveRDS(
  go_annotation,
  paste0(go_imed, "sgd_go_annotation.rds")
)

####################
# HELPER FUNCTIONS #
####################

parse_topgo_p = function(x) {

  x = as.character(x)
  x = gsub("^<\\s*", "", x)

  suppressWarnings(
    as.numeric(x)
  )
}

run_one_go = function(
  interesting_genes,
  set_name,
  ontology,
  universe,
  gene2GO,
  node_size = 3L
) {

  if (!all(interesting_genes %in% universe)) {
    stop(
      "Gene set ",
      set_name,
      " contains genes outside the universe."
    )
  }

  gene_score = as.integer(
    universe %in% interesting_genes
  )

  names(gene_score) = universe

  go_data = new(
    "topGOdata",
    ontology = ontology,
    allGenes = gene_score,
    geneSel = function(x) x == 1L,
    annot = annFUN.gene2GO,
    gene2GO = gene2GO,
    nodeSize = node_size
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

      members = intersect(
        term_members[[this_go]],
        interesting_genes
      )

      paste(
        sort(members),
        collapse = ";"
      )
    },
    character(1)
  )

  go_table =
    go_table %>%
    mutate(
      gene_set = set_name,
      ontology = ontology,
      n_gene_set = length(interesting_genes),
      n_universe = length(universe),
      fisher_p = parse_topgo_p(classicFisher),
      fisher_q = p.adjust(
        fisher_p,
        method = "BH"
      ),
      Annotated = as.numeric(Annotated),
      Significant = as.numeric(Significant),
      Expected = as.numeric(Expected),
      fold_enrichment = ifelse(
        Expected > 0,
        Significant / Expected,
        NA_real_
      )
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
      member_genes,
      classicFisher
    ) %>%
    arrange(
      fisher_p,
      GO.ID
    )

  list(
    data = go_data,
    test = go_test,
    table = go_table
  )
}

#####################
# RUN GO ENRICHMENT #
#####################

ontologies = c(
  "BP",
  "CC",
  "MF"
)

go_results = list()
go_tables = list()

for (set_name in names(gene_sets)) {

  interesting_genes = gene_sets[[set_name]]

  message(
    "Running GO enrichment for ",
    set_name,
    " [n = ",
    length(interesting_genes),
    "]"
  )

  for (ontology in ontologies) {

    message("  ", ontology)

    result = run_one_go(
      interesting_genes = interesting_genes,
      set_name = set_name,
      ontology = ontology,
      universe = supported_analyte,
      gene2GO = gene2GO,
      node_size = 3L
    )

    result_name = paste(
      set_name,
      ontology,
      sep = "_"
    )

    go_results[[result_name]] = result
    go_tables[[result_name]] = result$table
  }
}

all_go_results =
  bind_rows(go_tables) %>%
  arrange(
    gene_set,
    ontology,
    fisher_p,
    GO.ID
  )

#####################
# SUMMARIZE RESULTS #
#####################

go_summary =
  all_go_results %>%
  group_by(
    gene_set,
    ontology
  ) %>%
  summarise(
    n_gene_set = first(n_gene_set),
    n_terms_tested = n(),
    min_fisher_p = min(
      fisher_p,
      na.rm = TRUE
    ),
    n_nominal_0.05 = sum(
      fisher_p < 0.05,
      na.rm = TRUE
    ),
    n_fdr_0.05 = sum(
      fisher_q < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

go_report =
  all_go_results %>%
  filter(
    fisher_q < 0.05
  ) %>%
  arrange(
    gene_set,
    ontology,
    fisher_q,
    fisher_p
  )

##############################
# RIBOSOME CLASS CONTRAST    #
##############################

# Exploratory follow-up asking whether representative
# ribosome/translation annotations are concentrated in
# discordant genes among the 64 interaction genes.

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

discordant_genes =
  gene_sets$discordant

other_interaction_genes =
  setdiff(
    gene_sets$all_interaction,
    discordant_genes
  )

ribosome_contrast_rows =
  lapply(
    seq_len(
      nrow(
        ribosome_focus_terms
      )
    ),
    function(i) {

      this_term =
        ribosome_focus_terms[
          i,
          ,
          drop = FALSE
        ]

      result_name =
        paste(
          "all_interaction",
          this_term$ontology,
          sep = "_"
        )

      term_members =
        genesInTerm(
          go_results[[result_name]]$data,
          whichGO = this_term$GO.ID
        )[[this_term$GO.ID]]

      discordant_members =
        sort(
          intersect(
            discordant_genes,
            term_members
          )
        )

      other_members =
        sort(
          intersect(
            other_interaction_genes,
            term_members
          )
        )

      comparison_table =
        matrix(
          c(
            length(
              discordant_members
            ),
            length(
              discordant_genes
            ) -
              length(
                discordant_members
              ),
            length(
              other_members
            ),
            length(
              other_interaction_genes
            ) -
              length(
                other_members
              )
          ),
          nrow = 2,
          byrow = TRUE
        )

      fisher_result =
        fisher.test(
          comparison_table,
          alternative = "greater"
        )

      data.frame(
        ontology =
          this_term$ontology,
        GO.ID =
          this_term$GO.ID,
        Term =
          this_term$Term,
        n_discordant =
          length(
            discordant_genes
          ),
        n_discordant_member =
          length(
            discordant_members
          ),
        prop_discordant =
          length(
            discordant_members
          ) /
          length(
            discordant_genes
          ),
        n_other_interaction =
          length(
            other_interaction_genes
          ),
        n_other_interaction_member =
          length(
            other_members
          ),
        prop_other_interaction =
          length(
            other_members
          ) /
          length(
            other_interaction_genes
          ),
        odds_ratio =
          unname(
            fisher_result$estimate
          ),
        fisher_p =
          fisher_result$p.value,
        discordant_genes =
          paste(
            discordant_members,
            collapse = ";"
          ),
        other_interaction_genes =
          paste(
            other_members,
            collapse = ";"
          ),
        stringsAsFactors = FALSE
      )
    }
  )

ribosome_class_contrast =
  do.call(
    rbind,
    ribosome_contrast_rows
  )

ribosome_class_contrast$fisher_q =
  p.adjust(
    ribosome_class_contrast$fisher_p,
    method = "BH"
  )

gene_set_manifest =
  bind_rows(
    lapply(
      names(gene_sets),
      function(set_name) {
        data.frame(
          gene_set = set_name,
          gene = gene_sets[[set_name]],
          stringsAsFactors = FALSE
        )
      }
    )
  ) %>%
  arrange(
    gene_set,
    gene
  )

go_annotation_summary = data.frame(
  n_universe = length(supported_analyte),
  n_universe_go_annotated = length(annotated_genes),
  frac_universe_go_annotated =
    length(annotated_genes) /
    length(supported_analyte),
  stringsAsFactors = FALSE
)

################
# SAVE OUTPUTS #
################

write.table(
  all_go_results,
  file = paste0(
    go_asys,
    "analyte_interaction_GO_full.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  go_summary,
  file = paste0(
    go_asys,
    "analyte_interaction_GO_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  go_report,
  file = paste0(
    go_asys,
    "analyte_interaction_GO_report.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_set_manifest,
  file = paste0(
    go_asys,
    "analyte_interaction_gene_sets.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  go_annotation_summary,
  file = paste0(
    go_asys,
    "GO_annotation_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ribosome_class_contrast,
  file = paste0(
    go_asys,
    "analyte_interaction_ribosome_class_contrast.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

go_workbook = list(
  summary = go_summary,
  significant_terms = go_report,
  ribosome_focus = ribosome_class_contrast,
  gene_sets = gene_set_manifest,
  annotation_summary = go_annotation_summary,
  all_terms = all_go_results
)

openxlsx::write.xlsx(
  go_workbook,
  file = paste0(
    go_asys,
    "analyte_interaction_GO_review.xlsx"
  ),
  asTable = TRUE,
  overwrite = TRUE
)

saveRDS(
  go_results,
  paste0(
    go_imed,
    "analyte_interaction_GO_results.rds"
  )
)

message("")
message("GO enrichment complete.")
message("Output directory: ", go_asys)
message("")
print(go_summary)
