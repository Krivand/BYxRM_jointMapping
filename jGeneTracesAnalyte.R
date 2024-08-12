############################################
# SOURCE IMPORTANT FUNCTIONS AND CONSTANTS #
############################################

source("/home/albertf/vandy147/p05/remap_functions.R")

###################
# LOAD IN OBJECTS #
###################

# LOAD IN UNFILTERED ASSOCATIONS AND PREPARE
load(paste0(imed,"raAnalyte.rda"))

##################
# ANALYTE TRACES #
##################

raAnalyte %>% 
  rename(marker = geno.marker, gene = geno.gene) %>%
  select(-matches("\\.marker$"), 
         -matches("\\.gene$"), 
         -matches("\\.flavor$")) %>%
  mutate(pairID = paste(gene, marker, sep = "_")) %>% 
  select(pairID, gene, marker, everything()) %>%
  rename_with(
    .fn = ~ sub("\\.effect$", ".allelic_effect", .x),
    .cols = ends_with(".effect")
  ) %>%
  select(-starts_with("geno")) %>%
  mutate(across(
    .cols = all_of(paste0(c("analyte", "ixn", "protein", "mRNA"), ".pval")),
    .fns = ~ -log10(.) * 
      sign(get(sub(".pval$", ".allelic_effect", cur_column()))),
    .names = "{sub('.pval$', '', .col)}.smltp"
  )) %>%
  rename_with(~ str_replace_all(., "ixn", "interaction")) %>%
  rename_with(~ str_remove_all(., '.smltp')) %>%
  select(-matches("\\.")) %>%
  tracePlotter(
    fn = paste0(asys, 'trace_plots'), 
    plotWidth = 7.5, 
    plotHeight = 2.5
  )