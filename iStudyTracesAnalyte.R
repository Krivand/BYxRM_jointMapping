############################################
# SOURCE IMPORTANT FUNCTIONS AND CONSTANTS #
############################################

source("/PHShome/kv1026/p05/remap_functions.R")  # defines 'imed', 'asys', etc.

###################
# LOAD IN OBJECTS #
###################

studies = lapply(studyList, function(x) {
  load(paste0(imed, x, '.rda'))
  get(x)
})
names(studies) = studyList

classicMetadata$genes = lapply(classicMetadata$study,
                               function(x) colnames(studies[[x]]$pheno))
classicMetadata$segs  = lapply(classicMetadata$study,
                               function(x) rownames(studies[[x]]$pheno))

# NOTE: this line referenced scaMetadata inside its own creation.
# Keeping structure, but if you intended a predefined subset, adjust accordingly.
scaMetadata = classicMetadata %>% filter(study %in% scaMetadata$study)

doubleSupp = suppTwo(scaMetadata, type = "analyte")
print(paste0('Within the SC media with genotype A subset there are ',
             length(doubleSupp),
             ' genes supported by two mRNA and two protein datasets'))

# naming cleanups
studies = studies %>%
  {.[!names(.) %in% "Albert_2018"]} %>% 
  {setNames(., str_replace_all(names(.), "AlbertAlt_2018", "Albert_2018"))} %>% 
  {setNames(., str_replace_all(names(.), "SmithEth_2008", "Smith_2008_Ethanol"))} %>% 
  {setNames(., str_replace_all(names(.), "SmithGlc_2008", "Smith_2008_Glucose"))} %>% 
  {setNames(., str_replace_all(names(.), "GrossbachRNA_2022", "Grossbach_2022_mRNA"))} %>% 
  {setNames(., str_replace_all(names(.), "GrossbachPRT_2022", "Grossbach_2022_protein"))} %>% 
  {setNames(., str_replace_all(names(.), "_", " "))}

classicMetadata = classicMetadata %>%
  filter(study != "Albert_2018") %>%
  mutate(study = str_replace_all(study, "AlbertAlt_2018", "Albert_2018")) %>%
  mutate(study = str_replace_all(study, "SmithEth_2008", "Smith_2008_Ethanol")) %>%
  mutate(study = str_replace_all(study, "SmithGlc_2008", "Smith_2008_Glucose")) %>%
  mutate(study = str_replace_all(study, "GrossbachRNA_2022", "Grossbach_2022_mRNA")) %>%
  mutate(study = str_replace_all(study, "GrossbachPRT_2022", "Grossbach_2022_protein")) %>%
  mutate(study = str_replace_all(study, "_", " "))

###############
# TRACE PLOTS #
###############

exclusionMetadata = classicMetadata %>% mutate(study = gsub("_", " ", study)) %>%
  mutate(analyte = ifelse(
    study %in% c("Smith 2008 Ethanol", "Smith 2008 Glucose", "Albert 2018"),
    "excluded", analyte)) %>%
  mutate(analyte = factor(analyte, levels = c("mRNA", "protein", "excluded")))

studyTracer(studies = studies, 
            smeta = exclusionMetadata,
            suppGenes = doubleSupp, 
            fn = paste0(asys, 'individual_study'),
            plotHeight = 5,
            plotWidth = 7.5,
            traceType = "b",
            fixef = "analyte")
