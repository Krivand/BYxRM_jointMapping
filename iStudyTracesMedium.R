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

# NOTE: this line referenced medMeta inside its own creation.
# Keeping structure; set the intended subset if needed.
medMeta = classicMetadata %>% filter(study %in% medMeta$study)

doubleSupp = suppTwo(medMeta, type = "medium")
print(paste0('Within the RNA-only subset there are ',
             length(doubleSupp),
             ' genes supported by two SC and two YNB2G datasets'))

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

plotMeta = classicMetadata %>% mutate(study = gsub("_", " ", study)) %>%
  mutate(medium = ifelse(study == "Smith 2008 Ethanol", "YNB2E", medium)) %>%
  mutate(medium = factor(medium, levels = c("SC", "YNB2G", "YNB2E")))

studyTracer(studies = studies, 
            smeta = plotMeta,
            suppGenes = doubleSupp, 
            fn = paste0(asys, 'individual_study'),
            plotHeight = 5,
            plotWidth = 7.5,
            traceType = "b",
            fixef = "medium")
