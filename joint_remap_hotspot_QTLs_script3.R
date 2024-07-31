#############################
#############################
## JOINT DATASET REMAPPING ##
## FDR 5% ASSOCIATIONS     ##
## HOTSPOT MARKERS         ##
#############################
#############################

#########################################
# RETRIEVE AND FORMAT NECESSARY OBJECTS #
#########################################

# SOURCE FXNS AND CONSTANTS
source(paste0("/home/albertf/vandy147/p05/",
              "remap_functions.R"))

# create an object for easier joint mapping
studies = lapply(scaMetadata$study, function(x) {
  load(paste0(imed, x, '.rda'))
  get(x)
}) %>% `names<-`(scaMetadata$study) 

# retrieve background QTLs across studies
bgQTLs = lapply(scaMetadata$study, function(x) {
  load(paste0(imed, x, '_Background.rda'))
  get(paste0(x,'_Background'))
}) %>% `names<-`(scaMetadata$study) 

# add gene support to metadata
scaMetadata$genes = lapply(scaMetadata$study,function(x)
  colnames(studies[[x]]$pheno))

# add segregant support to metadata
scaMetadata$segs = lapply(scaMetadata$study,function(x)
  rownames(studies[[x]]$pheno))

# within the synthetic complete media (sc) datasets
# find which genes are supported by at least two datasets 
doubleSupp = suppTwo(scaMetadata)
print(paste0('Within the SC media with genotype A subset there are ',
             length(doubleSupp),
             ' genes supported by two mRNA and two protein datasets'))

# take the union of all background QTLs discovered across all datasets
if (!file.exists(paste0(imed, "unionBg.rda"))) {
  
  unionBg = collapseBgQTLs(bgQTLs)
  save(unionBg,file=paste0(imed,"unionBg.rda"))
  
} else {
  
  print('Union of background QTLs already processed and loading now...')
  load(file=paste0(imed,"unionBg.rda"))
  
}

# how many strong QTLs are most genes experiencing from other chroms?
# currently so far the worst is 48 for 2 genes
# not worth switching to random effects due to covars approaching sample num
sapply(unionBg,function(x)sapply(x,function(y)length(y))) %>% 
  unlist() %>% table()

# PREPARE THE HOTSPOT MARKERS FOR ASSOCATION

# hotspot data
hotspots = readxl::read_excel(
  paste0(raw,'SI_Data_08_hotspotTableForPaper_topGenes100_170420.xlsx')) %>%
  as.data.frame()

# get the indices of hotspots approximated in Brem-coordinates
gr_studies = markRanger(colnames(studies$Brem_2002$gdata))
gr_hotspots = markRanger(hotspots$hotspotMarker)
nearest_indices = nearest(gr_hotspots, gr_studies)
portedHotspots = gr_studies[nearest_indices] %>% grToJosh()
rm(gr_studies,gr_hotspots,nearest_indices)

#########################################################
# JOINTLY REMAP THE HOTSPOT MARKERS FOR GENO & IXN QTLS #
#########################################################

jointMarkerTest(
    studies = studies,
    targetMarkers = portedHotspots[53:77],
    unionBg = unionBg, 
    smetadata = scaMetadata,
    suppGenes = doubleSupp,
    nPerms = 1000,
    FDR = 0.05,
    fn = "joint_hotspot_QTLs",
    correctBG = TRUE,
    correctGRM = FALSE,
    cores = 25,
    chromosomes = NULL
)