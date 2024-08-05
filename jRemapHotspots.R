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

# retrieve args
args = commandArgs(trailingOnly = TRUE)
startIdx = as.numeric(args[1])
endIdx = as.numeric(args[2])

# source functions and constants
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
# not worth switching to random effects bc # covars not approaching # samples
sapply(unionBg,function(x)sapply(x,function(y)length(y))) %>% 
  unlist() %>% table()

# add new hand-curated hotspots
geneInfo = read.csv(paste0(
  raw,file='yeastmine_results_2024-06-18T13-15-41.csv'))
newHotspots = geneInfo[which(geneInfo$Gene...Standard.Name %in% 
                               c('URA3','LEU2','AMN1')),] %>%   
  select(
    `Gene...Chromosome...Identifier`,
    `Gene...Chromosome.Location...Start`,
    `Gene...Chromosome.Location...End`
  ) %>% setNames(c('chr','start','end')) %>%
  bind_rows(
    data.frame(
      chr = c("chrX", "chrXV"),
      start = c(295000, 1000000),
      end = c(295000, 1000000)
    )
  ) %>% makeGRangesFromDataFrame

print('Please note the KRE33 hotspot chrXIV:372376_G/A is dropped')  
print('Here are additional new hotspots')
print(newHotspots)

# get the indices of hotspots approximated in Grossbach-coordinates
gr_studies = markRanger(colnames(studies$Brem_2002$gdata))
hotspots = readxl::read_excel(
  paste0(raw,'SI_Data_08_hotspotTableForPaper_topGenes100_170420.xlsx')) %>%
  as.data.frame()
gr_hotspots = markRanger(hotspots$hotspotMarker[
  hotspots$hotspotMarker != 'chrXIV:372376_G/A']) %>%
  (function(x)c(x,newHotspots)) %>% sort()
mcols(gr_hotspots)$josh = grToJosh(gr_hotspots)
nearest_indices = nearest(gr_hotspots, gr_studies)
portedHotspots = gr_studies[nearest_indices] %>% grToJosh()
rm(gr_studies,gr_hotspots,nearest_indices,newHotspots)

#########################################################
# JOINTLY REMAP THE HOTSPOT MARKERS FOR GENO & IXN QTLS #
#########################################################

# Detect number of cores and run
jointMarkerTest(
    studies = studies,
    targetMarkers = portedHotspots[startIdx:endIdx],
    unionBg = unionBg, 
    smetadata = scaMetadata,
    suppGenes = doubleSupp,
    nPerms = 1000,
    FDR = 0.05,
    fn = "joint_hotspot_QTLs",
    correctBG = TRUE,
    correctGRM = FALSE,
    cores = detectCores(),
    chromosomes = NULL
)
