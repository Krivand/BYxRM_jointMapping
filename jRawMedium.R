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
studies = lapply(medMeta$study, function(x) {
  load(paste0(imed, x, '.rda'))
  get(x)
}) %>% `names<-`(medMeta$study) 

# retrieve background QTLs across studies
medQTLs = lapply(medMeta$study, function(x) {
  load(paste0(imed, x, '_Background.rda'))
  get(paste0(x,'_Background'))
}) %>% `names<-`(medMeta$study) 

# add gene support to metadata
medMeta$genes = lapply(medMeta$study,function(x)
  colnames(studies[[x]]$pheno))

# add segregant support to metadata
medMeta$segs = lapply(medMeta$study,function(x)
  rownames(studies[[x]]$pheno))

# within the mRNA-only datasets
# find which genes are supported by at least two datasets 
medSupp = suppTwo(medMeta,type = "medium")
print(paste0('Within the RNA data subset there are ',
             length(medSupp),
             ' genes supported by two SC and two YNB2G datasets'))

# take the union of all background QTLs discovered across all datasets
if (!file.exists(paste0(imed, "medBg.rda"))) {
  
  medBg = collapseBgQTLs(medQTLs)
  save(medBg,file=paste0(imed,"medBg.rda"))
  
} else {
  
  print('Union of background QTLs already processed and loading now...')
  load(file=paste0(imed,"medBg.rda"))
  
}

# how many strong QTLs are most genes experiencing from other chroms?
# not worth switching to random effects due to covars approaching sample num
sapply(medBg,function(x)sapply(x,function(y)length(y))) %>% 
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
rawAssocs(
  studies = studies,
  targetMarkers = colnames(studies$Brem_2002$gdata)[startIdx:endIdx],
  unionBg = medBg,
  smetadata = medMeta,
  suppGenes = medSupp,
  fixedEffect = 'medium',
  fn = 'unfiltered_media',
  correctBG = TRUE,
  correctGRM = FALSE,
  cores = detectCores(),
  chromosomes = NULL
)