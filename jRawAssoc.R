#############################
#############################
## JOINT DATASET REMAPPING ##
## UNFILTERED ASSOCIATIONS ##
## ALL MARKERS             ##
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

###########################################
# CRUDELY ASSOCIATE ALL GENE-MARKER PAIRS #
###########################################

# Detect number of cores and run
rawAssocs(
    studies = studies,
    targetMarkers = colnames(studies$Brem_2002$gdata)[startIdx:endIdx],
    unionBg = unionBg,
    smetadata = scaMetadata,
    suppGenes = doubleSupp,
    fn = 'unfiltered_assocs',
    correctBG = TRUE,
    correctGRM = FALSE,
    cores = detectCores(),
    chromosomes = NULL
)