#############################
# SOURCE FXNS AND CONSTANTS #
#############################
source(paste0("/home/albertf/vandy147/p05/",
              "remap_functions.R"))

#########################################
# RETRIEVE AND FORMAT NECESSARY OBJECTS #
#########################################

# RECREATE THE JOINT MAPPING OBJECT
studies = lapply(studyList, function(x) {
  load(paste0(imed, x, '.rda'))
  get(x)
})
names(studies) = studyList

# RETRIEVE BACKGROUND QTL INFORMATION
bgQTLs = lapply(studyList, function(x) {
  load(paste0(imed, x, '_Background.rda'))
  get(paste0(x,'_Background'))
}) %>% `names<-`(studyList) 

# TAILOR THE METADATA OBJECT FOR THIS MAPPING

# remove bulk segregant analyses
studyMetadata = studyMetadata[studyMetadata$mapping_approach == 'QTL',]

# add gene support to study metadata
studyMetadata$genes = lapply(studyMetadata$study,
                             function(x)colnames(studies[[x]]$pheno))

# add segregant support to medata
studyMetadata$segs = lapply(studyMetadata$study,
                            function(x)rownames(studies[[x]]$pheno))


# scaMetadata is metadata in SC using genotype A
scaMetadata = studyMetadata %>%
  dplyr::filter(medium == "SC", genetic_background == "A")

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

# FIND THE NEAREST MARKER TO EVERY GENE

# read in the gene metadata table from yeastmine
geneInfo = read.csv(paste0(
  raw,file='yeastmine_results_2024-06-18T13-15-41.csv'))

# ensure that all supported genes are in the yeastmine metadata
all(doubleSupp %in% geneInfo$Gene...Systematic.Name)

# subset the gene metadata to the relevant elements
subInfo = geneInfo[geneInfo$Gene...Systematic.Name %in% doubleSupp,]

# convert the list to a data frame
gene_df = list(subInfo$Gene...Systematic.Name,
                subInfo$Gene...Chromosome...Identifier,
                subInfo$Gene...Chromosome.Location...Start,
                subInfo$Gene...Chromosome.Location...End) %>% 
  as.data.frame() %>%
  `names<-`(c('gene','chr','start','end'))
rm(subInfo)

# create a GRanges object
gr_genes = GRanges(seqnames = gene_df$chr,
              ranges = IRanges(start = gene_df$start, end = gene_df$end),
              gene = gene_df$gene)
rm(gene_df)

# get the indices of hotspots approximated in Brem-coordinates
gr_studies = markRanger(colnames(studies$Brem_2002$gdata))
nearest_indices = nearest(gr_genes, gr_studies)
nearest_gene_markers = gr_studies[nearest_indices] %>% grToJosh()
rm(gr_studies,gr_genes,nearest_indices)

###############################################################
# ITERATE OVER CLOSEST GENE-MARKER PAIRS FOR THE FOUR FLAVORS #
###############################################################

pairAssocs(studies=studies[names(studies) %in% scaMetadata$study],
             targetMarkers = nearest_gene_markers,
             unionBg = unionBg[scaMetadata$name],
             smetadata = scaMetadata,
             suppGenes = doubleSupp)

##########
## KIND OF UNNECESSARY AND DISCONTINUED AFTER THE ALL-ASSOCIATION PARADIGM... THESE CAN JUST BE PULLED FROM THOSE