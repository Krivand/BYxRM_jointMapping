#################################
#################################
## JOINT ASSOCIATION ANALYSIS  ##
## HOTSPOT MARKERS             ##
#################################
#################################


############################################
# SOURCE IMPORTANT FUNCTIONS AND CONSTANTS #
############################################

source(paste0("/home/albertf/vandy147/p05/",
              "remap_functions.R"))

###################
# LOAD IN OBJECTS #
###################

# Get the script set type from command-line arguments
args = commandArgs(trailingOnly = TRUE)

# Check if scriptSet argument is provided
if (length(args) == 0) {
  stop("No scriptSet argument provided. Please specify 'analyte' or 'medium'.")
}

# get the script set type
scriptSet = args[1]

# Check if scriptSet is valid
if(!(scriptSet %in% c('analyte', 'medium'))){
  stop("Invalid scriptSet option. Please choose either 'analyte' or 'medium'.")
}

# load QTLs and raw assocations
if(scriptSet == 'analyte'){
  
  # QTLs
  load(paste0(imed,"jahQTLs.rda"))
  jointHotspotQTLs = jahQTLs
  rm(jahQTLs)
  
  # raw assocations
  load(paste0(imed,"raAnalyte.rda"))
  raAssocs = raAnalyte
  rm(raAnalyte)
  
} else {
  
  # QTLs
  load(paste0(imed, "jmhQTLs.rda"))
  jointHotspotQTLs = jmhQTLs %>%
    filter(!grepl("^chr(II|III|V|XIV):", marker))
  rm(jmhQTLs)
  
  # Raw associations
  load(paste0(imed, "raMedium.rda"))
  raAssocs = raMedium %>%
    filter(!grepl("^chr(II|III|V|XIV):", geno.marker))
  rm(raMedium)
  
}

# drop QTLs passing 5% FDR but having a weak nominal p-value, then format nicely
jointHotspotQTLs = jointHotspotQTLs %>%
  filter(p_value < 0.05) %>%
  mutate(pairID = paste(gene, marker, sep = "_")) %>%
  select(pairID, gene, marker, flavor, 
         samplesize, p_value, mlog10p, effect) %>%
  mutate(flavor = case_when(
    scriptSet == "analyte" & flavor == "analyte:geno" ~ "ixn",
    scriptSet == "medium" & flavor == "medium:geno" ~ "ixn",
    TRUE ~ flavor
  ))
ensure_directory(paste0(asys,"info_",scriptSet,"/"))
write.csv(jointHotspotQTLs,row.names=FALSE,file=paste0(
  asys,"info_",scriptSet,"/","joint_hotspot_QTLs_",scriptSet,".csv"))

# remove duplicate marker and gene columns
raAssocs = raAssocs %>% 
  rename(marker = geno.marker, gene = geno.gene) %>%
  select(-matches("\\.marker$"), 
         -matches("\\.gene$"), 
         -matches("\\.flavor$")) %>%
  mutate(pairID = paste(gene, marker, sep = "_")) %>% 
  select(pairID, gene, marker, everything()) %>%
  rename_with(
    .fn = ~ sub("\\.effect$", ".allelic_effect", .x),
    .cols = ends_with(".effect")
  )

# read in the gene metadata table from yeastmine
geneInfo = read.csv(paste0(
  raw,file='yeastmine_results_2024-06-18T13-15-41.csv'))

# List all files in the directory that start with "joint_QTLs_"
files = list.files(imed, pattern = "^joint_QTLs_.*\\.csv$", full.names = TRUE)

# GET THE STUDY OBJECT AND EXAMINE GENE SUPPORT
studies = lapply(studyList, function(x) {
  load(paste0(imed, x, '.rda'))
  get(x)
})
names(studies) = studyList

# get studies objects without the awkward AlbertAlt naming convention
studies = studies %>%
  {.[!names(.) %in% "Albert_2018"]} %>% 
  {setNames(., str_replace_all(
    names(.),"AlbertAlt_2018", "Albert_2018"))} %>% 
  {setNames(., str_replace_all(
    names(.),"SmithEth_2008", "Smith_2008_Ethanol"))} %>% 
  {setNames(., str_replace_all(
    names(.),"SmithGlc_2008", "Smith_2008_Glucose"))} %>% 
  {setNames(., str_replace_all(
    names(.),"GrossbachRNA_2022", "Grossbach_2022_mRNA"))} %>% 
  {setNames(., str_replace_all(
    names(.),"GrossbachPRT_2022", "Grossbach_2022_protein"))} %>% 
  {setNames(., str_replace_all(names(.), "_", " "))} 

# load the study metadata object and add gene support
studyMetadata = read.delim(file=paste0(raw,'studyMetadata.csv'),
                           header=TRUE,sep=',') %>%
  as.data.frame()

studyMetadata = studyMetadata %>%
  filter(study != "Albert_2018") %>%
  mutate(study = str_replace_all(
    study, "AlbertAlt_2018", "Albert_2018")) %>%
  mutate(study = str_replace_all(
    study, "SmithEth_2008", "Smith_2008_Ethanol")) %>%
  mutate(study = str_replace_all(
    study, "SmithGlc_2008", "Smith_2008_Glucose")) %>%
  mutate(study = str_replace_all(
    study, "GrossbachRNA_2022", "Grossbach_2022_mRNA")) %>%
  mutate(study = str_replace_all(
    study, "GrossbachPRT_2022", "Grossbach_2022_protein")) %>%
  mutate(study = str_replace_all(
    study, "_", " "))

# remove bulk segregant analyses
studyMetadata = studyMetadata[studyMetadata$mapping_approach == 'QTL',]

# add gene support to study metadata
studyMetadata$genes = lapply(studyMetadata$study,
                             function(x)colnames(studies[[x]]$pheno))

# add segregant support to medata
studyMetadata$segs = lapply(studyMetadata$study,
                            function(x)rownames(studies[[x]]$pheno))


if(scriptSet == 'analyte'){
  
  # metadata in SC using genotype A (not Albert 2018)
  subsetMeta = studyMetadata %>%
    dplyr::filter(medium == "SC", genetic_background == "A")
  
} else {
  
  # metadata must be mRNA only and exclude YNB2E
  subsetMeta = studyMetadata %>%
    filter(medium %in% c('SC', 'YNB2G')) %>%
    filter(analyte == 'mRNA')
  
}

# find which genes are supported by at least two datasets
doubleSupp = suppTwo(subsetMeta,type=scriptSet)
print(paste('Within the',scriptSet,"analyses there are",length(doubleSupp),
            'genes supported by two different studies'))
write.table(doubleSupp, file = 
              paste0(asys,"info_",scriptSet,"/",
                     "supported_genes_",scriptSet,".txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.csv(subsetMeta %>% select(-genes, -segs), file = paste0(
  asys,"info_",scriptSet,"/","supported_studies_",scriptSet,".csv"), 
  row.names = FALSE, quote = FALSE)

# hotspot data
hotspots = readxl::read_excel(
  paste0(raw,'SI_Data_08_hotspotTableForPaper_topGenes100_170420.xlsx')) %>%
  as.data.frame()

# prepare to add new hand-curated hotspots
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
# removes chromsomes with engineering diffs or mutations in the medium sets
gr_studies = markRanger(colnames(studies[[1]]$gdata))
gr_hotspots = markRanger(hotspots$hotspotMarker[
  hotspots$hotspotMarker != 'chrXIV:372376_G/A']) %>%
  (function(x)c(x,newHotspots)) %>% sort()
mcols(gr_hotspots)$josh = grToJosh(gr_hotspots)
nearest_indices = nearest(gr_hotspots, gr_studies)
portedHotspots = gr_studies[nearest_indices] %>% 
  grToJosh() %>%
  {if (scriptSet == "medium") {
    magrittr::extract(., !grepl("^chr(II|III|V|XIV):", .))} else {.}}
rm(gr_studies,gr_hotspots,nearest_indices)

# GET GRANGES THAT HOLD THE CLOSET (ANY) MARKER TO EACH GENE

# read in the gene metadata table from yeastmine
geneInfo = read.csv(paste0(
  raw,file='yeastmine_results_2024-06-18T13-15-41.csv'))

# ensure that all supported genes are in the yeastmine metadata
all(doubleSupp %in% geneInfo$Gene...Systematic.Name)

# subset to information about the closest genes
suppInfo = geneInfo[geneInfo$Gene...Systematic.Name %in% doubleSupp,]

# convert the list to a data frame
dfSupp = list(
  suppInfo$Gene...Systematic.Name,
  suppInfo$Gene...Chromosome...Identifier,
  suppInfo$Gene...Chromosome.Location...Start,
  suppInfo$Gene...Chromosome.Location...End
) %>% 
  as.data.frame() %>%
  `names<-`(c('gene', 'chr', 'start', 'end'))
rm(suppInfo)

# create GRanges objects
grSupp = GRanges(seqnames = dfSupp$chr,
                 ranges = IRanges(start = dfSupp$start, end = dfSupp$end),
                 gene = dfSupp$gene)

grAllMarkers = markRanger(unique(raAssocs$marker)) %>% sort()

# Find the nearest markers and create the closest genes and markers dataframe
nearest_indices = nearest(grSupp, grAllMarkers)
closeMarkers = data.frame(
  gene = mcols(grSupp)$gene[!is.na(nearest_indices)],
  marker = mcols(grAllMarkers[nearest_indices[!is.na(nearest_indices)]])$josh
)
rm(nearest_indices)

# get hotspot objects with pval adjustments
unfilteredHotspots = raAssocs[
  raAssocs$marker %in% portedHotspots,] %>%
  mutate(across(
    ends_with(".allelic_effect"),
    ~ abs(.),
    .names = "{gsub('\\\\.allelic_effect$', '', .col)}.absolute_allelic_effect"
  ))

##############
##############
## ANALYSIS ##
##############
##############

###################################################################
# COMPARE THE PERMUTATION-BASED FDR APPROACH WITH PVAL ADJUSTMENT #
###################################################################

# get the name and build the subfolder
subsetDir = paste0(asys,'association_summaries_',scriptSet,'/')
ensure_directory(subsetDir)

# do BY-adjustment on p-values to discover hotpsots and then save
testFrame = data.frame(
  pairID = rep(unfilteredHotspots$pairID, 2),
  flavor = rep(c("genotype", "interaction"), each = nrow(unfilteredHotspots)),
  pval = c(unfilteredHotspots$geno.pval, unfilteredHotspots$ixn.pval)
) %>% mutate(pBY = p.adjust(pval, method = "BY"))
qHotspots = filter(testFrame, pBY < 0.05)
qHotspots %>%
  separate(pairID, into = c("gene", "marker"), sep = "_", remove = FALSE) %>%
  write.csv(,file=paste0(
    subsetDir,'BY_adjusted_significant_hotspots.csv'),
    row.names = FALSE)

# compare permutational discovery to p-value adjustment discovery
compareMethodsFDR = as.data.frame(matrix(NA,ncol=2,nrow=3)) %>%
  `colnames<-`(c('genotype','interaction')) %>%
  `rownames<-`(c('permutation','BY-adjustment','intersection'))

compareMethodsFDR['permutation', 'genotype'] =
  length(jointHotspotQTLs$pairID[jointHotspotQTLs$flavor == 'geno'])
compareMethodsFDR['BY-adjustment', 'genotype'] =
  length(qHotspots$pairID[qHotspots$flavor == 'genotype'])
compareMethodsFDR['intersection', 'genotype'] =
  intersect(jointHotspotQTLs$pairID[jointHotspotQTLs$flavor == 'geno'],
            qHotspots$pairID[qHotspots$flavor == 'genotype']) %>% length()

compareMethodsFDR['permutation', 'interaction'] =
  length(jointHotspotQTLs$pairID[jointHotspotQTLs$flavor == 'ixn'])
compareMethodsFDR['BY-adjustment', 'interaction'] =
  length(qHotspots$pairID[qHotspots$flavor == 'interaction'])
compareMethodsFDR['intersection', 'interaction'] =
  intersect(jointHotspotQTLs$pairID[jointHotspotQTLs$flavor == 'ixn'],
            qHotspots$pairID[qHotspots$flavor == 'interaction']) %>% length()

write.csv(compareMethodsFDR,file=paste0(
  subsetDir,"QTL_discovery_method_comparison.csv"))

rm(compareMethodsFDR, qHotspots, testFrame)

#####################
# QTL ONLY ANALYSES #
#####################

top10ixn = jointHotspotQTLs %>%
  filter(flavor == 'ixn') %>% select(-pairID, -flavor) %>%
  arrange(desc(mlog10p)) %>% (function(x)x[1:10,]) %>%
  rename(`allelic effect` = effect) %>% rename(`p value` = p_value) %>%
  rename(`minus log10(p)` = mlog10p)

write.csv(top10ixn,paste0(
  subsetDir,'top10_interaction_QTLs.csv'),row.names = FALSE)

top10geno = jointHotspotQTLs %>%
  filter(flavor == 'geno') %>% select(-pairID, -flavor) %>%
  arrange(desc(mlog10p)) %>% (function(x)x[1:10,]) %>%
  rename(`allelic effect` = effect) %>% rename(`p value` = p_value) %>%
  rename(`minus log10(p)` = mlog10p)

write.csv(top10geno,paste0(
  subsetDir,'top10_genotype_QTLs.csv'),row.names = FALSE)

# create a barplot for the number of ixn QTLs discovered
p = jointHotspotQTLs %>%
  filter(flavor == 'ixn') %>%
  count(marker) %>%
  mutate(marker = gsub("chr", "", marker)) %>%
  mutate(marker = factor(marker, levels = rev(gsub("chr", "", marker)))) %>%
  ggplot(aes(y = marker, x = n)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
  geom_text(aes(label = n),
            color = "black",
            size = 2.5,
            fontface = "bold",
            hjust = -0.2) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  theme_bw() +
  scale_x_continuous() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black", size = 8),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black")) +
  labs(x = "Number of interaction QTLs observed", y = "Marker")

# Save the plot with adjusted dimensions
ggsave(paste0(subsetDir,"interaction_QTL_barplot.pdf"),
       plot = p, width = 7.5, height = 10)


rm(top10geno,top10ixn,p)

############################################
# EXPLORE DIFFERENCES IN INTERACTION CLASS #
############################################

subsetDir = paste0(asys,"difference_summaries_",scriptSet,"/")

ixnPairQTLs = jointHotspotQTLs$pairID[jointHotspotQTLs$flavor == 'ixn']
genoPairQTLs = jointHotspotQTLs$pairID[jointHotspotQTLs$flavor == 'geno']

# split unfiltered associations into those that are ixn and geno QTLs
isIxnQTL = raAssocs[raAssocs$pairID %in% ixnPairQTLs,]
isGenoQTL = raAssocs[raAssocs$pairID %in% genoPairQTLs,]

# get summaries of mRNA and protein differences using different slices
list(
  raAssocs = raAssocs,
  isIxnQTL = isIxnQTL,
  isGenoQTL = isGenoQTL
) %>% `names<-`(c("All Associations", "Interaction QTLs", "Genotype QTLs")) %>%
  mpSummarizer(directory=subsetDir,mode=scriptSet)

# show differences within QTL classes using funfetti plots
mpLinkPlotter(isGenoQTL, paste0(
  subsetDir, "interaction_class_differences_genotype_QTLs.pdf"),
  "geno", mode = scriptSet)
mpLinkPlotter(isIxnQTL, paste0(
  subsetDir, "interaction_class_differences_interaction_QTLs.pdf"),
  "ixn", mode = scriptSet)

# understand sign-flipping
sign_flips_table = data.frame(
  Category = c("Sign Flips", "Total", "Fraction"),
  `Genotype QTLs` = c(
    sum((isGenoQTL$mRNA.allelic_effect *
           isGenoQTL$protein.allelic_effect) < 0),
    nrow(isGenoQTL),
    sum((isGenoQTL$mRNA.allelic_effect *
           isGenoQTL$protein.allelic_effect) < 0) / nrow(isGenoQTL)
  ),
  `Interaction QTLs` = c(
    sum((isIxnQTL$mRNA.allelic_effect *
           isIxnQTL$protein.allelic_effect) < 0),
    nrow(isIxnQTL),
    sum((isIxnQTL$mRNA.allelic_effect *
           isIxnQTL$protein.allelic_effect) < 0) / nrow(isIxnQTL)
  )
) %>% write.csv(
  file=paste0(subsetDir,"sign_flip_counts.csv"),
  row.names = FALSE)

#################################
# CLOSEST MARKER (CIS) ANALYSIS #
#################################

# Define directories based on scriptSet
subsetDir = paste0(asys, "proximal_analyses_", scriptSet, "/")
ensure_directory(subsetDir)

# Define variables to swap between mRNA/protein and SC/YNB2G
if (scriptSet == "analyte") {
  qval1 = "mRNA.qval"
  qval2 = "protein.qval"
  effect1 = "mRNA.allelic_effect"
  effect2 = "protein.allelic_effect"
  flavor1 = "mRNA"
  flavor2 = "protein"
  pval1 = "mRNA.pval"
  pval2 = "protein.pval"
} else {
  qval1 = "SC.qval"
  qval2 = "YNB2G.qval"
  effect1 = "SC.allelic_effect"
  effect2 = "YNB2G.allelic_effect"
  flavor1 = "SC"
  flavor2 = "YNB2G"
  pval1 = "SC.pval"
  pval2 = "YNB2G.pval"
}

# Subset the raAssocs dataframe by creating a pairID of close pairs
closeAssocs = raAssocs[
  raAssocs$pairID %in% paste(closeMarkers$gene, closeMarkers$marker, sep = "_"), 
] %>%
  select(-matches("geno\\."))

# Get BY-adjusted qvals
closeAssocs = closeAssocs %>%
  mutate(
    ixn.qval = p.adjust(ixn.pval, method = "BY"),
    !!qval1 := p.adjust(get(pval1), method = "BY"),
    !!qval2 := p.adjust(get(pval2), method = "BY")
  )

# Filter and combine the data into one dataframe
qNearQTLs = list(
  closeAssocs %>% filter(ixn.qval < 0.05) %>% mutate(flavor = "ixn"),
  closeAssocs %>% filter(!!sym(qval1) < 0.05) %>% mutate(flavor = flavor1),
  closeAssocs %>% filter(!!sym(qval2) < 0.05) %>% mutate(flavor = flavor2)
) %>%
  bind_rows() %>%
  select(-pairID) %>%
  select(
    gene, marker,
    ixn.allelic_effect, !!sym(effect1), !!sym(effect2),
    ixn.qval, !!sym(qval1), !!sym(qval2),
    flavor
  ) %>%
  rename(
    effect_ixn = ixn.allelic_effect,
    effect1 = !!sym(effect1),
    effect2 = !!sym(effect2),
    qval_ixn = ixn.qval,
    qval1 = !!sym(qval1),
    qval2 = !!sym(qval2)
  )

# Replace variable names in rename function
names(qNearQTLs) = sub("effect1", paste0("effect_", flavor1), names(qNearQTLs))
names(qNearQTLs) = sub("effect2", paste0("effect_", flavor2), names(qNearQTLs))
names(qNearQTLs) = sub("qval1", paste0("qval_", flavor1), names(qNearQTLs))
names(qNearQTLs) = sub("qval2", paste0("qval_", flavor2), names(qNearQTLs))

write.csv(qNearQTLs,row.names = FALSE,
          file=paste0(subsetDir,'nearest_gene_marker_QTLs.csv'))

padists = c(0, 5000, 10000, 100000, 1000000)

# Loop over each flavor: "ixn" and "geno"
for (flavor_type in c("ixn", "geno")) {
  
  # Filter jointHotspotQTLs to only include genes of the current flavor
  filteredHotspotQTLs = jointHotspotQTLs %>% filter(flavor == flavor_type)
  grHotspotQTLs = markRanger(filteredHotspotQTLs$marker)
  mcols(grHotspotQTLs)$gene = filteredHotspotQTLs$gene
  seqlevels(grHotspotQTLs) = seqlevels(grSupp)
  
  # Create the GRanges object for genes in the same order as filteredHotspotQTLs
  grGeneQTLs = grSupp[match(filteredHotspotQTLs$gene, mcols(grSupp)$gene)]
  
  # Filter out pairs that are not on the same chromosome
  same_chr_indices = as.logical(seqnames(grHotspotQTLs) == seqnames(grGeneQTLs))
  grHotspotQTLs = grHotspotQTLs[same_chr_indices]
  filteredHotspotQTLs = filteredHotspotQTLs[same_chr_indices,]
  grGeneQTLs = grGeneQTLs[same_chr_indices]
  
  # Get the number of unique QTLs within each padding distance
  padists %>%
    sapply(function(x) {
      hits = padded_pintersect(grHotspotQTLs, grGeneQTLs, x)
      sum(hits)
    }) %>%
    enframe(name = "padding_distance", value = "cis_QTLs") %>%
    mutate(padding_distance = padists) %>%
    write.csv(file = paste0(subsetDir, "cisQTLs_by_padding_distance_", 
                            flavor_type, ".csv"),
              row.names = FALSE)
  
  filteredHotspotQTLs[
    padded_pintersect(grHotspotQTLs, grGeneQTLs, 5000), ] %>% 
    select(-pairID) %>% 
    write.csv(file = paste0(subsetDir, 'cis_QTLs_by_5K_padding_',
                            flavor_type, '.csv'), row.names = FALSE)
  
}

rm(filteredHotspotQTLs, grHotspotQTLs, same_chr_indices, padists)

########################
# GENETIC CORRELATIONS #
########################

# Define directories based on scriptSet
subsetDir = paste0(asys, "genetic_correlations_", scriptSet, "/")
ensure_directory(subsetDir)

# Define variables for allelic effect and p-values based on scriptSet
if (scriptSet == "analyte") {
  effect1 = "mRNA.allelic_effect"
  effect2 = "protein.allelic_effect"
  pval1 = "mRNA.pval"
  pval2 = "protein.pval"
} else {
  effect1 = "SC.allelic_effect"
  effect2 = "YNB2G.allelic_effect"
  pval1 = "SC.pval"
  pval2 = "YNB2G.pval"
}

# Calculate genetic correlations for allelic effects
hCorEffect = lapply(sort(unique(unfilteredHotspots$gene)), function(gene) {
  cor.test(
    unfilteredHotspots[unfilteredHotspots$gene == gene, ][[effect1]],
    unfilteredHotspots[unfilteredHotspots$gene == gene, ][[effect2]],
    method = 'pearson'
  )
}) %>% `names<-`(sort(unique(unfilteredHotspots$gene)))

# Calculate genetic correlations for -log10(p-values)
hCorLogP = lapply(sort(unique(unfilteredHotspots$gene)), function(gene) {
  cor.test(
    -log10(unfilteredHotspots[unfilteredHotspots$gene == gene, ][[pval1]]),
    -log10(unfilteredHotspots[unfilteredHotspots$gene == gene, ][[pval2]]),
    method = 'pearson'
  )
}) %>% `names<-`(sort(unique(unfilteredHotspots$gene)))

# Process effect estimates
hCE = list(
  p = sapply(hCorEffect, function(x) x$p.value),
  r = sapply(hCorEffect, function(x) x$estimate)
) %>% as.data.frame() %>%
  mutate(
    pBY = p.adjust(p, method = "BY"),
    value = "effect_estimate"
  ) %>%
  rownames_to_column(var = "genes")

# Process -log10(p-values)
hCP = list(
  p = sapply(hCorLogP, function(x) x$p.value),
  r = sapply(hCorLogP, function(x) x$estimate)
) %>% as.data.frame() %>%
  mutate(
    pBY = p.adjust(p, method = "BY"),
    value = "minus_log10_p"
  ) %>%
  rownames_to_column(var = "genes")

# Combine and sort the results
hA = rbind(hCE, hCP) %>% arrange(genes)

# Filter and save results for significant negative correlations (if any)

# Effect values
hA %>% 
  filter(value == "effect_estimate", r < 0, pBY < 0.05) %>% 
  (function(df) {
    write.csv(
      df, 
      paste0(subsetDir, 'significant_effect_anticorrelated_genes.csv'), 
      row.names = FALSE
    )
  })

# Locus significance
hA %>% 
  filter(value == "minus_log10_p", r < 0, pBY < 0.05) %>% 
  (function(df) {
    write.csv(
      df, 
      paste0(subsetDir, 'significant_basis_anticorrelated_genes.csv'), 
      row.names = FALSE
    )
  })

# Save full results to CSV
write.csv(hA, paste0(subsetDir, 'genetic_correlations_full.csv'),
          row.names = FALSE)

# Calculate and save summary statistics
hA %>%
  group_by(value) %>%
  summarize(
    mean_r = mean(r),
    median_r = median(r),
    sd_r = sd(r),
    min_r = min(r),
    max_r = max(r),
    mean_p = mean(p),
    median_p = median(p),
    sd_p = sd(p),
    min_p = min(p),
    max_p = max(p)
  ) %>% 
  write.csv(row.names = FALSE, file = 
              paste0(subsetDir, 'genetic_correlations_summaries.csv'))

# Create and save histogram plot
p = ggplot(hA, aes(x = r, fill = value)) +
  geom_histogram(
    aes(y = after_stat(density)), 
    position = "identity", 
    bins = 30, 
    color = "black",
    alpha = 0.5
  ) +
  scale_fill_viridis(
    discrete = TRUE,
    labels = c("Effect estimate", bquote(-log[10](italic(p))))
  ) +
  theme_bw() +
  theme(
    text = element_text(color = "black"), 
    axis.title = element_text(color = "black"), 
    axis.text = element_text(color = "black"), 
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(color = "black") 
  ) +
  labs(
    title = bquote(
      "Genetic correlations of effect estimates and" 
      ~ -log[10] * (italic(p)) ~ "values"),
    x = "R-value",
    y = "Density"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave(paste0(subsetDir, "genetic_correlations_histogram.pdf"),
       plot = p, width = 7.5, height = 5)

# Clean up variables
rm(hCorEffect, hCorLogP, hCE, hCP, hA, p)