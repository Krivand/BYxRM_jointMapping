############################################
# SOURCE IMPORTANT FUNCTIONS AND CONSTANTS #
############################################

source(paste0("/home/albertf/vandy147/p05/",
              "remap_functions.R"))

###################
# LOAD IN OBJECTS #
###################

# LOAD IN UNFILTERED ASSOCATIONS AND PREPARE
load(paste0(imed,"unfilteredAssocs.rda"))

# Remove duplicate marker and gene columns
unfilteredAssocs = unfilteredAssocs %>% 
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

# load in joint QTL objects
load(paste0(imed,"jhQTLs.rda"))

jhQTLs = jhQTLs %>% mutate(pairID = paste(gene, marker, sep = "_"))

jointHotspotQTLs = jhQTLs %>%
  filter(mlog10p_threshold > 0) %>%
  rename(threshold = mlog10p_threshold) %>%
  select(pairID, gene, marker, flavor, samplesize, p_value,
         mlog10p, effect, absolute_effect, threshold) %>%
  mutate(flavor = ifelse(flavor == "analyte:geno", "ixn", flavor))

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


# scaMetadata is metadata in SC using genotype A (not Albert 2018)
scaMetadata = studyMetadata %>%
  dplyr::filter(medium == "SC", genetic_background == "A")

# within the synthetic complete media (sc) datasets
# find which genes are supported by at least two datasets 
doubleSupp = suppTwo(scaMetadata)
print(paste0('Within the SC media with genotype A subset there are ',
             length(doubleSupp),
             ' genes supported by two mRNA and two protein datasets'))

# hotspot data
hotspots = readxl::read_excel(
  paste0(raw,'SI_Data_08_hotspotTableForPaper_topGenes100_170420.xlsx')) %>%
  as.data.frame()

# prepare to add new hand-curated hotspots
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
gr_studies = markRanger(colnames(studies[[1]]$gdata))
gr_hotspots = markRanger(hotspots$hotspotMarker[
  hotspots$hotspotMarker != 'chrXIV:372376_G/A']) %>%
  (function(x)c(x,newHotspots)) %>% sort()
mcols(gr_hotspots)$josh = grToJosh(gr_hotspots)
nearest_indices = nearest(gr_hotspots, gr_studies)
portedHotspots = gr_studies[nearest_indices] %>% grToJosh()
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
dfSupp = list(suppInfo$Gene...Systematic.Name,
               suppInfo$Gene...Chromosome...Identifier,
               suppInfo$Gene...Chromosome.Location...Start,
               suppInfo$Gene...Chromosome.Location...End) %>% 
  as.data.frame() %>%
  `names<-`(c('gene','chr','start','end'))
rm(suppInfo)


# create GRanges objects
grSupp = GRanges(seqnames = dfSupp$chr,
                   ranges = IRanges(start = dfSupp$start, end = dfSupp$end),
                   gene = dfSupp$gene)

grAllMarkers = markRanger(unique(unfilteredAssocs$marker)) %>% sort()

# Find the nearest markers
nearest_indices = nearest(grSupp, grAllMarkers)

# Create the dataframe with genes and closest markers
closeMarkers = data.frame(
  gene = mcols(grSupp)$gene,
  marker = mcols(grAllMarkers[nearest_indices])$josh
)

rm(nearest_indices)

# get hotspot objects with pval adjustments
unfilteredHotspots = unfilteredAssocs[
  unfilteredAssocs$marker %in% portedHotspots,] %>%
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

#########################################################
# MAKE SOME CRITICAL SUMMARY TABLES FROM LOADED OBJECTS #
#########################################################

# create a useful subfolder
ensure_directory(paste0(asys,'association_summaries'))

# get gene pair associations saved for the permutation discovered hotspots
file_path = paste0(asys, "association_summaries/joint_hotspot_QTLs.csv")
if (!file.exists(file_path)) {
  write.csv(jointHotspotQTLs, file = file_path, row.names = FALSE)
} else {
  message("File already exists: ", file_path)
}

# get gene pair associations saved for all gene-marker associations
file_path = paste0(
  asys, 'association_summaries/all_gene_marker_assocations.csv')
if (!file.exists(file_path)) {
  unfilteredAssocs %>%
    select(-pairID) %>%
    write.csv(file = file_path, row.names = FALSE)
} else {
  message("File already exists: ", file_path)
}

# get gene pair associations saved for all gene-hotspot associations
file_path = paste0(
  asys, 'association_summaries/all_gene_hotspots_assocations.csv')
if (!file.exists(file_path)) {
unfilteredHotspots = unfilteredAssocs[
  unfilteredAssocs$marker %in% portedHotspots,] %>% select(-pairID) %>%
  write.csv(file = file_path, row.names = FALSE)
} else {
  message("File already exists: ", file_path)
}

# remove file path and move on
rm(file_path)

###################################################################
# COMPARE THE PERMUTATION-BASED FDR APPROACH WITH PVAL ADJUSTMENT #
###################################################################

# do BY-adjustment on p-values to discover hotpsots and then save
testFrame = data.frame(
  pairID = rep(unfilteredHotspots$pairID, 2),
  flavor = rep(c("genotype", "interaction"), each = nrow(unfilteredHotspots)),
  pval = c(unfilteredHotspots$geno.pval, unfilteredHotspots$ixn.pval)
) %>% mutate(pBY = p.adjust(pval, method = "BY"))
qHotspots = filter(testFrame, pBY < 0.05)
qHotspots %>%
  separate(pairID, into = c("gene", "marker"), sep = "_", remove = FALSE) %>%
  write.csv(,file=paste0(asys,'association_summaries/',
                         'BY_adjusted_significant_hotspots.csv'),
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

write.csv(compareMethodsFDR,
          file=paste0(asys,'association_summaries/',
                      "QTL_discovery_method_comparison.csv"))

rm(compareMethodsFDR, qHotspots, testFrame)

#####################
# QTL ONLY ANALYSES #
#####################

top10ixn = jointHotspotQTLs %>%
  filter(flavor == 'ixn') %>% select(-pairID, -flavor, -absolute_effect) %>%
  arrange(desc(mlog10p)) %>% (function(x)x[1:10,]) %>% 
  rename(`allelic effect` = effect) %>% rename(`p value` = p_value) %>%
  rename(`minus log10(p)` = mlog10p)

write.csv(top10ixn,paste0(asys,'association_summaries/',
                          'top10_interaction_QTLs.csv'),row.names = FALSE)

top10geno = jointHotspotQTLs %>%
  filter(flavor == 'geno') %>% select(-pairID, -flavor, -absolute_effect) %>%
  arrange(desc(mlog10p)) %>% (function(x)x[1:10,]) %>%
  rename(`allelic effect` = effect) %>% rename(`p value` = p_value) %>%
  rename(`minus log10(p)` = mlog10p)

write.csv(top10geno,paste0(asys,'association_summaries/',
                           'top10_genotype_QTLs.csv'),row.names = FALSE)

# create a disjoint barplot for the number of ixn QTLs discovered
p = table(jointHotspotQTLs$marker[jointHotspotQTLs$flavor == 'ixn']) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  rename(marker = Var1, count = Freq) %>%
  mutate(marker = gsub("chr", "", marker)) %>%
  mutate(marker = factor(marker, levels = gsub("chr", "", names(
    table(jointHotspotQTLs$marker[jointHotspotQTLs$flavor == 'ixn']))))) %>%
  ggplot(aes(x = marker, y = count)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
  geom_text(aes(label = count), 
            vjust = -0.2, 
            color = "black", 
            size = 4, 
            fontface = "bold") +  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) + 
  theme_bw() +
  scale_y_continuous(breaks = c(seq(0, 25, by = 5), 134)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black")) +
  labs(x = "Marker", y = "Number of interaction QTLs observed")

ggsave(paste0(asys,'association_summaries/',
              "interaction_QTL_barplot.pdf"), 
       plot = p, width = 7.5, height = 5)

rm(top10geno,top10ixn,p)

########################################
# EXPLORE PROTEIN AND MRNA DIFFERENCES #
########################################

ixnPairQTLs = jointHotspotQTLs$pairID[jointHotspotQTLs$flavor == 'ixn']
genoPairQTLs = jointHotspotQTLs$pairID[jointHotspotQTLs$flavor == 'geno']

# split unfiltered associations into those that are ixn and geno QTLs
isIxnQTL = unfilteredAssocs[unfilteredAssocs$pairID %in% ixnPairQTLs,] 
isGenoQTL = unfilteredAssocs[unfilteredAssocs$pairID %in% genoPairQTLs,] 

# get summaries of mRNA and protein differences using different slices
list(
  unfilteredAssocs = unfilteredAssocs,
  isIxnQTL = isIxnQTL,
  isGenoQTL = isGenoQTL
) %>% `names<-`(c("All Associations", "Interaction QTLs", "Genotype QTLs")) %>% 
  mpSummarizer(directory=paste0(asys,'mRNA_vs_protein_difference_summaries/'))

# show differences within QTL classes using funfetti plots
mpLinkPlotter(isGenoQTL, paste0(asys,"mRNA_vs_protein_difference_summaries/",
                                "mRNA_vs_protein_genotype_QTLs.pdf"),"geno")
mpLinkPlotter(isIxnQTL, paste0(asys,"mRNA_vs_protein_difference_summaries/",
                               "mRNA_vs_protein_iteraction_QTLs.pdf"),"ixn")
mpLinkPlotter(isIxnQTL %>% filter(marker != "chrXIV:465187"),
              paste0(asys,"mRNA_vs_protein_difference_summaries/",
                     "mRNA_vs_protein_iteraction_QTLs_no_MKT1.pdf"),"ixn")

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
  file=paste0(asys,"mRNA_vs_protein_difference_summaries/",
              "sign_flip_counts.csv"),
  row.names = FALSE)

#################################
# CLOSEST MARKER (CIS) ANALYSIS #
#################################

ensure_directory(paste0(asys,'proximal_analyses'))

# Subset the unfilteredAssocs dataframe by creating a pairID of closet pairs
closeAssocs = unfilteredAssocs[
  unfilteredAssocs$pairID %in% 
    paste(closeMarkers$gene,
          closeMarkers$marker, sep = "_"), ] %>%
  select(-matches("geno\\."))

# get BY-adjusted qvals
closeAssocs = closeAssocs %>%
  mutate(
    ixn.qval = p.adjust(ixn.pval, method = "BY"),
    protein.qval = p.adjust(protein.pval, method = "BY"),
    mRNA.qval = p.adjust(mRNA.pval, method = "BY")
  )


# Filter and combine the data into one dataframe
qNearQTLs = list(
  closeAssocs %>% filter(ixn.qval < 0.05) %>% mutate(flavor = "ixn"),
  closeAssocs %>% filter(protein.qval < 0.05) %>% mutate(flavor = "protein"),
  closeAssocs %>% filter(mRNA.qval < 0.05) %>% mutate(flavor = "mRNA")
) %>%
  bind_rows() %>%
  select(-pairID) %>%
  select(
    gene, marker,
    ixn.allelic_effect, protein.allelic_effect, mRNA.allelic_effect,
    ixn.qval, protein.qval, mRNA.qval,
    flavor
  ) %>%
  rename(
    effect_ixn = ixn.allelic_effect,
    effect_protein = protein.allelic_effect,
    effect_mRNA = mRNA.allelic_effect,
    qval_ixn = ixn.qval,
    qval_protein = protein.qval,
    qval_mRNA = mRNA.qval
  )

write.csv(qNearQTLs,row.names = FALSE,
          file=paste0(asys,'proximal_analyses/',
                      'nearest_gene_marker_QTLs.csv'))

  
# GET GRANGES THAT HOLD THE CLOSET (HOTSPOT) MARKER TO EACH GENE
grHotspotQTLs = markRanger(jointHotspotQTLs$marker)

# Order the subset to match the order of genes in jointHotspotQTLs$gene
grGeneQTLs = grSupp[
  match(jointHotspotQTLs$gene[jointHotspotQTLs$gene %in% grSupp$gene],
        grSupp$gene)]

# get number of QTLs within each padding
c(0, 5e3, 1e4, 1e5, 1e6) %>%
  setNames(nm = .) %>%  # Set the names of the vector to the distances
  sapply(function(x) sum(padded_pintersect(grHotspotQTLs, grGeneQTLs, x)$hit)) %>%
  enframe(name = "Distance (bp)", value = "Num proximal QTLs") %>%
  mutate(`Distance (bp)` = as.numeric(`Distance (bp)`)) %>% 
  write.csv(file = paste0(asys,'proximal_analyses/',
                  'cis_QTLs_by_padding_distance.csv'),
    row.names = FALSE)
  
# get QTLs within 5k padding distance
write.csv(jointHotspotQTLs[
  padded_pintersect(grHotspotQTLs,grGeneQTLs,5e3)$hit,] %>%
    select(-pairID, -absolute_effect),
          file = paste0(asys,'proximal_analyses/',
                        'cis_QTLs_by_5k_padding.csv'),
  row.names = FALSE)

########################
# GENETIC CORRELATIONS #
########################

ensure_directory(paste0(asys,'genetic_correlations/'))

#CHECK FOR PER-GENE SYSTEMATIC NEGATIVE CORRELATION ACROSS HOTSPOTS

# CALCULATE CORRELATIONS ACROSS HOTSPOTS PER GENE

hCorEffect = lapply(sort(unique(unfilteredHotspots$gene)), function(gene) {
  cor.test(unfilteredHotspots[
    unfilteredHotspots$gene == gene,]$mRNA.allelic_effect,
           unfilteredHotspots[
             unfilteredHotspots$gene == gene,]$protein.allelic_effect,
           method = 'pearson')
}) %>% `names<-`(sort(unique(unfilteredHotspots$gene)))

hCorLogP = lapply(sort(unique(unfilteredHotspots$gene)), function(gene) {
  cor.test(-log10(unfilteredHotspots[
    unfilteredHotspots$gene == gene,]$mRNA.pval),
           -log10(unfilteredHotspots[
             unfilteredHotspots$gene == gene,]$protein.pval),
           method = 'pearson')
}) %>% `names<-`(sort(unique(unfilteredHotspots$gene)))

# SHAPE THE DATA FOR EASY REPORTING
hCE = list(
  p = sapply(hCorEffect,function(x)x$p.value),
  r = sapply(hCorEffect,function(x)x$estimate)
) %>% as.data.frame() %>%
  mutate(
    pBY = p.adjust(p, method = "BY"),
    value = "effect_estimate"
  ) %>%
  rownames_to_column(var = "genes")

hCP = list(
  p = sapply(hCorLogP,function(x)x$p.value),
  r = sapply(hCorLogP,function(x)x$estimate)
) %>% as.data.frame() %>%
  mutate(
    pBY = p.adjust(p, method = "BY"),
    value = "minus_log10_p"
  ) %>%
  rownames_to_column(var = "genes")

hA = rbind(hCE, hCP) %>% arrange(genes)

print('Are there any genes with significant negative correlatons...')
print('Across hotspots for...')
print('Effect values:')
hA %>% filter(value == "effect_estimate", r < 0, pBY < 0.05)
print('Locus significance:')
hA %>% filter(value == "minus_log_10_p", r < 0, pBY < 0.05)

write.csv(hA,paste0(asys,'genetic_correlations/',
  'genetic_correlations_full.csv'))

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
  write.csv(file=paste0(
    asys,'genetic_correlations/',
    'genetic_correlations_summaries.csv'))

p = ggplot(hA, aes(x = r, fill = value)) +
  geom_histogram(
    aes(y = ..density..), 
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
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  labs(
    title = bquote(
      "Genetic correlations of effect estimates and" 
      ~ -log[10] * (italic(p)) ~ "values"),
    x = "R-value",
    y = "Density"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# Save the plot to a PDF file with specified dimensions
ggsave(paste0(asys,'genetic_correlations/',
              "genetic_correlations_histogram.pdf"),
       plot = p, width = 7.5, height = 5)

rm(hCorEffect,hCorLogP,hCE,hCP,hA,p)

###############
# TRACE PLOTS #
###############

# get signed -log10(p) traces for joint mappings
unfilteredAssocs %>%
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
  tracePlotter(valueType = 'smltp',
               rootDir = paste0(asys,'joint_signed_-log10P_traces/'),
               rootName = "joint_signed_-log10P_trace_for_",
               plotWidth=7.5, plotHeight=2.5,
               r = c("geno", "ixn"), 
               n = c("genotype", "interaction"))

# get allelic effect traces for joint mappings
unfilteredAssocs %>%
  select(-starts_with("geno")) %>%
  rename_with(~ str_replace_all(., "ixn", "interaction")) %>%
  rename_with(~ str_remove_all(., '.allelic_effect')) %>%
  select(-matches("\\.")) %>% 
  tracePlotter(valueType = 'Allelic effect',
               rootDir = paste0(asys,'joint_effect_size_traces/'),
               rootName = "joint_effect_trace_for_",
               plotWidth=7.5, plotHeight=2.5,
               r = c("geno", "ixn"), 
               n = c("genotype", "interaction"))

# prepare metdata for nice trace plotting
exclusionMetadata = studyMetadata %>% mutate(study = gsub("_", " ", study)) %>%
  mutate(analyte = ifelse(
    study %in% c("Smith 2008 Ethanol", "Smith 2008 Glucose", "Albert 2018"),
    "excluded", analyte)) %>%
  mutate(analyte = factor(analyte, levels = c("mRNA", "protein", "excluded")))

# run signed LOD traces on every supported gene
studyTracer(studies = studies, 
            smeta = exclusionMetadata,
            suppGenes = doubleSupp, 
            fn = paste0(asys, 'study_signed_LOD_traces'),
            plotHeight = 5,
            plotWidth = 7.5,
            traceType = "s")

# run signed R-value traces on every supported gene
studyTracer(studies = studies, 
            smeta = exclusionMetadata,
            suppGenes = doubleSupp, 
            fn = paste0(asys, 'study_R-value_traces'),
            plotHeight = 5,
            plotWidth = 7.5,
            traceType = "r")

rm(exclusionMetadata)