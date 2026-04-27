##################################
##################################
## INDIVIDUAL DATASET REMAPPING ##
## MAIN SCRIPT                  ##
## ORIGINAL 2018 PARAMS         ##
##################################
##################################

##########################################
# DIRECTORIES, LIBRARIES, AND SOURCEFILE #
##########################################

# SOURCE FXNS AND CONSTANTS
source("/PHShome/kv1026/p05/remap_functions.R")

###############################################
# LOAD CEREVISIAE GENOMIC FEATURE INFORMATION #
###############################################

# yeast Chromosomal features
# from: "https://downloads.yeastgenome.org
# subdirectory /curation/chromosomal_feature/SGD_features.tab"

load(paste0(raw,'chrFeatures.rda'))

# ensure each gene with a name maps to a feature name
chrFeatures[!is.na(chrFeatures$Standard_Gene_Name),] %>% 
  (function(x)x$Feature_Name) %>% 
  (function(x)all(!is.na(x)))

###############################
# PREPARE GROSSBACH GENOTYPES #
###############################

# build a directory for genoT quality assurance
ensure_directory(paste0(asys,'genotype_QA/'))

# GROSSBACH 2022 DEDUPLICATED GENOTYPE MATRIX AND INTERVALS

# load table 
genoG = readxl::read_excel(
  paste0(raw,'msb202110712-sup-0004-datasetev2.xlsx'),
  sheet = 'mappingGenotype') %>% as.data.frame()

# create interval object
iG = genoG[,1:3] %>% dplyr::rename(chr = CHR)

# remove intervals from genotype matrix
genoG = genoG[,-c(1:3)]

# format genoG
genoG = genoG %>%
  (function(m)apply(m,2,function(x)as.numeric(x))) %>%
  as.matrix %>% t() %>% (function(x)ifelse(x == 0, -1, x))
colnames(genoG) = mapply(function(x,y)paste0(x,":",y),
                         iG$chr,iG$start)

# get the BY and RM names into an appropriate form
rownames(genoG) = unifyChain(rownames(genoG),c('BY','RM'))

# check to make sure the NA run sum matches the raw number of NAs in the data
row_NA_scan(genoG,iG$chr) %>% unlist() %>% 
  (function(test)sum(test) == sum(is.na(genoG)))

# get the indices of missing genotypes
msngGeno = is.na(genoG)

# use simple imputation to deal with the spattering of NA values
genoG = impute_geno(genoG,iG$chr)

# Per segregant, see how these were imputed
mapply(function(geno_row, msng_row) table(geno_row[msng_row]), 
       as.data.frame(t(genoG)), 
       as.data.frame(t(msngGeno)), 
       SIMPLIFY = FALSE) %>%
  (function(x) do.call(rbind, lapply(x, as.data.frame))) %>%
  rownames_to_column(var = "segregant") %>%
  spread(key = Var1, value = Freq, fill = 0) %>%
  gather(key = "imputation", value = "frequency", -segregant) %>%
  spread(key = "imputation", value = frequency, fill = 0) %>%
  mutate(prefix = gsub("\\.\\d+$", "", segregant)) %>%
  group_by(prefix) %>%
  summarise(across(c(`-1`, `0`, `1`), sum, .names = "{col}")) %>%
  ungroup() %>%
  rename(segregant = prefix) %>%
  mutate(
    num_imputed = `-1` + `0` + `1`, 
    tot_markers = ncol(genoG),
    `BY(-1)` = `-1`,
    `Avg(0)` = `0`,
    `RM(1)` = `1`
  ) %>%
  select(-`-1`, -`0`, -`1`) %>%
  write.csv(paste0(asys,'genotype_QA/',
                   'Grossbach_genotype_imputation.csv'),row.names=FALSE)
rm(msngGeno)

# CREATER A MARKER KEY FOR THE GROSSBACH GENOTYPES
cG = colnames(genoG) %>%
  
  # Split and extract the chromosome part
  map_chr(~str_split(.x, ":")[[1]][1]) %>%
  
  # Create a table of chromosome counts
  
  table() %>%
  # Convert to a vector and perform the shift operation
  
  { c(1, head(cumsum(.), -1) + 1) } %>%
  # Convert to a named vector
  
  { setNames(
    as.integer(.),
    names(table(colnames(genoG) %>% map_chr(~str_split(.x, ":")[[1]][1])))) }


# create genomic range from Grossbach intervals
gG = GenomicRanges::makeGRangesFromDataFrame(iG) 

#############################################################
# LOAD ALBERT COUNTS OBJECT AND PREPARE THE PORTED GENOTYPE #
#############################################################

# GET THE ALBERT 2018 OBJECT FIRST

# load the gdata and raw counts
load(paste0(raw,'counts.RData'))

# create the genomic range from gA
gA = markRanger(colnames(counts$gdata))

# PORT THE ALBERT GENOTYPE TO GROSSBACH BY AVERAGING

# find the nearest element of gG for each element of gA
rangeBond = GenomicRanges::distanceToNearest(gA, gG)

# prepare an empty matrix with the correct dimensions
genoT = matrix(NA,nrow=nrow(counts$gdata),ncol=length(gG))

# find the unique elements of gG that elements of gA map to
uniB = unique(subjectHits(rangeBond))

# average these groups of gA mapping to the same gG
for(i in uniB){
  meanIdx = queryHits(rangeBond)[which(subjectHits(rangeBond) == i)]
  if(length(meanIdx) == 1){
    genoT[,i] = counts$gdata[,meanIdx]
  } else {
    genoT[,i] = rowMeans(counts$gdata[,meanIdx])
  }
}

# genoT information object
infoGenoT = c(rep(NA,4)) %>% 
  `names<-`(c('pass1NA_num','num_elements','pass1NA_frac','isFractional'))

# get information about genoT transformation
infoGenoT['pass1NA_num'] = sum(is.na(genoT))
infoGenoT['num_elements'] = length(genoT)
infoGenoT['pass1NA_frac'] = infoGenoT['pass1NA_num']/infoGenoT['num_elements']
infoGenoT['isFractional'] = sum(abs(genoT) != 1,na.rm=T)/length(genoT)

# figure out which elements of gA had no gG mapped to them
missingCols = which(!1:length(gG) %in% uniB)

# take these missing elements and find the nearest elements of gG to them
rangeGaps = distanceToNearest(gG[missingCols], gA)

# Correcting the filling of missing values in genoT
for(i in seq_along(missingCols)) {
  
  # get the index if the nearest
  nearestIdx = subjectHits(rangeGaps)[i]
  
  # assign the value
  genoT[, missingCols[i]] = counts$gdata[, nearestIdx]
  
}

# ensure we have the right names for downstream processing
rownames(genoT) = rownames(counts$gdata)
colnames(genoT) = colnames(genoG)

# finish and save summary csv
infoGenoT = c(infoGenoT,(summary(genoT[abs(genoT) != 1])) %>% 
                `names<-`(c('frac_min','frac_1stQ','frac_median',
                            'frac_mean','frac_3rdQ','frac_max'))) %>%
  enframe(name = "measure", value = "value")

write.csv(infoGenoT, row.names = FALSE,
          file=paste0(
            asys,'genotype_QA/',
            'Albert_2018_to_Grossbach_2022_geno_coord_move_stats.csv'))

# plot and save
p = ggplot(data.frame(genoT = genoT[abs(genoT) != 1]), aes(x = genoT)) +
  geom_histogram(bins = 30, color = "black") +
  labs(
    title = 
      "Albert 2018 ported to Grossbach 2022 fractional genotype distribution",
    x = "Fractional genotype value",
    y = "Frequency"
  ) +
  theme_bw()
ggsave(plot = p, width = 7.5, height = 5,
       filename = 
         paste0(asys,'genotype_QA/',
                'Albert_2018_to_Grossbach_2022_geno_coord_move_stats.pdf'))

# remove QA objects and move on
rm(infoGenoT,p)

########################
# PREPROCESS BREM 2002 #
########################


if (!file.exists(paste0(imed, "Brem_2002.rda"))) {
  
  print('Brem_2002 not processed and creating now...')
  
  # Load dye channels, remove rows with too much missing, and impute by averages
  GDS91 = cleanGDS(raw,'GDS91.soft',"BY-",threshNA,chrFeatures,old='-',new='.')
  GDS92 = cleanGDS(raw,'GDS92.soft',"-BY",threshNA,chrFeatures,old='-',new='.')
  
  # Merge the dye channels and attach genotype
  Brem_2002 = list(-dyeMerge(GDS91,GDS92),genoG)
  names(Brem_2002) = c('pheno','gdata')
  rm(GDS91,GDS92) 
  
  # There are no covariates so add the mean center filler and pclean is scaled
  Brem_2002$covars = matrix(1, nrow = nrow(Brem_2002$gdata)) %>%
    `rownames<-`(rownames(Brem_2002$gdata))
  Brem_2002$pclean = scale(Brem_2002$pheno)
  
  # Harmonize and save
  Brem_2002 = harmonize(Brem_2002)
  save(Brem_2002,file=paste0(imed,"Brem_2002.rda"))
  
} else {
  
  print('Brem_2002 already processed and loading now...')
  load(file=paste0(imed,"Brem_2002.rda"))
  
}

########################
# PREPROCESS BREM 2005 #
########################

if (!file.exists(paste0(imed, "Brem_2005.rda"))) {
  
  print('Brem_2005 not processed and creating now...')
  
  # Load dye channels, remove rows with too much missing, and impute by averages
  GDS1115 = cleanGDS(raw,'GDS1115.soft','Cy5',
                     threshNA,chrFeatures,old='-',new='.')
  GDS1116 = cleanGDS(raw,'GDS1116.soft','Cy3',
                     threshNA,chrFeatures,old='-',new='.')
  
  # Merge the dye channels and attach genotype
  Brem_2005 = list(-dyeMerge(GDS1115,GDS1116),genoG)
  names(Brem_2005) = c('pheno','gdata')
  rm(GDS1115,GDS1116)
  
  # There are no covariates so add the mean center filler and pclean is scaled
  Brem_2005$covars = matrix(1, nrow = nrow(Brem_2005$gdata)) %>%
    `rownames<-`(rownames(Brem_2005$gdata))
  Brem_2005$pclean = scale(Brem_2005$pheno)
  
  # Harmonize and save
  Brem_2005 = harmonize(Brem_2005)
  save(Brem_2005,file=paste0(imed,"Brem_2005.rda"))
  
} else {
  
  print('Brem_2005 already processed and loading now...')
  load(file=paste0(imed,"Brem_2005.rda"))
  
}

#########################
# PREPROCESS SMITH 2008 #
#########################

if (!file.exists(paste0(imed, "SmithEth_2008.rda")) ||
    !file.exists(paste0(imed, "SmithGlc_2008.rda"))) {
  
  print('SmithEth_2008 or SmithGlc_2008 not processed and creating now...')
  
  # Prepare phenotype data for splitting
  Smith_2008_pheno = readxl::read_excel(
    paste0(raw,'pbio.0060083.sd002.xls')) %>%
    as.data.frame()
  rownames(Smith_2008_pheno) = sapply(Smith_2008_pheno$NAME,
                                      function(x)strsplit(x,("\\|\\|"))) %>%
    (function(x)sapply(x,function(y)y[4])) %>% 
    (function(x)sapply(x,function(y)gsub(" ", "",y)))
  Smith_2008_pheno = Smith_2008_pheno[,-c(1:2)]
  nuCol = rownames(Smith_2008_pheno)
  nuRow = colnames(Smith_2008_pheno)
  Smith_2008_pheno = apply(Smith_2008_pheno,c(1,2),as.numeric) %>% as.matrix()
  Smith_2008_pheno = t(Smith_2008_pheno)
  colnames(Smith_2008_pheno) = nuCol
  rownames(Smith_2008_pheno) = nuRow
  
  # Generate metadata for the split
  SmithMeta = sapply(nuRow,function(x)strsplit(x,(" "))) %>% 
    (function(x)do.call(rbind,x)) %>%
    `colnames<-`(c('id','segregant','condition')) %>% as.data.frame()
  rm(nuCol,nuRow)
  
  # Split the ethanol and glucose phenotype data
  SmithEth_2008_pheno = Smith_2008_pheno[
    which(SmithMeta$condition  == 'ethanol'),]
  rownames(SmithEth_2008_pheno) = SmithMeta$segregant[
    which(SmithMeta$condition  == 'ethanol')]
  SmithGlc_2008_pheno = Smith_2008_pheno[
    which(SmithMeta$condition  == 'glucose'),]
  rownames(SmithGlc_2008_pheno) = SmithMeta$segregant[
    which(SmithMeta$condition  == 'glucose')]
  rm(Smith_2008_pheno,SmithMeta)
  
  # Fix names and impute missing values by average (EtOH)
  SmithEth_2008_pheno_NAs = is.na(SmithEth_2008_pheno)
  print(paste('There are', sum(SmithEth_2008_pheno_NAs),
              'NAs of', length(SmithEth_2008_pheno_NAs),
              'datapoints in SmithEth_2008_pheno'))
  SmithEth_2008_pheno = SmithEth_2008_pheno %>% 
    rmIterNA(t=threshNA) %>% averageNA()
  print('Summary of non-imputed values in SmithEth_2008_pheno')
  print(summary(SmithEth_2008_pheno[!SmithEth_2008_pheno_NAs]))
  print('Summary of imputed values in SmithEth_2008_pheno')
  print(summary(SmithEth_2008_pheno[SmithEth_2008_pheno_NAs]))
  rm(SmithEth_2008_pheno_NAs)
  rownames(SmithEth_2008_pheno) = rownames(SmithEth_2008_pheno) %>%
    (function(x) sapply(x, function(y) gsub("_", ".", y)))
  
  # Fix names and impute missing values by average (Glucose)
  SmithGlc_2008_pheno_NAs = is.na(SmithGlc_2008_pheno)
  print(paste('There are', sum(SmithGlc_2008_pheno_NAs), 
              'NAs of', length(SmithGlc_2008_pheno_NAs), 
              'datapoints in SmithGlc_2008_pheno'))
  SmithGlc_2008_pheno = SmithGlc_2008_pheno %>% 
    rmIterNA(t=threshNA) %>% averageNA()
  print('Summary of non-imputed values in SmithGlc_2008_pheno')
  print(summary(SmithGlc_2008_pheno[!SmithGlc_2008_pheno_NAs]))
  print('Summary of imputed values in SmithGlc_2008_pheno')
  print(summary(SmithGlc_2008_pheno[SmithGlc_2008_pheno_NAs]))
  rm(SmithGlc_2008_pheno_NAs)
  
  rownames(SmithGlc_2008_pheno) = rownames(SmithGlc_2008_pheno) %>%
    (function(x) sapply(x, function(y) gsub("_", ".", y)))
  
  # finish ethanol list and save
  SmithEth_2008 = list(SmithEth_2008_pheno,genoG)
  names(SmithEth_2008) = c('pheno','gdata')
  rm(SmithEth_2008_pheno)
  SmithEth_2008$pclean = scale(SmithEth_2008$pheno)
  SmithEth_2008$covars = matrix(1, nrow = nrow(SmithEth_2008$gdata)) %>%
    `rownames<-`(rownames(SmithEth_2008$gdata))
  SmithEth_2008 = harmonize(SmithEth_2008)
  save(SmithEth_2008, file = paste0(imed, "SmithEth_2008.rda"))
  
  # finish gluclose list and save
  SmithGlc_2008 = list(SmithGlc_2008_pheno,genoG)
  names(SmithGlc_2008) = c('pheno','gdata')
  rm(SmithGlc_2008_pheno)
  SmithGlc_2008$pclean = scale(SmithGlc_2008$pheno)
  SmithGlc_2008$covars = matrix(1, nrow = nrow(SmithGlc_2008$gdata)) %>%
    `rownames<-`(rownames(SmithGlc_2008$gdata))
  SmithGlc_2008 = harmonize(SmithGlc_2008)
  save(SmithGlc_2008, file = paste0(imed, "SmithGlc_2008.rda"))
  
} else {
  
  print('SmithEth_2008 and SmithGlc_2008 already processed and loading now...')
  load(file = paste0(imed, "SmithEth_2008.rda"))
  load(file = paste0(imed, "SmithGlc_2008.rda"))
  
}

########################
# PREPROCESS KHAN 2009 #
########################

if (!file.exists(paste0(imed, "Khan_2009.rda"))) {
  
  print('Khan_2009 not processed and creating now...')
  
  # Load in data and impute by average
  Khan_2009_pheno = read.delim(paste0(raw,'Khan_2009_protein.tsv'),
                               sep='\t',row.names=1,check.names = F) 
  Khan_2009_pheno = as.matrix(Khan_2009_pheno) %>% rmIterNA(t=threshNA) 
  Khan_2009_pheno_NAs = is.na(Khan_2009_pheno)
  print(paste('There are', sum(Khan_2009_pheno_NAs), 'NAs of', 
              length(Khan_2009_pheno_NAs), 'datapoints in Khan_2009_pheno'))
  Khan_2009_pheno = averageNA(Khan_2009_pheno)
  print('Summary of non-imputed values in Khan_2009_pheno')
  print(summary(Khan_2009_pheno[!Khan_2009_pheno_NAs]))
  print('Summary of imputed values in Khan_2009_pheno')
  print(summary(Khan_2009_pheno[Khan_2009_pheno_NAs]))
  rm(Khan_2009_pheno_NAs)
  
  # Fix column names
  colnames(Khan_2009_pheno) = colnames(Khan_2009_pheno) %>% 
    (function(x)sapply(x,function(y)sub("^X", "", y))) %>%
    (function(x)sapply(x,function(y)gsub("_", ".", y)))
  
  # Harmonize and save
  Khan_2009 = list(t(Khan_2009_pheno),genoG)
  names(Khan_2009) = c('pheno','gdata')
  rm(Khan_2009_pheno)
  Khan_2009$pclean = scale(Khan_2009$pheno)
  Khan_2009$covars = matrix(1, nrow = nrow(Khan_2009$gdata)) %>%
    `rownames<-`(rownames(Khan_2009$gdata))
  Khan_2009 = harmonize(Khan_2009)
  save(Khan_2009,file=paste0(imed,"Khan_2009.rda"))
  
} else {
  
  print('Khan_2009 already processed and loading now...')
  load(file=paste0(imed,"Khan_2009.rda"))
  
}

########################
# PREPROCESS FOSS 2011 #
########################

if (!file.exists(paste0(imed, "Foss_2011.rda"))) {
  
  print('Foss_2011 not processed and creating now...')
  
  # Read Foss data into dataframe
  Foss_2011_pheno = readxl::read_excel(
    paste0(raw,'foss_11_protein_measures.xls')) %>%
    as.data.frame()
  
  # Fix dimnames and make appropriate columns numeric
  rownames(Foss_2011_pheno) = Foss_2011_pheno$`systematic name`
  holder = list(colnames(Foss_2011_pheno),rownames(Foss_2011_pheno))
  Foss_2011_pheno = Foss_2011_pheno[,-c(1:3)] %>% mutate_all(as.numeric) 
  colnames(Foss_2011_pheno) = colnames(Foss_2011_pheno) %>%
    (function(x)sapply(x,function(y)gsub("_", ".", y)))
  
  # Impute by averages and log
  Foss_2011_pheno = as.matrix(Foss_2011_pheno) %>% rmIterNA(t=threshNA)
  Foss_2011_pheno_NAs = is.na(Foss_2011_pheno)
  print(paste('There are', sum(Foss_2011_pheno_NAs), 'NAs of', 
              length(Foss_2011_pheno_NAs), 'datapoints in Foss_2011_pheno'))
  Foss_2011_pheno = averageNA(Foss_2011_pheno)
  print('Summary of non-imputed values in Foss_2011_pheno')
  print(summary(Foss_2011_pheno[!Foss_2011_pheno_NAs]))
  print('Summary of imputed values in Foss_2011_pheno')
  print(summary(Foss_2011_pheno[Foss_2011_pheno_NAs]))
  rm(Foss_2011_pheno_NAs)
  Foss_2011_pheno = t(Foss_2011_pheno) %>% log2()
  
  # Harmonize and save
  Foss_2011 = list(Foss_2011_pheno,genoG)
  names(Foss_2011) = c('pheno','gdata')
  rm(Foss_2011_pheno)
  Foss_2011$pclean = scale(Foss_2011$pheno)
  Foss_2011$covars = matrix(1, nrow = nrow(Foss_2011$gdata)) %>%
    `rownames<-`(rownames(Foss_2011$gdata))
  Foss_2011 = harmonize(Foss_2011)
  save(Foss_2011,file=paste0(imed,"Foss_2011.rda"))
  
} else {
  
  print('Foss_2011 already processed and loading now...')
  load(file=paste0(imed,"Foss_2011.rda"))
  
}

###########################
# PREPROCESS PICOTTI 2013 #
###########################

if (!file.exists(paste0(imed, "Picotti_2013.rda"))) {
  
  print('Picotti_2013 not processed and creating now...')
  
  # Load SRM dataset
  Picotti_2013_pheno = read.csv(
    paste0(raw,"Picotti_2013_SRM_dataset.csv"),check.names = F) %>%
    dplyr::rename(protein = `protein /strain`)
  
  # Go from common gene names to ORF names
  Picotti_2013_pheno$protein = toupper(
    Picotti_2013_pheno$protein) %>% 
    (function(v)ifelse(v == "ARG56", "ARG5,6", v)) %>% 
    (function(x)sapply(x,function(y)sub("\\..*", "", y))) %>% 
    dictReplace(chrFeatures[c('Standard_Gene_Name','Feature_Name')])
  
  # Check if anything exceeds NA threshold in Picotti
  print("For Picotti 2013:")
  Picotti_2013_pheno %>%
    { 
      row_na_percent = apply(., 1, function(row) 
        sum(is.na(row)) / length(row) * 100)
      col_na_percent = apply(., 2, function(col) 
        sum(is.na(col)) / length(col) * 100)
      
      rows_exceeding = which(row_na_percent > threshNA * 100)
      cols_exceeding = which(col_na_percent > threshNA * 100)
      
      if (length(rows_exceeding) == 0 && length(cols_exceeding) == 0) {
        cat("No rows or columns exceed missingness threshold\n")
      } else {
        if (length(rows_exceeding) > 0) {
          cat("Rows exceeding", threshNA * 100, "% NA values:\n")
          sapply(rows_exceeding, function(i) 
            cat(sprintf("Row %d: %.2f%%\n", i, row_na_percent[i])))
        } else {
          cat("No rows exceed missingness threshold\n")
        }
        
        if (length(cols_exceeding) > 0) {
          cat("Columns exceeding", threshNA * 100, "% NA values:\n")
          sapply(cols_exceeding, function(j) 
            cat(sprintf("Column %d: %.2f%%\n", j, col_na_percent[j])))
        } else {
          cat("No columns exceed missingness threshold\n")
        }
      }
    }
  
  print('Number of Picotti phenotypes with replicate measurements:')
  Picotti_2013_pheno %>%
    group_by(protein) %>%
    filter(n() > 1) %>%
    summarise(count = n()) %>%
    nrow() %>%
    print()
  
  print('Number of Picotti phenotypes without replicate measurements:')
  Picotti_2013_pheno$protein %>%
    unique() %>%
    length() %>%
    print()
  
  print('Number of Picotti strains with strain replicates:')
  colnames(Picotti_2013_pheno[,-1]) %>%
    { .[!grepl("BY_|RM_", .)] } %>%
    sub("_.*", "", .) %>%
    { sum(duplicated(.)) } %>%
    print()
  
  print('Number of Picotti strains without strain replicates:')
  colnames(Picotti_2013_pheno[,-1]) %>%
    { .[!grepl("BY_|RM_", .)] } %>%
    sub("_.*", "", .) %>%
    unique() %>%
    { length(.) } %>%
    print()
  
  # Average separate measurements of proteins
  Picotti_2013_pheno = meanDedupe(as.matrix(Picotti_2013_pheno[,-1]),
                                  Picotti_2013_pheno$protein) %>% t()
  
  # Remove all BY and RM replicates
  Picotti_2013_pheno = Picotti_2013_pheno %>%
    { .[!grepl("BY_|RM_", rownames(.)), , drop = FALSE] }
  
  # Average separate measurements of strains
  Picotti_2013_pheno = meanDedupe(
    Picotti_2013_pheno,rownames(Picotti_2013_pheno) %>% sub("_.*", "", .))

  # Impute by average
  Picotti_2013_pheno_NAs = is.na(Picotti_2013_pheno)
  print(paste('There are', sum(Picotti_2013_pheno_NAs), 
              'NAs of', length(Picotti_2013_pheno_NAs), 
              'datapoints in Picotti_2013_pheno'))
  Picotti_2013_pheno = averageNA(Picotti_2013_pheno)
  print('Summary of non-imputed values in Picotti_2013_pheno')
  print(summary(Picotti_2013_pheno[!Picotti_2013_pheno_NAs]))
  print('Summary of imputed values in Picotti_2013_pheno')
  print(summary(Picotti_2013_pheno[Picotti_2013_pheno_NAs]))
  rm(Picotti_2013_pheno_NAs)
  
  # Harmonize and save
  Picotti_2013 = list(Picotti_2013_pheno,genoG)
  names(Picotti_2013) = c('pheno','gdata')
  rm(Picotti_2013_pheno)
  Picotti_2013$pclean = scale(Picotti_2013$pheno)
  Picotti_2013$covars = matrix(1, nrow = nrow(Picotti_2013$gdata)) %>%
    `rownames<-`(rownames(Picotti_2013$gdata))
  Picotti_2013 = harmonize(Picotti_2013)
  save(Picotti_2013,file=paste0(imed,"Picotti_2013.rda"))
  
} else {
  
  print('Picotti_2013 already processed and loading now...')
  load(file=paste0(imed,"Picotti_2013.rda"))
  
}

##########################
# PREPROCESS ALBERT 2018 #
##########################

# This is the version with the original 2018 marker coordinates
if (!file.exists(paste0(imed, "Albert_2018.rda"))) {
  
  print('Albert_2018 not processed and creating now...')
  
  # Load the log2 transformed tpm matrix 
  load(paste0(raw,'log2_t.tpm.matrix.RData')) # obj name: t.tpm.matrix
  
  # Load growth batch and OD covariates
  load(paste0(raw,"covariates.OD.RData")) # obj name: covariates.OD
  
  # Prepare mappable object and save (no harmonizaton necessary)
  Albert_2018 = list()
  Albert_2018$pheno = t.tpm.matrix
  Albert_2018$gdata = counts$gdata[,!duplicated(counts$gdata, MARGIN=2)] #dedup
  Albert_2018$pclean = scale(residuals(lm(t.tpm.matrix ~ covariates.OD)))
  Albert_2018$covars = covariates.OD
  save(Albert_2018,file=paste0(imed,"Albert_2018.rda"))
  
} else {
  
  print('Albert_2018 already processed and loading now...')
  load(file=paste0(imed,"Albert_2018.rda"))
  
}

##############################
# PREPROCESS ALBERT ALT 2018 #
##############################

# This is the version with the markers ported into Grossbach space
if (!file.exists(paste0(imed, "AlbertAlt_2018.rda"))) {
  
  # Affix the ported markers to create AlbertAlt from Albert
  print('AlbertAlt_2018 not processed and creating now...')
  AlbertAlt_2018 = Albert_2018
  AlbertAlt_2018$gdata = genoT
  save(AlbertAlt_2018,file=paste0(imed,"AlbertAlt_2018.rda"))
  
} else {
  
  print('AlbertAlt_2018 already processed and loading now...')
  load(file=paste0(imed,"AlbertAlt_2018.rda"))
  
}

#####################################
# PREPROCESS GROSSBACH 2022 GENERAL #
#####################################

if (!file.exists(paste0(imed, "GrossbachPRT_2022.rda")) ||
    !file.exists(paste0(imed, "GrossbachRNA_2022.rda"))) {
  
  # process metadata
  Grossbach_2022_meta = readxl::read_excel(
    paste0(raw,"msb202110712-sup-0005-datasetev3.xlsx"),
    sheet='metadata') %>% as.data.frame() %>%
    mutate(platePosition = as.factor(platePosition),
           cultureBatch = as.factor(cultureBatch),
           OD = as.numeric(OD),
           RNAseqRun = as.factor(RNAseqRun),
           riboAmpBatch = as.factor(riboAmpBatch),
           MSbatch_unmod,
           LCcolumn_unmod)

  Grossbach_2022_meta$strain = Grossbach_2022_meta$strain %>% 
    unifyChain(c('BY','RM'))
  
}

#################################
# PREPROCESS GROSSBACH 2022 RNA #
#################################

if (!file.exists(paste0(imed, "GrossbachRNA_2022.rda"))) {
  
  print(' GrossbachRNA_2022 not processed and creating now...')
  
  # Initialize the Grossbach 2022 RNA object
  GrossbachRNA_2022 = list()
  GrossbachRNA_2022$pheno = readxl::read_excel(
    paste0(raw,"msb202110712-sup-0003-datasetev1.xlsx"),
    sheet='transcriptome') %>% as.data.frame()
  tempGeneNames = GrossbachRNA_2022$pheno$gene
  GrossbachRNA_2022$pheno = t(GrossbachRNA_2022$pheno[,-1])
  colnames(GrossbachRNA_2022$pheno) = tempGeneNames
  
  # Initialize the residualizer object and get the Josh dedupe vector
  cleanerRNA = Grossbach_2022_meta[match(rownames(GrossbachRNA_2022$pheno),
                                        Grossbach_2022_meta$culture),] %>%
    (function(x)x[,c('OD', 'cultureBatch', 'riboAmpBatch',"strain")])

  # Print out strain replicate information
  print('Number of replicated strains')
  cleanerRNA %>% filter(!strain %in% c("BY", "RM")) %>% 
    filter(duplicated(strain) | duplicated(strain, fromLast = TRUE)) %>% 
    nrow() %>% print()
  print('Number of strains after removing replicates')
  cleanerRNA %>%   filter(!strain %in% c("BY", "RM")) %>% 
    distinct(strain) %>% nrow()
  
  # Attach the other three necessary list elements
  GrossbachRNA_2022$gdata = genoG
  GrossbachRNA_2022$pclean = matrix(NA,nrow=nrow(GrossbachRNA_2022$pheno),
                                    ncol=ncol(GrossbachRNA_2022$pheno))
  GrossbachRNA_2022$covars = cleanerRNA[
    ,c('OD', 'cultureBatch', 'riboAmpBatch')]
  
  # Regress out covariates
  print('REGRESSING COVARIATES OUT OF GROSSBACH RNA DATA')
  pb = txtProgressBar(min = 0, max = ncol(GrossbachRNA_2022$pheno), style = 3)
  for (i in 1:ncol(GrossbachRNA_2022$pheno)) {
    cleanerRNA$trait = GrossbachRNA_2022$pheno[, i]
    GrossbachRNA_2022$pclean[, i] = 
      residuals(lm(trait ~ OD + cultureBatch + riboAmpBatch,
                   data = cleanerRNA))
    setTxtProgressBar(pb, i)
  }
  
  # Deduplicate strains for the joint mapping pipeline by averaging
  GrossbachRNA_2022$pclean = meanDedupe(GrossbachRNA_2022$pclean,
                                        cleanerRNA$strain) %>% scale()
  colnames(GrossbachRNA_2022$pclean) = colnames(GrossbachRNA_2022$pheno)
  
  # Ensure one-hot encoding of covariates
  GrossbachRNA_2022$covars = oneHotCo(
    GrossbachRNA_2022$covars)
  
  # Deduplicate the individual phenos and corresponding covars by averaging
  GrossbachRNA_2022$covars = meanDedupe(
    GrossbachRNA_2022$covars,cleanerRNA$strain)
  GrossbachRNA_2022$pheno = meanDedupe(
    GrossbachRNA_2022$pheno,cleanerRNA$strain)
  
  # Remove unnecessary items, harmonize unnecessary elements, and save
  rm(cleanerRNA)
  GrossbachRNA_2022 = harmonize(GrossbachRNA_2022)
  save(GrossbachRNA_2022,file=paste0(imed,'GrossbachRNA_2022.rda'))
  
} else {
  
  print('GrossbachRNA_2022 already processed and loading now...')
  load(paste0(imed,'GrossbachRNA_2022.rda'))
  
}

#####################################
# PREPROCESS GROSSBACH 2022 PROTEIN #
#####################################

if (!file.exists(paste0(imed, "GrossbachPRT_2022.rda"))) {
  
  print('GrossbachPRT_2022 not processed and creating now...')
  
  # Begin assembling Grossbach protein
  GrossbachPRT_2022 = list()
  GrossbachPRT_2022$pheno = readxl::read_excel(
    paste0(raw, "msb202110712-sup-0003-datasetev1.xlsx"),
    sheet = 'proteome',
    na = "NA"
  ) %>% as.data.frame()
  tempGeneNames = GrossbachPRT_2022$pheno$gene
  GrossbachPRT_2022$pheno = t(GrossbachPRT_2022$pheno[,-1])
  colnames(GrossbachPRT_2022$pheno) = tempGeneNames
  
  # Initialize the residualizer object and get the Josh dedupe vector
  cleanerPRT = Grossbach_2022_meta[match(rownames(GrossbachPRT_2022$pheno),
                                         Grossbach_2022_meta$culture), ] %>% 
    (function(x)x[,c("culture","OD","cultureBatch","strain",
                     "MSbatch_unmod","proteinConc_BCA",
                     "LCcolumn_unmod", "MSmethod_unmod", "TIC_unmod")])

  # Print out strain replicate information
  print('Number of replicated strains')
  cleanerPRT %>% filter(!strain %in% c("BY", "RM")) %>% 
    filter(duplicated(strain) | duplicated(strain, fromLast = TRUE)) %>% 
    nrow() %>% print()
  print('Number of strains after removing replicates')
  cleanerPRT %>%   filter(!strain %in% c("BY", "RM")) %>% 
    distinct(strain) %>% nrow()
  
  # Attach the other three necessary list elements
  GrossbachPRT_2022$gdata = genoG
  GrossbachPRT_2022$pclean = matrix(NA,nrow=nrow(GrossbachPRT_2022$pheno),
                                    ncol=ncol(GrossbachPRT_2022$pheno))
  GrossbachPRT_2022$covars = cleanerPRT[
    ,c("OD","cultureBatch","MSbatch_unmod","proteinConc_BCA",
       "LCcolumn_unmod", "MSmethod_unmod", "TIC_unmod")]
  
  
  # Impute missing protein values
  GrossbachPRT_2022$pheno = rmIterNA(GrossbachPRT_2022$pheno, t=threshNA)
  GrossbachPRT_2022_pheno_NAs = is.na(GrossbachPRT_2022$pheno)
  print('IMPUTING MISSING PROTEIN VALUES')
  pb = txtProgressBar(
    min = 0, max = sum(apply(GrossbachPRT_2022$pheno, 2, function(x)
      any(is.na(x)))), style = 3)
  counter = 0
  for (i in which(apply(GrossbachPRT_2022$pheno, 2, function(x)
    any(is.na(x))))) {
    msng = which(is.na(GrossbachPRT_2022$pheno[, i]))
    for (j in msng) {
      result = tryCatch({
        # Impute by creating predictions using covariates
        cleanerPRT$trait = GrossbachPRT_2022$pheno[, i]
        mdl = lm(trait ~ OD + proteinConc_BCA + cultureBatch + MSbatch_unmod + 
                   TIC_unmod + LCcolumn_unmod + MSmethod_unmod, 
                 data = cleanerPRT, na.action = na.omit)
        GrossbachPRT_2022$pheno[j, i] = predict(
          mdl, cleanerPRT[j, attr(terms(mdl), "term.labels")])
      }, error = function(e) {
        # If a factor cannot be estimated, drop it and impute without it
        error_factor = sub(".*factor (\\S+) has new level.*", "\\1", e$message)
        terms_list = attr(terms(mdl), "term.labels")
        terms_list = terms_list[terms_list != error_factor]
        formula = as.formula(paste("trait ~", 
                                   paste(terms_list, collapse = " + ")))
        errMdl = lm(formula, data = cleanerPRT, na.action = na.omit)
        GrossbachPRT_2022$pheno[j, i] <<- predict(
          errMdl, cleanerPRT[j, attr(terms(errMdl), "term.labels")])
      })
    }
    counter = counter + 1
    setTxtProgressBar(pb, counter)
  }
  
  print(paste('There are', sum(GrossbachPRT_2022_pheno_NAs), 
              'NAs of', length(GrossbachPRT_2022_pheno_NAs), 
              'datapoints in GrossbachPRT_2022$pheno'))
  print('Summary of non-imputed values in GrossbachPRT_2022$pheno')
  print(summary(GrossbachPRT_2022$pheno[!GrossbachPRT_2022_pheno_NAs]))
  print('Summary of imputed values in GrossbachPRT_2022$pheno')
  print(summary(GrossbachPRT_2022$pheno[GrossbachPRT_2022_pheno_NAs]))
  rm(GrossbachPRT_2022_pheno_NAs)
  
  # Regress out covariates
  print('REGRESSING COVARIATES OUT OF GROSSBACH PRT DATA')
  pb = txtProgressBar(min = 0, max = ncol(GrossbachPRT_2022$pheno), style = 3)
  for (i in 1:ncol(GrossbachPRT_2022$pheno)) {
    cleanerPRT$trait = GrossbachPRT_2022$pheno[, i]
    GrossbachPRT_2022$pclean[, i] = residuals(
      lm(trait ~ OD + proteinConc_BCA + cultureBatch + MSbatch_unmod +
           TIC_unmod + LCcolumn_unmod + MSmethod_unmod,
         data = cleanerPRT
      ))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Deduplicate strains for the joint mapping pipeline by averaging
  GrossbachPRT_2022$pclean = meanDedupe(GrossbachPRT_2022$pclean,
                                        cleanerPRT$strain) %>% scale()
  colnames(GrossbachPRT_2022$pclean) = colnames(GrossbachPRT_2022$pheno)
  
  # Ensure one-hot encoding of covariates
  GrossbachPRT_2022$covars = oneHotCo(
    GrossbachPRT_2022$covars)
  
  # Deduplicate the individual phenos and corresponding covars by averaging
  GrossbachPRT_2022$covars = meanDedupe(
    GrossbachPRT_2022$covars,cleanerPRT$strain)
  GrossbachPRT_2022$pheno = meanDedupe(
    GrossbachPRT_2022$pheno,cleanerPRT$strain)
  
  # Remove unnecessary items, harmonize unnecessary elements, and save
  rm(cleanerPRT)
  GrossbachPRT_2022 = harmonize(GrossbachPRT_2022)
  save(GrossbachPRT_2022,file=paste0(imed,'GrossbachPRT_2022.rda'))
  
} else {
  
  print('GrossbachPRT_2022 already processed and loading now...')
  load(paste0(imed,'GrossbachPRT_2022.rda'))
  
}

####################################################
# GET HISTOGRAMS OF THE RAW AND CLEANED PHENOTYPES #
####################################################

# Get the master study object
studies = lapply(classicMetadata$study, function(x) {
  load(paste0(imed, x, '.rda'))
  get(x)
})
names(studies) = classicMetadata$study

# Plot and save histograms for raw pheno data
# Plot and save histograms for raw and clean phenotypes
phenoHister(studies = studies, dir = paste0(asys,'phenotype_histograms/'), 
            trueNames = c(
              "Brem 2002", 
              "Brem 2005", 
              "Smith 2008 ethanol", 
              "Smith 2008 glucose", 
              "Khan 2009", 
              "Foss 2011", 
              "Picotti 2013", 
              "Albert 2018 in original marker coordinates", 
              "Albert 2018 in Grossbach 2022 marker coordinates", 
              "Grossbach 2022 mRNA", 
              "Grossbach 2022 protein"
            )) 

# Final readout of segregants by genes per dataset
lapply(studies,function(x)dim(x$pheno))

#################################################################
# REMAP ALL DATASETS USING THE 2018 METHOD WITH THE SAME PARAMS #
#################################################################

# Find background QTLs (sequential wrapper; internal parallelism handled inside)
bgQTLs <- findBackgroundQTLs(studies, imed)

# Use 100-column batches for permutation chunks
CHUNK <- 100L

# Map QTLs with population structure and background QTL correction
message("Mapping FULL (pop structure + background correction)...")
t_full <- proc.time()[["elapsed"]]
qtlResults_full <- fullMap(studies, bgQTLs, imed, method = "f", chunk = CHUNK)
message(sprintf("FULL mapping finished in %.2f s", proc.time()[["elapsed"]] - t_full))

# Collapse full QTLs into a nice object
indivFullQTLs <- qtlResults_full %>%
  map_dfr(~map_dfr(.x, ~.x, .id = "chr"), .id = "study")
save(indivFullQTLs, file = paste0(imed, "indivFullQTLs.rda"))