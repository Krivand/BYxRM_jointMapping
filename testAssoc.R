#############################
#############################
## JOINT DATASET REMAPPING ##
## TESTING MAPPING PARAMS  ##
#############################
#############################

#########################################
# RETRIEVE AND FORMAT NECESSARY OBJECTS #
#########################################

# SOURCE FXNS AND CONSTANTS
source("/home/albertf/vandy147/p05/remap_functions.R")

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

################################################################
# CRUDELY ASSOCIATE ALL GENE-MARKER PAIRS FOR THE FOUR FLAVORS #
################################################################

generate_suffix = function(correctBG, correctGRM, rankTransform) {
  paste0(
    "B", ifelse(correctBG, "1", "0"),
    "G", ifelse(correctGRM, "1", "0"),
    "R", ifelse(rankTransform, "1", "0")
  )
}

# Loop over all combinations of the flags
for (correctBG in c(FALSE, TRUE)) {
  for (correctGRM in c(FALSE, TRUE)) {
    for (rankTransform in c(FALSE, TRUE)) {

      # Generate the suffix for the current combination
      suffix = generate_suffix(correctBG, correctGRM, rankTransform)

      print(paste0('DEALING WITH CASE ', suffix))

      # Run rawAssocs with the current combination of flags
      testAssocs = rawAssocs(
        studies = studies,
        targetMarkers = colnames(studies$Brem_2002$gdata),
        unionBg = unionBg,
        smetadata = scaMetadata,
        suppGenes = c("YGL009C", "YEL034W"),
        fn = paste0('testAssocs_', suffix),
        correctBG = correctBG,
        correctGRM = correctGRM,
        rankTransform = rankTransform,
        cores = 1
      )

      # Read and combine all CSV files into a single dataframe
      testAssocs = list.files(
        path = paste0(imed, '/testAssocs_', suffix),
        pattern = "*.csv", full.names = TRUE
      ) %>%
        map_df(~ read.csv(.)) %>%
        select(-ends_with(".flavor")) %>%
        rename(
          gene = geno.gene,
          marker = geno.marker
        ) %>%
        select(-ends_with(".gene"), -ends_with(".marker")) %>%
        mutate(pairID = paste(gene, marker, sep = "_"))

      # Dynamically assign the combined dataframe to a variable
      assign(paste0("testAssocs", suffix), testAssocs)

      # Save the combined dataframe as an RDA file
      save(list = paste0("testAssocs", suffix),
           file = paste0(imed, 'testAssocs', suffix, '.rda'))

      # Generate trace plots
      tracePlotter(testAssocs, c("YGL009C", "YEL034W"),
                   rootDir = paste0(asys, 'testTraces_', suffix, '/'),
                   rootName = paste0("testTrace_", suffix, "_"))

      # Delete the folder and its contents
      unlink(paste0(imed, '/testAssocs_', suffix), recursive = TRUE)

    }
  }
}