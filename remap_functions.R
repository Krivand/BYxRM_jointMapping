###########################
###########################
## ALL DATASET REMAPPING ##
## FUNCTION BANK         ##
###########################
###########################

#######################################################
# ENSURE LIBRARIES ARE CORRECTLY INSTALLED AND LOADED #
#######################################################

# Define the library path
personal_lib_path = "~/R/x86_64-pc-linux-gnu-library/4.4"

# Create the directory if it doesn't exist
if (!dir.exists(personal_lib_path)) {
  dir.create(personal_lib_path, recursive = TRUE)
}

# Add the personal library path to the library paths
.libPaths(c(personal_lib_path, .libPaths()))

# Set the default CRAN mirror
options(repos = c(CRAN = "https://lib.stat.cmu.edu/R/CRAN/"))

# Install and use BiocManager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = personal_lib_path, quietly = TRUE)
}

# All package names
pkgs = c("broom.mixed", "GEOquery", "GenomicRanges", "gridExtra", "leaps", 
         "lme4", "lmerTest", "parallel", "patchwork", "pbapply", 
         "plyr", "progress", "readxl", "rrBLUP", "tidyverse", "viridis")

# Install using correct package installation tools
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message("Installing ", pkg, "...")
    if (pkg %in% c("GEOquery", "GenomicRanges")) {
      suppressMessages(BiocManager::install(pkg, lib = personal_lib_path, quietly = TRUE))
    } else {
      suppressMessages(install.packages(pkg, lib = personal_lib_path, quietly = TRUE))
    }
  } else {
    message(pkg, " is already installed.")
  }
}

# Load packages and print messages
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop("Package ", pkg, " could not be loaded.")
  } else {
    message("Successfully loaded ", pkg)
  }
}

# Remove the list of packages
rm(pkgs)

# broom.mixed    - for extraction of p-values and estimates from mixed models
# GEOquery       - for dealing with SOFT data
# GenomicRanges  - for gene interval operations
# gridExtra      - for nice multi-panel plots
# leaps          - for which.max and others
# lme4           - for the mixed models
# lmerTest       - to test the outcomes of mixed models
# parallel	 - for the mclapply function
# patchwork      - for the nice three-panel plots
# plyr           - needed for Brem 02 processing
# progress       - for progress bars
# readxl         - for the MS Excel files
# rrBLUP         - to deal with population structure
# tidyverse      - for the good good
# viridis        - for the nicely colored plots

##################################################
# DECIDE ON CONSTANTS AND LOAD UNIVERSAL OBJECTS #
##################################################

# PREPARE DIRECTORIES AND SET WORKING DIRECTORY
wd = "/home/albertf/vandy147/p05/"
raw = paste0(wd,"raw/")
imed = paste0(wd,'imed/')
asys = paste0(wd,'asys/')
setwd(wd)

# Create directories if they do not exist
if (!dir.exists(imed)) {
  dir.create(imed)
}

if (!dir.exists(asys)) {
  dir.create(asys)
}

# choose a random seed for the script
set.seed(464788913)

# load the universally used metadata object
studyMetadata = read.delim(file=paste0(raw,'studyMetadata.csv'),
                           header=TRUE,sep=',') %>%
  as.data.frame()

# create the non-bulk subset
classicMetadata = studyMetadata[
  studyMetadata$mapping_approach == 'QTL',]

# further subset to background A in synthetic complete
scaMetadata = classicMetadata %>%
  dplyr::filter(medium == "SC", genetic_background == "A")

# study root names
studyList = c("Brem_2002", "Brem_2005", "SmithEth_2008", "SmithGlc_2008",
              "Khan_2009", "Foss_2011", "Picotti_2013", "Albert_2018", 
              "AlbertAlt_2018", "GrossbachPRT_2022", "GrossbachRNA_2022")

studyTypes = c("mRNA","mRNA","special","special",
               "protein","protein","protein",
               "mRNA","mRNA","protein","mRNA")

# threshold for NA removal
threshNA = 0.25

#########################
# GENERAL FUNCTION BANK #
#########################

# negation of in
`%nin%` = Negate(`%in%`)

# scale without a bunch of extra data
rawscale = function(m){
  scale(m) %>%
    {attr(., "scaled:center") <- NULL; .} %>%
    {attr(., "scaled:scale") <- NULL; .}%>%
    return()
  }

roundupCSV = function(inputPath, outputPath, rName, cleanUp = FALSE){
  
  # gather CSVs and save combined object
  assign(x=rName,lapply(
    list.files(path = inputPath, pattern = "\\.csv$", full.names = TRUE),
    function(x)read.delim(x,sep=',')) %>% do.call(what=rbind))
  save(list = rName, file = file.path(outputPath, paste0(rName, ".rda")))
  
  # Remove the input path if cleanUp is TRUE
  if (cleanUp) {
    unlink(inputPath, recursive = TRUE)
  }
  
}

# replace anything containing a specific pattern with that pattern
unifyPattern = function(v,p,s=F) {
  v[grepl(p, v, ignore.case = !s)] = p
  return(v)
}

# p is a vector of replacement patterns in this case
unifyChain = function(v,p,s=F){
  for(i in p){
    v = unifyPattern(v,i,s)
  }
  return(v)
}

# CURRENTLY DEPRECATED?
enumPattern = function(vec) {

  # Find which elements are not unique
  dupes = names(which(table(vec) > 1))
  
  # For each non-unique element, append a unique identifier suffix
  for (dupe in dupes) {
    # Get indices of the non-unique element
    indices = which(vec == dupe)
    # Create new names with unique identifier suffixes
    vec[indices] = paste(dupe, ".r", 1:length(indices), sep = "")
  }
  
  vec
}

# replaces elements of a vector in the keys of a dictionary with values
dictReplace = function(vec, df){
  plyr::mapvalues(vec, from = df[[1]], to = df[[2]], warn_missing = FALSE)}

# randomly selects a representitives from a repeated class
randomDedupe = function(v) {
  
  u = unique(v)
  q = integer(length(u))
  
  for(i in seq_along(u)){
    matching_indices = which(v == u[i])
    if (length(matching_indices) == 1) {
      q[i] = matching_indices
    } else {
      q[i] = sample(matching_indices, 1)
    }
  }
  
  return(q)
}

# averages replicates based on a key system
meanDedupe = function(m, v) {
  
  u = unique(v)
  
  averaged_matrix = matrix(NA, nrow = length(u), ncol = ncol(m))
  
  for (i in seq_along(u)) {
    indices = which(v == u[i])
    averaged_matrix[i, ] = colMeans(m[indices, , drop = FALSE], na.rm = TRUE)
  }
  
  rownames(averaged_matrix) = u
  colnames(averaged_matrix) = colnames(m)
  
  return(averaged_matrix)
}

# removes rows or columns too heavy on missing values
rmHeavyNA = function(d, m, t) {
  
  if (!(t > 0 && t <= 1)) {
    stop("Invalid threshold. t must be in the range (0, 1].")
  }
  
  if (is.data.frame(d)) {
    numerix = sapply(d, is.numeric)
    numerNAs = is.na(d[, numerix])
    if (m == 1) {
      rk = rowSums(numerNAs) / sum(numerix) < t
      d = d[rk, ]
    } else if (m == 2) {
      ck = colSums(numerNAs) / nrow(d) < t
      allKeep = numerix
      allKeep[numerix] = ck
      allKeep = allKeep | !numerix
      d = d[, allKeep]
    } else {
      stop("Invalid margin. Specify margin as 1 for rows or 2 for columns.")
    }
  } else if (is.matrix(d)) {
    if (!is.numeric(d)) {
      stop("Matrix d must be numeric.")
    }
    if (m == 1) {
      d = d[rowSums(is.na(d)) / ncol(d) < t, ]
    } else if (m == 2) {
      d = d[, colSums(is.na(d)) / nrow(d) < t]
    } else {
      stop("Invalid margin. Specify margin as 1 for rows or 2 for columns.")
    }
  } else {
    stop("Invalid data type. d must be a data frame or a numeric matrix.")
  }
  return(d)
}

# alternate over rows and columns until rmHeavyNA is met on both margins
rmIterNA = function(d, t) {
  
  if (!(is.matrix(d) || is.data.frame(d))) {
    stop("Invalid input. Please provide a matrix or a data frame.")
  }
  
  if (!(t > 0 && t <= 1)) {
    stop("Invalid threshold. t must be in the range (0, 1].")
  }
  
  # Make a clean pass through rmHeavyNA on both the row and column margins
  d = rmHeavyNA(d, 1, t)
  d = rmHeavyNA(d, 2, t)
  
  keepGoing = TRUE
  
  while (keepGoing) {
    
    # flag to break out of the loop
    keepGoing = FALSE
    
    # Check and remove heavy NA rows and columns
    if(is.data.frame(d)){
      numerix = sapply(d, is.numeric)
      numerNAs = is.na(d[, numerix])
      if (nrow(d) > 0 && any(rowSums(numerNAs) >= ncol(d) * t)) {
        d = rmHeavyNA(d, 1, t)
        keepGoing = TRUE
      }
      
      if (ncol(d) > 0 && any(colSums(numerNAs) >= nrow(d) * t)) {
        d = rmHeavyNA(d, 2, t)
        keepGoing = TRUE
      }
    }else if(is.matrix(d)){
      numerNAs = is.na(d)
      if (nrow(d) > 0 && any(rowSums(numerNAs) >= ncol(d) * t)) {
        d = rmHeavyNA(d, 1, t)
        keepGoing = TRUE
      }
      
      if (ncol(d) > 0 && any(colSums(numerNAs) >= nrow(d) * t)) {
        d = rmHeavyNA(d, 2, t)
        keepGoing = TRUE
      }
    }
    
    # Check if either the row or column dimension reaches 0
    if (nrow(d) == 0 || ncol(d) == 0) {
      keepGoing = FALSE
    }
  }
  
  return(d)
}

# gets strain information from Brem .soft files
getStrainInfo = function(filename) {
  
  # Read in the file line by line
  lines = base::readLines(filename)
  
  # Filter the lines that start with a GSM identifier
  gsm_lines = lines[stringr::str_detect(lines, "^#GSM")]
  
  # Extract the GSM identifier and the strain information
  gsmID = stringr::str_extract(gsm_lines, "(?<=#)(GSM\\d+)")
  strain = stringr::str_extract(gsm_lines, "(?<=: )[a-zA-Z0-9\\-\\._ ]+(?=;)")
  
  # Return a data frame with the GSM ids and strain names
  return(data.frame(gsmID = gsmID, strain = strain, stringsAsFactors = FALSE))
}

# processes the SOFT files for Brem 2002 and 2005
cleanGDS = function(dir, filename, pat, naThreshold, feats, old = F, new = F) {
  
  strainKey = getStrainInfo(paste0(dir, filename))
  
  strainKey$strain = strainKey$strain %>%
    (function(x) sapply(x, function(y) sub(pat, "", y))) %>%
    unifyChain(c('BY', 'RM'))
  
  # replace any characters if necessary 
  if (is.character(old) & is.character(new)) {
    strainKey$strain = gsub(old, new, strainKey$strain)
  }
  
  datGDS = GEOquery::getGEO(filename = paste0(dir, filename)) %>%
    GEOquery::Table()
  
  # normalize using the blanks
  datGDS[, 3:ncol(datGDS)] = apply(as.matrix(datGDS[, -c(1:2)]), 2, function(x) {
    (x - mean(x[which(datGDS$IDENTIFIER == 'blank')], na.rm = T))
  })
  
  # replace gene names with ORF names
  datGDS$IDENTIFIER = dictReplace(
    datGDS$IDENTIFIER, feats[c('Standard_Gene_Name', 'Feature_Name')]
  )
  
  # get rid of any identifiers not associated with ORFs
  datGDS = datGDS[datGDS$IDENTIFIER %in% feats$Feature_Name,]
  
  # remove rows with too many NAs
  datGDS = rmIterNA(datGDS, naThreshold)
  
  # merge duplicates by average and then impute missing data by averages
  dKey = datGDS$IDENTIFIER
  theNAs = is.na(datGDS[,-c(1:2)])
  print(paste('There are',sum(theNAs),'NAs of',length(theNAs),'datapoints'))
  datGDS = averageNA(datGDS[, -c(1:2)])
  print('Summary of non-imputed values')
  print(summary(datGDS[!theNAs]))
  print('Summary of imputed values')
  print(summary(datGDS[theNAs]))
  
  # remove gene duplicates by averaging relevant values
  datGDS = meanDedupe(datGDS, dKey)
  
  # ensure strainKey is in order andremove strain duplicates by averaging 
  common_elements = intersect(strainKey$gsmID, colnames(datGDS))
  strainKey = strainKey[strainKey$gsmID %in% common_elements, ]
  strainKey = strainKey[match(colnames(datGDS), strainKey$gsmID), ]
  datGDS = meanDedupe(t(datGDS),strainKey$strain)
  
  # return the numeric elements as a transposed matrix
  return(datGDS)
  
}

# merges two matrices
dyeMerge = function(m, n) {
  
  # Check if both m and n are matrices
  if (!is.matrix(m) || !is.matrix(n)) {
    stop("Both inputs must be matrices.")
  }
  
  # Find the intersection of row and column names
  common_row_names = intersect(rownames(m), rownames(n))
  common_col_names = intersect(colnames(m), colnames(n))
  
  # If row or column names are not the same, throw a warning
  if (!(all(rownames(m) == rownames(n)) && all(colnames(m) == colnames(n)))) {
    warning("Row and column names of the two matrices are not the same.
            Proceeding by using the intersection of row and column names.")
    
    # Remove rows and columns not in the intersection
    m = m[common_row_names, common_col_names]
    n = n[common_row_names, common_col_names]
  }
  
  # Compute the element-wise average, returning NA if one of the elements is NA
  r = mapply(function(x, y) {
    if (is.na(x) || is.na(y)) {
      return(NA)
    } else {
      return((x + y) / 2)
    }
  }, m, n)
  
  # Ensure r is a matrix with appropriate row and column names
  dim(r) = dim(m)
  rownames(r) = rownames(m)
  colnames(r) = colnames(m)
  
  return(r)
}


# Scans for NA runs and figures out what the distribution of NA runs is
row_NA_scan = function(geno, meta) {
  # get the number of rows that may contain NA runs
  num_rows = nrow(geno)
  # holds NA chain lengths for rows
  na_chain_lengths = vector("list", num_rows)
  
  for (row in 1:num_rows) {
    # Initialize the NA chain length for the row
    chain_length = 0  
    
    for (col in 1:ncol(geno)) {
      if (is.na(geno[row, col])) {
        # Increment the chain length if the value is NA
        chain_length = chain_length + 1
      } else {
        # Store the chain length if it's not empty and reset the chain length
        if (chain_length > 0) {
          na_chain_lengths[[row]] = c(na_chain_lengths[[row]], chain_length)
          chain_length = 0
        }
      }
      
      # Check for end of chromosome or row or a non-NA value
      if (col == ncol(geno) || 
          meta[row + 1] != meta[row] || 
          !is.na(geno[row, col])) {
        # Store the chain length for the current chromosome (if it's not empty)
        if (chain_length > 0) {
          na_chain_lengths[[row]] = c(na_chain_lengths[[row]], chain_length)
          chain_length = 0
        }
      }
    }
  }
  
  return(na_chain_lengths)
}

# simple end matching genome imputation
impute_geno = function(geno, meta) {

  for (i in 1:nrow(geno)) {
    # Get the indices of NA values in the current row
    na_indices = which(is.na(geno[i,]))
    
    # Check if NA values exist in the current row
    if (length(na_indices) > 0) {
      # Get the chromosome value from the meta vector for the current row
      chromosome = meta[i]
      
      # Find the boundaries of the current chromosome
      chromosome_indices = which(meta == chromosome)
      
      # Check if the chromosome has any NA values
      if (length(chromosome_indices) > 0) {
        min_index = min(chromosome_indices)
        max_index = max(chromosome_indices)
        
        # Fill NA values at boundaries with values from the same chromosome
        boundary_nas = intersect(na_indices, c(min_index, max_index))
        geno[i, boundary_nas] = geno[i, boundary_nas][1]
        
        # Get the remaining NA indices within the row
        inner_nas = setdiff(na_indices, boundary_nas)
        
        # Iterate over inner NA indices and apply the new imputation approach
        for (j in inner_nas) {
          # Search leftwards for a non-NA value with absolute value of 1
          left_value = NA
          k = j - 1
          while (k >= 1 && is.na(left_value)) {
            left_value = ifelse(abs(geno[i, k]) == 1, geno[i, k], NA)
            k = k - 1
          }
          
          # Search rightward for a non-NA value with absolute value of 1
          right_value = NA
          k = j + 1
          while (k <= ncol(geno) && is.na(right_value)) {
            right_value = ifelse(abs(geno[i, k]) == 1, geno[i, k], NA)
            k = k + 1
          }
          
          # Assign 0 if left and right values are different, 
          # else assign the matching value
          geno[i, j] = ifelse(left_value != right_value, 0, left_value)
        }
      }
    }
  }
  
  # Return the modified geno matrix
  return(geno)
}

# replaces NA values with row or column averages
averageNA = function(m, g = 1) {
  
  # Check if input is a numeric matrix or dataframe
  if (!is.numeric(m) && !is.data.frame(m)) {
    stop("Input must be a numeric matrix or a dataframe.")
  }
  
  # If a dataframe is received, convert it into a numeric matrix
  if (is.data.frame(m)) {
    m = as.matrix(m)
  }
  
  # Iterate over each column or row based on the value of g
  if (g == 1) {
    margin = 1  # Margin for row-wise averaging
    dim_names = rownames(m)
  } else if (g == 2) {
    margin = 2  # Margin for column-wise averaging
    dim_names = colnames(m)
  } else {
    stop("Invalid value for 'g'. Specify 1 for row or 2 for column.")
  }
  
  # Iterate over the specified margin of the matrix
  for (i in 1:dim(m)[margin]) {
    # Check if the margin contains any NA values
    if (g == 1) {
      if (any(is.na(m[i, ]))) {
        # Calculate the average of non-NA values in the margin
        avg = mean(m[i, !is.na(m[i, ])], na.rm = TRUE)
        
        # Replace NA values with the average in the margin
        m[i, is.na(m[i, ])] = avg
      }
    } else if (g == 2) {
      if (any(is.na(m[, i]))) {
        # Calculate the average of non-NA values in the margin
        avg = mean(m[!is.na(m[, i]), i], na.rm = TRUE)
        
        # Replace NA values with the average in the margin
        m[is.na(m[, i]), i] = avg
      }
    }
  }
  
  # Update the dimension names of the matrix
  if (g == 1) {
    rownames(m) = dim_names
  } else if (g == 2) {
    colnames(m) = dim_names
  }
  
  # Return the updated matrix
  return(m)
}

# plots histograms of the raw and clean datasets
phenoHister = function(studies, file_name, asys) {
  raw_plots = list()
  cleaned_plots = list()
  
  for (study_name in names(studies)) {
    study = studies[[study_name]]
    
    df_pheno = data.frame(mean_values = colMeans(study$pheno))
    p_pheno = ggplot(df_pheno, aes(x = mean_values)) +
      geom_histogram(fill = "blue", color = "black") +
      labs(title = paste("Raw", study_name), x = "Mean", y = "Frequency") +
      theme_minimal()
    raw_plots[[length(raw_plots) + 1]] = p_pheno
    
    df_pclean = data.frame(mean_values = colMeans(study$pclean))
    p_pclean = ggplot(df_pclean, aes(x = mean_values)) +
      geom_histogram(fill = "blue", color = "black") +
      labs(title = paste("Cleaned", study_name), x = "Mean", y = "Frequency") +
      theme_minimal()
    cleaned_plots[[length(cleaned_plots) + 1]] = p_pclean
  }
  
  # Arrange plots: all raw plots on the left, all cleaned plots on the right
  plots = c(raw_plots, cleaned_plots)
  
  # Calculate number of rows for the grid
  nCol = 2
  nRow = length(studies)
  
  # Arrange plots and save to file
  g = marrangeGrob(plots, nrow = nRow, ncol = nCol, top = NULL)
  full_path = file.path(asys, paste0(file_name, ".png"))
  ggsave(filename = full_path, plot = g, width = 16, height = 9 * nRow, limitsize = FALSE)
}

# splits one of Josh's lists into three parts
split_marker = function(marker) {
  # Split by colon and underscore
  parts = str_split(marker, "[:_]", simplify = TRUE)
  
  # Extract components
  chromosome = parts[, 1]
  index = parts[, 2]
  substitution = parts[, 3]
  
  # Return a list of components
  list(chromosome = chromosome, 
       index = as.numeric(index), 
       substitution = substitution)
}

# turns a list of markers in Josh's format to granges
# can handle if there is variant underscore separation or not
markRanger = function(v) {

  # store the original names (Josh's format)
  original_names = v
  
  # transform into gRange and add Josh names as metadata column
  sapply(v, function(x) strsplit(x, "[:_]+")) %>%
    map_dfr(~{
      data = t(.)
      if (ncol(data) == 3) {
        data.frame(data, stringsAsFactors = FALSE)
      } else {
        data.frame(cbind(data, NA), stringsAsFactors = FALSE)
      }
    }, .id = "id") %>%
    setNames(c("id", "Chromosome", "Coordinate", "Variant")) %>%
    mutate(Variant = ifelse(is.na(Variant), "", Variant)) %>%
    select(Chromosome, Coordinate, Variant) %>%
    `colnames<-`(c('chr', 'start', 'variant')) %>%
    mutate(end = start, across(c("start", "end"), as.integer)) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
    { .@elementMetadata$josh = original_names; . }
}

# Function to transform GRanges back to Josh's format
grToJosh = function(gr){
  gr_df = as.data.frame(gr)
  gr_df %>%
    rowwise() %>%
    mutate(marker = paste(seqnames, start, sep=":")) %>%
    pull(marker)
}

# wrapper for the pairwise intersection from genomic ranges
# allows for the quick addition of padding on the first grange
# e.g. bigger flanks to reduce stringency
padded_pintersect = function(gr1, gr2, pad=0) {
  
  # Add padding to gr1
  gr1_padded = resize(gr1, width(gr1) + 2 * pad, fix = "center")
  
  # Perform pairwise intersection with padded gr1
  intersections = pintersect(gr1_padded, gr2)
  
  return(intersections)
}

# (CURRENTLY NOT IN USE) changes from peaks format to genomic range
peakRanger = function(m){
  
  k = paste0('chr',as.roman(1:16))
  n = sapply(m[,1],function(x)strsplit(x,"[:_]+")) %>% 
    (function(x)do.call(rbind,x)) %>% as.data.frame()
  n[,1] = sapply(n[,1],function(x)as.numeric(which(k == x)))
  n[,3] = sapply(m[,2],function(x)strsplit(x,"[:_]+")) %>% 
    (function(x)do.call(rbind,x)) %>% (function(x)x[,2]) %>% as.numeric()
  colnames(n) = c('chr','start','end')
  n$orf = m[,3]
  n$r = m[,4]
  n = n[order(n$orf,n$chr,n$start),]
  return(makeGRangesFromDataFrame(n,keep.extra.columns = T))
  
}

# orders a matrix by its row names
rowOrd = function(m){
  return(m[order(rownames(m)), ])
}

# converts from R,B and other symbols to -1,1, and NA
rb2one = function(x){
  if(x == 'R'){1}else if(x =='B'){-1}else{NA}
}

# switches names from underscore format to Josh format
formatMarks = function(v){
  sapply(v,function(x)strsplit(x,'_')) %>%
    (function(x)
      sapply(x,function(y)paste0('chr',as.roman(y[1]),':',y[2])))
  
}

# harmonizes matrix elements of a list by segregant representation
harmonize = function(matList) {

  # Check if required items are in matList
  required_items = c("pheno", "gdata", "pclean", "covars")
  if (!all(required_items %in% names(matList))) {
    stop("The matList must contain 'pheno', 'gdata', 'pclean', and 'covars'.")
  }
  
  # Ensure the order of items in matList
  matList = matList[required_items]
  
  # Function to filter out rows containing "BY" or "RM" in the rownames
  filter_rows = function(mat) {
    rownames_to_keep = !grepl("BY|RM", rownames(mat))
    return(mat[rownames_to_keep, , drop = FALSE])
  }
  
  # Apply the filter to each matrix
  matList$pheno = filter_rows(matList$pheno)
  matList$gdata = filter_rows(matList$gdata)
  matList$pclean = filter_rows(matList$pclean)
  matList$covars = filter_rows(matList$covars)
  
  # Get the rownames of the elements pheno, gdata, pclean, and covars
  rownames_pheno = rownames(matList$pheno)
  rownames_gdata = rownames(matList$gdata)
  rownames_pclean = rownames(matList$pclean)
  rownames_covars = rownames(matList$covars)
  
  # Intersect these rownames
  common_rownames = Reduce(intersect, list(
    rownames_pheno, rownames_gdata, rownames_pclean, rownames_covars))
  
  # Subset each matrix so that they only contain the intersection
  matList$pheno = matList$pheno[common_rownames, , drop = FALSE]
  matList$gdata = matList$gdata[common_rownames, , drop = FALSE]
  matList$pclean = matList$pclean[common_rownames, , drop = FALSE]
  matList$covars = matList$covars[common_rownames, , drop = FALSE]
  
  # Reorder gdata, pclean, and covars to match the order of pheno
  matList$gdata = matList$gdata[rownames(matList$pheno), , drop = FALSE]
  matList$pclean = matList$pclean[rownames(matList$pheno), , drop = FALSE]
  matList$covars = matList$covars[rownames(matList$pheno), , drop = FALSE]
  
  return(matList)
}

oneHotCo = function(covars) {
  # Identify columns that are either factors or characters
  categorical_columns = sapply(covars, function(x) is.factor(x) || is.character(x))
  covars_encoded = covars
  
  for (col_name in names(covars)[categorical_columns]) {
    covar_factor = as.factor(covars[[col_name]])
    dummy_vars = model.matrix(~ covar_factor - 1)
    colnames(dummy_vars) = paste(col_name, levels(covar_factor), sep = "_")
    covars_encoded = cbind(covars_encoded, dummy_vars)
    covars_encoded[[col_name]] = NULL
  }
  
  # Convert to matrix
  covars_encoded = as.matrix(covars_encoded)
  
  return(covars_encoded)
}

# QTL mapping dependency (J. Bloom)
calc.BLUPS= function(G,Z,Vinv,y,X,B){G%*%crossprod(Z,Vinv)%*%(y-X%*%B)}

# QTL mapping dependency (J. Bloom)
doEigenA_forMM=function(pheno.scaled,A ,X=NULL ) {
  n=nrow(pheno.scaled)
  if(is.null(X) ) {  X = matrix(rep(1, n), n, 1); p=1} else {p=ncol(X)}
  XtX = crossprod(X, X)
  XtXinv = solve(XtX)
  S = diag(n) - tcrossprod(X %*% XtXinv, X)
  SHbS = S %*% A %*% S
  SHbS.system = eigen(SHbS, symmetric = TRUE)
  theta = SHbS.system$values[1:(n - p)] 
  Q = SHbS.system$vectors[, 1:(n - p)]
  return(list(theta=theta, Q=Q))
}

# QTL mapping dependency (J. Bloom)
# borrowed and modified from from mixed.solve() in rrBLUP 08/02/15
# super optimized for one VC and fixed effects 
# ~1000X speedup by precomputing eigen decomp
m.S=function (y, K = NULL, bounds = c(1e-09, 1e+09),theta=NULL,Q=NULL,X=NULL) 
{
  n = length(y)
  y = matrix(y, n, 1)
  if(is.null(X) ) {  p = 1    } else { p = ncol(X) }
  Z = diag(n)
  m = ncol(Z)
  
  omega = crossprod(Q, y)
  omega.sq = omega^2
  
  f.REML = function(lambda, n.p, theta, omega.sq) {
    n.p * log(sum(omega.sq/(theta + lambda))) + sum(log(theta + lambda))
  }
  soln = optimize(f.REML, interval = bounds, n - p, theta,  omega.sq)
  lambda.opt = soln$minimum
  
  df = n - p
  Vu.opt = sum(omega.sq/(theta + lambda.opt))/df
  Ve.opt = lambda.opt * Vu.opt
  VCs=c(Vu.opt, Ve.opt)
  return(VCs)
}

# Faster LOD calculation (J. Bloom)
# no NA handling but very fast
# input is scaled phenotypes and genotypes
fasterLOD = function(n.pheno, pheno.s,gdata.s, betas=FALSE, sdx=1, pheno=NULL){
  r=crossprod(pheno.s, gdata.s)/(n.pheno-1)
  LOD=(-n.pheno*log(1-r^2))/(2*log(10))
  if(betas==FALSE) {
    return(LOD)
  } else {
    return(list(r=r, LOD=LOD))
  }
}

# Finds background QTLs (J. Bloom)
# covariates is numeric and segregant by number of covariates
# phenotype.raw is numeric with dim: segregants x number of traits 
# phenotype.scaled will be the same type as phenotype.raw
# gdata is numeric with dim: segregants x number of markers
  # the gdata colnames must be in a special format
  # chr:position e.g. chrXIII:591795 for the bp position on ChrXIII 
# gdata.scaled is just a scaled version of gdata
# bcl.thresh is a numeric LOD cutoff
# chromosomes is a vector of strings 
  # these should be present in the gdata names
  # lets you specify those other than the typical 16 in S. cerevisiae
# OUTPUT: a list of length equal to the number of chromosomes
  # each element (e.g. ChrI) contains a list of genes
  # genes have markers from OTHER chromosomes that influence them
  # this is used to remove genetic influence coming from other chromosomes
  # if you were scanning ChrI you would remove QTLs from ChrII through ChrXVI

find.background.QTL = function(covariates, phenotype.raw, pheno.scaled, gdata,
                               gdata.scaled, bck.thresh=3.5, chromosomes=NULL){
  
  # if no chromosome is provided, default to the 16 yeast chromosomes
  if(is.null(chromosomes)) { chromosomes=paste0('chr', as.roman(1:16)) } 
  
  # this vector tells you which chromosome each marker belongs to
  cvec=(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])
  
  # ITERATION 1  (find background QTL effects)
  LODS=fasterLOD(nrow(pheno.scaled),pheno.scaled,gdata.scaled)
  
  # lists to hold LODs, markers, and scaled markers by chromosome
  LODS.by.chr=list()
  gdata.s.by.chr=list()
  gdata.by.chr=list()
  for(cc in chromosomes) {   
    LODS.by.chr[[cc]]=LODS[,which(cvec %in% cc)]   
    gdata.s.by.chr[[cc]]=gdata.scaled[,which(cvec %in% cc)]   
    gdata.by.chr[[cc]]=gdata[,which(cvec %in% cc)]   
  }
  
  # max LOD per chromosome
  max.lod.chr=lapply(LODS.by.chr, function(x) { apply(x,1,max) })
  
  # index of max LOD per chromosome
  max.lod.chr.ind=lapply(LODS.by.chr, function(x) { apply(x,1,which.max) })
  
  # pulls out the marker index (by chromosome) of the peak marker 
  pmarker.chr.ind =  mapply(function(x,y,z) {
    
    # index gets included if LOD exceeds background threshold
    h=y[x>bck.thresh]
    
    # Retrieve marker names using indices from h; i holds the actual names of 
      # the markers at indices where LOD scores are significant.
    i=colnames(z)[h]
    
    # Assign names to i using the original indices in h. This step links each
      # marker name with its corresponding LOD score index to trace it
      # makes clear which LOD scores each marker name corresponds to.
    names(i)=names(h)
    
    return(i)
    
  },  max.lod.chr, max.lod.chr.ind,gdata.by.chr)
  
  # Convert markers by chromosome to a data frame
  pmarker.chr.ind.df=data.frame(
    chr=rep(names(pmarker.chr.ind),sapply(pmarker.chr.ind, length)),
    gene=as.vector(do.call('c', lapply(pmarker.chr.ind, function(x) names(x)))),
    marker=as.vector(do.call('c', lapply(pmarker.chr.ind, function(x) (x)))), 
    stringsAsFactors=F)
  pbg=split(pmarker.chr.ind.df$marker, pmarker.chr.ind.df$gene)
  
  # Adjust phenotypes for identified background QTL effects, iteration 1
  plist=lapply(colnames(phenotype.raw), function(i) {
    
    # if the gene not present in the QTL list return unchanged 
    if(!(i %in% names(pbg)) ) {return(pheno.scaled[,i]) } 
    
    # if gene present remove effect of gdata
    else{ residuals(lm(phenotype.raw[,i]~covariates+gdata[,pbg[[i]]]) )  }
  })
  
  # Combine residuals of all traits into a single matrix 
  # this matrix has the same dimensions as pheno.scaled
  p.resid1=do.call('cbind', plist)
  
  # Assign original phenotype names to the columns of the residuals matrix
  colnames(p.resid1)=colnames(pheno.scaled)
  
  # Normalize residuals to have zero mean and unit variance for further analysis
  pheno.scaled=scale(p.resid1)
  
  # ITERATION 2 (find background QTL effects)
  LODS1=fasterLOD(nrow(pheno.scaled),pheno.scaled,gdata.scaled)
  LODS1.by.chr=list()
  for(cc in chromosomes) { LODS1.by.chr[[cc]]=LODS1[,which(cvec %in% cc)]}
  max.lod1.chr = lapply(LODS1.by.chr, function(x) { apply(x,1,max) })
  max.lod1.chr.ind = lapply(LODS1.by.chr, function(x) { apply(x,1,which.max) })
  pmarker1.chr.ind = mapply(function(x,y,z) { 
    h=y[x>bck.thresh]
    i=colnames(z)[h]
    names(i)=names(h)
    return(i)
  },  max.lod1.chr, max.lod1.chr.ind,gdata.by.chr)
  pmarker1.chr.ind.df=data.frame(
    chr=rep(names(pmarker1.chr.ind), sapply(pmarker1.chr.ind, length)),
    gene=as.vector(do.call('c', lapply(pmarker1.chr.ind,function(x)names(x)))),
    marker=as.vector(do.call('c', lapply(pmarker1.chr.ind,function(x)(x)))),
    stringsAsFactors=F)
  pbg1=split(pmarker1.chr.ind.df$marker, pmarker1.chr.ind.df$gene)
  pbg.merge=pbg
  for(i in names(pbg1)){pbg.merge[[i]]=c(pbg.merge[[i]], pbg1[[i]])}
  
  # ITERATION 3 (find background QTL effects)
  plist=lapply(colnames(phenotype.raw),function(i) { 
    if(!(i %in% names(pbg.merge)) ) {  return(pheno.scaled[,i])} 
    else{ residuals(lm(phenotype.raw[,i]~covariates+gdata[,pbg.merge[[i]]]))}
  })
  p.resid2=do.call('cbind', plist)
  colnames(p.resid2)=colnames(pheno.scaled)
  pheno.scaled=scale(p.resid2)
  LODS2=fasterLOD(nrow(pheno.scaled),pheno.scaled,gdata.scaled)
  LODS2.by.chr=list()
  for(cc in chromosomes) {    LODS2.by.chr[[cc]]=LODS2[,which(cvec %in% cc)]}
  max.lod2.chr=lapply(LODS2.by.chr, function(x) { apply(x,1,max) })
  max.lod2.chr.ind=lapply(LODS2.by.chr, function(x) { apply(x,1,which.max) })
  pmarker2.chr.ind =  mapply(function(x,y,z) { 
    h=y[x>bck.thresh]
    i=colnames(z)[h]
    names(i)=names(h)
    return(i)
  },  max.lod2.chr, max.lod2.chr.ind,gdata.by.chr)
  pmarker2.chr.ind.df=data.frame(
    chr=rep(names(pmarker2.chr.ind), sapply(pmarker2.chr.ind, length)),
    gene=as.vector(do.call('c', lapply(pmarker2.chr.ind,function(x)names(x)))),
    marker=as.vector(do.call('c', lapply(pmarker2.chr.ind,function(x)(x)))), 
    stringsAsFactors=F)
  pbg2=split(pmarker2.chr.ind.df$marker, pmarker2.chr.ind.df$gene)
  
  pbg.final=pbg.merge
  for(i in names(pbg2)){pbg.final[[i]]=c(pbg.final[[i]], pbg2[[i]])}
  
  pm.chr.ind.df=rbind(pmarker.chr.ind.df, pmarker1.chr.ind.df, 
                      pmarker2.chr.ind.df)
  background.QTL=list()
  for(cc in chromosomes){
    background.QTL[[cc]]=split(pm.chr.ind.df$marker[pm.chr.ind.df$chr!=cc], 
                               pm.chr.ind.df$gene[pm.chr.ind.df$chr!=cc]) 
  }
  return(background.QTL)
}

# maps QTLs (Josh Bloom)
# covariates is a matrix of the covaraites
# background QTLs are 
# phenotype.raw is an unsacaled phenotypes matrix (segregant x trait)
  # elements are numeric measures
# pheno.scaled is the scaled pheno matrix
# gdata is an unscaled genotypes matrix (segregant x marker)
  # elements are typically binary -1 or 1 coded genotypes but can be numeric
# gdata.scaled is the scaled genotypes matrix
mapQTL=function(covariates, background.QTL, phenotype.raw, pheno.scaled, 
                gdata, gdata.scaled, 
                chromosomes=NULL,
                n.perm=1000, FDR.thresh=.05, Arandom=TRUE) {
  
  # if no chromosome is provided, default to the 16 yeast chromosomes
  if(is.null(chromosomes)) { chromosomes=paste0('chr', as.roman(1:16)) } 
  
  # this vector tells you which chromosome each marker belongs to
  cvec=(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])
  
  # get gdata and scaled gdata by chromosome
  gdata.s.by.chr=list()
  gdata.by.chr=list()
  for(cc in chromosomes) {   
    gdata.s.by.chr[[cc]]=gdata.scaled[,which(cvec %in% cc)]   
    gdata.by.chr[[cc]]=gdata[,which(cvec %in% cc)]   
  }
  
  # initialize peak holding list
  peakList=list()
  
  # iterate over the chromosomes in the phenotype matrix
  for(cc in chromosomes) {  
    
    # tell us what chromosome we are on
    print(cc) 
    
    # iterate over the genes in the phenotype matrix
    plist=lapply(colnames(phenotype.raw),function(i){ 
      
      # if there is no QTL for that gene
      if(!(i %in% names(background.QTL[[cc]]))){
        
        # then return only the covariate correction
        return(residuals(lm(phenotype.raw[,i]~covariates)))}
      
      else{
        
        # else return correction after 
        residuals(lm(
          phenotype.raw[,i]~covariates+gdata[,background.QTL[[cc]][[i]]]))}
    })
    
    # preal will be the same size as the phenotype matrix
    # it combines the residual slices
    preal=do.call('cbind', plist)
    
    # colnames for preal are lifted directly from the phenotype matrix
    colnames(preal)=colnames(phenotype.raw)
    
    # the presid close is made of preal
    presid = preal
    
    # remove residual additive effect
    if(Arandom) {
      Aloco=rrBLUP::A.mat(
        do.call('cbind',(gdata.by.chr[-match(cc, names(gdata.by.chr))])))/2
      eigA=doEigenA_forMM(presid,Aloco)
      svdAloco=svd(Aloco)
      pb=txtProgressBar(min=1, max=ncol(presid), style=3)
      for(tp in 1:ncol(presid)){
        setTxtProgressBar(pb,tp)
        rr=m.S(presid[,tp], K=Aloco,  theta=eigA$theta, Q=eigA$Q)
        W=svdAloco$u %*% tcrossprod(diag(1/((svdAloco$d*rr[1])+(rr[2]))),
                                    svdAloco$v)
        if(rr[1]>0) {
          blups=calc.BLUPS(rr[1]*Aloco,diag(nrow(presid)),
                           W,presid[,tp],matrix(1,nrow(presid),1),0 )[,1]
          presid[,tp]=as.vector(presid[,tp] - blups)
        }
      }
      rm(W)
      close(pb)
    }
    
    preal=presid
    presid=scale(presid)
    
    obsLOD = fasterLOD(nrow(preal),
                       presid,gdata.s.by.chr[[cc]],betas=TRUE,pheno=preal)
    
    # this flag forms the basis of a while loop
    findingPeaks=TRUE
    jump=1
    
    # while we remain in peakfinding mode we scan
    while(findingPeaks==TRUE) {
      
      # we find the maximum LODs
      max.obsLOD=apply(obsLOD$LOD,1,max)
      
      # we retrieve the indices (ind) for the maximum LODs
      max.ind.obsLOD=apply(obsLOD$LOD,1,which.max)   
      
      # we prepare the LOD permutation matrix
      permLOD=matrix(0, length(max.obsLOD), n.perm)
      
      # announce permutations and start progress bar
      print('doing permutations')
      pb = txtProgressBar(min = 1, max =n.perm, style = 3)
      
      # loop for the number of permutations
      for(i in 1:n.perm) {
        setTxtProgressBar(pb, i)
        permLOD[,i]=apply(fasterLOD(nrow(presid),presid[sample(nrow(presid)),],
                                    gdata.s.by.chr[[cc]]),1,max)
      }
      close(pb)
      
      # check how many are above the iterated LOD threshold
      obsPcnt = sapply(seq(1.5,9,.05),
                       function(thresh){sum(max.obsLOD>thresh)})
      names(obsPcnt) = seq(1.5, 9, .05)
      
      # if there are no more than 2 QTLs discovered, break the loop
      # we need to discover at least three to say that it's real
      if(sum(obsPcnt)<3) {break}
      
      # take the means of the permuted LODs to get the expected number of FDs
      expPcnt = sapply(seq(1.5, 9, .05),  
                       function(thresh) { 
                         mean(apply(permLOD,2,function(ll){sum(ll>thresh)}))
                       } )
      names(expPcnt) = seq(1.5, 9, .05)
      
      # FDR is then just the expected percents divided by the observed
      pFDR = expPcnt/obsPcnt
      
      # any numerical weirdness due to div by 0 or missingness is zeroed here
      pFDR[is.na(pFDR)]=0
      pFDR[!is.finite(pFDR)]=0
      
      # ensure monotonicity in FDR curve by cumulative max on double reverse
      pFDR = rev(cummax(rev(pFDR)))
      
      # linear interpolation of function
      fdrFX=approxfun(pFDR, seq(1.5,9,.05))
      thresh=fdrFX(FDR.thresh)
      print(paste('FDR Thresh 5%', thresh)) 
      peak.gene.ind=max.ind.obsLOD[which(max.obsLOD>thresh)]
      
      # continue 
      if(length(peak.gene.ind)>0) { 
        findingPeaks=TRUE 
        pldf=data.frame(gene=names(peak.gene.ind),
                        pmarker=colnames(gdata.by.chr[[cc]])[peak.gene.ind], 
                        pcind=as.vector(peak.gene.ind),
                        stringsAsFactors=F)
        pldf$r=obsLOD$r[cbind(pldf$gene, pldf$pmarker)]
        pldf$LOD=obsLOD$LOD[cbind(pldf$gene, pldf$pmarker)]
        lodInt=matrix(NA,nrow(pldf),2)
        
        for(peak in 1:nrow(pldf)){
          int.range=range(
            which(obsLOD$LOD[pldf[peak,]$gene,]>pldf[peak,]$LOD-1.5)) 
          lodInt[peak,]=colnames(gdata.by.chr[[cc]])[int.range]
        }
        
        pldf$CI.l=lodInt[,1]
        pldf$CI.r=lodInt[,2]
        peakList[[cc]][[jump]]=pldf
        
        presid.tmp=matrix(0, nrow(phenotype.raw), nrow(pldf))
        rownames(presid.tmp)=rownames(phenotype.raw)
        colnames(presid.tmp)=pldf$gene
        
        pldfm = do.call('rbind', peakList[[cc]])
        found.q.covs=split(pldfm$pmarker, pldfm$gene)
        
        for(i in pldf$gene){
          presid.tmp[,i]=scale(
            residuals(lm(presid[,i]~gdata[,found.q.covs[[i]]])))
        }
        obsLOD = fasterLOD(nrow(presid.tmp),presid.tmp,
                           gdata.s.by.chr[[cc]],betas=TRUE, pheno=presid.tmp)
        presid=presid.tmp
        jump=jump+1 
      } else {findingPeaks=FALSE} # set peakfinding to false if we find none
    }
  }
  return(peakList)
}

# LOD from n and R
# inputs: sample size (n) and correlation coefficient (r)
# outputs: LOD score
nrLOD = function(n, r){
  return(-n*log(1-r^2))/(2*log(10))
}

# wrapper for getting background QTLs
findBackgroundQTLs = function(studies, fn, cores = 1) {
  bgQTLs = mclapply(names(studies), function(study_name) {
    study = studies[[study_name]]
    
    if (!file.exists(paste0(fn, study_name, '_Background.rda'))) {
      print(paste0("MAPPING BACKGROUND QTLs for ", study_name))
      bgQTL = find.background.QTL(study$covars, 
                                  study$pheno, 
                                  scale(study$pheno),
                                  study$gdata,
                                  scale(study$gdata))
      # Save the object with the correct name
      assign(paste0(study_name, '_Background'), bgQTL)
      save(list = paste0(study_name, '_Background'), 
           file = paste0(fn, study_name, '_Background.rda'))
    } else {
      print(paste0("LOADING BACKGROUND QTLs for ", study_name))
      load(paste0(fn, study_name, '_Background.rda'))
      bgQTL = get(paste0(study_name, '_Background'))
    }
    
    return(bgQTL)
  }, mc.cores = cores)
  
  names(bgQTLs) = names(studies)
  return(bgQTLs)
}

# wrapper for mapping QTLs using Josh's pipeline over a list
fullMap = function(studies, bgQTLs, fn, method = 'f', cores = 1) {
  
  # EXPLANATION OF MAPPING TYPES
  # full map (f): includes empirical null, pop struc crxn, background crxn
  # simple map (s): includes empirical null
  
  # Alter mapping type based on method
  if (method == 'f') {
    typeRoot = '_fullPeaks'
    popCrxn = TRUE
  } else if (method == 's') {
    typeRoot = '_simplePeaks'
    popCrxn = FALSE
  } else {
    warning('Invalid method selected. Choose f or s. Returning NULL')
    return(NULL)
  }
  
  # Perform QTL mapping for each study
  qtlResults = mclapply(names(studies), function(study_name) {
    study = studies[[study_name]]
    background = if (method == 'f') bgQTLs[[study_name]] else list()
    
    # Check if results file exists
    if (file.exists(paste0(fn, study_name, typeRoot, '.rda'))) {
      print(paste0('Files already exist for ', study_name, ' using method ', method))
      print('Loading files now...')
      load(paste0(fn, study_name, typeRoot, '.rda'))
      qtlResult = get(paste0(study_name, typeRoot))
    } else {
      print(paste0('Attempting to map QTLs for ', study_name, ' using method ', method))
      
      # Scan for QTLs
      print("SCANNING FOR QTLS")
      qtlResult = mapQTL(study$covars, 
                         background,
                         study$pheno, scale(study$pheno), 
                         study$gdata, scale(study$gdata),
                         n.perm = 1000, FDR.thresh = .05,
                         Arandom = popCrxn)
      
      # Save results to file using list trick
      assign(paste0(study_name, typeRoot), qtlResult)
      save(list = paste0(study_name, typeRoot), 
           file = paste0(fn, study_name, typeRoot, '.rda'))
    }
    
    return(qtlResult)
  }, mc.cores = cores)
  
  # Assign names to results and return
  names(qtlResults) = names(studies)
  return(qtlResults)
}

# transforms peaks into a format usable by Matthew's code
easyTableQTL = function(peaks, ivals) {
  # First, we separate the chromosome and coordinate from pmarker, CI.l and CI.r
  peaks = peaks %>%
    separate(CI.l, into = c("chr", "CI.l"), sep = ":", convert = TRUE) %>%
    separate(CI.r, into = c("CI.r_chr", "CI.r"), sep = ":", remove = FALSE) %>%
    select(-CI.r_chr,-pcind) %>% mutate(CI.r = as.numeric(CI.r))
  
  # Now, we lookup the 'end' for CI.r_coord in the ivals table
  peaks = peaks %>%
    left_join(ivals, by = c("chr" = "chr", "CI.r" = "start")) %>%
    mutate(pmarker = round((CI.l + CI.r) / 2))
  
  # Format the output as desired
  peaks = peaks %>%
    select(gene, chr, pmarker, r, CI.l, CI.r,LOD)
  
  return(peaks)
}

# takes a nice table and creates a genomic range out of it
qtlTableToRange = function(df) {
  
  gr_df = df %>% 
    dplyr::rename(seqnames = chr, start = CI.l, end = CI.r) %>%
    dplyr::mutate(start = as.integer(start),
                  end = as.integer(end),
                  gene = as.character(gene),
                  pmarker = as.integer(pmarker),
                  r = as.numeric(r)) %>%
    dplyr::select(seqnames, start, end, gene, pmarker, r, LOD)
  
  gr = makeGRangesFromDataFrame(gr_df, keep.extra.columns = TRUE)
  
  return(gr)
}

# creates a simple gRange with gene name and marker metadata
geneRanger = function(df) {
  
  gr_df = df %>% 
    dplyr::rename(seqnames = chr, start = CI.l, end = CI.r) %>%
    dplyr::mutate(start = as.integer(start),
                  end = as.integer(end),
                  gene = as.character(gene),
                  pmarker = as.integer(pmarker)) %>%
    dplyr::select(seqnames, start, end, gene, pmarker)
  
  gr = makeGRangesFromDataFrame(gr_df, keep.extra.columns = TRUE)
  
  return(gr)
}

# compares a metadata measure between 'same' QTLs of two QTL granges
compr = function(old,new){
  
  # figure out which items are the same
  q = GenomicRanges::findOverlaps(old,new)
  k = old$gene[S4Vectors::queryHits(q)] == new$gene[S4Vectors::subjectHits(q)]
  return(q[k])
  
}

# extracts columns based on a overlap key
olapExtract = function(key,old,new,oCols,oNames,nCols,nNames){
  
  old = old[queryHits(key),oCols]
  new = new[subjectHits(key),nCols]
  colnames(old) = oNames
  colnames(new) = nNames
  return(cbind(old, new))
  
}

# wrapped for quick scatterplot creation
scatterer = function(df, colA, colB, title) {
  ggplot(df, aes(x = !!sym(colA), y = !!sym(colB))) +
    geom_point(alpha = 0.2, pch = 19) +
    theme_minimal() +
    labs(title = title, x = colA, y = colB)
}

############################
# JOINT ANALYSIS FUNCTIONS #
############################

# Collapse background QTLs across studies for joint mapping
collapseBgQTLs = function(bgQTLs, chroms = NULL) {
  
  # Define default 16 S. cerevisiae chromosomes if none are provided
  if (is.null(chroms)) {
    chroms =  paste0("chr", as.roman(1:16))
  }
  
  # Initialize an empty list to store the combined results
  combinedQTLs = list()
  
  # Loop through each study
  for (study in bgQTLs) {
    
    # Loop through each chromosome in the study
    for (chr in chroms) {
      
      if (!is.null(study[[chr]])) {
        
        # Loop through each gene in the chromosome
        for (gene in names(study[[chr]])) {
          
          # Gather all markers for this gene
          markers = study[[chr]][[gene]]
          
          # Store the markers in the combined list by gene first
          if (!gene %in% names(combinedQTLs)) {
            combinedQTLs[[gene]] = list()
          }
          
          if (!chr %in% names(combinedQTLs[[gene]])) {
            combinedQTLs[[gene]][[chr]] = markers
          } else {
            combinedQTLs[[gene]][[chr]] = unique(c(combinedQTLs[[gene]][[chr]], markers))
          }
        }
      }
    }
  }
  
  return(combinedQTLs)
}

# check which genes are represented in two mRNA and two protein datasets
suppTwo = function(studyMetadata, type = "analyte") {
  
  if (type == "analyte") {
    typeA = "mRNA"
    typeB = "protein"
  } else if (type == "medium") {
    typeA = "SC"
    typeB = "YNB2G"
  } else {
    stop("Invalid type specified. Use 'analyte' or 'medium'.")
  }
  
  # Filter data for typeA and typeB
  typeA_data = studyMetadata[studyMetadata[[type]] == typeA, ]
  typeB_data = studyMetadata[studyMetadata[[type]] == typeB, ]
  
  # Extract and count genes
  typeA_genes = table(unlist(typeA_data$genes))
  typeB_genes = table(unlist(typeB_data$genes))
  
  # Identify genes appearing in at least two datasets
  typeA_valid_genes = names(typeA_genes[typeA_genes >= 2])
  typeB_valid_genes = names(typeB_genes[typeB_genes >= 2])
  
  # Find common genes that meet the criteria for both typeA and typeB
  common_genes = intersect(typeA_valid_genes, typeB_valid_genes)
  
  return(common_genes)
}

# STUDY-STRATIFIED PERMUTATIONS
# takes as argument the group boundaries and number of perms to request
# critical for the functionality of the jointMarkerTest
stratPerm = function(v, n) {
  
  unique_values = unique(v)
  partitions = lapply(unique_values, function(x) which(v == x))
  
  perm_matrix = matrix(nrow = n, ncol = length(v))
  
  for (i in 1:n) {
    permuted_indices = 1:length(v)
    for (part in partitions) {
      permuted_indices[part] = sample(part, length(part), replace = FALSE)
    }
    perm_matrix[i, ] = permuted_indices
  }
  return(perm_matrix)
}

# CORRECT FOR PREVIOUSLY DISCOVERED BACKGROUND QTLS
# Function to adjust jointPhenoScaled for background QTLs
adjust_for_background_QTLs = function(
    geno_matrix, pheno_matrix, bgQTLs, chromosomes) {
  adjusted_pheno = pheno_matrix
  for (cc in chromosomes) {
    bg_markers = bgQTLs[[cc]]
    if (length(bg_markers) > 0) {
      for (gene in colnames(pheno_matrix)) {
        if (gene %in% names(bg_markers)) {
          adjusted_pheno[, gene] = residuals(
            lm(pheno_matrix[, gene] ~ geno_matrix[, bg_markers[[gene]]]))
        }
      }
    }
  }
  return(adjusted_pheno)
}

# CORRECT FOR GENETIC RELATEDNESS ACROSS SEGREGANTS
# Function to calculate BLUPs and correct relatedness out of phenotype
# NOTE: breaks when same segregant is in there multiple times (100% related)
calculate_BLUPs = function(geno_matrix, pheno_matrix) {
  Aloco = rrBLUP::A.mat(geno_matrix)/2
  eigA = doEigenA_forMM(pheno_matrix, Aloco)
  svdAloco = svd(Aloco)
  adjusted_residuals = pheno_matrix
  for (tp in 1:ncol(pheno_matrix)) {
    rr = m.S(pheno_matrix[,tp], K = Aloco, theta = eigA$theta, Q = eigA$Q)
    W = svdAloco$u %*% tcrossprod(diag(1/((svdAloco$d * rr[1]) + (rr[2]))),
                                  svdAloco$v)
    if (rr[1] > 0) {
      blups = calc.BLUPS(rr[1] * Aloco, diag(
        nrow(pheno_matrix)), W, pheno_matrix[,tp],
        matrix(1, nrow(pheno_matrix), 1), 0)[,1]
      adjusted_residuals[,tp] = as.vector(pheno_matrix[,tp] - blups)
    }
  }
  return(adjusted_residuals)
}

# MAP ASSOCIATION BETWEEN SPECIFIC MARKERS AND EXPRESSION JOINLY ACROSS DATASETS
# studies are the study objects used for indivudal mapping
# targetMarkers are the markers you are testing associations on
# these contain important information such as genotype and phenotype matrices
# bgQTLs are background QTLs mapped for each dataset
# smetadata is study metada (e.g. the analyte and the name of the study)
# the column gene is actually a list of pointers to vectors containing...
# the names of the genes which are quantified in that study
# supp are the genes supported by this specific analyses
# nPerm is the number of permutations to perform
# FDR is the false discovery rate threshold
# we are interested in finding QTLs of two flavors
# genotype (geno) and genotype-analyte interaction (ixn)
# mapping allows either winners via an "effect-size" or a "-log10(p)" hit
# function return includes columns for effect-size and -log10(p) hits...
  # these are negative if a QTL does not pass one or the other
  # so none of the reported and returned QTLs can be negative in both
  # if no QTL passes one of the measurement criteria...
    # the negative value returned is the most stringent FDR cutoff tested
jointMarkerTest = function(studies,targetMarkers,unionBg,smetadata,suppGenes,
                           nPerms,FDR,fn,correctBG=FALSE,correctGRM=FALSE,
                           cores=1,chromosomes=NULL){
  
  print('PREPARING MAPPING OBJECTS AND DEALING WITH GENE SUPPORT')
  
  # Create the directory if it does not already exist
  outPath = paste0(imed,fn,'/')
  if (!dir.exists(outPath)) {
    dir.create(outPath, showWarnings = TRUE, recursive = FALSE)
  }
  
  # get true marker names
  trueMarkers = colnames(studies[[1]]$gdata)
  
  # Define default 16 S. cervisiae chromosomes if none are provided
  if (is.null(chromosomes)) {chromosomes = paste0("chr",as.roman(1:16))}
  
  # get a list of which studies support every gene
  suppList = lapply(suppGenes, function(gene) {
    studies_containing_gene = sapply(smetadata$genes,function(gene_list){
      gene %in% gene_list})
    smetadata$study[studies_containing_gene]
  }) %>% `names<-`(suppGenes)
  
  # create per-gene metadata object (gene-long list with segregants as rows)
  longMeta = mapply(function(x, y) {
    rep(y, nrow(x$gdata))
  }, studies, names(studies)) %>%
    unlist() %>%
    unname() %>%
    as.data.frame() %>%
    `names<-`(c('study')) %>%
    left_join(smetadata %>% select(study, analyte), by = "study") %>%
    mutate(geno = NA)
  
  # create gene-specific metadata using longMeta for subsetting
  longMetaByGene = lapply(names(suppList), function(gene) {
    longMeta %>% filter(study %in% suppList[[gene]])
  }) %>% `names<-`(names(suppList))
  
  # create a joint genotype object across all studies
  jointGeno = do.call(rbind, lapply(studies[
    names(studies) %in% smetadata$study], function(x) x$gdata))
  
  # subset the joint genotype object by gene support
  jointGenoByGene = lapply(suppList, function(supp) {
    return(jointGeno[which(longMeta$study %in% supp),])
  }) %>% `names<-`(names(suppList))
  
  # chrom indexing key for all markers (e.g. "ChrI")
  cc = do.call('rbind', strsplit(colnames(jointGeno), ':'))[,1]
  
  # create phenotypes object (gene long list of Seg x Chrom matrices)
  jointPhenoScaled = lapply(names(suppList), function(gene) {
    pheno_values = unlist(lapply(suppList[[gene]], function(studyName) {
      if (studyName %in% names(studies)) {
        study = studies[[studyName]]
        if (gene %in% colnames(study$pclean)) {
          return(study$pclean[, gene])
        }
      }
      return(NULL)
    }))
    
    if (!is.null(pheno_values)) {
      pheno_matrix = matrix(rep(pheno_values, length(chromosomes)), 
                            ncol = length(chromosomes), byrow = FALSE)
      colnames(pheno_matrix) = chromosomes
      rownames(pheno_matrix) = names(pheno_values)
      return(pheno_matrix)
    } else {
      NULL
    }
  }) %>% `names<-`(names(suppList))
  
  if(correctBG){
    print('CORRECTING FOR STRONG GENETIC EFFECTS ON OTHER CHROMOSOMES')
    pb = txtProgressBar(min = 0, max = length(jointGenoByGene), style = 3)
    jointPhenoScaled = mapply(function(geno, pheno, index) {
      result = adjust_for_background_QTLs(geno, pheno, unionBg, chromosomes)
      setTxtProgressBar(pb, index)
      return(result)
    }, jointGenoByGene, jointPhenoScaled,
    seq_along(jointGenoByGene), SIMPLIFY = FALSE)
    cat("\n")
  }
  
  if(correctGRM){
    print('CORRECTING FOR RELATEDNESS ON OTHER CHROMOSOMES')
    pb = txtProgressBar(min = 0, max = length(jointGenoByGene), style = 3)
    jointPhenoScaled = mapply(function(geno, pheno, index) {
      result = calculate_BLUPs(geno, pheno)
      setTxtProgressBar(pb, index)
      return(result)
    }, jointGenoByGene, jointPhenoScaled,
    seq_along(jointGenoByGene), SIMPLIFY = FALSE)
    cat("\n")
  }
  
  # scale genotypes to control for variance at each locus gene-specifically
  jointGenoByGene = lapply(jointGenoByGene,scale)
  
  print('PREPARING PERMUTATION MATICES AND FINALIZAING MAPPING PREPARATIONS')
  
  # chrom indexing key for targeted markers (e.g. "ChrI")
  ct = do.call('rbind', strsplit(targetMarkers,':'))[,1]
  
  # get study-stratified permutation matrices per gene (Perm x Seg)
  permKey = lapply(longMetaByGene,function(x)stratPerm(x$study,nPerms))
  
  # important model terms for later listing
  mterms = c('geno', 'analyte:geno')
  
  # test associations at each marker
  print('BEGINNING JOINT QTL MAPPING')
  mclapply(targetMarkers,function(marker){
    
    # get index of current marker and true marker index
    mI = which(targetMarkers == marker)
    mT = which(trueMarkers == marker)
    print(paste0('mapping associations of marker #',
                 mI,' of ',length(targetMarkers),' at position ',marker))
    
    # real associations (terms many matrices)
    rassocs = rep(list(
      list(
        gene = suppGenes,
        marker = as.character(rep(marker,length(suppGenes))),
        flavor = as.character(rep(NA,length(suppGenes))),
        effect = as.numeric(rep(NA,length(suppGenes))),
        absolute_effect = as.numeric(rep(NA,length(suppGenes))),
        p_value = as.numeric(rep(NA,length(suppGenes)))
      ) %>% as.data.frame()
    ),length(mterms)) %>% `names<-`(mterms) %>%
      imap(function(df, name) {
        df %>% mutate(flavor = name)
      })
    
    # signifiance measures
    sgMeasures = c('absolute_effect','mlog10p')
    
    # permutation assocations (mterms long Genes x Perms mats)
    passocs = rep(list(rep(list(
      matrix(NA, nrow = length(suppGenes), ncol = nPerms)), 
      length(mterms)) %>% `names<-`(mterms)),2) %>% 
      `names<-`(sgMeasures)
    
    # sweep over genes
    for(gene in seq_along(suppGenes)) {
      
      # create the object to be used for association testing
      perMeta = data.frame(
        study = longMetaByGene[[suppGenes[gene]]]$study,
        analyte = longMetaByGene[[suppGenes[gene]]]$analyte,
        geno = jointGenoByGene[[suppGenes[gene]]][,marker],
        pheno = jointPhenoScaled[[suppGenes[gene]]][,ct[mI]]
      )
      
      # run the mixed model on the real data and get n and fixed effects
      datModel = lmer(
        pheno ~ (1|study) + analyte + geno + geno:analyte, data = perMeta,
        control = lmerControl(check.conv.singular = "ignore"))
      nGI = nrow(perMeta)
      fixies = fixef(datModel) 
      names(fixies)[3:4] = mterms
      datANOVA = anova(datModel)
      
      # assign term values
      for(mterm in mterms){
        rassocs[[mterm]][gene,'effect'] = fixies[mterm]
        rassocs[[mterm]][gene,'absolute_effect'] = abs(rassocs[[
          mterm]][gene,'effect'])
        rassocs[[mterm]][gene,'samplesize'] = nGI
        rassocs[[mterm]][gene,'p_value'] = datANOVA[mterm,"Pr(>F)"]
        rassocs[[mterm]][gene,'mlog10p'] = 
          -log10(rassocs[[mterm]][gene,'p_value'])
      }
      
      # permutations for FDR
      for(p in 1:nPerms){
        
        # break genotype association using the permutation key
        perMeta$geno = perMeta$geno[permKey[[gene]][p,]]
        
        # run the mixed model on the permutation data and get fixed effects
        datModel = lmer(
          pheno ~ (1|study) + analyte + geno + geno:analyte,data = perMeta,
          control = lmerControl(check.conv.singular = "ignore"))
        fixies = fixef(datModel) 
        names(fixies)[3:4] = mterms
        
        # do ANOVA for p-values
        datANOVA = anova(datModel)
        
        # assign permutation geno and 'analyte:geno' LOD values for FDR step
        for(measure in sgMeasures){
          for(mterm in mterms){
            if(measure == 'absolute_effect'){
              mv = abs(fixies[mterm])
            } else {
              mv = -log10(datANOVA[mterm,"Pr(>F)"])
            }
            passocs[[measure]][[mterm]][gene,p] = mv
          }
        }
      }
    }
    
    # use empirical FDRs to determine QTLs
    lapply(mterms,function(mterm){
      
      # test FDR thresholds over the data for both significance measures
      tt = sapply(sgMeasures,function(measure){
        
        # eFDRs are mean perm FDs over real discoveries at data-driven threshes
        lt = sort(rassocs[[mterm]][,measure])
        rPass = sapply(lt,function(t){sum(rassocs[[mterm]][,measure] >= t)})
        pPass = sapply(lt,function(t)
          mean(apply(passocs[[measure]][[mterm]],2,function(x)sum(x >= t))))
        eFDR = pPass/rPass
        
        # remove NaNess from div by zero etc. and make FDR curve monotone
        eFDR[is.na(eFDR)] = 0
        eFDR[!is.finite(eFDR)] = 0
        eFDR = rev(cummax(rev(eFDR)))
        
        # use minimum threshold passing eFDR cutoff to write QTLs (if any)
        cuts = eFDR < FDR
        if(any(cuts)){return(lt[min(which(cuts))])} else {return(-max(lt))}
        
      }) %>% `names<-`(sgMeasures)
      
      # if any QTL will pass any measure we figure out which
      if(any(tt > 0)){
        test = sapply(sgMeasures,function(measure){
          passed = rassocs[[mterm]][,measure] >= tt[measure]
          return(vapply(passed, function(x)
            if (x > 0) tt[measure] else -tt[measure], numeric(1)))
        }) %>% as.data.frame() %>% `names<-`(paste0(sgMeasures,'_threshold'))
        
        # filter to only QTLs passing FDR for at least one measure and save
        truQTLs = cbind(rassocs[[mterm]],test) %>% 
          (function(x)x[apply(test,1,function(y)any(y > 0)),])
        write.csv(truQTLs,file = paste0(
          outPath,'QTLs_for_marker_',mT,"_of_",length(trueMarkers),
          "_at_FDR_",FDR,"_and_type_",gsub(":", "-", mterm),
          ".csv"),row.names = FALSE)
        
      }
    })
    
    return(paste('Finished',marker))
    
  },mc.cores=cores)
}

# performs assocation tests between all markers and genes 
# does not include any form of FDR correction
rawAssocs = function(studies, targetMarkers, unionBg, smetadata,
                     suppGenes, fn, correctBG=FALSE, correctGRM=FALSE,
                     cores=1, chromosomes=NULL){
  
  print('PREPARING MAPPING OBJECTS AND DEALING WITH GENE SUPPORT')
  
  # Create the directory if it does not already exist
  outPath = paste0(imed,fn,'/')
  if (!dir.exists(outPath)) {
    dir.create(outPath, showWarnings = TRUE, recursive = FALSE)
  }
  
  # get true marker positions
  trueMarkers = colnames(studies[[1]]$gdata)
  
  # Define default 16 S. cervisiae chromosomes if none are provided
  if (is.null(chromosomes)) {chromosomes = paste0("chr",as.roman(1:16))}
  
  # get a list of which studies support every gene
  suppList = lapply(suppGenes, function(gene) {
    studies_containing_gene = sapply(smetadata$genes,function(gene_list){
      gene %in% gene_list})
    smetadata$study[studies_containing_gene]
  }) %>% `names<-`(suppGenes)
  
  # create per-gene metadata object (gene-long list with segregants as rows)
  longMeta = mapply(function(x, y) {
    rep(y, nrow(x$gdata))
  }, studies, names(studies)) %>%
    unlist() %>%
    unname() %>%
    as.data.frame() %>%
    `names<-`(c('study')) %>%
    left_join(smetadata %>% select(study, analyte), by = "study") %>%
    mutate(geno = NA)
  
  # create gene-specific metadata using longMeta for subsetting
  longMetaByGene = lapply(names(suppList), function(gene) {
    longMeta %>% filter(study %in% suppList[[gene]])
  }) %>% `names<-`(names(suppList))
  
  # create a joint genotype object across all studies
  jointGeno = do.call(rbind, lapply(studies[
    names(studies) %in% smetadata$study], function(x) x$gdata))
  
  # subset the joint genotype object by gene support
  jointGenoByGene = lapply(suppList, function(supp) {
    return(jointGeno[which(longMeta$study %in% supp),])
  }) %>% `names<-`(names(suppList))
  
  # chrom indexing key for all markers (e.g. "ChrI")
  cc = do.call('rbind', strsplit(colnames(jointGeno), ':'))[,1]
  
  # create phenotypes object (gene long list of Seg x Chrom matrices)
  jointPhenoScaled = lapply(names(suppList), function(gene) {
    pheno_values = unlist(lapply(suppList[[gene]], function(studyName) {
      if (studyName %in% names(studies)) {
        study = studies[[studyName]]
        if (gene %in% colnames(study$pclean)) {
          return(study$pclean[, gene])
        }
      }
      return(NULL)
    }))
    
    if (!is.null(pheno_values)) {
      pheno_matrix = matrix(rep(pheno_values, length(chromosomes)), 
                            ncol = length(chromosomes), byrow = FALSE)
      colnames(pheno_matrix) = chromosomes
      rownames(pheno_matrix) = names(pheno_values)
      return(pheno_matrix)
    } else {
      NULL
    }
  }) %>% `names<-`(names(suppList))
  
  if(correctBG){
    print('CORRECTING FOR STRONG GENETIC EFFECTS ON OTHER CHROMOSOMES')
    pb = txtProgressBar(min = 0, max = length(jointGenoByGene), style = 3)
    jointPhenoScaled = mapply(function(geno, pheno, index) {
      result = adjust_for_background_QTLs(geno, pheno, unionBg, chromosomes)
      setTxtProgressBar(pb, index)
      return(result)
    }, jointGenoByGene, jointPhenoScaled,
    seq_along(jointGenoByGene), SIMPLIFY = FALSE)
    cat("\n")
  }
  
  if(correctGRM){
    print('CORRECTING FOR RELATEDNESS ON OTHER CHROMOSOMES')
    pb = txtProgressBar(min = 0, max = length(jointGenoByGene), style = 3)
    jointPhenoScaled = mapply(function(geno, pheno, index) {
      result = calculate_BLUPs(geno, pheno)
      setTxtProgressBar(pb, index)
      return(result)
    }, jointGenoByGene, jointPhenoScaled,
    seq_along(jointGenoByGene), SIMPLIFY = FALSE)
    cat("\n")
  }
  
  # scale genotypes to control for variance at each locus gene-specifically
  jointGenoByGene = lapply(jointGenoByGene,scale)
  
  print('FINALIZING ASSOCIATION PREPARATIONS')
  
  # chrom indexing key for targeted markers (e.g. "ChrI")
  ct = do.call('rbind', strsplit(targetMarkers,':'))[,1]
  
  # CORE ASSOCIATION LOOP - OUTER LOOP IS MARKERS
  print('BEGINNING UNFILTERED ASSOCIATIONS AT SELECT MARKERS')
  mclapply(targetMarkers,function(marker){
    
    # get index of current marker and true marker index
    mI = which(targetMarkers == marker)
    mT = which(trueMarkers == marker)
    print(paste0('mapping associations of marker #',
                 mI,' of ',length(targetMarkers),' at position ',marker))
    
    # real associations (terms many matrices)
    rassocs = rep(list(
      list(
        gene = suppGenes,
        marker = as.character(rep(marker,length(suppGenes))),
        pval = as.numeric(rep(NA,length(suppGenes))),
        effect = as.numeric(rep(NA,length(suppGenes))),
        flavor = as.character(rep(NA,length(suppGenes))),
        samplesize = as.numeric(rep(NA,length(suppGenes)))
      ) %>% as.data.frame()
    ),5) %>% `names<-`(c('analyte','geno', 'ixn','protein','mRNA')) %>%
      imap(function(df, name) {
        df %>%
          mutate(flavor = name)
      })
    
    # sweep over the genes
    for (gene in seq_along(suppGenes)) {
      
      # create the object to be used for association testing
      perMeta = data.frame(
        study = longMetaByGene[[suppGenes[gene]]]$study,
        analyte = longMetaByGene[[suppGenes[gene]]]$analyte,
        geno = jointGenoByGene[[suppGenes[gene]]][,marker],
        pheno = jointPhenoScaled[[suppGenes[gene]]][,ct[mI]]
      )
      
      # run the full model and ANOVA
      datModel = lmer(
        pheno ~ (1|study) + analyte + geno + geno:analyte, data = perMeta,
        control = lmerControl(check.conv.singular = "ignore"))
      datANOVA = anova(datModel)
      
      # run the mRNA-only model and ANOVA
      mMeta = perMeta[perMeta$analyte == "mRNA",]
      mModel = lmer(
        pheno ~ (1|study) + geno, data = mMeta,
        control = lmerControl(check.conv.singular = "ignore"))
      mANOVA = anova(mModel)
      
      # run the protein model and ANOVA
      pMeta = perMeta[perMeta$analyte == "protein",]
      pModel = lmer(
        pheno ~ (1|study) + geno, data = pMeta,
        control = lmerControl(check.conv.singular = "ignore"))
      pANOVA = anova(pModel)
      
      # assign analyte values
      rassocs$analyte[gene,'pval'] = datANOVA["analyte","Pr(>F)"]
      rassocs$analyte[gene,'effect'] = fixef(datModel)['analyteprotein']
      rassocs$analyte[gene,'samplesize'] = nrow(perMeta)
      
      # assign geno values
      rassocs$geno[gene,'pval'] = datANOVA["geno","Pr(>F)"]
      rassocs$geno[gene,'effect'] = fixef(datModel)['geno']
      rassocs$geno[gene,'samplesize'] = nrow(perMeta)
      
      # assign ixn values
      rassocs$ixn[gene,'pval'] = datANOVA["analyte:geno","Pr(>F)"]
      rassocs$ixn[gene,'effect'] = fixef(datModel)['analyteprotein:geno']
      rassocs$ixn[gene,'samplesize'] = nrow(perMeta)
      
      # assign mRNA model
      rassocs$mRNA[gene,'pval'] = mANOVA["geno","Pr(>F)"]
      rassocs$mRNA[gene,'effect'] = fixef(mModel)['geno']
      rassocs$mRNA[gene,'samplesize'] = nrow(mMeta)
      
      # assign protein model
      rassocs$protein[gene,'pval'] = pANOVA["geno","Pr(>F)"]
      rassocs$protein[gene,'effect'] = fixef(pModel)['geno']
      rassocs$protein[gene,'samplesize'] = nrow(pMeta)
      
    }
    
    # update the list of all associations
    write.csv(rassocs, file = paste0(
      outPath,"marker_",mT,"_of_",length(trueMarkers),".csv"),
      row.names = FALSE)
    
    return(paste('completed assocations for',marker))
    
  },mc.cores=cores)
  
}

# TESTS INDIVIDUAL GENE-MARKER ASSOCATIONS FOR ALL DATASETS - DEPRECATED
# takes the study object and pairIDs as arguments
pairInterrogator = function(studies, pairs) {
  
  # study names
  studyNames = names(studies)
  
  # Precompute the table size
  n_pairs = length(pairs)
  n_studies = length(studies)
  n_rows = n_pairs * n_studies
  
  # Initialize the output data frame
  results = data.frame(
    pairID = character(n_rows),
    gene = character(n_rows),
    marker = character(n_rows),
    study = character(n_rows),
    LOD = numeric(n_rows),
    estimate = numeric(n_rows),
    stringsAsFactors = FALSE
  )
  
  row_idx = 1
  
  # announce  and initialize progress bar
  print('Starting investigation of select pairs across datasets')
  pb = txtProgressBar(min = 0, max = n_pairs, style = 3)
  
  # Loop over each pair
  for (i in seq_along(pairs)) {
    pair = pairs[i]
    gene_marker = strsplit(pair, split = '_')[[1]]
    gene = gene_marker[1]
    marker = gene_marker[2]
    
    # Loop over each study
    for (study_idx in seq_along(studies)) {
      
      # Check if the study contains the gene and marker
      if (gene %in% colnames(studies[[study_idx]]$pclean) &&
          marker %in% colnames(studies[[study_idx]]$gdata)) {
        
        # Compute the correlation test
        summary = cor.test(studies[[study_idx]]$pclean[, gene],
                           studies[[study_idx]]$gdata[, marker])
        esti = summary$estimate
        LOD = -log10(summary$p.value)
        
        # Store the results in the data frame
        results[row_idx, ] = list(
          pair, gene, marker, studyNames[study_idx], LOD, esti)
        row_idx = row_idx + 1
      }
    }
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close progress bar
  close(pb)
  
  # Remove unused rows if any
  results = results[1:(row_idx - 1), ]
  
  return(results)
}

# FIGURE OUT WHICH DATASETS SUPPORT A GIVEN GENE
# takes as argument a list of genes and special metadata object
datSupportGenes = function(genes, metadata) {

  # Initialize an empty list to store the results
  result_list <- list()
  
  # Loop through each gene in the input vector
  for (gene in genes) {
    # Find which rows in the metadata contain the gene
    supporting_datasets <- metadata %>%
      filter(str_detect(genes, gene)) %>%
      pull(study)
    
    # Add the supporting datasets to the result list
    result_list[[gene]] <- supporting_datasets
  }
  
  # Return the result list
  return(result_list)
}

########################################
# ANALYSIS AND VISUALIZATION FUNCTIONS #
########################################

titler = function(s, r = NULL, n = NULL, u = TRUE) {
  # Check if replacement vectors are valid
  if (!is.null(r) & !is.null(n)) {
    if (length(r) != length(n)) {
      stop("Replacement vectors 'r' and 'n' must be of the same length.")
    }
    for (i in seq_along(r)) {
      s = str_replace_all(s, fixed(r[i]), n[i])
    }
  }
  
  # Replace underscores and periods with spaces
  s = str_replace_all(s, "[._]", " ")
  
  # Capitalize the first letter of each sentence if u is TRUE
  if (u) {
    s = str_to_sentence(s)
  }
  
  return(s)
}

# PLOT KEY TRACES FROM THE UNFILTERED ASSOCIATIONS
tracePlotter = function(dat, column, rootDir = "plots/", rootName = "plot_",
                        mc.cores = 1, r = NULL, n = NULL) {
  # Create directory if it doesn't exist
  if (!dir.exists(rootDir)) {
    dir.create(rootDir, recursive = TRUE)
  }
  
  # Create GRanges object to sort markers
  grMarkers = markRanger(unique(dat$marker)) %>% sort()
  
  # Get the start indices for each chromosome
  chr_table = table(seqnames(grMarkers))
  chrStarts = cumsum(c(0, chr_table[-length(chr_table)])) + 1
  
  # Add the final index to chrStarts
  chrStarts = c(chrStarts, length(grMarkers) + 1)
  
  # Calculate chromosome label positions
  chr_labels_pos = round(
    (chrStarts[-length(chrStarts)] + chrStarts[-1] - 1) / 2)
  chr_labels = as.roman(1:(length(chrStarts)-1))
  
  # Determine which columns are present in the data
  available_prefixes = c("geno", "ixn", "protein", "mRNA")[
    sapply(c("geno", "ixn", "protein", "mRNA"), function(x) 
      any(grepl(paste0("^", x, "\\."), colnames(dat))))]
  
  # Create column names based on available prefixes
  column_names = paste0(available_prefixes, ".", column)
  
  # Ensure correct marker ordering
  dat$marker = factor(dat$marker, levels = grMarkers$josh)
  dat = dat %>% arrange(marker)
  
  # Gather the data for plotting
  dat = dat %>%
    select(gene, marker, all_of(column_names)) %>%
    gather(key = "flavor", value = "value", -gene, -marker)
  
  # Split for plotting
  dat = split(dat, dat$gene)
  
  # Apply the function to each gene using mclapply
  mclapply(dat, function(plot_data) {
    
    # Create the plot title and axis labels using titler
    column_cap = titler(column, r, n)
    plot_data$flavor = sapply(plot_data$flavor, titler, r = r, n = n, u = FALSE)
    
    # Create the plot
    p = ggplot(plot_data, aes(x = as.numeric(marker), 
                              y = value, color = flavor, group = flavor)) +
      geom_vline(xintercept = chrStarts, 
                 linetype = "dashed", color = "black") +
      geom_line() +
      scale_color_viridis_d(option = "viridis") +
      labs(title = paste(column_cap, "trace for", plot_data$gene[1]),
           x = "Chromosome", y = column_cap, color = "Flavor") +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.text = element_text(color = "black"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            plot.title = element_text(color = "black"),
            strip.text = element_text(color = "black")) +
      scale_x_continuous(breaks = chr_labels_pos, labels = chr_labels)
    
    # Save the plot as a PDF
    ggsave(filename = paste0(rootDir, rootName, plot_data$gene[1], ".pdf"), 
           plot = p, width = 8, height = 4.5, units = "in")
  }, mc.cores = mc.cores)
}

# GET GENETIC CORRELATION SUMMARY STATS AND HISTOGRAMS ACROSS GENES
genCorFxn = function(targetAssocs, file_prefix, output_folder, infix) {
  
  targetAssocs$mRNA.AE = abs(targetAssocs$mRNA.effect)
  targetAssocs$protein.AE = abs(targetAssocs$protein.effect)
  
  # Calculate correlation coefficients for effects
  genCorsEffect = lapply(sort(unique(targetAssocs$gene)), function(gene) {
    cor.test(targetAssocs[targetAssocs$gene == gene,]$mRNA.effect,
             targetAssocs[targetAssocs$gene == gene,]$protein.effect,
             method = 'pearson')
  }) %>% `names<-`(sort(unique(targetAssocs$gene)))
  
  # Calculate correlation coefficients for AEs
  genCorsAE = lapply(sort(unique(targetAssocs$gene)), function(gene) {
    cor.test(targetAssocs[targetAssocs$gene == gene,]$mRNA.AE,
             targetAssocs[targetAssocs$gene == gene,]$protein.AE,
             method = 'pearson')
  }) %>% `names<-`(sort(unique(targetAssocs$gene)))
  
  # Extract the correlation estimates
  effect_corrs = sapply(genCorsEffect, function(x) x$estimate)
  lod_corrs = sapply(genCorsAE, function(x) x$estimate)
  
  # Create the summary table
  summary_table = data.frame(
    Statistic = c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max"),
    Effect = as.numeric(summary(effect_corrs)),
    AE = as.numeric(summary(lod_corrs))
  )
  
  # Transpose the summary table
  summary_table_t = t(summary_table)
  colnames(summary_table_t) = summary_table$Statistic
  summary_table_t = as.data.frame(summary_table_t[-1,])
  rownames(summary_table_t) = c("Effect", "AE")
  
  # Save the summary table
  write.csv(summary_table_t, file.path(
    output_folder, paste0(file_prefix, "_genet_cor_summary.csv")),
    row.names = TRUE)
  
  # Combine correlation data for plotting
  corr_data = data.frame(
    type = rep(c("Effect", "AE"), each = length(effect_corrs)),
    correlation = c(effect_corrs, lod_corrs)
  )
  
  # Plot normalized histogram
  p = ggplot(corr_data, aes(x = correlation, fill = type)) +
    geom_histogram(aes(y = ..density..), bins = 30,
                   position = "identity", alpha = 0.5) +
    scale_fill_manual(
      values = c("Effect" = "red", "AE" = "blue"), name = "Type") +
    theme_minimal() +
    labs(x = "Correlation Coefficient", y = "Density",
         title = paste("Normalized Histogram of Correlation Coefficients (",
                       infix, ")", sep = "")) +
    theme(legend.position = "bottom")
  
  # Save the plot
  ggsave(file.path(output_folder, "_genet_cor_hist.pdf"),
         plot = p, width = 7.5, height = 5)
}

# CREATE A DENSITY NORMALIZED HISTOGRAM FOR QTL VERSUS CONTROL
create_normalized_histogram = 
  function(data, variable, title_suffix, root_filename) {
    # Calculate the weights for normalization
    data = data %>%
      group_by(group) %>%
      mutate(weight = 1 / n())
    
    # Plot the normalized histograms on the same axes
    p = ggplot(data, aes(x = !!sym(variable), fill = group, weight = weight)) +
      geom_histogram(alpha = 0.6, position = 'identity',
                     bins = 30, aes(y = after_stat(density))) +
      labs(title = paste("Normalized Histogram of", title_suffix),
           x = title_suffix,
           y = "Density") +
      scale_fill_manual(values = c("Control" = "blue", "QTL" = "red"),
                        name = "Group") +
      theme_minimal()
    
    # Save the plot with a 3:2 width-to-height ratio
    plot_filename = paste0(root_filename, "_",
                           tolower(title_suffix), "_hist.png")
    ggsave(plot_filename, plot = p, width = 9, height = 6)
  }

# POINT BY POINT LINKED BOXPLOT OVERLAY PLOT
mpLinkPlotter = function(data, filename, infix, r = NULL, n = NULL) {
  data_long = data %>%
    dplyr::select(pairID, mRNA.effect, protein.effect) %>%
    tidyr::pivot_longer(cols = c(mRNA.effect, protein.effect), 
                        names_to = "effect_type", 
                        values_to = "effect_value")
  
  # Truncate effect_type labels after the period
  data_long$effect_type = sub("\\..*", "", data_long$effect_type)
  
  # Use the titler function to format the x-axis labels
  data_long$effect_type = sapply(data_long$effect_type, titler, r = r, n = n, u = FALSE)
  
  # Determine the title based on infix
  if (infix == "ixn") {
    title_infix = "of interaction QTLs"
  } else if (infix == "geno") {
    title_infix = "of genotype QTLs"
  } else {
    title_infix = infix
  }
  
  p = ggplot(data_long, aes(x = effect_type, y = effect_value, color = pairID)) +
    geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
    geom_line(aes(group = pairID), linewidth = 0.5) +
    scale_color_viridis(discrete = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "black") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title = element_text(color = "black")
    ) +
    labs(x = "Effect Type", y = "Allelic Effect", 
         title = paste("mRNA and protein effects", title_infix))
  
  # Add the infix before the file extension
  file_parts = tools::file_path_sans_ext(filename)
  new_filename = paste0(file_parts, "_", infix, ".pdf")
  
  ggsave(new_filename, plot = p, width = 7.5, height = 5)
}

# GET GENETIC CORRELATION SUMMARY STATS AND HISTOGRAMS ACROSS GENES
genCorFxn = function(targetAssocs, file_prefix, output_folder, infix) {
  
  targetAssocs$mRNA.AE = abs(targetAssocs$mRNA.effect)
  targetAssocs$protein.AE = abs(targetAssocs$protein.effect)
  
  # Calculate correlation coefficients for effects
  genCorsEffect = lapply(sort(unique(targetAssocs$gene)), function(gene) {
    cor.test(targetAssocs[targetAssocs$gene == gene,]$mRNA.effect,
             targetAssocs[targetAssocs$gene == gene,]$protein.effect,
             method = 'pearson')
  }) %>% `names<-`(sort(unique(targetAssocs$gene)))
  
  # Calculate correlation coefficients for AEs
  genCorsAE = lapply(sort(unique(targetAssocs$gene)), function(gene) {
    cor.test(targetAssocs[targetAssocs$gene == gene,]$mRNA.AE,
             targetAssocs[targetAssocs$gene == gene,]$protein.AE,
             method = 'pearson')
  }) %>% `names<-`(sort(unique(targetAssocs$gene)))
  
  # Extract the correlation estimates
  effect_corrs = sapply(genCorsEffect, function(x) x$estimate)
  lod_corrs = sapply(genCorsAE, function(x) x$estimate)
  
  # Create the summary table
  summary_table = data.frame(
    Statistic = c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max"),
    Effect = as.numeric(summary(effect_corrs)),
    AE = as.numeric(summary(lod_corrs))
  )
  
  # Transpose the summary table
  summary_table_t = t(summary_table)
  colnames(summary_table_t) = summary_table$Statistic
  summary_table_t = as.data.frame(summary_table_t[-1,])
  rownames(summary_table_t) = c("Effect", "AE")
  
  # Save the summary table
  write.csv(summary_table_t, file.path(
    output_folder, paste0(file_prefix, "_genet_cor_summary.csv")),
    row.names = TRUE)
  
  # Combine correlation data for plotting
  corr_data = data.frame(
    type = rep(c("Effect", "AE"), each = length(effect_corrs)),
    correlation = c(effect_corrs, lod_corrs)
  )
  
  # Plot normalized histogram
  p = ggplot(corr_data, aes(x = correlation, fill = type)) +
    geom_histogram(aes(y = ..density..), bins = 30,
                   position = "identity", alpha = 0.5) +
    scale_fill_manual(
      values = c("Effect" = "red", "AE" = "blue"), name = "Type") +
    theme_minimal() +
    labs(x = "Correlation Coefficient", y = "Density",
         title = paste("Normalized Histogram of Correlation Coefficients (",
                       infix, ")", sep = "")) +
    theme(legend.position = "bottom")
  
  # Save the plot
  ggsave(file.path(output_folder, "_genet_cor_hist.png"),
    plot = p, width = 9, height = 6)
}

# SUMMARY MRNA AND PROTEIN DIFFERENCES
mpSummarizer = function(data, root_filename) {
  # Create a combined dataframe in long format
  effects_list = list(
    mRNA_effect = lapply(data, function(x) x$mRNA.effect),
    protein_effect = lapply(data, function(x) x$protein.effect),
    raw_difference = lapply(data, function(x) x$mRNA.effect - x$protein.effect),
    abs_difference = lapply(data, function(x) abs(x$mRNA.effect - x$protein.effect))
  )
  
  # Convert list to long format dataframe
  combined_data = do.call(rbind, lapply(names(effects_list), function(effect_type) {
    do.call(rbind, lapply(names(effects_list[[effect_type]]), function(name) {
      data.frame(
        measure = effect_type,
        effect = effects_list[[effect_type]][[name]],
        set = name
      )
    }))
  }))
  
  # Generate colors using the viridis palette
  types = unique(combined_data$set)
  colors = viridis(length(types))
  names(colors) = types
  
  # Plot with facet_wrap
  p = ggplot(combined_data, aes(x = effect, fill = set)) +
    geom_histogram(alpha = 0.6, position = 'identity', bins = 30, aes(y = after_stat(density))) +
    labs(x = "Effect", y = "Density") +
    scale_fill_manual(values = colors, name = "QTL Type") +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title = element_text(color = "black"),
      strip.text = element_text(color = "black")
    ) +
    facet_wrap(~ measure, scales = "free", ncol = 1, labeller = labeller(measure = c(
      mRNA_effect = "Normalized Histogram of mRNA Effect",
      protein_effect = "Normalized Histogram of Protein Effect",
      raw_difference = "Normalized Histogram of Raw Differences",
      abs_difference = "Normalized Histogram of Absolute Differences"
    )))

  # Save the plot with a 3:2 width-to-height ratio
  plot_filename = paste0(root_filename, "_combined_hist.pdf")
  ggsave(plot_filename, plot = p, width = 8.5, height = 11)
  
  # Generate summary statistics
  summary_stats = combined_data %>%
    group_by(set, measure) %>%
    summarise(
      Min = min(effect),
      `1st Qu.` = quantile(effect, 0.25),
      Median = median(effect),
      Mean = mean(effect),
      `3rd Qu.` = quantile(effect, 0.75),
      Max = max(effect),
      .groups = 'drop'
    ) %>%
    pivot_wider(names_from = measure, values_from = c(Min, `1st Qu.`, Median, Mean, `3rd Qu.`, Max))
  
  # Save summary statistics
  summary_filename = paste0(root_filename, "_summary.csv")
  write.table(summary_stats, file = summary_filename, sep = ",", row.names = FALSE, quote = FALSE)
}

# TRACE LOD OR EFFECT SIZE ACROSS STUDIES INDIVIDUALLY
studyTracer = function(studies, smeta, suppGenes, fn, 
                       plotHeight, plotWidth, traceType) {
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(fn)) {
    dir.create(fn, recursive = TRUE)
  }
  
  # Get the start indices for each chromosome
  chrStarts = sapply(colnames(studies[[1]]$gdata), function(x)
    sub(":.*", "", x)) %>%
    table() %>%
    as.vector() %>%
    { x = .; x[2:length(x)] = x[1:(length(x) - 1)]; x[1] = 1; x } %>%
    cumsum()
  
  # Add the final index to chrStarts
  chrStarts = c(chrStarts, length(unique(unlist(lapply(studies, function(study) colnames(study$gdata))))) + 1)
  
  # Calculate chromosome label positions
  chr_labels_pos = round((chrStarts[-length(chrStarts)] + chrStarts[-1] - 1) / 2)
  chr_labels = as.roman(1:(length(chrStarts) - 1))
  
  # Get supporting data matrix
  suppMat = sapply(suppGenes, function(gene)
    sapply(smeta$genes, function(set)
      gene %in% set)) %>% t()
  
  colnames(suppMat) = smeta$study
  
  # Ensure you are only keeping studies in your metadata
  studies = studies[names(studies) %in% smeta$study]
  
  # Get matrices based on traceType
  sLODs = lapply(studies, function(study) {
    lod = fasterLOD(nrow(study$pheno), scale(study$pclean),
                    scale(study$gdata), betas = TRUE) 
    if (traceType == "s") {
      return(lod$LOD * sign(lod$r))
    } else if (traceType == "r") {
      return(lod$r)
    }
  })
  
  # Define a fixed set of colors for the datasets
  dataset_colors = setNames(viridis::viridis(length(unique(smeta$study))),
                            unique(smeta$study))
  
  # Initialize the progress bar
  pb = pbapply::startpb(min = 0, max = nrow(suppMat))
  
  for(i in 1:nrow(suppMat)){
    
    # Update the progress bar
    pbapply::setpb(pb, i)
    
    # Get the supported datasets for the current gene
    supported = colnames(suppMat)[suppMat[i, ]]
    
    # Prepare the plot data
    plotDatList = lapply(supported, function(dataset) {
      sLOD = sLODs[[dataset]]
      gene_row = rownames(sLOD) == rownames(suppMat)[i]
      if(any(gene_row)){
      return(data.frame(
        marker = colnames(sLOD),
        value = sLOD[gene_row, ],
        gene = rownames(suppMat)[i],
        dataset = dataset))
      }
    })
    
    # Filter out NULL data frames
    plotDatList = Filter(Negate(is.null), plotDatList)
    
    # If plotDatList is empty, continue to the next iteration
    if (length(plotDatList) == 0) {
      next
    }
    
    # Combine the non-empty data frames
    plotDat = do.call(rbind, plotDatList)
    
    # Merge with metadata to get the analyte information
    plotDat = merge(plotDat, smeta[, c("study", "analyte")],
                    by.x = "dataset", by.y = "study")
    
    # Lock in the factor levels for the dataset
    plotDat$dataset = factor(plotDat$dataset, levels = smeta$study)
    
    # Plot the results using ggplot2
    plotDat$marker = factor(plotDat$marker, levels = unique(plotDat$marker))
    plotDat$gene = factor(plotDat$gene, levels = suppGenes)
    
    plotTitle = if (traceType == "s") "Signed LOD Scores for" else 
      "R-values for"
    yAxisLabel = if (traceType == "s") "Signed LOD" else 
      "R-value"
    
    p = ggplot(plotDat, aes(
      x = as.numeric(marker), 
      y = value, color = dataset, group = dataset)) +
      geom_vline(xintercept = chrStarts, linetype = "dashed", color = "black") +
      geom_line() +
      facet_wrap(~ analyte, scales = "free_y", ncol = 1) +
      scale_color_manual(values = dataset_colors) +
      labs(title = paste(plotTitle, rownames(suppMat)[i]),
           x = "Chromosome",
           y = yAxisLabel,
           color = "Dataset") +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.box = "horizontal",
            legend.key.size = unit(0.5, "lines"),
            legend.text = element_text(size = 8, color = "black"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            plot.title = element_text(color = "black"),
            strip.text = element_text(color = "black")) +
      scale_x_continuous(breaks = chr_labels_pos, labels = chr_labels)
    
    ggsave(filename = paste0(fn, "/", rownames(suppMat)[i], ".pdf"), 
           plot = p, width = plotWidth, height = plotHeight, units = "in")
  }
  
  pbapply::closepb(pb)
}

print('SOURCE SUCCESSFUL')