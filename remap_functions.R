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
      suppressMessages(BiocManager::install(
        pkg, lib = personal_lib_path, quietly = TRUE))
    } else {
      suppressMessages(install.packages(
        pkg, lib = personal_lib_path, quietly = TRUE))
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
# parallel       - for parallelization tools (e.g. mclapply)
# patchwork      - for the nice three-panel plots
# pbapply        - for nice progress bars on apply functions
# plyr           - needed for Brem 02 processing
# progress       - for progress bars
# readxl         - for the MS Excel files
# rrBLUP         - to deal with population structure
# tidyverse      - for the good good (data wrangling, plots, pipes, etc.)
# viridis        - for the nicely colored plots

##################################################
# DECIDE ON CONSTANTS AND LOAD UNIVERSAL OBJECTS #
##################################################

# Prepare directories and set the working directory
wd = "/home/albertf/vandy147/p05/"
raw = paste0(wd,"raw/")
imed = paste0(wd,'imed/')
asys = paste0(wd,'asys/')
logs = paste0(wd,'logs/')
setwd(wd)

# Create directories if they do not exist
for (dir in c(imed,asys,logs)) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}

# Choose a random seed for the script
set.seed(464788913)

# Load the universally used metadata object
studyMetadata = read.delim(file=paste0(raw,'studyMetadata.csv'),
                           header=TRUE,sep=',') %>%
  as.data.frame()

# Create the non-bulk subset
classicMetadata = studyMetadata[
  studyMetadata$mapping_approach == 'QTL',]

# Further subset to background A in synthetic complete
scaMetadata = classicMetadata %>%
  dplyr::filter(medium == "SC", genetic_background == "A")

# Study root names
studyList = c("Brem_2002", "Brem_2005", "SmithEth_2008", "SmithGlc_2008",
              "Khan_2009", "Foss_2011", "Picotti_2013", "Albert_2018", 
              "AlbertAlt_2018", "GrossbachPRT_2022", "GrossbachRNA_2022")

studyTypes = c("mRNA","mRNA","special","special",
               "protein","protein","protein",
               "mRNA","mRNA","protein","mRNA")

# Threshold for NA removal
threshNA = 0.25

#########################
# GENERAL FUNCTION BANK #
#########################

# Negation of in
`%nin%` = Negate(`%in%`)

# Ensures a directory exists by creating it if it doesn't
ensure_directory = function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}


# Scale without a bunch of extra data
rawscale = function(m){
  scale(m) %>%
    {attr(., "scaled:center") <- NULL; .} %>%
    {attr(., "scaled:scale") <- NULL; .}%>%
    return()
}

# Round up CSVs into a nice R object
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

# Replace anything containing a specific pattern with that pattern
unifyPattern = function(v,p,s=F) {
  v[grepl(p, v, ignore.case = !s)] = p
  return(v)
}

# Where p is a vector of replacement patterns in this case
unifyChain = function(v,p,s=F){
  for(i in p){
    v = unifyPattern(v,i,s)
  }
  return(v)
}

# Replaces elements of a vector in the keys of a dictionary with values
dictReplace = function(vec, df){
  plyr::mapvalues(vec, from = df[[1]], to = df[[2]], warn_missing = FALSE)}

# Averages replicates based on a key system
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

# Removes rows or columns too heavy on missing values
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

# Alternate over rows and columns until rmHeavyNA is met on both margins
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
    
    # Flag to break out of the loop
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

# Get strain information from Brem .soft files
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

# Processes the SOFT files for Brem 2002 and 2005
cleanGDS = function(dir, filename, pat, naThreshold, feats, old = F, new = F) {
  
  strainKey = getStrainInfo(paste0(dir, filename))
  
  strainKey$strain = strainKey$strain %>%
    (function(x) sapply(x, function(y) sub(pat, "", y))) %>%
    unifyChain(c('BY', 'RM'))
  
  # Replace any characters if necessary 
  if (is.character(old) & is.character(new)) {
    strainKey$strain = gsub(old, new, strainKey$strain)
  }
  
  datGDS = GEOquery::getGEO(filename = paste0(dir, filename)) %>%
    GEOquery::Table()
  
  # Normalize using the blanks
  datGDS[, 3:ncol(datGDS)] = apply(
    as.matrix(datGDS[, -c(1:2)]), 2, function(x) {
      (x - mean(x[which(datGDS$IDENTIFIER == 'blank')], na.rm = T))
    })
  
  # Replace gene names with ORF names
  datGDS$IDENTIFIER = dictReplace(
    datGDS$IDENTIFIER, feats[c('Standard_Gene_Name', 'Feature_Name')]
  )
  
  # Get rid of any identifiers not associated with ORFs
  datGDS = datGDS[datGDS$IDENTIFIER %in% feats$Feature_Name,]
  
  # Remove rows with too many NAs
  datGDS = rmIterNA(datGDS, naThreshold)
  
  # Merge duplicates by average and then impute missing data by averages
  dKey = datGDS$IDENTIFIER
  theNAs = is.na(datGDS[,-c(1:2)])
  print(paste('There are',sum(theNAs),'NAs of',length(theNAs),'datapoints'))
  datGDS = averageNA(datGDS[, -c(1:2)])
  print('Summary of non-imputed values')
  print(summary(datGDS[!theNAs]))
  print('Summary of imputed values')
  print(summary(datGDS[theNAs]))
  
  # Remove gene duplicates by averaging relevant values
  datGDS = meanDedupe(datGDS, dKey)
  
  # Ensure strainKey is in order andremove strain duplicates by averaging 
  common_elements = intersect(strainKey$gsmID, colnames(datGDS))
  strainKey = strainKey[strainKey$gsmID %in% common_elements, ]
  strainKey = strainKey[match(colnames(datGDS), strainKey$gsmID), ]
  datGDS = meanDedupe(t(datGDS),strainKey$strain)
  
  # Return the numeric elements as a transposed matrix
  return(datGDS)
  
}

# Merges the two dye-swap matrices
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
  
  # Get the number of rows that may contain NA runs
  num_rows = nrow(geno)
  
  # Holds NA chain lengths for rows
  na_chain_lengths = vector("list", num_rows)
  
  # Check each row of NAs
  for (row in 1:num_rows) {
    # Initialize the NA chain length for the row and walk
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

# Simple end-matching genome imputation by average
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
          # Else assign the matching value
          geno[i, j] = ifelse(left_value != right_value, 0, left_value)
        }
      }
    }
  }
  
  # Return the modified geno matrix
  return(geno)
}

# Replaces NA values with row or column averages
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

# Turns a list of markers in Josh's format to granges
# Tan handle if there is variant underscore separation or not
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

# Wrapper for the pairwise intersection from genomic ranges
# Allows for the quick addition of padding on the first grange
# E.g. bigger flanks to reduce stringency
# Add padding to gr1, pairwise intersect, and return
padded_pintersect = function(gr1, gr2, pad=0) {
  gr1_padded = resize(gr1, width(gr1) + 2 * pad, fix = "center")
  intersections = pintersect(gr1_padded, gr2)
  return(intersections)
}

# Harmonizes matrix elements of a list by segregant representation
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

# Create a one-hot encoding for any colummns that are factors or characters
oneHotCo = function(covars) {
  
  # Identify columns that are either factors or characters
  categorical_columns = sapply(covars, function(x)
    is.factor(x) || is.character(x))
  covars_encoded = covars
  
  # Encode
  for (col_name in names(covars)[categorical_columns]) {
    covar_factor = as.factor(covars[[col_name]])
    dummy_vars = model.matrix(~ covar_factor - 1)
    colnames(dummy_vars) = paste(col_name, levels(covar_factor), sep = "_")
    covars_encoded = cbind(covars_encoded, dummy_vars)
    covars_encoded[[col_name]] = NULL
  }
  
  # Convert to matrix and return
  covars_encoded = as.matrix(covars_encoded)
  return(covars_encoded)
}

###############################################
# INDIVIDUAL MAPPING LEGACY CODE (JOSH BLOOM) #
###############################################

# QTL mapping dependency
calc.BLUPS= function(G,Z,Vinv,y,X,B){G%*%crossprod(Z,Vinv)%*%(y-X%*%B)}

# QTL mapping dependency
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

# QTL mapping dependency
# Borrowed and modified from from mixed.solve() in rrBLUP 08/02/15
# Super optimized for one VC and fixed effects 
# Roughly 1000X speedup by precomputing eigen decomp
m.S=function (y, K = NULL, bounds = c(1e-09, 1e+09),theta=NULL,Q=NULL,X=NULL) 
{
  n = length(y)
  y = matrix(y, n, 1)
  if(is.null(X)) {p = 1} else {p = ncol(X)}
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

# Faster LOD calculation
# No NA handling but very fast
# Input is scaled phenotypes and genotypes
fasterLOD = function(n.pheno, pheno.s,gdata.s, betas=FALSE, sdx=1, pheno=NULL){
  r=crossprod(pheno.s, gdata.s)/(n.pheno-1)
  LOD=(-n.pheno*log(1-r^2))/(2*log(10))
  if(betas==FALSE) {
    return(LOD)
  } else {
    return(list(r=r, LOD=LOD))
  }
}

# Finds background QTLs
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
  
  # If no chromosome is provided, default to the 16 yeast chromosomes
  if(is.null(chromosomes)) { chromosomes=paste0('chr', as.roman(1:16)) } 
  
  # This vector tells you which chromosome each marker belongs to
  cvec=(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])
  
  # Iteration 1  (find background QTL effects)
  LODS=fasterLOD(nrow(pheno.scaled),pheno.scaled,gdata.scaled)
  
  # Lists to hold LODs, markers, and scaled markers by chromosome
  LODS.by.chr=list()
  gdata.s.by.chr=list()
  gdata.by.chr=list()
  for(cc in chromosomes) {   
    LODS.by.chr[[cc]]=LODS[,which(cvec %in% cc)]   
    gdata.s.by.chr[[cc]]=gdata.scaled[,which(cvec %in% cc)]   
    gdata.by.chr[[cc]]=gdata[,which(cvec %in% cc)]   
  }
  
  # Max LOD per chromosome and then indices for them
  max.lod.chr=lapply(LODS.by.chr, function(x) { apply(x,1,max) })
  max.lod.chr.ind=lapply(LODS.by.chr, function(x) { apply(x,1,which.max) })
  
  # Pulls out the marker index (by chromosome) of the peak marker 
  pmarker.chr.ind =  mapply(function(x,y,z) {
    
    # Index gets included if LOD exceeds background threshold
    h=y[x>bck.thresh]
    
    # Retrieve marker names using indices from h; i holds the actual names of 
    # The markers at indices where LOD scores are significant.
    i=colnames(z)[h]
    
    # Assign names to i using the original indices in h. This step links each
    # Marker name with its corresponding LOD score index to trace it
    # Makes clear which LOD scores each marker name corresponds to.
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
    
    # If the gene not present in the QTL list return unchanged 
    if(!(i %in% names(pbg)) ) {return(pheno.scaled[,i]) } 
    
    # If gene present remove effect of gdata
    else{ residuals(lm(phenotype.raw[,i]~covariates+gdata[,pbg[[i]]]) )  }
  })
  
  # Combine residuals of all traits into a single matrix 
  # This matrix has the same dimensions as pheno.scaled
  p.resid1=do.call('cbind', plist)
  
  # Assign original phenotype names to the columns of the residuals matrix
  colnames(p.resid1)=colnames(pheno.scaled)
  
  # Normalize residuals to have zero mean and unit variance for further analysis
  pheno.scaled=scale(p.resid1)
  
  # Iteration 2 (find background QTL effects)
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
  
  # Iteration 3 (find background QTL effects)
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

# Maps QTLs optional background QTL and GRM corrections
# Covariates is a matrix of the covaraites
# Background QTLs are an input for background QTL correction
# phenotype.raw is an unsacaled phenotypes matrix (segregant x trait)
# Elements are numeric measures
# Note pheno.scaled is the scaled pheno matrix
# gdata is an unscaled genotypes matrix (segregant x marker)
# Elements are typically binary -1 or 1 coded genotypes but can be numeric
# gdata.scaled is the scaled genotypes matrix
mapQTL=function(covariates, background.QTL, phenotype.raw, pheno.scaled, 
                gdata, gdata.scaled, 
                chromosomes=NULL,
                n.perm=1000, FDR.thresh=.05, Arandom=TRUE) {
  
  # If no chromosome is provided, default to the 16 yeast chromosomes
  if(is.null(chromosomes)) { chromosomes=paste0('chr', as.roman(1:16)) } 
  
  # This vector tells you which chromosome each marker belongs to
  cvec=(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])
  
  # Get gdata and scaled gdata by chromosome
  gdata.s.by.chr=list()
  gdata.by.chr=list()
  for(cc in chromosomes) {   
    gdata.s.by.chr[[cc]]=gdata.scaled[,which(cvec %in% cc)]   
    gdata.by.chr[[cc]]=gdata[,which(cvec %in% cc)]   
  }
  
  # Initialize peak holding list
  peakList=list()
  
  # Iterate over the chromosomes in the phenotype matrix
  for(cc in chromosomes) {  
    
    # Tell us what chromosome we are on
    print(cc) 
    
    # Iterate over the genes in the phenotype matrix
    plist=lapply(colnames(phenotype.raw),function(i){ 
      
      # If there is no QTL for that gene
      if(!(i %in% names(background.QTL[[cc]]))){
        
        # Then return only the covariate correction
        return(residuals(lm(phenotype.raw[,i]~covariates)))}
      
      else{
        
        # Else return correction after 
        residuals(lm(
          phenotype.raw[,i]~covariates+gdata[,background.QTL[[cc]][[i]]]))}
    })
    
    # preal will be the same size as the phenotype matrix
    # preal combines the residual slices
    preal=do.call('cbind', plist)
    
    # Colnames for preal are lifted directly from the phenotype matrix
    colnames(preal)=colnames(phenotype.raw)
    
    # presid is made of preal
    presid = preal
    
    # Remove residual additive effect of genetics
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
    
    # This flag forms the basis of a while loop
    findingPeaks=TRUE
    jump=1
    
    # While we remain in peakfinding mode we scan
    while(findingPeaks==TRUE) {
      
      # We find the maximum LODs
      max.obsLOD=apply(obsLOD$LOD,1,max)
      
      # We retrieve the indices (ind) for the maximum LODs
      max.ind.obsLOD=apply(obsLOD$LOD,1,which.max)   
      
      # We prepare the LOD permutation matrix
      permLOD=matrix(0, length(max.obsLOD), n.perm)
      
      # Announce permutations and start progress bar
      print('doing permutations')
      pb = txtProgressBar(min = 1, max =n.perm, style = 3)
      
      # Loop for the number of permutations
      for(i in 1:n.perm) {
        setTxtProgressBar(pb, i)
        permLOD[,i]=apply(fasterLOD(nrow(presid),presid[sample(nrow(presid)),],
                                    gdata.s.by.chr[[cc]]),1,max)
      }
      close(pb)
      
      # Check how many are above the iterated LOD threshold
      obsPcnt = sapply(seq(1.5,9,.05),
                       function(thresh){sum(max.obsLOD>thresh)})
      names(obsPcnt) = seq(1.5, 9, .05)
      
      # If there are no more than 2 QTLs discovered, break the loop
      # We need to discover at least three to say that it's real
      if(sum(obsPcnt)<3) {break}
      
      # Take the means of the permuted LODs to get the expected number of FDs
      expPcnt = sapply(seq(1.5, 9, .05),  
                       function(thresh) { 
                         mean(apply(permLOD,2,function(ll){sum(ll>thresh)}))
                       } )
      names(expPcnt) = seq(1.5, 9, .05)
      
      # FDR is then just the expected percents divided by the observed
      pFDR = expPcnt/obsPcnt
      
      # Any numerical weirdness due to div by 0 or missingness is zeroed here
      pFDR[is.na(pFDR)]=0
      pFDR[!is.finite(pFDR)]=0
      
      # Ensure monotonicity in FDR curve by cumulative max on double reverse
      pFDR = rev(cummax(rev(pFDR)))
      
      # Linear interpolation of function
      fdrFX=approxfun(pFDR, seq(1.5,9,.05))
      thresh=fdrFX(FDR.thresh)
      print(paste('FDR Thresh 5%', thresh)) 
      peak.gene.ind=max.ind.obsLOD[which(max.obsLOD>thresh)]
      
      # Continue 
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
      } else {findingPeaks=FALSE} # stop peakfinding if we find none
    }
  }
  return(peakList)
}

# LOD from n and R
# Inputs: sample size (n) and correlation coefficient (r)
# Outputs: LOD score
nrLOD = function(n, r){
  return(-n*log(1-r^2))/(2*log(10))
}

# Wrapper for getting background QTLs
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

# Wrapper for mapping QTLs using Josh's pipeline over a list
fullMap = function(studies, bgQTLs, fn, method = 'f', cores = 1) {
  
  # EXPLANATION OF MAPPING TYPES
  # Full map (f): includes empirical null, pop struc crxn, background crxn
  # Simple map (s): includes empirical null
  
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
      print(paste0('Files already exist for ', study_name,
                   ' using method ', method))
      print('Loading files now...')
      load(paste0(fn, study_name, typeRoot, '.rda'))
      qtlResult = get(paste0(study_name, typeRoot))
    } else {
      print(paste0('Attempting to map QTLs for ', study_name,
                   ' using method ', method))
      
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
            combinedQTLs[[gene]][[chr]] = unique(c(
              combinedQTLs[[gene]][[chr]], markers))
          }
        }
      }
    }
  }
  
  return(combinedQTLs)
}

# Check which genes are represented in two mRNA and two protein datasets
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

# Do statified permutations (e.g. for study)
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

# Map association between markers and expression jointly across datasets
# Studies are the study objects used for indivudal mapping
# targetMarkers are the markers you are testing associations on
# bgQTLs are background QTLs mapped for each dataset
# smetadata is study metada (e.g. the analyte and the name of the study)
# The column gene is actually a list of pointers to vectors containing...
# The names of the genes which are quantified in that study
# supp are the genes supported by this specific analyses
# nPerm is the number of permutations to perform
# FDR is the false discovery rate threshold
# We are interested in finding QTLs of two flavors
# Genotype (geno) and genotype-analyte interaction (ixn)
# Mapping allows either winners via an "effect-size" or a "-log10(p)" hit
# Function return includes columns for effect-size and -log10(p) hits...
# These are negative if a QTL does not pass one or the other
# So none of the reported and returned QTLs can be negative in both
# If no QTL passes one of the measurement criteria...
# The negative value returned is the most stringent FDR cutoff tested
# Bit of a vile memory hog that I did not optimize for memory-usage
jointMarkerTest = function(studies,targetMarkers,unionBg,smetadata,suppGenes,
                           nPerms,FDR,fn,correctBG=FALSE,correctGRM=FALSE,
                           cores=1,chromosomes=NULL){
  
  print('PREPARING MAPPING OBJECTS AND DEALING WITH GENE SUPPORT')
  
  # Create the directory if it does not already exist
  outPath = paste0(imed,fn,'/')
  if (!dir.exists(outPath)) {
    dir.create(outPath, showWarnings = TRUE, recursive = FALSE)
  }
  
  # Get true marker names
  trueMarkers = colnames(studies[[1]]$gdata)
  
  # Define default 16 S. cervisiae chromosomes if none are provided
  if (is.null(chromosomes)) {chromosomes = paste0("chr",as.roman(1:16))}
  
  # Get a list of which studies support every gene
  suppList = lapply(suppGenes, function(gene) {
    studies_containing_gene = sapply(smetadata$genes,function(gene_list){
      gene %in% gene_list})
    smetadata$study[studies_containing_gene]
  }) %>% `names<-`(suppGenes)
  
  # Create per-gene metadata object (gene-long list with segregants as rows)
  longMeta = mapply(function(x, y) {
    rep(y, nrow(x$gdata))
  }, studies, names(studies)) %>%
    unlist() %>%
    unname() %>%
    as.data.frame() %>%
    `names<-`(c('study')) %>%
    left_join(smetadata %>% select(study, analyte), by = "study") %>%
    mutate(geno = NA)
  
  # Create gene-specific metadata using longMeta for subsetting
  longMetaByGene = lapply(names(suppList), function(gene) {
    longMeta %>% filter(study %in% suppList[[gene]])
  }) %>% `names<-`(names(suppList))
  
  # Create a joint genotype object across all studies
  jointGeno = do.call(rbind, lapply(studies[
    names(studies) %in% smetadata$study], function(x) x$gdata))
  
  # Subset the joint genotype object by gene support
  jointGenoByGene = lapply(suppList, function(supp) {
    return(jointGeno[which(longMeta$study %in% supp),])
  }) %>% `names<-`(names(suppList))
  
  # Chrom indexing key for all markers (e.g. "ChrI")
  cc = do.call('rbind', strsplit(colnames(jointGeno), ':'))[,1]
  
  # Create phenotypes object (gene long list of Seg x Chrom matrices)
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
  
  # Scale genotypes to control for variance at each locus gene-specifically
  jointGenoByGene = lapply(jointGenoByGene,scale)
  
  print('PREPARING PERMUTATION MATICES AND FINALIZAING MAPPING PREPARATIONS')
  
  # Chrom indexing key for targeted markers (e.g. "ChrI")
  ct = do.call('rbind', strsplit(targetMarkers,':'))[,1]
  
  # Get study-stratified permutation matrices per gene (Perm x Seg)
  permKey = lapply(longMetaByGene,function(x)stratPerm(x$study,nPerms))
  
  # Important model terms for later listing
  mterms = c('geno', 'analyte:geno')
  
  # Test associations at each marker
  print('BEGINNING JOINT QTL MAPPING')
  mclapply(targetMarkers,function(marker){
    
    # Get index of current marker and true marker index
    mI = which(targetMarkers == marker)
    mT = which(trueMarkers == marker)
    print(paste0('mapping associations of marker #',
                 mI,' of ',length(targetMarkers),' at position ',marker))
    
    # Real associations (terms many matrices)
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
    
    # Signifiance measures
    sgMeasures = c('absolute_effect','mlog10p')
    
    # Permutation assocations (mterms long Genes x Perms mats)
    passocs = rep(list(rep(list(
      matrix(NA, nrow = length(suppGenes), ncol = nPerms)), 
      length(mterms)) %>% `names<-`(mterms)),2) %>% 
      `names<-`(sgMeasures)
    
    # Sweep over genes
    for(gene in seq_along(suppGenes)) {
      
      # Create the object to be used for association testing
      perMeta = data.frame(
        study = longMetaByGene[[suppGenes[gene]]]$study,
        analyte = longMetaByGene[[suppGenes[gene]]]$analyte,
        geno = jointGenoByGene[[suppGenes[gene]]][,marker],
        pheno = jointPhenoScaled[[suppGenes[gene]]][,ct[mI]]
      )
      
      # Run the mixed model on the real data and get n and fixed effects
      datModel = lmer(
        pheno ~ (1|study) + analyte + geno + geno:analyte, data = perMeta,
        control = lmerControl(check.conv.singular = "ignore"))
      nGI = nrow(perMeta)
      fixies = fixef(datModel) 
      names(fixies)[3:4] = mterms
      datANOVA = anova(datModel)
      
      # Assign term values
      for(mterm in mterms){
        rassocs[[mterm]][gene,'effect'] = fixies[mterm]
        rassocs[[mterm]][gene,'absolute_effect'] = abs(rassocs[[
          mterm]][gene,'effect'])
        rassocs[[mterm]][gene,'samplesize'] = nGI
        rassocs[[mterm]][gene,'p_value'] = datANOVA[mterm,"Pr(>F)"]
        rassocs[[mterm]][gene,'mlog10p'] = 
          -log10(rassocs[[mterm]][gene,'p_value'])
      }
      
      # Permutations for FDR
      for(p in 1:nPerms){
        
        # Break genotype association using the permutation key
        perMeta$geno = perMeta$geno[permKey[[gene]][p,]]
        
        # Run the mixed model on the permutation data and get fixed effects
        datModel = lmer(
          pheno ~ (1|study) + analyte + geno + geno:analyte,data = perMeta,
          control = lmerControl(check.conv.singular = "ignore"))
        fixies = fixef(datModel) 
        names(fixies)[3:4] = mterms
        
        # Do ANOVA for p-values
        datANOVA = anova(datModel)
        
        # Assign permutation geno and 'analyte:geno' LOD values for FDR step
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
    
    # Use empirical FDRs to determine QTLs
    lapply(mterms,function(mterm){
      
      # Test FDR thresholds over the data for both significance measures
      tt = sapply(sgMeasures,function(measure){
        
        # eFDRs are mean perm FDs over real discoveries at data-driven threshes
        lt = sort(rassocs[[mterm]][,measure])
        rPass = sapply(lt,function(t){sum(rassocs[[mterm]][,measure] >= t)})
        pPass = sapply(lt,function(t)
          mean(apply(passocs[[measure]][[mterm]],2,function(x)sum(x >= t))))
        eFDR = pPass/rPass
        
        # Remove NaNess from div by zero etc. and make FDR curve monotone
        eFDR[is.na(eFDR)] = 0
        eFDR[!is.finite(eFDR)] = 0
        eFDR = rev(cummax(rev(eFDR)))
        
        # Use minimum threshold passing eFDR cutoff to write QTLs (if any)
        cuts = eFDR < FDR
        if(any(cuts)){return(lt[min(which(cuts))])} else {return(-max(lt))}
        
      }) %>% `names<-`(sgMeasures)
      
      # If any QTL will pass any measure we figure out which
      if(any(tt > 0)){
        test = sapply(sgMeasures,function(measure){
          passed = rassocs[[mterm]][,measure] >= tt[measure]
          return(vapply(passed, function(x)
            if (x > 0) tt[measure] else -tt[measure], numeric(1)))
        }) %>% as.data.frame() %>% `names<-`(paste0(sgMeasures,'_threshold'))
        
        # Filter to only QTLs passing FDR for at least one measure and save
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

# Performs assocation tests between all markers and genes 
# essentially a truncated and reformatted jointMarkerTest
# Does not include any form of FDR correction
# Also a memory-unoptimized hog
rawAssocs = function(studies, targetMarkers, unionBg, smetadata,
                     suppGenes, fn, correctBG=FALSE, correctGRM=FALSE,
                     cores=1, chromosomes=NULL){
  
  print('PREPARING MAPPING OBJECTS AND DEALING WITH GENE SUPPORT')
  
  # Create the directory if it does not already exist
  outPath = paste0(imed,fn,'/')
  if (!dir.exists(outPath)) {
    dir.create(outPath, showWarnings = TRUE, recursive = FALSE)
  }
  
  # Get true marker positions
  trueMarkers = colnames(studies[[1]]$gdata)
  
  # Define default 16 S. cervisiae chromosomes if none are provided
  if (is.null(chromosomes)) {chromosomes = paste0("chr",as.roman(1:16))}
  
  # Get a list of which studies support every gene
  suppList = lapply(suppGenes, function(gene) {
    studies_containing_gene = sapply(smetadata$genes,function(gene_list){
      gene %in% gene_list})
    smetadata$study[studies_containing_gene]
  }) %>% `names<-`(suppGenes)
  
  # Create per-gene metadata object (gene-long list with segregants as rows)
  longMeta = mapply(function(x, y) {
    rep(y, nrow(x$gdata))
  }, studies, names(studies)) %>%
    unlist() %>%
    unname() %>%
    as.data.frame() %>%
    `names<-`(c('study')) %>%
    left_join(smetadata %>% select(study, analyte), by = "study") %>%
    mutate(geno = NA)
  
  # Create gene-specific metadata using longMeta for subsetting
  longMetaByGene = lapply(names(suppList), function(gene) {
    longMeta %>% filter(study %in% suppList[[gene]])
  }) %>% `names<-`(names(suppList))
  
  # Create a joint genotype object across all studies
  jointGeno = do.call(rbind, lapply(studies[
    names(studies) %in% smetadata$study], function(x) x$gdata))
  
  # Subset the joint genotype object by gene support
  jointGenoByGene = lapply(suppList, function(supp) {
    return(jointGeno[which(longMeta$study %in% supp),])
  }) %>% `names<-`(names(suppList))
  
  # Chrom indexing key for all markers (e.g. "ChrI")
  cc = do.call('rbind', strsplit(colnames(jointGeno), ':'))[,1]
  
  # Create phenotypes object (gene long list of Seg x Chrom matrices)
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
  
  # Scale genotypes to control for variance at each locus gene-specifically
  jointGenoByGene = lapply(jointGenoByGene,scale)
  
  print('FINALIZING ASSOCIATION PREPARATIONS')
  
  # Chrom indexing key for targeted markers (e.g. "ChrI")
  ct = do.call('rbind', strsplit(targetMarkers,':'))[,1]
  
  # oop for all markers
  print('BEGINNING UNFILTERED ASSOCIATIONS AT SELECT MARKERS')
  mclapply(targetMarkers,function(marker){
    
    # Get index of current marker and true marker index
    mI = which(targetMarkers == marker)
    mT = which(trueMarkers == marker)
    print(paste0('mapping associations of marker #',
                 mI,' of ',length(targetMarkers),' at position ',marker))
    
    # Real associations (terms many matrices)
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
    
    # Sweep over the genes
    for (gene in seq_along(suppGenes)) {
      
      # Create the object to be used for association testing
      perMeta = data.frame(
        study = longMetaByGene[[suppGenes[gene]]]$study,
        analyte = longMetaByGene[[suppGenes[gene]]]$analyte,
        geno = jointGenoByGene[[suppGenes[gene]]][,marker],
        pheno = jointPhenoScaled[[suppGenes[gene]]][,ct[mI]]
      )
      
      # Run the full model and ANOVA
      datModel = lmer(
        pheno ~ (1|study) + analyte + geno + geno:analyte, data = perMeta,
        control = lmerControl(check.conv.singular = "ignore"))
      datANOVA = anova(datModel)
      
      # Run the mRNA-only model and ANOVA
      mMeta = perMeta[perMeta$analyte == "mRNA",]
      mModel = lmer(
        pheno ~ (1|study) + geno, data = mMeta,
        control = lmerControl(check.conv.singular = "ignore"))
      mANOVA = anova(mModel)
      
      # Run the protein model and ANOVA
      pMeta = perMeta[perMeta$analyte == "protein",]
      pModel = lmer(
        pheno ~ (1|study) + geno, data = pMeta,
        control = lmerControl(check.conv.singular = "ignore"))
      pANOVA = anova(pModel)
      
      # Assign analyte values
      rassocs$analyte[gene,'pval'] = datANOVA["analyte","Pr(>F)"]
      rassocs$analyte[gene,'effect'] = fixef(datModel)['analyteprotein']
      rassocs$analyte[gene,'samplesize'] = nrow(perMeta)
      
      # Assign geno values
      rassocs$geno[gene,'pval'] = datANOVA["geno","Pr(>F)"]
      rassocs$geno[gene,'effect'] = fixef(datModel)['geno']
      rassocs$geno[gene,'samplesize'] = nrow(perMeta)
      
      # Assign ixn values
      rassocs$ixn[gene,'pval'] = datANOVA["analyte:geno","Pr(>F)"]
      rassocs$ixn[gene,'effect'] = fixef(datModel)['analyteprotein:geno']
      rassocs$ixn[gene,'samplesize'] = nrow(perMeta)
      
      # Assign mRNA model
      rassocs$mRNA[gene,'pval'] = mANOVA["geno","Pr(>F)"]
      rassocs$mRNA[gene,'effect'] = fixef(mModel)['geno']
      rassocs$mRNA[gene,'samplesize'] = nrow(mMeta)
      
      # Assign protein model
      rassocs$protein[gene,'pval'] = pANOVA["geno","Pr(>F)"]
      rassocs$protein[gene,'effect'] = fixef(pModel)['geno']
      rassocs$protein[gene,'samplesize'] = nrow(pMeta)
      
    }
    
    # Update the list of all associations
    write.csv(rassocs, file = paste0(
      outPath,"marker_",mT,"_of_",length(trueMarkers),".csv"),
      row.names = FALSE)
    
    return(paste('completed assocations for',marker))
    
  },mc.cores=cores)
  
}

# Figures out which datasets support specific genes
# Takes as argument a list of genes and special metadata object
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

# Formats string nicely for titles, plot legends, etc.
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

# Plots histograms of the raw and clean datasets on a study-by-study basis
phenoHister = function(studies, dir, trueNames = NULL) {
  
  # Use the names of studies if trueNames is not provided
  if (is.null(trueNames)) {
    trueNames = names(studies)
  }
  
  # Check that trueNames is the same length as studies
  if (length(trueNames) != length(studies)) {
    stop("trueNames must be the same length as studies")
  }
  
  # Create the directory if it doesn't exist
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  # Iterate over each study and corresponding true name
  for (i in seq_along(studies)) {
    study = studies[[i]]
    true_name = trueNames[i]
    
    # Create data frames for raw and processed data with a 'type' column
    df_combined = rbind(
      data.frame(mean_values = colMeans(study$pheno),
                 type = factor("Raw", levels = c("Raw", "Processed"))),
      data.frame(mean_values = colMeans(study$pclean),
                 type = factor("Processed", levels = c("Raw", "Processed"))))
    
    # Create the plot with facet wrap
    p_combined = ggplot(df_combined, aes(x = mean_values)) +
      geom_histogram(bins = 30, color = "black") + 
      facet_wrap(~type, scales = "free", ncol = 1) +  
      labs(
        title = true_name,
        x = "Mean expression of gene across segregants",
        y = "Frequency"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)  
      ) +
      scale_x_continuous(labels = function(x) sprintf("%.2f", x)) 
    
    file_name = gsub(" ", "_", true_name) 
    full_path = file.path(dir, paste0(file_name, "_phenotype_histograms.pdf"))
    ggsave(filename = full_path, plot = p_combined, width = 7.5,
           height = 5, limitsize = FALSE)
  }
}



# Plots key traces from unfiltered assocations
tracePlotter = function(dat, valueType, rootDir = "plots/", 
                        rootName = "plot_",plotWidth=7.5,plotHeight=5,
                        r = NULL, n = NULL) {
  
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
  chr_labels = as.roman(1:(length(chrStarts) - 1))
  
  # Define prefixes and order them
  all_prefixes = c("mRNA", "protein", "genotype", "interaction", "analyte")
  
  # Determine which prefixes are present in the data
  present_prefixes = all_prefixes[sapply(all_prefixes, function(x)
    any(x == colnames(dat)))]
  
  # Ensure correct marker ordering
  dat$marker = factor(dat$marker, levels = grMarkers$josh)
  dat = dat %>% arrange(marker)
  
  # Gather the data for plotting
  dat = dat %>%
    select(gene, marker, all_of(present_prefixes)) %>%
    gather(key = "flavor", value = "value", -gene, -marker)
  
  # Set the order of the flavor factor for plotting and legend
  # reverse plotting order but normal legend order
  dat$flavor_plot = factor(dat$flavor, levels = rev(present_prefixes))
  dat$flavor_legend = factor(dat$flavor, levels = present_prefixes)
  
  # Split for plotting
  dat = split(dat, dat$gene)
  
  # Apply the function to each gene using mclapply
  lapply(dat, function(plot_data) {
    
    # Determine the plot title and axis labels
    if (valueType == "smltp") {
      y_label = bquote("Signed " ~ -log[10](italic(p)))
      plot_title = 
        bquote(
          "Signed " ~ -log[10](italic(p)) ~ "traces for" ~ .(plot_data$gene[1]))
    } else {
      y_label = valueType
      plot_title = paste(valueType, "traces for", plot_data$gene[1])
    }
    
    # Create the plot and save
    p = ggplot(plot_data, aes(x = as.numeric(marker), 
                              y = value, color = flavor_legend,
                              group = flavor_plot)) +
      geom_vline(xintercept = chrStarts, 
                 linetype = "dashed", color = "black") +
      geom_line() +
      scale_color_viridis_d(option = "viridis") +
      labs(title = plot_title,
           x = "Chromosome", y = y_label, color = NULL) + 
      theme_bw() +
      theme(legend.position = "bottom",
            legend.text = element_text(color = "black"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1,
                                       hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            plot.title = element_text(color = "black"),
            strip.text = element_text(color = "black")) +
      scale_x_continuous(breaks = chr_labels_pos, labels = chr_labels)
    ggsave(filename = paste0(rootDir, rootName, plot_data$gene[1], ".pdf"), 
           plot = p, width = plotWidth, height = plotHeight, units = "in")
  })
}

# Gets genetic correlation summary stats and histograms across genes 
mpSummarizer = function(data, directory) {
  
  # Ensure the directory exists
  ensure_directory(directory)
  
  # Create a list of effect types with corresponding transformations
  effects_list = list(
    mRNA_effects = lapply(data, function(x)
      x$mRNA.allelic_effect),
    Protein_effects = lapply(data, function(x)
      x$protein.allelic_effect),
    Difference_in_mRNA_and_protein_raw_effects = lapply(data, function(x)
      x$mRNA.allelic_effect - x$protein.allelic_effect),
    Difference_in_mRNA_and_protein_absolute_effects = lapply(data, function(x)
      abs(x$mRNA.allelic_effect - x$protein.allelic_effect)),
    Log10P_mRNA = lapply(data, function(x)
      -log10(x$mRNA.pval)),
    Log10P_Protein = lapply(data, function(x)
      -log10(x$protein.pval)),
    Log10P_mRNA_Protein = lapply(data, function(x)
      -log10(x$mRNA.pval / x$protein.pval))
  )
  
  # Loop through each effect type to create plots and summary statistics
  for (effect_type in names(effects_list)) {
    
    # Convert list to long format dataframe for each effect type
    combined_data = do.call(rbind, lapply(
      names(effects_list[[effect_type]]), function(name) {
        data.frame(
          effect = effects_list[[effect_type]][[name]],
          set = name
        )
      }))
    
    # Generate colors using the viridis palette
    types = unique(combined_data$set)
    colors = viridis::viridis(length(types))
    names(colors) = types
    
    # Determine plot title, file naming, and x-label based on effect type
    if (effect_type == "Log10P_mRNA_Protein") {
      plot_title = bquote(
        "-log"[10] ~ "(" ~ italic(p) ~ ") ratio of mRNA to protein")
      effect_name = "Log10_P_Value_Ratio"
      x_label = bquote("-log"[10] ~ "(" ~ italic(p) ~ ") ratio")
    } else if (effect_type == "Log10P_mRNA") {
      plot_title = bquote("-log"[10] ~ "mRNA" ~ italic(p) ~ "-values")
      effect_name = "Log10_mRNA_P_values"
      x_label = bquote("-log"[10] ~ "(" ~ italic(p) ~ ")")
    } else if (effect_type == "Log10P_Protein") {
      plot_title = bquote("-log"[10] ~ "protein" ~ italic(p) ~ "-values")
      effect_name = "Log10_Protein_P_Values"
      x_label = bquote("-log"[10] ~ "(" ~ italic(p) ~ ")")
    } else if (effect_type == "mRNA_effects"){
      plot_title = gsub("_", " ", effect_type)
      effect_name = effect_type
      x_label ='Allelic effect'
    } else if (effect_type == "protein_effects") {
      plot_title = gsub("_", " ", effect_type)
      effect_name = effect_type
      x_label ='Allelic effect'
    } else {
      plot_title = gsub("_", " ", effect_type)
      effect_name = effect_type
      x_label = "Difference in allelic effect"
    }
    
    # Correct file naming
    plot_filename = paste0(directory, effect_name, ".pdf")
    summary_filename = paste0(directory, effect_name, ".csv")

    # Plot for each effect type
    p = ggplot(combined_data, aes(x = effect, fill = set)) +
      geom_histogram(alpha = 0.6, position = 'identity', 
                     bins = 30, aes(y = after_stat(density)), 
                     color = "black") + # Add black outlines to the bars
      labs(x = x_label, y = "Density", title = plot_title) +
      scale_fill_manual(values = colors, name = NULL) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        legend.text = element_text(color = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"),
        plot.title = element_text(color = "black"),
        strip.text = element_text(color = "black")
      )
    
    # Save the plot with a specified width and height
    ggsave(plot_filename, plot = p, width = 7.5, height = 5)
    
    # Generate summary statistics for each effect type
    summary_stats = combined_data %>%
      group_by(set) %>%
      summarise(
        Min = min(effect),
        `1st Qu.` = quantile(effect, 0.25),
        Median = median(effect),
        Mean = mean(effect),
        `3rd Qu.` = quantile(effect, 0.75),
        Max = max(effect),
        .groups = 'drop'
      )
    
    # Save summary statistics
    write.table(summary_stats, file = summary_filename, sep = ",", 
                row.names = FALSE, quote = FALSE)
  }
}

# Point by point linked boxplot (the famous funfetti plot)
mpLinkPlotter = function(data, filename, infix, r = NULL, n = NULL) {
  data_long = data %>%
    dplyr::select(pairID, mRNA.allelic_effect, protein.allelic_effect) %>%
    tidyr::pivot_longer(cols = c(mRNA.allelic_effect, protein.allelic_effect), 
                        names_to = "effect_type", 
                        values_to = "effect_value")
  
  # Truncate effect_type labels after the period
  data_long$effect_type = sub("\\..*", "", data_long$effect_type)
  
  # Use the titler function to format the x-axis labels
  data_long$effect_type = sapply(data_long$effect_type,
                                 titler, r = r, n = n, u = FALSE)
  
  # Determine the title based on infix
  if (infix == "ixn") {
    title_infix = "of interaction QTLs"
  } else if (infix == "geno") {
    title_infix = "of genotype QTLs"
  } else {
    title_infix = infix
  }
  
  p = ggplot(data_long, aes(x = effect_type, 
                            y = effect_value, color = pairID)) +
    geom_boxplot(outlier.shape = 1, color = "black") +
    geom_line(aes(group = pairID), linewidth = 0.5) +
    scale_color_viridis(discrete = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 1, color = "black") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title = element_text(color = "black")
    ) +
    labs(x = "Effect type", y = "Allelic effect", 
         title = paste("mRNA and protein effects", title_infix))
  
  ggsave(filename, plot = p, width = 7.5, height = 5)
  
}

# Summarize mRNA and protein differences
mpSummarizer = function(data, directory) {
  
  # Ensure the directory exists
  ensure_directory(directory)
  
  # Create a list of effect types with corresponding transformations
  effects_list = list(
    mRNA_effect = lapply(data, function(x) 
      x$mRNA.allelic_effect),
    Protein_effect = lapply(data, function(x) 
      x$protein.allelic_effect),
    Genotype_effect = lapply(data, function(x) 
      x$mRNA.allelic_effect - x$protein.allelic_effect),
    Absolute_effect = lapply(data, function(x) 
      abs(x$mRNA.allelic_effect - x$protein.allelic_effect)),
    MinusLog10P_difference = lapply(data, function(x) 
      abs(x$mRNA.pval) - abs(x$protein.pval))
  )
  
  # Function to lowercase only the first letter of a string
  tolower_first = function(x) {
    paste0(tolower(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  }
  
  # Loop through each effect type to create plots and summary statistics
  for (effect_type in names(effects_list)) {
    
    # Convert list to long format dataframe for each effect type
    combined_data = do.call(rbind, lapply(
      names(effects_list[[effect_type]]), function(name) {
        data.frame(
          effect = effects_list[[effect_type]][[name]],
          set = name
        )
      }))
    
    # Generate colors using the viridis palette
    types = unique(combined_data$set)
    colors = viridis::viridis(length(types))
    names(colors) = types
    
    # Create plot title with special handling for MinusLog10P_difference
    if (effect_type == "MinusLog10P_difference") {
      plot_title = 
        bquote("mRNA and protein differences in" ~ -log[10](italic(p)))
    } else {
      plot_title =
        paste("mRNA and protein differences in", 
              tolower_first(gsub("_", "", effect_type)))
    }
    
    # Plot for each effect type
    p = ggplot(combined_data, aes(x = effect, fill = set)) +
      geom_histogram(alpha = 0.6, position = 'identity', 
                     bins = 30, aes(y = after_stat(density))) +
      labs(x = "Effect", y = "Density", title = plot_title) +
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
      )
    
    # Save the plot with a 3:2 width-to-height ratio
    plot_filename = paste0(directory, "/mRNA_and_protein_differences_in_",
                           tolower_first(gsub("_", "", effect_type)), ".pdf")
    ggsave(plot_filename, plot = p, width = 7.5, height = 5)
    
    # Generate summary statistics for each effect type
    summary_stats = combined_data %>%
      group_by(set) %>%
      summarise(
        Min = min(effect),
        `1st Qu.` = quantile(effect, 0.25),
        Median = median(effect),
        Mean = mean(effect),
        `3rd Qu.` = quantile(effect, 0.75),
        Max = max(effect),
        .groups = 'drop'
      )
    
    # Save summary statistics
    summary_filename = paste0(directory, "/mRNA_and_protein_differences_in_",
                              tolower_first(gsub("_", "", effect_type)), ".csv")
    write.table(summary_stats, file = summary_filename, sep = ",", 
                row.names = FALSE, quote = FALSE)
  }
}


# Trace LOD or effect size across studies individually
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
  chrStarts = c(chrStarts, length(unique(unlist(lapply(studies, function(study)
    colnames(study$gdata))))) + 1)
  
  # Calculate chromosome label positions
  chr_labels_pos = ((chrStarts[-length(chrStarts)] + chrStarts[-1] - 1)/2) %>%
    round()
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
    
    plotTitle = if (traceType == "s") "Signed LOD score traces for" else 
      "R-value traces for"
    yAxisLabel = if (traceType == "s") "Signed LOD score" else 
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
            axis.text.x = element_text(angle = 45, vjust = 1,
                                       hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            plot.title = element_text(color = "black"),
            strip.text = element_text(color = "black")) +
      scale_x_continuous(breaks = chr_labels_pos, labels = chr_labels)
    
    # Set the file prefix based on traceType
    file_prefix = if (traceType == "s") {
      "Signed_LOD_traces_for"
    } else {
      "R_value_traces_for"
    }
    
    # Save the plot with the appropriate file name
    ggsave(filename = paste0(fn, "/", file_prefix, "_", 
                             rownames(suppMat)[i], ".pdf"), 
           plot = p, width = plotWidth, height = plotHeight, units = "in")
  }
  
  pbapply::closepb(pb)
}

print('SOURCE OF CONSTANTS AND FUNCTIONS WAS SUCCESSFUL')