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
personal_lib_path = "/data/gusev/ADQTL/envs/adqtl/lib/R/library"

# Create the directory if it doesn't exist
if (!dir.exists(personal_lib_path)) {
  dir.create(personal_lib_path, recursive = TRUE)
}

# Add the personal library path to the library paths
.libPaths(personal_lib_path)

# Set the default CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install and use BiocManager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = personal_lib_path, quietly = TRUE)
}

# # ANCHOR
# # All package names
# pkgs = c("broom.mixed", "GEOquery", "GenomicRanges", "gridExtra", "leaps", 
#          "lme4", "lmerTest", "parallel", "patchwork", "pbapply", 
#          "plyr", "progress", "readxl", "rrBLUP", "tidyverse", "viridis")

# # Install using correct package installation tools
# for (pkg in pkgs) {
#   if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
#     message("Installing ", pkg, "...")
#     if (pkg %in% c("GEOquery", "GenomicRanges")) {
#       suppressMessages(BiocManager::install(
#         pkg, lib = personal_lib_path, quietly = TRUE))
#     } else {
#       suppressMessages(install.packages(
#         pkg, lib = personal_lib_path, quietly = TRUE))
#     }
#   } else {
#     message(pkg, " is already installed.")
#   }
# }

# All package names
pkgs = c("broom.mixed", "GenomicRanges", "gridExtra", "leaps", "lme4", "lmerTest",
         "mashr", "parallel", "patchwork", "pbapply", "pbmcapply", "plyr", 
         "progress", "RColorBrewer", "readxl", "rrBLUP", "tidyverse", "viridis")

# Install using correct package installation tools
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message("Installing ", pkg, "...")
    if (pkg %in% c("GenomicRanges")) {
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
wd = "/PHShome/kv1026/p05/"
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

# subset of RNA-only data used to understand genotype-medium interactions
medMeta = classicMetadata[classicMetadata$medium %in% c('SC','YNB2G') &
                            classicMetadata$study != 'Albert_2018' &
                            classicMetadata$analyte == 'mRNA',]

# Study root names
studyList = c("Brem_2002", "Brem_2005", "SmithEth_2008", "SmithGlc_2008",
              "Khan_2009", "Foss_2011", "Picotti_2013", "Albert_2018", 
              "AlbertAlt_2018", "GrossbachPRT_2022", "GrossbachRNA_2022")

studyTypes = c("mRNA","mRNA","special","special",
               "protein","protein","protein",
               "mRNA","mRNA","protein","mRNA")

studyNiceNames = c(
  "Brem_2002" = "Brem 2002",
  "Brem_2005" = "Brem 2005",
  "SmithEth_2008" = "Smith 2008 (EtOH)",
  "SmithGlc_2008" = "Smith 2008 (Glucose)",
  "Khan_2009" = "Khan 2009",
  "Foss_2011" = "Foss 2011",
  "Picotti_2013" = "Picotti 2013",
  "Albert_2018" = "Albert 2018",
  "AlbertAlt_2018" = "Albert 2018 (Alt)",
  "Albert_2014" = "Albert 2014",
  "BrionRNA_2020" = "Brion 2020 (RNA)",
  "BrionPRT_2020" = "Brion 2020 (Protein)",
  "GrossbachRNA_2022" = "Grossbach 2022 (RNA)",
  "GrossbachPRT_2022" = "Grossbach 2022 (Protein)"
)

# Threshold for NA removal
threshNA = 0.25

# Store original thread settings
old_threads = Sys.getenv(
  c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS","BLIS_NUM_THREADS"),
  unset = NA
)

# Anlyte two color scheme
colAnalyte = c(mRNA = "#191970", protein = "#EC5800")

#########################
# GENERAL FUNCTION BANK #
#########################

# create an integer seed from one or more strings
seed_from = function(...) {
  s = enc2utf8(paste(..., collapse = "|"))
  if (requireNamespace("digest", quietly = TRUE)) {
    hx = digest::digest(s, algo = "xxhash32", serialize = FALSE)   # 8 hex chars
    out = suppressWarnings(strtoi(hx, base = 16L))
    out = as.integer(out %% .Machine$integer.max)
  } else {
    # base-R fallback: rolling hash (no overflow), mod 2^31-1
    u = utf8ToInt(s)
    h = 0
    for (i in seq_along(u)) {
      h = (h * 131 + u[i]) %% 2147483647  # 2^31-1
    }
    out = as.integer(h)
  }
  if (is.na(out) || out == 0L) out = 1L   # set.seed() needs 1..2^31-1
  out
}

# Respect explicit override via env or option; fall back to SLURM; then detectCores()
get_cores = function(max_tasks = NULL) {
  # explicit override first
  opt = getOption("mapqtl.cores", NULL)
  env = suppressWarnings(as.integer(Sys.getenv("MAPQTL_CORES", NA)))
  if (!is.null(opt))           cores = as.integer(opt)
  else if (!is.na(env))        cores = env
  else {
    slurm = suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
    job   = suppressWarnings(as.integer(sub("\\D.*$", "", Sys.getenv("SLURM_JOB_CPUS_PER_NODE", ""))))
    det   = parallel::detectCores()
    cores = if (!is.na(slurm)) slurm else if (!is.na(job)) job else det
  }
  if (is.na(cores) || cores < 1L) cores = 1L
  if (!is.null(max_tasks)) cores = min(cores, as.integer(max_tasks))
  cores
}

# Functions to set single-threaded operation
set_single_thread = function() {
  Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
             MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
             BLIS_NUM_THREADS="1")
}

# Function to restore original thread settings
restore_threads = function() {
  nm = names(old_threads)
  for (j in seq_along(nm)) {
    if (is.na(old_threads[j])) Sys.unsetenv(nm[j])
    else Sys.setenv(structure(old_threads[j], names = nm[j]))
  }
}

# Negation of in
`%nin%` = Negate(`%in%`)

# Szudzik's pairing function
# used to generate unique seeds for multithreading
szudzik_pair = function(a, b) {
  if (a >= b) {
    return(a * a + a + b)
  } else {
    return(b * b + a)
  }
}

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
  
  # ANCHOR
  # datGDS = GEOquery::getGEO(filename = paste0(dir, filename)) %>%
  #   GEOquery::Table()
  
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

# Function for calculating intersections between two GRanges objects
# add optional padding
# adds padding to the first GRanges object (gr1) to expand its ranges
# This expansion reduces stringency 
# Resizes gr1 with the specified padding
# finds overlaps of rezied gr1 with gr2
# returns the intersected regions
padded_pintersect = function(gr1, gr2, pad = 0) {
  # Ensure gr1 and gr2 are of the same length for pairwise comparison
  if (length(gr1) == 0 || length(gr2) == 0 || length(gr1) != length(gr2)) {
    # Return an empty logical vector if inputs are empty or lengths do not match
    return(logical(length(gr1)))
  }
  
  # Adjust gr1 to include padding
  gr1_padded = resize(gr1, width(gr1) + 2 * pad, fix = "center")
  
  # Ensure no negative start positions
  start(gr1_padded) = pmax(start(gr1_padded), 1)
  
  # Perform pairwise intersection and return a logical vector
  overlaps = pintersect(gr1_padded, gr2)
  hits = width(overlaps) > 0  # Determine if each intersection is non-empty
  
  return(hits)
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

# faster p-value calculation
# Adds non-numeric p-handling
fasterP = function(n.pheno, pheno.s, gdata.s) {
  r  = crossprod(pheno.s, gdata.s) / (n.pheno - 1)
  df = n.pheno - 2
  t  = r * sqrt(df / (1 - r^2))
  p  = 2 * pt(-abs(t), df = df)
  p[!is.finite(p)] = 1
  return(p)
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
mapQTL = function(covariates, background.QTL, phenotype.raw, pheno.scaled, 
                  gdata, gdata.scaled, 
                  chromosomes = NULL,
                  n.perm = 1000, 
                  FDR.thresh = 0.05, 
                  Arandom = TRUE,
                  chunk = 100L,
                  seed = NULL) {

  # time the entire function
  t_all = proc.time()[["elapsed"]]
  on.exit({
    cat(sprintf("mapQTL: total elapsed = %.2f s\n",
                proc.time()[["elapsed"]] - t_all))
  }, add = TRUE)
                  
  # set the seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
    cat(sprintf("mapQTL: seed set to %d\n", seed))
  }

  # figure out how many cores you have
  cores_total = get_cores()

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
    cat("Mapping:", cc, "\n")
    
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
    if (Arandom) {

      Aloco = rrBLUP::A.mat(do.call('cbind', gdata.by.chr[-match(cc, names(gdata.by.chr))]))/2
      eigA  = doEigenA_forMM(presid, Aloco)

      # SVD/eigen basis for K
      svdA = svd(Aloco)        # symmetric => eigenvectors in svd$u
      U    = svdA$u
      lam  = svdA$d

      # Helper: (r1*K + r2*I)^{-1} y without forming W
      solve_rr = function(y, r1, r2) {
        Uy = crossprod(U, y)                 # t(U) %*% y
        Uy = Uy / (r1 * lam + r2)            # elementwise
        U %*% Uy                              # back to original space
      }

      # Parallel over columns (cap by work size), keep workers single-threaded
      ncoresMM = min(cores_total, ncol(presid))
      mm_chunk = getOption("mapqtl.mm_chunk", 8L * ncoresMM)
      mm_chunk = max(1L, min(mm_chunk, ncol(presid)))

      set_single_thread()
      on.exit(restore_threads(), add = TRUE)

      p      = ncol(presid)
      nbatch = ceiling(p / mm_chunk)
      res_list = vector("list", p)

      t0 = proc.time()[["elapsed"]]
      for (b in seq_len(nbatch)) {
        lo = (b - 1L) * mm_chunk + 1L
        hi = min(b * mm_chunk, p)
        cat(sprintf("%s: MM residualization batch %d/%d (cols %d–%d) using %d cores\n",
                    cc, b, nbatch, lo, hi, ncoresMM))

        tb = proc.time()[["elapsed"]]
        res_b = parallel::mclapply(
          lo:hi,
          function(tp) {
            y  = presid[, tp]
            rr = m.S(y, K = Aloco, theta = eigA$theta, Q = eigA$Q)  # c(r1, r2)
            if (rr[1] <= 0) return(y)
            inv_y = solve_rr(y, r1 = rr[1], r2 = rr[2])
            Uiny  = crossprod(U, inv_y)
            Kiny  = U %*% (lam * Uiny)
            as.vector(y - rr[1] * Kiny)
          },
          mc.cores = ncoresMM,
          mc.set.seed = TRUE
        )
        res_list[lo:hi] = res_b
        rm(res_b); gc()

        # batch timing + ETA
        elapsed_b = proc.time()[["elapsed"]] - tb
        done = hi
        elapsed_total = proc.time()[["elapsed"]] - t0
        rate = done / max(1e-9, elapsed_total)           # cols/sec so far
        remaining = p - done
        eta = if (rate > 0) remaining / rate else NA_real_

        cat(sprintf("%s: batch %d done in %.2f s; %.1f cols/sec; ETA %.1f s\n",
                    cc, b, elapsed_b, rate, eta))
      }

      presid = do.call(cbind, res_list)
      colnames(presid) = colnames(preal)
      rownames(presid) = rownames(preal)
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
      
      # Announce permutations and time them
      chunk = max(10L, as.integer(chunk))  # enforce lower bound first
      cat(sprintf("%s: running %d permutations (FDR target = %.3f, chunk = %d)\n",
            cc, n.perm, FDR.thresh, chunk))

      t_perm = proc.time()

      # --- SCOPE BLAS/OMP TO 1 THREAD ---
      set_single_thread()
      perm_elapsed = NA_real_

      tryCatch({

        # announce the number of cores we are using
        cat(sprintf("%s: using %d cores\n", cc, cores_total))

        permLOD = matrix(0, nrow = nrow(obsLOD$LOD), ncol = n.perm)

        k = 1L
        while (k <= n.perm) {
          m = min(chunk, n.perm - k + 1L)

          repeat {
            # try parallel with current cores_total
            res = try(
              parallel::mclapply(
                X = seq_len(m),
                FUN = function(j) {
                  ridx = sample.int(nrow(presid))
                  apply(
                    fasterLOD(nrow(presid), presid[ridx, , drop = FALSE], gdata.s.by.chr[[cc]]),
                    1, max
                  )
                },
                mc.cores    = cores_total,
                mc.set.seed = TRUE
              ),
              silent = TRUE
            )

            if (!inherits(res, "try-error")) break  # success

            # adaptive backoff: first reduce cores, then chunk
            if (cores_total > 1L) {
              cores_total = max(1L, cores_total %/% 2L)
              message(sprintf("  batch failed; reducing cores to %d and retrying", cores_total))
              next
            }

            if (m > 10L) {
              # shrink chunk and recompute m for this batch
              chunk = max(10L, m %/% 2L)
              m = min(chunk, n.perm - k + 1L)
              message(sprintf("  still failing at 1 core; reducing chunk to %d and retrying", chunk))
              next
            }

            stop("Permutation chunk failed even at 1 core and minimal chunk size.")
          }

          permLOD[, k:(k + m - 1L)] = do.call(cbind, res)
          cat(sprintf("  %d / %d permutations done (cores=%d, chunk=%d)\n", k + m - 1L, n.perm, cores_total, chunk))
          k = k + m
          rm(res); gc()
        }

        perm_elapsed = (proc.time() - t_perm)[["elapsed"]]
        cat(sprintf("permutations finished in %.2f s (%.1f perms/sec)\n",
            perm_elapsed, n.perm / perm_elapsed))

      }, finally = {
        # --- ALWAYS RESTORE THREADING RIGHT AFTER PERMS ---
        restore_threads()
      })
      
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
findBackgroundQTLs = function(studies, fn) {

  bgQTLs = lapply(names(studies), function(study_name) {
    study = studies[[study_name]]
    rda_path = paste0(fn, study_name, "_Background.rda")
    obj_name = paste0(study_name, "_Background")

    if (!file.exists(rda_path)) {
      cat(sprintf("MAPPING BACKGROUND QTLs for %s\n", study_name))
      t0 = proc.time()[["elapsed"]]

      bgQTL = find.background.QTL(study$covars,
                                  study$pheno,
                                  scale(study$pheno),
                                  study$gdata,
                                  scale(study$gdata))

      # Save
      assign(obj_name, bgQTL)
      save(list = obj_name, file = rda_path)

      cat(sprintf("Finished mapping background QTLs for %s in %.2f s\n",
                  study_name, proc.time()[["elapsed"]] - t0))

    } else {
      cat(sprintf("LOADING BACKGROUND QTLs for %s\n", study_name))
      load(rda_path)
      bgQTL = get(obj_name)
    }

    bgQTL
  })

  names(bgQTLs) = names(studies)
  bgQTLs
}

# Wrapper for mapping QTLs using Josh's pipeline over a list
fullMap = function(studies, bgQTLs, fn, method = 'f',
                   n.perm = 1000, FDR.thresh = 0.05, chunk = 100L) {

  # Mapping type
  if (method == 'f') {
    typeRoot = '_fullPeaks';  popCrxn = TRUE
  } else if (method == 's') {
    typeRoot = '_simplePeaks'; popCrxn = FALSE
  } else {
    warning('Invalid method selected. Choose f or s. Returning NULL')
    return(NULL)
  }

  # Sequential over studies; mapQTL handles internal parallelism
  qtlResults = lapply(names(studies), function(study_name) {
    study = studies[[study_name]]
    background = if (method == 'f') bgQTLs[[study_name]] else list()
    rda_path = paste0(fn, study_name, typeRoot, '.rda')
    obj_name = paste0(study_name, typeRoot)

    if (file.exists(rda_path)) {
      cat(sprintf("Files already exist for %s (%s). Loading...\n", study_name, method))
      load(rda_path)
      return(get(obj_name))
    }

    cat(sprintf("Mapping QTLs for %s (%s)...\n", study_name, method))
    t0 = proc.time()[["elapsed"]]

    qtlResult = mapQTL(study$covars,
                       background,
                       study$pheno, scale(study$pheno),
                       study$gdata, scale(study$gdata),
                       n.perm = n.perm,
                       FDR.thresh = FDR.thresh,
                       Arandom = popCrxn,
                       chunk = chunk,
                       seed = seed_from(study_name))   # <-- forward chunk

    cat(sprintf("Finished mapping QTLs for %s in %.2f s\n", 
                study_name, proc.time()[["elapsed"]] - t0))


    assign(obj_name, qtlResult)
    save(list = obj_name, file = rda_path)
    qtlResult
  })

  names(qtlResults) = names(studies)
  qtlResults
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

# does a stratified permutation on a vector using the groups (used in JMT)
stratPerm = function(s, v) {

  if (length(s) != length(v)) {
    stop("Strata and vector must have the same length.")
  }
  
  unique_strata = unique(s)
  
  permuted_v = v
  
  for (stratum in unique_strata) {
    indices = which(s == stratum)
    permuted_v[indices] = sample(v[indices], length(indices), replace = FALSE)
  }
  
  return(permuted_v)
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
jointMarkerTest = function(studies, targetMarkers, unionBg, smetadata, 
                           suppGenes, nPerms, FDR, fixedEffect, fn, 
                           correctBG = FALSE, correctGRM = FALSE,
                           cores = 1, chromosomes = NULL) {
  
  print('PREPARING MAPPING OBJECTS AND DEALING WITH GENE SUPPORT')
  
  # Create the directory if it does not already exist
  outPath = paste0(imed, fn, '/')
  if (!dir.exists(outPath)) {
    dir.create(outPath, showWarnings = TRUE, recursive = FALSE)
  }
  
  # Get true marker names
  trueMarkers = colnames(studies[[1]]$gdata)
  
  # Define default 16 S. cervisiae chromosomes if none are provided
  if (is.null(chromosomes)) {
    chromosomes = paste0("chr", as.roman(1:16))
  }
  
  # Get a list of which studies support every gene
  suppList = lapply(suppGenes, function(gene) {
    studies_containing_gene = sapply(smetadata$genes, function(gene_list) {
      gene %in% gene_list
    })
    smetadata$study[studies_containing_gene]
  }) %>% `names<-`(suppGenes)
  
  # get unique study combinations
  uqCombos = unique(suppList)
  uqKey = match(suppList, uqCombos)
  
  # Create per-gene metadata object (gene-long list with segregants as rows)
  longMeta = mapply(function(x, y) {
    rep(y, nrow(x$gdata))
  }, studies, names(studies)) %>%
    unlist() %>%
    unname() %>%
    as.data.frame() %>%
    `names<-`(c('study')) %>%
    left_join(smetadata %>% select(study, all_of(fixedEffect)),
              by = "study") %>% mutate(geno = NA)
  
  # Create gene-specific metadata using longMeta for subsetting
  longMeta = lapply(uqCombos, function(combo){
    longMeta %>% filter(study %in% combo)})
  
  # Create a joint genotype object across all studies
  jointGeno = lapply(uqCombos, function(combo) {
    names(studies) %>%
      .[. %in% combo] %>%
      lapply(function(study_name) studies[[study_name]]$gdata) %>%
      do.call(rbind, .) %>%
      return()
  })
  
  # Chrom indexing key for all markers (e.g. "ChrI")
  cc = do.call('rbind', strsplit(colnames(jointGeno[[1]]), ':'))[,1]
  
  # Create phenotypes object (gene long list of Seg x Chrom matrices)
  jointPhenoScaled = lapply(suppGenes, function(gene) {
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
  
  if (correctBG) {
    print('CORRECTING FOR STRONG GENETIC EFFECTS ON OTHER CHROMOSOMES')
    pb = txtProgressBar(min = 0, max = length(uqKey), style = 3)
    jointPhenoScaled = mapply(function(key, pheno, index) {
      result = adjust_for_background_QTLs(
        jointGeno[[key]], pheno, unionBg, chromosomes)
      setTxtProgressBar(pb, index)
      return(result)
    }, uqKey, jointPhenoScaled,
    seq_along(uqKey), SIMPLIFY = FALSE)
    cat("\n")
  }
  
  if (correctGRM) {
    print('CORRECTING FOR RELATEDNESS ON OTHER CHROMOSOMES')
    pb = txtProgressBar(min = 0, max = length(uqKey), style = 3)
    jointPhenoScaled = mapply(function(key, pheno, index) {
      result = calculate_BLUPs(jointGeno[[key]], pheno)
      setTxtProgressBar(pb, index)
      return(result)
    }, uqKey, jointPhenoScaled,
    seq_along(uqKey), SIMPLIFY = FALSE)
    cat("\n")
  }
  
  # Scale genotypes to control for variance at each locus gene-specifically
  jointGeno = lapply(jointGeno, scale)
  
  print('PREPARING PERMUTATION MATICES AND FINALIZAING MAPPING PREPARATIONS')
  
  # Chrom indexing key for targeted markers (e.g. "ChrI")
  ct = do.call('rbind', strsplit(targetMarkers, ':'))[,1]
  
  # Important model terms for later listing
  mterms = c('geno', paste0(fixedEffect, ':geno'))
  
  # Test associations at each marker
  print('BEGINNING JOINT QTL MAPPING')
  lapply(targetMarkers, function(marker) {
    
    # Get index of current marker and true marker index
    mI = which(targetMarkers == marker)
    mT = which(trueMarkers == marker)
    print(paste0('mapping associations of marker #',
                 mI, ' of ', length(targetMarkers), ' at position ', marker))
    
    # Real associations (terms many matrices)
    rassocs = rep(list(
      list(
        gene = suppGenes,
        marker = as.character(rep(marker, length(suppGenes))),
        flavor = as.character(rep(NA, length(suppGenes))),
        effect = as.numeric(rep(NA, length(suppGenes))),
        p_value = as.numeric(rep(NA, length(suppGenes)))
      ) %>% as.data.frame()
    ), length(mterms)) %>% `names<-`(mterms) %>%
      imap(function(df, name) {
        df %>% mutate(flavor = name)
      })
    
    # Permutation assocations (mterms long Genes x Perms mats)
    passocs = lapply(mterms, function(mterm) {
      matrix(NA, nrow = length(suppGenes), ncol = nPerms)
    }) %>% `names<-`(mterms)
    
    # create the correct metadata objects
    perMeta = mapply(function(g, m) {
      
      df = data.frame(
        study = m$study,
        geno = g[,marker],
        pheno = rep(NA, length(m$study))
      )
      df[[fixedEffect]] = m[[fixedEffect]]
      return(df)
    }, jointGeno, longMeta, SIMPLIFY = FALSE)
    
    # Sweep over genes
    for (gene in seq_along(suppGenes)) {
      
      perMeta[[uqKey[gene]]]$pheno = 
        jointPhenoScaled[[gene]][,ct[mI]]
      
      # Run the mixed model on the real data and get n and fixed effects
      formula = as.formula(paste("pheno ~ (1|study) +", 
                                 fixedEffect, "+ geno + geno:", fixedEffect))
      datModel = lmer(
        formula, data = perMeta[[uqKey[gene]]],
        control = lmerControl(check.conv.singular = "ignore"))
      nGI = nrow(perMeta[[uqKey[gene]]])
      fixies = fixef(datModel) 
      names(fixies)[3:4] = mterms
      datANOVA = anova(datModel)
      
      # Assign term values
      for (mterm in mterms) {
        rassocs[[mterm]][gene, 'effect'] = fixies[mterm]
        rassocs[[mterm]][gene, 'samplesize'] = nGI
        rassocs[[mterm]][gene, 'p_value'] = datANOVA[mterm, "Pr(>F)"]
        rassocs[[mterm]][gene, 'mlog10p'] = 
          -log10(rassocs[[mterm]][gene, 'p_value'])
      }
      
      
      permObj = mclapply(1:nPerms,function(perm){
        
        # set the seed using the szudzik pair for uniqueness
        set.seed(szudzik_pair(gene, perm))
        
        # Break genotype association using the permutation key
        perMeta[[uqKey[gene]]]$geno = stratPerm(
          s = perMeta[[uqKey[gene]]]$study,
          v = perMeta[[uqKey[gene]]]$geno)
        
        # Run the mixed model on the permutation data and get fixed effects
        datModel = lmer(
          formula, data = perMeta[[uqKey[gene]]],
          control = lmerControl(check.conv.singular = "ignore"))
        fixies = fixef(datModel) 
        
        names(fixies)[3:4] = mterms
        
        # Do ANOVA for p-values
        datANOVA = anova(datModel)
        
        # Assign permutation geno and 'analyte:geno' LOD values for FDR step
        return(sapply(mterms,function(mterm)-log10(datANOVA[mterm, "Pr(>F)"])))
        
      },mc.cores=cores) %>% do.call(what=rbind)
      
      
      for (p in seq_along(mterms)) {
        passocs[[mterms[p]]][gene,] = permObj[, p]
      }
      
    }
    
    # Use empirical FDRs to determine QTLs
    lapply(mterms, function(mterm) {
      
      # eFDRs are mean perm FDs over real discoveries at data-driven threshes
      lt = sort(rassocs[[mterm]]$mlog10p)
      rPass = sapply(lt, function(t) {sum(rassocs[[mterm]] >= t)})
      pPass = sapply(lt, function(t)
        mean(apply(passocs[[mterm]], 2, function(x) sum(x >= t))))
      eFDR = pPass / rPass
      
      # Remove NaNess from div by zero etc. and make FDR curve monotone
      eFDR[is.na(eFDR)] = 0
      eFDR[!is.finite(eFDR)] = 0
      eFDR = rev(cummax(rev(eFDR)))
      
      # Use minimum threshold passing eFDR cutoff to write QTLs (if any)
      cuts = eFDR < FDR
      if(any(cuts)){
        rassocs[[mterm]][rassocs[[mterm]]$mlog10p >= lt[min(which(cuts))],] %>% 
          write.csv(file = paste0(
            outPath, 'QTLs_for_marker_', mT, "_of_", length(trueMarkers),
            "_at_FDR_", FDR, "_and_type_", gsub(":", "-", mterm),
            ".csv"), row.names = FALSE)}
      
    })
    
    return(paste('Finished', marker))
    
  })
}

# jointMarker FDR
jointMarkerFDR = function(studies, targetMarkers, unionBg, smetadata,
                          suppGenes, FDR, fixedEffect = NULL, fn,
                          dropLOD = 1.5,
                          correctBG = FALSE, correctGRM = FALSE,
                          cores = 1, chromosomes = NULL) {

  message("PREPARING MAPPING OBJECTS AND DEALING WITH GENE SUPPORT")

  outPath = paste0(imed, fn, "/")
  if (!dir.exists(outPath)) dir.create(outPath, showWarnings = TRUE, recursive = FALSE)

  trueMarkers = colnames(studies[[1]]$gdata)
  if (is.null(chromosomes)) chromosomes = paste0("chr", as.roman(1:16))

  # which studies support each gene
  suppList = lapply(suppGenes, function(gene) {
    studies_containing_gene = sapply(smetadata$genes, function(gene_list) gene %in% gene_list)
    smetadata$study[studies_containing_gene]
  }) %>% `names<-`(suppGenes)

  uqCombos = unique(suppList)
  uqKey = match(suppList, uqCombos)

  # base long metadata
  baseLongMeta =
    mapply(function(x, y) rep(y, nrow(x$gdata)), studies, names(studies)) %>%
    unlist() %>% unname() %>% as.data.frame() %>%
    `names<-`("study")

  if (!is.null(fixedEffect)) {
    baseLongMeta = baseLongMeta %>%
      left_join(smetadata %>% dplyr::select(study, all_of(fixedEffect)), by = "study")
  }

  baseLongMeta = baseLongMeta %>% mutate(geno = NA)

  longMeta = lapply(uqCombos, function(combo) baseLongMeta %>% dplyr::filter(study %in% combo))

  jointGeno = lapply(uqCombos, function(combo) {
    names(studies) %>%
      .[. %in% combo] %>%
      lapply(function(study_name) studies[[study_name]]$gdata) %>%
      do.call(rbind, .)
  })

  # chrom key for markers
  ct = do.call("rbind", strsplit(targetMarkers, ":"))[, 1]

  # phenotype matrices, replicated across chromosomes
  jointPhenoScaled = lapply(suppGenes, function(gene) {
    pheno_values = unlist(lapply(suppList[[gene]], function(studyName) {
      if (studyName %in% names(studies)) {
        st = studies[[studyName]]
        if (gene %in% colnames(st$pclean)) return(st$pclean[, gene])
      }
      NULL
    }))

    if (!is.null(pheno_values)) {
      pheno_matrix = matrix(rep(pheno_values, length(chromosomes)),
                            ncol = length(chromosomes), byrow = FALSE)
      colnames(pheno_matrix) = chromosomes
      rownames(pheno_matrix) = names(pheno_values)
      pheno_matrix
    } else {
      NULL
    }
  }) %>% `names<-`(names(suppList))

  if (correctBG) {
    message("CORRECTING FOR STRONG GENETIC EFFECTS ON OTHER CHROMOSOMES")
    pb = txtProgressBar(min = 0, max = length(uqKey), style = 3)
    jointPhenoScaled = mapply(function(key, pheno, index) {
      out = adjust_for_background_QTLs(jointGeno[[key]], pheno, unionBg, chromosomes)
      setTxtProgressBar(pb, index)
      out
    }, uqKey, jointPhenoScaled, seq_along(uqKey), SIMPLIFY = FALSE)
    close(pb); cat("\n")
  }

  if (correctGRM) {
    message("CORRECTING FOR RELATEDNESS ON OTHER CHROMOSOMES")
    pb = txtProgressBar(min = 0, max = length(uqKey), style = 3)
    jointPhenoScaled = mapply(function(key, pheno, index) {
      out = calculate_BLUPs(jointGeno[[key]], pheno)
      setTxtProgressBar(pb, index)
      out
    }, uqKey, jointPhenoScaled, seq_along(uqKey), SIMPLIFY = FALSE)
    close(pb); cat("\n")
  }

  jointGeno = lapply(jointGeno, scale)

  mapMode = if (is.null(fixedEffect)) "genotype" else "interaction"

  # final assoc container (single flavor per run)
  fassocs = data.frame(
    gene       = character(),
    marker     = character(),  # peak marker
    flavor     = character(),  # "genotype" or "interaction"
    effect     = numeric(),
    p_value    = numeric(),
    samplesize = numeric(),
    mlog10p    = numeric(),
    LOD        = numeric(),
    LEFT       = character(),
    RIGHT      = character(),
    BH         = numeric(),
    pass       = integer(),
    stringsAsFactors = FALSE
  )

  message("BEGINNING JOINT QTL MAPPING WITH BH-FDR AND MULTIPLE PASSES")

  # helper: compute wings for one gene using LOD curve of THIS PASS
  compute_wings = function(gene, peak_marker, peak_LOD, pass_df, ct_vec) {
    idx_peak = match(peak_marker, targetMarkers)
    chr_peak = ct_vec[idx_peak]

    chr_idx = which(ct_vec == chr_peak)
    chr_idx = chr_idx[order(chr_idx)]
    pos = match(idx_peak, chr_idx)

    thr = peak_LOD - dropLOD

    # left
    left_pos = pos
    while (left_pos > 1) {
      cand_idx = chr_idx[left_pos - 1]
      cand_marker = targetMarkers[cand_idx]
      cand_LOD = pass_df$LOD[pass_df$gene == gene & pass_df$marker == cand_marker][1]
      if (is.na(cand_LOD) || cand_LOD < thr) break
      left_pos = left_pos - 1
    }
    LEFT = targetMarkers[chr_idx[left_pos]]

    # right
    right_pos = pos
    while (right_pos < length(chr_idx)) {
      cand_idx = chr_idx[right_pos + 1]
      cand_marker = targetMarkers[cand_idx]
      cand_LOD = pass_df$LOD[pass_df$gene == gene & pass_df$marker == cand_marker][1]
      if (is.na(cand_LOD) || cand_LOD < thr) break
      right_pos = right_pos + 1
    }
    RIGHT = targetMarkers[chr_idx[right_pos]]

    list(LEFT = LEFT, RIGHT = RIGHT)
  }

  nGenes   = length(suppGenes)
  nMarkers = length(targetMarkers)

  pass = 1

  repeat {

    message(paste0("--- PASS ", pass, " (", mapMode, ") ---"))

    # -----------------------------
    # PREALLOCATE PER-PASS RESULTS
    # -----------------------------
    nRows = nGenes * nMarkers

    gene_vec   = rep(suppGenes, times = nMarkers)
    marker_vec = rep(targetMarkers, each  = nGenes)
    pass_vec   = rep(pass, nRows)
    flavor_vec = rep(mapMode, nRows)

    effect_vec     = rep(NA_real_, nRows)
    pval_vec       = rep(NA_real_, nRows)
    n_vec          = rep(NA_real_, nRows)
    mlog10p_vec    = rep(NA_real_, nRows)
    lod_vec        = rep(NA_real_, nRows)

    pb_pass = txtProgressBar(min = 0, max = nMarkers, style = 3)

    for (mI in seq_along(targetMarkers)) {

      marker = targetMarkers[mI]

      # per marker, per unique combo metadata (geno changes with marker; study/fixedEffect constant)
      perMeta = mapply(function(g, m) {
        df = data.frame(
          study = m$study,
          geno  = g[, marker],
          pheno = rep(NA, length(m$study))
        )
        if (!is.null(fixedEffect)) df[[fixedEffect]] = m[[fixedEffect]]
        df
      }, jointGeno, longMeta, SIMPLIFY = FALSE)

      # build formulas once per marker
      if (is.null(fixedEffect)) {
        f0 = as.formula("pheno ~ (1|study)")
        f1 = as.formula("pheno ~ (1|study) + geno")
      } else {
        f0 = as.formula(paste0("pheno ~ (1|study) + ", fixedEffect, " + geno"))
        f1 = as.formula(paste0("pheno ~ (1|study) + ", fixedEffect, " + geno + geno:", fixedEffect))
      }

      ctrl = lmerControl(check.conv.singular = "ignore")

      run_one_gene = function(g) {
        perMeta_g = perMeta[[uqKey[g]]]
        perMeta_g$pheno = jointPhenoScaled[[g]][, ct[mI]]

        m0 = lmer(f0, data = perMeta_g, control = ctrl, REML = FALSE)
        m1 = lmer(f1, data = perMeta_g, control = ctrl, REML = FALSE)

        ll0 = as.numeric(logLik(m0))
        ll1 = as.numeric(logLik(m1))
        lod = (ll1 - ll0) / log(10)

        # LRT p-value for m1 vs m0
        pval = anova(m0, m1)[2, "Pr(>Chisq)"]

        fx = fixef(m1)
        eff =
          if (is.null(fixedEffect)) {
            unname(fx["geno"])
          } else {
            ixn_terms = grep("geno:", names(fx), value = TRUE)
            if (length(ixn_terms) == 0) ixn_terms = grep(":geno", names(fx), value = TRUE)
            ixn_terms = unique(ixn_terms)

            if (length(ixn_terms) == 1) {
              unname(fx[ixn_terms])
            } else if (length(ixn_terms) > 1) {
              max(abs(unname(fx[ixn_terms])))
            } else {
              NA_real_
            }
          }

        list(
          effect     = eff,
          p_value    = pval,
          samplesize = nrow(perMeta_g),
          mlog10p    = -log10(pval),
          LOD        = lod
        )
      }

      gene_out =
        if (cores > 1) parallel::mclapply(seq_len(nGenes), run_one_gene, mc.cores = cores)
        else lapply(seq_len(nGenes), run_one_gene)

      # fill preallocated vectors for this marker block
      block_start = (mI - 1) * nGenes + 1
      block_end   = mI * nGenes
      idx_block   = block_start:block_end

      effect_vec[idx_block]  = vapply(gene_out, `[[`, numeric(1), "effect")
      pval_vec[idx_block]    = vapply(gene_out, `[[`, numeric(1), "p_value")
      n_vec[idx_block]       = vapply(gene_out, `[[`, numeric(1), "samplesize")
      mlog10p_vec[idx_block] = vapply(gene_out, `[[`, numeric(1), "mlog10p")
      lod_vec[idx_block]     = vapply(gene_out, `[[`, numeric(1), "LOD")

      setTxtProgressBar(pb_pass, mI)
    }

    close(pb_pass); cat("\n")

    # build per-pass df once
    passRassocs = data.frame(
      gene       = gene_vec,
      marker     = marker_vec,
      flavor     = flavor_vec,
      effect     = effect_vec,
      p_value    = pval_vec,
      samplesize = n_vec,
      mlog10p    = mlog10p_vec,
      LOD        = lod_vec,
      BH         = NA_real_,
      pass       = pass_vec,
      stringsAsFactors = FALSE
    )

    passRassocs$BH = p.adjust(passRassocs$p_value, method = "BH")

    sig_hits = passRassocs %>% dplyr::filter(!is.na(BH) & BH <= FDR)

    if (nrow(sig_hits) == 0) {
      message(paste0("No BH-significant QTLs at FDR ", FDR, " in pass ", pass, "; stopping."))
      break
    }

    # one QTL per gene per pass: choose smallest p-value (peak)
    best_hits = sig_hits %>%
      dplyr::group_by(gene) %>%
      dplyr::slice_min(order_by = p_value, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()

    wings = lapply(seq_len(nrow(best_hits)), function(i) {
      g = best_hits$gene[i]
      peak = best_hits$marker[i]
      peakLOD = best_hits$LOD[i]
      compute_wings(gene = g, peak_marker = peak, peak_LOD = peakLOD,
                    pass_df = passRassocs, ct_vec = ct)
    })

    best_hits$LEFT  = vapply(wings, function(x) x$LEFT,  character(1))
    best_hits$RIGHT = vapply(wings, function(x) x$RIGHT, character(1))

    message(paste0("PASS ", pass, ": discovered ",
                   nrow(best_hits), " QTL(s) across ",
                   length(unique(best_hits$gene)), " gene(s)."))

    fassocs = rbind(fassocs, best_hits %>% dplyr::select(names(fassocs)))

    # residualize winners using the m1 model for this mode
    for (i in seq_len(nrow(best_hits))) {

      gene_i   = best_hits$gene[i]
      marker_i = best_hits$marker[i]

      g_index = which(suppGenes == gene_i)
      key     = uqKey[g_index]

      base_meta = longMeta[[key]]
      geno_vec  = jointGeno[[key]][, marker_i]
      pheno_vec = jointPhenoScaled[[g_index]][, 1]

      df = data.frame(
        study = base_meta$study,
        geno  = geno_vec,
        pheno = pheno_vec
      )
      if (!is.null(fixedEffect)) df[[fixedEffect]] = base_meta[[fixedEffect]]

      if (is.null(fixedEffect)) {
        f1 = as.formula("pheno ~ (1|study) + geno")
      } else {
        f1 = as.formula(paste0("pheno ~ (1|study) + ", fixedEffect, " + geno + geno:", fixedEffect))
      }

      qtlModel = lmer(f1, data = df,
                      control = lmerControl(check.conv.singular = "ignore"),
                      REML = FALSE)

      resid_vec = resid(qtlModel)

      rn = rownames(jointPhenoScaled[[g_index]])
      jointPhenoScaled[[g_index]] =
        matrix(rep(resid_vec, length(chromosomes)),
               ncol = length(chromosomes), byrow = FALSE,
               dimnames = list(rn, chromosomes))
    }

    pass = pass + 1
  }

  out = list()
  out[[mapMode]] = fassocs
  out
}

# Performs assocation tests between all markers and genes 
# essentially a truncated and reformatted jointMarkerTest
# Does not include any form of FDR correction
# Also a memory-unoptimized hog
rawAssocs = function(studies, targetMarkers, unionBg, smetadata, suppGenes, 
                     fixedEffect, fn, correctBG=FALSE, correctGRM=FALSE,
                     cores=1, chromosomes=NULL) {
  
  print('PREPARING MAPPING OBJECTS AND DEALING WITH GENE SUPPORT')
  
  # Create the directory if it does not already exist
  outPath = paste0(imed, fn, '/')
  if (!dir.exists(outPath)) {
    dir.create(outPath, showWarnings = TRUE, recursive = FALSE)
  }
  
  # Get true marker positions
  trueMarkers = colnames(studies[[1]]$gdata)
  
  # Define default 16 S. cerevisiae chromosomes if none are provided
  if (is.null(chromosomes)) {
    chromosomes = paste0("chr", as.roman(1:16))
  }
  
  # Get a list of which studies support every gene
  suppList = lapply(suppGenes, function(gene) {
    studies_containing_gene = sapply(smetadata$genes, function(gene_list) {
      gene %in% gene_list
    })
    smetadata$study[studies_containing_gene]
  }) %>% `names<-`(suppGenes)
  
  # Get unique study combinations and their keys
  uqCombos = unique(suppList)
  uqKey = match(suppList, uqCombos)
  
  # Create per-gene metadata object (gene-long list with segregants as rows)
  longMeta = mapply(function(x, y) {
    rep(y, nrow(x$gdata))
  }, studies, names(studies)) %>%
    unlist() %>%
    unname() %>%
    as.data.frame() %>%
    `names<-`(c('study')) %>%
    left_join(smetadata %>% 
                select(study, all_of(fixedEffect)), by = "study") %>%
    mutate(geno = NA)
  
  # Create gene-specific metadata using longMeta for subsetting
  longMeta = lapply(uqCombos, function(combo) {
    longMeta %>% filter(study %in% combo)
  })
  
  # Create a joint genotype object across all studies
  jointGeno = lapply(uqCombos, function(combo) {
    names(studies) %>%
      .[. %in% combo] %>%
      lapply(function(study_name) studies[[study_name]]$gdata) %>%
      do.call(rbind, .) %>%
      return()
  })
  
  # Chrom indexing key for all markers (e.g., "ChrI")
  cc = do.call('rbind', strsplit(colnames(jointGeno[[1]]), ':'))[,1]
  
  # Create phenotypes object (gene long list of Seg x Chrom matrices)
  jointPhenoScaled = lapply(suppGenes, function(gene) {
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
  
  if (correctBG) {
    print('CORRECTING FOR STRONG GENETIC EFFECTS ON OTHER CHROMOSOMES')
    pb = txtProgressBar(min = 0, max = length(uqKey), style = 3)
    jointPhenoScaled = mapply(function(key, pheno, index) {
      result = adjust_for_background_QTLs(
        jointGeno[[key]], pheno, unionBg, chromosomes)
      setTxtProgressBar(pb, index)
      return(result)
    }, uqKey, jointPhenoScaled,
    seq_along(uqKey), SIMPLIFY = FALSE)
    cat("\n")
  }
  
  if (correctGRM) {
    print('CORRECTING FOR RELATEDNESS ON OTHER CHROMOSOMES')
    pb = txtProgressBar(min = 0, max = length(uqKey), style = 3)
    jointPhenoScaled = mapply(function(key, pheno, index) {
      result = calculate_BLUPs(jointGeno[[key]], pheno)
      setTxtProgressBar(pb, index)
      return(result)
    }, uqKey, jointPhenoScaled,
    seq_along(uqKey), SIMPLIFY = FALSE)
    cat("\n")
  }
  
  
  # Scale genotypes to control for variance at each locus gene-specifically
  jointGeno = lapply(jointGeno, scale)
  
  print('FINALIZING ASSOCIATION PREPARATIONS')
  
  # Chrom indexing key for targeted markers (e.g., "ChrI")
  ct = do.call('rbind', strsplit(targetMarkers, ':'), TRUE)[,1]
  
  # Determine model terms based on fixedEffect
  model_terms = if (fixedEffect == "analyte") {
    c(fixedEffect, 'geno', 'ixn', 'mRNA', 'protein')
  } else {
    c(fixedEffect, 'geno', 'ixn', 'SC', 'YNB2G')
  }
  
  # Loop for all markers
  print('BEGINNING UNFILTERED ASSOCIATIONS AT SELECT MARKERS')
  mclapply(targetMarkers, function(marker) {
    
    # Get index of current marker and true marker index
    mI = which(targetMarkers == marker)
    mT = which(trueMarkers == marker)
    
    # Real associations (terms many matrices)
    rassocs = rep(list(
      list(
        gene = suppGenes,
        marker = as.character(rep(marker, length(suppGenes))),
        pval = as.numeric(rep(NA, length(suppGenes))),
        effect = as.numeric(rep(NA, length(suppGenes))),
        flavor = as.character(rep(NA, length(suppGenes))),
        samplesize = as.numeric(rep(NA, length(suppGenes)))
      ) %>% as.data.frame()
    ), length(model_terms)) %>% `names<-`(model_terms) %>%
      imap(function(df, name) {
        df %>% mutate(flavor = name)
      })
    
    # Sweep over the genes
    for (gene in seq_along(suppGenes)) {
      
      # Create the correct metadata objects
      perMeta = data.frame(
        study = longMeta[[uqKey[gene]]]$study,
        geno = jointGeno[[uqKey[gene]]][, marker],
        pheno = jointPhenoScaled[[gene]][, ct[mI]]
      )
      perMeta[[fixedEffect]] = longMeta[[uqKey[gene]]][[fixedEffect]]
      
      # Run the full model and ANOVA
      formula = as.formula(paste("pheno ~ (1|study) +", 
                                 fixedEffect, "+ geno + geno:", fixedEffect))
      datModel = lmer(formula, data = perMeta, 
                      control = lmerControl(check.conv.singular = "ignore"))
      datANOVA = anova(datModel)
      
      # Run the mRNA-only model and ANOVA
      mMeta = perMeta[perMeta[[fixedEffect]] == 
                        ifelse(fixedEffect == "analyte", "mRNA", "SC"),]
      mModel = lmer(pheno ~ (1|study) + geno, data = mMeta, control = 
                      lmerControl(check.conv.singular = "ignore"))
      mANOVA = anova(mModel)
      
      # Run the protein model and ANOVA
      pMeta = perMeta[perMeta[[fixedEffect]] == 
                        ifelse(fixedEffect == "analyte", "protein", "YNB2G"),]
      pModel = lmer(pheno ~ (1|study) + geno, data = pMeta, control =
                      lmerControl(check.conv.singular = "ignore"))
      pANOVA = anova(pModel)
      
      # Assign fixedEffect values
      rassocs[[fixedEffect]][gene, 'pval'] = datANOVA[fixedEffect, "Pr(>F)"]
      rassocs[[fixedEffect]][gene, 'effect'] = 
        fixef(datModel)[paste0(fixedEffect, ifelse(
          fixedEffect == "analyte", "protein", "YNB2G"))]
      rassocs[[fixedEffect]][gene, 'samplesize'] = nrow(perMeta)
      
      # Assign geno values
      rassocs$geno[gene, 'pval'] = datANOVA["geno", "Pr(>F)"]
      rassocs$geno[gene, 'effect'] = fixef(datModel)['geno']
      rassocs$geno[gene, 'samplesize'] = nrow(perMeta)
      
      # Assign ixn values
      rassocs$ixn[gene, 'pval'] = datANOVA[
        paste0(fixedEffect, ":geno"), "Pr(>F)"]
      rassocs$ixn[gene, 'effect'] = 
        fixef(datModel)[paste0(fixedEffect, ifelse(
          fixedEffect == "analyte", "protein:geno", "YNB2G:geno"))]
      rassocs$ixn[gene, 'samplesize'] = nrow(perMeta)
      
      # Assign mRNA or SC model
      rassocs[[ifelse(fixedEffect == "analyte", "mRNA", "SC")]][
        gene, 'pval'] = mANOVA["geno", "Pr(>F)"]
      rassocs[[ifelse(fixedEffect == "analyte", "mRNA", "SC")]][
        gene, 'effect'] = fixef(mModel)['geno']
      rassocs[[ifelse(fixedEffect == "analyte", "mRNA", "SC")]][
        gene, 'samplesize'] = nrow(mMeta)
      
      # Assign protein or YNB2G model
      rassocs[[ifelse(fixedEffect == "analyte", "protein", "YNB2G")]][
        gene, 'pval'] = pANOVA["geno", "Pr(>F)"]
      rassocs[[ifelse(fixedEffect == "analyte", "protein", "YNB2G")]][
        gene, 'effect'] = fixef(pModel)['geno']
      rassocs[[ifelse(fixedEffect == "analyte", "protein", "YNB2G")]][
        gene, 'samplesize'] = nrow(pMeta)
      
    }
    
    # Update the list of all associations
    write.csv(rassocs, file = paste0(
      outPath, "marker_", mT, "_of_", length(trueMarkers), ".csv"),
      row.names = FALSE)
    
    return(paste('completed associations for', marker))
    
  }, mc.cores = cores)
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

# Function to trace LOD and/or effect size across studies individually
# traceType = s for signed LOD, r for R-value, or b for both
studyTracer = function(studies, smeta, suppGenes, fn, 
                       plotHeight, plotWidth, traceType, fixef = "analyte") {
  
  # Define directory names for signed LOD and R-value
  signedLOD_dir = paste0(fn, "_signed_LOD_traces_", fixef)
  rValue_dir = paste0(fn, "_R_value_traces_", fixef)
  
  # Create the output directories if they don't exist
  if (!dir.exists(signedLOD_dir)) {
    dir.create(signedLOD_dir, recursive = TRUE)
  }
  if (!dir.exists(rValue_dir)) {
    dir.create(rValue_dir, recursive = TRUE)
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
  
  # Prepare data matrices for signed LOD scores and R-values
  sLODs = lapply(studies, function(study) {
    lod = fasterLOD(nrow(study$pheno), scale(study$pclean),
                    scale(study$gdata), betas = TRUE) 
    list(signedLOD = lod$LOD * sign(lod$r), rValue = lod$r)
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
    
    # Prepare the plot data for both signed LOD and R-value
    plotDatList = lapply(supported, function(dataset) {
      sLOD = sLODs[[dataset]]
      gene_row = rownames(sLOD$signedLOD) == rownames(suppMat)[i]
      if (any(gene_row)) {
        data.frame(
          marker = colnames(sLOD$signedLOD),
          signedLOD = sLOD$signedLOD[gene_row, ],
          rValue = sLOD$rValue[gene_row, ],
          gene = rownames(suppMat)[i],
          dataset = dataset
        )
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
    
    # Merge with metadata to get the analyte or medium information
    plotDat = merge(plotDat, smeta[, c("study", fixef)],
                    by.x = "dataset", by.y = "study")
    
    # Lock in the factor levels for the dataset
    plotDat$dataset = factor(plotDat$dataset, levels = smeta$study)
    
    # Plotting based on traceType
    plotTypes = if (traceType == "b") c("s", "r") else traceType
    
    for (type in plotTypes) {
      # Select the correct data based on the trace type
      plotDat$value = if (type == "s") plotDat$signedLOD else plotDat$rValue
      
      # Plot the results using ggplot2
      plotDat$marker = factor(plotDat$marker, levels = unique(plotDat$marker))
      plotDat$gene = factor(plotDat$gene, levels = suppGenes)
      
      plotTitle = if (type == "s") "Signed LOD score traces for" else 
        "R-value traces for"
      yAxisLabel = if (type == "s") "Signed LOD score" else 
        "R-value"
      
      p = ggplot(plotDat, aes(
        x = as.numeric(marker), 
        y = value, color = dataset, group = dataset)) +
        geom_vline(xintercept = chrStarts, 
                   linetype = "dashed", color = "black") +
        geom_line() +
        facet_wrap(as.formula(paste("~", fixef)), scales = "free_y", ncol = 1) +
        scale_color_manual(values = dataset_colors) +
        labs(title = paste(plotTitle, rownames(suppMat)[i]),
             x = "Genomic position",
             y = yAxisLabel,
             color = NULL) +  # Remove legend title
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
      
      # Set the file prefix and directory based on traceType
      if (type == "s") {
        file_prefix = "Signed_LOD_traces_for"
        dir_to_use = signedLOD_dir
      } else {
        file_prefix = "R_value_traces_for"
        dir_to_use = rValue_dir
      }
      
      # Save the plot with the appropriate file name
      ggsave(filename = paste0(dir_to_use, "/", file_prefix, "_", 
                               rownames(suppMat)[i], ".pdf"), 
             plot = p, width = plotWidth, height = plotHeight, units = "in")
    }
  }
  
  pbapply::closepb(pb)
}

# Plots key traces from unfiltered associations
tracePlotter = function(dat, fn = "plots", plotWidth = 7.5,
                        plotHeight = 5, r = NULL, n = NULL) {
  
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
  
  # Determine prefixes based on data content
  if ("mRNA" %in% colnames(dat) || "protein" %in% colnames(dat)) {
    all_prefixes = c("mRNA", "protein", "genotype", "interaction", "analyte")
    fixef = "analyte"
  } else if ("SC" %in% colnames(dat) || "YNB2G" %in% colnames(dat)) {
    all_prefixes = c("SC", "YNB2G", "genotype", "interaction", "medium")
    fixef = "medium"
  } else {
    stop("Unknown prefixes in the data columns")
  }
  
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
  
  # Define trace types to iterate over
  trace_types = c("smltp", "Allelic effect")
  
  # Apply the function to each gene using lapply
  lapply(dat, function(plot_data) {
    for (trace_type in trace_types) {
      
      formatted_trace_type = gsub(" ", "_", tolower(trace_type))
      plot_dir = paste0(fn, "_", formatted_trace_type, "_for_", fixef, "/")
      
      # Create directory if it doesn't exist
      if (!dir.exists(plot_dir)) {
        dir.create(plot_dir, recursive = TRUE)
      }
      
      # Determine the plot title and axis labels
      if (trace_type == "smltp") {
        y_label = bquote("Signed " ~ -log[10](italic(p)))
        plot_title = bquote("Signed " ~ -log[10](italic(p)) ~ 
                              "traces for" ~ .(plot_data$gene[1]))
      } else {
        y_label = trace_type
        plot_title = paste(trace_type, "traces for", plot_data$gene[1])
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
             x = "Genomic position", y = y_label, color = NULL) + 
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
      
      # Construct the filename with gene name and suffix
      plot_filename = paste0(plot_dir, plot_data$gene[1], "_", 
                             formatted_trace_type, "_for_", fixef, ".pdf")
      
      # Save the plot
      ggsave(filename = plot_filename, 
             plot = p, width = plotWidth, height = plotHeight, units = "in")
    }
  })
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


# Gets genetic correlation summary stats and histograms across genes 
mpSummarizer = function(data, directory, mode) {
  
  # Ensure the directory exists
  ensure_directory(directory)
  
  # Create a list of effect types with corresponding transformations
  if (mode == "analyte") {
    effects_list = list(
      mRNA_effects = lapply(data, function(x) x$mRNA.allelic_effect),
      Protein_effects = lapply(data, function(x) x$protein.allelic_effect),
      Difference_in_mRNA_and_protein_raw_effects = lapply(data, function(x) 
        x$mRNA.allelic_effect - x$protein.allelic_effect),
      Difference_in_mRNA_and_protein_absolute_effects = lapply(data, function(x) 
        abs(x$mRNA.allelic_effect - x$protein.allelic_effect)),
      Log10P_mRNA = lapply(data, function(x) -log10(x$mRNA.pval)),
      Log10P_Protein = lapply(data, function(x) -log10(x$protein.pval)),
      Log10P_mRNA_Protein = lapply(data, function(x)
        -log10(x$mRNA.pval / x$protein.pval))
    )
  } else if (mode == "medium") {
    effects_list = list(
      SC_effects = lapply(data, function(x)
        x$SC.allelic_effect),
      YNB2G_effects = lapply(data, function(x)
        x$YNB2G.allelic_effect),
      Difference_in_SC_and_YNB2G_raw_effects = lapply(data, function(x)
        x$SC.allelic_effect - x$YNB2G.allelic_effect),
      Difference_in_SC_and_YNB2G_absolute_effects = lapply(data, function(x)
        abs(x$SC.allelic_effect - x$YNB2G.allelic_effect)),
      Log10P_SC = lapply(data, function(x) -log10(x$SC.pval)),
      Log10P_YNB2G = lapply(data, function(x) -log10(x$YNB2G.pval)),
      Log10P_SC_YNB2G = lapply(data, function(x)
        -log10(x$SC.pval / x$YNB2G.pval))
    )
  } else {
    stop("Invalid mode specified. Please use 'analyte' or 'medium'.")
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
    
    # Determine plot title, file naming, and x-label based on effect type
    if (grepl("Log10P_", effect_type)) {
      if (grepl("ratio", effect_type)) {
        plot_title = bquote("-log"[10] ~ "(" ~ italic(p) ~ ") ratio")
        effect_name = "Log10_P_Value_Ratio"
        x_label = bquote("-log"[10] ~ "(" ~ italic(p) ~ ") ratio")
      } else {
        entity = if (mode == "analyte") {
          if (grepl("mRNA", effect_type)) "mRNA" else "protein"
        } else {
          if (grepl("SC", effect_type)) "SC" else "YNB2G"
        }
        plot_title = bquote("-log"[10] ~ .(entity) ~ italic(p) ~ "-values")
        effect_name = paste0("Log10_", entity, "_P_Values")
        x_label = bquote("-log"[10] ~ "(" ~ italic(p) ~ ")")
      }
    } else {
      plot_title = gsub("_", " ", effect_type)
      effect_name = effect_type
      if (grepl("Difference", effect_type)) {
        x_label = "Difference in allelic effect"
      } else {
        x_label = 'Allelic effect'
      }
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
mpLinkPlotter = function(data, filename, infix, mode, r = NULL, n = NULL) {
  # Create subdirectory for plots based on filename
  plot_subdir = gsub("\\.pdf$", "", filename)
  ensure_directory(plot_subdir)
  
  # Select and pivot data based on the mode
  if (mode == "analyte") {
    data_long = data %>%
      dplyr::select(
        pairID, marker, mRNA.allelic_effect, protein.allelic_effect) %>%
      tidyr::pivot_longer(cols = c(mRNA.allelic_effect, protein.allelic_effect), 
                          names_to = "effect_type", 
                          values_to = "effect_value")
  } else if (mode == "medium") {
    data_long = data %>%
      dplyr::select(
        pairID, marker, SC.allelic_effect, YNB2G.allelic_effect) %>%
      tidyr::pivot_longer(cols = c(SC.allelic_effect, YNB2G.allelic_effect), 
                          names_to = "effect_type", 
                          values_to = "effect_value")
  } else {
    stop("Invalid mode specified. Please use 'analyte' or 'medium'.")
  }
  
  # Truncate effect_type labels after the period
  data_long$effect_type = sub("\\..*", "", data_long$effect_type)
  
  # Use the titler function to format the x-axis labels
  data_long$effect_type = sapply(data_long$effect_type, 
                                 titler, r = r, n = n, u = FALSE)
  
  # Determine the title prefix based on mode
  title_prefix = if (mode == "analyte") "mRNA and protein effects" else 
    "SC and YNB2G effects"
  
  # Determine the title infix based on the infix argument
  title_infix = if (infix == "ixn") {
    "of interaction QTLs"
  } else if (infix == "geno") {
    "of genotype QTLs"
  } else {
    infix
  }
  
  # Loop over unique markers
  unique_markers = unique(data_long$marker)
  
  for (marker in unique_markers) {
    marker_data = data_long %>% filter(marker == !!marker)
    
    plot_title = paste(title_prefix, "for", marker)
    
    # Create the plot
    p = ggplot(marker_data, aes(x = effect_type, 
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
      labs(x = NULL, y = "Allelic effect", title = plot_title)
    
    # Save the plot
    marker_filename = gsub(":", "_", marker)  # Replace ':' with '_'
    ggsave(paste0(plot_subdir, "/", marker_filename, ".pdf"), 
           plot = p, width = 7.5, height = 5)
  }
  
  # Create a combined plot for all markers
  combined_plot_title = paste(title_prefix, title_infix, "for all markers")
  
  p_combined = ggplot(data_long, aes(x = effect_type, 
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
    labs(x = "Effect type", y = "Allelic effect", title = combined_plot_title)
  
  # Save the combined plot
  ggsave(paste0(plot_subdir, "/all_markers.pdf"), 
         plot = p_combined, width = 7.5, height = 5)
}

# --- DRILL: now also stores p-values (LRT like your mapper) -------------------
drillQTL = function(studies,
                    dfQTL,
                    unionBg,
                    smetadata,
                    suppGenes,
                    fixedEffect = NULL,
                    correctBG = FALSE,
                    correctGRM = FALSE,
                    cores = 1,
                    chromosomes = NULL) {

  stopifnot(is.list(studies), length(studies) >= 1)
  stopifnot(is.data.frame(dfQTL))
  stopifnot(all(c("gene", "marker") %in% colnames(dfQTL)))

  if (!requireNamespace("lme4", quietly = TRUE)) stop("Package 'lme4' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")

  dfQTL$gene   = as.character(dfQTL$gene)
  dfQTL$marker = as.character(dfQTL$marker)

  if (is.null(chromosomes)) chromosomes = paste0("chr", as.roman(1:16))

  if (!is.null(suppGenes)) {
    suppGenes = as.character(suppGenes)
    dfQTL = dfQTL[dfQTL$gene %in% suppGenes, , drop = FALSE]
  }
  if (nrow(dfQTL) == 0) return(list())

  trueMarkers = colnames(studies[[1]]$gdata)
  keep_pairs = dfQTL$marker %in% trueMarkers
  if (!all(keep_pairs)) {
    missing_markers = unique(dfQTL$marker[!keep_pairs])
    warning("Dropping ", sum(!keep_pairs),
            " QTL rows whose marker is not in gdata. Example: ",
            paste(head(missing_markers, 3), collapse = ", "))
    dfQTL = dfQTL[keep_pairs, , drop = FALSE]
  }
  if (nrow(dfQTL) == 0) return(list())

  suppGenes_use = unique(dfQTL$gene)

  suppList = lapply(suppGenes_use, function(gene) {
    studies_containing_gene =
      sapply(smetadata$genes, function(gene_list) gene %in% gene_list)
    smetadata$study[studies_containing_gene]
  })
  names(suppList) = suppGenes_use

  uqCombos = unique(suppList)
  uqKey = match(suppList, uqCombos)

  baseLongMeta =
    mapply(function(x, y) rep(y, nrow(x$gdata)), studies, names(studies)) |>
    unlist(use.names = FALSE) |>
    data.frame(study = _, stringsAsFactors = FALSE)

  if (!is.null(fixedEffect)) {
    if (!fixedEffect %in% colnames(smetadata)) {
      stop("fixedEffect '", fixedEffect, "' is not a column in smetadata.")
    }
    baseLongMeta =
      dplyr::left_join(
        baseLongMeta,
        smetadata |> dplyr::select(study, dplyr::all_of(fixedEffect)),
        by = "study"
      )
  }

  baseLongMeta$geno  = NA_real_
  baseLongMeta$pheno = NA_real_

  longMeta = lapply(uqCombos, function(combo) {
    baseLongMeta |> dplyr::filter(study %in% combo)
  })

  jointGeno = lapply(uqCombos, function(combo) {
    st_names = names(studies)[names(studies) %in% combo]
    mats = lapply(st_names, function(study_name) studies[[study_name]]$gdata)
    do.call(rbind, mats)
  })

  jointPhenoScaled = lapply(suppGenes_use, function(gene) {
    pheno_values = unlist(lapply(suppList[[gene]], function(studyName) {
      if (!studyName %in% names(studies)) return(NULL)
      st = studies[[studyName]]
      if (!gene %in% colnames(st$pclean)) return(NULL)
      st$pclean[, gene]
    }))
    if (length(pheno_values) == 0) return(NULL)

    mat = matrix(rep(pheno_values, length(chromosomes)),
                 ncol = length(chromosomes), byrow = FALSE)
    colnames(mat) = chromosomes
    rownames(mat) = names(pheno_values)
    mat
  })
  names(jointPhenoScaled) = suppGenes_use

  ok_gene = !vapply(jointPhenoScaled, is.null, logical(1))
  if (!all(ok_gene)) {
    dropped = names(jointPhenoScaled)[!ok_gene]
    warning("Dropping ", sum(!ok_gene),
            " genes with no phenotype values. Example: ",
            paste(head(dropped, 3), collapse = ", "))
    jointPhenoScaled = jointPhenoScaled[ok_gene]
    suppGenes_use = names(jointPhenoScaled)
    suppList = suppList[suppGenes_use]
    uqCombos = unique(suppList)
    uqKey = match(suppList, uqCombos)

    longMeta = lapply(uqCombos, function(combo) {
      baseLongMeta |> dplyr::filter(study %in% combo)
    })
    jointGeno = lapply(uqCombos, function(combo) {
      st_names = names(studies)[names(studies) %in% combo]
      mats = lapply(st_names, function(study_name) studies[[study_name]]$gdata)
      do.call(rbind, mats)
    })
  }

  if (correctBG) {
    if (!exists("adjust_for_background_QTLs", mode = "function")) {
      stop("correctBG=TRUE but adjust_for_background_QTLs() is not found.")
    }
    message("CORRECTING FOR STRONG GENETIC EFFECTS ON OTHER CHROMOSOMES")
    pb = txtProgressBar(min = 0, max = length(uqKey), style = 3)
    jointPhenoScaled = mapply(function(key, pheno, index) {
      out = adjust_for_background_QTLs(jointGeno[[key]], pheno, unionBg, chromosomes)
      out = as.matrix(out)
      if (ncol(out) == length(chromosomes)) colnames(out) = chromosomes
      setTxtProgressBar(pb, index)
      out
    }, uqKey, jointPhenoScaled, seq_along(uqKey), SIMPLIFY = FALSE)
    close(pb); cat("\n")
  }

  if (correctGRM) {
    if (!exists("calculate_BLUPs", mode = "function")) {
      stop("correctGRM=TRUE but calculate_BLUPs() is not found.")
    }
    message("CORRECTING FOR RELATEDNESS ON OTHER CHROMOSOMES")
    pb = txtProgressBar(min = 0, max = length(uqKey), style = 3)
    jointPhenoScaled = mapply(function(key, pheno, index) {
      out = calculate_BLUPs(jointGeno[[key]], pheno)
      out = as.matrix(out)
      if (ncol(out) == length(chromosomes)) colnames(out) = chromosomes
      setTxtProgressBar(pb, index)
      out
    }, uqKey, jointPhenoScaled, seq_along(uqKey), SIMPLIFY = FALSE)
    close(pb); cat("\n")
  }

  jointGeno = lapply(jointGeno, scale)

  fixed_levels = NULL
  if (!is.null(fixedEffect)) {
    fixed_levels = sort(unique(as.character(baseLongMeta[[fixedEffect]])))
    fixed_levels = fixed_levels[!is.na(fixed_levels)]
  }

  ctrl = lme4::lmerControl(check.conv.singular = "ignore")

  f0_geno = as.formula("pheno ~ (1|study)")
  f1_geno = as.formula("pheno ~ (1|study) + geno")

  f0_int = NULL
  f1_int = NULL
  if (!is.null(fixedEffect)) {
    f0_int = as.formula(paste0("pheno ~ (1|study) + ", fixedEffect, " + geno"))
    f1_int = as.formula(paste0("pheno ~ (1|study) + ", fixedEffect, " + geno + geno:", fixedEffect))
  }

  dfPairs = unique(dfQTL[dfQTL$gene %in% suppGenes_use, c("gene", "marker"), drop = FALSE])
  if (nrow(dfPairs) == 0) return(list())

  fit_one_pair = function(gene, marker) {
    g_index = match(gene, suppGenes_use)
    if (is.na(g_index)) return(list(gene = gene, marker = marker, error = "gene_not_in_suppGenes_use"))
    key = uqKey[g_index]

    base_meta = longMeta[[key]]

    if (!marker %in% colnames(jointGeno[[key]])) {
      return(list(gene = gene, marker = marker, error = "marker_not_in_jointGeno"))
    }

    chr = strsplit(marker, ":", fixed = TRUE)[[1]][1]
    chr_idx = match(chr, chromosomes)
    if (is.na(chr_idx)) return(list(gene = gene, marker = marker, error = paste0("bad_chr: ", chr)))

    ph = jointPhenoScaled[[g_index]]
    if (is.null(ph)) return(list(gene = gene, marker = marker, error = "pheno_is_NULL"))

    ph = as.matrix(ph)
    if (ncol(ph) < chr_idx) {
      return(list(gene = gene, marker = marker, error = paste0("pheno_has_only_", ncol(ph), "_cols")))
    }

    df = base_meta
    df$geno  = jointGeno[[key]][, marker]
    df$pheno = ph[, chr_idx]

    out = list(gene = gene, marker = marker)

    # GENO model + pval (LRT)
    m0 = lme4::lmer(f0_geno, data = df, control = ctrl, REML = FALSE)
    m1 = lme4::lmer(f1_geno, data = df, control = ctrl, REML = FALSE)
    out$geno = m1
    out$geno_pval = tryCatch(stats::anova(m0, m1)[2, "Pr(>Chisq)"], error = function(e) NA_real_)

    # INTERACTION model + pval (LRT)
    if (!is.null(fixedEffect)) {
      m_main = lme4::lmer(f0_int, data = df, control = ctrl, REML = FALSE)
      m_int  = lme4::lmer(f1_int, data = df, control = ctrl, REML = FALSE)
      out$interaction = m_int
      out$interaction_pval = tryCatch(stats::anova(m_main, m_int)[2, "Pr(>Chisq)"], error = function(e) NA_real_)

      # per-level subset geno models + pvals (LRT within subset)
      for (lv in fixed_levels) {
        df_sub = df[as.character(df[[fixedEffect]]) == lv, , drop = FALSE]
        if (nrow(df_sub) < 5) {
          out[[lv]] = NULL
          out[[paste0(lv, "_pval")]] = NA_real_
        } else {
          m0s = lme4::lmer(f0_geno, data = df_sub, control = ctrl, REML = FALSE)
          m1s = lme4::lmer(f1_geno, data = df_sub, control = ctrl, REML = FALSE)
          out[[lv]] = m1s
          out[[paste0(lv, "_pval")]] = tryCatch(stats::anova(m0s, m1s)[2, "Pr(>Chisq)"], error = function(e) NA_real_)
        }
      }
    }

    out
  }

  message("FITTING LMER MODELS FOR PROVIDED (gene, marker) PAIRS")
  nPairs = nrow(dfPairs)
  idx = seq_len(nPairs)
  pb = txtProgressBar(min = 0, max = nPairs, style = 3)

  if (cores > 1) {
    res = parallel::mclapply(idx, function(i) fit_one_pair(dfPairs$gene[i], dfPairs$marker[i]),
                             mc.cores = cores)
    setTxtProgressBar(pb, nPairs)
  } else {
    res = lapply(idx, function(i) {
      ans = fit_one_pair(dfPairs$gene[i], dfPairs$marker[i])
      setTxtProgressBar(pb, i)
      ans
    })
  }

  close(pb); cat("\n")
  names(res) = paste0(dfPairs$gene, "_", dfPairs$marker)
  res
}


# summarize the QTL models from the targeted search into a df
modelQTL_to_df = function(modelQTL) {
  if (!requireNamespace("lme4", quietly = TRUE)) stop("Package 'lme4' is required.")

  is_lmer = function(m) inherits(m, "lmerMod")

  get_geno_effect = function(m) {
    if (!is_lmer(m)) return(NA_real_)
    fx = lme4::fixef(m)
    if (!"geno" %in% names(fx)) return(NA_real_)
    unname(fx["geno"])
  }

  get_interaction_effect = function(m) {
    if (!is_lmer(m)) return(NA_real_)
    fx = lme4::fixef(m)
    ixn_names = grep("(^geno:|:geno$)", names(fx), value = TRUE)
    if (length(ixn_names) == 0) return(NA_real_)
    vals = fx[ixn_names]
    if (length(vals) == 1) return(unname(vals))
    unname(vals[which.max(abs(vals))])
  }

  infer_levels = function(x) {
    nms = names(x)
    # keep only subset models (exclude pval slots)
    levs = setdiff(nms, c("gene", "marker", "geno", "interaction", "geno_pval", "interaction_pval"))
    levs = levs[!grepl("_pval$", levs)]
    levs
  }

  subset_levels = character(0)
  for (x in modelQTL) {
    lv = infer_levels(x)
    if (length(lv) > 0) { subset_levels = lv; break }
  }

  rows = vector("list", length(modelQTL))

  for (i in seq_along(modelQTL)) {
    x = modelQTL[[i]]

    row = list(
      gene = x$gene,
      marker = x$marker,

      geno_effect = get_geno_effect(x$geno),
      geno_pval = if (!is.null(x$geno_pval)) as.numeric(x$geno_pval) else NA_real_,

      interaction_effect = get_interaction_effect(x$interaction),
      interaction_pval = if (!is.null(x$interaction_pval)) as.numeric(x$interaction_pval) else NA_real_
    )

    for (lv in subset_levels) {
      row[[paste0(lv, "_effect")]] = get_geno_effect(x[[lv]])
      row[[paste0(lv, "_pval")]] =
        if (!is.null(x[[paste0(lv, "_pval")]])) as.numeric(x[[paste0(lv, "_pval")]]) else NA_real_
    }

    rows[[i]] = row
  }

  df = as.data.frame(do.call(rbind, lapply(rows, as.data.frame)), stringsAsFactors = FALSE)

  num_cols = setdiff(names(df), c("gene", "marker"))
  for (cc in num_cols) df[[cc]] = suppressWarnings(as.numeric(df[[cc]]))
  df
}

binTwoLevel = function(df, level1, level2, alpha = 0.05,
  pval_suffix = "_pval", effect_suffix = "_effect") {

  l1_p = paste0(level1, pval_suffix)
  l2_p = paste0(level2, pval_suffix)
  l1_e = paste0(level1, effect_suffix)
  l2_e = paste0(level2, effect_suffix)

  sig1 = !is.na(df[[l1_p]]) & df[[l1_p]] < alpha
  sig2 = !is.na(df[[l2_p]]) & df[[l2_p]] < alpha

  e1 = df[[l1_e]]
  e2 = df[[l2_e]]

  # binning direction

  same_dir =
    !is.na(e1) & !is.na(e2) &
    e1 != 0 & e2 != 0 &
    sign(e1) == sign(e2)

  opp_dir =
    !is.na(e1) & !is.na(e2) &
    e1 != 0 & e2 != 0 &
    sign(e1) != sign(e2)

  # all bins
  bin = rep(NA_character_, nrow(df))

  bin[!sig1 & !sig2] = "nd in singles"
  bin[ sig1 & !sig2] = paste0(level1, " specific")
  bin[!sig1 &  sig2] = paste0(level2, " specific")
  bin[ sig1 &  sig2 & same_dir] = "same direction"
  bin[ sig1 &  sig2 & opp_dir]  = "line cross"

  # edge case: both significant but one/both effects are 0 or NA
  bin[sig1 & sig2 & is.na(bin)] = "Both significant (direction NA)"

  df$bin = bin

  # quick count table
  counts = as.data.frame(table(df$bin, useNA = "ifany"), stringsAsFactors = FALSE)
  names(counts) = c("bin", "count")
  counts = counts[order(-counts$count), , drop = FALSE]
  rownames(counts) = NULL

  list(df = df, counts = counts)
}

aggregate_chunks = function(inp_dir, out_rda) {
  files = list.files(inp_dir, pattern = "^chunk[0-9]+\\.rda$", full.names = TRUE)
  if (length(files) == 0) stop("No chunk*.rda in: ", inp_dir)

  # numeric sort by chunk index
  idx = as.integer(sub("^chunk([0-9]+)\\.rda$", "\\1", basename(files)))
  files = files[order(idx)]

  agg = NULL
  for (f in files) {
    e = new.env(parent = emptyenv())
    load(f, envir = e)

    x = e$chunkAssocs
    if (is.null(x)) stop("chunkAssocs missing in: ", f)

    # NEW: if chunkAssocs is a length-1 list (your [[1]] wrapper), unwrap it
    if (is.list(x) && length(x) == 1L) x = x[[1]]

    if (!is.data.frame(x)) {
      stop("chunkAssocs in ", f, " is not a data.frame after unwrapping; class=",
           paste(class(x), collapse = "/"))
    }

    if (is.null(agg)) agg = x else agg = dplyr::bind_rows(agg, x)
  }

  out_name = tools::file_path_sans_ext(basename(out_rda))
  assign(out_name, agg, envir = .GlobalEnv)

  out_dir = dirname(out_rda)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  save(list = out_name, file = out_rda)
  message("Wrote ", out_rda, " (object: ", out_name, ", rows: ", nrow(agg), ")")
}

# Build augmented QTL df (original + added cols + attached modelQTL) and save it
attach_drill_info_and_save = function(
  dfQTL,
  studies,
  unionBg,
  smetadata,
  suppGenes = NULL,
  fixedEffect = NULL,
  correctBG = FALSE,
  correctGRM = FALSE,
  cores = 1,
  chromosomes = NULL,
  do_bin = FALSE,
  bin_level1 = NULL,
  bin_level2 = NULL,
  bin_alpha = 0.05,
  out_name,
  out_file
) {
  df_in = dfQTL
  df_in$gene = as.character(df_in$gene)
  df_in$marker = as.character(df_in$marker)

  modelQTL = drillQTL(
    studies = studies,
    dfQTL = df_in,
    unionBg = unionBg,
    smetadata = smetadata,
    suppGenes = suppGenes,
    fixedEffect = fixedEffect,
    correctBG = correctBG,
    correctGRM = correctGRM,
    cores = cores,
    chromosomes = chromosomes
  )

  df_out = df_in

  if (length(modelQTL) > 0) {
    df_add = modelQTL_to_df(modelQTL)

    df_out = merge(
      x = df_in,
      y = df_add,
      by = c("gene", "marker"),
      all.x = TRUE,
      sort = FALSE
    )
  }

  if (isTRUE(do_bin)) {
    if (is.null(bin_level1) || is.null(bin_level2)) {
      stop("do_bin=TRUE requires bin_level1 and bin_level2.")
    }

    b = binTwoLevel(
      df_out,
      level1 = bin_level1,
      level2 = bin_level2,
      alpha = bin_alpha
    )

    df_out = b$df
    attr(df_out, "bin_counts") = b$counts
  }

  assign(out_name, df_out, envir = .GlobalEnv)
  save(list = out_name, file = out_file)

  rm(modelQTL)
  invisible(gc())

  df_out
}

##########################
# mashr analysis section #
##########################

load_rda_object = function(path, object_name = NULL) {
  env = new.env(parent = emptyenv())
  nm = load(path, envir = env)

  if (!is.null(object_name)) {
    if (!(object_name %in% nm)) {
      stop("Object '", object_name, "' not found in ", path)
    }
    return(env[[object_name]])
  }

  if (length(nm) != 1L) {
    stop(
      "load_rda_object() expected exactly one object in ",
      path,
      " when object_name is NULL, but found: ",
      paste(nm, collapse = ", ")
    )
  }

  env[[nm]]
}

write_tsv_if_needed = function(dt, path, overwrite = FALSE, sep = "\t") {
  if (!file.exists(path) || isTRUE(overwrite)) {
    data.table::fwrite(dt, path, sep = sep)
  }
  invisible(path)
}

write_lines_if_needed = function(x, path, overwrite = FALSE) {
  if (!file.exists(path) || isTRUE(overwrite)) {
    writeLines(x, path)
  }
  invisible(path)
}

prepare_output_dirs = function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

fit_gene_in_study = function(study, gene) {
  if (!(gene %in% colnames(study$pclean))) {
    m = ncol(study$gscaled)
    return(list(
      beta = rep(NA_real_, m),
      se = rep(NA_real_, m),
      n = NA_integer_
    ))
  }

  y = as.numeric(study$pclean[, gene])
  G = study$gscaled

  ok = is.finite(y)
  y = y[ok]
  G = G[ok, , drop = FALSE]

  n = length(y)
  m = ncol(G)

  if (n < 3L) {
    return(list(
      beta = rep(NA_real_, m),
      se = rep(NA_real_, m),
      n = n
    ))
  }

  numer = as.numeric(crossprod(y, G))
  denom = colSums(G^2)
  beta = numer / denom

  ssy = sum(y^2)
  rss = ssy - (numer^2 / denom)
  se = sqrt((rss / (n - 2)) / denom)

  beta[!is.finite(beta)] = NA_real_
  se[!is.finite(se)] = NA_real_

  list(beta = beta, se = se, n = n)
}

build_gene_mats = function(gene, gene_fits, study_names, markers) {
  fits = gene_fits[[gene]]

  Bhat = do.call(cbind, lapply(fits, `[[`, "beta"))
  Shat = do.call(cbind, lapply(fits, `[[`, "se"))

  colnames(Bhat) = study_names
  colnames(Shat) = study_names
  rownames(Bhat) = markers
  rownames(Shat) = markers

  ns = vapply(fits, `[[`, integer(1), "n")

  list(Bhat = Bhat, Shat = Shat, n = ns)
}

build_shared_ld_reference = function(overwrite = FALSE) {
  path = outfile("shared_ld_reference.rds")

  load_or_compute_rds(
    path,
    expr = {
      message("Loading studies once to compute shared LD reference...")
      studies_file = infile("studies.rda")
      studies_all = load_rda_object(studies_file, "studies")

      if ("Albert_2018" %in% names(studies_all)) {
        studies_all = studies_all[names(studies_all) != "Albert_2018"]
      }

      studies_all = lapply(names(studies_all), function(x) {
        harmonize_study_rowkeys(studies_all[[x]], study_name = x)
      }) |> `names<-`(names(studies_all))

      ref_study = choose_ld_reference_study(studies_all)
      message("Using ", ref_study, " as the LD reference study.")

      G = studies_all[[ref_study]]$gdata
      if (is.null(G)) {
        stop("Reference study ", ref_study, " does not contain gdata.")
      }

      ld_r = stats::cor(G, use = "pairwise.complete.obs")

      out = list(
        study = ref_study,
        markers = colnames(G),
        ld_r = ld_r
      )

      rm(studies_all, G, ld_r)
      gc()

      out
    },
    overwrite = overwrite,
    label = "shared LD reference"
  )
}

align_ld_matrix = function(ld_reference, markers) {
  if (!all(markers %in% ld_reference$markers)) {
    missing_markers = setdiff(markers, ld_reference$markers)
    stop(
      "Some markers are missing from the shared LD matrix: ",
      paste(head(missing_markers, 10L), collapse = ", ")
    )
  }

  ld_reference$ld_r[markers, markers, drop = FALSE]
}

ld_prune_gene_rows = function(
  Bhat,
  Shat,
  ld_submat,
  abs_cor_threshold = 0.95,
  max_keep = Inf
) {
  stopifnot(identical(rownames(Bhat), rownames(Shat)))
  stopifnot(identical(rownames(Bhat), rownames(ld_submat)))
  stopifnot(identical(colnames(ld_submat), rownames(ld_submat)))

  Z = Bhat / Shat
  row_score = apply(abs(Z), 1, max, na.rm = TRUE)
  row_score[!is.finite(row_score)] = -Inf

  ord = order(row_score, decreasing = TRUE)
  marker_names = rownames(Bhat)
  keep = integer(0)
  blocked = rep(FALSE, length(marker_names))
  names(blocked) = marker_names

  for (idx in ord) {
    marker_i = marker_names[idx]

    if (isTRUE(blocked[marker_i])) {
      next
    }

    keep = c(keep, idx)

    if (length(keep) >= max_keep) {
      break
    }

    linked = abs(ld_submat[marker_i, marker_names]) >= abs_cor_threshold
    linked[!is.finite(linked)] = FALSE
    blocked[marker_names[linked]] = TRUE
  }

  sort(unique(keep))
}

print('SOURCE OF CONSTANTS AND FUNCTIONS WAS SUCCESSFUL')
