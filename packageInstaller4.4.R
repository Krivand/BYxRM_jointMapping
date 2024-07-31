# Define the library path
personal_lib_path <- "~/R/x86_64-pc-linux-gnu-library/4.4"

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
         "lme4", "lmerTest", "patchwork", "plyr", "progress", 
         "readxl", "rrBLUP", "tidyverse", "viridis")

# Install using correct package installation tools
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message("Installing ", pkg, "...")
    if (pkg %in% c("GEOquery", "GenomicRanges")) {
      BiocManager::install(pkg, lib = personal_lib_path, quietly = TRUE)
    } else {
      install.packages(pkg, lib = personal_lib_path, quietly = TRUE)
    }
  }
}

# Load packages and remove list
invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE))
rm(pkgs)
