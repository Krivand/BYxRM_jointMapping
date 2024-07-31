# SOURCE FXNS AND CONSTANTS
source(paste0("/home/albertf/vandy147/p05/",
              "remap_functions.R"))

# List all files in the directory that start with "joint_QTLs_"
files = list.files(imed, pattern = "^joint_QTLs_.*\\.csv$", full.names = TRUE)

# Initialize an empty dataframe
jointQTLs = data.frame()

# Loop through each file and read its contents to create a single dataframe
for (file in files) {
  data = read.csv(file)
  jointQTLs = bind_rows(jointQTLs, data)
}

# Display the combined dataframe
print(jointQTLs)
unique(jointQTLs$gene) %>% length()
unique(jointQTLs$marker) %>% length()
table(jointQTLs$flavor)
unique_geno = jointQTLs[jointQTLs$flavor == 'geno',]
unique_ixn = jointQTLs[jointQTLs$flavor == 'ixn',]

# Initialize counters
geno_only_count = 0
ixn_only_count = 0
both_count = 0

# Get unique markers
unique_markers = unique(jointQTLs$marker)

# Iterate over each unique marker
for (marker in unique_markers) {
  # Subset data for the current marker
  subset_data = jointQTLs %>% filter(marker == !!marker)
  
  # Get unique genes for the current marker
  unique_genes = unique(subset_data$gene)
  
  # Iterate over each unique gene
  for (gene in unique_genes) {
    # Subset data for the current gene and marker
    gene_data = subset_data %>% filter(gene == !!gene)
    
    # Get counts for geno and ixn QTLs
    geno_count = nrow(gene_data %>% filter(flavor == "geno"))
    ixn_count = nrow(gene_data %>% filter(flavor == "ixn"))
    
    # Update counters based on counts
    if (geno_count > 0 & ixn_count == 0) {
      geno_only_count = geno_only_count + 1
    } else if (geno_count == 0 & ixn_count > 0) {
      ixn_only_count = ixn_only_count + 1
    } else if (geno_count > 0 & ixn_count > 0) {
      both_count = both_count + 1
    }
  }
}

# Print the counts
cat("Number of geno QTLs without an ixn QTL:", geno_only_count, "\n")
cat("Number of ixn QTLs without a geno QTL:", ixn_only_count, "\n")
cat("Number of QTLs with both geno and ixn QTLs:", both_count, "\n")
