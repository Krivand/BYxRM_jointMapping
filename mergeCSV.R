# SOURCE FXNS AND CONSTANTS
source(paste0("/home/albertf/vandy147/p05/",
              "remap_functions.R"))

# Save the combined dataframe as an RDA file
roundupCSV(inputPath = paste0(imed, "joint_hotspot_QTLs"),
           outputPath = imed, rName = "jhQTLs", cleanUp = FALSE)
roundupCSV(inputPath = paste0(imed, "unfiltered_assocs"),
           outputPath = imed, rName = "unfilteredAssocs", cleanUp = FALSE)
