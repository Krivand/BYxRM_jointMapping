# Source functions and constants
source("/home/albertf/vandy147/p05/remap_functions.R")

# Parse command-line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Four arguments must be supplied (inputPath, outputPath, rName, cleanUp).", call. = FALSE)
}
inputPath = args[1]
outputPath = args[2]
rName = args[3]
cleanUp = as.logical(args[4])

# Call the roundupCSV function
roundupCSV(inputPath = inputPath, outputPath = outputPath, rName = rName, cleanUp = cleanUp)
