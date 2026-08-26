########################################
########################################
## GENETICS REVIEW FINAL OUTPUT BUNDLE ##
## WORKBOOK PACKAGING ENTRY POINT      ##
########################################
########################################

Sys.setenv(REVIEW_POWER_NVALID = "1000")

source("/PHShome/kv1026/p05/iReviewWorkbookBundle.R")

allpairs_env = new.env(parent = globalenv())
sys.source(
  "/PHShome/kv1026/p05/iReviewPowerAllPairs.R",
  envir = allpairs_env
)
allpairs_validation_env = new.env(parent = globalenv())
sys.source(
  "/PHShome/kv1026/p05/iReviewPowerAllPairsValidate.R",
  envir = allpairs_validation_env
)

source("/PHShome/kv1026/p05/iReviewWorkbookFinalize.R")
source("/PHShome/kv1026/p05/iReviewPrePowerSupportFinalize.R")
source("/PHShome/kv1026/p05/iReviewPowerAllPairsFinalize.R")
