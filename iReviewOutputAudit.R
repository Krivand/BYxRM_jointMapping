########################################
########################################
## GENETICS REVIEW FINAL OUTPUT BUNDLE ##
## WORKBOOK PACKAGING ENTRY POINT      ##
########################################
########################################

Sys.setenv(REVIEW_POWER_NVALID = "1000")

source("/PHShome/kv1026/p05/iReviewWorkbookBundle.R")

# iReviewWorkbookBundle.R intentionally removes nonfinal reviewer staging
# outputs after it builds the first-pass workbooks. The all-pair power
# extension validates against the compact Albert-reference and Figure 3
# calibration tables, so regenerate only those inexpensive power prerequisites
# before running the all-pair analysis. Keep them in isolated environments so
# the objects retained from the workbook bundle are not overwritten.
allpairs_power_prereq_env = new.env(parent = globalenv())
sys.source(
  "/PHShome/kv1026/p05/iReviewPowerEffectCurve.R",
  envir = allpairs_power_prereq_env
)
allpairs_figure3_prereq_env = new.env(parent = globalenv())
sys.source(
  "/PHShome/kv1026/p05/iReviewPowerFigure3Compare.R",
  envir = allpairs_figure3_prereq_env
)

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

allpairs_attenuation_env = new.env(parent = globalenv())
sys.source(
  "/PHShome/kv1026/p05/iReviewPowerAllPairsAttenuation.R",
  envir = allpairs_attenuation_env
)
source("/PHShome/kv1026/p05/iReviewPowerAllPairsAttenuationFinalize.R")
