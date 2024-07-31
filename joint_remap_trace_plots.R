############################################
# SOURCE IMPORTANT FUNCTIONS AND CONSTANTS #
############################################

source(paste0("/home/albertf/vandy147/p05/",
              "remap_functions.R"))

#########################
# LOAD IN CRITICAL DATA #
#########################

# LOAD IN UNFILTERED ASSOCATIONS AND PREPARE
load(paste0(imed,"unfilteredAssocs.rda"))

# Remove duplicate marker and gene columns
unfilteredAssocs = unfilteredAssocs %>% 
  rename(marker = geno.marker, gene = geno.gene) %>%
  select(-matches("\\.marker$"), 
         -matches("\\.gene$"), 
         -matches("\\.flavor$")) %>%
  mutate(pairID = paste(gene, marker, sep = "_")) %>% 
  select(pairID, everything())

#########################
# DETECT CORES AND PLOT #
#########################

# detect cores
nCores = detectCores()

# change the pval name to something we actually want to plot
uaPlot = unfilteredAssocs %>% select(-starts_with("geno")) %>% 
  mutate(across(ends_with("pval"), ~ -log10(.),
                .names = "{sub('pval$', '', .col)}-log10(p-value)"))

# plot the -log10(p) traceplots
tracePlotter(dat = uaPlot, 
             column = '-log10(p-value)',
             rootDir = paste0(asys,'ThreeFlavorTracePlotsLogP/'),
             rootName = "ThreeTraceLogP_",
             mc.cores = nCores,
             r = c("geno", "ixn"), 
             n = c("genotype", "interaction"))

# plot the effect-size traceplots
tracePlotter(dat = uaPlot, 
             column = 'effect',
             rootDir = paste0(asys,'ThreeFlavorTracePlotsEffect/'),
             rootName = "ThreeTraceEffect_",
             mc.cores = nCores,
             r = c("geno", "ixn"), 
             n = c("genotype", "interaction"))