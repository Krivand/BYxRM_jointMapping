# Joint mapping of mRNA and protein abundance in the BY × RM yeast cross

This repository accompanies [High concordance between genetic effects on mRNA and protein abundance](https://doi.org/10.64898/2026.05.26.727960). It combines published expression datasets from the BY × RM cross in *Saccharomyces cerevisiae* to estimate how often QTLs affect mRNA and protein in the same direction.

## Inputs and environment

The manuscript reports R 4.5.2 for most analyses. The project scripts use the ERIS root `/PHShome/kv1026/p05`. Running them elsewhere requires adapting absolute paths in the entry points, the shared `remap_functions.R`, and any selected runners. The shared file sets the working directory, selects the R library, and loads or installs packages including rrBLUP, lme4, lmerTest, mashr and tidyverse. The GO analysis additionally uses topGO and the SGD annotation database. Some review workflows use Bash or Python 3.

The manuscript's Data sources section identifies the original studies and accessions. The scripts also expect project-prepared files outside Git, including `raw/studyMetadata.csv`, `raw/chrFeatures.rda`, study-specific expression and genotype inputs, and saved study objects for downstream analyses. The clone alone does not supply these inputs. Prepared objects and chunk results go under `imed/`, figures and analysis tables under `asys/`, and mapping logs under `logs/`.

## Analysis entry points

| Analysis | How to run it |
| --- | --- |
| Individual studies | `iRemapAll.R` prepares the study objects. Its final commented block records the calls to `findBackgroundQTLs` and `fullMap` for background, full and simple maps. Running the file as supplied does not execute that block. |
| Joint mapping | With the study objects and background QTLs available, run `Rscript jBH_analyte_chunk.R 1 50` for the first 50-gene chunk of the genotype scan. Add `analyte` as the third argument for the interaction scan. Repeat for each chunk. `jBH_medium_chunk.R` handles the medium comparison. `jAnalysis.R` reads the saved chunks from both comparisons. |
| Cross-study effect estimates | `jMashrTrain.R` builds the mash fits under `imed/mashr/`. With the required arm-specific inputs in place, use `jMashrCov.R` for covariance summaries, `jMashrShare.R` for effect sharing and `jMashrCross.R` for cross-study comparisons. |
| Review analyses | Run `bash iReviewFrankRun.sh /PHShome/kv1026/p05 /path/to/new/output` for supported-gene regulation, segregant dependence, threshold recovery and power. The runner uses the ERIS R installation and records completion for each analysis. See the [GO methods](GO_REVIEW_METHODS.md) and [rediscovery export and figure instructions](REDISCOVERY_EXPORT.md) for their separate entry points and required inputs. |

## Reproducing the review comparisons

- The supported-gene analysis uses the Albert 2018 QTL table in `imed/Frank/qtls_list_260204.rda` and the tested genes in `imed/Albert_2018.rda`. The project table's identity with the original published table was not established in the inspection.
- The dependence analysis targets the union of 868 joint and 64 strict interaction gene–marker pairs. The model crosses study and segregant random intercepts. Treating normalized segregant labels as the same individual across studies is an assumption of this fit.
- The power analysis recovers 160 first-pass FULL-map study–chromosome thresholds from the historical log. Relating them to Figure 3 assumes the log corresponds to its QTL tables. The calculation estimates detection at the first-pass cutoff and does not replay conditional mapping or localization. Compare thresholds within the same query set. Scaling source effects to 80% is a sensitivity analysis.
- `iReviewPower.R` contains the simulation, and `iReviewFrankSimulationFigure.R` produces its calibration figure. The rediscovery export uses predictions and observed counts from the Figure 3 matched query set.
