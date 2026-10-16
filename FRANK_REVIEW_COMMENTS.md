# Frank reviewer comments 15 October 2026

Based on the comments in KJV_reviewresponse and revision source at
3bc11caa95e61c1a48e9e8e86f674f57e446ed43. This branch stages work;
it does not claim scientific results.

No conventions.md or AGENTS.md appears in either inspected GitHub tree.
The p05 inspection also checks the project and immediate home-level copies
before remote analysis execution.

## Comment coverage

| Comment | Requested action | Implementation and remaining evidence |
| --- | --- | --- |
| 0 | Compare local and trans regulation of the 385 genes | jReviewGeneRegulation.R: Fisher tests for presence of cis and trans eQTLs, odds ratios and confidence intervals; supplementary per-gene count tests. Requires a verified published AB18 QTL annotation adapter and full tested-gene list. Never infer the background from significant QTL genes. |
| 1 | Define targeted | Existing iReviewInteractionDependence.R refits the unique union of joint-QTL and strict interaction-QTL gene-marker pairs, not a new genome-wide scan. |
| 2 | Segregant or panel random effect? | Existing code adds (1 \| segregant), crossed with (1 \| study). Segregant labels come from genotype row names after normalization. This is an individual-segregant intercept, not a panel intercept; label equivalence across datasets remains a scientific assumption to verify. |
| 3 | Why 928 fits? | Union of 868 joint and 64 strict gene-marker pairs. If the saved union is 928, the intersection is four exact pairs. iReviewFrankDependence.R verifies this from stored fit rows and summarizes both random-intercept variances. Do not infer which component causes singularity without the variance evidence. |
| 4 | Explain simulation plot and axes | Existing iReviewPower.R uses balanced two-allele genotypes and Gaussian phenotype noise, recomputes correlations and LODs, and compares empirical detection frequencies to the noncentral-t calculation. The validation table separates analytic_power, simulated_power and same-direction quantities. Verify the actual embedded plot before specifying which quantity is on each axis. |
| 5 | Use actual LOD for Albert projections | iReviewFrankPower.R includes the nine Albert-to-target projections using recovered FULL-map first-pass target-study/chromosome thresholds. |
| 6 | Use actual LOD for all 90 pairs | Same script outputs all 90 directed pairs, actual-threshold versus fixed-4.5 predictions, original and 0.8-scaled effects, correlation and absolute-error summaries. |

## Interpretation and execution

The threshold audit already exists on the parent revision branch.
The new power script requires all 160 first-pass study/chromosome thresholds,
completed selected FULL calls, and exact agreement with their saved log lines.
Missing or ambiguous coverage stops execution; there is no pooled-threshold fallback.

Before interpretation, verify that these mapping-log calls correspond to the QTL
tables used by Figure 3. This provenance check remains open. First-pass thresholds
are actual recovered cutoffs, but a fixed-effect power calculation does not replay
later conditional mapping passes or localization. Describe it accordingly.

Both historical gene-overlap and Figure 3 matched query sets are retained.
Compare actual and fixed thresholds within each denominator, so changes in
technical-marker exclusions are not mistaken for changes caused by thresholds.
The 0.8 effect factor remains a sensitivity analysis, not a fitted correction.

All three analysis scripts require a new output directory and refuse overwrite.
They do not source remap_functions.R, install packages, delete outputs or rerun
mapping. No manuscript edits have been made.

## Current execution record

- Existing bigi allocation: 5904729, four CPUs and 16 GiB, scheduled end
  2026-10-15 22:47:11 ERIS scheduler time.
- Inspection card card_a0a531ad3c6c02c5c60d07ea71a701f4 ran as step
  5904729.4 and failed before report creation because git was absent from PATH.
  The script had already verified that /PHShome/kv1026/p05 is a directory.
  No scientific analysis ran.
- iReviewFrankInspectV2.py handles unavailable Git and uses an explicit R PATH.
  It adds bounded AB18 schema inspection. Its new card requires separate approval.
- Python inspection source parses locally. R is unavailable in the local runtime;
  R parse checks, real-data execution and numerical validation remain pending.
- The AB18 input adapter will be completed after the inspection identifies the
  canonical local/trans annotation source and tested-gene universe.
