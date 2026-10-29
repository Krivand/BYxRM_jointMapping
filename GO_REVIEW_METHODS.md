# Interaction GO analysis

- Background: 385 supported genes; five sets of 64, 37, 14, 9 and 4 genes.
- Annotations: org.Sc.sgd.db via AnnotationDbi; topGO propagates membership.
- Primary tests: topGO classic one-sided Fisher, separately for BP/CC/MF.
- Eligibility: at least three annotated background genes per term (nodeSize=3).
- Primary correction: BH over all tested terms within each set x ontology.
- Full-precision topGO scores enter BH; displayed/rounded P values do not.
- Follow-up: 37 discordant versus 27 other interaction genes, one-sided Fisher.
- Follow-up eligibility: the same complete term list as the primary discordant
  analysis, without filtering by primary significance or membership in the 64.
- Follow-up correction: BH separately within BP/CC/MF over the complete family.
- Empty/full membership among all 64: retain P=1 in correction; OR is NA and
  estimable=FALSE. An informative zero cell may yield Inf, displayed as text.
- Summary: every primary discordant hit (BH P<0.05), with its follow-up P/Q and
  a significant_both flag. This is descriptive overlap, not an independently
  validated discovery set or a claim of a separate FDR guarantee for overlap.
- Full follow-up results remain available, including follow-up-only hits.
- MKT1/BUL2: descriptive class counts and effect directions only. The bespoke
  six-term ribosome union and targeted hotspot enrichment tests are removed.

## GO-only execution

Run `Rscript --vanilla iReviewGOStandalone.R /absolute/new/output/directory`
with jReviewGO.R and iReviewWorkbookFinalize.R beside the entry point.
The destination must not already exist. Existing mapping inputs and annotation
packages are read; no power simulation, remapping or final-directory replacement
is performed. Execution still requires the normal reviewed ERIS card.

Outputs: interaction_GO.xlsx (10 sheets), full primary/follow-up TSVs, summary
TSVs, GO_results.rds (annotations, gene sets and results), and sessionInfo.txt.
The workbook Summary includes all first-pass discordant hits. Discordant vs other
replaces Ribosome comparison. The five primary sheets, Significant terms,
Gene sets and descriptive MKT1 BUL2 summary remain.

The entry point validates full term coverage, 37/27 denominators, per-ontology
BH, Fisher against an independent hypergeometric tail, zero-cell boundaries,
and workbook round-trip values. The previous Inf export fix passed on ERIS with
R 4.4.3/openxlsx 4.2.8.1 (card_cf831dab8499d99996c9d6132395f8fd).
The full-family GO run completed on 2026-10-29 using source commit
f09751a7af9aad2492312c0d8c1af2e508356822 (Slurm step 6268410.1;
card_64c19d3c9a2f864a6d2f3211894785fa). The analysis process exited 0 and
reported that all statistical and workbook checks passed. Input hashes were
unchanged. Results comprise 1,279 follow-up terms (842 BP, 176 CC, 261 MF),
29 primary discordant hits (12 BP, 13 CC, 4 MF), and 11 significant follow-ups
(2 BP, 7 CC, 2 MF), all also primary hits.

Outputs are in /PHShome/kv1026/p05/asys/review/go_full_family_f09751a_20261029.
The interaction_GO.xlsx SHA-256 is
d7ab130bae22a5caea82a0956cf82eea9ec67d206cdd820bba2bb22aa5bb0549.
The Observer card remains FAILED: after successful analysis, its output
collector raised FileNotFoundError for 'p05' because output paths were declared
home-relative rather than card-relative. This collection failure does not
indicate an analysis or validation failure; outputs and logs were inspected.
