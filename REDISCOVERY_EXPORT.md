# Rediscovery comparison export

Generate the upstream predictions with the committed `iReviewFrankPower.R`.
It computes same-direction detection using the upper noncentral-t tail,
with recovered target-study/chromosome first-pass LOD thresholds. Its
`pairwise_actual_threshold_power.tsv` contains two denominators and two
effect scales. Do not substitute the older fixed-LOD pairwise files.

Export the Figure 3 matched denominator (including the technical-marker
exclusions) at full and 80% source-effect magnitude:

```bash
python3 iReviewRediscoveryExport.py /path/to/power/pairwise_actual_threshold_power.tsv /path/to/new/export
```

The destination must be new. The exporter does not rerun mapping or simulations.
It preserves the original ten columns and appends explicitly named same-direction
prediction columns, a direction flag, and the denominator name. The explicit
prediction columns are aliases of the original prediction columns, not a second
analysis. All probabilities are fractions, not percentages. Observed counts are
recovered from the matched denominator and observed fraction, with integer checks.
The provenance JSON records input/output hashes. Direction must be verified from
the generating run: the aggregated input alone cannot establish it.

Use `--reference previous.tsv` to require all original columns and all 90 keys
to agree with a prior export. Omit this option when intentionally exporting
changed upstream results; review those changes separately.

## Recovered original export

The user-supplied `rediscovery_predictions_comparisons (1).tsv` contains all
90 comparisons and matches the `figure3_matched` rows from
`card_1e080f28d4b1c883cb4f66157701e1a7/power/pairwise_actual_threshold_power.tsv`.
Both prediction columns match exactly as parsed floating-point numbers.
All ten original columns pass the exporter reference comparison.
The archived generating source in
`bundle_b303a517e324d839795eb37bbe370f4cf3ce2d7575b8d477a729e59d5f812b10`
already contains `same_direction_detection` (source SHA-256
`9bb8781df99341b4737c4db7e3c10e138c3ffab912f4067b142950e8404be7f5`).
Thus the existing table already incorporates direction; the new columns make
that definition explicit. The exact original ad hoc export script was not found
among the committed Frank review/export scripts; this script records the mapping.

The response-document scatterplots use `historical_gene_overlap` predictions,
whereas this table uses `figure3_matched`. Their predictions can differ for
comparisons affected by technical-marker exclusions. Do not silently substitute
one denominator for the other or describe their values as identical.
