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
It preserves exactly the original ten columns, their names and their order.
Both prediction columns already require matching effect direction; no duplicate
columns or definition fields are added. All probabilities are fractions, not
percentages. Observed counts are
recovered from the matched denominator and observed fraction, with integer checks.
The provenance JSON records input/output hashes. Direction must be verified from
the generating run: the aggregated input alone cannot establish it.

Use `--reference previous.tsv` to require all original columns and all 90 keys
to agree with a prior export. Omit this option when intentionally exporting
changed upstream results; review those changes separately.

## Recovered original export

All ten columns of the archived comparison table agree across all 90 directed
study pairs with the export from the Figure 3 matched rows. Both prediction
columns match at full floating-point precision after parsing. The archived
generating source implements same-direction detection (SHA-256
`9bb8781df99341b4737c4db7e3c10e138c3ffab912f4067b142950e8404be7f5`).

## Matched manuscript figures

Both manuscript scatterplots now use `figure3_matched` predictions, like the
comparison table and the paired-panel figures. The plotting script verifies that
prediction denominators, observed denominators and point-size counts agree.
It preserves the original ggplot layouts and writes both PNG and PDF files.

```bash
Rscript --vanilla iReviewFrankOriginalFigures.R /path/to/output /path/to/old/all_pairs /path/to/power/pairwise_actual_threshold_power.tsv
```

The optional input arguments remove dependence on archived absolute paths;
the original one-argument invocation remains supported. The old all-pairs file
supplies comparison labels and observed counts, not prediction values.
No mapping or simulation rerun is required. The original ten-column TSV is
unchanged. All-90 rounded summaries remain 65.0% / 46.8% predicted and 26.3%
observed; the Albert full-effect mean is 6.05% (range 1.02–8.92%) and its
80%-effect mean is 3.37%.
