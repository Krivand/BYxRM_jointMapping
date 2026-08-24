#!/usr/bin/env python3
from pathlib import Path

path = Path('iReviewPower.R')
text = path.read_text()

start_marker = '''############################\n# EMPIRICAL LOD CALIBRATION #\n############################\n'''
end_marker = '''################################\n# EMPIRICAL EFFECT CALIBRATION #\n################################\n'''
if start_marker not in text or end_marker not in text:
    raise SystemExit('Expected LOD calibration markers were not found; refusing to patch.')

start = text.index(start_marker)
end = text.index(end_marker, start)
new_block = '''#############################
# LOD THRESHOLD SENSITIVITY #
#############################

# Primary inference is made across a common LOD sweep, not from a pooled or
# study-specific average of historical permutation thresholds. In each sweep
# run every study receives the same candidate LOD, so differences arise from
# study sample size and query conditioning. The recovered FULL first-pass
# thresholds span LOD 2.42-6.35; 2-7 deliberately brackets that range.
lod_sweep = seq(2, 7, by = 0.25)
threshold_sensitivity_lod = lod_sweep
threshold_sensitivity_label = sprintf("common_lod_%0.2f", lod_sweep)

# LOD 4.5 is retained only as a display/diagnostic anchor for existing plots
# and validation objects. It is not the primary inferential threshold.
lod_reference = 4.5
empirical_thresholds = numeric(0)
threshold_summary = data.frame(
  statistic = rep("common_sensitivity_sweep", length(lod_sweep)),
  lod = lod_sweep,
  stringsAsFactors = FALSE
)

message(
  "Primary LOD sensitivity sweep: ",
  min(lod_sweep), " to ", max(lod_sweep),
  " by ", unique(diff(lod_sweep)),
  "; no historical threshold averaging."
)
message("LOD 4.5 is a display/diagnostic anchor only.")

study_design$display_lod = lod_reference
study_design$display_rcrit = lod_to_rcrit(
  lod = lod_reference,
  n = study_design$n
)

'''
text = text[:start] + new_block + text[end:]

old_validation = '''validation_grid = expand.grid(
  n = validation_n,
  rho = validation_rho,
  stringsAsFactors = FALSE
)
'''
new_validation = '''validation_grid = expand.grid(
  n = validation_n,
  rho = validation_rho,
  lod = c(min(lod_sweep), lod_reference, max(lod_sweep)),
  stringsAsFactors = FALSE
)
'''
if old_validation not in text:
    raise SystemExit('Validation grid block not found; refusing to patch.')
text = text.replace(old_validation, new_validation, 1)

old_call = '''        n = validation_grid$n[[i]],
        rho = validation_grid$rho[[i]],
        lod = lod_reference,
'''
new_call = '''        n = validation_grid$n[[i]],
        rho = validation_grid$rho[[i]],
        lod = validation_grid$lod[[i]],
'''
if old_call not in text:
    raise SystemExit('Validation call block not found; refusing to patch.')
text = text.replace(old_call, new_call, 1)

anchor = '''observed_threshold_sensitivity_summary =
  observed_threshold_sensitivity %>%
  dplyr::filter(
    is.finite(observed_rediscovery),
    is.finite(expected_rediscovery)
  ) %>%
  dplyr::group_by(threshold, lod) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    mean_observed = mean(observed_rediscovery),
    mean_expected = mean(expected_rediscovery),
    mean_difference = mean(difference),
    .groups = "drop"
  )
'''
if anchor not in text:
    raise SystemExit('Threshold-sensitivity summary block not found; refusing to patch.')
addition = anchor + '''
# Primary threshold-robust result: for each ordered study pair, summarize the
# full expected rediscovery range over LOD 2-7 rather than privileging one LOD.
observed_threshold_envelope = observed_threshold_sensitivity %>%
  dplyr::group_by(query, target) %>%
  dplyr::summarise(
    observed_rediscovery = dplyr::first(observed_rediscovery),
    expected_minimum = min(expected_rediscovery, na.rm = TRUE),
    expected_maximum = max(expected_rediscovery, na.rm = TRUE),
    lod_at_expected_minimum = lod[which.min(expected_rediscovery)][1],
    lod_at_expected_maximum = lod[which.max(expected_rediscovery)][1],
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    relative_to_expected_range = dplyr::case_when(
      observed_rediscovery < expected_minimum ~ "below",
      observed_rediscovery > expected_maximum ~ "above",
      TRUE ~ "within"
    )
  )
'''
text = text.replace(anchor, addition, 1)

write_anchor = '''write_table(
  observed_threshold_sensitivity_summary,
  "observed_threshold_sensitivity_summary.tsv"
)
'''
if write_anchor not in text:
    raise SystemExit('Threshold-sensitivity write block not found; refusing to patch.')
text = text.replace(
    write_anchor,
    write_anchor + '''write_table(
  observed_threshold_envelope,
  "observed_threshold_envelope.tsv"
)
''',
    1
)

list_anchor = '''  observed_threshold_sensitivity_summary =
    observed_threshold_sensitivity_summary,
'''
if list_anchor not in text:
    raise SystemExit('review_power threshold list block not found; refusing to patch.')
text = text.replace(
    list_anchor,
    list_anchor + '''  observed_threshold_envelope = observed_threshold_envelope,
''',
    1
)

text = text.replace(
    '"pooled median historical 5% FDR LOD calibration = ",',
    '"display/diagnostic LOD anchor only = ",',
    1
)
text = text.replace(
    '"threshold_reference",',
    '"threshold_display_anchor",',
    1
)
text = text.replace(
    '"threshold_sensitivity",',
    '"threshold_primary_sensitivity",',
    1
)

old_final = '''message(
  "Reference LOD: ",
  signif(lod_reference, 4)
)
'''
new_final = '''message(
  "Primary LOD sweep: ",
  min(lod_sweep), " to ", max(lod_sweep),
  " by ", unique(diff(lod_sweep))
)
message("Historical LOD thresholds averaged in primary model: no")
message("LOD 4.5 retained only as display/diagnostic anchor")
'''
if old_final not in text:
    raise SystemExit('Final reference-LOD message not found; refusing to patch.')
text = text.replace(old_final, new_final, 1)

path.write_text(text)
print('Patched iReviewPower.R for common LOD 2-7 sensitivity analysis.')
print('Primary output: asys/review/power/observed_threshold_envelope.tsv')
