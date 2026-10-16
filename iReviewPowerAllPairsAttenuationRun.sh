#!/usr/bin/env bash

set -euo pipefail

ROOT="/PHShome/kv1026/p05"
cd "$ROOT"

echo "=== CONVENTIONS ==="
CONVENTIONS="$(find .. -name conventions.md -print -quit)"
if [[ -n "$CONVENTIONS" ]]; then
    echo "Found: $CONVENTIONS"
    cat "$CONVENTIONS"
else
    echo "No conventions.md found"
fi

echo
echo "=== INPUT CHECKS ==="
test -s asys/review/power_effect_curve/all_pairs/pairwise_observed_vs_power.tsv
test -s asys/review/final/pairwise_power_predictions.tsv
test -s asys/review/final/power_analysis.xlsx
test -s asys/review/review_output_manifest.tsv
echo "Required existing all-pair outputs present"

export REVIEW_POWER_CORES="${REVIEW_POWER_CORES:-8}"
echo
echo "=== ATTENUATION SENSITIVITY ==="
echo "REVIEW_POWER_CORES=$REVIEW_POWER_CORES"
Rscript iReviewPowerAllPairsAttenuation.R

echo
echo "=== ATTENUATION PLOTS ==="
Rscript iReviewPowerAllPairsAttenuationPlots.R

echo
echo "=== FINAL SUPPORT UPDATE ==="
Rscript iReviewPowerAllPairsAttenuationFinalize.R

echo
echo "=== PLOT FILES ==="
ls -lh \
  asys/review/power_effect_curve/all_pairs/all_pairs_observed_vs_power_lod4_5_scale0_8.png \
  asys/review/power_effect_curve/all_pairs/all_pairs_observed_vs_power_lod4_5_scale0_8.pdf \
  asys/review/power_effect_curve/all_pairs/all_pairs_observed_vs_power_lod4_5_scale0_6.png \
  asys/review/power_effect_curve/all_pairs/all_pairs_observed_vs_power_lod4_5_scale0_6.pdf

echo
echo "=== SUMMARY ==="
cat asys/review/power_effect_curve/all_pairs/power_attenuation_summary.tsv

echo
echo "=== GIT STATUS ==="
git status --short

echo
echo "ATTENUATION PLOT RUN COMPLETE"
