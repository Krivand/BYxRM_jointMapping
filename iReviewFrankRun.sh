#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "Usage: bash iReviewFrankRun.sh P05_ROOT NEW_OUTPUT_DIRECTORY" >&2
  exit 2
fi
project_root="$1"
review_output="$2"
source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rscript="/data/gusev/ADQTL/envs/adqtl/bin/Rscript"
export PATH="/data/gusev/ADQTL/envs/adqtl/bin:/usr/bin:/bin"
export R_LIBS_USER="/data/gusev/ADQTL/envs/adqtl/lib/R/library"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
[[ -d "$project_root" && -x "$rscript" ]]
if [[ -e "$review_output" ]]; then
  echo "Output exists; refusing overwrite: $review_output" >&2
  exit 2
fi
# Parse every script before creating outputs or executing analyses.
for script in jReviewGeneRegulation.R iReviewFrankDependence.R iReviewFrankThresholdAudit.R iReviewFrankPower.R; do
  "$rscript" --vanilla -e 'invisible(parse(commandArgs(TRUE)[[1]]))' "$source_dir/$script"
done
mkdir "$review_output"
failed=0
run_analysis() {
  local script="$1"
  local label="$2"
  shift 2
  if "$rscript" --vanilla "$source_dir/$script" "$project_root" "$review_output/$label" "$@" >"$review_output/$label.log" 2>&1; then
    printf '%s\tcompleted\n' "$label" >>"$review_output/status.tsv"
  else
    printf '%s\tfailed\n' "$label" >>"$review_output/status.tsv"
    failed=1
    return 1
  fi
}
printf 'analysis\tstatus\n' >"$review_output/status.tsv"
run_analysis jReviewGeneRegulation.R gene_regulation || true
run_analysis iReviewFrankDependence.R dependence || true
if run_analysis iReviewFrankThresholdAudit.R thresholds; then
  run_analysis iReviewFrankPower.R power "$review_output/thresholds" || true
else
  printf 'power\tblocked_missing_verified_thresholds\n' >>"$review_output/status.tsv"
fi
cat "$review_output/status.tsv"
exit "$failed"
