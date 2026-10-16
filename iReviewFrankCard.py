"""Run the staged review analyses and retain bounded summaries even on failure."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

source = Path(__file__).resolve().parent
output = Path(os.environ["ERIS_CARD_OUTPUT"]).resolve()
root = Path("/PHShome/kv1026/p05")
result = {"started": time.time(), "project": str(root), "analyses": {}, "source_sha256": {}}
for name in ("iReviewFrankRun.sh", "iReviewFrankPower.R", "iReviewFrankDependence.R",
             "iReviewFrankThresholdAudit.R", "jReviewGeneRegulation.R"):
    result["source_sha256"][name] = hashlib.sha256((source / name).read_bytes()).hexdigest()
command = ["/bin/bash", str(source / "iReviewFrankRun.sh"), str(root), str(output / "redux")]
env = os.environ.copy()
env["PATH"] = "/data/gusev/ADQTL/envs/adqtl/bin:/usr/bin:/bin"
try:
    with (output / "redux_execution.log").open("x") as log:
        proc = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
    result["exit_code"] = proc.returncode
except OSError as exc:
    result["exit_code"] = 1
    result["error"] = str(exc)
names = ["status.tsv", "gene_regulation.log", "dependence.log", "thresholds.log", "power.log",
         "gene_regulation/supported_gene_cis_trans_presence.tsv",
         "gene_regulation/supported_gene_cis_trans_counts.tsv",
         "gene_regulation/regulation_denominators.tsv",
         "dependence/targeted_fit_counts.tsv",
         "dependence/random_intercept_variances.tsv",
         "power/actual_threshold_fit_summary.tsv",
         "power/pairwise_actual_threshold_power.tsv",
         "thresholds/threshold_full_first_pass_by_study_chr.tsv",
         "thresholds/threshold_selected_full_calls.tsv"]
for name in names:
    path = output / "redux" / name
    if path.is_file():
        with path.open("rb") as handle:
            data = handle.read(64001)
        result["analyses"][name] = {"text": data[:64000].decode("utf-8", errors="replace"),
                                    "truncated": len(data) > 64000}
    else:
        result["analyses"][name] = {"present": False}
log = output / "redux_execution.log"
if log.is_file():
    with log.open("rb") as handle:
        result["execution_log"] = handle.read(16000).decode("utf-8", errors="replace")
result["finished"] = time.time()
with (output / "redux_report.json").open("x") as handle:
    json.dump(result, handle, indent=2)
print(json.dumps({"exit_code": result["exit_code"], "report": "redux_report.json"}))
raise SystemExit(result["exit_code"])
