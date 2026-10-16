"""Read bounded p05 source and input metadata for Frank's revision analyses."""
import hashlib
import json
import os
from pathlib import Path
import subprocess

root = Path("/PHShome/kv1026/p05")
out = Path(os.environ["ERIS_CARD_OUTPUT"])
if not root.is_dir():
    raise SystemExit("p05 is not accessible")
report = {"root": str(root.resolve()), "files": {}, "directories": {}, "git": {}}
for args in (["rev-parse", "HEAD"], ["status", "--short"], ["remote", "get-url", "origin"]):
    p = subprocess.run(["git", "-C", str(root)] + args, capture_output=True, text=True, timeout=15)
    report["git"][" ".join(args)] = {"returncode": p.returncode, "stdout": p.stdout[:20000], "stderr": p.stderr[:1000]}
paths = ["conventions.md", "AGENTS.md", "raw/studyMetadata.csv",
         "iReviewPowerThresholdAudit.R", "iReviewPowerAllPairs.R",
         "jReviewGenes.R", "iReviewInteractionDependence.R",
         "asys/review/power/threshold_selected_full_calls.tsv",
         "asys/review/power/threshold_full_first_pass_by_study_chr.tsv",
         "asys/review/power/threshold_full_all_pass_by_study.tsv",
         "asys/review/interaction_dependence/segregant_dependence_summary.tsv",
         "asys/review/power_effect_curve/all_pairs/pairwise_power_fit_summary.tsv"]
for rel in paths:
    p = root / rel
    if not p.is_file():
        report["files"][rel] = {"exists": False}
        continue
    size = p.stat().st_size
    if size > 200000:
        report["files"][rel] = {"size": size, "read": False}
        continue
    data = p.read_bytes()
    report["files"][rel] = {"size": size, "sha256": hashlib.sha256(data).hexdigest(), "text": data.decode("utf-8")}
# No raw matrices or individual identifiers are exported.
for rel in ("", "raw", "imed", "imed/Frank", "asys/review", "asys/review/power"):
    p = root / rel
    if p.is_dir():
        entries = sorted(p.iterdir(), key=lambda x: x.name)
        report["directories"][rel] = {"truncated": len(entries) > 200,
            "entries": [{"name": x.name, "is_dir": x.is_dir(), "size": x.stat().st_size}
                        for x in entries[:200] if not x.name.startswith(".")]}
for p in (root.parent / "conventions.md", root.parent / "AGENTS.md"):
    if p.is_file() and p.stat().st_size <= 50000:
        report["files"][str(p)] = {"text": p.read_text()}
rscript = Path("/data/gusev/ADQTL/envs/adqtl/bin/Rscript")
report["rscript_exists"] = rscript.is_file()
if rscript.is_file():
    p = subprocess.run([str(rscript), "--vanilla", "-e",
        'cat(R.version.string, "\\n"); for (p in c("readxl", "GenomicRanges", "IRanges", "ggplot2")) cat(p, requireNamespace(p, quietly = TRUE), "\\n")'],
        capture_output=True, text=True, timeout=90)
    report["r_runtime"] = {"returncode": p.returncode, "stdout": p.stdout, "stderr": p.stderr[-4000:]}
with (out / "p05_inspection.json").open("x") as handle:
    json.dump(report, handle, indent=2)
print(json.dumps({"p05_accessible": True, "output": "p05_inspection.json", "rscript_exists": report["rscript_exists"]}))
