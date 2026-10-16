"""Render original manuscript plots and package lossless PNG downloads."""
import base64, hashlib, json, os, subprocess
from pathlib import Path
source = Path(__file__).resolve().parent
output = Path(os.environ["ERIS_CARD_OUTPUT"]).resolve()
env = os.environ.copy()
env.update(PATH="/data/gusev/ADQTL/envs/adqtl/bin:/usr/bin:/bin",
           R_LIBS_USER="/data/gusev/ADQTL/envs/adqtl/lib/R/library",
           OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
with (output / "figures.log").open("x") as log:
    p = subprocess.run(["/data/gusev/ADQTL/envs/adqtl/bin/Rscript", "--vanilla",
        str(source / "iReviewFrankOriginalFigures.R"), str(output)],
        env=env, stdout=log, stderr=subprocess.STDOUT)
report = {"exit_code": p.returncode, "log": (output / "figures.log").read_text()}
with (output / "figure_report.json").open("x") as f:
    json.dump(report, f)
if p.returncode:
    raise SystemExit(p.returncode)
exports = {}
for name in ["all_pairs_actual_thresholds.png", "all_pairs_actual_thresholds_scale0_8.png"]:
    data = (output / name).read_bytes()
    assert data[:8] == bytes([137,80,78,71,13,10,26,10]) and len(data) < 2097152
    exports[name] = {"sha256": hashlib.sha256(data).hexdigest(),
                     "size_bytes": len(data), "base64": base64.b64encode(data).decode()}
with (output / "plots_export.json").open("x") as f:
    json.dump(exports, f)
