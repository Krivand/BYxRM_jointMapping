"""Run the corrected simulation figure and export its PNG and PDF in the same card."""
import base64, hashlib, json, os, subprocess, sys
from pathlib import Path
source = Path(__file__).resolve().parent
output = Path(os.environ["ERIS_CARD_OUTPUT"]).resolve()
env = os.environ.copy()
env.update(PATH="/data/gusev/ADQTL/envs/adqtl/bin:/usr/bin:/bin",
           R_LIBS_USER="/data/gusev/ADQTL/envs/adqtl/lib/R/library",
           OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
with (output / "simulation.log").open("x") as log:
    p = subprocess.run(["/data/gusev/ADQTL/envs/adqtl/bin/Rscript", "--vanilla",
        str(source / "iReviewFrankSimulationFigure.R"), str(output), *sys.argv[1:]],
        env=env, stdout=log, stderr=subprocess.STDOUT)
with (output / "simulation_report.json").open("x") as f:
    json.dump({"exit_code": p.returncode, "log": (output / "simulation.log").read_text()}, f)
if p.returncode:
    raise SystemExit(p.returncode)
exports = {}
for extension, signature in (("png", b"\x89PNG\r\n\x1a\n"), ("pdf", b"%PDF-")):
    name = "simulation_actual_thresholds." + extension
    data = (output / name).read_bytes()
    if not data.startswith(signature) or len(data) >= 2097152:
        raise ValueError("Invalid or oversized figure: " + name)
    exports[name] = {"sha256": hashlib.sha256(data).hexdigest(),
                     "size_bytes": len(data), "base64": base64.b64encode(data).decode()}
with (output / "plots_export.json").open("x") as f:
    json.dump(exports, f)
