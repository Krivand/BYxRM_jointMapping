"""Run threshold recovery, power and figures in a new approved card directory."""
import hashlib
import json
import os
from pathlib import Path
import subprocess

source = Path(__file__).resolve().parent
output = Path(os.environ['ERIS_CARD_OUTPUT']).resolve()
root = '/PHShome/kv1026/p05'
rscript = '/data/gusev/ADQTL/envs/adqtl/bin/Rscript'
env = os.environ.copy()
env.update(PATH='/data/gusev/ADQTL/envs/adqtl/bin:/usr/bin:/bin',
           R_LIBS_USER='/data/gusev/ADQTL/envs/adqtl/lib/R/library',
           OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
scripts = ['iReviewFrankThresholdAudit.R', 'iReviewFrankPower.R', 'iReviewFrankPowerFigures.R']
report = {'source_sha256': {s: hashlib.sha256((source / s).read_bytes()).hexdigest() for s in scripts}, 'steps': []}
steps = [('syntax', [rscript, '--vanilla', '-e',
          'invisible(lapply(commandArgs(TRUE), parse))', *[str(source / s) for s in scripts]]),
         ('thresholds', [rscript, '--vanilla', str(source / scripts[0]), root, str(output / 'thresholds')]),
         ('power', [rscript, '--vanilla', str(source / scripts[1]), root, str(output / 'power'), str(output / 'thresholds')]),
         ('figures', [rscript, '--vanilla', str(source / scripts[2]), str(output / 'power')])]
exit_code = 0
try:
    for name, command in steps:
        with (output / (name + '.log')).open('x') as log:
            proc = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        report['steps'].append({'step': name, 'exit_code': proc.returncode})
        if proc.returncode:
            exit_code = proc.returncode
            break
except Exception as exc:
    report['error'] = str(exc)
    exit_code = 1
report['summaries'] = {}
for name in ['syntax.log', 'thresholds.log', 'power.log', 'figures.log',
             'thresholds/threshold_mapqtl_calls.tsv', 'thresholds/threshold_identity_evidence.tsv',
             'thresholds/threshold_full_first_pass_by_study_chr.tsv',
             'power/actual_threshold_fit_summary.tsv', 'power/figure_captions.txt']:
    path = output / name
    if path.is_file():
        data = path.read_bytes()
        report['summaries'][name] = {'text': data[:100000].decode('utf-8', errors='replace'), 'truncated': len(data) > 100000}
report['exit_code'] = exit_code
with (output / 'power_recovery_report.json').open('x') as handle:
    json.dump(report, handle, indent=2)
raise SystemExit(exit_code)
