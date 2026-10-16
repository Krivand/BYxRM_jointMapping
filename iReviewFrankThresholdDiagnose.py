"""Collect bounded, read-only evidence for unresolved historical threshold IDs.

Run only through an approved ERIS card. No mapping, fitting, source changes,
network calls or job submissions. Outputs are created exclusively in card output.
"""
import hashlib
import json
import os
from pathlib import Path

root = Path('/PHShome/kv1026/p05')
previous = Path('/PHShome_actual/k/kv1026/tools/eris-jobcards/actions/card_36fb24393d12769900830e35e4625315/redux/thresholds')
output = Path(os.environ['ERIS_CARD_OUTPUT'])
sources = {
    'mapping_log': root / 'logs/iRemapAll.out',
    'mapping_driver': root / 'iRemapAll.R',
    'mapping_functions': root / 'remap_functions.R',
    'metadata': root / 'raw/studyMetadata.csv',
    'calls': previous / 'threshold_mapqtl_calls.tsv',
    'records': previous / 'threshold_records_all.tsv',
    'context': previous / 'threshold_context.tsv',
}
report = {'purpose': 'Resolve threshold study identities without guessing call order', 'files': {}}
for name, path in sources.items():
    try:
        with path.open('rb') as handle:
            data = handle.read(8 * 1024 * 1024 + 1)
        if len(data) > 8 * 1024 * 1024:
            report['files'][name] = {'path': str(path), 'error': 'exceeds 8 MiB bound'}
            continue
        report['files'][name] = {
            'path': str(path), 'bytes': len(data),
            'sha256': hashlib.sha256(data).hexdigest(),
            'text': data.decode('utf-8', errors='replace'),
        }
    except OSError as exc:
        report['files'][name] = {'path': str(path), 'error': str(exc)}
with (output / 'threshold_diagnosis.json').open('x') as handle:
    json.dump(report, handle, indent=2)
print(json.dumps({name: {k: v for k, v in item.items() if k != 'text'}
                  for name, item in report['files'].items()}))
raise SystemExit(int(any('error' in item for item in report['files'].values())))
