#!/usr/bin/env python3
"""Re-scan identical raw inputs with the app and a separately validated reference.

Development/release check, not a self-baseline test. No outputs from the app are
accepted as expected answers. Requires ROOT/PyROOT and NumPy. Leaves only a JSON
receipt; raw inputs are read-only and temporary scans/builds are removed.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
import ROOT

HERE = Path(__file__).resolve().parents[1]
KEYS = ('detector_order rate_asymmetry_sums rate_asymmetry_mc_covariance '
        'reference_asymmetry_ppb rates F F_mc_sigma H B C correlation M '
        'experimental_covariance mc_covariance combined_covariance').split()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def array(obj):
    if not obj:
        raise ValueError('Missing required ROOT result')
    if obj.InheritsFrom('TMatrixTBase<double>'):
        return np.array([[obj[i][j] for j in range(obj.GetNcols())]
                         for i in range(obj.GetNrows())])
    return np.array([obj[i] for i in range(obj.GetNrows())])


def compare(actual, expected, key):
    if actual.shape != expected.shape:
        return {'passed': False, 'reason': 'shape mismatch'}
    if not (np.isfinite(actual).all() and np.isfinite(expected).all()):
        return {'passed': False, 'reason': 'nonfinite result'}
    if key == 'detector_order':
        return {'passed': bool(np.array_equal(actual, expected))}
    # Normalize covariance errors by the associated diagonal standard deviations,
    # not by the single largest entry (which could hide a small channel's error).
    if key == 'rate_asymmetry_mc_covariance':
        sd = np.sqrt(np.maximum(np.diag(expected), 0))
        scale = np.outer(sd, sd)
        tolerance = 1e-8 * scale + 1e-13 * np.max(scale, initial=0)
    else:
        tolerance = 1e-8 * np.abs(expected) + 1e-13 * np.max(np.abs(expected), initial=0)
    error = np.abs(actual - expected)
    failed = error > tolerance
    first = np.argwhere(failed)
    result = {'passed': not bool(failed.any()), 'mismatched_entries': int(failed.sum()),
              'max_absolute_difference': float(error.max(initial=0))}
    if len(first):
        idx = tuple(first[0])
        result.update(first_index=list(map(int, idx)), actual=float(actual[idx]), expected=float(expected[idx]))
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference-macro', type=Path, required=True)
    ap.add_argument('--reference-validation', type=Path, required=True,
                    help='Passed validation.json binding the independent macro SHA256')
    ap.add_argument('--input-receipt', type=Path, required=True,
                    help='JSON inputs: [{path, sha256, events}]; must cover all three generators')
    ap.add_argument('--reference-main-tables', type=Path, required=True)
    ap.add_argument('--reference-sm-tables', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    args.output = args.output.resolve()
    if args.output.exists():
        ap.error('Receipt already exists; preserve previous evidence and choose another path.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    reference = args.reference_macro.resolve()
    reference_tables = {'main_pe': args.reference_main_tables.resolve(),
                        'response': args.reference_sm_tables.resolve()}
    for folder in reference_tables.values():
        if not folder.is_dir() or not list(folder.glob('*.csv')):
            ap.error('Missing independent response tables: ' + str(folder))
    validation = json.loads(args.reference_validation.read_text())
    if validation['status'] != 'passed' or sha(reference) != validation['source']['sha256']:
        ap.error('Reference macro does not match its passed validation receipt.')
    inputs = json.loads(args.input_receipt.read_text())['inputs']
    records = []
    for item in inputs:
        path = Path(item['path']).resolve()
        if sha(path) != item['sha256']:
            ap.error('Raw input identity changed: ' + str(path))
        records.append({**item, 'path': str(path), 'size': path.stat().st_size,
                        'mtime_ns': path.stat().st_mtime_ns})
    sources = [HERE / 'Deconvolution.C', *sorted((HERE / 'root').glob('*')),
               *sorted((HERE / 'vendor').rglob('*.C')), reference, Path(__file__).resolve()]
    sources += sorted((HERE / 'data').rglob('*.csv')) + [HERE / 'data/channels.json', HERE / 'data/geometry.md5']
    for folder in reference_tables.values():
        sources += sorted(folder.glob('*.csv'))
    fingerprints = {str(p): sha(p) for p in sources if p.is_file()}
    receipt = {'status': 'running', 'scope': 'Raw ROOT to complete correlated fit; both weighting modes',
               'histories': sum(x['events'] for x in records), 'inputs': records,
               'reference_validation_sha256': sha(args.reference_validation),
               'source_sha256': fingerprints, 'comparisons': {},
               'settings': {'minimum_kinetic_MeV': 1, 'beam_days': 344, 'polarization': .8, 'current_uA': 65},
               'analysis_contract': {
                   'detectors': '224 physical quartz tiles + 28 ShowerMax entrance planes, compared in recorded order',
                   'reference': 'Physical Ring 5 Open tiles; rate-weighted charged crossings K>=1 MeV, both directions, for H in both modes',
                   'main_rate': 'Charged crossings K>=1 MeV, both directions',
                   'main_pe': 'Forward e+/- K>1.1 MeV; total-energy PE lookup with polar cap 4 degrees',
                   'showermax': 'Forward entries, total-energy expected PE in the recorded module frame',
                   'weights': 'Event rate once in signal and physical B; squared in finite-MC products; current and per-process file normalization',
                   'components': 'Moller, elastic, inelastic W=[1,1.4),[1.4,2.5),[2.5,6) GeV',
                   'units': 'Main crossings/s or PE/s; SM PE/s; reference ppb; C ppb^2; M ppb^-2',
                   'reference_tables': {k: str(v) for k, v in reference_tables.items()}},
               'limitation': 'Checks consistency on these inputs, not all possible inputs or independent validity of the shared estimator.'}
    try:
        with tempfile.TemporaryDirectory(prefix='consistency-', dir=args.output.parent) as tmp:
            tmp = Path(tmp)
            input_list = tmp / 'inputs.txt'
            input_list.write_text('\n'.join(x['path'] for x in records) + '\n')
            for name, folder in reference_tables.items():
                target = tmp / 'reference-tables' / name
                target.mkdir(parents=True)
                for table in folder.glob('*.csv'):
                    shutil.copyfile(table, target / table.name)
            for label, source in [('reference', reference), ('app', HERE / 'Deconvolution.C')]:
                driver = tmp / (label + '.C')
                calls = []
                for mode in (False, True):
                    output = tmp / f'{label}-{int(mode)}.root'
                    if label == 'reference':
                        calls.append('fnms={' + ','.join(json.dumps(x['path']) for x in records) + '};')
                        calls.append(f'deconvolution_correlated({str(mode).lower()}, {json.dumps(str(tmp / "reference-tables"))}, {json.dumps(str(output))});')
                    else:
                        calls.append(f'Deconvolution({json.dumps(str(input_list))}, {json.dumps(str(output))}, {str(mode).lower()});')
                driver.write_text('#include ' + json.dumps(str(source)) + '\nvoid ' + label + '(){\n' + '\n'.join(calls) + '\n}\n')
                command = ['root', '-n', '-l', '-b', '-e', 'gSystem->SetBuildDir(' + json.dumps(str(tmp / ('build-' + label))) + ',true);', '-q', str(driver) + '+']
                print('RUN', label, 'both weighting modes', flush=True)
                result = subprocess.run(command, cwd=HERE, capture_output=True, text=True, timeout=1200,
                                        env={**os.environ, 'TMPDIR': str(tmp)})
                if result.returncode or not all((tmp / f'{label}-{i}.root').is_file() for i in (0, 1)):
                    raise RuntimeError(label + ' failed:\n' + (result.stdout + result.stderr)[-6000:])
            for mode in (0, 1):
                actual = ROOT.TFile.Open(str(tmp / f'app-{mode}.root'))
                expected = ROOT.TFile.Open(str(tmp / f'reference-{mode}.root'))
                checks = {key: compare(array(actual.Get(key)), array(expected.Get(key)), key) for key in KEYS}
                receipt['comparisons']['pe' if mode else 'rate'] = checks
                for key, check in checks.items():
                    print('PASS' if check['passed'] else 'FAIL', mode, key, flush=True)
                actual.Close(); expected.Close()
            for path, digest in fingerprints.items():
                if sha(path) != digest:
                    raise RuntimeError('Source changed during verification: ' + path)
            for item in records:
                stat = Path(item['path']).stat()
                if (stat.st_size, stat.st_mtime_ns) != (item['size'], item['mtime_ns']):
                    raise RuntimeError('Input changed during verification: ' + item['path'])
        receipt['status'] = 'passed' if all(x['passed'] for c in receipt['comparisons'].values() for x in c.values()) else 'failed'
    except Exception as error:
        receipt.update(status='failed', error=str(error))
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(receipt['status'].upper(), args.output)
    if 'error' in receipt:
        print(receipt['error'])
    return 0 if receipt['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
