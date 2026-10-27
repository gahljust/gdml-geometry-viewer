#!/usr/bin/env python3
"""Run known-answer and failure-handling checks; delete generated fixtures afterward."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parents[1]
CHECKS = ('solver', 'root-inputs', 'responses', 'cache', 'launch', 'package', 'pe', 'hits', 'saved-view')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--remoll', type=Path, help='remoll source directory; also REMOLL_DIR')
    parser.add_argument('--work-dir', type=Path, default=HERE / 'build', help='Temporary test directory parent')
    parser.add_argument('--verbose', action='store_true', help='Show successful test details too')
    parser.add_argument('--only', choices=CHECKS, help='Run one check instead of the full suite')
    args = parser.parse_args()
    if not shutil.which('root'):
        parser.error('ROOT is required. Run in the terminal where ROOT is configured.')
    remoll = args.remoll or (Path(os.environ['REMOLL_DIR']) if os.environ.get('REMOLL_DIR') else
                            next((p for p in HERE.parents if (p / 'geometry/mollerMother.gdml').is_file()), None))
    if remoll:
        remoll = remoll.expanduser().resolve()
        if not (remoll / 'geometry/mollerMother.gdml').is_file():
            parser.error('--remoll must point to a remoll source directory.')
    parent = args.work_dir.expanduser().resolve()
    if any(c.isspace() for c in str(parent)):
        parser.error('ROOT/ACLiC needs a test build path without spaces. Set --work-dir accordingly.')
    parent.mkdir(parents=True, exist_ok=True)
    passed = failed = skipped = 0
    with tempfile.TemporaryDirectory(prefix='checks-', dir=parent) as temporary:
        scratch = Path(temporary)
        print(f'Temporary fixtures: {scratch} (removed when finished)', flush=True)
        environment = {**os.environ, 'PYTHONDONTWRITEBYTECODE': '1',
                       'MOLLER_OUTPUT': str(scratch / 'output'),
                       'TMPDIR': str(scratch), 'TMP': str(scratch), 'TEMP': str(scratch)}
        if remoll:
            environment['REMOLL_DIR'] = str(remoll)
        (scratch / 'output').mkdir()
        root_driver = scratch / 'root_inputs.C'
        root_driver.write_text('#include ' + json.dumps((HERE / 'tests/test_root_inputs.C').as_posix()) +
                              '\nvoid root_inputs(){test_root_inputs(' + json.dumps((scratch / 'fixtures').as_posix()) + ');}\n')
        root_command = ['root', '-n', '-l', '-b', '-e',
                        'gSystem->SetBuildDir(' + json.dumps((scratch / 'aclic').as_posix()) + ',true);', '-q']
        python = lambda name: [sys.executable, '-B', str(HERE / 'tests' / name)]
        checks = [
            ('Known-answer deconvolution', root_command + [str(HERE / 'tests/test_root_solver.C')], 'PASS: exact statistical sigmas'),
            ('ROOT metadata and PE table validation', root_command + [str(root_driver) + '+'], 'PASS: real ROOT metadata'),
            ('Response resource validation', python('test_response_sources.py'), None),
            ('Cache invalidation', python('test_result_cache.py'), None),
            ('Launcher', python('test_launch.py'), None),
            ('Source package', python('test_package.py'), None),
            ('Known-answer PE rates from ROOT files', python('test_weighted_maps.py') if remoll else None, 'PASS all-step PE weighting rejected'),
            ('ROOT cuts, histograms and malformed files', python('test_hit_loader.py') if remoll and shutil.which('node') else None, 'PASS reload after errors'),
            ('Saved views and display units', ['node', str(HERE / 'tests/test_save_view.js')] if shutil.which('node') else None, None),
        ]
        for name, (label, command, marker) in zip(CHECKS, checks):
            if args.only and args.only != name:
                continue
            if command is None:
                skipped += 1
                print(f'SKIP {label} — needs ' + ('--remoll' if 'PE rates' in label else 'Node.js' + (' and --remoll' if 'ROOT cuts' in label else '')), flush=True)
                continue
            print(f'RUN  {label}', flush=True)
            try:
                result = subprocess.run(command, cwd=HERE, env=environment, capture_output=True, text=True, timeout=600)
                output = result.stdout + result.stderr
                ok = result.returncode == 0 and (marker is None or marker in output)
            except (OSError, subprocess.TimeoutExpired) as error:
                ok, output = False, str(error)
            passed += int(ok)
            failed += int(not ok)
            print(f'{"PASS" if ok else "FAIL"} {label}', flush=True)
            if not ok or args.verbose:
                print(output[-16000:], flush=True)
    print(f'\n{passed} passed · {failed} failed · {skipped} skipped. Temporary files removed.')
    if skipped:
        print('Skipped checks are untested, not passes. Supply --remoll and Node.js to run the skipped checks.')
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
