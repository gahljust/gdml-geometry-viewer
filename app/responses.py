"""Explicit PE resource selection, validation and provenance. No silent fallback."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parents[1]
SPECIES = ('e-', 'e+', 'gamma', 'mu-', 'mu+', 'pi-', 'pi+', 'neutron')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_config(path):
    return json.loads(path.read_text()) if path.is_file() else {'mode': 'bundled'}


def table(path, columns):
    if not path.is_file():
        raise ValueError(f'Missing response file: {path}')
    with path.open(newline='', encoding='utf-8-sig') as stream:
        rows = csv.reader(stream)
        header = next(rows, [])
        if len(header) != columns:
            raise ValueError(f'{path.name}: expected {columns} columns.')
        values = []
        for line, row in enumerate(rows, 2):
            if not row:
                continue
            try:
                nums = [float(x) for x in row]
                if len(nums) != columns or not all(math.isfinite(x) for x in nums):
                    raise ValueError()
            except ValueError:
                raise ValueError(f'{path.name}: invalid numeric row {line}.') from None
            values.append(nums)
    if not values:
        raise ValueError(f'{path.name}: empty response table.')
    return values


def resolve(config, base=HERE):
    """Local means local files only; missing resources are actionable errors."""
    if not isinstance(config, dict) or config.get('mode') not in ('bundled', 'local'):
        raise ValueError('Choose bundled or local PE response.')
    mode = config['mode']
    if mode == 'bundled':
        main = base / 'data/main_pe'
        sm = base / 'data/response'
        code = base / 'vendor/main_response'
    else:
        for key in ('remoll', 'main_tables', 'showermax_tables'):
            if not isinstance(config.get(key, ''), str):
                raise ValueError(f'{key}: expected a directory path.')
        if not config.get('remoll', '').strip():
            raise ValueError('Choose your remoll directory.')
        remoll = Path(config['remoll']).expanduser().resolve()
        code = remoll / 'analysis'
        main = (Path(config['main_tables']).expanduser().resolve() if config.get('main_tables', '').strip()
                else code / 'table')
        sm = (Path(config['showermax_tables']).expanduser().resolve() if config.get('showermax_tables', '').strip()
              else code / 'showermax/response_lookup/data')
    files = {}
    for name in ('lookuptable.C', 'convert_labframe_to_local_thinquartz.C'):
        path = code / name
        if not path.is_file():
            raise ValueError(f'Missing response code: {path}')
        files['code/' + name] = path
    for ring in range(1, 7):
        # All three names occur in upstream examples / DocDB releases.
        choices = [main / f'r{ring}_{suffix}.csv' for suffix in
                   ('NewLookupTable', 'PositionDepTable', 'LookupTable')]
        available = [p for p in choices if p.is_file()]
        if not available:
            raise ValueError(f'Missing ring {ring} position coefficients in {main}. Select the extracted DocDB table folder.')
        if len({sha(p) for p in available}) > 1:
            raise ValueError(f'Conflicting ring {ring} position tables in {main}. Keep one version.')
        pos = available[0]
        values = table(pos, 13)
        if any(row[0] >= row[1] or row[2] >= row[3] for row in values):
            raise ValueError(f'{pos.name}: invalid position bin bounds.')
        energy = main / f'r{ring}_EnergyDepTable.csv'
        rows = table(energy, 2)
        # The current evaluator selects within ±1 MeV and assumes this grid.
        if len(rows) != 500 or any(abs(row[0] - (2.511 + 2*i)) > .001 for i, row in enumerate(rows)):
            raise ValueError(f'{energy.name}: unsupported energy grid; expected 500 rows at 2.511 + 2*i MeV.')
        if any(row[1] < 0 for row in rows):
            raise ValueError(f'{energy.name}: negative energy scaling.')
        files[f'main_pe/r{ring}_NewLookupTable.csv'] = pos
        files[f'main_pe/{energy.name}'] = energy
    for species in SPECIES:
        path = sm / f'fit_param_xy_{species}_ifarm.csv'
        rows = table(path, 6)
        seen = {}
        for row in rows:
            if row[0] <= 0 or (row[0] in seen and seen[row[0]] != row[1:]):
                raise ValueError(f'{path.name}: invalid energy or conflicting coefficients at one energy.')
            seen[row[0]] = row[1:]  # Upstream includes an identical repeated 1000 MeV row.
        files['response/' + path.name] = path
    manifest = json.loads((base / 'data/response-manifest.json').read_text())
    hashes = {name: sha(path) for name, path in files.items()}
    if mode == 'bundled':
        expected = manifest['sha256'] if 'sha256' in manifest else manifest['files']
        for name, digest in hashes.items():
            key = 'vendor/main_response/' + name[5:] if name.startswith('code/') else name
            target = manifest['code_sha256'].get(key) if name.startswith('code/') else expected.get(key)
            if digest != target:
                raise ValueError(f'Bundled response changed: {name}. Choose Local to use modified resources.')
    fingerprint = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    return {'mode': mode, 'files': {k: str(v.resolve()) for k, v in files.items()},
            'sha256': hashes, 'fingerprint': fingerprint,
            'source': manifest['main_source'] if mode == 'bundled' else {'origin': 'User selected local files'},
            'evaluator': 'Selected main fit code; fixed detector ID based ShowerMax adapter.'}


def stage(resolved, output):
    """Small immutable snapshot per content hash; aliases do not alter CSV bytes."""
    folder = output / resolved['fingerprint']
    for name, source in resolved['files'].items():
        data = Path(source).read_bytes()
        if hashlib.sha256(data).hexdigest() != resolved['sha256'][name]:
            raise ValueError('Response files changed during selection; select them again.')
        dest = folder / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.is_file() or sha(dest) != resolved['sha256'][name]:
            dest.write_bytes(data)
    (folder / 'provenance.json').write_text(json.dumps(resolved, indent=2) + '\n')
    return folder


def save_config(path, config):
    resolved = resolve(config)  # Validate before replacing a working selection.
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, suffix='.json', delete=False) as stream:
        temp = Path(stream.name)
        json.dump(config, stream, indent=2)
    temp.replace(path)
    return resolved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', choices=('bundled', 'local'), default='bundled')
    parser.add_argument('--remoll', default='')
    parser.add_argument('--main-tables', default='')
    parser.add_argument('--showermax-tables', default='')
    parser.add_argument('--config', type=Path, help='Read a saved source selection')
    parser.add_argument('--inputs', type=Path, help='Run the existing ROOT deconvolution on this file list')
    parser.add_argument('--output', type=Path, default=Path('deconvolution.root'))
    args = parser.parse_args()
    config = read_config(args.config) if args.config else dict(mode=args.source, remoll=args.remoll,
        main_tables=args.main_tables, showermax_tables=args.showermax_tables)
    resolved = resolve(config)
    if not args.inputs:
        print(json.dumps(resolved, indent=2))
        return
    output = args.output.resolve()
    if output.exists():
        raise ValueError('Output already exists; choose another file.')
    output.parent.mkdir(parents=True, exist_ok=True)
    folder = stage(resolved, output.parent / 'response-resources')
    # Includes are resolved before ROOT compilation, so Local uses local fit code too.
    quote = lambda p: json.dumps(str(p).replace('\\', '/'))
    driver = folder / 'ResponseDeconvolution.C'
    driver.write_text('#define MOLLER_MAIN_LOOKUP_HEADER ' + quote(folder / 'code/lookuptable.C') + '\n'
        '#define MOLLER_MAIN_COORDINATE_HEADER ' + quote(folder / 'code/convert_labframe_to_local_thinquartz.C') + '\n'
        '#include ' + quote(HERE / 'Deconvolution.C') + '\n'
        'void ResponseDeconvolution(){Deconvolution(' + quote(args.inputs.resolve()) + ',' + quote(output)
        + ',true,1,344,.8,65,' + quote(folder / 'main_pe') + ',' + quote(folder / 'response') + ');}\n')
    result = subprocess.run(['root', '-n', '-l', '-b', '-q', str(driver)+'+'])
    if result.returncode or not output.is_file():
        raise RuntimeError('ROOT deconvolution failed; see its diagnostic above.')
    output.with_suffix('.response.json').write_text(json.dumps(resolved, indent=2) + '\n')


if __name__ == '__main__':
    main()
