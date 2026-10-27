#!/usr/bin/env python3
"""Export the current app only; no research data, builds, or retired launchers."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

HERE = Path(__file__).resolve().parents[1]
FILES = (
    'README.md', 'AGENTS.md', 'THIRD_PARTY.md', '.gitignore',
    'launch', 'launch.py', 'Start-Mac.command', 'Start-Windows.cmd', 'Deconvolution.C',
    'app/responses.py', 'docs/ANALYSIS.md', 'docs/GUIDES.txt', 'tests/README.md', 'tests/run.py',
    'tests/test_hit_loader.py', 'tests/test_weighted_maps.py', 'tests/test_geometry.py',
    'tests/test_result_cache.py', 'tests/test_launch.py', 'tests/test_save_view.js', 'tests/test_response_sources.py',
    'tests/test_root_solver.C', 'tests/test_root_inputs.C',
    'tests/test_package.py', 'tools/build_package.py',
)
TREES = {'geometry': {'.py', '.js', '.html', '.cc', '.C', '.txt'},
         'root': {'.C', '.hh', '.h'}, 'data': {'.json', '.csv', '.md5'},
         'vendor': {'.C', '.hh'},
         'docs': {'.pdf'}}


def build(destination):
    destination = Path(destination).resolve()
    if destination == HERE or HERE in destination.parents:
        raise ValueError('Export outside the source directory.')
    if destination.exists():
        raise ValueError('Choose a new package directory.')
    sources = {HERE / name for name in FILES}
    for directory, suffixes in TREES.items():
        tree = HERE / directory
        candidates = tree.glob('*') if directory in {'geometry', 'docs'} else tree.rglob('*')
        sources.update(p for p in candidates
                       if p.is_file() and p.suffix in suffixes and '__pycache__' not in p.parts)
    missing = sorted(str(p) for p in sources if not p.is_file())
    if missing:
        raise ValueError('Missing package sources: ' + ', '.join(missing))
    copied = {}
    destination.mkdir(parents=True)
    for source in sorted(sources):
        relative = source.relative_to(HERE)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        copied[relative.as_posix()] = hashlib.sha256(target.read_bytes()).hexdigest()
    (destination / 'bundle.json').write_text(json.dumps(
        {'schema': 2, 'copied_sha256': copied}, indent=2) + '\n')
    return destination


def archive(destination):
    destination = Path(destination).resolve()
    archive_path = destination.with_name(destination.name + '.zip')
    if archive_path.exists():
        raise ValueError('Archive already exists: ' + str(archive_path))
    manifest = json.loads((destination / 'bundle.json').read_text())
    for name, expected in manifest['copied_sha256'].items():
        path = (destination / name).resolve()
        if destination not in path.parents or not path.is_file():
            raise ValueError('Invalid package entry: ' + name)
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError('Package changed after export: ' + name)
    with zipfile.ZipFile(archive_path, 'x', compression=zipfile.ZIP_DEFLATED) as output:
        for name in sorted(set(manifest['copied_sha256']) | {'bundle.json'}):
            output.write(destination / name, Path(destination.name) / name)
    return archive_path


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='New directory outside source')
    parser.add_argument('--archive-only', action='store_true')
    args = parser.parse_args()
    if not args.archive_only:
        build(args.output)
    print(archive(args.output))
