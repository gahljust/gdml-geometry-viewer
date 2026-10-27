"""Check exported entry points, resources, hashes and exclusion of local artifacts."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

HERE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('package_builder', HERE / 'tools/build_package.py')
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


class PackageTests(unittest.TestCase):
    def test_export(self):
        # Exports must sit outside the source, even when the runner sets TMPDIR
        # to build/checks-* inside it. Remove the sibling directory on exit.
        with tempfile.TemporaryDirectory(prefix='detector-package-', dir=HERE.parent) as temporary:
            folder = builder.build(Path(temporary) / 'workspace')
            manifest = json.loads((folder / 'bundle.json').read_text())
            for name, expected in manifest['copied_sha256'].items():
                self.assertEqual(hashlib.sha256((folder / name).read_bytes()).hexdigest(), expected)
            self.assertFalse((folder / 'native').exists())
            self.assertFalse((folder / 'web').exists())
            self.assertFalse((folder / 'app/server.py').exists())
            self.assertFalse((folder / 'experiment.json').exists())
            self.assertFalse(list(folder.rglob('*.root')))
            self.assertFalse(list(folder.rglob('__pycache__')))
            self.assertFalse((folder / 'docs/visual').exists())
            self.assertTrue((folder / 'docs/GUIDES.txt').is_file())
            visuals = list((folder / 'docs').glob('*.pdf'))
            expected_pdfs = {p.relative_to(HERE) for p in (HERE / 'docs').glob('*.pdf')}
            self.assertTrue(expected_pdfs)
            self.assertEqual({p.relative_to(folder) for p in visuals}, expected_pdfs)
            environment = {**os.environ, 'PYTHONPATH': '', 'PYTHONDONTWRITEBYTECODE': '1',
                           'MOLLER_OUTPUT': str(Path(temporary) / 'results')}
            run = subprocess.run([sys.executable, '-B', str(folder / 'launch.py'), '--help'],
                                 cwd=temporary, env=environment, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn('--port', run.stdout)
            self.assertNotIn('--cache-results', run.stdout)
            archive = builder.archive(folder)
            with zipfile.ZipFile(archive) as package:
                self.assertEqual(len(package.namelist()), len(manifest['copied_sha256']) + 1)
            with self.assertRaises(ValueError):
                builder.build(folder)
            with self.assertRaises(ValueError):
                builder.archive(folder)
            # A user's local CMake build must not enter a later export.
            artifacts = folder / 'geometry/build'
            artifacts.mkdir()
            (artifacts / 'CMakeCache.txt').write_text('local paths')
            (artifacts / 'generated.C').write_text('local build')
            source = builder.HERE
            try:
                builder.HERE = folder
                repackaged = builder.build(Path(temporary) / 'repackaged')
                self.assertFalse((repackaged / 'geometry/build').exists())
            finally:
                builder.HERE = source
            archive.unlink()
            (folder / 'launch.py').write_text('modified')
            with self.assertRaisesRegex(ValueError, 'changed after export'):
                builder.archive(folder)


if __name__ == '__main__':
    unittest.main()
