"""Launch/discovery contracts; no ROOT, Geant4 or platform emulation required."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch, MagicMock

HERE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('launcher_test', HERE / 'geometry/run.py')
app = importlib.util.module_from_spec(spec)
spec.loader.exec_module(app)


class LaunchTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='detector-launch-')
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name).resolve()
        self.remoll = self.base / 'remoll with spaces'
        (self.remoll / 'geometry').mkdir(parents=True)
        (self.remoll / 'geometry/mollerMother.gdml').write_text('<gdml/>')
        self.g4 = self.base / 'Geant4 config'
        self.g4.mkdir()
        (self.g4 / 'Geant4Config.cmake').write_text('# test config')
        self.patches = patch.multiple(app, OUTPUT=self.base / 'build', DEFAULT='', GEANT4_DIR='')
        self.patches.start()
        self.addCleanup(self.patches.stop)

    def test_discovery(self):
        cache = self.remoll / 'build-release/CMakeCache.txt'
        cache.parent.mkdir()
        cache.write_text(f'Geant4_DIR:PATH={self.g4}\n')
        self.assertEqual(app.geant4_directory(self.remoll), self.g4)
        self.assertIsNone(app.geant4_directory(None))
        with patch.object(app, 'GEANT4_DIR', str(self.g4)):
            self.assertEqual(app.geant4_directory(None), self.g4)
        with patch.object(app, 'GEANT4_DIR', str(self.base / 'missing')):
            with self.assertRaisesRegex(ValueError, 'Geant4Config.cmake'):
                app.geant4_directory(self.remoll)
        with patch.object(app.shutil, 'which', return_value=None):
            with self.assertRaisesRegex(ValueError, 'not on PATH'):
                app.require_tool('root')

    def test_geometry_build_commands(self):
        # Both generator layouts use bin/, including a native Windows .exe.
        for platform, executable in [('posix', 'export_geometry'), ('nt', 'export_geometry.exe')]:
            commands = []

            def execute(command, log, env=None):
                commands.append([str(v) for v in command])
                if str(command[0]).endswith(executable):
                    command[-1].write_text(json.dumps({'meshes': [], 'mesh_count': 0}))

            with patch.object(app.os, 'name', platform), \
                 patch.object(app, 'require_tool', return_value='cmake'), \
                 patch.object(app, 'geant4_directory', return_value=self.g4), \
                 patch.object(app, 'execute', side_effect=execute):
                app.load(self.remoll)
            self.assertTrue(app.state['ready'], app.state['message'])
            self.assertIn('-DCMAKE_BUILD_TYPE=Release', commands[0])
            self.assertIn('-DGeant4_DIR=' + str(self.g4), commands[0])
            self.assertEqual(commands[1][-4:], ['--config', 'Release', '--parallel', '2'])
            self.assertEqual(commands[2][0], str(app.OUTPUT / 'bin' / executable))
            self.assertTrue((app.OUTPUT / 'surface-faces.tsv').is_file())

    def test_cli_failures(self):
        with self.assertRaisesRegex(ValueError, 'requires --remoll'):
            app.main(['--geometry-only'])
        with self.assertRaisesRegex(ValueError, 'mollerMother.gdml'):
            app.main(['--remoll', str(self.base / 'absent'), '--check'])
        with self.assertRaisesRegex(ValueError, 'Port must'):
            app.main(['--port', '0'])
        with self.assertRaisesRegex(ValueError, 'separate build directory'):
            app.main(['--output', str(HERE)])

    def test_aclic_build_path(self):
        with patch.object(app, 'OUTPUT', self.base / 'build with spaces'):
            with self.assertRaisesRegex(ValueError, 'ROOT input filenames may contain spaces'):
                app.root_compiler()

    def test_explicit_remoll_starts_geometry_without_browser(self):
        server = MagicMock()
        with patch.object(app, 'ThreadingHTTPServer', return_value=server), \
             patch.object(app.threading, 'Thread') as thread, \
             patch.object(app.webbrowser, 'open') as browser, patch.object(app, 'stop'):
            self.assertEqual(app.main(['--remoll', str(self.remoll), '--no-open']), 0)
        self.assertEqual(thread.call_args.kwargs['args'], (self.remoll,))
        thread.return_value.start.assert_called_once()
        browser.assert_not_called()
        server.server_close.assert_called_once()


if __name__ == '__main__':
    unittest.main()
