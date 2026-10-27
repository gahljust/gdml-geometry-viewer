"""Bundled/local identity, malformed inputs and actual ROOT evaluator checks."""
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / 'app'))
import responses


class Sources(unittest.TestCase):
    def setUp(self):
        # Locate the research data root without placing large artifacts in /tmp.
        out = next((p / 'rootfiles/local_work/experiments/detector_launcher'
                    for p in HERE.parents if (p / 'rootfiles').is_dir()), HERE / 'results')
        out.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=out, prefix='response-test-')
        self.base = Path(self.temp.name)
        self.remoll = self.base / 'remoll'
        self.code = self.remoll / 'analysis'
        shutil.copytree(HERE / 'vendor/main_response', self.code)
        self.tables = self.base / 'docdb'
        shutil.copytree(HERE / 'data/main_pe', self.tables)
        shutil.copytree(HERE / 'data/response', self.code / 'showermax/response_lookup/data')
        self.config = dict(mode='local', remoll=str(self.remoll), main_tables=str(self.tables))

    def tearDown(self):
        self.temp.cleanup()

    def test_same_bytes_and_docdb_aliases(self):
        for path in self.tables.glob('*_NewLookupTable.csv'):
            path.rename(path.with_name(path.name.replace('NewLookupTable', 'PositionDepTable')))
        bundled, local = responses.resolve({'mode': 'bundled'}), responses.resolve(self.config)
        self.assertEqual(bundled['fingerprint'], local['fingerprint'])
        staged = responses.stage(local, self.base / 'staged')
        self.assertEqual((staged / 'main_pe/r5_NewLookupTable.csv').read_bytes(),
                         (self.tables / 'r5_PositionDepTable.csv').read_bytes())

    def test_errors_preserve_configuration(self):
        saved = self.base / 'config.json'
        responses.save_config(saved, {'mode': 'bundled'})
        original = saved.read_bytes()
        path = self.tables / 'r1_NewLookupTable.csv'
        path.unlink()
        with self.assertRaisesRegex(ValueError, 'Missing ring 1'):
            responses.save_config(saved, self.config)
        self.assertEqual(saved.read_bytes(), original)
        with self.assertRaisesRegex(ValueError, 'Choose bundled or local'):
            responses.resolve([])

    def test_invalid_numbers_and_conflicting_aliases(self):
        path = self.tables / 'r2_NewLookupTable.csv'
        original = path.read_text()
        path.write_text(original.replace('272', 'nan', 1) if '272' in original else original + 'nan\n')
        with self.assertRaisesRegex(ValueError, 'invalid numeric'):
            responses.resolve(self.config)
        path.write_text(original)
        path.with_name('r2_PositionDepTable.csv').write_text(original + '\n')
        with self.assertRaisesRegex(ValueError, 'Conflicting ring 2'):
            responses.resolve(self.config)

    def test_duplicate_showermax_energy(self):
        path = self.code / 'showermax/response_lookup/data/fit_param_xy_e-_ifarm.csv'
        with path.open('a') as stream:
            stream.write('1000,1,2,3,4,5\n')
        with self.assertRaisesRegex(ValueError, 'conflicting coefficients'):
            responses.resolve(self.config)

    def test_detect_changes_during_selection(self):
        result = responses.resolve(self.config)
        with (self.tables / 'r1_NewLookupTable.csv').open('a') as stream:
            stream.write('\n')
        with self.assertRaisesRegex(ValueError, 'changed during selection'):
            responses.stage(result, self.base / 'staged')

    def test_bundled_integrity(self):
        bundle = self.base / 'bundle'
        shutil.copytree(HERE / 'data', bundle / 'data')
        shutil.copytree(HERE / 'vendor', bundle / 'vendor')
        path = bundle / 'data/main_pe/r1_NewLookupTable.csv'
        path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'Bundled response changed'):
            responses.resolve({'mode': 'bundled'}, base=bundle)

    @unittest.skipUnless(shutil.which('root'), 'ROOT unavailable')
    def test_root_evaluators_and_local_code_override(self):
        def evaluate(config, label):
            selected = responses.resolve(config)
            folder = responses.stage(selected, self.base / 'staged')
            q = lambda p: json.dumps(str(p))
            driver = self.base / (label + '.C')
            driver.write_text('#define MOLLER_MAIN_LOOKUP_HEADER ' + q(folder / 'code/lookuptable.C') + '\n'
                '#define MOLLER_MAIN_COORDINATE_HEADER ' + q(folder / 'code/convert_labframe_to_local_thinquartz.C') + '\n'
                '#include ' + q(HERE / 'Deconvolution.C') + '\n'
                'void ' + label + '(){\nroot_deconv::MainModel model(' + q(folder / 'main_pe') + ');\n'
                'std::cout<<std::setprecision(17); for(int r=0;r<6;r++){auto c=model.position[r][0];'
                'double v=root_main_pe::cal_PE((c.hmin+c.hmax)/2,(c.vmin+c.vmax)/2,1.,30.,model.position[r]);'
                'v*=root_main_pe::cal_PE_energy(123.,model.energy[r]);std::cout<<"CHECK "<<v<<"\\n";}\n'
                'auto sm=portable_deconv::load(' + q(folder / 'response') + ');'
                'std::cout<<"CHECK "<<portable_deconv::pe(sm,73001,11,1000,1100,0)<<"\\n";}\n')
            proc = subprocess.run(['root', '-n', '-l', '-b', '-q', str(driver)], text=True, capture_output=True)
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            values = [float(line.split()[1]) for line in proc.stdout.splitlines() if line.startswith('CHECK ')]
            self.assertEqual(len(values), 7, proc.stdout + proc.stderr)
            return values
        bundled = evaluate({'mode': 'bundled'}, 'BundledCheck')
        local = evaluate(self.config, 'LocalCheck')
        self.assertEqual(bundled, local)
        # Deliberately altered local evaluator must actually alter the result.
        code = self.code / 'lookuptable.C'
        text = code.read_text()
        self.assertIn('return PE;', text)
        code.write_text(text.replace('return PE;', 'return 2*PE;'))
        changed = evaluate(self.config, 'ChangedCheck')
        for a, b in zip(bundled[:6], changed[:6]):
            self.assertGreater(a, 0)
            self.assertAlmostEqual(2*a, b)
        self.assertEqual(bundled[6], changed[6])


if __name__ == '__main__':
    unittest.main()
