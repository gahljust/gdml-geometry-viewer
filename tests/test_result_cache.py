"""Bounded latest-result cache; no ROOT jobs or persistent test files."""
import importlib.util
import json
import os
from pathlib import Path
import tempfile
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('geometry_server_cache_test', HERE/'geometry/run.py')
app = importlib.util.module_from_spec(spec)
spec.loader.exec_module(app)
with tempfile.TemporaryDirectory(prefix='result-cache-test-') as folder:
    app.OUTPUT = Path(folder)
    source = app.OUTPUT/'input.root'; source.write_bytes(b'ROOT test identity')
    app.state['ready'] = True
    face = app.OUTPUT/'surface-faces.tsv'; face.write_text('geometry A')
    key = lambda cut='1', response=None: app.result_key('hits', [source], [cut, {}, False], response)
    first = key()
    assert first == key()
    assert first != key('hit.pid==11')
    assert key(response={'fingerprint':'a'}) != key(response={'fingerprint':'b'})
    face.write_text('geometry B'); assert first != key()
    previous = key(); old = source.stat(); source.write_bytes(b'ROOT changed data!')
    os.utime(source, ns=(old.st_atime_ns, old.st_mtime_ns))
    assert previous != key(), 'Content replacement with restored mtime must invalidate'
    previous = key(); original_read = Path.read_bytes
    with patch.object(Path, 'read_bytes', autospec=True, side_effect=lambda p: b'changed analyzer' if p == app.HERE/'load_hits.C' else original_read(p)):
        assert previous != key()
    value = {'events':12, 'cells':[[1,2,3,4]]}
    app.remember_result('hits', 'a', value)
    assert app.cached_result('hits','a') is value
    assert app.cached_result('hits','b') is None
    assert not (app.OUTPUT/'analysis-cache.json').exists(), 'Default must not save cache'
    app.remember_result('fit','c',{'matrices':{'M':[[2]]}})
    assert app.cached_result('fit','c')['matrices']['M'] == [[2]]
    for i in range(12): app.remember_result('hits',str(i),{'events':i})
    assert len(app.result_cache['hits']) == 8
    assert app.cached_result('hits','0') is None
    assert app.cached_result('hits','11') == {'events':11}
    app.CACHE_MEMORY_LIMIT=100
    app.remember_result('hits','large',{'cells':list(range(1000))})
    assert app.cached_result('hits','large') is None
    app.CACHE_MEMORY_LIMIT=32*1024*1024; app.CACHE_JSON_LIMIT=20
    assert not app.small_result({'value':'x'*30})
    assert not (app.OUTPUT/'analysis-cache.json').exists()
print('PASS: input/geometry/analyzer/response invalidation; bounded memory; no automatic disk cache.')
