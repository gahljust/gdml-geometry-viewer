"""Geometry slice integration check. Use a separate test server, not an active user's load.

python3 tests/test_geometry.py http://127.0.0.1:8786 /path/to/remoll
"""
import json
import math
import sys
import time
from urllib.error import HTTPError
from urllib.request import Request, urlopen

url, remoll = sys.argv[1:]


def api(path, data=None):
    req = Request(url + path, data=None if data is None else json.dumps(data).encode(),
                  headers={'Content-Type': 'application/json'})
    with urlopen(req, timeout=10) as response:
        return json.load(response)


def load():
    api('/load', {'root': remoll})
    deadline = time.monotonic() + 120
    while time.monotonic() < deadline:
        state = api('/status')
        if not state['busy']:
            assert state['ready'], state
            return api('/geometry')
        time.sleep(.1)
    raise AssertionError('Geometry load timed out')


first = load()
assert first['units'] == 'mm'
assert first['mesh_count'] == len(first['meshes']) > 0
assert first['max_transform_error_mm'] < 1e-7
assert set(m['system'] for m in first['meshes']) == {*(f'Ring {i}' for i in range(1, 7)), 'ShowerMax'}
for mesh in first['meshes']:
    assert all(len(v) == 3 and all(math.isfinite(x) for x in v) for v in mesh['vertices'])
    assert all(len(f) >= 3 and all(0 <= i < len(mesh['vertices']) for i in f) for f in mesh['faces'])
assert first['vertex_count'] == sum(len(m['vertices']) for m in first['meshes'])
surfaces=[m for m in first['meshes'] if m['hit_surface']]
assert len(surfaces)==252 and len({m['geometry_detector'] for m in surfaces})==252
assert {m['geometry_detector'] for m in surfaces if m['system']=='ShowerMax'}==set(range(73001,73029))
print('PASS: 224 main tiles + 28 ShowerMax virtual planes; quartz layers are context only')
print('PASS: six rings + ShowerMax, finite vertices, valid facets, Geant4 transform check')
try:
    api('/load', {'root': remoll + '/missing-directory'})
    raise AssertionError('Invalid directory accepted')
except HTTPError as error:
    assert error.code == 400
assert api('/geometry') == first
print('PASS: invalid directory reports an error and preserves loaded geometry')
assert load() == first
print('PASS: repeated load produces identical geometry')
