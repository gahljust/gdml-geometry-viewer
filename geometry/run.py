#!/usr/bin/env python3
"""Serve detector geometry, ROOT hit/PE analysis and ROOT deconvolution locally."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'app'))
import responses
REPO = next((p for p in HERE.parents if (p / 'geometry/mollerMother.gdml').is_file()), None)
DEFAULT = os.environ.get('REMOLL_DIR', str(REPO or ''))
GEANT4_DIR = os.environ.get('Geant4_DIR', '')
OUTPUT = Path(os.environ.get('MOLLER_OUTPUT', str(
    REPO / 'rootfiles/local_work/experiments/detector_launcher/geometry-build'
    if REPO else HERE.parent / 'build'))).expanduser().resolve()
RESPONSE_CONFIG = OUTPUT / 'response-config.json'
lock = threading.RLock()
state = {'busy': False, 'message': 'Choose a remoll directory.', 'revision': 0, 'ready': False}
process = None
closing = threading.Event()
fit_state = {'busy': False, 'message': 'Select ROOT files, then check inputs.', 'revision': 0}
fit_data = None

# Recent selections stay in bounded memory only; no automatic saved-analysis files.
CACHE_JSON_LIMIT = 8 * 1024 * 1024  # per result; memory only
CACHE_MEMORY_LIMIT = 32 * 1024 * 1024  # Python objects per result
result_cache = {}


def file_identity(path):
    path = Path(path).resolve()
    st = path.stat()
    return [str(path), st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns]


def result_key(kind, paths, settings, response=None, *, content_state=False):
    """No ROOT scan: file identity + exact cuts/settings + analyzer/resource hashes."""
    source = {HERE / 'run.py', HERE / 'load_hits.C', HERE / 'deconvolution_bridge.C', HERE.parent / 'Deconvolution.C'}
    for folder in ('root', 'app', 'vendor'):
        source.update(p for p in (HERE.parent / folder).rglob('*') if p.suffix in ('.C', '.cc', '.hh', '.h', '.py'))
    portable = HERE.parent / 'vendor/portable'
    source.update(portable / name for name in ('scan.C', 'detector_selection.hh'))
    source.update(HERE.parent / 'data' / name for name in ('channels.json', 'geometry.md5'))
    digest = hashlib.sha256()
    for path in sorted(source):
        digest.update(os.path.relpath(path, HERE.parent).encode()); digest.update(path.read_bytes())
    faces = OUTPUT / 'surface-faces.tsv'
    geometry = hashlib.sha256(faces.read_bytes()).hexdigest() if kind == 'hits' and state['ready'] and faces.is_file() else None
    root = shutil.which('root')
    sidecars = [Path(p).with_suffix('.neyman.tsv') for p in paths]
    identities = [file_identity(p) for p in paths]
    # ROOT may change metadata timestamps while reading. The cache lookup uses
    # ctime; the before/after read guard uses inode, size and mtime instead.
    if content_state: identities = [entry[:-1] for entry in identities]
    payload = [1, kind, identities,
               [file_identity(p) if p.exists() else None for p in sidecars], settings, digest.hexdigest(),
               response if response else None, geometry, file_identity(root) if root else None]
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def small_result(value):
    seen = set(); total = 0; pending = [value]
    while pending:
        item = pending.pop()
        if id(item) in seen: continue
        seen.add(id(item)); total += sys.getsizeof(item)
        if total > CACHE_MEMORY_LIMIT: return False
        if isinstance(item, dict): pending.extend(item.keys()); pending.extend(item.values())
        elif isinstance(item, (list, tuple)): pending.extend(item)
    size = 0
    for chunk in json.JSONEncoder(separators=(',', ':'), allow_nan=False).iterencode(value):
        size += len(chunk.encode())
        if size > CACHE_JSON_LIMIT: return False
    return True


def remember_result(kind, key, value):
    entries = [e for e in result_cache.get(kind, []) if e['key'] != key]
    if small_result(value):
        entries.append(dict(key=key, value=value))
        entries = entries[-(8 if kind == 'hits' else 1):]
        while entries and not small_result(entries): entries.pop(0)
    result_cache[kind] = entries


def cached_result(kind, key):
    entries = result_cache.get(kind, [])
    for index, entry in enumerate(entries):
        if entry['key'] == key:
            entries.append(entries.pop(index))
            return entry['value']
    return None



def choose_directory(prompt="Choose remoll"):
    """Open the native folder picker; nothing is copied or uploaded."""
    if sys.platform == 'darwin':
        command = ['osascript', '-e', 'POSIX path of (choose folder with prompt '+json.dumps(prompt)+')']
    elif os.name == 'nt':
        command = ['powershell', '-NoProfile', '-STA', '-Command',
                   'Add-Type -AssemblyName System.Windows.Forms; $d=New-Object System.Windows.Forms.FolderBrowserDialog; if($d.ShowDialog() -eq "OK"){$d.SelectedPath}']
    else:
        command = (['zenity', '--file-selection', '--directory', '--title='+prompt]
                   if shutil.which('zenity') else [sys.executable, '-c',
                   'import tkinter as t; from tkinter.filedialog import askdirectory; r=t.Tk();r.withdraw();print(askdirectory());r.destroy()'])
    result = subprocess.run(command, capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else ''


def require_tool(name):
    path = shutil.which(name)
    if not path:
        raise ValueError(f'{name} is not on PATH. Activate your existing installation in this terminal.')
    return path


def root_compiler():
    # ACLiC can emit an unquoted linker rpath even when SetBuildDir is relative.
    if any(c.isspace() for c in str(OUTPUT)):
        raise ValueError('Use a build directory without spaces for ROOT/ACLiC: relaunch with --output /path/without/spaces. ROOT input filenames may contain spaces.')
    return require_tool('root')


def remoll_directory(value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('Choose your remoll source directory.')
    root = Path(value).expanduser().resolve()
    if not (root / 'geometry/mollerMother.gdml').is_file():
        raise ValueError('Choose the remoll directory containing geometry/mollerMother.gdml.')
    return root


def geant4_directory(root):
    """Honor explicit CMake configuration, then reuse remoll's existing build cache."""
    if GEANT4_DIR:
        path = Path(GEANT4_DIR).expanduser().resolve()
        if not (path / 'Geant4Config.cmake').is_file():
            raise ValueError('Geant4_DIR must name the directory containing Geant4Config.cmake.')
        return path
    if root:
        for cache in sorted(root.glob('build*/CMakeCache.txt')):
            match = re.search(r'^Geant4_DIR:[^=]+=(.+)$', cache.read_text(errors='replace'), re.M)
            if match and (Path(match[1]) / 'Geant4Config.cmake').is_file():
                return Path(match[1]).resolve()
    return None  # CMake searches its standard locations and CMAKE_PREFIX_PATH.


def check_environment(root):
    print(f'Python {sys.version.split()[0]} · {sys.executable}')
    for name in ('cmake', 'root'):
        tool = root_compiler() if name == 'root' else require_tool(name)
        result = subprocess.run([tool, '--version'], capture_output=True, text=True, timeout=20)
        if result.returncode:
            raise ValueError(f'{name} could not start. Activate its runtime libraries.\n{result.stderr.strip()}')
        print(f'{name}: {tool}')
    g4 = geant4_directory(root)
    print(f'Geant4: {g4 or "CMake search; set --geant4-dir if needed"}')
    print(f'remoll: {root or "select in the app"}\nOutput: {OUTPUT}')
    print('The geometry build checks the C++ compiler and GDML support.')


def execute(command, log, env=None):
    """Track one child process so Quit/Ctrl+C also stops compilation/loading."""
    global process
    with lock:
        if closing.is_set():
            raise RuntimeError('Stopped')
        process = subprocess.Popen([str(x) for x in command], stdout=log, stderr=log,
                                   start_new_session=os.name != 'nt', env=env)
        child = process
    code = child.wait()
    with lock:
        process = None
    if code:
        raise RuntimeError(f'{Path(str(command[0])).name} failed; see {log.name}')


def surface_faces(meshes):
    """Campaign face coordinates: outward normal, local u/v, clipped 10 mm cells."""
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    sub=lambda a,b:[x-y for x,y in zip(a,b)]
    cross=lambda a,b:[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
    unit=lambda v:[x/math.sqrt(dot(v,v)) for x in v]
    faces=[]
    for mesh in meshes:
        if not mesh['hit_surface']: continue
        points=mesh['vertices'];center=[sum(p[i] for p in points)/len(points) for i in range(3)]
        for indices in mesh['faces']:
            poly=[points[i] for i in indices];origin=poly[0];u=unit(sub(poly[1],origin))
            normal=unit(cross(sub(poly[1],origin),sub(poly[2],origin)))
            if dot(normal,sub(origin,center))<0: normal=[-v for v in normal]
            v=cross(normal,u)
            faces.append(dict(detector=mesh['geometry_detector'],system='showermax' if mesh['system']=='ShowerMax' else 'ring'+mesh['system'][-1],
                origin=origin,u=u,v=v,normal=normal,polygon=[[dot(sub(p,origin),u),dot(sub(p,origin),v)] for p in poly]))
    return faces


def load(root):
    """Build the helper, then re-read the selected checkout on every Load."""
    try:
        OUTPUT.mkdir(parents=True, exist_ok=True)
        gdml = root / 'geometry/mollerMother.gdml'
        cmake = require_tool('cmake')
        g4 = geant4_directory(root)
        with (OUTPUT / 'load.log').open('w') as log:
            execute([cmake, '-S', HERE, '-B', OUTPUT, '-DCMAKE_BUILD_TYPE=Release',
                     *(['-DGeant4_DIR=' + str(g4)] if g4 else [])], log)
            execute([cmake, '--build', OUTPUT, '--config', 'Release', '--parallel', '2'], log)
            with lock:
                state['message'] = 'Reading GDML…'
            executable = OUTPUT / 'bin' / ('export_geometry.exe' if os.name == 'nt' else 'export_geometry')
            execute([executable, gdml, OUTPUT / 'geometry.next.json'], log)
        data = json.loads((OUTPUT / 'geometry.next.json').read_text())
        data['faces']=surface_faces(data['meshes'])
        (OUTPUT / 'geometry.next.json').write_text(json.dumps(data))
        (OUTPUT / 'surface-faces.tsv').write_text(''.join(' '.join(map(str,[i,f['detector'],*f['origin'],*f['u'],*f['v'],*f['normal'],len(f['polygon']),*[x for uv in f['polygon'] for x in uv]]))+'\n' for i,f in enumerate(data['faces'])))
        (OUTPUT / 'geometry.next.json').replace(OUTPUT / 'geometry.json')
        with lock:
            state.update(ready=True, revision=state['revision'] + 1,
                         message=f'{data["mesh_count"]} detector meshes', root=str(root))
    except Exception as error:
        with lock:
            state.update(ready=False, message=str(error))
    finally:
        with lock:
            state['busy'] = False


def stop():
    closing.set()
    with lock:
        child = process
    if child and child.poll() is None:
        if os.name == 'nt':
            subprocess.run(['taskkill', '/PID', str(child.pid), '/T', '/F'], capture_output=True)
        else:
            try:
                os.killpg(child.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *_):
        pass

    def reply(self, value, code=200, content_type='application/json'):
        body = value if isinstance(value, bytes) else json.dumps(value, allow_nan=False).encode()
        self.send_response(code)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('X-Frame-Options', 'DENY')
        self.send_header('Content-Security-Policy', "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; object-src 'none'; frame-ancestors 'none'; base-uri 'none'")
        self.end_headers()
        self.wfile.write(body)

    def local_request(self):
        """Reject remote origins and forged Host headers before exposing local paths."""
        hosts = {f'127.0.0.1:{self.server.server_port}', f'localhost:{self.server.server_port}'}
        supplied = self.headers.get_all('Host', [])
        origin = self.headers.get('Origin')
        if len(supplied) != 1 or supplied[0] not in hosts:
            self.reply({'error': 'Local host only'}, 403)
            return False
        if origin and origin != 'http://' + supplied[0]:
            self.reply({'error': 'Local page only'}, 403)
            return False
        if self.headers.get('Sec-Fetch-Site') == 'cross-site' and self.path != '/':
            self.reply({'error': 'Local page only'}, 403)
            return False
        return True

    def do_GET(self):
        if not self.local_request(): return
        if self.path == '/responses':
            try:
                config = responses.read_config(RESPONSE_CONFIG)
                self.reply({'config': config, 'resources': responses.resolve(config)})
            except (ValueError, OSError) as error:
                self.reply({'error': str(error)}, 400)
        elif self.path == '/status':
            with lock:
                result = state['hits']['result']
                compact = {k:v for k,v in result.items() if k not in ('cells','histograms','histograms2d','particles')} if result else None
                self.reply({**state, 'fit': fit_state, 'hits': {**state['hits'], 'result': compact}, 'default_root': DEFAULT})
        elif self.path == '/deconvolution':
            with lock:
                self.reply(fit_data)
        elif self.path == '/hits':
            with lock:
                self.reply(state['hits']['result'])
        elif self.path == '/geometry':
            with lock:
                if not state['ready']:
                    self.reply({'error': 'Load geometry first.'}, 409)
                    return
                self.reply((OUTPUT / 'geometry.json').read_bytes())
        elif self.path in ['/', '/view.js', '/hit-map-3d.js', '/source/run.py', '/source/export_geometry.cc',
                           '/source/CMakeLists.txt', '/source/view.js', '/source/index.html', '/deconvolution.js', '/save-view.js']:
            name = 'index.html' if self.path == '/' else self.path.rsplit('/', 1)[-1]
            kind = 'text/html' if self.path == '/' else 'text/javascript' if self.path in ['/view.js','/hit-map-3d.js','/deconvolution.js','/save-view.js'] else 'text/plain'
            self.reply((HERE / name).read_bytes(), content_type=kind + '; charset=utf-8')
        else:
            self.reply({'error': 'Not found'}, 404)

    def do_POST(self):
        if not self.local_request(): return
        try:
            if self.headers.get('Transfer-Encoding') or len(self.headers.get_all('Content-Length', [])) != 1:
                raise ValueError('Provide one Content-Length header.')
            size = int(self.headers.get('Content-Length', '0'))
            if not 0 <= size <= 131072:
                raise ValueError('Request size must be 0–131072 bytes.')
            self.connection.settimeout(10)
            def nonfinite(value):
                raise ValueError('Nonfinite JSON number: '+value)
            data = json.loads(self.rfile.read(size) or '{}', parse_constant=nonfinite)
            if not isinstance(data, dict):
                raise ValueError('Request must be a JSON object.')
            if self.path == '/responses':
                with lock:
                    if state['busy']:
                        self.reply({'error': 'Wait for the current analysis to finish.'}, 409)
                        return
                    resolved = responses.save_config(RESPONSE_CONFIG, data)
                self.reply({'config': data, 'resources': resolved})
            elif self.path in ('/deconvolution/check', '/deconvolution/run'):
                paths = validate_paths(data.get('paths'))
                extra = data.get('extra', [])
                if not isinstance(extra, list) or len(extra)>8 or any(not isinstance(x,str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',x) for x in extra):
                    raise ValueError('Use up to eight additional generator names.')
                if len(set(extra))!=len(extra) or set(extra)&{'moller','elastic','inelastic'}:
                    raise ValueError('Background names must be unique.')
                forecast={k:data.get(k,default) for k,default in [('days',344),('polarization',.8),('current',65)]}
                if any(type(v) not in (int,float) or not 0<v<1e9 for v in forecast.values()) or forecast['polarization']>1:
                    raise ValueError('Use positive beam days/current and polarization between 0 and 1.')
                resolved=responses.resolve(responses.read_config(RESPONSE_CONFIG))
                with lock:
                    if state['busy']:
                        self.reply({'error':'Wait for the current calculation.'},409);return
                    state['busy']=True
                    fit_state.update(busy=True,message='Checking ROOT files…' if self.path.endswith('check') else 'Running ROOT deconvolution…')
                    threading.Thread(target=run_fit,args=(paths,extra,forecast,resolved,self.path.endswith('check')),daemon=True).start()
                self.reply({'ok':True})
            elif self.path == '/choose':
                self.reply({'path': choose_directory()})
            elif self.path == '/load':
                root = remoll_directory(data.get('root', ''))
                with lock:
                    if state['busy']:
                        self.reply({'error': 'Geometry is already loading.'}, 409)
                        return
                    state.update(busy=True, ready=False, message='Building geometry reader…')
                    state['hits'].update(result=None,message='Geometry changed. Apply the ROOT cut again.')
                    threading.Thread(target=load, args=(root,), daemon=True).start()
                self.reply({'ok': True})
            elif self.path == '/quit':
                stop()
                self.reply({'ok': True})
                threading.Thread(target=self.server.shutdown, daemon=True).start()
            elif self.path == '/choose-hits':
                self.reply({'paths': choose_hit_file()})
            elif self.path == '/load-hits':
                paths = validate_paths(data.get('paths', [data.get('path', '')]))
                if not isinstance(data.get('cut', ''), str):
                    raise ValueError('The ROOT cut must be text.')
                cut = data.get('cut', '').strip() or '1'
                if len(cut) > 4096 or '\0' in cut:
                    raise ValueError('Invalid cut expression.')
                binning = validate_binning(data.get('binning', {}))
                weighting = data.get('weighting', 'hits')
                if weighting not in ('hits', 'pe'):
                    raise ValueError('Choose hit counts or PE weighting.')
                resolved = responses.resolve(responses.read_config(RESPONSE_CONFIG)) if weighting == 'pe' else None
                with lock:
                    if state['busy']:
                        self.reply({'error': 'A file is already loading.'}, 409)
                        return
                    state['busy'] = True
                    state['hits'] = dict(loading=True, path=str(paths[0]), paths=list(map(str, paths)), cut=cut, binning=binning, weighting=weighting, message='Reading hits…', result=None)
                    threading.Thread(target=load_hit_file, args=(paths, cut, binning, resolved), daemon=True).start()
                self.reply({'ok': True})
            else:
                self.reply({'error': 'Not found'}, 404)
        except (ValueError, OSError, RecursionError) as error:
            self.reply({'error': str(error)}, 400)


state['hits'] = dict(loading=False, path='', message='All stored hits · No cuts', result=None)


def choose_hit_file():
    """Use the standard OS file picker; no upload or file copy."""
    if sys.platform == 'darwin':
        # Use AppKit's native dialog directly; avoid AppleScript alias coercion.
        command = ['osascript', '-l', 'JavaScript', '-e', '''
function run() {
    ObjC.import('AppKit');
    const app = $.NSApplication.sharedApplication;
    app.setActivationPolicy($.NSApplicationActivationPolicyAccessory);
    app.activateIgnoringOtherApps(true);
    const panel = $.NSOpenPanel.openPanel;
    panel.title = 'Choose remoll ROOT files';
    panel.canChooseFiles = true;
    panel.canChooseDirectories = false;
    panel.allowsMultipleSelection = true;
    if (Number(panel.runModal) !== Number($.NSModalResponseOK)) return '[]';
    const paths = [];
    for (let i=0; i<Number(panel.URLs.count); i++) paths.push(ObjC.unwrap(panel.URLs.objectAtIndex(i).path));
    return JSON.stringify(paths);
}''']
    elif os.name == 'nt':
        command = ['powershell', '-NoProfile', '-STA', '-Command',
                   'Add-Type -AssemblyName System.Windows.Forms; $d=New-Object System.Windows.Forms.OpenFileDialog; $d.Filter="ROOT files (*.root)|*.root"; $d.Multiselect=$true; if($d.ShowDialog() -eq "OK"){ConvertTo-Json -InputObject @($d.FileNames) -Compress}else{"[]"}']
    else:
        command = (['zenity', '--file-selection', '--multiple', '--separator=\n', '--title=Choose remoll ROOT files']
                   if shutil.which('zenity') else [sys.executable, '-c',
                   'import tkinter as t; from tkinter.filedialog import askopenfilenames; import json; r=t.Tk();r.withdraw();print(json.dumps(list(askopenfilenames(filetypes=[("ROOT files","*.root"),("All files","*")]))));r.destroy()'])
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        if sys.platform.startswith('linux') and result.returncode == 1 and not result.stderr.strip():
            return []  # Zenity cancellation.
        raise ValueError(result.stderr.strip() or 'Could not open the file picker.')
    value = result.stdout.strip()
    return value.splitlines() if command[0] == 'zenity' else json.loads(value or '[]')


# Only named histogram quantities are configurable; no arbitrary ROOT code here.
PLOT_KEYS = set('energy total_energy edep radius vertex vertex_r x y z time momentum px py pz phi theta spin_x spin_y spin_z asymmetry q2 w rate beam_p'.split())
def validate_binning(value):
    if not isinstance(value, dict) or len(value)>len(PLOT_KEYS):
        raise ValueError('Invalid histogram bin settings.')
    result = {}
    for key, item in value.items():
        if key not in PLOT_KEYS or not isinstance(item, dict):
            raise ValueError('Unknown histogram quantity.')
        n, low, high, spacing = (item.get(k) for k in ('bins','min','max','spacing'))
        if type(n) is not int or not 1 <= n <= 256:
            raise ValueError('Use 1–256 bins per axis.')
        if any(type(v) not in (int,float) or abs(v)>1e15 or not math.isfinite(v) for v in (low,high)) or low>=high:
            raise ValueError('Histogram limits must be finite and increasing.')
        if spacing not in ('linear','log') or (spacing=='log' and low<=0):
            raise ValueError('Logarithmic bins require positive limits.')
        if spacing=='log' and not math.isfinite(high/low):
            raise ValueError('Histogram range is too large.')
        if low+(high-low)/n == low:
            raise ValueError('Histogram bins are too narrow at this scale.')
        result[key] = dict(bins=n,min=low,max=high,spacing=spacing)
    return result


def merge_hits(total, part):
    """Sum ROOT-produced counts; never mix or renormalize event rate weights."""
    if total is None:
        return part
    for key in ('events', 'hits', 'stored_hits', 'selected_events', 'invalid_positions',
                'unmatched_positions', 'other_detectors', 'invalid_particles'):
        total[key] += part[key]
    for key in ('detectors', 'particles', 'detector_signal', 'particle_signal'):
        for item, count in part.get(key, {}).items():
            total[key][item] = total[key].get(item, 0) + count
    for key in ('signal_total','unsupported_pe'):
        total[key] = total.get(key,0)+part.get(key,0)
    total['fields'] = sorted(set(total['fields']) & set(part['fields']))
    for key in ('positions_available', 'surfaces_available', 'particles_available'):
        total[key] = total[key] and part[key]
    for name, hist in total['histograms'].items():
        other = part['histograms'][name]
        if hist['edges'] != other['edges']:
            raise ValueError('Histogram bins differ between files.')
        hist['available'] = hist['available'] and other['available']
        hist['counts'] = [a+b for a,b in zip(hist['counts'], other['counts'])]
        for key in ('underflow','overflow','invalid'):
            hist[key] += other[key]
    for name, hist in total['histograms2d'].items():
        other = part['histograms2d'][name]
        if (hist['x'], hist['y']) != (other['x'], other['y']):
            raise ValueError('2D histogram axes differ between files.')
        hist['available'] = hist['available'] and other['available']
        for key in ('invalid', 'outside'):
            hist[key] += other[key]
        bins = {(x,y):n for x,y,n in hist['cells']}
        for x,y,n in other['cells']:
            bins[x,y] = bins.get((x,y),0) + n
        hist['cells'] = [[x,y,n] for (x,y),n in sorted(bins.items())]
    cells = {tuple(c[:3]):c[3] for c in total['cells']}
    for *cell, count in part['cells']:
        key = tuple(cell)
        cells[key] = cells.get(key, 0) + count
    if len(cells) > 200000:
        raise ValueError('Too many occupied surface bins. Narrow the cut.')
    total['cells'] = [list(key)+[count] for key,count in sorted(cells.items())]
    return total


def load_hit_file(paths, cut='1', binning=None, response=None):
    paths = [paths] if isinstance(paths, Path) else paths
    binning = validate_binning(binning or {})
    bin_text = '\n'.join(f'{k} {v["bins"]} {v["min"]:.17g} {v["max"]:.17g} {int(v["spacing"]=="log")}' for k,v in binning.items())
    receipt = OUTPUT / 'hit-counts.next.json'
    try:
        key = result_key('hits', paths, [cut, binning, bool(response)], response)
        cached = cached_result('hits', key)
        if cached is not None:
            with lock: state['hits'].update(result=cached, message=f'Reused analysis · {cached["events"]:,} events')
            return
        unchanged = result_key('hits', paths, [cut, binning, bool(response)], response, content_state=True)
        root = root_compiler()
        OUTPUT.mkdir(parents=True, exist_ok=True)
        total = None
        groups = {}
        if response:
            folder = responses.stage(response, OUTPUT / 'response-resources')
            driver = response_driver(folder, 'WeightedHits', '#define MOLLER_WEIGHTED_HITS\n#include '+cpp_path(HERE.parent/'Deconvolution.C')+'\n#include '+cpp_path(HERE/'load_hits.C'), 'load_hits();')
        for index, path in enumerate(paths):
            if closing.is_set():
                raise ValueError('Stopped')
            with lock:
                state['hits']['message'] = f'Reading file {index+1} / {len(paths)} · {path.name}'
            receipt.unlink(missing_ok=True)
            log_path = OUTPUT / 'hits.log'
            try:
                with log_path.open('w') as log:
                    command = [root, '-n', '-l', '-b']
                    if not response:
                        command += ['-e', 'gSystem->SetBuildDir('+json.dumps(str(OUTPUT))+',true);']
                    command += ['-q', str(driver if response else HERE / 'load_hits.C')+'+']
                    execute(command, log,
                            env={**os.environ, 'MOLLER_HIT_INPUT': str(path), 'MOLLER_HIT_OUTPUT': str(receipt), 'MOLLER_HIT_CUT': cut, 'MOLLER_HIT_BINS': bin_text, 'MOLLER_HIT_FACES': str(OUTPUT / 'surface-faces.tsv') if state['ready'] else '', **({'MOLLER_MAIN_PE':str(folder/'main_pe'),'MOLLER_SM_PE':str(folder/'response')} if response else {})})
            except RuntimeError:
                errors = [line.split('HIT_ERROR ', 1)[1] for line in log_path.read_text(errors='replace').splitlines() if 'HIT_ERROR ' in line]
                raise ValueError(f'{path.name}: ' + (errors[-1] if errors else f'ROOT could not read the file. See {log_path}.')) from None
            part = json.loads(receipt.read_text())
            if response:
                if any(g['part']['transport'] != part['transport'] for g in groups.values()):
                    raise ValueError('Incompatible transport: source, geometry, magnetic fields or beam/physics settings differ. Analyze separately.')
                cohort=part['normalization']
                if any(g['part']['generator']==part['generator'] and k!=cohort for k,g in groups.items()):
                    raise ValueError('The same generator has different run settings or geometry. Analyze those files separately.')
                if cohort not in groups:groups[cohort]={'seeds':set(), 'part':None}
                if part['seed'] in groups[cohort]['seeds']:raise ValueError('Repeated seed in PE input cohort.')
                groups[cohort]['seeds'].add(part['seed'])
                groups[cohort]['part']=merge_hits(groups[cohort]['part'],part)
            else:total = merge_hits(total, part)
        if response:
            for group in groups.values():
                part=group['part'];scale_signal(part,1/part['events']);total=merge_hits(total,part)
            total['response_source']=response
            total['unit']='PE/s'
            total['weighting']='pe'
            total['beam_current_uA']=65
        total['file_count'] = len(paths)
        # Do not display a partial histogram/map as though it covered every file.
        if not total['positions_available']:
            total['cells'] = []
        if result_key('hits', paths, [cut, binning, bool(response)], response, content_state=True) != unchanged:
            raise ValueError('Inputs or analysis code changed during analysis. Load again.')
        remember_result('hits', result_key('hits', paths, [cut, binning, bool(response)], response), total)
        with lock:
            state['hits'].update(result=total, message=f'{len(paths)} file(s) · {total["hits"]:,} / {total["stored_hits"]:,} hits · {total["selected_events"]:,} / {total["events"]:,} events')
    except Exception as error:
        with lock:
            state['hits'].update(result=None, message=str(error))
    finally:
        receipt.unlink(missing_ok=True)
        with lock:
            state['hits']['loading'] = False
            state['busy'] = False



def validate_paths(values):
    """The same local-file contract for maps and fits; never split a filename."""
    if not isinstance(values, list) or not 1 <= len(values) <= 256:
        raise ValueError('Choose 1–256 ROOT files.')
    paths = []
    for value in values:
        if not isinstance(value, str) or not value.strip():
            raise ValueError('Each file path must be text.')
        if any(c in value for c in ('\0', '\r', '\n')):
            raise ValueError('ROOT filenames cannot contain line breaks or null characters.')
        path = Path(value).expanduser()
        if path.is_symlink() and not path.exists():
            raise ValueError('This file link is unavailable. Reconnect its drive or choose the original file.')
        path = path.resolve()
        if not path.is_file():
            raise ValueError(f'File not found: {path}')
        if any(path.samefile(other) for other in paths):
            raise ValueError(f'File selected twice: {path.name}')
        paths.append(path)
    return paths


def cpp_path(path):
    return json.dumps(str(path).replace('\\', '/'))


def response_driver(folder,name,includes,body):
    driver=folder/(name+'.C')
    code = ('#define MOLLER_MAIN_LOOKUP_HEADER '+cpp_path(folder/'code/lookuptable.C')+'\n'
        '#define MOLLER_MAIN_COORDINATE_HEADER '+cpp_path(folder/'code/convert_labframe_to_local_thinquartz.C')+'\n'
        +includes+'\nvoid '+name+'(){'+body+'}\n')
    if not driver.is_file() or driver.read_text() != code: driver.write_text(code)
    return driver


def scale_signal(part,factor):
    part['signal_total']*=factor
    for key in ('detector_signal','particle_signal'):
        part[key]={k:v*factor for k,v in part[key].items()}
    for cell in part['cells']:cell[3]*=factor
    for hist in part['histograms'].values():
        hist['counts']=[x*factor for x in hist['counts']]
        for key in ('underflow','overflow'):hist[key]*=factor
    for hist in part['histograms2d'].values():
        for cell in hist['cells']:cell[2]*=factor


def run_fit(paths,extra,forecast,response,check):
    global fit_data
    result=OUTPUT/'deconvolution.next.json';root_result=OUTPUT/'deconvolution.next.root'
    try:
        root = root_compiler()
        key = result_key('fit', paths, [extra, forecast], response)
        cached = cached_result('fit', key) if not check else None
        if cached is not None:
            with lock:
                fit_data = cached
                fit_state.update(revision=fit_state['revision']+1, message='Reused last deconvolution result.')
            return
        if not check: result_cache.pop('fit', None)
        unchanged = result_key('fit', paths, [extra, forecast], response, content_state=True)
        folder=responses.stage(response,OUTPUT/'response-resources')
        inputs=OUTPUT/'deconvolution-inputs.txt';inputs.write_text('\n'.join(map(str,paths))+'\n')
        driver=response_driver(folder,'FitDriver','#include '+cpp_path(HERE/'deconvolution_bridge.C'),'deconvolution_bridge();')
        result.unlink(missing_ok=True);root_result.unlink(missing_ok=True)
        env={**os.environ,'MOLLER_INPUT_LIST':str(inputs),'MOLLER_FIT_JSON':str(result),'MOLLER_FIT_ROOT':str(root_result),
             'MOLLER_FIT_ACTION':'check' if check else 'run','MOLLER_EXTRA':' '.join(extra),
             'MOLLER_MAIN_PE':str(folder/'main_pe'),'MOLLER_SM_PE':str(folder/'response'),
             'MOLLER_DAYS':str(forecast['days']),'MOLLER_POLARIZATION':str(forecast['polarization']),'MOLLER_CURRENT':str(forecast['current'])}
        log_path=OUTPUT/'deconvolution.log'
        try:
            with log_path.open('w') as log:
                if not check:
                    execute([root,'-n','-l','-b','-q',str(driver)+'+'],log,{**env,'MOLLER_FIT_ACTION':'check'})
                    if not json.loads(result.read_text())['ready']:check=True
                if not check or not result.is_file():
                    execute([root,'-n','-l','-b','-q',str(driver)+'+'],log,env)
        except RuntimeError:
            errors=[x.split('FIT_ERROR ',1)[1] for x in log_path.read_text(errors='replace').splitlines() if 'FIT_ERROR ' in x]
            raise ValueError(errors[-1] if errors else f'ROOT failed. See {log_path}') from None
        data=json.loads(result.read_text());data.update(response_source=response,forecast=forecast,paths=list(map(str,paths)),extra=extra)
        if result_key('fit', paths, [extra, forecast], response, content_state=True) != unchanged:
            raise ValueError('Inputs or analysis code changed during analysis. Run again.')
        if not check: remember_result('fit', result_key('fit', paths, [extra, forecast], response), data)
        if not check:root_result.replace(OUTPUT/'deconvolution.root')
        result.write_text(json.dumps(data,separators=(',',':'),allow_nan=False))
        result.replace(OUTPUT/'deconvolution.json')
        with lock:
            fit_data=data
            fit_state.update(revision=fit_state['revision']+1,message='Ready to calculate.' if check and data['ready'] else 'Resolve the listed inputs.' if check else 'ROOT calculation complete.')
    except Exception as error:
        with lock:
            fit_data=None
            fit_state.update(revision=fit_state['revision']+1,message=str(error))
    finally:
        result.unlink(missing_ok=True);root_result.unlink(missing_ok=True)
        with lock:state['busy']=False;fit_state['busy']=False

def main(argv=None):
    global DEFAULT, OUTPUT, RESPONSE_CONFIG, GEANT4_DIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--remoll', help='remoll source directory; also REMOLL_DIR')
    parser.add_argument('--output', type=Path, help='Build/result directory; also MOLLER_OUTPUT')
    parser.add_argument('--geant4-dir', default=GEANT4_DIR, help='Directory containing Geant4Config.cmake')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--check', action='store_true', help='Check existing tools without installing anything')
    modes.add_argument('--geometry-only', action='store_true', help='Build and export geometry JSON, then exit')
    parser.add_argument('--port', type=int, default=8784)
    parser.add_argument('--no-open', action='store_true')
    args = parser.parse_args(argv)
    selected = args.remoll if args.remoll is not None else DEFAULT
    root = remoll_directory(selected) if selected else None
    DEFAULT = str(root or '')
    GEANT4_DIR = args.geant4_dir
    if args.output:
        OUTPUT = args.output.expanduser().resolve()
    if OUTPUT in (HERE, HERE.parent):
        raise ValueError('Choose a separate build directory, for example --output build.')
    RESPONSE_CONFIG = OUTPUT / 'response-config.json'
    if args.check:
        check_environment(root)
        return 0
    if args.geometry_only:
        if not root:
            raise ValueError('--geometry-only requires --remoll or REMOLL_DIR.')
        load(root)
        if not state['ready']:
            raise ValueError(state['message'])
        print(f'{OUTPUT / "geometry.json"} · {state["message"]}')
        return 0
    if not 1 <= args.port <= 65535:
        raise ValueError('Port must be between 1 and 65535.')
    server = ThreadingHTTPServer(('127.0.0.1', args.port), Handler)
    address = f'http://127.0.0.1:{args.port}/'
    print(f'{address} · Quit or Ctrl+C to stop. Build files: {OUTPUT}', flush=True)
    if root and args.remoll is not None:
        state.update(busy=True, message='Building geometry reader…')
        threading.Thread(target=load, args=(root,), daemon=True).start()
    if not args.no_open:
        webbrowser.open(address)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        stop()
        server.server_close()
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, subprocess.TimeoutExpired) as error:
        sys.exit(str(error))
