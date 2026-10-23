"""Tiny develop-schema fixtures; no Geant4 simulation. Replaces the demo load.
Run: python3 test_root.py [http://127.0.0.1:8785]
Fixtures/build products are created visibly in build/ and deleted after testing.
"""

import collections
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from urllib.request import urlopen, Request

HERE = Path(__file__).resolve().parent
from run import DEFAULT_REMOLL

if not DEFAULT_REMOLL:
    raise SystemExit("Set REMOLL_DIR to your remoll checkout.")
ROOT = Path(DEFAULT_REMOLL)
URL = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8785"


def api(path, data=None):
    with urlopen(
        Request(
            URL + path,
            data=None if data is None else json.dumps(data).encode(),
            headers={"Content-Type": "application/json"},
        ),
        timeout=20,
    ) as r:
        return json.load(r)


def wait():
    end = time.monotonic() + 60
    while time.monotonic() < end:
        s = api("/status")
        if not s["busy"]:
            return s
        time.sleep(0.03)
    raise AssertionError("Reader timeout")


api("/stop", {})
wait()
api("/load", {"root": str(ROOT), "delay": 0.01})
api("/finish", {})
assert wait()["stage"] == "ready"
geometry = api("/geometry")["data"]
main = [r["detector"] for r in geometry["rows"] if r["system"] != "showermax"]
planes = [r["plane_detector"] for r in geometry["rows"] if r["system"] == "showermax"]
quartz = [
    i
    for r in geometry["rows"]
    if r["system"] == "showermax"
    for i in r["gdml_detectors"]
]
cases = []


def case(name, events, rates=None, kind="normal", error=None):
    cases.append(dict(name=name, events=events, rates=rates, kind=kind, error=error))


for mode in ["moller", "elastic", "inelastic", "beam", "pion", "external"]:
    # Generator/kinematic branches differ; the reader must only depend on det/rate.
    case(
        mode,
        [[150252, 150252, 73003, 999], [150252, 150240, 73003], []],
        [2.0, 3.0, 11.0],
    )
for name, ids in [
    ("all_detectors", main + planes),
    ("main_only", main),
    ("showermax_only", planes),
    ("physical_quartz", quartz),
    ("both_sm_types", quartz + planes),
    ("all_disabled", []),
    ("unmatched", [-1, 0, 999, 2147483647]),
    ("one_tile", [150252]),
    ("many_hits", [150252] * 20000),
]:
    case(name, [ids], [7.0])
for ring in range(1, 7):
    case(
        "ring" + str(ring),
        [
            [
                r["detector"]
                for r in geometry["rows"]
                if r["system"] == "ring" + str(ring)
            ]
        ],
        [1.25],
    )
case("empty_tree", [], [])
case("zero_rate", [[150252], [150252]], [0.0, 0.0])
case("tiny_rate", [[150252], [150252]], [1e-250, 2e-250])
case("vector_det", [[150252, 73003], [150252]], [2.0, 3.0], "vector_det")
case(
    "leaf_array_det",
    [[150252, 73003], [150252]],
    [2.0, 3.0],
    "leaf_array_det",
    "Unsupported leaf-list hit.det",
)
case("missing_rate", [[150252], [73003]], None, "missing_rate")
for name, val in [
    ("nan", float("nan")),
    ("infinity", float("inf")),
    ("negative", -1.0),
]:
    case(name, [[150252], [73003]], [2.0, val])
case("overflow", [[150252, 150252]], [1.7e308])
case("float_rate", [[150252]], [2.0], "float_rate")
case("array_rate", [[150252]], [2.0], "array_rate")
case("short_rate", [[150252], [73003]], [2.0, 3.0], "short_rate")
case("missing_tree", [[150252]], [2.0], "missing_tree", "Missing tree T")
case("missing_det", [[150252]], [2.0], "missing_det", "Missing field hit.det")
case("wrong_det", [[150252]], [2.0], "wrong_det", "hit.det must contain integers")
case(
    "short_det",
    [[150252], [73003]],
    [2.0, 3.0],
    "short_det",
    "Incomplete hit.det branch",
)
case("unsplit_hit", [[150252]], [2.0], "unsplit_hit", "Missing field hit.det")
# literal header from the local develop ref, not our modified struct
if os.environ.get("REMOLL_TEST_SCHEMA_HEADER"):
    header = Path(os.environ["REMOLL_TEST_SCHEMA_HEADER"]).read_text()
    commit = os.environ["REMOLL_TEST_SCHEMA_REVISION"]
else:
    header = subprocess.check_output(
        ["git", "-C", str(ROOT), "show", "develop:include/remolltypes.hh"], text=True
    )
    commit = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "develop"], text=True
    ).strip()
with tempfile.TemporaryDirectory(dir=HERE / "build", prefix="root-check-") as temp:
    d = Path(temp)
    (d / "remolltypes.hh").write_text(header)
    code = [
        "#include <vector>",
        "#include <limits>",
        "#include <TFile.h>",
        "#include <TTree.h>",
        "#include <TNamed.h>",
        "#include <TInterpreter.h>",
        '#include "remolltypes.hh"',
        "void fixtures(){",
        'gInterpreter->GenerateDictionary("vector<remollGenericDetectorHit_t>;vector<remollEventParticle_t>","remolltypes.hh;vector");',
    ]
    for n, c in enumerate(cases):
        code += [
            "{",
            f'TFile f("case {n}.root","RECREATE");',
            f'TTree t("'
            + ("wrong" if c["kind"] == "missing_tree" else "T")
            + '","fixture");',
            "std::vector<remollGenericDetectorHit_t> hits;std::vector<remollEventParticle_t> part;remollEvent_t ev{};",
            "std::vector<double> wrong;std::vector<int> integer;int n=0,a[10];double rate=0,ar[2]={1,2};float fr=1;",
        ]
        if c["kind"] == "wrong_det":
            code += ['t.Branch("hit.det",&wrong);']
        elif c["kind"] == "vector_det":
            code += ['t.Branch("hit.det",&integer);']
        elif c["kind"] == "leaf_array_det":
            code += ['t.Branch("n",&n);t.Branch("hit.det",a,"det[n]/I");']
        elif c["kind"] != "missing_det":
            code += [
                f't.Branch("hit",&hits,32000,{0 if c["kind"]=="unsplit_hit" else 99});'
            ]
        code += ['t.Branch("part",&part);t.Branch("ev",&ev);']
        kind = c["kind"]
        if kind == "float_rate":
            code += ['t.Branch("rate",&fr,"rate/F");']
        elif kind == "array_rate":
            code += ['t.Branch("rate",ar,"rate[2]/D");']
        elif kind != "missing_rate":
            code += ['t.Branch("rate",&rate,"rate/D");']
        code += ["part.resize(" + ("2" if c["name"] == "moller" else "1") + ");"]
        for j, ids in enumerate(c["events"]):
            val = c["rates"][j] if c["rates"] else 0
            v = (
                "std::numeric_limits<double>::quiet_NaN()"
                if math.isnan(val)
                else (
                    "std::numeric_limits<double>::infinity()"
                    if math.isinf(val)
                    else repr(val)
                )
            )
            if len(ids) > 1000:
                fill = f"for(int i=0;i<{len(ids)};i++){{remollGenericDetectorHit_t h{{}};h.det={ids[0]};hits.push_back(h);}}"
            else:
                fill = (
                    "for(int id:std::vector<int>{"
                    + ",".join(map(str, ids))
                    + "}){remollGenericDetectorHit_t h{};h.det=id;h.pid=(id%2?22:11);h.k=.5;h.pz=-1;hits.push_back(h);wrong.push_back(id);integer.push_back(id);}"
                )
            code += [
                "hits.clear();wrong.clear();integer.clear();"
                + fill
                + "n=integer.size();for(int i=0;i<n&&i<10;i++)a[i]=integer[i];",
                f"rate={v};t.Fill();",
            ]
        if kind == "short_det":
            code += ['t.GetBranch("hit.det")->SetEntries(1);']
        if kind == "short_rate":
            code += ['t.GetBranch("rate")->SetEntries(1);']
        code += [
            't.Write();TNamed mode("fixture_generator","'
            + c["name"]
            + '");mode.Write();}'
        ]
    code += ["}"]
    macro = d / "fixtures.C"
    macro.write_text("\n".join(code))
    proc = subprocess.run(
        ["root", "-l", "-b", "-q", str(macro)], cwd=d, capture_output=True, text=True
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    size = sum(p.stat().st_size for p in d.glob("*.root"))

    def analyze(n, delay=0.01):
        api("/analyze", {"path": str(d / f"case {n}.root"), "delay": delay})

    for n, c in enumerate(cases):
        analyze(n)
        s = wait()
        if c["error"]:
            assert s["stage"] == "error" and c["error"] in s["detail"], (c["name"], s)
        else:
            h = s["hits"]
            assert s["stage"] == "ready" and h["complete"], (c["name"], s)
            expected = dict(collections.Counter(str(i) for e in c["events"] for i in e))
            assert h["counts"] == expected and h["events"] == len(c["events"]), (
                c["name"],
                h,
            )
            if (
                c["rates"] is None
                or c["kind"]
                not in ("normal", "short_det", "vector_det", "leaf_array_det")
                or any(not math.isfinite(v) or v < 0 for v in c["rates"])
                or c["name"] == "overflow"
            ):
                assert h["rates"] is None and h["rate_error"], (c["name"], h)
            else:
                expected_rate = collections.defaultdict(float)
                for e, w in zip(c["events"], c["rates"]):
                    for i in e:
                        expected_rate[str(i)] += w
                assert h["rates"].keys() == expected_rate.keys(), (c["name"], h)
                assert all(
                    math.isclose(h["rates"][i], v, rel_tol=1e-12, abs_tol=0)
                    for i, v in expected_rate.items()
                ), (c["name"], h)
        print("PASS", c["name"], flush=True)
    # A renamed text file must not look like an empty successful ROOT scan.
    bad = d / "bad.root"
    bad.write_text("not a ROOT file")
    api("/analyze", {"path": str(bad)})
    assert wait()["stage"] == "error"
    print("PASS corrupt file")
    analyze(0, 1.6)
    end = time.monotonic() + 15
    while not api("/status")["hits"]:
        assert time.monotonic() < end
        time.sleep(0.03)
    api("/pause", {})
    s = api("/status")
    time.sleep(0.2)
    assert api("/status")["hits"] == s["hits"]
    api("/stop", {})
    s = wait()
    assert (
        not s["hits"]["complete"]
        and s["child_pid"] is None
        and api("/geometry")["data"]
    )
    analyze(0)
    assert wait()["hits"]["complete"]
    print("PASS pause, stop, rerun")
    print(
        json.dumps(
            {
                "schema_commit": commit,
                "fixtures": len(cases) + 1,
                "total_root_bytes": size,
            }
        )
    )
# Leave no test file selected in the product.
api("/load", {"root": str(ROOT), "delay": 0.01})
api("/finish", {})
assert wait()["stage"] == "ready"
