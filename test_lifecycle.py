#!/usr/bin/env python3
"""Exercise a running demo with real GDML. No third-party test packages.

Start run.py, then: python3 test_lifecycle.py http://127.0.0.1:8785
This replaces the current demo load. It does not shut down the server.
"""

import concurrent.futures
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError

URL = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8785"
from run import DEFAULT_REMOLL

if not DEFAULT_REMOLL:
    raise SystemExit("Set REMOLL_DIR to your remoll checkout.")
ROOT = DEFAULT_REMOLL
checks = []


def request(path, data=None):
    req = Request(
        URL + path,
        data=json.dumps(data).encode() if data is not None else None,
        headers={"Content-Type": "application/json"},
    )
    try:
        with urlopen(req, timeout=10) as r:
            return r.status, json.load(r)
    except HTTPError as e:
        return e.code, json.load(e)


def status():
    return request("/status")[1]


def wait(predicate, timeout=45):
    end = time.monotonic() + timeout
    while time.monotonic() < end:
        s = status()
        if predicate(s):
            return s
        time.sleep(0.03)
    raise AssertionError("Timed out: " + str(status()))


def check(name, condition=True):
    assert condition, name
    checks.append(name)
    print("PASS", name, flush=True)


def start(delay=0.8):
    code, d = request("/load", {"root": ROOT, "delay": delay})
    assert code == 200, d
    return d["run"]


def stop():
    pid = status()["child_pid"]
    request("/stop", {})
    s = wait(lambda s: not s["busy"])
    assert s["stage"] == "stopped" and s["child_pid"] is None
    if pid:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            pass
        else:
            raise AssertionError("Reader survived Stop")


request("/stop", {})
wait(lambda s: not s["busy"])
for data in (
    {},
    {"root": 42},
    {"root": ROOT, "delay": -1},
    {"root": ROOT, "delay": float("nan")},
):
    assert request("/load", data)[0] == 400
check("Malformed load requests rejected")
request("/load", {"root": "/not/remoll"})
wait(lambda s: s["stage"] == "error")
check("Bad directory reports error without killing app")
for i in range(3):
    start()
    stop()
check("Three immediate start/stop cycles; no surviving reader")

start()
s = wait(lambda s: s["stage"] == "placement", 90)
revision = s["revision"]
t = time.monotonic()
s = wait(lambda s: s["revision"] > revision)
check("Default stage remains visible for human reading", time.monotonic() - t >= 0.65)
request("/pause", {})
s = status()
rev = s["revision"]
count = s["mesh_count"]
time.sleep(1.1)
check(
    "Pause freezes the native reader",
    status()["revision"] == rev and status()["mesh_count"] == count,
)
for i in range(2):
    request("/step", {})
    s = wait(lambda s: s["revision"] > rev)
    check("Step advances exactly one trace " + str(i + 1), s["revision"] == rev + 1)
    rev = s["revision"]
    time.sleep(0.25)
    assert status()["revision"] == rev
check("Overlapping load rejected", request("/load", {"root": ROOT})[0] == 409)
stop()
check("Stop while paused terminates reader and keeps server alive")

hashes = []
for i in range(2):
    start()
    request("/finish", {})
    wait(lambda s: s["stage"] == "ready", 90)
    d = request("/geometry")[1]["data"]
    check(
        "Complete load " + str(i + 1), len(d["rows"]) == 252 and len(d["meshes"]) == 336
    )
    hashes.append(hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest())
check("Repeated full loads produce identical geometry", len(set(hashes)) == 1)
stop()
request("/stop", {})
check("Repeated Stop is harmless", status()["stage"] == "stopped")

with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
    answers = list(pool.map(lambda _: request("/load", {"root": ROOT})[0], range(5)))
check(
    "Rapid duplicate loads start exactly one reader",
    answers.count(200) == 1 and answers.count(409) == 4,
)
stop()
check("Recovery after concurrent clicks")
print(json.dumps({"passed": len(checks), "geometry_sha256": hashes[0]}, indent=2))
