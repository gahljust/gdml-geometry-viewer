#!/usr/bin/env python3
"""One local walkthrough: choose remoll, read GDML with Geant4, inspect the mesh."""

import argparse
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
import webbrowser
from urllib.request import urlopen
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_REMOLL = os.environ.get("REMOLL_DIR") or next(
    (str(p) for p in HERE.parents if (p / "geometry/mollerMother.gdml").is_file()), ""
)
# The same parser, polyhedron and placement operations as export_geometry.cc.
# TRACE messages are emitted by the actual executing code, not a replay.
CPP = r"""
#include "G4GDMLParser.hh"
#include "G4Polyhedron.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NavigationHistory.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <filesystem>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <regex>
#include <map>
#include <vector>
// Geometry traversal: accumulate local-to-world placements and export quartz.
G4NavigationHistory history;
double error=0; int vertices=0, meshes=0;
void point(std::ostream& o,const G4ThreeVector& p){o<<'['<<p.x()<<','<<p.y()<<','<<p.z()<<']';}
std::string jsonText(const std::string& text){std::string out="\"";for(char c:text){if(c=='\n')out+="\\n";else if(c=='\r')out+="\\r";else if(c=='\t')out+="\\t";else{if(c=='\\'||c=='\"')out+='\\';out+=c;}}return out+'\"';}
std::string placement="null";
std::vector<std::string> ancestry;
std::map<G4LogicalVolume*,bool> detectorBranches;
bool quartz(G4LogicalVolume* lv){
 const auto& name=lv->GetName();
 return std::regex_search(name,std::regex("R[1-6]-.*Quartz"))||std::regex_match(name,std::regex("logic_quartz_[0-9]+_[0-9]+"));
}
bool detectorBranch(G4LogicalVolume* lv){
 auto found=detectorBranches.find(lv);if(found!=detectorBranches.end())return found->second;
 bool keep=quartz(lv);for(int i=0;i<lv->GetNoDaughters();i++)keep=detectorBranch(lv->GetDaughter(i)->GetLogicalVolume())||keep;
 return detectorBranches[lv]=keep;
}
void rotation(std::ostream& o,const G4RotationMatrix& r){
 o<<'[';for(int i=0;i<3;i++){if(i)o<<',';o<<'[';for(int j=0;j<3;j++){if(j)o<<',';o<<r(i,j);}o<<']';}o<<']';
}
void trace(const char* stage,const std::string& detail=""){
 std::cout<<"TRACE {\"stage\":"<<std::quoted(stage)<<",\"detail\":"<<jsonText(detail)<<",\"placement\":"<<placement<<"}"<<std::endl;
 std::string resume; if(!std::getline(std::cin,resume))throw std::runtime_error("Reader stopped");
}
void walk(G4VPhysicalVolume* pv,const G4RotationMatrix& parentR,const G4ThreeVector& parentT,G4GDMLParser& parser,int depth=0){
 if(!detectorBranch(pv->GetLogicalVolume()))return;
 if(depth==0)history.SetFirstEntry(pv);else history.NewLevel(pv,kNormal,pv->GetCopyNo());
 auto R=parentR*pv->GetObjectRotationValue();
 auto T=parentR*pv->GetObjectTranslation()+parentT;
 auto* lv=pv->GetLogicalVolume();std::string name=lv->GetName();
 std::ostringstream info;info<<std::setprecision(12)<<"{\"volume\":"<<jsonText(pv->GetName())<<",\"logical\":"<<jsonText(name)
  <<",\"parent\":"<<(ancestry.empty()?"null":jsonText(ancestry.back()))<<",\"copy\":"<<pv->GetCopyNo()<<",\"quartz\":"<<(quartz(lv)?"true":"false")<<",\"path\":[";
 for(size_t i=0;i<ancestry.size();i++){if(i)info<<',';info<<jsonText(ancestry[i]);}
 info<<"],\"local_position\":";point(info,pv->GetObjectTranslation());info<<",\"local_rotation\":";rotation(info,pv->GetObjectRotationValue());
 info<<",\"world_position\":";point(info,T);info<<",\"world_rotation\":";rotation(info,R);info<<'}';placement=info.str();
 ancestry.push_back(pv->GetName());
 if(!quartz(lv))trace("volume",name);
 if(quartz(lv)){
  trace("placement",name);
  auto* poly=lv->GetSolid()->GetPolyhedron();
  if(!poly)throw std::runtime_error("No polyhedron: "+name);
  trace("mesh",name);
  int detector=0;for(auto& a:parser.GetVolumeAuxiliaryInformation(lv))if(a.type=="DetNo")detector=std::stoi(a.value);
  std::ostringstream out;out<<std::setprecision(12)<<"{\"name\":"<<std::quoted(name)<<",\"detector\":"<<detector<<",\"center\":";point(out,T);out<<",\"vertices\":[";
  for(int i=1;i<=poly->GetNoVertices();i++){
   auto v=poly->GetVertex(i);G4ThreeVector local(v.x(),v.y(),v.z());
   auto world=R*local+T;
   auto check=history.GetTopTransform().Inverse().TransformPoint(local);
   error=std::max(error,(world-check).mag());vertices++;
   if(i>1)out<<',';point(out,world);
   if(i==1){std::ostringstream detail;detail<<name<<"\nlocal [mm] = ";point(detail,local);detail<<"\nworld [mm] = ";point(detail,world);trace("vertices",detail.str());}
  }
  out<<"],\"faces\":[";
  for(int i=1;i<=poly->GetNoFacets();i++){
   int n,ix[4];poly->GetFacet(i,n,ix);
   if(i>1)out<<',';out<<'[';for(int j=0;j<n;j++){if(j)out<<',';out<<ix[j]-1;}out<<']';
  }
  out<<"]}";trace("faces",name);
  std::cout<<"MESH "<<out.str()<<std::endl;meshes++;
  trace("placed",name);
 }
 for(int i=0;i<lv->GetNoDaughters();i++)walk(lv->GetDaughter(i),R,T,parser,depth+1);
 ancestry.pop_back();
 if(depth>0)history.BackLevel();
}
int main(int argc,char** argv){try{
 if(argc!=2)return 2;
 G4UnitDefinition::BuildUnitsTable();
 new G4UnitDefinition("inch","in","Length",2.54*cm);new G4UnitDefinition("mil","mil","Length",.001*2.54*cm);
 auto gdml=std::filesystem::absolute(argv[1]);std::filesystem::current_path(gdml.parent_path());
 G4GDMLParser parser;parser.SetOverlapCheck(false);
 trace("parse",gdml.string());
 parser.Read(gdml.filename().string(),false);
 trace("parsed","Geant4 finished reading the GDML tree.");
 walk(parser.GetWorldVolume(),G4RotationMatrix(),G4ThreeVector(),parser);
 placement="null";
 std::ostringstream result;result<<"Meshes: "<<meshes<<"\nVertices: "<<vertices<<"\nPlacement check: "<<error<<" mm";trace("check",result.str());
 if(error>1e-7)throw std::runtime_error("Placement verification failed");
 trace("done");return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
"""
SNIPPETS = {
    "volume": [
        "auto* lv=pv->GetLogicalVolume();std::string name=lv->GetName();",
        "auto R=parentR*pv->GetObjectRotationValue();",
        "auto T=parentR*pv->GetObjectTranslation()+parentT;",
    ],
    "placed": ['std::cout<<"MESH "<<out.str()<<std::endl;meshes++;'],
    "placement": [
        "auto R=parentR*pv->GetObjectRotationValue();",
        "auto T=parentR*pv->GetObjectTranslation()+parentT;",
    ],
    "mesh": ["auto* poly=lv->GetSolid()->GetPolyhedron();"],
    "vertices": [
        "auto v=poly->GetVertex(i);G4ThreeVector local(v.x(),v.y(),v.z());",
        "auto world=R*local+T;",
    ],
    "faces": ["int n,ix[4];poly->GetFacet(i,n,ix);"],
    "parse": [
        "G4GDMLParser parser;parser.SetOverlapCheck(false);",
        "parser.Read(gdml.filename().string(),false);",
    ],
    "parsed": [
        "walk(parser.GetWorldVolume(),G4RotationMatrix(),G4ThreeVector(),parser);"
    ],
    "check": [
        "auto check=history.GetTopTransform().Inverse().TransformPoint(local);",
        "error=std::max(error,(world-check).mag());vertices++;",
    ],
}

# A plain ROOT macro. Only T and hit.det are needed; no remoll library.
ROOT_MACRO = r"""
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TLeaf.h>
#include <TSystem.h>
#include <cmath>
#include <iomanip>
#include <map>
#include <iostream>
#include <memory>
void count_hits(const char* path){try{
 std::unique_ptr<TFile> file(TFile::Open(path,"READ"));
 if(!file||file->IsZombie())throw std::runtime_error("Cannot open ROOT file");
 auto* tree=file->Get<TTree>("T");
 if(!tree)throw std::runtime_error("Missing tree T");
 if(!tree->GetBranch("hit.det")&&!tree->GetLeaf("hit.det"))throw std::runtime_error("Missing field hit.det");
 if(auto* b=tree->GetBranch("hit.det")){if(b->GetEntries()!=tree->GetEntries())throw std::runtime_error("Incomplete hit.det branch");}
 if(auto* l=tree->GetLeaf("hit.det")){
  std::string type=l->GetTypeName();
  if(type!="Int_t"&&type!="int"&&type!="vector<int>"&&type!="std::vector<int>")throw std::runtime_error("hit.det must contain integers");
 }
 if(auto* b=tree->GetBranch("hit.det")){
  if(!b->InheritsFrom("TBranchElement"))throw std::runtime_error("Unsupported leaf-list hit.det; use a split remoll hit branch or vector<int>");
 }
 TTreeReader reader(tree);TTreeReaderArray<int> det(reader,"hit.det");
 std::unique_ptr<TTreeReaderValue<double>> rate;
 std::string rateError;
 auto* leaf=tree->GetLeaf("rate");
 if(!leaf)rateError="Missing rate branch";
 else if(std::string(leaf->GetTypeName())!="Double_t"||leaf->GetLeafCount()||leaf->GetLenStatic()!=1)rateError="rate must be a scalar double";
 else if(leaf->GetBranch()->GetEntries()!=tree->GetEntries())rateError="Incomplete rate branch";
 else rate=std::make_unique<TTreeReaderValue<double>>(reader,"rate");
 std::map<int,long long> counts;long long events=0,records=0;
 std::map<int,double> rates;std::cout<<std::setprecision(17);
 auto report=[&](){
  std::cout<<"COUNTS {\"events\":"<<events<<",\"total_events\":"<<tree->GetEntries()<<",\"records\":"<<records<<",\"counts\":{";
  bool first=true;for(auto [id,n]:counts){if(!first)std::cout<<',';first=false;std::cout<<'"'<<id<<"\":"<<n;}
  std::cout<<"},\"rates\":";
  if(!rateError.empty())std::cout<<"null";
  else {std::cout<<'{';first=true;for(auto [id,v]:rates){if(!first)std::cout<<',';first=false;std::cout<<'"'<<id<<"\":"<<v;}std::cout<<'}';}
  std::cout<<",\"rate_error\":\""<<rateError<<"\"}"<<std::endl;
  std::cout<<"TRACE {\"stage\":\"root_count\",\"detail\":\"Count each stored hit once\"}"<<std::endl;
  std::string resume;if(!std::getline(std::cin,resume))throw std::runtime_error("Reader stopped");
 };
 while(reader.Next()){
  if(det.GetSetupStatus()<0)throw std::runtime_error("Cannot read hit.det as an integer array");
  double w=0;
  if(rate&&rateError.empty()){
   if(rate->GetSetupStatus()<0)rateError="Cannot read rate";
   else {w=**rate;if(!std::isfinite(w)||w<0)rateError="Nonfinite or negative event rate";}
  }
  for(int id:det){
   ++counts[id];++records;
   if(rateError.empty()){rates[id]+=w;if(!std::isfinite(rates[id]))rateError="Accumulated rate overflow";}
  }
  ++events;if(events<=3||events%10000==0)report();
 }
 if(events!=tree->GetEntries())throw std::runtime_error("ROOT reader stopped before all events were read");
 report();
 }catch(const std::exception& e){std::cerr<<"ERROR "<<e.what()<<std::endl;gSystem->Exit(1);}}
"""
SNIPPETS["root_count"] = (
    ROOT_MACRO.split("  for(int id:det){", 1)[1].split("  ++events", 1)[0].splitlines()
)
SNIPPETS["root_count"].insert(0, "for(int id:det){")


def code_location(code):
    """Resolve displayed statements to the actual source, including generated C++."""
    if not code:
        return dict(file="", lines=[])
    sources = [
        ("build/export_geometry.cc", CPP),
        ("build/count_hits.C", ROOT_MACRO),
        ("run.py", Path(__file__).read_text()),
    ]
    if code.startswith("find_package("):
        sources.insert(
            0,
            (
                "build/CMakeLists.txt",
                "cmake_minimum_required(VERSION 3.20)\nproject(geometry_walkthrough LANGUAGES CXX)\nset(CMAKE_CXX_STANDARD 17)\nfind_package(Geant4 REQUIRED gdml)\n",
            ),
        )
    for name, source in sources:
        lines = source.splitlines()
        numbers = []
        after = 0
        for statement in code.splitlines():
            if not statement.strip():
                numbers.append(None)
                continue
            matches = [
                i for i, line in enumerate(lines) if line.strip() == statement.strip()
            ]
            if not matches:
                break
            i = next((i for i in matches if i >= after), matches[0])
            numbers.append(i + 1)
            after = i + 1
        else:
            return dict(file=name, lines=numbers)
    return dict(file="", lines=[None] * len(code.splitlines()))


def choose_root_file():
    """Use the OS picker: return a local path without uploading/copying the file."""
    if sys.platform == "darwin":
        command = [
            "osascript",
            "-e",
            "try",
            "-e",
            'POSIX path of (choose file with prompt "Choose a remoll ROOT file")',
            "-e",
            "on error number -128",
            "-e",
            'return ""',
            "-e",
            "end try",
        ]
    elif os.name == "nt":
        command = [
            "powershell",
            "-NoProfile",
            "-STA",
            "-Command",
            'Add-Type -AssemblyName System.Windows.Forms; $d=New-Object System.Windows.Forms.OpenFileDialog; $d.Filter="ROOT files (*.root)|*.root"; if($d.ShowDialog() -eq "OK"){$d.FileName}',
        ]
    elif shutil.which("zenity"):
        command = [
            "zenity",
            "--file-selection",
            "--title=Choose a remoll ROOT file",
            "--file-filter=*.root",
        ]
    else:
        command = [
            sys.executable,
            "-c",
            'import tkinter as t;from tkinter.filedialog import askopenfilename;r=t.Tk();r.withdraw();print(askopenfilename(filetypes=[("ROOT files","*.root")]));r.destroy()',
        ]
    p = subprocess.run(command, capture_output=True, text=True)
    cancelled = command[0] == "zenity" and p.returncode == 1
    if p.returncode and not cancelled:
        message = (
            p.stderr.strip().splitlines()[-1]
            if p.stderr.strip()
            else "File picker could not open."
        )
        if command[0] == sys.executable:
            message = "Linux file picker needs Zenity or Python Tk. " + message
        raise ValueError(message)
    path = p.stdout.strip()
    if path and not Path(path).is_file():
        raise ValueError("The selected file is unavailable.")
    return path


lock = threading.RLock()
current = None
serial = 0


class Cancelled(Exception):
    pass


class Run:
    def __init__(self, ident, root, delay):
        self.id = ident
        self.root = root
        self.delay = delay
        self.cancel = threading.Event()
        self.condition = threading.Condition(lock)
        self.proc = None
        self.paused = False
        self.credits = 0
        self.fast = False
        self.data = dict(
            kind="geometry",
            stage="locate",
            detail="",
            code="",
            log=[],
            meshes=[],
            geometry=None,
            busy=True,
            revision=0,
            hits=None,
            placement=None,
        )
        self.thread = threading.Thread(target=self.work, daemon=True)

    def check(self):
        if self.cancel.is_set():
            raise Cancelled()

    def status(self, stage, detail="", code=""):
        with lock:
            self.check()
            self.data.update(
                stage=stage,
                detail=detail,
                code=code or "\n".join(SNIPPETS.get(stage, [])),
            )
            self.data["code_location"] = code_location(self.data["code"])
            self.data["revision"] += 1

    def gate(self):
        # The C++ reader blocks on stdin until this exact step is released.
        elapsed = 0.0
        previous = time.monotonic()
        with self.condition:
            while True:
                self.check()
                now = time.monotonic()
                if not self.paused:
                    elapsed += now - previous
                previous = now
                if self.credits:
                    self.credits -= 1
                    return
                if not self.paused and (self.fast or elapsed >= self.delay):
                    return
                self.condition.wait(0.025)

    def stop(self):
        with self.condition:
            self.cancel.set()
            self.condition.notify_all()
            proc = self.proc
            if self.data["busy"]:
                self.data.update(stage="stopping", detail="Stopping the reader…")
            else:
                self.data.update(
                    stage="stopped", detail="Stopped. Load again when ready."
                )
            if self.data["kind"] == "geometry":
                self.data["geometry"] = None
        if proc is not None:
            self.kill(proc)

    @staticmethod
    def kill(proc):
        if proc.poll() is not None:
            return
        try:
            if os.name == "posix":
                os.killpg(proc.pid, signal.SIGTERM)
            else:
                proc.terminate()
            try:
                proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                if os.name == "posix":
                    os.killpg(proc.pid, signal.SIGKILL)
                else:
                    proc.kill()
                proc.wait(timeout=2)
        except ProcessLookupError:
            pass

    def execute(self, args):
        with lock:
            self.check()
            proc = subprocess.Popen(
                [str(a) for a in args],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                start_new_session=os.name == "posix",
            )
            self.proc = proc
        try:
            for line in proc.stdout:
                self.check()
                line = line.rstrip()
                if line.startswith("TRACE "):
                    item = json.loads(line[6:])
                    with lock:
                        self.data["placement"] = item.get("placement")
                    self.status(item["stage"], item["detail"])
                    self.gate()
                    proc.stdin.write("\n")
                    proc.stdin.flush()
                elif line.startswith("MESH "):
                    with lock:
                        self.data["meshes"].append(json.loads(line[5:]))
                elif line.startswith("COUNTS "):
                    with lock:
                        self.data["hits"] = dict(
                            json.loads(line[7:]), complete=False, path=self.path
                        )
                else:
                    with lock:
                        self.data["log"] = (self.data["log"] + [line])[-30:]
            self.check()
            if proc.wait():
                errors = [s[6:] for s in self.data["log"] if s.startswith("ERROR ")]
                raise RuntimeError(
                    errors[-1] if errors else "\n".join(self.data["log"][-12:])
                )
        finally:
            try:
                self.kill(proc)
                for stream in (proc.stdout, proc.stdin):
                    try:
                        stream.close()
                    except (BrokenPipeError, OSError):
                        pass
            finally:
                with lock:
                    if self.proc is proc:
                        self.proc = None

    def work(self):
        try:
            root = Path(self.root).expanduser().resolve()
            gdml = root / "geometry/mollerMother.gdml"
            if not gdml.is_file():
                raise ValueError(
                    "Missing geometry/mollerMother.gdml in that directory."
                )
            self.status(
                "locate", str(gdml), 'gdml = root / "geometry/mollerMother.gdml"'
            )
            g4 = os.environ.get("Geant4_DIR")
            for cache in [
                root / "build/CMakeCache.txt",
                *root.glob("build*/CMakeCache.txt"),
            ]:
                if not g4 and cache.is_file():
                    match = re.search(
                        r"^Geant4_DIR:[^=]+=(.+)$", cache.read_text(), re.M
                    )
                    if match and Path(match[1], "Geant4Config.cmake").is_file():
                        g4 = match[1]
            if not g4:
                candidate = root.parent / "geant4-install/lib/cmake/Geant4"
                if (candidate / "Geant4Config.cmake").is_file():
                    g4 = str(candidate)
            build = HERE / "build"
            build.mkdir(exist_ok=True)
            stamp = hashlib.sha256((CPP + str(g4)).encode()).hexdigest()
            exe = build / (
                "bin/Release/export_geometry.exe"
                if os.name == "nt"
                else "bin/export_geometry"
            )
            if (
                not exe.is_file()
                or not (build / "stamp").is_file()
                or (build / "stamp").read_text() != stamp
            ):
                self.status(
                    "compile",
                    "Compiling the small Geant4 geometry reader.",
                    "find_package(Geant4 REQUIRED gdml)",
                )
                (build / "export_geometry.cc").write_text(CPP)
                (build / "CMakeLists.txt").write_text(
                    "cmake_minimum_required(VERSION 3.20)\nproject(geometry_walkthrough LANGUAGES CXX)\nset(CMAKE_CXX_STANDARD 17)\nfind_package(Geant4 REQUIRED gdml)\ninclude(${Geant4_USE_FILE})\nadd_executable(export_geometry export_geometry.cc)\ntarget_link_libraries(export_geometry ${Geant4_LIBRARIES})\n"
                )
                self.execute(
                    [
                        "cmake",
                        "-S",
                        build,
                        "-B",
                        build / "bin",
                        *(["-DGeant4_DIR=" + g4] if g4 else []),
                    ]
                )
                self.execute(
                    [
                        "cmake",
                        "--build",
                        build / "bin",
                        "--config",
                        "Release",
                        "-j",
                        "2",
                    ]
                )
                (build / "stamp").write_text(stamp)
            self.execute([exe, gdml])
            self.check()
            self.status(
                "package",
                "Assign labels and assemble the browser mesh.",
                'faces = [[vertices[j] for j in face] for face in v["faces"]]',
            )
            geometry = self.package(str(gdml))
            with lock:
                self.check()
                self.data["geometry"] = geometry
            self.status(
                "ready",
                f'{len(geometry["rows"])} channels · {len(geometry["meshes"])} quartz meshes',
            )
        except Cancelled:
            pass
        except Exception as error:
            if not self.cancel.is_set():
                self.status("error", str(error))
        finally:
            with lock:
                self.data["busy"] = False
                if self.cancel.is_set():
                    self.data.update(
                        stage="stopped",
                        detail="Stopped. Load again when ready.",
                        geometry=None,
                    )

    def package(self, gdml):
        rows = []
        meshes = []
        keys = {}
        for v in self.data["meshes"]:
            self.check()
            sm = re.fullmatch(r"logic_quartz_(\d+)_(\d+)", v["name"])
            ring = re.search(r"R([1-6])-", v["name"])
            if sm:
                module = int(sm[1])
                key = ("showermax", module)
                system = "showermax"
                detector = f"SM {module}"
                region = {1: "closed", 3: "open"}.get(module % 4, "transition")
            else:
                system = "ring" + ring[1]
                detector = v["detector"]
                key = (system, detector)
                if not detector:
                    raise ValueError("Missing DetNo: " + v["name"])
                seg = re.match(r"(FF|BF)Seg(\d+)", v["name"])
                region = (
                    "transition"
                    if seg[1] == "BF"
                    else "closed" if int(seg[2]) % 2 else "open"
                )
            if key not in keys:
                keys[key] = len(rows)
                rows.append(
                    dict(
                        index=len(rows),
                        detector=detector,
                        system=system,
                        region=region,
                        gdml_detectors=[],
                    )
                )
            i = keys[key]
            rows[i]["gdml_detectors"].append(v["detector"])
            vertices = v["vertices"]
            faces = [[vertices[j] for j in face] for face in v["faces"]]
            meshes.append(dict(index=i, center=v["center"], faces=faces))
        if {r["system"] for r in rows} != {
            *[f"ring{i}" for i in range(1, 7)],
            "showermax",
        }:
            raise ValueError(
                "This GDML does not contain all six quartz rings and ShowerMax."
            )
        order = sorted(range(len(rows)), key=lambda i: (rows[i]["system"], i))
        remap = {old: new for new, old in enumerate(order)}
        rows = [dict(rows[i], index=remap[i]) for i in order]
        for m in meshes:
            m["index"] = remap[m["index"]]
        parallel = Path(gdml).with_name("mollerParallel.gdml")
        if parallel.is_file():
            planes = {
                int(module): int(det)
                for module, det in re.findall(
                    r'<volume name="showerMaxStackVirtualPlane_(\d+)_log">.*?<auxiliary auxtype="DetNo" auxvalue="(\d+)"',
                    parallel.read_text(),
                    re.S,
                )
            }
            for row in rows:
                if row["system"] == "showermax":
                    row["plane_detector"] = planes.get(int(row["detector"].split()[1]))
        return {"rows": rows, "meshes": meshes, "source": gdml}


class HitRun(Run):
    def __init__(self, ident, previous, path, delay):
        super().__init__(ident, previous.root, delay)
        self.path = path
        self.data.update(
            kind="hits",
            geometry=previous.data["geometry"],
            stage="root_open",
            detail=path,
        )

    def work(self):
        try:
            path = Path(self.path).expanduser().resolve()
            if not path.is_file():
                raise ValueError("The selected ROOT file no longer exists.")
            root = shutil.which("root")
            if not root:
                raise ValueError(
                    "ROOT is not on PATH. Activate your ROOT installation and relaunch."
                )
            macro = HERE / "build/count_hits.C"
            macro.parent.mkdir(exist_ok=True)
            macro.write_text(ROOT_MACRO)
            self.status(
                "root_open",
                str(path),
                'std::unique_ptr<TFile> file(TFile::Open(path,"READ"));',
            )
            self.execute(
                [root, "-l", "-b", "-q", str(macro) + "(" + json.dumps(str(path)) + ")"]
            )
            self.check()
            with lock:
                if self.data["hits"] is None:
                    raise ValueError("ROOT produced no hit counts.")
                self.data["hits"]["complete"] = True
            self.status(
                "ready", "Hit counts complete", "\n".join(SNIPPETS["root_count"])
            )
        except Cancelled:
            pass
        except Exception as error:
            if not self.cancel.is_set():
                self.status("error", str(error))
        finally:
            with lock:
                self.data["busy"] = False
                if self.cancel.is_set():
                    self.data.update(
                        stage="stopped", detail="Analysis stopped. Counts are partial."
                    )


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def send(self, body, kind="application/json", code=200):
        data = body.encode()
        self.send_response(code)
        self.send_header("Content-Type", kind)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        try:
            self.wfile.write(data)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def do_GET(self):
        if self.path == "/":
            self.send((HERE / "index.html").read_text(), "text/html; charset=utf-8")
            return
        with lock:
            job = current
            if self.path == "/status":
                d = (
                    job.data
                    if job
                    else dict(
                        stage="choose",
                        detail="",
                        code="",
                        log=[],
                        busy=False,
                        revision=0,
                    )
                )
                result = {k: v for k, v in d.items() if k not in ("geometry", "meshes")}
                result.update(
                    app="geometry_walkthrough_v2",
                    run=job.id if job else 0,
                    mesh_count=len(d.get("meshes", [])),
                    paused=job.paused if job else False,
                    delay=job.delay if job else 0.8,
                    default_root=DEFAULT_REMOLL,
                    root=job.root if job else "",
                    child_pid=job.proc.pid if job and job.proc else None,
                )
            elif self.path in ("/geometry", "/preview"):
                result = {
                    "run": job.id if job else 0,
                    "data": (
                        job.data["geometry" if self.path == "/geometry" else "meshes"]
                        if job
                        else None
                    ),
                }
            else:
                self.send("{}", code=404)
                return
            body = json.dumps(result)
        self.send(body)

    def do_POST(self):
        global current, serial
        if self.headers.get("Origin") not in (
            None,
            "http://" + self.headers.get("Host", ""),
        ):
            self.send("{}", code=403)
            return
        try:
            size = int(self.headers.get("Content-Length", "0"))
            if size < 0 or size > 8192:
                raise ValueError("Request too large")
            args = json.loads(self.rfile.read(size) or b"{}")
            if self.path == "/choose-root":
                self.send(json.dumps({"path": choose_root_file()}))
                return
            if self.path == "/analyze":
                path = args.get("path")
                delay = float(args.get("delay", 0.8))
                if not isinstance(path, str) or not path or not 0.01 <= delay <= 5:
                    raise ValueError("Choose a ROOT file first.")
                with lock:
                    if not current or not current.data["geometry"]:
                        raise ValueError("Load the geometry first.")
                    if current.data["busy"]:
                        self.send(
                            '{"error":"Stop the current operation first."}', code=409
                        )
                        return
                    serial += 1
                    current = HitRun(serial, current, path, delay)
                    current.thread.start()
                    ident = current.id
                self.send(json.dumps({"run": ident}))
                return
            if self.path == "/load":
                root = args.get("root")
                if not isinstance(root, str) or not root.strip():
                    raise ValueError("Enter the remoll directory.")
                delay = float(args.get("delay", 0.8))
                if not 0.01 <= delay <= 5:
                    raise ValueError("Invalid step duration.")
                with lock:
                    if current and current.data["busy"]:
                        self.send('{"error":"Stop the current load first."}', code=409)
                        return
                    serial += 1
                    current = Run(serial, root, delay)
                    current.thread.start()
                    ident = current.id
                self.send(json.dumps({"run": ident}))
                return
            with lock:
                job = current
            if self.path == "/quit":
                if job:
                    job.stop()
                self.send("{}")
                threading.Thread(target=self.server.shutdown, daemon=True).start()
                return
            if self.path == "/stop":
                if job:
                    job.stop()
                self.send("{}")
                return
            if self.path not in ("/pause", "/step", "/pace", "/finish"):
                self.send("{}", code=404)
                return
            if job:
                with job.condition:
                    if self.path == "/pause":
                        job.paused = not job.paused
                    elif self.path == "/step":
                        job.paused = True
                        job.credits += 1
                    elif self.path == "/finish":
                        job.fast = True
                        job.paused = False
                    else:
                        delay = float(args["delay"])
                        if not 0.01 <= delay <= 5:
                            raise ValueError("Invalid step duration.")
                        job.delay = delay
                        job.fast = False
                    job.condition.notify_all()
            self.send("{}")
        except (ValueError, TypeError, KeyError, OSError) as e:
            self.send(json.dumps({"error": str(e)}), code=400)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8785)
    ap.add_argument("--no-open", action="store_true")
    args = ap.parse_args()
    try:
        server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    except OSError:
        try:
            with urlopen(f"http://127.0.0.1:{args.port}/status", timeout=2) as response:
                existing = json.load(response)
            if existing.get("app") != "geometry_walkthrough_v2":
                raise ValueError("Different app")
        except Exception:
            raise SystemExit(
                f"Port {args.port} is unavailable. Use --port with another port."
            )
        print(f"Already running: http://127.0.0.1:{args.port}/", flush=True)
        if not args.no_open:
            webbrowser.open(f"http://127.0.0.1:{args.port}/")
        sys.exit(0)
    print(
        f"Open http://127.0.0.1:{server.server_port}/ · Stop cancels loading; Quit or Ctrl+C closes the app.",
        flush=True,
    )
    if not args.no_open:
        webbrowser.open(f"http://127.0.0.1:{server.server_port}/")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        if current:
            current.stop()
            current.thread.join(timeout=4)
        server.server_close()
