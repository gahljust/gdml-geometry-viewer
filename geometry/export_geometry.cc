#include "../root/json.hh"
// Geometry slice: Geant4 reads GDML; export placed detector surfaces in mm.
// Adapted from showermax_live/export_geometry.cc and geometry_walkthrough.
#include "G4GDMLParser.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Polyhedron.hh"
#include "G4VSolid.hh"
#include "G4NavigationHistory.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>

G4NavigationHistory history;
double maxError = 0;
int meshCount = 0, vertexCount = 0;
bool parallelPlanes = false;

void point(std::ostream& out, const G4ThreeVector& p) {
    out << '[' << p.x() << ',' << p.y() << ',' << p.z() << ']';
}

// Display selection only: six quartz rings and ShowerMax quartz geometry.
// No hit selection or signal calculation occurs in this program.
std::string systemFor(const std::string& name) {
    static const std::regex main("R([1-6])-.*Quartz");
    static const std::regex showermax("logic_quartz_[0-9]+_[0-9]+");
    std::smatch match;
    if (std::regex_search(name, match, main)) return "Ring " + match[1].str();
    if (std::regex_match(name, showermax)) return "ShowerMax";
    return "";
}

void walk(G4VPhysicalVolume* placed, const G4RotationMatrix& parentRotation,
          const G4ThreeVector& parentPosition, G4GDMLParser& parser,
          std::ostream& out, int depth = 0) {
    if (depth == 0) history.SetFirstEntry(placed);
    else history.NewLevel(placed, kNormal, placed->GetCopyNo());

    // Compose each placement with its ancestors, including nested modules.
    const auto rotation = parentRotation * placed->GetObjectRotationValue();
    const auto position = parentRotation * placed->GetObjectTranslation() + parentPosition;
    auto* volume = placed->GetLogicalVolume();
    const std::string name = volume->GetName();
    int detector = 0;
    for (const auto& aux : parser.GetVolumeAuxiliaryInformation(volume))
        if (aux.type == "DetNo") detector = std::stoi(aux.value);
    const auto system = parallelPlanes ? (detector>=73001 && detector<=73028 ? std::string("ShowerMax") : std::string()) : systemFor(name);
    if (!system.empty()) {
        auto* mesh = volume->GetSolid()->GetPolyhedron();
        if (!mesh) throw std::runtime_error("Cannot mesh " + name);
        if (meshCount++) out << ',';
        out << "{\"name\":" << viewer_json::quote(name)
            << ",\"system\":" << viewer_json::quote(system)
            << ",\"geometry_detector\":" << detector << ",\"hit_surface\":" << (parallelPlanes || system!="ShowerMax" ? "true" : "false") << ",\"vertices\":[";
        for (int i = 1; i <= mesh->GetNoVertices(); ++i) {
            const auto vertex = mesh->GetVertex(i);
            const G4ThreeVector local(vertex.x(), vertex.y(), vertex.z());
            const auto world = rotation * local + position;
            // Cross-check against Geant4's own navigation transform.
            const auto check = history.GetTopTransform().Inverse().TransformPoint(local);
            maxError = std::max(maxError, (world - check).mag());
            ++vertexCount;
            if (i > 1) out << ',';
            point(out, world);
        }
        out << "],\"faces\":[";
        for (int i = 1; i <= mesh->GetNoFacets(); ++i) {
            int count, indices[4];
            mesh->GetFacet(i, count, indices);
            if (i > 1) out << ',';
            out << '[';
            for (int j = 0; j < count; ++j) {
                if (j) out << ',';
                out << indices[j] - 1;
            }
            out << ']';
        }
        out << "]}";
    }
    for (int i = 0; i < volume->GetNoDaughters(); ++i)
        walk(volume->GetDaughter(i), rotation, position, parser, out, depth + 1);
    if (depth > 0) history.BackLevel();
}

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::runtime_error("Usage: export_geometry mollerMother.gdml geometry.json");
        G4UnitDefinition::BuildUnitsTable();
        new G4UnitDefinition("inch", "in", "Length", 2.54 * cm);
        new G4UnitDefinition("mil", "mil", "Length", .001 * 2.54 * cm);
        const auto gdml = std::filesystem::absolute(argv[1]);
        const auto output = std::filesystem::absolute(argv[2]);
        std::filesystem::current_path(gdml.parent_path());
        G4GDMLParser parser;
        parser.SetOverlapCheck(false);
        parser.Read(gdml.filename().string(), false);
        std::ofstream out(output);
        if (!out) throw std::runtime_error("Cannot write geometry JSON");
        out << std::setprecision(12) << "{\"units\":\"mm\",\"meshes\":[";
        walk(parser.GetWorldVolume(), G4RotationMatrix(), G4ThreeVector(), parser, out);
        const auto parallel = gdml.parent_path() / "mollerParallel.gdml";
        if (std::filesystem::exists(parallel)) {
            parallelPlanes = true;
            G4GDMLParser planeParser;
            planeParser.SetOverlapCheck(false);
            planeParser.Read(parallel.string(), false);
            walk(planeParser.GetWorldVolume(), G4RotationMatrix(), G4ThreeVector(), planeParser, out);
        }
        if (!meshCount) throw std::runtime_error("No supported remoll detector geometry found");
        if (maxError > 1e-7) throw std::runtime_error("Placement transform check failed");
        out << "],\"mesh_count\":" << meshCount << ",\"vertex_count\":" << vertexCount
            << ",\"max_transform_error_mm\":" << maxError << "}\n";
        out.close();
        if (!out) throw std::runtime_error("Geometry JSON write failed");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
