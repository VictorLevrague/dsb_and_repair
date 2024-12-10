#include <map>
#include "G4Types.hh"
#include <array>

std::map<std::array<G4int,3>,G4long> g_map_voxel; // map de coordonnées spatiale pour le voxel (Mathieu)
std::map<G4long, std::string> g_map_voxel2; //G4String: voxel_type
