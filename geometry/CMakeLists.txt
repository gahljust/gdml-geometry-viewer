cmake_minimum_required(VERSION 3.20)
project(detector_geometry LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Geant4 REQUIRED gdml)
if(Geant4_CXX_STANDARD GREATER CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD ${Geant4_CXX_STANDARD})
endif()
add_executable(export_geometry export_geometry.cc)
target_link_libraries(export_geometry PRIVATE ${Geant4_LIBRARIES})
# Same location for Make/Ninja and multi-configuration Visual Studio/Xcode builds.
set_target_properties(export_geometry PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/bin")
