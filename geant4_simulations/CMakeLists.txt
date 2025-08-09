cmake_minimum_required(VERSION 3.16)
project(MuonThroughRock)

find_package(Geant4 REQUIRED ui_all vis_all)

# Header path
include_directories(${PROJECT_SOURCE_DIR}/include)

set(SRCS
    src/main.cc
    src/DetectorConstruction.cc
    src/PrimaryGeneratorAction.cc
    src/RunAction.cc
    src/SteppingAction.cc
)

add_executable(MuonThroughRock ${SRCS})
target_link_libraries(MuonThroughRock ${Geant4_LIBRARIES})

# Nice out-of-tree build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

