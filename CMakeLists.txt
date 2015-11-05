cmake_minimum_required(VERSION 2.8)
project(GDSCAD CXX)
set(GDSCAD_MAJOR_VERSION 0)
set(GDSCAD_MINOR_VERSION 1)
set(GDSCAD_PATCH_VERSION 0)

# We use delegated constructors so we must use c++11 or greater
add_definitions("-std=c++11")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/GDSCADConfig.h.in"
  "${PROJECT_BINARY_DIR}/GDSCADConfig.h"
  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

add_library(gdscad SHARED containers.cpp shapes.cpp gdsfile.cpp)
install(TARGETS gdscad
  LIBRARY DESTINATION lib)