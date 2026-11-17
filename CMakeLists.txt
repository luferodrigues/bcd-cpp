cmake_minimum_required(VERSION 3.16)
project(patchy_bcd LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(NATIVE "Optimize for this machine's CPU (-march=native)" ON)

# Everything except main(): shared by the simulation program and the test tool.
add_library(bcd_core STATIC
  src/params.cpp
  src/cells.cpp
  src/system.cpp
  src/clusters.cpp
  src/steps.cpp
  src/transitions.cpp
  src/output.cpp
)
target_include_directories(bcd_core PUBLIC src)

# PUBLIC: these settings also apply to every program that links bcd_core.
target_compile_options(bcd_core PUBLIC -Wall -Wextra -Wpedantic)

# Never merge a*b+c into one fused instruction: keeps arithmetic identical to NumPy's, so
# deterministic results can be compared with the Python code digit for digit.
target_compile_options(bcd_core PUBLIC -ffp-contract=off)

# Release: full optimization for production runs.
target_compile_options(bcd_core PUBLIC $<$<CONFIG:Release>:-O3>)
if(NATIVE)
  target_compile_options(bcd_core PUBLIC $<$<CONFIG:Release>:-march=native>)
endif()

# Debug: address + undefined-behaviour sanitizers catch out-of-bounds and similar bugs.
target_compile_options(bcd_core PUBLIC
  $<$<CONFIG:Debug>:-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined>)
target_link_options(bcd_core PUBLIC
  $<$<CONFIG:Debug>:-fsanitize=address,undefined>)

# The simulation program.
add_executable(patchy_bcd src/main.cpp)
target_link_libraries(patchy_bcd PRIVATE bcd_core)

# Small helper used by tools/test_numerics.py to compare C++ functions with Python/SciPy.
add_executable(check_tool tests/check_tool.cpp)
target_link_libraries(check_tool PRIVATE bcd_core)
