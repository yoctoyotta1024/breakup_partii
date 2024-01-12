# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
  # cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler 
endif()

# set project name and source directory of .cpp files
project("buii")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")
message("CLEO LIBS SOURCE DIR: ${CLEOLIBS_SOURCE_DIR}")

# Kokkos library stuff
set(CMAKE_CXX_FLAGS_B4 "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_STANDARD_B4 "${CMAKE_CXX_STANDARD}")
include(FetchContent)
FetchContent_Declare(
  Kokkos
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/kokkos/kokkos/archive/refs/tags/4.2.00.tar.gz
  GIT_TAG 71a9bcae52543bd065522bf3e41b5bfa467d8015
)
set(CMAKE_CXX_STANDARD "20")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
FetchContent_MakeAvailable(Kokkos)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_B4}")
set(CMAKE_CXX_STANDARD "${CMAKE_CXX_STANDARD_B4}")

# Set libraries and compiler warnings to use
set(SDMLIBS gridboxes initialise observers runcleo superdrops zarr)

# set(compileops -Wall -pedantic -g -gdwarf-4) # less strict without -Werror
# set compiler options
if(NOT DEFINED compileops)
  set(compileops ${defaultcompileops})
endif()

# primary executable for CLEO sdm coalbure
set(coalbure "buii_coalbure")
add_executable(${coalbure} "main_coalbure.cpp")
target_link_libraries(${coalbure} PRIVATE coupldyn_fromfile cartesiandomain "${SDMLIBS}")
target_link_libraries(${coalbure} PUBLIC Kokkos::kokkos)
target_include_directories(${coalbure} PRIVATE "${CLEOLIBS_SOURCE_DIR}")
target_compile_options(${coalbure} PRIVATE ${compileops})

set(coalonly "buii_coalonly")
add_executable(${coalonly} "main_coalonly.cpp")
target_link_libraries(${coalonly} PRIVATE coupldyn_fromfile cartesiandomain "${SDMLIBS}")
target_link_libraries(${coalonly} PUBLIC Kokkos::kokkos)
target_include_directories(${coalonly} PRIVATE "${CLEOLIBS_SOURCE_DIR}")
target_compile_options(${coalonly} PRIVATE ${compileops})

# set compiler properties for target(s)
set_target_properties(${coalbure} ${coalonly} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)