cmake_minimum_required(VERSION 3.12)
project(PHES_Searching VERSION 1.0.0)

message("-- Building RE100 PHES Searching software --")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message("-- CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

enable_language(CXX)
set(CMAKE_CXX_STANDARD 20)

set(PHES_SOURCE_PATH "${PROJECT_SOURCE_DIR}/src")

# ##################################################################################################
# setup general C++ compilation flags
# ##################################################################################################
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -g")

find_package(GDAL REQUIRED)
find_package(
  Boost
  COMPONENTS system filesystem
  REQUIRED)

# if("${MPI_ROOT}" STREQUAL "") message("-- Searching for MPI") find_package(MPI REQUIRED)
# elseif(EXISTS ${MPI_ROOT}) if(EXISTS ${MPI_ROOT}/include64) set(MPI_CXX_INCLUDE_PATH
# ${MPI_ROOT}/include64) else() set(MPI_CXX_INCLUDE_PATH ${MPI_ROOT}/include) endif() if(EXISTS
# ${MPI_ROOT}/lib64) set(MPI_CXX_LIBRARIES ${MPI_ROOT}/lib64) else() set(MPI_CXX_LIBRARIES
# ${MPI_ROOT}/lib) endif() else() message( FATAL_ERROR "MPI not found! Please use the -DMPI_ROOT
# flag to select an MPI for installation.") endif()
#
# if(DEFINED MPI_CXX_INCLUDE_PATH) option(HAVE_MPI "enable mpi" ON) message("-- Found MPI")
# message("   ** MPI include directory: ${MPI_CXX_INCLUDE_PATH}") endif()

# ##################################################################################################
# set JSON
# ##################################################################################################

# if("${JSON_ROOT}" STREQUAL "") message("-- Searching for JSON parser") find_package(nlohmann_json
# 3.2.0 REQUIRED) elseif(EXISTS ${JSON_ROOT}) message("-- Using ${JSON_ROOT}/include as the JSON
# install --") set(JSON_INCLUDE_PATH ${JSON_ROOT}/include) else() message(FATAL_ERROR "JSON parser
# not found!" "Please use the -DJSON_ROOT flag to select a JSON parser for installation.") endif()

# list(LENGTH MPI_CXX_LIBRARIES NUM_MPI_LIBS) if(NUM_MPI_LIBS GREATER 1) list(JOIN MPI_CXX_LIBRARIES
# " -L" MPI_LIBRARIES) set(MPI_LIBRARIES "-L${MPI_LIBRARIES} -lmpi") if(${MPI_LIBRARIES} MATCHES
# ".*libmpi_cxx.so.*") set(MPI_LIBRARIES "${MPI_LIBRARIES} -lmpi_cxx") endif() elseif(EXISTS
# "${MPI_CXX_LIBRARIES}/libmpi_ibm.so") set(MPI_LIBRARIES "${MPI_LIBRARIES}-lmpi_ibm") else()
# set(MPI_LIBRARIES "-L${MPI_CXX_LIBRARIES} -lmpi") endif() if(EXISTS
# "${MPI_CXX_LIBRARIES}/libmpi_cxx.so") set(MPI_LIBRARIES "${MPI_LIBRARIES} -lmpi_cxx") endif()

add_subdirectory(src)

foreach(target screening pairing pretty_set constructor search_driver shapefile_tiling
    reservoir_constructor depression_volume_finding gpkg_tiling)
  set(TARGETS $<TARGET_OBJECTS:util_objects> $<TARGET_OBJECTS:${target}_objects>)
  add_executable(${target} ${TARGETS})

  target_include_directories(
    ${target}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC SYSTEM ${GDAL_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH} ${JSON_INCLUDE_PATH})

  target_link_libraries(${target} PUBLIC ${MPI_LIBRARIES} ${GDAL_LIBRARY} shp
                                         ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
  set_target_properties(${target} PROPERTIES LINKER_LANGUAGE CXX)

  install(TARGETS ${target} DESTINATION ${PROJECT_SOURCE_DIR}/bin)

endforeach()
