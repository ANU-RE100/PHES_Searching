set(UTIL_SOURCES
    reservoir.cpp
    polygons.cpp
    csv.cpp
    kml.cpp
    coordinates.cpp
    phes_base.cpp
    variable_parser.cpp
    constructor_helpers.cpp
    model2D.cpp
    search_config.cpp
    fsolve.cpp)

include_directories(${MPI_CXX_INCLUDE_PATH} ${JSON_INCLUDE_PATH})

add_library(screening_objects OBJECT screening.cpp)
add_library(pairing_objects OBJECT pairing.cpp)
add_library(pretty_set_objects OBJECT pretty_set.cpp)
add_library(constructor_objects OBJECT constructor.cpp)
add_library(search_driver_objects OBJECT search_driver.cpp)
add_library(shapefile_tiling_objects OBJECT shapefile_tiling.cpp)
add_library(reservoir_constructor_objects OBJECT reservoir_constructor.cpp)
add_library(depression_volume_finding_objects OBJECT depression_volume_finding.cpp)
add_library(util_objects OBJECT ${UTIL_SOURCES})
