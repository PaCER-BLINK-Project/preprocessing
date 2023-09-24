cmake_minimum_required(VERSION 3.10.0)
project(blink-preprocessing)
set(CMAKE_CXX_STANDARD 17)

include(CTest)

option(USE_OPENMP "Compile with OpenMP enabled." ON)
option(USE_HIP "Compile the HIP implementation." OFF)

if(USE_HIP AND (NOT CMAKE_CXX_COMPILER MATCHES "hipcc"))
    message(FATAL_ERROR "HIP code must be compiled with hipcc.")
endif()

# External dependencies
set(BLINK_TEST_DATADIR $ENV{BLINK_TEST_DATADIR})
if(NOT BLINK_TEST_DATADIR)
message(FATAL_ERROR "The 'BLINK_TEST_DATADIR' environment variable was not set and is needed for testing.")
endif()

find_library(ASTROIO_LIB blink_astroio HINTS ENV LD_LIBRARY_PATH)
if( NOT ASTROIO_LIB )
    message(FATAL_ERROR "ASTROIO not found.")
endif()

file(GLOB preprocessing_sources "src/*.cpp")
file(GLOB preprocessing_headers "src/*.hpp")

if(USE_HIP)
	add_definitions("-DUSE_HIP")
	file(GLOB hip_sources "src/hip/*.cpp")
	file(GLOB hip_headers "src/hip/*.hpp")
	list(APPEND preprocessing_sources ${hip_sources})
	list(APPEND preprocessing_headers ${hip_headers})
endif()

add_library(blink_preprocessing SHARED ${preprocessing_sources})
set_target_properties(blink_preprocessing PROPERTIES PUBLIC_HEADER "${preprocessing_headers}")
target_link_libraries(blink_preprocessing ${ASTROIO_LIB})


install(TARGETS blink_preprocessing
    LIBRARY DESTINATION "lib"
    PUBLIC_HEADER DESTINATION "include"
)


# TESTS
add_executable(mapping_test tests/mapping_test.cpp)
target_link_libraries(mapping_test blink_preprocessing)
add_test(NAME mapping_test COMMAND mapping_test)

add_executable(calibration_test tests/calibration_test.cpp)
target_link_libraries(calibration_test blink_preprocessing)
add_test(NAME calibration_test COMMAND calibration_test)
