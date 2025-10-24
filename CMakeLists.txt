cmake_minimum_required(VERSION 3.10.0)
project(blink-preprocessing)
set(CMAKE_CXX_STANDARD 17)
include(CTest)

option(USE_OPENMP "Compile with OpenMP enabled." ON)
option(USE_CUDA "Compile with NVIDIA CUDA support." OFF)
option(USE_HIP "Compile with AMD HIP support." OFF)

if(USE_CUDA)
enable_language(CUDA CXX)
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

file(GLOB sources 
    "src/calibration.cpp"
    "src/mapping.cpp"
)

file(GLOB headers
    "src/calibration.hpp"
    "src/mapping.hpp"
)

file(GLOB apps "apps/*.cpp")
file(GLOB tests "tests/*.cpp")


if(CMAKE_CXX_COMPILER MATCHES "hipcc" OR USE_CUDA)
	file(GLOB gpu_sources
        "src/mapping_gpu.cpp"
        "src/calibration_gpu.cpp"
    )
	file(GLOB gpu_headers 
        "src/mapping_gpu.hpp"
        "src/calibration_gpu.hpp"
    )
	list(APPEND sources ${gpu_sources})
	list(APPEND headers ${gpu_headers})
endif()

if(USE_CUDA)
set_source_files_properties( ${sources} ${headers} ${apps} ${tests} PROPERTIES LANGUAGE CUDA)
add_compile_options("--expt-relaxed-constexpr")
add_compile_options("-G")
endif()

add_library(blink_preprocessing SHARED ${sources})
set_target_properties(blink_preprocessing PROPERTIES PUBLIC_HEADER "${headers}")
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
