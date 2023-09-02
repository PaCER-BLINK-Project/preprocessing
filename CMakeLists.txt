cmake_minimum_required(VERSION 3.10.0)
project(blink-preprocessing)
set(CMAKE_CXX_STANDARD 17)

include(CTest)


# External dependencies
find_library(ASTROIO_LIB blink_astroio HINTS ENV LD_LIBRARY_PATH)
if( NOT ASTROIO_LIB )
    message(FATAL_ERROR "ASTROIO not found.")
endif()

file(GLOB preprocessing_sources "src/*.cpp")
file(GLOB preprocessing_headers "src/*.hpp")

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
