cmake_minimum_required(VERSION 3.23)
project(labview_thermodynamics LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES
    src/impl/lvt.cc
    src/impl/backend/lvt_naive.cc
)

add_library(labview_thermodynamics SHARED ${SOURCES})
target_include_directories(labview_thermodynamics PUBLIC src PRIVATE src/impl)
target_compile_definitions(labview_thermodynamics PRIVATE LVT_BUILD)

add_executable(test_rod_1d tests/test_rod_1d.cc)
target_link_libraries(test_rod_1d PRIVATE labview_thermodynamics)
target_include_directories(test_rod_1d PRIVATE src)

add_executable(test_plate_2d tests/test_plate_2d.cc)
target_link_libraries(test_plate_2d PRIVATE labview_thermodynamics)
target_include_directories(test_plate_2d PRIVATE src)

add_executable(test_plate_2d_radiation tests/test_plate_2d_radiation.cc)
target_link_libraries(test_plate_2d_radiation PRIVATE labview_thermodynamics)
target_include_directories(test_plate_2d_radiation PRIVATE src)
