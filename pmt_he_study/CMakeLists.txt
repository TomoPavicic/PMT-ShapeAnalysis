cmake_minimum_required(VERSION 3.8)
project(PMT_Study VERSION 0.1.0)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wfatal-errors -Wshadow -pedantic")

find_package(Boost REQUIRED program_options)

find_package(ROOT REQUIRED Tree)

set_source_files_properties(xml_reader.cpp PROPERTIES LANGUAGE CXX)
set_source_files_properties(pmt_charge_analysis.cpp PROPERTIES LANGUAGE CXX)

add_executable(xml_reader xml_reader.cpp)
add_executable(pmt_charge_analysis pmt_charge_analysis.cpp)

target_link_libraries(xml_reader PRIVATE Boost::program_options)
target_link_libraries(pmt_charge_analysis PRIVATE Boost::program_options)

target_link_libraries(xml_reader PRIVATE ROOT::Core)
target_link_libraries(pmt_charge_analysis PRIVATE ROOT::Core)

target_link_libraries(xml_reader PRIVATE ROOT::Tree)
target_link_libraries(pmt_charge_analysis PRIVATE ROOT::Tree)

target_link_libraries(xml_reader PRIVATE ROOT::MathCore)
target_link_libraries(pmt_charge_analysis PRIVATE ROOT::MathCore)

target_link_libraries(xml_reader PRIVATE ROOT::Hist)
target_link_libraries(pmt_charge_analysis PRIVATE ROOT::Hist)

target_link_libraries(xml_reader PRIVATE ROOT::Gpad)
target_link_libraries(pmt_charge_analysis PRIVATE ROOT::Gpad)


target_compile_features(xml_reader
        PRIVATE cxx_auto_type cxx_range_for cxx_uniform_initialization
        )

target_compile_features(pmt_charge_analysis
        PRIVATE cxx_auto_type cxx_range_for cxx_uniform_initialization
        )

