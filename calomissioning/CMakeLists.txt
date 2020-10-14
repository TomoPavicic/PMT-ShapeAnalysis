# - Example project for accessing SuperNEMO RTD data (Raw Trigger Data)

cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(snfee_read_rtd LANGUAGES CXX)

# - Load Custom Modules
list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)
message(STATUS "[info] CMAKE_MODULE_PATH = '${CMAKE_MODULE_PATH}'")

find_package(SNFrontEndElectronics)
message(STATUS  "SNFrontEndElectronics_INCLUDE_DIR=${SNFrontEndElectronics_INCLUDE_DIR}")
include_directories(${SNFrontEndElectronics_INCLUDE_DIRS})

# - Executable:
add_executable(read_rtd read_rtd.cxx)

target_link_libraries(read_rtd SNFrontEndElectronics::snfee)

set_target_properties(read_rtd
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )

install(FILES
    ${CMAKE_BINARY_DIR}/snfee-tutorial_lv0
    DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
    )


# end
