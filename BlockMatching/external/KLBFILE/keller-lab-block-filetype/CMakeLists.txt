project(klb)
cmake_minimum_required(VERSION 3.0)

# check git version (when cloning from bitbucket.org/fernandoamat/keller-lab-block-filetype.git)
#
# include(cmake/git-versioning.cmake)

# for export through KLBUse.cmake
set(KLB_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_INSTALL_PREFIX}/include/nifti)
set(KLB_LIBRARY_DIRS 
    ${PROJECT_BINARY_DIR}/src
    ${PROJECT_BINARY_DIR}/src/external/zlib-1.2.8
    ${PROJECT_BINARY_DIR}/src/external/bzip2-1.0.6
    ${CMAKE_INSTALL_PREFIX}/lib)

set(KLB_USE_FILE     ${PROJECT_BINARY_DIR}/KLBUse.cmake)
set(KLB_CMAKE_DIR    ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/KLBConfig.cmake.in)
   configure_file( ## Build tree configure file
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/KLBConfig.cmake.in
    ${PROJECT_BINARY_DIR}/KLBConfig.cmake
    @ONLY IMMEDIATE)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/KLBUse.cmake.in)
   configure_file( ## Common use file
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/KLBUse.cmake.in
    ${PROJECT_BINARY_DIR}/KLBUse.cmake
    @ONLY IMMEDIATE)
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_C_STANDARD   99)
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(src)

option(build_tests "Build tests" OFF)
if (build_tests)
    enable_testing()
    add_subdirectory(test)
endif()

include(InstallRequiredSystemLibraries)
include(CPack)
