### vtPack.cmake --- 
## 
## Author: Gregoire Malandain
## Copyright (C) 2013 - Gregoire Malandain, Inria.
## Version: $Id$
##
##
##
##
######################################################################
## 
### Commentary: 
## 
######################################################################
## 
### Change log:
## 
######################################################################

include (InstallRequiredSystemLibraries)

## #################################################################
## Global settings
## #################################################################

set(CPACK_PACKAGE_NAME ${PROJECT_NAME})

#set(CPACK_COMPONENTS_ALL apps devel libs)

if("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
  execute_process(COMMAND lsb_release -irs
    COMMAND sed "s/ //"
    COMMAND sed "s/Fedora/fc/"
    COMMAND tr -d '\n'
    OUTPUT_VARIABLE DISTRIB
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND arch 
    OUTPUT_VARIABLE ARCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-${${PROJECT_NAME}_VERSION}-1.${DISTRIB}-${ARCH}")

  set(CPACK_PACKAGING_INSTALL_PREFIX  "usr" CACHE PATH "Where you want to install your package")
  mark_as_advanced(CPACK_PACKAGING_INSTALL_PREFIX)

else("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
  set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-${${PROJECT_NAME}_VERSION}.${CMAKE_SYSTEM_PROCESSOR}")
  set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")

set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PROJECT_NAME}-${${PROJECT_NAME}_VERSION}-src")

set(CPACK_PACKAGE_VENDOR "INRIA")
set(CPACK_PACKAGE_CONTACT "Gregoire Malandain <gregoire.malandain@inria.fr>")
# set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.txt)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_NAME})
set(CPACK_PACKAGE_VERSION_MAJOR ${${PROJECT_NAME}_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${${PROJECT_NAME}_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${${PROJECT_NAME}_VERSION_BUILD})
#set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/COPYING.txt)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")

## #################################################################
## Generator list
## #################################################################

set(CPACK_BINARY_STGZ OFF)
set(CPACK_BINARY_TBZ2 OFF)
set(CPACK_BINARY_TGZ OFF)
set(CPACK_SOURCE_TBZ2 OFF)
set(CPACK_SOURCE_TGZ OFF)
set(CPACK_SOURCE_TZ OFF)
set(CPACK_SOURCE_ZIP OFF)

if(APPLE AND NOT UNIX)
  set(CPACK_GENERATOR "PackageMaker")
   set(CPACK_BUNDLE_NAME ${PROJECT_NAME})
endif(APPLE AND NOT UNIX)

if(WIN32)
  set(CPACK_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
  set(CPACK_PACKAGING_INSTALL_PREFIX "")
  set(CPACK_BINARY_NSIS ON CACHE BOOL "Enable to build NSIS packages")
  set(CPACK_PACKAGE_EXECUTABLES "vt-platform" "STracking")
  set(CPACK_NSIS_MODIFY_PATH ON)
  set(CPACK_CREATE_DESKTOP_LINKS "STracking")
  set(CPACK_NSIS_PACKAGE_NAME "dtk")
  set(CPACK_NSIS_MUI_FINISHPAGE_RUN "vt-platform.exe")
  set(CPACK_GENERATOR "NSIS")
  set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
  set(CPACK_NSIS_EXECUTABLES_DIRECTORY "bin")
endif(WIN32)

## #############################################################################
## Add postinst and prerm script
## #############################################################################

if(UNIX AND NOT APPLE)
  if(${DISTRIB} MATCHES fc|fedora|Fedora|Centos|centos|SUSE|Suse|suse)
    ## #################################################################
    ## RPM generator settings
    ## #################################################################


    set(CPACK_RPM_USER_BINARY_SPECFILE "${CMAKE_CURRENT_BINARY_DIR}/vt.spec")

    set(CPACK_GENERATOR "RPM")

#    set(CPACK_RPM_PACKAGE_DEBUG 1)

    set(CPACK_RPM_PACKAGE_REQUIRES "qt")

    set(CPACK_RPM_PACKAGE_REQUIRES "dtk")

    set(CPACK_RPM_PACKAGE_LICENSE "BSD")
    set(CPACK_RPM_PACKAGE_GROUP "System Environment/Libraries")
    set(CPACK_INSTALL_STUFF  ${CMAKE_BINARY_DIR})

    # set(CPACK_RPM_COMPONENT_INSTALL ON)
    # set(CPACK_COMPONENTS_ALL apps devel libs)

  else()
    ## #################################################################
    ## DEB generator settings
    ## #################################################################

    set(CPACK_GENERATOR "DEB")
  endif()
endif(UNIX AND NOT APPLE)

## #################################################################

include(CPack)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/vt.spec.in" "${CMAKE_CURRENT_BINARY_DIR}/vt.spec")