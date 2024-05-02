# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.22)

# version needed to generate CPPTEMPLATEConfig.cmake
project( CPPTEMPLATE VERSION 1.0.0 )

if(MSVC) # Windows build
  message(STATUS "CMAKE_HOST_SYSTEM_NAME ${CMAKE_HOST_SYSTEM_NAME}")
  message(STATUS "CMAKE_HOST_SYSTEM_VERSION ${CMAKE_HOST_SYSTEM_VERSION}")
  add_compile_options(/permissive- /w14242 /w14254 /w14263 /w14265 /w14287 /we4289 /w14296 /w14311 /w14545 /w14546 /w14547 /w14549 /w14555 /w14640 /w14826 /w14905 /w14906 /w14928 /Zm2000 /MP)
else() # Linux build
  message(STATUS "CMAKE_HOST_SYSTEM_NAME ${CMAKE_HOST_SYSTEM_NAME}")
  message(STATUS "CMAKE_HOST_SYSTEM_VERSION ${CMAKE_HOST_SYSTEM_VERSION}")
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13)
    message(FATAL_ERROR "gcc version has to be greater or eqaul 13")
  endif()
  add_compile_options(-Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Wparentheses)
endif()

# check GCC version: it has to be greater or equal 13
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  add_definitions(-DGEOMDEBUG)
  # to make sure cgal goes through the same path during debugging.
  # CGAL has a lot of internal checks (some related to exact computation) which we
  # do not wanna get involved with.
  #
  # see page https://doc.cgal.org/latest/STL_Extension/index.html (section 10 Checks)
  add_definitions(-DCGAL_NDEBUG)
endif()

include(cmake/utilities.cmake)
include(GenerateExportHeader)


# require c++23
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSION OFF)
set(CMAKE_CXX_STANDARD 23)

begin_package()

if(ENABLE_TESTING)
  include(CTest)
  enable_testing()
endif()

option(BUILD_SHARED "Enable Building of Shared Libraries" ON)
message(STATUS "BUILD_SHARED: ${BUILD_SHARED}")



find_package(Eigen3 CONFIG REQUIRED)
find_package(CGAL REQUIRED)

# add individual libraries
add_subdirectory(utilities)

# add the directory including the apps
add_subdirectory(clis)



# add tests if testing is enabled
if(BUILD_TESTING)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        f8d7d77c06936315286eb55f8de22cd23c188571 # v1.14
    )
  FetchContent_MakeAvailable(googletest)
  add_subdirectory(tests)
endif()

end_package()
