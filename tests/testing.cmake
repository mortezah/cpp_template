# google tests
# Note: DISCOVERY_MODE PRE_TEST forces the build system to not call test post build
include(GoogleTest)
gtest_discover_tests(vector_algebra_utils_tests DISCOVERY_MODE PRE_TEST)

# # for slower tests
# if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RELEASE")
#
# endif()
