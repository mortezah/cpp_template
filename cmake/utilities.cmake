include(GNUInstallDirs)
include(InstallRequiredSystemLibraries)
function(always_full_rpath)
  # CMake RPATH "always full" configuration, see:
  # https://cmake.org/Wiki/CMake_RPATH_handling#Always_full_RPATH
  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH False PARENT_SCOPE)
  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH False PARENT_SCOPE)

  if(APPLE)
      set(base @loader_path)
  else()
      set(base $ORIGIN)
  endif()
  file(RELATIVE_PATH relDir
       ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
       ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
  )
  set(CMAKE_INSTALL_RPATH ${base} ${base}/${relDir} PARENT_SCOPE)
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True PARENT_SCOPE)
endfunction(always_full_rpath)

macro(begin_package)
  message(STATUS "CMAKE_VERSION: ${CMAKE_VERSION}")
  if (${PROJECT_NAME}_VERSION)
    message(STATUS "${PROJECT_NAME}_VERSION: ${${PROJECT_NAME}_VERSION}")
  endif()
  option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
  always_full_rpath()
  message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
endmacro(begin_package)

function(export_target tgt_name)
  install(TARGETS ${tgt_name} EXPORT ${tgt_name}-target
      RUNTIME DESTINATION bin
      ARCHIVE DESTINATION lib
      LIBRARY DESTINATION lib)
  install(EXPORT ${tgt_name}-target NAMESPACE ${PROJECT_NAME}::
          DESTINATION lib/cmake/${PROJECT_NAME})
  set(${PROJECT_NAME}_EXPORTED_TARGETS
      ${${PROJECT_NAME}_EXPORTED_TARGETS} ${tgt_name} PARENT_SCOPE)
endfunction(export_target)

macro(export_library target)
  export_target(${target})
  install(FILES ${HEADERS} DESTINATION include)
endmacro(export_library)

function(post_build_install target_name destination)
  install(TARGETS ${target_name} DESTINATION ${destination})
endfunction()

macro(export_object obj_name)
  set(${PROJECT_NAME}_EXPORTED_OBJECTS
    ${${PROJECT_NAME}_EXPORTED_OBJECTS} ${obj_name} PARENT_SCOPE)
  install(FILES ${HEADERS} DESTINATION include)
endmacro(export_object)

macro(end_subdir)
  set(${PROJECT_NAME}_EXPORTED_TARGETS
      ${${PROJECT_NAME}_EXPORTED_TARGETS} PARENT_SCOPE)
  set(${PROJECT_NAME}_DEPS ${${PROJECT_NAME}_DEPS} PARENT_SCOPE)
  set(${PROJECT_NAME}_DEP_PREFIXES ${${PROJECT_NAME}_DEP_PREFIXES} PARENT_SCOPE)
endmacro(end_subdir)

function(end_package)
  include(CMakePackageConfigHelpers)
  set(INCLUDE_INSTALL_DIR include)
  set(LIB_INSTALL_DIR lib)
  set(CONFIG_CONTENT "
set(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_VERSION})
include(CMakeFindDependencyMacro)
# we will use find_dependency, but we don't want to force
# our users to have to specify where all of our dependencies
# were installed; that defeats the whole point of automatically
# importing dependencies.
# since the documentation for find_dependency() doesn't mention
# a PATHS argument, we'll temporarily add the prefixes to
# CMAKE_PREFIX_PATH.
set(${PROJECT_NAME}_DEPS \"${${PROJECT_NAME}_DEPS}\")
set(${PROJECT_NAME}_DEP_PREFIXES \"${${PROJECT_NAME}_DEP_PREFIXES}\")
set(${PROJECT_NAME}_BACKUP_PREFIX_PATH \"\${CMAKE_PREFIX_PATH}\")
set(CMAKE_PREFIX_PATH \"\${${PROJECT_NAME}_DEP_PREFIXES};\${CMAKE_PREFIX_PATH}\")
foreach(dep IN LISTS ${PROJECT_NAME}_DEPS)
  find_dependency(\${dep})
endforeach()
set(CMAKE_PREFIX_PATH \"\${${PROJECT_NAME}_BACKUP_PREFIX_PATH}\")
set(${PROJECT_NAME}_EXPORTED_TARGETS \"${${PROJECT_NAME}_EXPORTED_TARGETS}\")
foreach(tgt IN LISTS ${PROJECT_NAME}_EXPORTED_TARGETS)
  include(\${CMAKE_CURRENT_LIST_DIR}/\${tgt}-target.cmake)
endforeach()
foreach(_comp \${${PROJECT_NAME}_FIND_COMPONENTS})
  if (NOT \";\${${PROJECT_NAME}_EXPORTED_TARGETS};\" MATCHES \${_comp})
    set(${PROJECT_NAME}_\${_comp}_FOUND False)
    if ( ${PROJECT_NAME}_FIND_REQUIRED_\${_comp} )
      MESSAGE(SEND_ERROR \"Required ${PROJECT_NAME} component not found: \${_comp}\")
    endif()
  else()
    set(${PROJECT_NAME}_\${_comp}_FOUND True)
  endif()
  MESSAGE(STATUS \"${PROJECT_NAME} component \${_comp} found: \${${PROJECT_NAME}_\${_comp}_FOUND}\")
endforeach()
set(${PROJECT_NAME}_COMPILER \"${CMAKE_CXX_COMPILER}\")
set(${PROJECT_NAME}_CXX_FLAGS \"${CMAKE_CXX_FLAGS}\")
set(${PROJECT_NAME}_INCLUDE_DIRS \"${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_DIR}\")
set(${PROJECT_NAME}_LIBRARIES \"${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}\")
")
  install(FILES
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    DESTINATION lib/cmake/${PROJECT_NAME})
  if(PROJECT_VERSION)
    file(WRITE
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        "${CONFIG_CONTENT}")
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)
    install(FILES
      "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
      DESTINATION lib/cmake/${PROJECT_NAME})
  endif()
endfunction(end_package)

function(add_test_nogtest execname srcname)
  add_executable(${execname} ${srcname})
  target_link_libraries(${execname} PRIVATE Eigen3::Eigen CGAL::CGAL)
  target_link_libraries(${execname} PRIVATE utilities)
endfunction(add_test_nogtest)

function(add_test_ execname srcname)
  add_test_nogtest(${execname} ${srcname})
  target_link_libraries(${execname} PRIVATE GTest::gtest_main GTest::gmock_main)
  # this is too copy dlls to where the tests are for cases were the test is called post build
  if (WIN32 OR MSVC)
    add_custom_command(
      TARGET ${execname} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:${execname}> $<TARGET_FILE_DIR:${execname}>
      COMMAND_EXPAND_LISTS)
  endif()
  export_target(${execname})
endfunction(add_test_)

function(add_exec execname srcname)
  add_test_nogtest(${execname} ${srcname})
  export_target(${execname})
  # manually disable compile options for examples
  set_target_properties(${execname} PROPERTIES COMPILE_OPTIONS "")
endfunction(add_exec)

function(copy_target_to_destination target destination)
  get_target_property(deps ${target} INTERFACE_LINK_LIBRARIES)
  message(STATUS "INTERFACE_LINK_LIBRARIES are ${INTERFACE_LINK_LIBRARIES}")
  # Copy the target itself.
  add_custom_command(TARGET ${target} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${target}> ${destination}/$<TARGET_FILE_NAME:${target}>
    COMMENT "Copying $<TARGET_FILE:${target}> and its dependencies to ${destination}")
endfunction()

function(declare_system_library target)
  message(STATUS "Declaring system library ${target}")
  get_target_property(target_aliased_name ${target} ALIASED_TARGET)
  if (target_aliased_name)
    set(target ${target_aliased_name})
  endif()
  set_target_properties(${target} PROPERTIES
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES $<TARGET_PROPERTY:${target},INTERFACE_INCLUDE_DIRECTORIES>)
endfunction()
