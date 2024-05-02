# sources
set(SOURCES
  vector_algebra_utils.cpp
)

# headers
set(HEADERS
  cgal_typedefs.h
  vector_algebra_utils.h
)

add_library(utilities ${SOURCES})
if (BUILD_SHARED)
  set_property(TARGET utilities PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

target_link_libraries(utilities PRIVATE CGAL::CGAL)
foreach(comp ${CGAL_REQUESTED_COMPONENTS})
  if(TARGET CGAL::CGAL_${comp})
    target_link_libraries(utilities PRIVATE CGAL::CGAL_${comp})
  endif()
endforeach()
if(CGAL_3RD_PARTY_LIBRARIES)
  target_link_libraries(utilities PRIVATE ${CGAL_3RD_PARTY_LIBRARIES})
endif()
# Include directories
target_include_directories(utilities PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

export_library(utilities)

end_subdir()
