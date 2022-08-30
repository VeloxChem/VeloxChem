set(Libxc_pinned "5.2.3")

find_package(Libxc ${Libxc_pinned} CONFIG QUIET COMPONENTS C)

if(TARGET Libxc::xc)
  get_property(_loc TARGET Libxc::xc PROPERTY LOCATION)
  message(STATUS "Found Libxc: ${_loc} (found version ${Libxc_VERSION})")
else()
  message(STATUS "Suitable Libxc could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(Libxc
    QUIET
    URL
      https://gitlab.com/libxc/libxc/-/archive/${Libxc_pinned}/libxc-${Libxc_pinned}.tar.gz
    )

  enable_language(C)

  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
  set(ENABLE_XHOST ${ENABLE_ARCH_FLAGS} CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(Libxc)
endif()
