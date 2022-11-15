set(Libxc_pinned "6.0.0")

find_package(Libxc ${Libxc_pinned} CONFIG QUIET COMPONENTS C)

if(TARGET Libxc::xc)
  get_property(_loc TARGET Libxc::xc PROPERTY LOCATION)
  message(STATUS "Found Libxc: ${_loc} (found version ${Libxc_VERSION})")
else()
  message(STATUS "Suitable Libxc could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(Libxc
    QUIET
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    URL
      https://gitlab.com/libxc/libxc/-/archive/${Libxc_pinned}/libxc-${Libxc_pinned}.tar.gz
    )

  enable_language(C)

  # compile 3rd derivative code
  set(DISABLE_KXC OFF CACHE BOOL "" FORCE)
  # compile 4th derivative code
  set(DISABLE_LXC OFF CACHE BOOL "" FORCE)
  # disable compilation of testing infrastructure
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  # only build shared libraries
  set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
  # enable arch-dependent flags if they are used for VeloxChem
  set(ENABLE_XHOST ${ENABLE_ARCH_FLAGS} CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(Libxc)

  # Provide an alias, so linking to Libxc looks the same regardless if it was
  # found on the system or if it was fetched at configuration
  add_library(Libxc::xc ALIAS xc)
endif()
