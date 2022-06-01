if(DEFINED VLX_USE_HIP OR DEFINED VLX_USE_CUDA)
  target_compile_definitions(vlxobjs
    PUBLIC
      VLX_USE_DEVICE
    )
endif()

if(DEFINED VLX_USE_HIP)
  # HIP-specific code
  add_subdirectory(hip)
elseif(DEFINED VLX_USE_CUDA)
  # CUDA-specific code
  add_subdirectory(cuda)
endif()

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
