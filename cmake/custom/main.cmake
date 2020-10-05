# Export compile commands for each file to JSON
# This is useful for static analysis tools and linters
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED COMPONENTS CXX)

include(${CMAKE_CURRENT_LIST_DIR}/FindMKL.cmake)

add_subdirectory(src)
