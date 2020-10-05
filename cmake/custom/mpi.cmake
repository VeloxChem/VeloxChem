#.rst:
#
# Enables MPI support.
# This was adapted from Autocmake
#
# Variables used::
#
#   USE_MPI
#
# autocmake.yml configuration::
#
#   docopt: "--mpi Enable MPI parallelization [default: False]."
#   define: "'-DUSE_MPI={0}'.format(arguments['--mpi'])"

option(ENABLE_MPI "Enable MPI parallelization" OFF)

if(USE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()
