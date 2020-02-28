#
#                           VELOXCHEM 1.0-RC
#      ---------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact

# Automatically generated settings

USE_MPI := true
USE_MKL := true

MATH_INC := -I/opt/intel/compilers_and_libraries_2019.0.117/linux/mkl/include
MATH_LIB := -L/opt/intel/compilers_and_libraries_2019.0.117/linux/mkl/lib/intel64
MATH_LIB += -Wl,-rpath,/opt/intel/compilers_and_libraries_2019.0.117/linux/mkl/lib/intel64
MATH_LIB += -lmkl_rt -lmkl_gnu_thread
MATH_LIB += -lgomp -lpthread -lm -ldl

PYTHON := python3
PYTHON_USER_INC := 

CXX := mpicxx

CXX_REL_FLG := -fopenmp
CXX_DEB_FLG := -fopenmp

MACLIBS := 

LTOFLAG := -fno-lto

# Generic settings

# Compiler flags

CXX_REL_FLG += -std=c++11 -O3 -Wall
CXX_DEB_FLG += -std=c++11 -O3 -Wall -g

# Preprocessor flags

CPP_REL_FLG := -DVLX_ALIGN=64
CPP_DEB_FLG := -DVLX_ALIGN=64 -DVAR_DEBUG=1

ifeq ($(USE_MPI),true)
  CPP_REL_FLG += -DENABLE_MPI=1
  CPP_DEB_FLG += -DENABLE_MPI=1
else
  CPP_REL_FLG += -DENABLE_MPI=0
  CPP_DEB_FLG += -DENABLE_MPI=0
endif

ifeq ($(USE_MKL),true)
  CPP_REL_FLG += -DENABLE_MKL
  CPP_DEB_FLG += -DENABLE_MKL
endif

# Math library

LIB_LIST := $(MATH_LIB)

# Extra flag for MacOS

LIB_LIST += $(MACLIBS)

# Python interface

PYTHON_INC := $(shell $(PYTHON)-config --includes)
PYTHON_INC += -I$(shell $(PYTHON) -c "import mpi4py; print(mpi4py.get_include())")
PYTHON_INC += -I$(shell $(PYTHON) -c "import numpy; print(numpy.get_include())")
PYTHON_INC += -I$(PYTHON_USER_INC) $(shell $(PYTHON) -m pybind11 --includes)
PYTHON_LD  := $(shell $(PYTHON)-config --ldflags) $(LTOFLAG)

CXX_REL_FLG += -fPIC $(PYTHON_INC)
CXX_DEB_FLG += -fPIC $(PYTHON_INC)

# Native shell commands

CP := cp -f
MV := mv -f
RM := rm -f

# Native libraries generation tools

AR      := ar
ARFLAGS := rc
