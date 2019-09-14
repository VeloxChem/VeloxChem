#
#                     V.E.L.O.X. C.H.E.M. MP
#      ---------------------------------------------------
#           An Electronic Structure Code for Nanoscale
#
#  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
#  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

# Automatically generated settings

USE_MPI := true
USE_GPU := false
USE_MKL := true

MATH_INC := -I/opt/intel/mkl/include
MATH_LIB := -L/opt/intel/mkl/lib/intel64
MATH_LIB += -Wl,-rpath,/opt/intel/mkl/lib/intel64
MATH_LIB += -lmkl_rt -lmkl_gnu_thread
MATH_LIB += -lgomp -lpthread -lm -ldl

PYTHON := python3
PYBIND11_ROOT := /home/xinli/software/pybind11

CXX := mpicxx

CXX_REL_FLG := -fopenmp
CXX_DEB_FLG := -fopenmp

MACLIBS := 

GST_ROOT := /home/xinli/software/googletest/googletest
GST_LIB := /home/xinli/software/googletest/googletest/lib/libgtest.a

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

ifeq ($(USE_GPU),true)
  CPP_REL_FLG += -DENABLE_GPU
  CPP_DEB_FLG += -DENABLE_GPU
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
PYTHON_INC += -I$(PYBIND11_ROOT)/include
PYTHON_LD  := $(shell $(PYTHON)-config --ldflags) -fno-lto

CXX_REL_FLG += -fPIC $(PYTHON_INC)
CXX_DEB_FLG += -fPIC $(PYTHON_INC)

# Native shell commands

CP := cp -f
MV := mv -f
RM := rm -f

# Native libraries generation tools

AR      := ar
ARFLAGS := rc
