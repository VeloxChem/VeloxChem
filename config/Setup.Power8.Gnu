#
#                     V.E.L.O.X. C.H.E.M. MP
#      ---------------------------------------------------
#           An Electronic Structure Code for Nanoscale
#
#  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
#  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

# NOTE: setup for monkey.theochem.kth.se

USE_MPI := true
USE_GPU := true
USE_MKL := false

# C++ compiler & linker

CXX := mpicxx

# Compiler flags for generic intel CPUs

CXX_REL_FLG := -std=c++11 -O3 -fopenmp -Wall
CXX_DEB_FLG := -std=c++11 -O3 -fopenmp -Wall -g

# Generic preprocessor flags

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

MATH_ROOT := /home/lxin/software/openblas
MATH_INC := -I$(MATH_ROOT)/include

LIB_LIST := -L$(MATH_ROOT)/lib -Wl,-rpath,$(MATH_ROOT)/lib -lopenblas
LIB_LIST += -lgomp -lpthread -lm -ldl

# Google Test root directory

GST_ROOT := /home/lxin/software/googletest/googletest

GST_LIB := $(GST_ROOT)/lib/libgtest.a

# CUDA compiler & linker

CUDA_ROOT := /usr/local/cuda-8.0
NVCC := $(CUDA_ROOT)/bin/nvcc

ifeq ($(USE_GPU),true)
  NVCC_REL_FLG := -std=c++11 -O3 --compiler-options "-fopenmp" --compiler-options "-Wall"
  NVCC_DEB_FLG := -std=c++11 -O3 --compiler-options "-fopenmp" --compiler-options "-Wall" --compiler-options "-g"

  NVCC_REL_FLG += --compiler-options "-fPIC"
  NVCC_DEB_FLG += --compiler-options "-fPIC"

  ifdef CUDA_ROOT
    LIB_LIST += -L$(CUDA_ROOT)/lib64
  endif

  LIB_LIST += -lcuda -lcudart
endif

# Python interface

PYTHON := python3
PYBIND11_ROOT := /home/lxin/software/pybind11

PYTHON_INC := $(shell $(PYTHON)-config --includes)
PYTHON_INC += -I$(shell $(PYTHON) -c "import mpi4py; print(mpi4py.get_include())")
PYTHON_INC += -I$(shell $(PYTHON) -c "import numpy; print(numpy.get_include())")
PYTHON_INC += -I$(PYBIND11_ROOT)/include

PYTHON_LD  := $(shell $(PYTHON)-config --ldflags)

CXX_REL_FLG += -fPIC $(PYTHON_INC)
CXX_DEB_FLG += -fPIC $(PYTHON_INC)

# Native shell commands

CP := cp -f
MV := mv -f
RM := rm -f

# Native libraries generation tools

AR      := ar
ARFLAGS := rc
