#
#                     V.E.L.O.X. C.H.E.M. MP
#      ---------------------------------------------------
#           An Electronic Structure Code for Nanoscale
#
#  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
#  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

# NOTE: setup for monkey.theochem.kth.se

USE_MPI    := true
USE_GPU    := true
USE_PYTHON := false

# C++ compiler & linker

CXX := mpiicpc

# Compiler flags for generic intel CPUs

CXX_REL_FLG := -std=c++11 -O3 -xHost -qopenmp -Wall -qoverride-limits
CXX_DEB_FLG := -std=c++11 -O3 -xHost -qopenmp -Wall -qoverride-limits -g
CXX_DEB_FLG += -qopt-report=5 -qopt-report-phase=vec

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

ifeq ($(USE_PYTHON),true)
  CPP_REL_FLG += -DENABLE_PYTHON -DBOOST_NO_AUTO_PTR
  CPP_DEB_FLG += -DENABLE_PYTHON -DBOOST_NO_AUTO_PTR
endif

# MKL header and library

MKL_INCLUDES := $(MKLROOT)/include

LIB_LIST := -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64
LIB_LIST += -lmkl_intel_thread -lmkl_core -lmkl_avx -liomp5 -lpthread -lm -ldl

# Google Test root directory

GST_ROOT := /home/rinkevic/Development/FrameWorks/GoogleTest/googletest

GST_LIB := $(GST_ROOT)/lib/gtest_main.a

# CUDA compiler & linker

CUDA_ROOT := /usr/local/cuda
NVCC := ${CUDA_ROOT}/bin/nvcc

ifeq ($(USE_GPU),true)
  NVCC_REL_FLG := -std=c++11 -O3 --compiler-options "-fopenmp" --compiler-options "-Wall"
  NVCC_DEB_FLG := -std=c++11 -O3 --compiler-options "-fopenmp" --compiler-options "-Wall" --compiler-options "-g"

  ifeq ($(USE_PYTHON),true)
    NVCC_REL_FLG += --compiler-options "-fPIC"
    NVCC_DEB_FLG += --compiler-options "-fPIC"
  endif

  ifdef CUDA_ROOT
    LIB_LIST += -L${CUDA_ROOT}/lib64
  endif

  LIB_LIST += -lcuda -lcudart
endif

# Boost header and library

PYTHON := python3
BOOST_PATH := /home/xinli/software/boost

ifeq ($(USE_PYTHON),true)
  PYTHON_INC := $(shell $(PYTHON)-config --includes)
  PYTHON_INC += -I$(shell $(PYTHON) -c "import mpi4py; print(mpi4py.get_include())")
  PYTHON_INC += -I$(shell $(PYTHON) -c "import numpy; print(numpy.get_include())")
  PYTHON_INC += -I$(BOOST_PATH)/include

  PYTHON_LD  := $(shell $(PYTHON)-config --ldflags)
  PYTHON_LD  += -L$(BOOST_PATH)/lib -lboost_python3 -lboost_numpy3

  CXX_REL_FLG += -fPIC ${PYTHON_INC}
  CXX_DEB_FLG += -fPIC ${PYTHON_INC}

  LIB_LIST += -shared ${PYTHON_LD}
endif

# Native shell commands

MV := mv -f
RM := rm -f

# Native libraries generation tools

AR      := ar
ARFLAGS := rc
