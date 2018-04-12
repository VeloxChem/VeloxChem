#
#                     V.E.L.O.X. C.H.E.M. MP
#      ---------------------------------------------------
#           An Electronic Structure Code for Nanoscale
#
#  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
#  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

# NOTE: setup for monkey.theochem.kth.se

# C++ compiler & linker

CXX := mpiicpc

# Compiler flags for generic intel CPUs

CXX_REL_FLG := -std=c++11 -Wall -O3 -xHost -qopenmp
CXX_DEB_FLG := -std=c++11 -Wall -O3 -xHost -qopenmp -g -qopt-report=5
CXX_DEB_FLG += -qopt-report-phase=vec

# Generic preprocessor flags

CPP_REL_FLG := -DVLX_ALIGN=64
CPP_DEB_FLG := -DVAR_DEBUG=1 -DVLX_ALIGN=64

# Google Test root directory

GST_ROOT := /home/rinkevic/Development/FrameWorks/GoogleTest/googletest

GST_LIB := $(GST_ROOT)/lib/gtest_main.a

# MKL includes

MKL_INCLUDES := $(MKLROOT)/include

# Library includes

LIB_LIST := -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64
LIB_LIST += -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

# Native shell commands

MV := mv -f
RM := rm -f

# Native libraries generation tools

AR      := ar
ARFLAGS := rc
