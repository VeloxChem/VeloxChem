#
#                      V.E.L.O.X. C.H.E.M. MP
#      ---------------------------------------------------
#           An Electronic Structure Code for Nanoscale
#
#  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
#  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

# Makefile for internal Velox Chem MP library

LOCAL_LIB := gpu.a
LOCAL_SCR := $(wildcard *.cpp)
LOCAL_OBJ := $(subst .cpp,.o,$(LOCAL_SCR))
LOCAL_DEP := $(subst .cpp,.d,$(LOCAL_SCR))
LOCAL_OPT := $(subst .cpp,.optrpt,$(LOCAL_SCR))

LOCAL_CU_SCR := $(wildcard *.cu)
LOCAL_CU_OBJ := $(subst .cu,.o,$(LOCAL_CU_SCR))
LOCAL_CU_DEP := $(subst .cu,.d,$(LOCAL_CU_SCR))
LOCAL_CU_OPT := $(subst .cu,.optrpt,$(LOCAL_CU_SCR))

# Read in compilers & shell commands info
include ../Makefile.setup

# Add cuda objects
ifeq ($(USE_GPU),true)
LOCAL_SCR += $(LOCAL_CU_SCR)
LOCAL_OBJ += $(LOCAL_CU_OBJ)
LOCAL_DEP += $(LOCAL_CU_DEP)
LOCAL_OPT += $(LOCAL_CU_OPT)
endif

# Select compiler flags
ifneq "$(MAKECMDGOALS)" "release"
NVCCFLAGS := $(NVCC_DEB_FLG)
CXXFLAGS := $(CXX_DEB_FLG)
CPPFLAGS := $(CPP_DEB_FLG)
else
NVCCFLAGS := $(NVCC_REL_FLG)
CXXFLAGS := $(CXX_DEB_FLG)
CPPFLAGS := $(CPP_REL_FLG)
endif

# Add includes directories
VLX_INCLUDES := ${shell find ../* -type d -print}

# Update CXX flags & includes path
NVCCFLAGS += $(addprefix -I $(CURDIR)/,$(VLX_INCLUDES))
vpath %.hpp $(VLX_INCLUDES)

# Add MKL includes
NVCCFLAGS += $(MATH_INC)

# Update CXX flags & includes path
CXXFLAGS += $(addprefix -I $(CURDIR)/,$(VLX_INCLUDES))
vpath %.hpp $(VLX_INCLUDES)

# Add MKL includes
CXXFLAGS += $(MATH_INC)

# Default rule for linking library
debug: $(LOCAL_OBJ)
	@echo ====== Linking $(LOCAL_LIB)...
	$(AR) $(ARFLAGS) $(LOCAL_LIB) $^

# Default rule for linking library
release: $(LOCAL_OBJ)
	@echo ====== Linking $(LOCAL_LIB)...
	$(AR) $(ARFLAGS) $(LOCAL_LIB) $^

# Default rule for cleaning library directory
clean:
	@echo ====== Removing temporary files...
	$(RM) $(LOCAL_LIB) $(LOCAL_DEP) $(LOCAL_OBJ) $(LOCAL_OPT)

# Include dependency files
ifneq "$(MAKECMDGOALS)" "clean"
  -include $(LOCAL_DEP)
endif

# Default rule for compiling cu files
%.o: %.cu
	@echo ====== Compiling $@...
	$(NVCC) -c $(NVCCFLAGS) $(CPPFLAGS) $(LIB_INC) $(MPI_INC) -o $@ $<

# Default rule for generating dependencies
%.d: %.cu
	@echo ====== Updating dependencies for $@
	$(NVCC) $(NVCCFLAGS) $(CPPFLAGS) -M $< > $@

# Default rule for compiling cpp files
%.o: %.cpp
	@echo ====== Compiling $@...
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(LIB_INC) -o $@ $<

# Default rule for generating dependencies
%.d: %.cpp
	@echo ====== Updating dependencies for $@
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -M $< > $@
