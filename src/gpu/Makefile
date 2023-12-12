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

LOCAL_DEV_SCR := $(wildcard *.hip)
LOCAL_DEV_OBJ := $(subst .hip,.o,$(LOCAL_DEV_SCR))
LOCAL_DEV_DEP := $(subst .hip,.d,$(LOCAL_DEV_SCR))
LOCAL_DEV_OPT := $(subst .hip,.optrpt,$(LOCAL_DEV_SCR))

# Read in compilers & shell commands info
include ../Makefile.setup

# Add device objects
ifeq ($(USE_GPU),true)
LOCAL_SCR += $(LOCAL_DEV_SCR)
LOCAL_OBJ += $(LOCAL_DEV_OBJ)
LOCAL_DEP += $(LOCAL_DEV_DEP)
LOCAL_OPT += $(LOCAL_DEV_OPT)
endif

# Select compiler flags
ifneq "$(MAKECMDGOALS)" "release"
DEVCCFLAGS := $(DEVCC_DEB_FLG)
CXXFLAGS := $(CXX_DEB_FLG)
CPPFLAGS := $(CPP_DEB_FLG)
else
DEVCCFLAGS := $(DEVCC_REL_FLG)
CXXFLAGS := $(CXX_DEB_FLG)
CPPFLAGS := $(CPP_REL_FLG)
endif

# Add includes directories
VLX_INCLUDES := ${shell find ../* -type d -print}

# Update CXX flags & includes path
DEVCCFLAGS += $(addprefix -I $(CURDIR)/,$(VLX_INCLUDES))
vpath %.hpp $(VLX_INCLUDES)

# Add MKL includes
DEVCCFLAGS += $(MATH_INC)

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
%.o: %.hip
	@echo ====== Compiling $@...
	$(DEVCC) -c $(DEVCCFLAGS) $(CPPFLAGS) $(LIB_INC) -o $@ $<

# Default rule for generating dependencies
%.d: %.hip
	@echo ====== Updating dependencies for $@
	$(DEVCC) $(DEVCCFLAGS) $(CPPFLAGS) -M $< > $@

# Default rule for compiling cpp files
%.o: %.cpp
	@echo ====== Compiling $@...
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(LIB_INC) -o $@ $<

# Default rule for generating dependencies
%.d: %.cpp
	@echo ====== Updating dependencies for $@
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -M $< > $@
