#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

# Makefile for internal VeloxChem MP library

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
