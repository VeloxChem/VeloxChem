#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
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

# Makefile for internal VeloxChem library

LOCAL_LIB := oneints_ovl.a
LOCAL_SCR := $(wildcard *.cpp)
LOCAL_OBJ := $(subst .cpp,.o,$(LOCAL_SCR))
LOCAL_DEP := $(subst .cpp,.d,$(LOCAL_SCR))
LOCAL_OPT := $(subst .cpp,.optrpt,$(LOCAL_SCR))

# Read in compilers & shell commands info
include ../Makefile.setup

# Select compiler flags
ifneq "$(MAKECMDGOALS)" "release"
CXXFLAGS := $(CXX_DEB_FLG)
CPPFLAGS := $(CPP_DEB_FLG)
else
CXXFLAGS := $(CXX_REL_FLG)
CPPFLAGS := $(CPP_REL_FLG)
endif

# Add includes directories
VLX_INCLUDES := ${shell find ../* -type d -print}

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

# Default rule for compiling cpp files
%.o: %.cpp
	@echo ====== Compiling $@...
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(LIB_INC) -o $@ $<

# Default rule for generating dependencies
%.d: %.cpp
	@echo ====== Updating dependencies for $@
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -M $< > $@
