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

# Include makefile setup

include Makefile.setup

# Set TARGET and linker flags for VeloxChem shared library
VLX_BAS := ../basis
VLX_MODULE := $(wildcard pymodule/*.py)

LIB_LIST += -shared $(LTOFLAG)

# Select compiler flags
ifneq "$(MAKECMDGOALS)" "release"
  CXXFLAGS := $(CXX_DEB_FLG)
  CPPFLAGS := $(CPP_DEB_FLG)
else
  CXXFLAGS := $(CXX_REL_FLG)
  CPPFLAGS := $(CPP_REL_FLG)
endif 

# Add includes directories
VLX_INCLUDES := ${shell find * -not -path pymodule -type d -print}

# Update CXX flags & includes path
CXXFLAGS += $(addprefix -I $(CURDIR)/,$(VLX_INCLUDES))
vpath %.hpp $(VLX_INCLUDES)

# Add MKL includes
CXXFLAGS += $(MATH_INC)

# Add list of internal libraries
# exclude `pymodule` and any folders that DO NOT contain any .cpp files
VLX_LIBS_DIR := $(shell find * -not -path pymodule -type d -exec sh -c 'ls -1 "{}" | egrep -i -q "*\.cpp$$"' ';' -print)

# Generate all internal libraries related files
VLX_LIBS = $(foreach d,$(VLX_LIBS_DIR),$(addsuffix .a,$(addprefix $d/,$d)))
VLX_LIBS_OBJS = $(foreach d,$(VLX_LIBS_DIR),$(addprefix $d/,*.o))
VLX_LIBS_DEPS = $(foreach d,$(VLX_LIBS_DIR),$(addprefix $d/,*.d))
VLX_LIBS_OPTS = $(foreach d,$(VLX_LIBS_DIR),$(addprefix $d/,*.optrpt))
VLX_LIBS_LIST = $(foreach d,$(VLX_LIBS_DIR),$(addprefix -L$d -l,$d))

# Create build folder (if it doesn't exist already)
build_folder:
	@echo ====== Creating build folder...
	mkdir -p $(BUILD_LIB)

# Copy basis set files
basis: build_folder
	@echo ====== Copying basis set files...
	$(CP) -r $(VLX_BAS) $(BUILD_LIB)

# Copy module files
pymodule: build_folder
	@echo ====== Copying module files...
	$(CP) -r $(VLX_MODULE) $(BUILD_LIB)

# Copy pytest files
pytests: build_folder
	@echo ====== Copying Python test files...
	rsync -a --force --exclude='Makefile' --exclude='__pycache__' ../python_tests/ $(BUILD_LIB)/tests

# Debug
debug: $(VLX_LIBS) basis pymodule pytests
	@echo ====== Linking debug version of $(VLX_TARGET)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(VLX_TARGET) \
	$(VLX_LIBS_OBJS) $(LIB_LIST)
	@echo ====== Building bcastdebug program...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) bcastdebug.cpp -o bcastdebug.x \
	$(VLX_LIBS) $(MATH_LIB) general/general.a

# Release
release: $(VLX_LIBS) basis pymodule pytests
	@echo ====== Linking release version of $(VLX_TARGET)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(VLX_TARGET) \
	$(VLX_LIBS_OBJS) $(LIB_LIST)

# Clean up
clean:
	@echo ====== Removing temporary files...
	$(RM) $(VLX_LIBS)
	$(RM) $(VLX_LIBS_OBJS)
	$(RM) $(VLX_LIBS_DEPS)
	$(RM) $(VLX_LIBS_OPTS)

# Set internal libraries generation rule
.PHONY: $(VLX_LIBS) $(VLX_LIBS_DIR)
$(VLX_LIBS): $(VLX_LIBS_DIR)
$(VLX_LIBS_DIR):
	$(MAKE) --directory=$@ $(MAKECMDGOALS)
