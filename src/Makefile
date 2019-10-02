#
#                      V.E.L.O.X. C.H.E.M. MP
#      ---------------------------------------------------
#           An Electronic Structure Code for Nanoscale
#
#  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
#  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

# Include makefile setup

include Makefile.setup

# Set TARGET and linker flags for VeloxChemMP shared library

VLX_BAS := ../basis
VLX_TARGET := veloxchemlib.so
VLX_MODULE := $(wildcard pymodule/*.py)
PYVLX_MAIN := vlx

LIB_LIST += -shared $(PYTHON_LD)

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
VLX_LIBS_DIR := ${shell find * -not -path pymodule -type d -print}

# Generate all internal libraries related files
VLX_LIBS = $(foreach d,$(VLX_LIBS_DIR),$(addsuffix .a,$(addprefix $d/,$d)))
VLX_LIBS_OBJS = $(foreach d,$(VLX_LIBS_DIR),$(addprefix $d/,*.o))
VLX_LIBS_DEPS = $(foreach d,$(VLX_LIBS_DIR),$(addprefix $d/,*.d))
VLX_LIBS_OPTS = $(foreach d,$(VLX_LIBS_DIR),$(addprefix $d/,*.optrpt))
VLX_LIBS_LIST = $(foreach d,$(VLX_LIBS_DIR),$(addprefix -L$d -l,$d))

# Debug
debug: $(VLX_LIBS)
	@echo ====== Linking debug version of $(VLX_TARGET)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(VLX_TARGET) \
	$(VLX_LIBS_OBJS) $(LIB_LIST)
	@echo ====== Moving $(VLX_TARGET) and copying $(PYVLX_MAIN)...
	mkdir -p ../build/python/veloxchem
	$(MV) $(VLX_TARGET) $(CURDIR)/../build/python/veloxchem
	$(CP) $(VLX_MODULE) $(CURDIR)/../build/python/veloxchem
	$(CP) -r $(VLX_BAS) $(CURDIR)/../build/python/veloxchem
	mkdir -p ../build/bin
	$(CP) $(PYVLX_MAIN) $(CURDIR)/../build/bin

# Release
release: $(VLX_LIBS)
	@echo ====== Linking release version of $(VLX_TARGET)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(VLX_TARGET) \
	$(VLX_LIBS_OBJS) $(LIB_LIST)
	@echo ====== Moving $(VLX_TARGET) and copying $(PYVLX_MAIN)...
	mkdir -p ../build/python/veloxchem
	$(MV) $(VLX_TARGET) $(CURDIR)/../build/python/veloxchem
	$(CP) $(VLX_MODULE) $(CURDIR)/../build/python/veloxchem
	$(CP) -r $(VLX_BAS) $(CURDIR)/../build/python/veloxchem
	mkdir -p ../build/bin
	$(CP) $(PYVLX_MAIN) $(CURDIR)/../build/bin

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
