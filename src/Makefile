#
#                      V.E.L.O.X. C.H.E.M. MP
#      ---------------------------------------------------
#           An Electronic Structure Code for Nanoscale
#
#  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
#  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

# Include makefile setup

include Makefile.setup

# Set TARGET for VeloxChemMP package

VLX_TARGET := VeloxChemMP.so

# Select compiler flags
ifneq "$(MAKECMDGOALS)" "release"
  CXXFLAGS := $(CXX_DEB_FLG)
  CPPFLAGS := $(CPP_DEB_FLG)
else
  CXXFLAGS := $(CXX_REL_FLG)
  CPPFLAGS := $(CPP_REL_FLG)
endif 

# Add includes directories
ifeq ($(USE_GPU),true)
  VLX_INCLUDES := ${shell find * -type d -print}
else
  VLX_INCLUDES := ${shell find * -not -path gpu -type d -print}
endif

# Update CXX flags & includes path
CXXFLAGS += $(addprefix -I $(CURDIR)/,$(VLX_INCLUDES))
vpath %.hpp $(VLX_INCLUDES)

# Add MKL includes
CXXFLAGS += -I $(MKL_INCLUDES)

# Add includes/ldflags for Boost Python
CXX_REL_FLG += -fPIC ${PYTHON_INC}
CXX_DEB_FLG += -fPIC ${PYTHON_INC}
LIB_LIST += -shared ${PYTHON_LD}

# Add list of internal libraries
ifeq ($(USE_GPU),true)
  VLX_LIBS_DIR := ${shell find * -type d -print}
else
  VLX_LIBS_DIR := ${shell find * -not -path gpu -type d -print}
endif

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
	@echo ====== Moving $(VLX_TARGET) into $(CURDIR)/../build directory...
	mkdir -p ../build
	$(MV) $(VLX_TARGET) $(CURDIR)/../build

# Release
release: $(VLX_LIBS)
	@echo ====== Linking release version of $(VLX_TARGET)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(VLX_TARGET) \
	$(VLX_LIBS_OBJS) $(LIB_LIST)
	@echo ====== Moving $(VLX_TARGET) into $(CURDIR)/../build directory...
	mkdir -p ../build
	$(MV) $(VLX_TARGET) $(CURDIR)/../build

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
	make --directory=$@ $(MAKECMDGOALS)
