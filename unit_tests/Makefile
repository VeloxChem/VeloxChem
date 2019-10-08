#
#                      V.E.L.O.X. C.H.E.M. MP
#      ---------------------------------------------------
#           An Electronic Structure Code for Nanoscale
#
#  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
#  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

# Makefile for unit tests suite

include ../src/Makefile.setup

# Set TARGET for unit tests suite

VLX_TARGET   := UnitTestVeloxChemMP.x
VLX_MAIN_SCR := $(subst .x,.cpp,$(VLX_TARGET))
VLX_MAIN_OBJ := $(subst .x,.o,$(VLX_TARGET))
VLX_MAIN_DEP := $(subst .x,.d,$(VLX_TARGET))
VLX_MAIN_OPT := $(subst .x,.optrpt,$(VLX_TARGET))
VLX_MAIN_OPT += ipo_out.optrpt

# Select compiler flags
ifneq "$(MAKECMDGOALS)" "release"
  CXXFLAGS := $(CXX_DEB_FLG)
  CPPFLAGS := $(CPP_DEB_FLG)
else
  CXXFLAGS := $(CXX_REL_FLG)
  CPPFLAGS := $(CPP_REL_FLG)
endif

# Add includes directories
VLX_INCLUDES := ${shell find * -type d -print}

# Update CXX flags & includes path
CXXFLAGS += $(addprefix -I $(CURDIR)/,$(VLX_INCLUDES))
vpath %.hpp $(VLX_INCLUDES)

# Add MKL includes
CXXFLAGS += $(MATH_INC)

# Add Google Test includes
CXXFLAGS += -I $(GST_ROOT)/include

# Add list of internal libraries
VST_LIBS_DIR := ${shell find * -type d -print}

# Generate all internal libraries related files
VST_LIBS = $(foreach d,$(VST_LIBS_DIR),$(addsuffix .a,$(addprefix $d/,$d)))
VST_LIBS_OBJS = $(foreach d,$(VST_LIBS_DIR),$(addprefix $d/,*.o))
VST_LIBS_DEPS = $(foreach d,$(VST_LIBS_DIR),$(addprefix $d/,*.d))
VST_LIBS_OPTS = $(foreach d,$(VST_LIBS_DIR),$(addprefix $d/,*.optrpt))
VST_LIBS_LIST = $(foreach d,$(VST_LIBS_DIR),$(addprefix -L$d -l,$d))

# Add list of internal libraries
VLX_LIBS_DIR := ${shell find ../src/*        \
                  -not -path ../src/pymodule \
                  -not -path ../src/python   \
                  -type d -print}

# Generate all internal libraries related files
VLX_LIBS = $(foreach d,$(VLX_LIBS_DIR),$(addsuffix .a,$(addprefix $d/,$d)))
VLX_LIBS_OBJS = $(foreach d,$(VLX_LIBS_DIR),$(addprefix $d/,*.o))
VLX_LIBS_DEPS = $(foreach d,$(VLX_LIBS_DIR),$(addprefix $d/,*.d))
VLX_LIBS_OPTS = $(foreach d,$(VLX_LIBS_DIR),$(addprefix $d/,*.optrpt))
VLX_LIBS_LIST = $(foreach d,$(VLX_LIBS_DIR),$(addprefix -L$d -l,$d))

# Debug
debug: $(VLX_MAIN_OBJ) $(VLX_LIBS) $(VST_LIBS)
	@echo ====== Linking debug version of $(VLX_TARGET)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(VLX_TARGET) $(VLX_MAIN_OBJ) \
		$(VLX_LIBS_OBJS) $(VST_LIBS_OBJS) $(LIB_LIST) $(GST_LIB)
	@echo ====== Moving $(VLX_TARGET) into $(CURDIR)/../build/bin directory...
	mkdir -p ../build/bin
	$(MV) $(VLX_TARGET) $(CURDIR)/../build/bin

# Release
release: $(VLX_MAIN_OBJ) $(VLX_LIBS) $(VST_LIBS)
	@echo ====== Linking release version of $(VLX_TARGET)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(VLX_TARGET) $(VLX_MAIN_OBJ) \
		$(VLX_LIBS_OBJS) $(VST_LIBS_OBJS) $(LIB_LIST) $(GST_LIB)
	@echo ====== Moving $(VLX_TARGET) into $(CURDIR)/../build/bin directory...
	mkdir -p ../build/bin
	$(MV) $(VLX_TARGET) $(CURDIR)/../build/bin

# Clean up
clean:
	@echo ====== Removing temporary files...
	$(RM) $(VLX_MAIN_DEP) $(VLX_MAIN_OBJ) $(VLX_MAIN_OPT) $(VLX_LIBS) $(VST_LIBS)
	$(RM) $(VLX_LIBS_OBJS)
	$(RM) $(VST_LIBS_OBJS)
	$(RM) $(VLX_LIBS_DEPS)
	$(RM) $(VST_LIBS_DEPS)
	$(RM) $(VLX_LIBS_OPTS)
	$(RM) $(VST_LIBS_OPTS)
	
# Request dependencies generation
ifneq "$(MAKECMDGOALS)" "clean"
  -include $(VLX_MAIN_DEP)
endif

# Set TARGET compilation rule
$(VLX_MAIN_OBJ): $(VLX_MAIN_SCR)
	@echo ====== Compiling $(VLX_MAIN_OBJ)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LIB_INC) -c $(VLX_MAIN_SCR) -o $(VLX_MAIN_OBJ)

# Set TARGET dependencies rule
$(VLX_MAIN_DEP): $(VLX_MAIN_SCR)
	@echo ====== Updating dependencies in $(VLX_MAIN_OBJ)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -M $(VLX_MAIN_SCR) > $(VLX_MAIN_DEP)

# Set internal libraries generation rule
.PHONY: $(VLX_LIBS) $(VLX_LIBS_DIR)
$(VLX_LIBS): $(VLX_LIBS_DIR)
$(VLX_LIBS_DIR):
	$(MAKE) --directory=$@ $(MAKECMDGOALS)

# Set internal libraries generation rule
.PHONY: $(VST_LIBS) $(VST_LIBS_DIR)
$(VST_LIBS): $(VST_LIBS_DIR)
$(VST_LIBS_DIR):
	$(MAKE) --directory=$@ $(MAKECMDGOALS)
