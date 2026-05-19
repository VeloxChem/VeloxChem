
# Makefile for the Tabula shared library

LOCAL_LIB := tabulalib.so
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

# Default rule for linking the Tabula shared library — VeloxChem and Python
# symbols are resolved at load time from the already-imported veloxchemlib
debug: $(LOCAL_OBJ)
	@echo ====== Linking $(LOCAL_LIB)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -shared -undefined dynamic_lookup -o $(LOCAL_LIB) $^
	mkdir -p $(CURDIR)/../../build/python/veloxchem
	$(CP) $(LOCAL_LIB) $(CURDIR)/../../build/python/veloxchem

# Default rule for linking the Tabula shared library
release: $(LOCAL_OBJ)
	@echo ====== Linking $(LOCAL_LIB)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -shared -undefined dynamic_lookup -o $(LOCAL_LIB) $^
	mkdir -p $(CURDIR)/../../build/python/veloxchem
	$(CP) $(LOCAL_LIB) $(CURDIR)/../../build/python/veloxchem

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
