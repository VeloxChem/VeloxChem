
# Makefile for internal Velox Chem library

LOCAL_LIB := t2c_overlap_geom_300.a
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
