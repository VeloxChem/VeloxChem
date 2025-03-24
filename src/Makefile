# Include makefile setup

include Makefile.setup

# Set TARGET and linker flags for VeloxChem shared library

VLX_BAS := ../basis
VLX_TARGET := veloxchemlib.so
VLX_MODULE := $(wildcard pymodule/*.py)

LIB_LIST += -shared $(PYTHON_LD) $(OPENMP_LD) $(MATH_LD) $(LIBXC_LD)

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
VLX_LIBS_DIR := general
VLX_LIBS_DIR += math
VLX_LIBS_DIR += moldata
VLX_LIBS_DIR += orbdata
VLX_LIBS_DIR += python
VLX_LIBS_DIR += t2c_func
VLX_LIBS_DIR += t2c_geometrical_derivatives
VLX_LIBS_DIR += t2c_overlap
VLX_LIBS_DIR += t2c_kinetic_energy
VLX_LIBS_DIR += t2c_nuclear_potential
VLX_LIBS_DIR += t2c_electric_dipole_momentum
VLX_LIBS_DIR += t2c_nuclear_potential_geom_0X0
VLX_LIBS_DIR += t2c_electron_repulsion
VLX_LIBS_DIR += t4c_func
VLX_LIBS_DIR += t4c_electron_repulsion
VLX_LIBS_DIR += t4c_geometrical_derivatives
VLX_LIBS_DIR += dft_func
VLX_LIBS_DIR += t3c_overlap
VLX_LIBS_DIR += t3c_electron_repulsion
VLX_LIBS_DIR += t3c_func
VLX_LIBS_DIR += t3c_electron_repulsion_geom_X00
VLX_LIBS_DIR += t3c_electron_repulsion_geom_0X0
VLX_LIBS_DIR += visualization
VLX_LIBS_DIR +=	onee_ints
VLX_LIBS_DIR +=	pdft_func

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
	$(CP) -r ../database $(CURDIR)/../build/python/veloxchem
	$(CP) pymodule/mklconf.py.in $(CURDIR)/../build/python/veloxchem/mklconf.py

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
	$(CP) -r ../database $(CURDIR)/../build/python/veloxchem
	$(CP) pymodule/mklconf.py.in $(CURDIR)/../build/python/veloxchem/mklconf.py

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
