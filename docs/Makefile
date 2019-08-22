# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

VLX_SOURCE = source
PYTHONPATH = $(VLX_SOURCE):mocks

VLX_EXCLUDE = $(VLX_SOURCE)/veloxchem/veloxchemlib.py
VLX_EXCLUDE += $(VLX_SOURCE)/veloxchem/aodensitymatrix.py
VLX_EXCLUDE += $(VLX_SOURCE)/veloxchem/aofockmatrix.py
VLX_EXCLUDE += $(VLX_SOURCE)/veloxchem/molecularbasis.py
VLX_EXCLUDE += $(VLX_SOURCE)/veloxchem/molecularorbitals.py
VLX_EXCLUDE += $(VLX_SOURCE)/veloxchem/molecule.py
VLX_EXCLUDE += $(VLX_SOURCE)/veloxchem/visualizationdriver.py

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@PYTHONPATH=$(PYTHONPATH) $(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

api: links
	PYTHONPATH=$(PYTHONPATH) sphinx-apidoc -o api $(VLX_SOURCE) $(VLX_EXCLUDE)
links:
	mkdir -p source/veloxchem
	for f in ../src/pymodule/*py; do ln -sf ../../$$f source/veloxchem; done
	ln -sf ../../mocks/veloxchemlib.py source/veloxchem
serve:
	python3 -c 'import webbrowser; webbrowser.open_new_tab("http://localhost:5000")' && python3 server.py
