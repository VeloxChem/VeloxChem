# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

VLX_SOURCE = "Reference Guide"
PYTHONPATH = $(VLX_SOURCE):mocks

VLX_EXCLUDE = $(VLX_SOURCE)/veloxchem/veloxchemlib.py
VLX_EXCLUDE += $(VLX_SOURCE)/veloxchem/aodensitymatrix.py
VLX_EXCLUDE += $(VLX_SOURCE)/veloxchem/aofockmatrix.py
VLX_EXCLUDE += $(VLX_SOURCE)/veloxchem/main.py
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

api: clean links log
	PYTHONPATH=$(PYTHONPATH) sphinx-apidoc -o api $(VLX_SOURCE) $(VLX_EXCLUDE)
	echo "   log" >> api/modules.rst

links:
	mkdir -p $(VLX_SOURCE)/veloxchem
	for f in ../src/pymodule/*py; do ln -sf ../../$$f $(VLX_SOURCE)/veloxchem; done
	ln -sf ../../mocks/veloxchemlib.py $(VLX_SOURCE)/veloxchem

log:
	mkdir -p api
	touch api/log.rst
	echo "Change log\n==========\n::\n" > api/log.rst
	git log --date=short --pretty=format:'    %h %cd %s' >> api/log.rst

serve:
	python3 -c 'import webbrowser; webbrowser.open_new_tab("http://localhost:5000")' && python3 server.py

clean:
	rm -rf _build api/*
