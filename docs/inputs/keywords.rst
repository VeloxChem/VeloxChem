Details of the input keywords
=============================

The jobs group
--------------

- **task**

  - the type of calculation

  .. csv-table::
    :widths: 1, 3

    "``scf``", "Run self-consistent field calculation."
    "``response``", "Run response calculation."
    "``pulses``", "Run response calculation for pulses."
    "``mp2``", "Run second-order Møller–Plesset perturbation theory."
    "``visualization``", "Generate cube files for visualization."
    "``loprop``", "Run LoProp calculation."
    "``exciton``", "Run ab initio exciton model calculation."
    "``optimize``", "Run geometry optimization."

The method settings group
-------------------------

- **basis**

  - the basis set

  .. csv-table::
    :widths: 1, 3

    "``cc-pVDZ``", "Correlation consistent Dunning basis set, double-zeta"
    "``aug-cc-pVDZ``", "Correlation consistent Dunning basis set, double-zeta
    with diffuse functions"
    "``def2-SVP``", "Karlsruhe basis set, split valence polarization"
    "``def2-SVPD``", "Karlsruhe basis set, split valence polarization with
    diffuse functions"
    "``...``", ""

- **xtb**

  - the xTB method, if xTB calculation is to be carried out. The ``basis``
    keyword will be ignored if ``xtb`` is specified.

  .. csv-table::
    :widths: 1, 3

    "``gfn0``", "The GFN0-xTB method"
    "``gfn1``", "The GFN1-xTB method"
    "``gfn2``", "The GFN2-xTB method"

- **xcfun**

  - the exchange-correlation functional, if DFT calculation is to be carried
    out

  .. csv-table::
    :widths: 1, 3

    "``SLDA``", "Local density exchange-correlation functional"
    "``BLYP``", "Becke Lee-Yang-Parr exchange-correlation functional"
    "``B3LYP``", "Becke three-parameter Lee-Yang-Parr hybrid
    exchange-correlation functional"
    "``...``", ""

- **grid_level**

  - the accuracy level of DFT grid
  - default: 4

- **potfile**

  - the name of the potential file, if polarizable embedding calculation is to
    be carried out

- **use_split_comm**

  - use split communicators for ERI/DFT/PE
  - default: no

The molecule group
------------------

- **charge**

  - net charge
  - default: 0

- **multiplicity**

  - spin multiplicity
  - default: 1

- **xyz**

  - xyz string, multiple lines

- **units**

  - unit of coordinates in xyz string

  .. csv-table::
    :widths: 1, 3

    "``angstrom``", "Angstrom (default)"
    "``au``", "Atomic unit"
    "``bohr``", "Atomic unit"

The scf group
-------------

- **max_iter**

  - maximum number of iterations
  - default: 50

- **conv_thresh**

  - convergence threshold for SCF
  - default: 1.0e-6

- **eri_thresh**

  - screening threshold for electron repulsion integrals
  - default: 1.0e-12

- **restart**

  - restart from checkpoint file if possible
  - default: yes

- **checkpoint_file**

  - name of the checkpiont file
  - default: <output_file_name>.scf.h5

- **timing**

  - prints timing for SCF iterations
  - default: no

The response group
------------------

- **property**

  - the response property to be calculated

  .. csv-table::
    :widths: 1, 3

    "``polarizability``", "Electric dipole polarizability"
    "``absorption``", "UV-Vis absorption spectrum"
    "``absorption (cpp)``", "Absorption spectrum using complex polarization
    propagator"
    "``ecd (cpp)``", "Electronic circular dichroism using complex polarization
    propagator"

- **frequencies**

  - frequencies for polarizability or CPP calculations, in atomic unit
  - format: start-stop (step)
  - example::

      frequencies: 0.1-0.2 (0.01)

  - default: 0

- **nstates**

  - number of excited states for UV-Vis absorption
  - default: 3

- **tamm_dancoff**

  - use Tamm--Dancoff approximation for UV-Vis absorption
  - default: no

- **max_iter**

  - maximum number of iterations
  - default: 150

- **conv_thresh**

  - convergence threshold for response calculation
  - default: 1.0e-4

- **restart**

  - restart from checkpoint file if possible
  - default: yes

- **checkpoint_file**

  - name of the checkpiont file
  - default: <output_file_name>.rsp.h5

- **timing**

  - prints timing for response iterations
  - default: no

The mp2 group
-------------

- **conventional**

  - use conventional O(N\ :sup:`5`) algorithm for integral transformation
  - default: no

The visualization group
-----------------------

- **grid**

  - number of grid points in three dimensions
  - default: 80,80,80

- **cubes**

  - densities or orbitals for cube files
  - example::

      cubes: density(alpha), mo(homo)

- **files**

  - name of the cube files to be generated
  - example::

      files: density.cube, homo.cube

The exciton group
-----------------------

- **fragments**

  - number of fragments with the same number of atoms
  - example (2 fragments with *m* atoms per fragment and 3 fragments with *n*
    atoms per fragment, where *m* and *n* are specified by the
    ``atoms_per_fragment`` keyword)::

      fragments: 2, 3

- **atoms_per_fragment**

  - number of atoms in each fragment
  - example (18 atoms per fragment in the first group of fragments, 26 atoms
    per fragment in the second group of fragments)::

      atoms_per_fragment: 18, 26

- **nstates**

  - number of locally excited (LE) states in each fragment
  - default: 3

- **ct_nocc**

  - number of occupied oribtals to be involved in charge-transfer (CT) excited
    states
  - default: 0

- **ct_nvir**

  - number of virtual oribtals to be involved in charge-transfer (CT) excited
    states
  - default: 0

The optimize group
-----------------------

- **coordsys**

  - the coordinate system

  .. csv-table::
    :widths: 1, 3

    "``tric``", "Translation-rotation internal coordinates (default)"
    "``cart``", "Cartesian coordinates"
    "``prim``", "Primitive (a.k.a redundant) coordinates"
    "``dlc``", "Delocalized internal coordinates"
    "``hdlc``", "Hybrid delocalized internal coordinates"
