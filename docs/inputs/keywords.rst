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

    "``angs``", "Angstroms (default)"
    "``au``", "Atomic unit"

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
  - format: start-end (step)
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
