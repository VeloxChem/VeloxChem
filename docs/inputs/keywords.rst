Details of the input keywords
=============================

The jobs group
--------------

The ``jobs`` group specifies the type of calculation with the following
keyword(s).

- task

  - The type of calculation
  - scf | response | pulses | mp2 | visualization | loprop

The method settings group
-------------------------

The ``method settings`` group specifies the method of calculation with the
following keyword(s).

- xcfun

  - The density functional
  - Case insensitive
  - SLDA | BLYP | B3LYP | BHandH | BHandHLYP

- basis

  - The basis set
  - Case insensitive
  - cc-pVDZ | aug-cc-pVDZ | def2-SVP | def2-SVPD | Sadlej-pVTZ | ...

The molecule group
------------------

The ``molecule`` group specifies the molecule with the following keyword(s):

- charge

  - The net charge of the molecule
  - Default: 0

- multiplicity

  - The spin multiplicity of the molecule
  - Default: 1

- units:

  - The unit of molecular geometry input
  - Default: angs (Angstroms)
  - [angs | au]

- xyz

  - the xyz string (multiple lines)

