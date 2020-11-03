The jobs group
--------------

Specifies type of calculation.

List of keywords:

- task

  - string
  - scf | response | pulses | mp2 | visualization | loprop

Example::

  @jobs
  task: scf
  @end

The method settings group
-------------------------

Specifies method of calculation.

List of keywords:

- xcfun

  - string
  - slda | blyp | b3lyp | bhandh | bhandhlyp

- basis

  - string
  - cc-pvdz | aug-cc-pvdz | def2-svp | def2-svpd | sadlej-pvtz | ...

Example::

  @method settings
  xcfun: b3lyp
  basis: def2-svp
  @end

The molecule group
------------------

Specifies the molecule.

List of keywords:

- charge

  - integer
  - default: 0

- multiplicity

  - integer
  - default: 1

- units

  - string
  - angs | au
  - default: angs

- xyz

  - multiple lines of string

Example::

  @molecule
  charge: 0
  multiplicity: 1
  units: au
  xyz:  
  O   0.0   0.0   0.0
  H   0.0   1.4   1.1
  H   0.0  -1.4   1.1
  @end 

