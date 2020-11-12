Details of the input keywords
=============================

The jobs group
--------------

- **task**

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

- **xcfun**

  .. csv-table::
    :widths: 1, 3

    "``SLDA``", "Local density exchange-correlation functional"
    "``BLYP``", "Becke Lee-Yang-Parr exchange-correlation functional"
    "``B3LYP``", "Becke three-parameter Lee-Yang-Parr hybrid exchange-correlation functional"
    "``...``", ""

- **basis**

  .. csv-table::
    :widths: 1, 3

    "``cc-pVDZ``", "Dunning basis set, double-zeta"
    "``aug-cc-pVDZ``", "Dunning basis set, double-zeta with diffuse functions"
    "``def2-SVP``", "Karlsruhe basis set, split valence polarization"
    "``def2-SVPD``", "Karlsruhe basis set, split valence polarization with diffuse functions"
    "``...``", ""

The molecule group
------------------

- **charge**

  - net charge, default 0

- **multiplicity**

  - spin multiplicity, default 1

- **units**

  .. csv-table::
    :widths: 1, 3

    "``angs``", "Angstroms (default)"
    "``au``", "Atomic unit"

- **xyz**

  - xyz string (multiple lines)
