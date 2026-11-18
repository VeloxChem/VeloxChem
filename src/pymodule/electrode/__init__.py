#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
"""Electrode interfaces: graphite slabs with solutes and electrolytes under a uniform field.

Opt-in subpackage (not imported by ``import veloxchem``)::

    from veloxchem.electrode import GraphiteBuilder, build, run_fep
    python -m veloxchem.electrode.build_system my_leg.json

* ``graphite``       carbon models (slab, flake, nanotube), carbon force-field sets, mobility
                     modes, surface features and groups, constant face charge.
* ``system_tools``   solute placement, slab-aware packing, electrolyte registry
                     (H3O+/HSO4-, KOH, ...), OpenMM assembly, uniform field, height restraints,
                     z-barostat, the electrode-aware FEP driver subclass, adsorption free energy.
* ``build_system``   one leg from a JSON config: solute FF, graphite, packing, MD, FEP.

Data files (surface-group templates, HSO4- force field, an example solute) are in
``<VLXDATAPATH>/electrode``.
"""

from .graphite import GraphiteBuilder, make_slab, nonbonded_force
from .system_tools import (adsorption_free_energy, build, build_bulk, electrode_fep_driver,
                           fep_driver_base, physical_fingerprint, run_fep, window_length)

__all__ = ['GraphiteBuilder', 'make_slab', 'nonbonded_force', 'adsorption_free_energy', 'build',
           'build_bulk', 'electrode_fep_driver', 'fep_driver_base', 'physical_fingerprint',
           'run_fep', 'window_length']
