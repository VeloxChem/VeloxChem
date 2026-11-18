"""Graphite builder: the carbon electrode only (geometry, surface features, force field).

builder.GraphiteBuilder      the class to use: slab / flake / nanotube, steps, defects,
                             oxygen groups, carbon force field, standalone export
graphite                     lattice generators (make_slab, make_flake, make_nanotube)
surface                      steps, vacancies, Stone-Wales, oxygen groups
surface_templates            QM parametrisation of the oxygen groups (run on the cluster)
forces                       carbon terms: mobility, group terms, NBFIX correction
export                       carbon-only OpenMM system and files
params                       carbon parameter registry with literature status

Solvent, solute, field and FEP are added on top: veloxchem.electrode.build_system.
"""
from .builder import GraphiteBuilder
from .export import carbon_system, write_carbon
from .forces import add_carbon_mobility, add_nbfix_correction, nonbonded_force
from .graphite import (CarbonModel, GraphiteSlab, graphite_forcefield_xml, make_flake,
                       make_nanotube, make_slab, residue_names)
from .params import CARBON_SETS, FORCE_GROUPS, SURFACE_GROUPS
from .surface import add_defects, add_groups, add_step

__all__ = ['GraphiteBuilder', 'CarbonModel', 'GraphiteSlab', 'make_slab', 'make_flake',
           'make_nanotube', 'add_step', 'add_defects', 'add_groups', 'carbon_system',
           'write_carbon', 'add_carbon_mobility', 'add_nbfix_correction', 'nonbonded_force',
           'graphite_forcefield_xml', 'residue_names', 'CARBON_SETS', 'FORCE_GROUPS',
           'SURFACE_GROUPS']
