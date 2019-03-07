from .veloxchemlib import Molecule
from .veloxchemlib import bohr_in_angstroms
from .errorhandler import assert_msg_critical

import numpy as np


@staticmethod
def _Molecule_read_str(xyzstr, units='angs'):

    labels = []
    coords = []

    for line in xyzstr.strip().split('\n'):
        if line:
            content = line.split()
            labels.append(content[0])
            coords.append([float(x) for x in content[1:4]])

    return Molecule(labels, coords, units)


@staticmethod
def _Molecule_read_xyz(xyzfile):

    xyzstr = ''

    with open(xyzfile, 'r') as f_xyz:
        natoms = int(f_xyz.readline().split()[0])
        f_xyz.readline()
        for a in range(natoms):
            xyzstr += f_xyz.readline().strip() + '\n'

    return Molecule.read_str(xyzstr, 'angs')


@staticmethod
def _Molecule_from_dict(mol_dict):

    xyzstr = mol_dict['xyzstr']

    units = 'angs'
    if 'units' in mol_dict:
        units = mol_dict['units'].lower()

    charge = 0.0
    if 'charge' in mol_dict:
        charge = float(mol_dict['charge'])

    multiplicity = 1
    if 'multiplicity' in mol_dict:
        multiplicity = int(mol_dict['multiplicity'])

    mol = Molecule.read_str(xyzstr, units)
    mol.set_charge(charge)
    mol.set_multiplicity(multiplicity)
    mol.check_multiplicity()
    mol.check_proximity(0.1)

    return mol


Molecule.read_str = _Molecule_read_str
Molecule.read_xyz = _Molecule_read_xyz
Molecule.from_dict = _Molecule_from_dict
