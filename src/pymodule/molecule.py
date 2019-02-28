from .veloxchemlib import Molecule
from .veloxchemlib import bohr_in_angstroms


@staticmethod
def _Molecule_read_str(xyzstr, units='angs', charge=0, spinmult=1):

    scale = 1.0 / bohr_in_angstroms()
    if units.lower() in ['au', 'bohr', 'bohrs']:
        scale = 1.0

    atom_labels = []
    x_coords = []
    y_coords = []
    z_coords = []

    for line in xyzstr.strip().split('\n'):
        content = line.split()
        atom_labels.append(content[0])
        x_coords.append(float(content[1]) * scale)
        y_coords.append(float(content[2]) * scale)
        z_coords.append(float(content[3]) * scale)

    mol = Molecule(atom_labels, x_coords, y_coords, z_coords)
    mol.set_charge(charge)
    mol.set_multiplicity(spinmult)
    mol.check_multiplicity()
    mol.check_proximity(0.1)

    return mol


@staticmethod
def _Molecule_read_xyz(xyzfile, charge=0, spinmult=1):

    xyzstr = ""

    with open(xyzfile, 'r') as f_xyz:
        natoms = int(f_xyz.readline().split()[0])
        f_xyz.readline()
        for a in range(natoms):
            xyzstr += f_xyz.readline().strip() + '\n'

    return Molecule.read_str(xyzstr)


@staticmethod
def _Molecule_from_dict(mol_dict):

    xyzstr = mol_dict['xyzstr']

    charge = 0
    if 'charge' in mol_dict.keys():
        charge = int(mol_dict['charge'])

    spinmult = 1
    if 'multiplicity' in mol_dict.keys():
        spinmult = int(mol_dict['multiplicity'])

    units = 'angs'
    if 'units' in mol_dict.keys():
        units = mol_dict['units'].lower()

    return Molecule.read_str(xyzstr, units, charge, spinmult)


Molecule.read_str = _Molecule_read_str
Molecule.read_xyz = _Molecule_read_xyz
Molecule.from_dict = _Molecule_from_dict
