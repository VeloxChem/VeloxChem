from .veloxchemlib import VisualizationDriver


@staticmethod
def _gen_grid(molecule, nx=80, ny=80, nz=80):

    # need implementation

    pass


@staticmethod
def _write_cube(molecule, basis, mol_orbs, mo_idx, mo_spin):

    # need implementation

    pass


VisualizationDriver.gen_grid = _gen_grid
VisualizationDriver.write_cube = _write_cube
