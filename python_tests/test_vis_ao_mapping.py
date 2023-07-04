import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.visualizationdriver import VisualizationDriver
from veloxchem.cubicgrid import CubicGrid


class TestVisAOMapping:

    def run_vis_ao_mapping(self, molecule, basis, ref_atom_to_ao, ref_ao_info,
                           ref_vals):

        vis_drv = VisualizationDriver()

        if is_mpi_master():

            atom_to_ao = vis_drv.map_atom_to_atomic_orbitals(molecule, basis)
            assert atom_to_ao == ref_atom_to_ao

            # get atomic orbital (AO) information, including:
            # 1) atomic number of the atom
            # 2) angular momentum of the AO
            # 3) spherical harmonic index of the AO
            # 4) basis function index of the AO
            ao_info = vis_drv.get_atomic_orbital_info(molecule, basis)
            assert ao_info == ref_ao_info

            # create cubic grid
            cube_origin = (-0.1, -0.1, -0.1)
            cube_stepsize = (0.1, 0.1, 0.1)
            cube_points = (3, 3, 3)
            grid = CubicGrid(cube_origin, cube_stepsize, cube_points)

            vis_drv.compute_atomic_orbital_for_grid(grid, basis, ao_info[2])
            vals = grid.values_to_numpy()
            assert np.max(np.abs(vals - ref_vals)) < 1.0e-6

    def test_vis_ao_mapping(self):

        mol_str = """
            O     0.000000    0.000000    0.000000
            H     0.000000    0.504284    0.758602
            H     0.000000   -0.504284    0.758602
        """
        molecule = Molecule.read_molecule_string(mol_str, units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        ref_atom_to_ao = [
            [0, 1, 2, 7, 8, 11, 12, 15, 16, 19, 20, 21, 22, 23],
            [3, 4, 9, 13, 17],
            [5, 6, 10, 14, 18],
        ]

        ref_ao_info = [
            [8, 0, 0, 0],
            [8, 0, 0, 1],
            [8, 0, 0, 2],
            [1, 0, 0, 0],
            [1, 0, 0, 1],
            [1, 0, 0, 0],
            [1, 0, 0, 1],
            [8, 1, 0, 0],
            [8, 1, 0, 1],
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [8, 1, 1, 0],
            [8, 1, 1, 1],
            [1, 1, 1, 0],
            [1, 1, 1, 0],
            [8, 1, 2, 0],
            [8, 1, 2, 1],
            [1, 1, 2, 0],
            [1, 1, 2, 0],
            [8, 2, 0, 0],
            [8, 2, 1, 0],
            [8, 2, 2, 0],
            [8, 2, 3, 0],
            [8, 2, 4, 0],
        ]

        ref_vals = np.array([
            [
                [0.25402806, 0.25467744, 0.25402806],
                [0.25467744, 0.25532848, 0.25467744],
                [0.25402806, 0.25467744, 0.25402806],
            ],
            [
                [0.25467744, 0.25532848, 0.25467744],
                [0.25532848, 0.25598119, 0.25532848],
                [0.25467744, 0.25532848, 0.25467744],
            ],
            [
                [0.25402806, 0.25467744, 0.25402806],
                [0.25467744, 0.25532848, 0.25467744],
                [0.25402806, 0.25467744, 0.25402806],
            ],
        ])

        self.run_vis_ao_mapping(molecule, basis, ref_atom_to_ao, ref_ao_info,
                                ref_vals)
