from pathlib import Path
import numpy as np
import h5py
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.orbitalviewer import OrbitalViewer
from veloxchem.visualizationdriver import VisualizationDriver
from veloxchem.cubicgrid import CubicGrid


class TestOrbitalViewer:

    def test_orbital_viewer(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2se.inp')
        h5file = str(here / 'data' / 'h2se.orbv.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis)

        if scf_drv.rank == mpi_master():
            mo_coefs = scf_drv.molecular_orbitals.alpha_to_numpy()
            homo = (task.molecule.number_of_alpha_occupied_orbitals(
                task.ao_basis) - 1)

            orb_viewer = OrbitalViewer()
            orb_viewer.initialize(task.molecule, task.ao_basis)
            mo_values = orb_viewer.compute_orbital(mo_coefs, homo)

            hf = h5py.File(h5file)
            ref_mo_values = np.array(hf.get('cube_values'))
            hf.close()

            if np.vdot(mo_values[0, 0, :], ref_mo_values[0, 0, :]) < 0.0:
                mo_values *= -1.0

            assert np.max(np.abs(mo_values - ref_mo_values)) < 1.0e-5

        task.finish()

    def test_orbital_viewer_interpolate(self):

        molstr = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        mol = Molecule.read_molecule_string(molstr, units='au')
        bas = MolecularBasis.read(mol, 'aug-cc-pvdz', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results_not_used = scf_drv.compute(mol, bas)

        mol_orbs = scf_drv.molecular_orbitals
        mol_orbs = mol_orbs.broadcast(scf_drv.comm, root=mpi_master())

        orbviewer = OrbitalViewer()
        orbviewer.interpolate = True
        orbviewer.initialize(mol, bas)
        orbviewer._i_orb = mol.number_of_alpha_occupied_orbitals(bas) - 1
        orbviewer._mo_coefs = mol_orbs.alpha_to_numpy()
        orbital = orbviewer.compute_orbital(orbviewer._mo_coefs,
                                            orbviewer._i_orb)

        vis_drv = VisualizationDriver()
        cubic_grid = CubicGrid(orbviewer.origin, orbviewer.stepsize,
                               orbviewer.npoints)
        vis_drv.compute(cubic_grid, mol, bas, orbviewer._mo_coefs,
                        orbviewer._i_orb)

        if scf_drv.rank == mpi_master():
            orbital_ref = cubic_grid.values_to_numpy()

            max_rel_diff = 0.0
            for i in range(orbital_ref.shape[0]):
                for j in range(orbital_ref.shape[1]):
                    for k in range(orbital_ref.shape[2]):
                        if abs(orbital_ref[i, j, k]) < 1e-6:
                            assert abs(orbital[i, j, k] -
                                       orbital_ref[i, j, k]) < 1e-6
                        else:
                            rel_diff = abs(orbital[i, j, k] /
                                           orbital_ref[i, j, k] - 1.0)
                            max_rel_diff = max(rel_diff, max_rel_diff)
            assert max_rel_diff < 0.05

    def test_orbital_viewer_initialize_with_atom_centers_and_ghost_atom(self):

        mol = Molecule.read_molecule_string("""
            Bq_N  0.0  0.0  0.0
            H     0.0  0.0  1.0
            H     0.0  0.0 -1.0
        """,
                                            units='au')
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        orbviewer = OrbitalViewer()
        orbviewer.atom_centers = [1, 2]
        orbviewer.initialize(mol, bas)

        assert tuple(orbviewer._atomnr) == (6, 0, 0)
        assert orbviewer.origin[2] == pytest.approx(-4.0)
        assert orbviewer.npoints[2] > 0

    def test_compute_orbital_handles_boundary_clipping_and_discard(self):

        orbviewer = OrbitalViewer()
        orbviewer.origin = np.array([0.0, 0.0, 0.0])
        orbviewer.stepsize = np.array([1.0, 1.0, 1.0])
        orbviewer.npoints = (3, 3, 3)
        orbviewer._atom_origin = np.array([-1.0, -1.0, -1.0])
        orbviewer._atom_npoints = (2, 2, 2)
        orbviewer._coords = np.array([[0.4, 0.4, 0.4], [5.0, 5.0, 5.0],
                                      [1.0, 1.0, 1.0]])
        orbviewer._atomnr = np.array([0, 0, 0])
        orbviewer._ao_to_atom = [(0, '', 0), (1, '', 0), (2, '', 0)]
        orbviewer._ao_dict = {(0, ''): [np.ones((2, 2, 2))]}

        orbital = np.array([[1.0], [2.0], [0.001]])

        values = orbviewer.compute_orbital(orbital, 0)

        assert values.shape == (3, 3, 3)
        assert values[0, 0, 0] == pytest.approx(1.0)
        assert np.count_nonzero(values) == 1

        orbviewer.interpolate = True
        orbviewer._atom_npoints = (1, 1, 1)
        orbviewer._atom_origin = np.array([-0.5, -0.5, -0.5])
        orbviewer._coords = np.array([[0.2, 0.2, 0.2], [4.0, 4.0, 4.0]])
        orbviewer._ao_to_atom = [(0, '', 0), (1, '', 0)]
        orbviewer._ao_dict = {(0, ''): [np.array([[[2.0]]])]}

        interp_values = orbviewer.compute_orbital(np.array([[1.0], [1.0]]), 0)

        assert interp_values[0, 0, 0] == pytest.approx(2.0 * 0.7**3)
        assert np.count_nonzero(interp_values) == 1

    def test_get_orbital_cube_str(self):

        mol = Molecule.read_molecule_string("""
            H  0.0  0.0  -0.7
            H  0.0  0.0   0.7
        """,
                                            units='au')
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.compute(mol, bas)

        mol_orbs = scf_drv.molecular_orbitals.broadcast(scf_drv.comm,
                                                        root=mpi_master())

        orbviewer = OrbitalViewer()
        orbviewer.initialize(mol, bas)
        orbviewer._i_orb = 0

        cube_str = orbviewer._get_orbital_cube_str(
            orbviewer.compute_orbital(mol_orbs.alpha_to_numpy(), 0), '(alpha)')

        cube_lines = cube_str.splitlines()

        assert cube_lines[0] == 'Cube file generated by VeloxChem'
        assert cube_lines[1] == 'Orbital 1 (alpha)'
        assert cube_lines[2].strip().endswith('1')
        assert cube_lines[3].split()[0] == str(orbviewer.npoints[0])
        assert cube_lines[4].split()[0] == str(orbviewer.npoints[1])
        assert cube_lines[5].split()[0] == str(orbviewer.npoints[2])

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_invalid_orbital_color_scheme_raises_assertion(self):

        orbviewer = OrbitalViewer()
        orbviewer.orbital_color_scheme = 'invalid'

        with pytest.raises(
                AssertionError,
                match="orbital_color_scheme must be 'default' or 'alternative'"
        ):
            orbviewer._get_orbital_colors()
