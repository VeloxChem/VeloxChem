from pathlib import Path
import numpy as np
import h5py

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
        inpfile = str(here / 'inputs' / 'h2se.inp')
        h5file = str(here / 'inputs' / 'h2se.orbv.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if scf_drv.rank == mpi_master():
            mo_coefs = scf_drv.molecular_orbitals.alpha_to_numpy()
            homo = task.molecule.number_of_alpha_electrons() - 1

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
        scf_results = scf_drv.compute(mol, bas)

        mol_orbs = scf_drv.molecular_orbitals
        mol_orbs.broadcast(scf_drv.rank, scf_drv.comm)

        orbviewer = OrbitalViewer()
        orbviewer.interpolate = True
        orbviewer.initialize(mol, bas)
        orbviewer.i_orb = mol.number_of_alpha_electrons() - 1
        orbviewer.mo_coefs = mol_orbs.alpha_to_numpy()
        orbital = orbviewer.compute_orbital(orbviewer.mo_coefs, orbviewer.i_orb)

        vis_drv = VisualizationDriver()
        cubic_grid = CubicGrid(orbviewer.origin, orbviewer.stepsize,
                               orbviewer.npoints)
        vis_drv.compute(cubic_grid, mol, bas, mol_orbs, orbviewer.i_orb,
                        'alpha')

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
