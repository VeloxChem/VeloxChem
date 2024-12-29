from pathlib import Path
import numpy as np
import h5py

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.densityviewer import DensityViewer


class TestDensityViewer:

    def test_density_viewer(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2se.inp')
        h5file = str(here / 'data' / 'h2se.denv.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_results = scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if scf_drv.rank == mpi_master():

            density_viewer = DensityViewer()
            density_viewer.initialize(task.molecule, task.ao_basis)
            density_viewer.use_visualization_driver = True
            density_matrix = 2 * scf_results['D_alpha']
            density = density_viewer.compute_density(density_matrix)

            hf = h5py.File(h5file, "r")
            ref_den_values = np.array(hf.get('cube-values'))
            hf.close()

            assert np.max(np.abs(density - ref_den_values)) < 1.0e-6

        task.finish()

    def test_density_viewer_compute_simple(self):
        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2se.inp')
        h5file = str(here / 'data' / 'h2se.denv.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_results = scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if scf_drv.rank == mpi_master():

            density_viewer = DensityViewer()
            density_viewer.initialize(task.molecule, task.ao_basis)
            density_viewer.use_visualization_driver = False
            density_matrix = 2 * scf_results['D_alpha']
            density = density_viewer.compute_density(density_matrix)

            hf = h5py.File(h5file, "r")
            ref_den_values = np.array(hf.get('cube-values-simple'))
            hf.close()

            assert np.max(np.abs(density - ref_den_values)) < 1.0e-6

        task.finish()

    def test_density_viewer_relative_difference(self):

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

        if scf_drv.rank == mpi_master():
            density_matrix = 2 * scf_results['D_alpha']

            density_viewer = DensityViewer()
            density_viewer.grid_density = 6
            density_viewer.atombox_radius = 9
            density_viewer.initialize(mol, bas)
            density_viewer.use_visualization_driver = False
            density = density_viewer.compute_density(density_matrix)

            density_viewer.use_visualization_driver = True
            density_ref = density_viewer.compute_density(density_matrix)

            max_rel_diff = 0.0
            for i in range(density_ref.shape[0]):
                for j in range(density_ref.shape[1]):
                    for k in range(density_ref.shape[2]):
                        if abs(density_ref[i, j, k]) < 1e-6:
                            assert abs(density[i, j, k] -
                                       density_ref[i, j, k]) < 1e-6
                        else:
                            rel_diff = abs(density[i, j, k] /
                                           density_ref[i, j, k] - 1.0)
                            max_rel_diff = max(rel_diff, max_rel_diff)
            assert max_rel_diff < 0.5
