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

        scf_results = scf_drv.compute(task.molecule, task.ao_basis,
                                      task.min_basis)

        density_viewer = DensityViewer()
        density_viewer.initialize(task.molecule, task.ao_basis)
        density_viewer.use_visualization_driver = True
        if scf_drv.rank == mpi_master():
            density_matrix = scf_results['D_alpha'] + scf_results['D_beta']
        else:
            density_matrix = None
        density_matrix = scf_drv.comm.bcast(density_matrix, root=mpi_master())
        density = density_viewer.compute_density(density_matrix)

        if scf_drv.rank == mpi_master():

            hf = h5py.File(h5file, "r")
            ref_den_values = np.array(hf.get('cube-values'))
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
            density_matrix_a = scf_results['D_alpha']
        else:
            density_matrix_a = None
        density_matrix_a = scf_drv.comm.bcast(density_matrix_a,
                                              root=mpi_master())

        density_viewer = DensityViewer()
        density_viewer.grid_density = 6
        density_viewer.atombox_radius = 9
        density_viewer.initialize(mol, bas)
        density_viewer.use_visualization_driver = False
        density = density_viewer.compute_density(density_matrix_a)

        density_viewer.use_visualization_driver = True
        density_ref = density_viewer.compute_density(density_matrix_a)

        if scf_drv.rank == mpi_master():
            max_rel_diff = 0.0
            for i in range(density_ref.shape[0]):
                for j in range(density_ref.shape[1]):
                    for k in range(density_ref.shape[2]):
                        # Note: skip large densities
                        if abs(density_ref[i, j, k]) > 1.0:
                            continue
                        elif abs(density_ref[i, j, k]) < 1e-6:
                            assert abs(density[i, j, k] -
                                       density_ref[i, j, k]) < 1e-6
                        else:
                            rel_diff = abs(density[i, j, k] /
                                           density_ref[i, j, k] - 1.0)
                            max_rel_diff = max(rel_diff, max_rel_diff)
            assert max_rel_diff < 0.1
