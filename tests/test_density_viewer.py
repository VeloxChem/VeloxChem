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
