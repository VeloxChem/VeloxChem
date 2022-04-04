from pathlib import Path
import numpy as np
import h5py

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.orbitalviewer import OrbitalViewer


class TestOrbitalViewer:

    def test_orbital_viewer(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'h2se.inp')
        h5file = str(here / 'inputs' / 'h2se.orbv.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if scf_drv.rank == mpi_master():
            mo_coefs = scf_drv.mol_orbs.alpha_to_numpy()
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
