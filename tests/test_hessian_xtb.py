from pathlib import Path
import numpy as np
import pytest
import h5py

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule
from veloxchem.xtbdriver import XtbDriver
from veloxchem.xtbhessiandriver import XtbHessianDriver


class TestXtbHessianDriver:

    @pytest.mark.skipif(not XtbDriver.is_available(),
                        reason='tblite-python not available')
    def test_xtb_hessian_driver(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_xtb.inp')
        h5file = str(here / 'data' / 'water_hessian_xtb.h5')

        task = MpiTask([inpfile, None])

        xtb_method = 'gfn2'

        xtb_drv = XtbDriver(task.mpi_comm, task.ostream)
        xtb_drv.ostream.mute()

        xtb_drv.set_method(xtb_method.lower())
        xtb_drv.compute(task.molecule)

        xtb_hessian_drv = XtbHessianDriver(xtb_drv)
        xtb_hessian_drv.ostream.mute()
        xtb_hessian_drv.compute(task.molecule)

        if task.mpi_rank == mpi_master():

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            hf.close()

            diff_hessian = np.max(np.abs(xtb_hessian_drv.hessian - ref_hessian))

            assert diff_hessian < 1.0e-5

        task.finish()

    @pytest.mark.skipif(not XtbDriver.is_available(),
                        reason='tblite-python not available')
    def test_dipole_gradient_of_charged_molecule(self):

        molecule = Molecule.read_xyz_string("""4
        hydronium
        O   0.000000   0.000000   0.100000
        H   0.000000   0.940000  -0.250000
        H   0.810000  -0.470000  -0.250000
        H  -0.810000  -0.470000  -0.250000
        """)
        molecule.set_charge(1)

        xtb_drv = XtbDriver()
        xtb_drv.ostream.mute()
        xtb_drv.compute(molecule)

        xtb_hessian_drv = XtbHessianDriver(xtb_drv)
        xtb_hessian_drv.ostream.mute()
        xtb_hessian_drv.compute(molecule)

        # Moving the whole molecule by t changes the dipole moment of a
        # molecule with charge q by q t when the origin is fixed, so the
        # dipole gradient summed over the atoms is q times the unit matrix.
        if xtb_drv.rank == mpi_master():
            dipole_gradient = xtb_hessian_drv.dipole_gradient.reshape(3, -1, 3)
            assert np.max(np.abs(dipole_gradient.sum(axis=1) -
                                 np.eye(3))) < 1.0e-4
