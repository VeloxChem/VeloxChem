import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.xtbdriver import XtbDriver


class TestXtbDriver:

    def run_xtb(self, xyz_string, charge, multiplicity):

        molecule = Molecule.read_xyz_string(xyz_string)
        molecule.set_charge(charge)
        molecule.set_multiplicity(multiplicity)

        xtb_drv = XtbDriver()
        xtb_drv.ostream.mute()
        xtb_drv.compute(molecule)

        return xtb_drv

    @pytest.mark.skipif(not XtbDriver.is_available(),
                        reason='tblite-python not available')
    def test_charged_molecule(self):

        xtb_drv = self.run_xtb(
            """4
            hydronium
            O   0.000000   0.000000   0.100000
            H   0.000000   0.940000  -0.250000
            H   0.810000  -0.470000  -0.250000
            H  -0.810000  -0.470000  -0.250000
            """, 1, 1)

        if xtb_drv.rank == mpi_master():
            partial_charges = xtb_drv._xtb_res.get('charges')
            assert abs(np.sum(partial_charges) - 1.0) < 1.0e-10
            assert abs(xtb_drv.get_energy() - (-5.0870560281)) < 1.0e-8

    @pytest.mark.skipif(not XtbDriver.is_available(),
                        reason='tblite-python not available')
    def test_multiplicity(self):

        xyz_string = """2
        oxygen
        O   0.000000   0.000000   0.000000
        O   0.000000   0.000000   1.210000
        """

        singlet_drv = self.run_xtb(xyz_string, 0, 1)
        triplet_drv = self.run_xtb(xyz_string, 0, 3)

        if singlet_drv.rank == mpi_master():
            assert abs(singlet_drv.get_energy() - (-7.9067523570)) < 1.0e-8
            assert abs(triplet_drv.get_energy() - (-7.9041182796)) < 1.0e-8
