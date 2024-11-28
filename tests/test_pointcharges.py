from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver


class TestPointCharges:

    @staticmethod
    def get_molecule_and_basis():

        xyz_string = """
        3
        xyz
        O    1.2361419   1.0137761  -0.0612424
        H    0.5104418   0.8944555   0.5514190
        H    1.9926927   1.1973129   0.4956931
        """
        mol = Molecule.read_xyz_string(xyz_string)
        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        return mol, bas

    def test_scf_with_point_charges(self):

        mol, bas = self.get_molecule_and_basis()

        here = Path(__file__).parent
        potfile = str(here / 'data' / 'pe_water.pot')
        vdwfile = str(here / 'data' / 'pe_water.qm_vdw_params.txt')

        for scf_drv in [ScfRestrictedDriver(), ScfUnrestrictedDriver()]:

            scf_drv._point_charges = potfile
            scf_drv._qm_vdw_params = vdwfile
            scf_drv.ostream.mute()

            scf_results = scf_drv.compute(mol, bas)

            if scf_drv.rank == mpi_master():
                assert abs(scf_results['scf_energy'] - (-75.9732334209)) < 1e-9

            grad_drv = ScfGradientDriver(scf_drv)
            grad_drv.compute(mol, bas, scf_results)
            grad = grad_drv.get_gradient()

            ref_grad = np.array([
                [-0.012426949439, -0.009449823622, -0.013955922853],
                [-0.003881905016, -0.002485992070, 0.002890530027],
                [0.010380233206, 0.003716170298, 0.009359249589],
            ])

            assert np.max(np.abs(grad - ref_grad)) < 1e-6
