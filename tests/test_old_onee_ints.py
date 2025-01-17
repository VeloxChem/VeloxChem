from pathlib import Path
import numpy as np
import h5py

from veloxchem import ElectricDipoleIntegralsDriver
from veloxchem.veloxchemlib import LinearMomentumIntegralsDriver
from veloxchem.veloxchemlib import AngularMomentumIntegralsDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis


class TestOldOneElectronIntegrals:

    def test_old_onee_ints(self):

        molstr = """
        N   -3.710    3.019   -0.037
        H   -3.702    4.942    0.059
        H   -4.704    2.415    1.497
        H   -4.780    2.569   -1.573
        C   -1.621   -5.080    0.444
        H   -0.819   -6.698   -0.465
        H   -3.412   -4.654   -0.393
        H   -0.381   -3.498    0.222
        H   -1.872   -5.468    2.413
        """
        mol = Molecule.read_molecule_string(molstr, units='bohr')

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        dipole_drv = ElectricDipoleIntegralsDriver()
        dipole_mats = dipole_drv.compute(mol, bas)

        linmom_drv = LinearMomentumIntegralsDriver()
        linmom_mats = linmom_drv.compute(mol, bas)

        angmom_drv = AngularMomentumIntegralsDriver()
        angmom_mats = angmom_drv.compute(mol, bas)

        data = {
            'edip_x': dipole_mats.x_to_numpy(),
            'edip_y': dipole_mats.y_to_numpy(),
            'edip_z': dipole_mats.z_to_numpy(),
            'lmom_x': linmom_mats.x_to_numpy(),
            'lmom_y': linmom_mats.y_to_numpy(),
            'lmom_z': linmom_mats.z_to_numpy(),
            'amom_x': angmom_mats.x_to_numpy(),
            'amom_y': angmom_mats.y_to_numpy(),
            'amom_z': angmom_mats.z_to_numpy(),
        }

        here = Path(__file__).parent
        onee_h5_file = str(here / 'data' / 'nh3_ch4_old_onee_ints.h5')

        hf = h5py.File(onee_h5_file, 'r')

        for key in data:
            ref_mat = np.array(hf.get(key))
            assert np.max(np.abs(data[key] - ref_mat)) < 1.0e-12

        hf.close()
