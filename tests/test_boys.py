import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.localizationdriver import LocalizationDriver


@pytest.mark.solvers
class TestBoys:

    def test_run_boys(self, tol=1e-6):

        xyz_string = """12
        D6h
        C      1.396792    0.000000    0.000000
        C      0.698396    1.209951    0.000000
        C     -0.698396    1.209951    0.000000
        C     -1.396792    0.000000    0.000000
        C     -0.698396   -1.209951    0.000000
        C      0.698396   -1.209951    0.000000
        H      2.479452    0.000000    0.000000
        H      1.239726    2.146842    0.000000
        H     -1.239726    2.146842    0.000000
        H     -2.479452    0.000000    0.000000
        H     -1.239726   -2.146842    0.000000
        H      1.239726   -2.146842    0.000000
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_res = scf_drv.compute(mol, bas)

        local_scf_res = scf_res.copy()
        localization_drv = LocalizationDriver()
        C_loc = localization_drv.compute(mol, bas, local_scf_res, list(range(6)))

        if scf_drv.rank == mpi_master():
            assert np.max(np.abs(C_loc.T @ scf_res['S'] @ C_loc - np.eye(C_loc.shape[1]))) < tol
            assert np.max(np.abs(C_loc[:, :6])) > 0.98
