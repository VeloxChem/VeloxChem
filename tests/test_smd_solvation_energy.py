import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver

try:
    import rdkit
except ImportError:
    pass


@pytest.mark.solvers
class TestSmdSolvation:

    def run_smd_solvation(self, mol, xcfun_label, smd_solvent, ref_epsilon,
                          ref_solv_energy):

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label

        scf_drv.solvation_model = 'smd'
        scf_drv.smd_solvent = smd_solvent
        scf_drv.ostream.mute()
        scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            assert scf_drv.cpcm_drv.epsilon == ref_epsilon

        if scf_drv.rank == mpi_master():
            assert abs(ref_solv_energy[0] - scf_drv.cpcm_drv.cpcm_epol) < 1.0e-4
            assert abs(ref_solv_energy[1] - scf_drv.smd_cds_energy) < 1.0e-3

    @pytest.mark.skipif("rdkit" not in sys.modules,
                        reason="rdkit not available")
    def test_pbe0(self):

        xyz_string = """9

        H  -1.9489660000  0.3684310000  0.0000000000
        O  -1.1900810000 -0.2164680000  0.0000000000
        C   0.0000000000  0.5547550000  0.0000000000
        C   1.1687830000 -0.4037250000  0.0000000000
        H   0.0394610000  1.1978290000  0.8849450000
        H   0.0394610000  1.1978290000 -0.8849450000
        H   2.1117090000  0.1417520000  0.0000000000
        H   1.1331390000 -1.0401400000 -0.8830790000
        H   1.1331390000 -1.0401400000  0.8830790000
        """
        mol = Molecule.read_xyz_string(xyz_string)

        self.run_smd_solvation(
            mol,
            xcfun_label='pbe0',
            smd_solvent='water',
            ref_epsilon=78.355,
            ref_solv_energy=[-0.01208379720722, 0.00383188957126])
