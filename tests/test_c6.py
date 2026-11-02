import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.c6driver import C6Driver


@pytest.mark.solvers
class TestC6:

    def run_c6(self, xcfun_label, ri, ref_c6, tol):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.acc_type = 'l2_diis'
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        lr_drv = C6Driver()
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert abs(lr_results['c6'] - ref_c6) < tol

    @pytest.mark.parametrize('ri', [False, True])
    def test_hf(self, ri):

        # vlxtag: RHF, C6, LR

        xcfun_label = 'hf'
        ref_c6 = {False: 16.852249, True: 16.8498302}[ri]

        self.run_c6(xcfun_label, ri, ref_c6, 2.0e-5)

    @pytest.mark.parametrize('ri', [False, True])
    def test_b3lyp(self, ri):

        # vlxtag: RKS, C6, LR

        xcfun_label = 'b3lyp'
        ref_c6 = {False: 17.421351, True: 17.422050}[ri]

        self.run_c6(xcfun_label, ri, ref_c6, 1.0e-4)
