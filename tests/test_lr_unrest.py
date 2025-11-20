import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.lrsolverunrest import LinearResponseUnrestrictedSolver


@pytest.mark.solvers
class TestUnrestrictedLR:

    def run_lr(self, xcfun_label, ref_rsp_func, tol):

        xyz_string = """3
        xyz
        O      0.000000   0.000000   0.117790
        H      0.000000   0.755453  -0.471161
        H      0.000000  -0.755453  -0.471161
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(3)

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseUnrestrictedSolver()
        lr_drv.ostream.mute()
        lr_drv.frequencies = [0.0, 0.01, 0.05]
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            max_diff = 0.0
            for key, val in lr_results['response_functions'].items():
                ref_val = ref_rsp_func[key]
                max_diff = max(max_diff, abs(val - ref_val))
            assert max_diff < tol

    def test_camb3lyp(self):

        ref_rsp_func = {}

        for w in [0.0, 0.01, 0.05]:
            for op_a in 'xyz':
                for op_b in 'xyz':
                    ref_rsp_func[(op_a, op_b, w)] = 0.0

        ref_rsp_func[('x', 'x', 0.0)] = -3.54375468
        ref_rsp_func[('y', 'y', 0.0)] = -66.53060889
        ref_rsp_func[('z', 'z', 0.0)] = -4.79288428

        ref_rsp_func[('x', 'x', 0.01)] = -3.54934293
        ref_rsp_func[('y', 'y', 0.01)] = -67.29003967
        ref_rsp_func[('z', 'z', 0.01)] = -4.79357607

        ref_rsp_func[('x', 'x', 0.05)] = -3.74295389
        ref_rsp_func[('y', 'y', 0.05)] = -93.80105016
        ref_rsp_func[('z', 'z', 0.05)] = -4.81028519

        self.run_lr('cam-b3lyp', ref_rsp_func, 1.0e-5)
