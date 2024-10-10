import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lrsolver import LinearResponseSolver


@pytest.mark.solvers
class TestLR:

    def run_lr(self, xcfun_label, ref_rsp_func, tol):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.check_multiplicity()

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseSolver()
        lr_drv.ostream.mute()
        lr_drv.frequencies = [0.05, 0.06]
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            max_diff = 0.0
            for key, val in lr_results['response_functions'].items():
                ref_val = ref_rsp_func[key]
                max_diff = max(max_diff, abs(val - ref_val))
            assert max_diff < tol

    def test_hf(self):

        ref_rsp_func = {
            ('x', 'x', 0.05): -3.05561028,
            ('y', 'x', 0.05): -0.35361889,
            ('z', 'x', 0.05): 0.05536406,
            ('x', 'x', 0.06): -3.06301301,
            ('y', 'x', 0.06): -0.35462235,
            ('z', 'x', 0.06): 0.05546917,
            ('x', 'y', 0.05): -0.35361889,
            ('y', 'y', 0.05): -5.16980438,
            ('z', 'y', 0.05): 0.60122556,
            ('x', 'y', 0.06): -0.35462235,
            ('y', 'y', 0.06): -5.18317763,
            ('z', 'y', 0.06): 0.60255167,
            ('x', 'z', 0.05): 0.05536406,
            ('y', 'z', 0.05): 0.60122556,
            ('z', 'z', 0.05): -6.60300613,
            ('x', 'z', 0.06): 0.05546917,
            ('y', 'z', 0.06): 0.60255167,
            ('z', 'z', 0.06): -6.61957112,
        }

        self.run_lr('hf', ref_rsp_func, 1.0e-6)

    def test_slda(self):

        ref_rsp_func = {
            ('x', 'x', 0.05): -3.22866675,
            ('y', 'x', 0.05): -0.38764460,
            ('z', 'x', 0.05): 0.04947888,
            ('x', 'x', 0.06): -3.24144637,
            ('y', 'x', 0.06): -0.38902139,
            ('z', 'x', 0.06): 0.04937325,
            ('x', 'y', 0.05): -0.38764460,
            ('y', 'y', 0.05): -5.54006288,
            ('z', 'y', 0.05): 0.57715700,
            ('x', 'y', 0.06): -0.38902139,
            ('y', 'y', 0.06): -5.56089547,
            ('z', 'y', 0.06): 0.57715078,
            ('x', 'z', 0.05): 0.04947888,
            ('y', 'z', 0.05): 0.57715700,
            ('z', 'z', 0.05): -6.92248532,
            ('x', 'z', 0.06): 0.04937325,
            ('y', 'z', 0.06): 0.57715078,
            ('z', 'z', 0.06): -6.94349198,
        }

        self.run_lr('slda', ref_rsp_func, 1.0e-5)

    def test_b3lyp(self):

        ref_rsp_func = {
            ('x', 'x', 0.05): -3.20369355,
            ('y', 'x', 0.05): -0.38823232,
            ('z', 'x', 0.05): 0.05404290,
            ('x', 'x', 0.06): -3.21577868,
            ('y', 'x', 0.06): -0.38952401,
            ('z', 'x', 0.06): 0.05400927,
            ('x', 'y', 0.05): -0.38823233,
            ('y', 'y', 0.05): -5.52109035,
            ('z', 'y', 0.05): 0.61083581,
            ('x', 'y', 0.06): -0.38952402,
            ('y', 'y', 0.06): -5.54076712,
            ('z', 'y', 0.06): 0.61130842,
            ('x', 'z', 0.05): 0.05404290,
            ('y', 'z', 0.05): 0.61083581,
            ('z', 'z', 0.05): -6.98116681,
            ('x', 'z', 0.06): 0.05400927,
            ('y', 'z', 0.06): 0.61130842,
            ('z', 'z', 0.06): -7.00210891,
        }

        self.run_lr('b3lyp', ref_rsp_func, 1.0e-5)

    def test_camb3lyp(self):

        ref_rsp_func = {
            ('x', 'x', 0.05): -3.17832038,
            ('y', 'x', 0.05): -0.38098629,
            ('z', 'x', 0.05): 0.05477747,
            ('x', 'x', 0.06): -3.18988932,
            ('y', 'x', 0.06): -0.38221874,
            ('z', 'x', 0.06): 0.05476353,
            ('x', 'y', 0.05): -0.38098629,
            ('y', 'y', 0.05): -5.45343301,
            ('z', 'y', 0.05): 0.61217358,
            ('x', 'y', 0.06): -0.38221875,
            ('y', 'y', 0.06): -5.47225553,
            ('z', 'y', 0.06): 0.61275712,
            ('x', 'z', 0.05): 0.05477748,
            ('y', 'z', 0.05): 0.61217358,
            ('z', 'z', 0.05): -6.91559859,
            ('x', 'z', 0.06): 0.05476353,
            ('y', 'z', 0.06): 0.61275712,
            ('z', 'z', 0.06): -6.93593385,
        }

        self.run_lr('cam-b3lyp', ref_rsp_func, 1.0e-5)

    def test_tpssh(self):

        ref_rsp_func = {
            ('x', 'x', 0.05): -3.22903999,
            ('y', 'x', 0.05): -0.39359565,
            ('z', 'x', 0.05): 0.05279942,
            ('x', 'x', 0.06): -3.24066854,
            ('y', 'x', 0.06): -0.39488579,
            ('z', 'x', 0.06): 0.05276587,
            ('x', 'y', 0.05): -0.39359565,
            ('y', 'y', 0.05): -5.57732792,
            ('z', 'y', 0.05): 0.60469118,
            ('x', 'y', 0.06): -0.39488579,
            ('y', 'y', 0.06): -5.59653871,
            ('z', 'y', 0.06): 0.60516294,
            ('x', 'z', 0.05): 0.05279942,
            ('y', 'z', 0.05): 0.60469118,
            ('z', 'z', 0.05): -7.02397815,
            ('x', 'z', 0.06): 0.05276587,
            ('y', 'z', 0.06): 0.60516294,
            ('z', 'z', 0.06): -7.04445203,
        }

        self.run_lr('tpssh', ref_rsp_func, 1.0e-5)
