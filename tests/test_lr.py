from mpi4py import MPI
from pathlib import Path
import h5py
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.errorhandler import VeloxChemError


@pytest.mark.solvers
class TestLR:

    def run_lr(self, xcfun_label, ref_rsp_func, tol, ri_coulomb=False):

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
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.ri_coulomb = ri_coulomb
        scf_drv.acc_type = 'l2_c2diis'
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

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_openshell_molecule(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(3)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        # empty scf results just for testing
        scf_results = {}

        lr_drv = LinearResponseSolver()
        lr_drv.ostream.mute()

        with pytest.raises(
                VeloxChemError,
                match="LinearResponseSolver: not implemented for unrestricted case"):
            lr_results_not_used = lr_drv.compute(mol, bas, scf_results)

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

    def test_ri_blyp(self):

        ref_rsp_func = {
            ('x', 'x', 0.05): -3.28252412,
            ('y', 'x', 0.05): -0.40522649,
            ('z', 'x', 0.05): 0.05390479,
            ('x', 'x', 0.06): -3.29747011,
            ('y', 'x', 0.06): -0.40664672,
            ('z', 'x', 0.06): 0.05379217,
            ('x', 'y', 0.05): -0.40522650,
            ('y', 'y', 0.05): -5.69997023,
            ('z', 'y', 0.05): 0.61928204,
            ('x', 'y', 0.06): -0.40664672,
            ('y', 'y', 0.06): -5.72322129,
            ('z', 'y', 0.06): 0.61924903,
            ('x', 'z', 0.05): 0.05390479,
            ('y', 'z', 0.05): 0.61928204,
            ('z', 'z', 0.05): -7.18182779,
            ('x', 'z', 0.06): 0.05379217,
            ('y', 'z', 0.06): 0.61924903,
            ('z', 'z', 0.06): -7.20519718,
        }

        self.run_lr('blyp', ref_rsp_func, 1.0e-5, ri_coulomb=True)

    def run_lr_with_ecp(self, ref_rsp_func, tol):

        xyz_string = """3
        xyz
        Hg   0.00  0.00  0.00
        Cl   0.00  0.00  2.35
        Cl   0.00  0.10 -2.40
        """
        mol = Molecule.read_xyz_string(xyz_string)

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseSolver()
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            max_diff = 0.0
            for key, val in lr_results['response_functions'].items():
                ref_val = ref_rsp_func[key]
                max_diff = max(max_diff, abs(val - ref_val))
            assert max_diff < tol

    def test_hf_with_ecp(self):

        ref_rsp_func = {
            ('x', 'x', 0.0): -27.35540784,
            ('y', 'x', 0.0): 0.00000000,
            ('z', 'x', 0.0): 0.00000000,
            ('x', 'y', 0.0): 0.00000000,
            ('y', 'y', 0.0): -27.40381961,
            ('z', 'y', 0.0): 1.06711071,
            ('x', 'z', 0.0): 0.00000000,
            ('y', 'z', 0.0): 1.06711071,
            ('z', 'z', 0.0): -76.97515361,
        }

        self.run_lr_with_ecp(ref_rsp_func, 1.0e-5)

    def test_checkpoint_restart_and_saved_solutions(self, tmp_path):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()

        filename = str(tmp_path / 'lr_restart')
        filename = scf_drv.comm.bcast(filename, root=mpi_master())

        scf_drv.filename = filename
        scf_drv.compute(mol, bas)
        scf_drv.restart = True
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseSolver()
        lr_drv.filename = filename
        lr_drv.ostream.mute()

        lr_drv.frequencies = [0.05]
        lr_drv.compute(mol, bas, scf_results)

        lr_drv.restart = True
        lr_drv.frequencies = [0.05, 0.06]
        restarted_results = lr_drv.compute(mol, bas, scf_results)
        assert lr_drv.restart is True

        lr_drv.restart = False
        lr_drv.frequencies = [0.05, 0.06]
        fresh_results = lr_drv.compute(mol, bas, scf_results)
        assert lr_drv.restart is False

        if lr_drv.rank == mpi_master():
            for key, value in restarted_results['response_functions'].items():
                assert value == pytest.approx(fresh_results['response_functions'][key],
                                              abs=1.0e-8)

            solution_file = Path(f'{filename}.h5')
            checkpoint_file = Path(f'{filename}_rsp.h5')

            assert solution_file.is_file()
            assert checkpoint_file.is_file()

            with h5py.File(solution_file, 'r') as h5file:
                assert 'rsp' in h5file
                assert 'x_x_0.05000000' in h5file['rsp']
                assert 'z_z_0.06000000' in h5file['rsp']

    def test_compute_with_external_rhs_returns_solutions_only(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseSolver()
        lr_drv.ostream.mute()
        lr_drv.frequencies = [0.0]

        b_grad = lr_drv.get_prop_grad('electric dipole', 'xyz', mol, bas,
                                      scf_results)
        v_grad = {(op, 0.0): vec for op, vec in zip('xyz', b_grad)}

        results = lr_drv.compute(mol, bas, scf_results, v_grad=v_grad)

        assert 'solutions' in results
        assert 'response_functions' not in results
