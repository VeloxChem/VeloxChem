from mpi4py import MPI
from pathlib import Path
import h5py
import numpy as np
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

    def run_lr_with_ecp(self, ref_rsp_func, tol):

        xyz_string = """2
        xyz
        Au     0.000000   0.000000   0.000000
        H      0.000000   0.000000   1.550000
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseUnrestrictedSolver()
        lr_drv.ostream.mute()
        lr_drv.frequencies = [0.05]
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            max_diff = 0.0
            for key, val in lr_results['response_functions'].items():
                ref_val = ref_rsp_func[key]
                max_diff = max(max_diff, abs(val - ref_val))
            assert max_diff < tol

    def test_hf_with_ecp(self):

        ref_rsp_func = {}

        for w in [0.0, 0.01, 0.05]:
            for op_a in 'xyz':
                for op_b in 'xyz':
                    ref_rsp_func[(op_a, op_b, w)] = 0.0

        ref_rsp_func[('x', 'x', 0.05)] = -14.00881203
        ref_rsp_func[('y', 'y', 0.05)] = -14.00881203
        ref_rsp_func[('z', 'z', 0.05)] = -24.65918702

        self.run_lr_with_ecp(ref_rsp_func, 1.0e-5)

    def test_checkpoint_restart_and_saved_solutions(self, tmp_path):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()

        filename = str(tmp_path / 'lr_unrest_restart')
        filename = scf_drv.comm.bcast(filename, root=mpi_master())

        scf_drv.filename = filename
        scf_drv.compute(mol, bas)
        scf_drv.restart = True
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseUnrestrictedSolver()
        lr_drv.filename = filename
        lr_drv.ostream.mute()

        lr_drv.frequencies = [0.0, 0.01]
        lr_drv.compute(mol, bas, scf_results)

        lr_drv.restart = True
        lr_drv.frequencies = [0.0, 0.01, 0.05]
        restarted_results = lr_drv.compute(mol, bas, scf_results)
        assert lr_drv.restart is True

        lr_drv.restart = False
        lr_drv.frequencies = [0.0, 0.01, 0.05]
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
                assert 'x_x_0.00000000' in h5file['rsp']
                assert 'z_z_0.05000000' in h5file['rsp']

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_external_rhs(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseUnrestrictedSolver()
        lr_drv.ostream.mute()
        lr_drv.frequencies = [0.0]

        v_grad = {
            ('x', 0.0): (np.zeros(10), np.zeros(9)),
        }

        with pytest.raises(
                AssertionError,
                match='LinearResponseUnrestrictedSolver: not implemented for external rhs'):
            lr_drv.compute(mol, bas, scf_results, v_grad=v_grad)
