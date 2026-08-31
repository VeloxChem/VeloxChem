from pathlib import Path
from copy import deepcopy

from mpi4py import MPI
import numpy as np
import h5py
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.optimizationengine import OptimizationEngine
from veloxchem.optimizationdriver import OptimizationDriver


class _FakeOStream:

    def print_header(self, *args, **kwargs):
        pass

    def print_blank(self, *args, **kwargs):
        pass

    def print_info(self, *args, **kwargs):
        pass

    def flush(self):
        pass


class _FakeGradientDriver:

    def __init__(self, comm):
        self.comm = comm
        self.rank = comm.Get_rank()
        self.ostream = _FakeOStream()

    def compute_energy(self, molecule, *args):
        energy = float(np.sum(molecule.get_coordinates_in_bohr()))
        return self.comm.bcast(energy, root=mpi_master())

    def compute(self, molecule, *args):
        self.comm.bcast(0, root=mpi_master())

    def get_gradient(self):
        return np.zeros((3, 3))


class _FailingGradientDriver(_FakeGradientDriver):

    def compute_energy(self, molecule, *args):
        raise AssertionError('SCF did not converge')


class TestOptimizeMiscellaneous:

    @staticmethod
    def get_ch3_molecule_and_basis():

        xyz_string = """
        4
        ch3
        C   -1.85334300   -0.63945100    1.29623300
        H   -2.40884500   -1.56570200    1.04276400
        H   -2.24160900   -0.22442700    2.24900500
        H   -1.98830700    0.10613700    0.48589200
        """
        molecule = Molecule.read_xyz_string(xyz_string)
        molecule.set_multiplicity(2)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        return molecule, basis

    @staticmethod
    def get_small_molecule():

        xyz_string = """
        3
        small
        O    0.00000000    0.00000000    0.00000000
        H    1.00000000    0.00000000    0.00000000
        H    0.00000000    1.00000000    0.00000000
        """
        return Molecule.read_xyz_string(xyz_string)

    @staticmethod
    def run_unrestricted_opt(molecule,
                             basis,
                             xcfun='hf',
                             filename=None,
                             constraints=None,
                             first_hessian=False,
                             last_hessian=False,
                             hessian_stop=False):

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun
        scf_drv.filename = filename

        opt_drv = OptimizationDriver(scf_drv)
        opt_drv.constraints = constraints
        if hessian_stop:
            opt_drv.hessian = 'stop'
        elif first_hessian and last_hessian:
            opt_drv.hessian = 'first+last'
        elif first_hessian:
            opt_drv.hessian = 'first'
        elif last_hessian:
            opt_drv.hessian = 'last'
        opt_results = opt_drv.compute(molecule, basis)

        return opt_drv, opt_results

    def test_optimization_engine_and_driver_deepcopy(self):

        molecule, basis = self.get_ch3_molecule_and_basis()

        opt_drv, opt_results = self.run_unrestricted_opt(
            molecule, basis, 'b3lyp')

        if opt_drv.rank == mpi_master():
            filename = opt_drv.grad_drv.scf_driver.filename
            for fname in [f'{filename}.h5', f'{filename}_scf.h5']:
                fpath = Path(fname)
                if fpath.is_file():
                    fpath.unlink()

        opt_engine = OptimizationEngine(opt_drv.grad_drv, molecule, basis)

        opt_engine_copy = deepcopy(opt_engine)

        assert opt_engine_copy.molecule == opt_engine.molecule
        assert opt_engine_copy.opt_unparsed_input == opt_engine.opt_unparsed_input
        assert opt_engine_copy.grad_drv.xcfun == opt_engine.grad_drv.xcfun
        assert opt_engine_copy.args[0] == opt_engine.args[0]

        opt_drv_copy = deepcopy(opt_drv)

        assert opt_drv_copy.grad_drv.xcfun == opt_drv.grad_drv.xcfun
        assert opt_drv_copy.grad_drv.ostream is opt_drv.grad_drv.ostream

    def test_read_settings(self, tmp_path):

        molecule, basis = self.get_ch3_molecule_and_basis()

        filename = str(tmp_path / "opt_import_settings")

        comm = MPI.COMM_WORLD
        filename = comm.bcast(filename, root=mpi_master())

        opt_drv, opt_results = self.run_unrestricted_opt(
            molecule,
            basis,
            'b3lyp',
            filename=filename,
            constraints=[
                'freeze distance 1 2',
                'set angle 2 1 3 115.0',
            ])

        checkpoint_file = f'{filename}_scf.h5'

        new_scf_drv = ScfUnrestrictedDriver()
        new_opt_drv = OptimizationDriver(new_scf_drv)
        new_opt_drv.read_settings(checkpoint_file)

        # these should be updated by read_settings
        assert new_opt_drv.grad_drv.scf_driver.xcfun == opt_drv.grad_drv.xcfun
        assert new_opt_drv.grad_drv.xcfun == opt_drv.grad_drv.xcfun
        assert new_opt_drv.constraints == opt_drv.constraints

        # these should not be updated by read_settings
        assert new_opt_drv.grad_drv.scf_driver.filename is None
        assert new_opt_drv.grad_drv.scf_driver.checkpoint_file is None
        assert new_opt_drv.grad_drv.gradient is None

    def test_opt_with_hessian(self, tmp_path):

        molecule, basis = self.get_ch3_molecule_and_basis()

        filename = str(tmp_path / "opt_with_hessian")

        comm = MPI.COMM_WORLD
        filename = comm.bcast(filename, root=mpi_master())

        opt_drv, opt_results = self.run_unrestricted_opt(molecule,
                                                         basis,
                                                         'hf',
                                                         filename=filename,
                                                         first_hessian=False,
                                                         last_hessian=True)

        final_h5_file = f'{filename}.h5'

        if opt_drv.rank == mpi_master():
            hf = h5py.File(final_h5_file)
            assert 'vib' in hf
            assert 'vib_frequencies' in hf['vib']
            assert 'ir_intensities' in hf['vib']
            ref_vib_freqs = np.array(
                [584.46, 1690.07, 1690.07, 3552.39, 3827.10, 3827.10])
            ref_ir_intens = np.array(
                [16.8701, 0.3816, 0.3816, 0.0144, 0.0357, 0.0357])
            assert np.max(
                np.abs(np.array(hf['vib/vib_frequencies']) -
                       ref_vib_freqs)) < 0.1
            assert np.max(
                np.abs(np.array(hf['vib/ir_intensities']) -
                       ref_ir_intens)) < 0.005
            hf.close()

    def test_engine_worker_protocol(self, tmp_path):
        """
        Master-driven engine protocol: the master evaluates a few geometries
        (including a repeated one served from the engine cache), then signals
        'stop'; non-master ranks mirror the evaluations via run_worker and
        return on 'stop'.
        """

        molecule = self.get_small_molecule()

        engine = OptimizationEngine(_FakeGradientDriver(MPI.COMM_WORLD),
                                    molecule)

        if engine.rank == mpi_master():
            dirname = str(tmp_path)
            geometries = [
                np.array([0.0, 0.0, 0.0, 1.0 + s, 0.0, 0.0, 0.0, 1.0 + s,
                          0.0]) for s in (0.0, 0.01, 0.02)
            ]
            results = [engine.calc(coords, dirname) for coords in geometries]
            assert len({r['energy'] for r in results}) == 3
            assert results[0]['gradient'].shape == (9,)

            # repeated geometry: engine cache hit, no collective evaluation
            # on any rank
            cached = engine.calc(geometries[0], dirname)
            assert cached['energy'] == results[0]['energy']

            engine.comm.bcast('stop', root=mpi_master())
            engine.comm.bcast(True, root=mpi_master())
        else:
            engine.run_worker()
            flag = engine.comm.bcast(None, root=mpi_master())
            assert flag is True

    def test_engine_worker_protocol_collective_failure(self, tmp_path):
        """
        A failure raised inside calc_new on all ranks propagates on every
        rank (no deadlock) with the original error message.
        """

        molecule = self.get_small_molecule()

        engine = OptimizationEngine(_FailingGradientDriver(MPI.COMM_WORLD),
                                    molecule)

        coords = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0])

        if engine.rank == mpi_master():
            with pytest.raises(AssertionError, match='SCF did not converge'):
                engine.calc(coords, str(tmp_path))
        else:
            with pytest.raises(AssertionError, match='SCF did not converge'):
                engine.run_worker()

    def test_opt_hessian_exit(self, tmp_path):
        """
        geomeTRIC raises HessianExit for hessian='stop': the driver takes the
        hessian_exit branch (no optimization steps, input geometry is the
        final geometry) and the 'stop'/'flag' broadcasts keep all ranks in
        step.
        """

        molecule, basis = self.get_ch3_molecule_and_basis()

        filename = str(tmp_path / "opt_hessian_exit")

        comm = MPI.COMM_WORLD
        filename = comm.bcast(filename, root=mpi_master())

        opt_drv, opt_results = self.run_unrestricted_opt(molecule,
                                                         basis,
                                                         'hf',
                                                         filename=filename,
                                                         hessian_stop=True)

        if opt_drv.rank == mpi_master():
            # HessianExit: no optimization step is taken, so the final
            # geometry is the input geometry
            opt_mol = Molecule.read_xyz_string(opt_results['final_geometry'])
            assert np.max(
                np.abs(opt_mol.get_coordinates_in_bohr() -
                       molecule.get_coordinates_in_bohr())) < 1.0e-6
            assert set(opt_results.keys()) == {'final_geometry'}
