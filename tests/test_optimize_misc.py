from pathlib import Path
from copy import deepcopy

from mpi4py import MPI
import numpy as np
import h5py

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.optimizationengine import OptimizationEngine
from veloxchem.optimizationdriver import OptimizationDriver


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
    def run_unrestricted_opt(molecule,
                             basis,
                             xcfun='hf',
                             filename=None,
                             constraints=None,
                             first_hessian=False,
                             last_hessian=False):

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun
        scf_drv.filename = filename

        opt_drv = OptimizationDriver(scf_drv)
        opt_drv.constraints = constraints
        if first_hessian and last_hessian:
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
                                                         first_hessian=True,
                                                         last_hessian=True)

        final_h5_file = f'{filename}.h5'

        if opt_drv.rank == mpi_master():
            hf = h5py.File(final_h5_file)
            assert 'vib' in hf
            assert 'vib_frequencies' in hf['vib']
            assert 'ir_intensities' in hf['vib']
            ref_vib_freqs = np.array(
                [862.98, 1752.54, 1752.54, 3368.44, 3536.17, 3536.18])
            ref_ir_intens = np.array(
                [163.2794, 0.0232, 0.0232, 2.6263, 0.6178, 0.6178])
            assert np.max(
                np.abs(np.array(hf['vib/vib_frequencies']) -
                       ref_vib_freqs)) < 0.1
            assert np.max(
                np.abs(np.array(hf['vib/ir_intensities']) -
                       ref_ir_intens)) < 0.005
            hf.close()
