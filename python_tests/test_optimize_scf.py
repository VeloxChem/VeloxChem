from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.veloxchemlib import bohr_in_angstroms
from veloxchem.mpitask import MpiTask
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.optimizationdriver import OptimizationDriver


class TestOptimizeSCF:

    def run_opt(self, inpfile, basis_label, ref_coords):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None
        task.input_dict['method_settings']['basis'] = basis_label.upper()

        if is_mpi_master(task.mpi_comm):
            task.ao_basis = MolecularBasis.read(task.molecule,
                                                basis_label.upper())
        else:
            task.ao_basis = MolecularBasis()
        task.ao_basis.broadcast(task.mpi_rank, task.mpi_comm)

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        grad_drv = ScfGradientDriver(scf_drv, task.mpi_comm, task.ostream)
        opt_drv = OptimizationDriver(grad_drv, 'SCF')
        opt_drv.update_settings({
            'coordsys': 'tric',
            'filename': task.input_dict['filename'],
        })
        opt_mol = opt_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if is_mpi_master(task.mpi_comm):
            opt_coords = opt_mol.get_coordinates()
            assert np.max(np.abs(opt_coords - ref_coords)) < 1.0e-6

            inpfile = Path(inpfile)
            optfile = Path(str(inpfile.with_name(inpfile.stem)) + '_optim.xyz')
            if optfile.is_file():
                optfile.unlink()

    def test_nh3(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'nh3.inp')

        basis_label = 'sto-3g'

        ref_coords = np.array([
            [-1.941114808233, 1.588257603470, -0.020632249491],
            [-1.965745778628, 2.619171893222, 0.031449914839],
            [-2.496768727202, 1.280798298868, 0.793505674044],
            [-2.536814090843, 1.362016810479, -0.833073376431],
        ]) / bohr_in_angstroms()

        self.run_opt(inpfile, basis_label, ref_coords)
